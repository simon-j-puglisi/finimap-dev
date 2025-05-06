use std::{cmp::min, io::BufRead, ops::DerefMut, path::{Path, PathBuf}, sync::{atomic::AtomicU8, Arc, Mutex}};

use clap::builder::styling::Color;
use crossbeam::channel::{Receiver, RecvError, Sender};
use sbwt::{self, BitPackedKmerSorting, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix, SubsetSeq};
use bitvec::prelude::*;
use rayon::prelude::*;
use simple_sds_sbwt::{self, raw_vector::AccessRaw};

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

pub(crate) fn is_dna_char(c: u8) -> bool {
    IS_DNA[c as usize]
}

fn has_dna_alphabet(seq: &[u8]) -> bool {
    for &c in seq {
        if !is_dna_char(c) {
            return false;
        }
    }
    true
}


// Returns a vec v of length |SBWT| such that v[i] is zero if k-mer i does not have a finimizer as a suffix, otherwise
// v[i] is the length of the finimizer that is a suffix of k-mer i 
pub fn identify_finimizers<P: AsRef<Path> + Send + Sync>(sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray, filenames: &[P], k: usize, n_threads: usize) -> Vec<u8> {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    let finimizer_lengths = thread_pool.install(|| {
        log::info!("Reading input sequences into memory");
        let mut input_stream = InputStream::new(filenames); // Read the data into memory again
        let dbs = &input_stream.dbs;

        let index = StreamingIndex::new(sbwt, lcs);
        let mut finimizer_lengths: Vec<AtomicU8> = (0..sbwt.n_sets()).map(|_| AtomicU8::new(0)).collect();

        while let Some(seq) = input_stream.stream_next() {
            let sus = index.shortest_freq_bound_suffixes(seq, 1); // Shortes unique suffixes
            for kmer_end in k..seq.len() { // Exclusive end
                if !has_dna_alphabet(&seq[(kmer_end-k)..kmer_end]) {
                    continue;
                }

                // The finimizer is the shortest unique suffix, with ties broken by colex
                let mut best: (usize, usize) = (usize::MAX, usize::MAX); // (length, colex)
                for i in (kmer_end - k)..kmer_end {
                    if let Some((len, range)) = sus[i].clone() {
                        best = min(best, (len, range.start));
                    }
                }

                if best == (usize::MAX, usize::MAX){
                    panic!("Finimizer not found for kmer {}", String::from_utf8_lossy(&seq[(kmer_end-k)..kmer_end]));
                }

                finimizer_lengths[best.1].fetch_min(best.0 as u8, std::sync::atomic::Ordering::Relaxed);
            }
        }
        finimizer_lengths
    });

    // Convert atomic integers to regular integers. This seems to compile into a basically no-op.
    let finimizer_lengths: Vec<u8> = finimizer_lengths.into_iter().map(|x| x.into_inner()).collect();

    finimizer_lengths

    //std::thread::scope(|scope| {
}

#[derive(Debug)]
pub struct ColoredFinimizers {
    kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    distinct_color_sets: BitVec, // Concatenation of distinct color sets
    colex_to_color_set_id: Vec<usize>,
    empty_set: BitVec, // So that we can return a bitslice to an empty set
    n_colors: usize,
}

impl ColoredFinimizers {

    pub fn n_colors(&self) -> usize {
        self.n_colors
    }

    pub fn n_kmers(&self) -> usize {
        self.kmers.n_kmers()
    }

    pub fn get_k(&self) -> usize {
        self.kmers.k()
    }

    pub fn get_color_set(&self, kmer: &[u8]) -> &BitSlice {
        match self.kmers.search(kmer){
            Some(range) => {
                let row = range.start;
                &self.distinct_color_sets[row*self.n_colors..(row+1)*self.n_colors]
            }
            None => &self.empty_set,
        }
    }

    pub fn get_all_color_sets(&self, seq: &[u8]) -> Vec<&BitSlice> {
        let mut sets = Vec::new();
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);

        // Start iterating from index k-1 because those before do not correspond to a full k-mer
        for (match_len, colex_range) in index.matching_statistics(seq)[self.kmers.k()-1 ..].iter() {
            if *match_len == self.kmers.k() {
                assert_eq!(colex_range.len(), 1);
                let id = self.colex_to_color_set_id[colex_range.start];
                sets.push(self.get_color_set_by_id(id));
            } else {
                sets.push(&self.empty_set);
            }
        }
        sets
    }

    fn get_color_set_by_id(&self, id: usize) -> &BitSlice {
        &self.distinct_color_sets[id*self.n_colors..(id+1)*self.n_colors]
    }

    pub fn serialize<W: std::io::Write>(&self, mut out: &mut W) {
        self.kmers.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();

        out.write_all(&(self.n_colors as u64).to_le_bytes()).unwrap();
        out.write_all(&(self.distinct_color_sets.len() as u64).to_le_bytes()).unwrap();

        bincode::serialize_into(&mut out, &self.distinct_color_sets).unwrap();
        bincode::serialize_into(&mut out, &self.colex_to_color_set_id).unwrap();
    }

    pub fn load<R: std::io::Read>(mut input: &mut R) -> Self {
        let kmers = SbwtIndex::<SubsetMatrix>::load(input).unwrap();
        let lcs = sbwt::LcsArray::load(input).unwrap();

        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf).unwrap();
        let n_colors = u64::from_le_bytes(buf);

        input.read_exact(&mut buf).unwrap();
        let _ = u64::from_le_bytes(buf); // Total length of distinct color sets

        let distinct_color_sets: BitVec = bincode::deserialize_from(&mut input).unwrap();
        let colex_to_color_set_id: Vec<usize> = bincode::deserialize_from(&mut input).unwrap(); // Todo: u64 instead of usize

        ColoredFinimizers {kmers, lcs, colex_to_color_set_id, n_colors: n_colors as usize, distinct_color_sets, empty_set: bitvec![0; n_colors as usize]}
    }
}

struct InputStream {
    dbs: Vec<jseqio::seq_db::SeqDB>,
    cur_db_idx: usize, // Index of the db currently being iterated over
    seq_idx_in_cur_db: usize,
}

impl InputStream {
    fn new<P: AsRef<Path>>(filenames: &[P]) -> InputStream {
        let mut dbs: Vec<jseqio::seq_db::SeqDB> = vec![];
        for path in filenames {
            let reader = jseqio::reader::DynamicFastXReader::from_file(path).unwrap();
            let (mut fw, rc) = reader.into_db_with_revcomp().unwrap();

            if fw.sequence_count() == 0 {
                panic!("No sequences found in file {}", path.as_ref().display());
            }

            // Append reverse complement records to the forward database
            for rec in rc.iter() {
                fw.push_record(rec);
            }
            dbs.push(fw);
        }
        Self {dbs, cur_db_idx: 0, seq_idx_in_cur_db: 0}
    }
}

impl SeqStream for InputStream {

    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.cur_db_idx == self.dbs.len() {
            return None; // Done
        }

        // Fetch the next sequence
        let db = &self.dbs[self.cur_db_idx];
        assert!(db.sequence_count() > 0);
        let seq = db.get(self.seq_idx_in_cur_db).seq;

        // Update the "cursor"
        self.seq_idx_in_cur_db += 1;
        if self.seq_idx_in_cur_db == db.sequence_count() {
            self.cur_db_idx += 1;
            self.seq_idx_in_cur_db = 0;
        }
        Some(seq)
    } 
}

fn mark_bits(bv: &mut BitVec, color: usize, num_colors: usize, to_mark: Vec<usize>) {
    for i in to_mark {
        bv.set(i*num_colors + color, true);
    }

}

fn mark_all_kmers_of_seq(bv: Arc<Mutex<BitVec>>, num_colors: usize, color: usize, seq: &[u8], k: usize, mark_buffer_size: usize, index: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>){
    // Search all k-mers
    let mut marking_buffer: Vec<usize> = Vec::new(); // These bits should be marked
    for (len, colex) in index.matching_statistics(seq) {
        if len == k {
            // Full kmer -> set the bit in the color set of the k-mer
            assert!(colex.len() == 1);
            marking_buffer.push(colex.start);
            if marking_buffer.len() == mark_buffer_size {
                mark_bits(bv.lock().unwrap().deref_mut(), color, num_colors, marking_buffer);
                marking_buffer = Vec::new();
            }
        }
    }

    if !marking_buffer.is_empty() { 
        // Mark the rest
        mark_bits(bv.lock().unwrap().deref_mut(), color, num_colors, marking_buffer);
    }
} 

impl ColoredFinimizers {

    #[allow(clippy::type_complexity)]
    pub fn new<P: AsRef<Path> + Send + Sync>(filenames: &[P], k: usize, n_threads: usize, temp_dir: &Path) -> Self {

        log::info!("Loading {} sequence files (colors) into memory", filenames.len());
        let input_stream = InputStream::new(filenames);
        let num_colors = input_stream.dbs.len();
        log::info!("Building SBWT");
        let sbwt_builder = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false) // Already added in the input stream
            .k(k)
            .build_lcs(true)
            .n_threads(n_threads)
            .precalc_length(8)
            .algorithm(BitPackedKmerSorting::new()
                .dedup_batches(true)
                .temp_dir(temp_dir)
        );
        let (sbwt, lcs) = sbwt_builder.run(input_stream);
        let lcs = lcs.unwrap(); // Ok since used build_lcs(true) above

        let finimizer_lengths = identify_finimizers(&sbwt, &lcs, filenames, k, n_threads);
        let mut n_finimizers = 0_usize;
        let mut total_length = 0_usize;

        for x in finimizer_lengths {
            if x > 0 {
                n_finimizers += 1;
                total_length += x as usize;
            }
        }

        eprintln!("Number of finimizers = {}", n_finimizers);
        eprintln!("Avg length = {}", total_length as f64 / n_finimizers as f64);

        /*

        let sbwt_len = sbwt.n_sets();
        let streaming_index_owned = StreamingIndex::new(&sbwt, &lcs);
        let streaming_index = &streaming_index_owned; // Pass by reference into the scope

        log::info!("Reading input sequences back into memory again");
        let input_stream = InputStream::new(filenames); // Read the data into memory again
        let dbs = &input_stream.dbs;

        let color_sets = std::thread::scope(|scope| {

            log::info!("Building colors");

            let work_input_queue: (Sender<(usize, &[u8])>, Receiver<(usize, &[u8])>) = crossbeam::channel::unbounded();

            // Push work to the input queue
            for color in 0..filenames.len() {
                for rec in dbs[color].iter() {
                    work_input_queue.0.send((color, rec.seq)).unwrap(); 
                }
            }
            drop(work_input_queue.0); // Close the channel

            // Spawn worker threads
            let mut worker_handles = Vec::new();
            let color_sets_lock = Arc::new(Mutex::new(bitvec![0; num_colors*sbwt_len])); // Concatenation of color sets
            for thread_id in 0..n_threads {
                let recv_clone = work_input_queue.1.clone();
                let color_sets_lock_clone = color_sets_lock.clone();
                let consumer_handle = scope.spawn(move || {
                    loop {
                        match recv_clone.recv() {
                            Ok((color, seq)) => {
                                mark_all_kmers_of_seq(color_sets_lock_clone.clone(), num_colors, color, seq, k, 100000, streaming_index);
                            },
                            Err(RecvError) => {
                                log::info!("Thread {} finished", thread_id);
                                break;
                            }
                        }
                    }
                });
                worker_handles.push(consumer_handle);
            }

            // Wait for all workers to finish
            for handle in worker_handles {
                handle.join().unwrap();
            }

            // Since we have joined the workers, there should be only one clone of the
            // Arc<Mutex> (the one owned by this thread), so we can consume the lock and return the data.
            Arc::try_unwrap(color_sets_lock).unwrap().into_inner().unwrap()

        }); // End of thread scope 

        // Todo: deduplicate color sets

        let colex_to_color_set_id: Vec<usize> = (0..sbwt_len).collect(); // Identity mapping

        ColoredFinimizers{kmers: sbwt, lcs, distinct_color_sets: color_sets, empty_set: BitVec::new(), colex_to_color_set_id, n_colors: filenames.len()}
        */

        todo!();
    }
}

#[cfg(test)]
mod tests {

}
