use std::{cmp::min, ops::DerefMut, path::{Path, PathBuf}, sync::{atomic::{AtomicU64, AtomicU8}, Arc, Mutex}};

use sbwt::{self, BitPackedKmerSorting, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use bitvec::prelude::*;
use rayon::prelude::*;

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

fn colex_compare(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    a.iter().rev().cmp(b.iter().rev())
}


// Returns (len, colex, position in sfs slice)
#[allow(clippy::manual_flatten)] // More readable
fn pick_finimizer(sfs_slice: &[Option<(usize, std::ops::Range<usize>)>]) -> (usize, usize, usize){
    // The finimizer is the shortest unique suffix, with ties broken by colex

    let mut best: (usize, usize) = (usize::MAX, usize::MAX); // (length, colex)
    let mut best_i: usize = usize::MAX;
    for (i,x) in sfs_slice.iter().enumerate() {
        if let Some((len, range)) = x {
            if (*len, range.start) < best {
                best = min(best, (*len, range.start));
                best_i = i;
            } 
        }
    }

    if best == (usize::MAX, usize::MAX){
        panic!("Finimizer not found for kmer");
    }

    (best.0, best.1, best_i)

}

// Returns (v, bitmap), where:
// v is a vec of length |SBWT| such that v[i] is zero if k-mer i does not have a finimizer as a suffix, otherwise
// v[i] is the length of the finimizer that is a suffix of k-mer i.
// bitmap is the concatenation of the bitmap representation of all finimizers in colex order.
pub fn identify_finimizers<P: AsRef<Path> + Send + Sync>(sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray, filenames: &[P], k: usize, n_threads: usize, reverse: bool) -> (Vec<u8>, BitVec<u64, Lsb0>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    let (finimizer_lengths, finimizer_matrix) = thread_pool.install(|| {
        log::info!("Reading input sequences into memory");
        let input_stream = InputStream::new(filenames, reverse); // Read the data into memory again
        let dbs = &input_stream.dbs;
        let sbwt_len = sbwt.n_sets();
        let n_colors = dbs.len();

        let index = StreamingIndex::new(sbwt, lcs);
        let finimizer_lengths: Vec<AtomicU8> = (0..sbwt_len).map(|_| AtomicU8::new(0)).collect();

        let per_color_bitmaps: Vec<BitVec> = dbs.par_iter().enumerate().map(|(db_idx, db)| {
            log::info!("Processing file {}/{}", db_idx + 1, dbs.len());
            let mut this_color_finimizer_marks = bitvec![usize, Lsb0; 0; sbwt_len];
            for rec in db.iter() {
                let seq = rec.seq;
                let sus = index.shortest_freq_bound_suffixes(seq, 1); // Shortest unique suffixes
                for kmer_end in k..seq.len() { // Exclusive end
                    if !has_dna_alphabet(&seq[(kmer_end-k)..kmer_end]) {
                        continue;
                    }

                    let (len, colex, _) = pick_finimizer(&sus[(kmer_end-k) .. kmer_end]);

                    finimizer_lengths[colex].store(len as u8, std::sync::atomic::Ordering::Relaxed);
                    this_color_finimizer_marks.set(colex, true);
                }
            }

            this_color_finimizer_marks
        }).collect();

        // Convert atomic integers to regular integers. This seems to compile into a basically no-op.
        let finimizer_lengths: Vec<u8> = finimizer_lengths.into_iter().map(|x| x.into_inner()).collect();

        // Figure out how many finimizers we have and their average length
        let mut n_finimizers = 0_usize;
        let mut total_length = 0_usize;

        for &x in finimizer_lengths.iter() {
            if x > 0 {
                n_finimizers += 1;
                total_length += x as usize;
            }
        }

        log::info!("Number of finimizers = {}", n_finimizers);
        log::info!("Avg length = {}", total_length as f64 / n_finimizers as f64);

        log::info!("Interleaving individual per-color marks");
        let total_n_bits = n_colors*n_finimizers;
        let total_n_words = total_n_bits.div_ceil(64);
        let finimizer_matrix: Vec<AtomicU64> = (0..total_n_words).map(|_| AtomicU64::new(0)).collect();

        let bar = indicatif::ProgressBar::new((n_colors) as u64);
        (0..n_colors).into_par_iter().for_each(|color| { 
            let mut finimizer_idx = 0_usize;
            for colex in 0..sbwt_len {
                if finimizer_lengths[colex] > 0 { // This position is a finimizer
                    let bit_idx = finimizer_idx*n_colors + color; // We want to set this bit
                    let word_idx = bit_idx / 64;
                    let word_offset = bit_idx % 64;
                    let bit = per_color_bitmaps[color][colex];
                    let mask = (1_u64 << word_offset) * bit as u64;

                    finimizer_matrix[word_idx].fetch_or(mask, std::sync::atomic::Ordering::Relaxed);
                    finimizer_idx += 1;
                }
            }
            bar.inc(1);
            assert_eq!(finimizer_idx, n_finimizers);

        });
        bar.finish();
        log::info!("Color matrix finished");

        // Convert atomic integers to regular integers. This seems to compile into a basically no-op.
        let finimizer_matrix: Vec<u64> = finimizer_matrix.into_iter().map(|x| x.into_inner()).collect();

        (finimizer_lengths, BitVec::<u64, Lsb0>::from_vec(finimizer_matrix))
    });

    (finimizer_lengths, finimizer_matrix)

}

pub fn verify_finimizer_shortest_uniqueness(cf: &ColoredFinimizers, sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, n_threads: usize) {
    // Look up all finimizer strings, verify that their SBWT frequency is 1, and the SBWT frequency of the
    // longest proper suffix is more than .
    let n_finimizers = cf.finimizer_lengths.len();

    // Finimizer length cumulative sum up to but not including the array index 
    let mut f_starts: Vec<usize> = vec![0; n_finimizers + 1];
    for i in 0..n_finimizers {
        f_starts[i+1] = f_starts[i] + cf.finimizer_lengths[i] as usize;
    }

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(|| {
        (0..n_finimizers).par_bridge().for_each(|i| {
            let s = &cf.finimizers_concat[f_starts[i]..f_starts[i+1]];
            assert!(sbwt.search(s).unwrap().len() == 1); // Is unique
            assert!(sbwt.search(&s[1..]).unwrap().len() > 1); // Is not unique
        });
    });
}

pub fn verify_finimizer_colors(cf: &ColoredFinimizers, sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &LcsArray, n_threads: usize, paths: &[PathBuf], reverse: bool) {
    let n_finimizers = cf.finimizer_lengths.len();
    let n_colors = cf.n_colors;
    assert!(cf.n_colors == paths.len());

    let mut finimizers: Vec<&[u8]> = Vec::with_capacity(n_finimizers);
    let mut cumul_finimizer_len = 0_usize;
    for i in 0..n_finimizers {
        let len = cf.finimizer_lengths[i] as usize;
        finimizers.push(&cf.finimizers_concat[cumul_finimizer_len..(cumul_finimizer_len + len)]);
        cumul_finimizer_len += len;
    }


    let k = sbwt.k();
    let index = StreamingIndex::new(sbwt, lcs);

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| {
        paths.par_iter().enumerate().for_each(|(color, path)| {
            log::info!("Processing color {}", path.display());
            let is = InputStream::new(&[path], reverse);
            let db = is.dbs.first().unwrap(); // Also has reverse complements

            // Mark all finimizers (in finimizer order, not sbwt colex order) in this color
            let mut finimizer_marks: BitVec = bitvec![0; n_finimizers];
            for rec in db.iter() {
                let sfs = index.shortest_freq_bound_suffixes(rec.seq, 1);
                for kmer_end in k..rec.seq.len() { // Exclusive end
                    if !has_dna_alphabet(&rec.seq[(kmer_end-k)..kmer_end]) {
                        continue;
                    }
                    let (len, _, mut end) = pick_finimizer(&sfs[(kmer_end-k)..kmer_end]);
                    end += 1; // Exclusive end
                    end += kmer_end-k; // Shift to where we are at rec.seq
                    let finimizer = &rec.seq[(end-len) .. end];

                    match finimizers.binary_search_by(|test| colex_compare(test, finimizer)) {
                        Ok(f_idx) => {
                            finimizer_marks.set(f_idx, true);
                        },
                        Err(_) => panic!("Finimizer {} not found with binary search", String::from_utf8_lossy(finimizer)),
                    }
                }
            }

            // Check that the marks in ColoredFinimizers match our marks
            for f_idx in 0..n_finimizers {
                assert!(cf.color_sets_bitmap[f_idx * n_colors + color] == finimizer_marks[f_idx]);
            }
            log::info!("Color {} OK", path.display());
        });
    });

}

//pub fn verify<P: AsRef<Path> + Send + Sync>(cf: ColoredFinimizers, sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray, filenames: &[P], k: usize, n_threads: usize) {
    // For each color:
    //   Search all finimizers, check that they have the color


// finimizer_lengths should have length |sbwt|, with nonzero entries indicating k-mers who have a finimizer as a suffix
// Returns the concatenation of finimizers in colex order
pub fn get_finimizer_strings(sbwt: &mut sbwt::SbwtIndex<sbwt::SubsetMatrix>, finimizer_lengths: &[u8], n_threads: usize) -> Vec<u8> {
    assert_eq!(sbwt.n_sets(), finimizer_lengths.len());

    log::info!("Building select suppport");
    sbwt.build_select();

    log::info!("Planning input and output ranges");
    let mut n_finimizers = 0_usize;
    let mut total_finimizer_length = 0_usize;
    for &x in finimizer_lengths.iter() {
        if x > 0 {
            n_finimizers += 1;
            total_finimizer_length += x as usize;
        }
    }

    let mut concat = vec![0_u8; total_finimizer_length];

    // Split into mutable segments, one for each thread
    let mut output_segments: Vec<&mut[u8]> = Vec::new();
    let mut input_sbwt_ranges: Vec<std::ops::Range<usize>> = Vec::new();
    let mut end_ptr = 0_usize;
    let n_finimizers_in_range = n_finimizers.div_ceil(n_threads);
    let mut tail = concat.as_mut_slice();
    for _ in 0..n_threads {
        let mut n_finimizers_found = 0;
        let mut total_finimizer_length_found = 0_usize;
        let start_ptr = end_ptr;
        while end_ptr < sbwt.n_sets() && n_finimizers_found < n_finimizers_in_range {
            if finimizer_lengths[end_ptr] > 0 {
                n_finimizers_found += 1;
                total_finimizer_length_found += finimizer_lengths[end_ptr] as usize;
            }
            end_ptr += 1;
        }
        input_sbwt_ranges.push(start_ptr..end_ptr);
        let (a,b) = tail.split_at_mut(total_finimizer_length_found);
        output_segments.push(a);
        tail = b;
    }
    assert!(tail.is_empty()); // Segmentation is complete

    log::info!("Extracting finimizer strings");
    // Fill in segments in parallel
    output_segments.iter_mut().enumerate().par_bridge().for_each(|(t, out_segment)| {
        let mut output_idx = 0_usize;
        for colex in input_sbwt_ranges[t].clone() {
            if finimizer_lengths[colex] > 0 {
                // Write the finimizer to the concatenation
                let kmer = sbwt.access_kmer(colex);
                let finimizer = &kmer[kmer.len() - finimizer_lengths[colex] as usize ..];
                for &c in finimizer {
                    out_segment[output_idx] = c;
                    output_idx += 1;
                }
            }
        }
        assert!(output_idx == out_segment.len());
    });

    concat
}

#[derive(Debug)]
pub struct ColoredFinimizers {
    finimizer_lengths: Vec<u8>, // Lengths of finimizers, in colex order
    finimizers_concat: Vec<u8>, // concatenation of finimizers in colex order in ASCII
    color_sets_bitmap: BitVec<u64, Lsb0>, // concatenation of finimizer color set bit maps, in colex order of finimizers 
    n_colors: usize,
}

impl ColoredFinimizers {

    pub fn serialize<W: std::io::Write>(&self, mut out: &mut W) {
        let n_finimizers = self.finimizer_lengths.len();
        let finimizer_total_length = self.finimizers_concat.len();

        // Write constants
        out.write_all(&(n_finimizers as u64).to_le_bytes()).unwrap();
        out.write_all(&(finimizer_total_length as u64).to_le_bytes()).unwrap();
        out.write_all(&(self.n_colors as u64).to_le_bytes()).unwrap();

        // Write data
        out.write_all(self.finimizer_lengths.as_slice()).unwrap();
        out.write_all(self.finimizers_concat.as_slice()).unwrap();

        let data_bytes: &[u8] = unsafe {std::mem::transmute(self.color_sets_bitmap.as_raw_slice()) };
        assert_eq!(data_bytes.len(), ((n_finimizers * self.n_colors + 63) / 64) * 8);
        out.write_all(data_bytes).unwrap();
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> Self {
        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf).unwrap();
        let n_finimizers = u64::from_le_bytes(buf);

        input.read_exact(&mut buf).unwrap();
        let finimizer_total_length = u64::from_le_bytes(buf);

        input.read_exact(&mut buf).unwrap();
        let n_colors = u64::from_le_bytes(buf);

        let mut finimizer_lengths = vec![0_u8; n_finimizers as usize]; 
        input.read_exact(&mut finimizer_lengths).unwrap();

        let mut finimizers_concat = vec![0_u8; finimizer_total_length as usize]; 
        input.read_exact(&mut finimizers_concat).unwrap();

        let mut bitmap_bytes = vec![0_u8; ((n_finimizers * n_colors + 63) as usize / 64) * 8];
        input.read_exact(bitmap_bytes.as_mut_slice()).unwrap();
        let bitmap_words = unsafe { std::mem::transmute::<&[u8], &[u64]>(&bitmap_bytes) };
        let color_sets_bitmap = BitVec::<u64, Lsb0>::from_slice(bitmap_words); // Todo: this makes a copy

        ColoredFinimizers {finimizer_lengths,  finimizers_concat, color_sets_bitmap, n_colors: n_colors as usize}
    }

    pub fn verify_colex_order(&self) {
        let mut cumul_len = 0_usize;
        for i in 0..self.finimizer_lengths.len() {
            if i > 0 {
                let a = &self.finimizers_concat[(cumul_len - self.finimizer_lengths[i-1] as usize) .. cumul_len];
                let b = &self.finimizers_concat[cumul_len .. cumul_len + self.finimizer_lengths[i] as usize];
                match a.iter().rev().cmp(b.iter().rev()) {
                    std::cmp::Ordering::Equal => panic!("Equal finimizers"),
                    std::cmp::Ordering::Greater => panic!("Non-colex-ordered finimizers"),
                    std::cmp::Ordering::Less => {}, // Ok
                }
            }
            cumul_len += self.finimizer_lengths[i] as usize;
        }
    }
}

struct InputStream {
    dbs: Vec<jseqio::seq_db::SeqDB>,
    cur_db_idx: usize, // Index of the db currently being iterated over
    seq_idx_in_cur_db: usize,
    reverse: bool, // Reverse sequences (not complement) 
}

impl InputStream {
    fn new<P: AsRef<Path>>(filenames: &[P], reverse: bool) -> InputStream {
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

        if reverse {
            log::info!("Reversing sequences");
            for db in dbs.iter_mut() {
                for i in 0..db.sequence_count() {
                    let rec = db.get_mut(i);
                    rec.seq.reverse();
                    rec.qual.map(|q| q.reverse());
                    // Header is not reversed
                }
            }
        }
        Self {dbs, cur_db_idx: 0, seq_idx_in_cur_db: 0, reverse}
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
    pub fn new<P: AsRef<Path> + Send + Sync>(filenames: &[P], k: usize, n_threads: usize, mem_gb: usize, temp_dir: &Path, reverse: bool) -> Self {

        log::info!("Loading {} sequence files (colors) into memory", filenames.len());
        let input_stream = InputStream::new(filenames, reverse);
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
                .mem_gb(mem_gb)
        );
        let (mut sbwt, lcs) = sbwt_builder.run(input_stream);
        let lcs = lcs.unwrap(); // Ok since used build_lcs(true) above

        //eprintln!("Building select (TODO REMOVE THIS");
        //sbwt.build_select(); // DEBUG TODO REMOVE THIS

        let (finimizer_lengths, interleaved_bitmap) = identify_finimizers(&sbwt, &lcs, filenames, k, n_threads, reverse);
        let concat = get_finimizer_strings(&mut sbwt, &finimizer_lengths, n_threads);

        let nonzero_lengths: Vec<u8> = finimizer_lengths.into_iter().filter(|x| *x > 0).collect();

        ColoredFinimizers { finimizer_lengths: nonzero_lengths, finimizers_concat: concat, color_sets_bitmap: interleaved_bitmap, n_colors: num_colors}
    }
}

#[cfg(test)]
mod tests {
    use bitvec::prelude::*;
    #[test] 
    fn test_bit_order(){
        let mut my_bv: Vec<u64> = vec![0; 2];
        my_bv[0] |= (1u64 << 3);
        my_bv[1] |= (1u64 << 4);

        let my_bv = BitVec::<u64, Lsb0>::from_vec(my_bv);
        let mut bv = bitvec![u64, Lsb0; 0; 128];
        bv.set(3, true);
        bv.set(64 + 4, true);

        assert_eq!(bv, my_bv);
    }
}
