use bytemuck::{cast_slice, Pod};
use std::{
    cmp::{max, min},
    collections::HashMap,
    io::Write,
    ops::DerefMut,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicU64, AtomicU8},
        Arc, Mutex,
    },
};

use bitvec::prelude::*;
use rayon::prelude::*;
use sbwt::{
    self, BitPackedKmerSorting, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix,
};

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
fn pick_finimizer(sfs_slice: &[Option<(usize, std::ops::Range<usize>)>]) -> (usize, usize, usize) {
    // The finimizer is the shortest unique suffix, with ties broken by colex

    // The full k-mer should have an existing unique match
    assert!(
        sfs_slice
            .last()
            .expect("Empty slice")
            .as_ref()
            .expect("Last SFS pos is None")
            .1
            .len()
            == 1
    );

    let mut best: (usize, usize) = (usize::MAX, usize::MAX); // (length, colex)
    let mut best_i: usize = usize::MAX;
    for (i, x) in sfs_slice.iter().enumerate() {
        if let Some((len, range)) = x {
            if i + 1 >= *len && (*len, range.start) < best {
                best = min(best, (*len, range.start));
                best_i = i;
            }
        }
    }

    if best == (usize::MAX, usize::MAX) {
        dbg!(sfs_slice);
        panic!("Finimizer not found for kmer");
    }

    (best.0, best.1, best_i)
}

pub fn identify_finimizers_sparse<P: AsRef<Path> + Send + Sync>(
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: &sbwt::LcsArray,
    filenames: &[P],
    k: usize,
    n_threads: usize,
    reverse: bool,
) -> (Vec<u8>, VecConcat, Vec<u8>) {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    let (finimizer_lengths, finimizer_matrix, global_length_counts) = thread_pool.install(|| {
        log::info!("Reading input sequences into memory");
        // Read the data into memory again
        let dbs: Vec<jseqio::seq_db::SeqDB> =
            filenames.iter().map(|p| load_db(p, reverse)).collect();
        let dbs = &dbs;
        let sbwt_len = sbwt.n_sets();
        let n_colors = dbs.len();

        let index = StreamingIndex::new(sbwt, lcs);
        let finimizer_lengths: Vec<AtomicU8> = (0..sbwt_len).map(|_| AtomicU8::new(0)).collect();

        let per_color_results_sparse: Vec<(Vec<usize>, Vec<u32>)> = dbs
            .par_iter()
            .enumerate()
            .map(|(db_idx, db)| {
                log::info!("Processing file {}/{}", db_idx + 1, dbs.len());
                let mut this_color_colex_ranks = Vec::<usize>::new();
                let mut unique_lengths: Vec<u32> = vec![0; k + 1];
                let mut prev_colex: Option<usize> = None;
                for rec in db.iter() {
                    let seq = rec.seq;
                    if seq.len() >= k {
                        let sus = index.shortest_freq_bound_suffixes(seq, 1); // Shortest unique suffixes
                        for kmer_end in k..seq.len() + 1 {
                            // Exclusive end
                            if !has_dna_alphabet(&seq[(kmer_end - k)..kmer_end]) {
                                continue;
                            }

                            let (len, colex, _) = pick_finimizer(&sus[(kmer_end - k)..kmer_end]);

                            finimizer_lengths[colex]
                                .store(len as u8, std::sync::atomic::Ordering::Relaxed);
                            this_color_colex_ranks.push(colex);

                            if Some(colex) != prev_colex {
                                unique_lengths[len] += 1;
                            }
                            prev_colex = Some(colex);
                        }
                    }
                }

                // Sort and deduplicate. Later code depends on the assumption that this is done.
                this_color_colex_ranks.sort();
                this_color_colex_ranks.dedup();
                (this_color_colex_ranks, unique_lengths)
            })
            .collect();

        let mut per_color_colex_ranks: Vec<Vec<usize>> = per_color_results_sparse
            .iter()
            .map(|(bv, _)| bv.clone())
            .collect();

        let mut global_length_counts: Vec<u32> = vec![0; k + 1];
        for (_, counts) in per_color_results_sparse {
            for (len, &count) in counts.iter().enumerate() {
                global_length_counts[len] += count;
            }
        }

        // Convert atomic integers to regular integers. This seems to compile into a basically no-op.
        let finimizer_lengths: Vec<u8> = finimizer_lengths
            .into_iter()
            .map(|x| x.into_inner())
            .collect();

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

        let color_set_sizes: Vec<AtomicU64> =
            (0..n_finimizers).map(|_| AtomicU64::new(0)).collect();

        log::info!("Computing color set sizes");
        let bar = indicatif::ProgressBar::new((n_colors) as u64);
        (0..n_colors).into_par_iter().for_each(|color| {
            let mut finimizer_idx = 0_usize; // Keep track of the current finimizer index
            let mut sbwt_pos = 0_usize;
            for &colex in per_color_colex_ranks[color].iter() {
                // Assumes the Vec is sorted
                while sbwt_pos < colex {
                    if finimizer_lengths[sbwt_pos] > 0 {
                        // This finimizer exists
                        finimizer_idx += 1;
                    }
                    sbwt_pos += 1;
                }

                if finimizer_lengths[colex] > 0 {
                    // This colex position is a finimizer
                    color_set_sizes[finimizer_idx]
                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
            }
            bar.inc(1);
        });
        bar.finish();

        // Convert from atomic to regular integers
        let color_set_sizes: Vec<u64> = color_set_sizes
            .into_iter()
            .map(|x| x.into_inner())
            .collect();

        // Figure out the starting points of color sets in the concatenation
        let mut color_set_starts: Vec<usize> = vec![0; color_set_sizes.len() + 1]; // Sentinel value at the end
        for i in 1..color_set_starts.len() {
            color_set_starts[i] = color_set_starts[i - 1] + color_set_sizes[i - 1] as usize;
        }
        drop(color_set_sizes); // Free memory

        // Allocate the space for the color set concatenation
        let mut color_set_concat = vec![0; *color_set_starts.last().unwrap()];

        log::info!("Constructing color sets");
        for (color, colex_ranks_for_color) in per_color_colex_ranks.into_iter().enumerate() {
            let mut finimizer_idx = 0_usize; // Keep track of the current finimizer index
            let mut sbwt_pos = 0_usize;
            for colex in colex_ranks_for_color {
                while sbwt_pos < colex {
                    if finimizer_lengths[sbwt_pos] > 0 {
                        // This finimizer exists
                        finimizer_idx += 1;
                    }
                    sbwt_pos += 1;
                }

                if finimizer_lengths[colex] > 0 {
                    // This colex position is a finimizer
                    color_set_concat[color_set_starts[finimizer_idx]] = color;
                    color_set_starts[finimizer_idx] += 1;
                }
            }
        }

        // color_set_starts are now ends
        color_set_starts.pop(); // Remove the sentinel at the end

        (
            finimizer_lengths,
            VecConcat {
                concat: color_set_concat,
                ends: color_set_starts,
            },
            global_length_counts,
        )
    });
    let mut lengths_by_freq: Vec<u8> = {
        let mut vec: Vec<(u8, u32)> = global_length_counts
            .iter()
            .enumerate()
            .filter(|&(len, &count)| len > 0 && count > 0)
            .map(|(len, &count)| (len as u8, count))
            .collect();
        vec.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        vec.into_iter().map(|(len, _)| len).collect()
    };

    (finimizer_lengths, finimizer_matrix, lengths_by_freq)
}

// Thread-seconds spent in the parts of the finimizer search, summed over all threads.
#[derive(Default)]
struct SearchTimes {
    load_ns: AtomicU64,   // Reading and parsing the input files
    suffix_ns: AtomicU64, // Shortest unique suffix search in the SBWT
    pick_ns: AtomicU64,   // Picking the finimizer of each k-mer
}

fn add_elapsed(counter: &AtomicU64, start: std::time::Instant) {
    counter.fetch_add(
        start.elapsed().as_nanos() as u64,
        std::sync::atomic::Ordering::Relaxed,
    );
}

// Calls f(len, colex) for the finimizer of every k-mer of the sequences in db, in sequence order.
//
// The finimizer of a k-mer is its shortest unique substring, with ties broken by colex rank
// (see pick_finimizer, which finds it by scanning the k-mer). Here the k-mers of a sequence are
// handled together with a sliding window minimum: a shortest unique suffix of length len that
// ends at position p lies inside the k-mers that end (exclusive) at p+1 ..= p+1-len+k, and
// `candidates` holds those that can still be the minimum of some k-mer, smallest first.
fn for_each_finimizer<F: FnMut(usize, usize)>(
    index: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>,
    db: &jseqio::seq_db::SeqDB,
    k: usize,
    times: &SearchTimes,
    mut f: F,
) {
    // (len, colex, last exclusive k-mer end that contains the candidate)
    let mut candidates = std::collections::VecDeque::<(usize, usize, usize)>::new();
    for rec in db.iter() {
        let seq = rec.seq;
        let start = std::time::Instant::now();
        let sus = index.shortest_freq_bound_suffixes(seq, 1); // Shortest unique suffixes
        add_elapsed(&times.suffix_ns, start);

        let start = std::time::Instant::now();
        candidates.clear();
        let mut dna_run = 0_usize; // Number of consecutive DNA characters ending here
        for (p, x) in sus.iter().enumerate() {
            dna_run = if is_dna_char(seq[p]) { dna_run + 1 } else { 0 };

            if let Some((len, range)) = x {
                // A new candidate replaces the earlier ones that are not smaller,
                // because it stays inside the k-mers for at least as long as they do.
                while let Some(&(back_len, back_colex, _)) = candidates.back() {
                    if (back_len, back_colex) >= (*len, range.start) {
                        candidates.pop_back();
                    } else {
                        break;
                    }
                }
                candidates.push_back((*len, range.start, p + 1 + k - len));
            }

            let kmer_end = p + 1; // Exclusive end
            if kmer_end < k || dna_run < k {
                continue;
            }

            // The full k-mer should have an existing unique match
            match x {
                Some((_, range)) if range.len() == 1 => (),
                _ => panic!("The k-mer ending at {} has no unique match", kmer_end),
            }

            while let Some(&(_, _, last_end)) = candidates.front() {
                if last_end < kmer_end {
                    candidates.pop_front();
                } else {
                    break;
                }
            }
            let &(len, colex, _) = candidates.front().expect("Finimizer not found for kmer");
            f(len, colex);
        }
        add_elapsed(&times.pick_ns, start);
    }
}

// Rank support over the SBWT positions that are finimizers: rank(i) is the number of
// finimizer positions before position i. Takes two bits per SBWT position.
struct FinimizerRank {
    words: Vec<u64>,       // Bit i is set iff SBWT position i is a finimizer
    block_ranks: Vec<u64>, // Number of set bits before each word
}

impl FinimizerRank {
    fn new(finimizer_lengths: &[u8]) -> Self {
        let mut words = vec![0_u64; finimizer_lengths.len().div_ceil(64)];
        for (i, &len) in finimizer_lengths.iter().enumerate() {
            if len > 0 {
                words[i / 64] |= 1_u64 << (i % 64);
            }
        }
        let mut block_ranks = Vec::with_capacity(words.len());
        let mut total = 0_u64;
        for w in words.iter() {
            block_ranks.push(total);
            total += w.count_ones() as u64;
        }
        Self { words, block_ranks }
    }

    fn rank(&self, i: usize) -> usize {
        let below = self.words[i / 64] & ((1_u64 << (i % 64)) - 1);
        self.block_ranks[i / 64] as usize + below.count_ones() as usize
    }
}

// Appends x to out as a variable-length integer (7 bits per byte, low bits first).
fn push_varint(out: &mut Vec<u8>, mut x: usize) {
    while x >= 0x80 {
        out.push((x as u8 & 0x7f) | 0x80);
        x >>= 7;
    }
    out.push(x as u8);
}

// Encodes a strictly increasing list of integers as the varints of its gaps.
fn encode_increasing(values: &[usize]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut prev = 0_usize;
    for &x in values {
        push_varint(&mut out, x - prev);
        prev = x;
    }
    out.shrink_to_fit();
    out
}

// Calls f on each integer of a list encoded with encode_increasing, in increasing order.
fn decode_increasing<F: FnMut(usize)>(bytes: &[u8], mut f: F) {
    let mut value = 0_usize;
    let mut gap = 0_usize;
    let mut shift = 0_u32;
    for &b in bytes {
        gap |= ((b & 0x7f) as usize) << shift;
        if b & 0x80 == 0 {
            value += gap;
            f(value);
            gap = 0;
            shift = 0;
        } else {
            shift += 7;
        }
    }
}

// Returns (v, bitmap), where:
// v is a vec of length |SBWT| such that v[i] is zero if k-mer i does not have a finimizer as a suffix, otherwise
// v[i] is the length of the finimizer that is a suffix of k-mer i.
// bitmap is the concatenation of the bitmap representation of all finimizers in colex order.
//
// Each file is loaded only while it is being processed. What is kept for each color is the
// sorted list of the SBWT positions of its finimizers, gap-encoded. The color matrix is filled
// from those lists once the full set of finimizers is known.
pub fn identify_finimizers<P: AsRef<Path> + Send + Sync>(
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: &sbwt::LcsArray,
    filenames: &[P],
    k: usize,
    n_threads: usize,
    reverse: bool,
) -> (Vec<u8>, BitVec<u64, Lsb0>, Vec<u8>) {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    let (finimizer_lengths, finimizer_matrix, global_length_counts) = thread_pool.install(|| {
        let sbwt_len = sbwt.n_sets();
        let n_colors = filenames.len();

        let index = StreamingIndex::new(sbwt, lcs);
        let finimizer_lengths: Vec<AtomicU8> = (0..sbwt_len).map(|_| AtomicU8::new(0)).collect();

        log::info!("Identifying finimizers");
        let times = SearchTimes::default();
        let per_color_results: Vec<(Vec<u8>, Vec<u32>)> = filenames
            .par_iter()
            .enumerate()
            .map(|(db_idx, path)| {
                log::info!("Processing file {}/{}", db_idx + 1, n_colors);
                let start = std::time::Instant::now();
                let db = load_db(path, reverse);
                add_elapsed(&times.load_ns, start);
                let mut this_color_colex_ranks = Vec::<usize>::new();
                let mut unique_lengths: Vec<u32> = vec![0; k + 1];
                let mut prev_colex: Option<usize> = None;
                for_each_finimizer(&index, &db, k, &times, |len, colex| {
                    // Same effect as always storing, without writing to shared memory needlessly
                    if finimizer_lengths[colex].load(std::sync::atomic::Ordering::Relaxed) != len as u8 {
                        finimizer_lengths[colex]
                            .store(len as u8, std::sync::atomic::Ordering::Relaxed);
                    }
                    if Some(colex) != prev_colex {
                        this_color_colex_ranks.push(colex);
                        unique_lengths[len] += 1;
                    }
                    prev_colex = Some(colex);
                });
                drop(db);

                this_color_colex_ranks.sort_unstable();
                this_color_colex_ranks.dedup();
                (encode_increasing(&this_color_colex_ranks), unique_lengths)
            })
            .collect();

        let mut global_length_counts: Vec<u32> = vec![0; k + 1];
        let mut per_color_lists: Vec<Vec<u8>> = Vec::with_capacity(n_colors);
        for (list, counts) in per_color_results {
            for (len, &count) in counts.iter().enumerate() {
                global_length_counts[len] += count;
            }
            per_color_lists.push(list);
        }

        // Convert atomic integers to regular integers. This seems to compile into a basically no-op.
        let finimizer_lengths: Vec<u8> = finimizer_lengths
            .into_iter()
            .map(|x| x.into_inner())
            .collect();

        // Figure out how many finimizers we have and their average length
        let mut n_finimizers = 0_usize;
        let mut total_length = 0_usize;

        for &x in finimizer_lengths.iter() {
            if x > 0 {
                n_finimizers += 1;
                total_length += x as usize;
            }
        }

        let seconds = |x: &AtomicU64| x.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9;
        log::info!(
            "Thread-seconds: reading input {:.0}, suffix search {:.0}, picking finimizers {:.0}",
            seconds(&times.load_ns),
            seconds(&times.suffix_ns),
            seconds(&times.pick_ns)
        );
        log::info!("Number of SBWT positions = {}", sbwt_len);
        log::info!("Number of finimizers = {}", n_finimizers);
        log::info!("Avg length = {}", total_length as f64 / n_finimizers as f64);
        log::info!(
            "Per-color finimizer lists take {:.3} GB",
            per_color_lists.iter().map(|l| l.len()).sum::<usize>() as f64 / 1e9
        );

        let finimizer_rank = FinimizerRank::new(&finimizer_lengths);

        let total_n_bits = n_colors * n_finimizers;
        let total_n_words = total_n_bits.div_ceil(64);
        log::info!(
            "Filling the color matrix ({} colors x {} finimizers, {:.3} GB)",
            n_colors,
            n_finimizers,
            (total_n_words * 8) as f64 / 1e9
        );
        let finimizer_matrix: Vec<AtomicU64> =
            (0..total_n_words).map(|_| AtomicU64::new(0)).collect();

        let bar = indicatif::ProgressBar::new((n_colors) as u64);
        per_color_lists
            .into_par_iter()
            .enumerate()
            .for_each(|(color, list)| {
                decode_increasing(&list, |colex| {
                    let bit_idx = finimizer_rank.rank(colex) * n_colors + color; // We want to set this bit
                    finimizer_matrix[bit_idx / 64]
                        .fetch_or(1_u64 << (bit_idx % 64), std::sync::atomic::Ordering::Relaxed);
                });
                bar.inc(1);
            });
        bar.finish();
        log::info!("Color matrix finished");

        // Convert atomic integers to regular integers. This seems to compile into a basically no-op.
        let finimizer_matrix: Vec<u64> = finimizer_matrix
            .into_iter()
            .map(|x| x.into_inner())
            .collect();

        (
            finimizer_lengths,
            BitVec::<u64, Lsb0>::from_vec(finimizer_matrix),
            global_length_counts,
        )
    });

    let mut lengths_by_freq: Vec<u8> = {
        let mut vec: Vec<(u8, u32)> = global_length_counts
            .iter()
            .enumerate()
            .filter(|&(len, &count)| len > 0 && count > 0)
            .map(|(len, &count)| (len as u8, count))
            .collect();
        vec.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        vec.into_iter().map(|(len, _)| len).collect()
    };

    (finimizer_lengths, finimizer_matrix, lengths_by_freq)
}

pub fn verify_finimizer_shortest_uniqueness(
    cf: &ColoredFinimizers,
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    n_threads: usize,
) {
    // Look up all finimizer strings, verify that their SBWT frequency is 1, and the SBWT frequency of the
    // longest proper suffix is more than .
    let n_finimizers = cf.finimizer_lengths.len();

    // Finimizer length cumulative sum up to but not including the array index
    let mut f_starts: Vec<usize> = vec![0; n_finimizers + 1];
    for i in 0..n_finimizers {
        f_starts[i + 1] = f_starts[i] + cf.finimizer_lengths[i] as usize;
    }

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        (0..n_finimizers).par_bridge().for_each(|i| {
            let s = &cf.finimizers_concat[f_starts[i]..f_starts[i + 1]];
            assert!(sbwt.search(s).unwrap().len() == 1); // Is unique
            assert!(sbwt.search(&s[1..]).unwrap().len() > 1); // Is not unique
        });
    });
}

pub fn verify_finimizer_colors(
    cf: &ColoredFinimizers,
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: &LcsArray,
    n_threads: usize,
    paths: &[PathBuf],
    reverse: bool,
) {
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

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();
    thread_pool.install(|| {
        paths.par_iter().enumerate().for_each(|(color, path)| {
            log::info!("Processing color {}", path.display());
            let db = load_db(path, reverse); // Also has reverse complements

            // Mark all finimizers (in finimizer order, not sbwt colex order) in this color
            let mut finimizer_marks: BitVec = bitvec![0; n_finimizers];
            for rec in db.iter() {
                let sfs = index.shortest_freq_bound_suffixes(rec.seq, 1);
                for kmer_end in k..rec.seq.len() + 1 {
                    // Exclusive end
                    if !has_dna_alphabet(&rec.seq[(kmer_end - k)..kmer_end]) {
                        continue;
                    }
                    let (len, _, mut end) = pick_finimizer(&sfs[(kmer_end - k)..kmer_end]);
                    end += 1; // Exclusive end
                    end += kmer_end - k; // Shift to where we are at rec.seq
                    let finimizer = &rec.seq[(end - len)..end];

                    match finimizers.binary_search_by(|test| colex_compare(test, finimizer)) {
                        Ok(f_idx) => {
                            finimizer_marks.set(f_idx, true);
                        }
                        Err(_) => panic!(
                            "Finimizer {} not found with binary search",
                            String::from_utf8_lossy(finimizer)
                        ),
                    }
                }
            }

            // Check that the marks in ColoredFinimizers match our marks
            for f_idx in 0..n_finimizers {
                assert!(cf.set_has_color(f_idx, color) == finimizer_marks[f_idx]);
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
pub fn get_finimizer_strings(
    sbwt: &mut sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    finimizer_lengths: &[u8],
    n_threads: usize,
) -> Vec<u8> {
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
    let mut output_segments: Vec<&mut [u8]> = Vec::new();
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
        let (a, b) = tail.split_at_mut(total_finimizer_length_found);
        output_segments.push(a);
        tail = b;
    }
    assert!(tail.is_empty()); // Segmentation is complete

    log::info!("Extracting finimizer strings");
    // Fill in segments in parallel
    output_segments
        .iter_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(t, out_segment)| {
            let mut output_idx = 0_usize;
            for colex in input_sbwt_ranges[t].clone() {
                if finimizer_lengths[colex] > 0 {
                    // Write the finimizer to the concatenation
                    let kmer = sbwt.access_kmer(colex);
                    let finimizer = &kmer[kmer.len() - finimizer_lengths[colex] as usize..];
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
struct VecConcat {
    concat: Vec<usize>,
    ends: Vec<usize>, // Exclusive ends
}

impl VecConcat {
    fn serialize<W: Write>(&self, out: &mut W) {
        let n_elements = self.concat.len() as u64;
        let n_sets = self.ends.len() as u64;

        out.write_all(&n_elements.to_le_bytes()).unwrap();
        out.write_all(&n_sets.to_le_bytes()).unwrap();

        let concat_bytes: &[u8] = bytemuck::cast_slice(self.concat.as_slice());
        out.write_all(concat_bytes).unwrap();

        let ends_bytes: &[u8] = bytemuck::cast_slice(self.ends.as_slice());
        out.write_all(ends_bytes).unwrap();
    }

    fn load<R: std::io::Read>(input: &mut R) -> Self {
        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf).unwrap();
        let n_elements = u64::from_le_bytes(buf);

        input.read_exact(&mut buf).unwrap();
        let n_sets = u64::from_le_bytes(buf);

        let mut concat_bytes = vec![0_u8; (n_elements as usize) * std::mem::size_of::<usize>()];
        input.read_exact(concat_bytes.as_mut_slice()).unwrap();
        let concat: Vec<usize> = bytemuck::cast_slice(&concat_bytes).to_vec();
        // Todo: cast_slice makes a copy

        let mut ends_bytes = vec![0_u8; (n_sets as usize) * std::mem::size_of::<usize>()];
        input.read_exact(ends_bytes.as_mut_slice()).unwrap();
        let ends: Vec<usize> = bytemuck::cast_slice(&ends_bytes).to_vec();
        // Todo: cast_slice makes a copy

        VecConcat { concat, ends }
    }
}

#[derive(Debug)]
enum ColorSetStorage {
    Dense(BitVec<u64, Lsb0>), // Concatenation of bitvectors
    Sparse(VecConcat),
}

#[derive(Debug)]
pub struct ColoredFinimizers {
    finimizer_lengths: Vec<u8>, // Lengths of finimizers, in colex order
    finimizers_concat: Vec<u8>, // concatenation of finimizers in colex order in ASCII
    //color_sets_bitmap: BitVec<u64, Lsb0>, // concatenation of finimizer color set bit maps, in colex order of finimizers
    color_sets: ColorSetStorage,
    n_colors: usize,
    pub lengths_by_freq: Vec<u8>,
}

impl ColoredFinimizers {
    pub fn serialize<W: std::io::Write>(&self, mut out: &mut W) {
        let n_finimizers = self.finimizer_lengths.len();
        let finimizer_total_length = self.finimizers_concat.len();

        // Write constants
        out.write_all(&(n_finimizers as u64).to_le_bytes()).unwrap();
        out.write_all(&(finimizer_total_length as u64).to_le_bytes())
            .unwrap();
        out.write_all(&(self.n_colors as u64).to_le_bytes())
            .unwrap();

        // Write data
        out.write_all(self.finimizer_lengths.as_slice()).unwrap();
        out.write_all(self.finimizers_concat.as_slice()).unwrap();

        match &self.color_sets {
            ColorSetStorage::Dense(bitmap) => {
                out.write_all(&[0u8]).unwrap(); // Byte indicating that this is a dense structure
                let data_bytes: &[u8] = bytemuck::cast_slice(bitmap.as_raw_slice());
                assert_eq!(
                    data_bytes.len(),
                    ((n_finimizers * self.n_colors + 63) / 64) * 8
                );
                out.write_all(data_bytes).unwrap();
            }
            ColorSetStorage::Sparse(vec_concat) => {
                out.write_all(&[1u8]).unwrap(); // Byte indicating that this is a sparse structure
                vec_concat.serialize(out);
            }
        }
        out.write_all(&(self.lengths_by_freq.len() as u64).to_le_bytes())
            .unwrap();
        out.write_all(self.lengths_by_freq.as_slice()).unwrap();
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

        let mut is_sparse = [0_u8; 1];
        input.read_exact(&mut is_sparse).unwrap();
        if is_sparse[0] == 1 {
            let vec_concat = VecConcat::load(input);

            // Read lengths_by_freq (as Vec<u8>)
            input.read_exact(&mut buf).unwrap();
            let n_lengths_by_freq = u64::from_le_bytes(buf);
            let mut lengths_by_freq = vec![0_u8; n_lengths_by_freq as usize];
            input.read_exact(&mut lengths_by_freq).unwrap();

            ColoredFinimizers {
                finimizer_lengths,
                finimizers_concat,
                color_sets: ColorSetStorage::Sparse(vec_concat),
                n_colors: n_colors as usize,
                lengths_by_freq,
            }
        } else {
            let mut bitmap_bytes = vec![0_u8; ((n_finimizers * n_colors + 63) as usize / 64) * 8];
            input.read_exact(bitmap_bytes.as_mut_slice()).unwrap();
            let bitmap_words: &[u64] = bytemuck::cast_slice(&bitmap_bytes);
            let color_sets_bitmap = BitVec::<u64, Lsb0>::from_slice(bitmap_words); // Todo: this makes a copy

            // Read lengths_by_freq (as Vec<u8>)
            input.read_exact(&mut buf).unwrap();
            let n_lengths_by_freq = u64::from_le_bytes(buf);
            let mut lengths_by_freq = vec![0_u8; n_lengths_by_freq as usize];
            input.read_exact(&mut lengths_by_freq).unwrap();

            ColoredFinimizers {
                finimizer_lengths,
                finimizers_concat,
                color_sets: ColorSetStorage::Dense(color_sets_bitmap),
                n_colors: n_colors as usize,
                lengths_by_freq,
            }
        }
    }

    pub fn verify_colex_order(&self) {
        let mut cumul_len = 0_usize;
        for i in 0..self.finimizer_lengths.len() {
            if i > 0 {
                let a = &self.finimizers_concat
                    [(cumul_len - self.finimizer_lengths[i - 1] as usize)..cumul_len];
                let b = &self.finimizers_concat
                    [cumul_len..cumul_len + self.finimizer_lengths[i] as usize];
                match a.iter().rev().cmp(b.iter().rev()) {
                    std::cmp::Ordering::Equal => panic!("Equal finimizers"),
                    std::cmp::Ordering::Greater => panic!("Non-colex-ordered finimizers"),
                    std::cmp::Ordering::Less => {} // Ok
                }
            }
            cumul_len += self.finimizer_lengths[i] as usize;
        }
    }

    // Returns the set with the given index as an indicator bit vector
    pub fn set_has_color(&self, set_idx: usize, color: usize) -> bool {
        match &self.color_sets {
            ColorSetStorage::Dense(bitmap) => bitmap[set_idx * self.n_colors + color],
            ColorSetStorage::Sparse(ds) => {
                let start = if set_idx == 0 {
                    0
                } else {
                    ds.ends[set_idx - 1]
                };
                let end = ds.ends[set_idx];
                ds.concat[start..end].iter().any(|&c| c == color)
            }
        }
    }
}

// Loads one sequence file into memory, with the reverse complements of its sequences
// appended after the forward sequences. If reverse is true, every sequence is reversed.
fn load_db<P: AsRef<Path>>(path: P, reverse: bool) -> jseqio::seq_db::SeqDB {
    let reader = jseqio::reader::DynamicFastXReader::from_file(&path).unwrap();
    let (mut fw, rc) = reader.into_db_with_revcomp().unwrap();

    if fw.sequence_count() == 0 {
        panic!("No sequences found in file {}", path.as_ref().display());
    }

    // Append reverse complement records to the forward database
    for rec in rc.iter() {
        fw.push_record(rec);
    }

    if reverse {
        for i in 0..fw.sequence_count() {
            let rec = fw.get_mut(i);
            rec.seq.reverse();
            rec.qual.map(|q| q.reverse());
            // Header is not reversed
        }
    }
    fw
}

// State shared between an InputStream and its loader threads.
struct LoaderState {
    next_to_load: usize,                           // Index of the next file a loader should take
    next_to_consume: usize,                        // Index of the next file the stream will return
    loaded: HashMap<usize, jseqio::seq_db::SeqDB>, // Files that are loaded and not yet consumed
    failure: Option<String>,                       // Set if a loader failed
    stop: bool,                                    // Set when the stream is dropped
}

struct Loader {
    paths: Vec<PathBuf>,
    reverse: bool,   // Reverse sequences (not complement)
    window: usize,   // Loaders stay less than this many files ahead of the consumer
    state: Mutex<LoaderState>,
    changed: std::sync::Condvar,
    load_ns: AtomicU64, // Thread-time spent reading and parsing files
}

impl Loader {
    // Body of a loader thread: load files in index order, staying inside the window.
    fn run(&self) {
        loop {
            let idx = {
                let mut state = self.state.lock().unwrap();
                loop {
                    if state.stop || state.failure.is_some() || state.next_to_load == self.paths.len() {
                        return;
                    }
                    if state.next_to_load < state.next_to_consume + self.window {
                        break;
                    }
                    state = self.changed.wait(state).unwrap();
                }
                state.next_to_load += 1;
                state.next_to_load - 1
            };

            let start = std::time::Instant::now();
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                load_db(&self.paths[idx], self.reverse)
            }));
            add_elapsed(&self.load_ns, start);

            let mut state = self.state.lock().unwrap();
            match result {
                Ok(db) => {
                    state.loaded.insert(idx, db);
                }
                Err(_) => {
                    state.failure = Some(format!("Could not load {}", self.paths[idx].display()));
                }
            }
            self.changed.notify_all();
        }
    }
}

// Streams the sequences of the input files in order. The files are read and parsed by
// background threads, a bounded number of files ahead of the one being streamed.
struct InputStream {
    loader: Arc<Loader>,
    cur_db: Option<jseqio::seq_db::SeqDB>, // The file currently being iterated over
    seq_idx_in_cur_db: usize,              // Index of the next sequence to return from cur_db
}

impl InputStream {
    fn new<P: AsRef<Path>>(filenames: &[P], reverse: bool, n_threads: usize) -> InputStream {
        let n_threads = max(n_threads, 1);
        let loader = Arc::new(Loader {
            paths: filenames.iter().map(|p| p.as_ref().to_path_buf()).collect(),
            reverse,
            window: 2 * n_threads,
            state: Mutex::new(LoaderState {
                next_to_load: 0,
                next_to_consume: 0,
                loaded: HashMap::new(),
                failure: None,
                stop: false,
            }),
            changed: std::sync::Condvar::new(),
            load_ns: AtomicU64::new(0),
        });
        for _ in 0..n_threads {
            let loader = Arc::clone(&loader);
            std::thread::spawn(move || loader.run());
        }
        Self {
            loader,
            cur_db: None,
            seq_idx_in_cur_db: 0,
        }
    }

    // Waits for the next file in order. Returns None when all files have been consumed.
    fn next_db(&mut self) -> Option<jseqio::seq_db::SeqDB> {
        let mut state = self.loader.state.lock().unwrap();
        let idx = state.next_to_consume;
        if idx == self.loader.paths.len() {
            return None;
        }
        loop {
            if let Some(message) = &state.failure {
                panic!("{}", message);
            }
            if let Some(db) = state.loaded.remove(&idx) {
                state.next_to_consume += 1;
                self.loader.changed.notify_all();
                return Some(db);
            }
            state = self.loader.changed.wait(state).unwrap();
        }
    }
}

impl Drop for InputStream {
    fn drop(&mut self) {
        self.loader.state.lock().unwrap().stop = true;
        self.loader.changed.notify_all();
    }
}

impl SeqStream for InputStream {
    fn stream_next(&mut self) -> Option<&[u8]> {
        // Move to the next file if the current one is exhausted
        let exhausted = match &self.cur_db {
            Some(db) => self.seq_idx_in_cur_db == db.sequence_count(),
            None => true,
        };
        if exhausted {
            let had_db = self.cur_db.take().is_some();
            match self.next_db() {
                Some(db) => self.cur_db = Some(db),
                None => {
                    if had_db {
                        log::info!(
                            "Reading the input files took {:.0} thread-seconds",
                            self.loader.load_ns.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9
                        );
                    }
                    return None; // Done
                }
            }
            self.seq_idx_in_cur_db = 0;
        }

        let seq = self.cur_db.as_ref().unwrap().get(self.seq_idx_in_cur_db).seq;
        self.seq_idx_in_cur_db += 1;
        Some(seq)
    }
}

fn mark_bits(bv: &mut BitVec, color: usize, num_colors: usize, to_mark: Vec<usize>) {
    for i in to_mark {
        bv.set(i * num_colors + color, true);
    }
}

fn mark_all_kmers_of_seq(
    bv: Arc<Mutex<BitVec>>,
    num_colors: usize,
    color: usize,
    seq: &[u8],
    k: usize,
    mark_buffer_size: usize,
    index: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>,
) {
    // Search all k-mers
    let mut marking_buffer: Vec<usize> = Vec::new(); // These bits should be marked
    for (len, colex) in index.matching_statistics(seq) {
        if len == k {
            // Full kmer -> set the bit in the color set of the k-mer
            assert!(colex.len() == 1);
            marking_buffer.push(colex.start);
            if marking_buffer.len() == mark_buffer_size {
                mark_bits(
                    bv.lock().unwrap().deref_mut(),
                    color,
                    num_colors,
                    marking_buffer,
                );
                marking_buffer = Vec::new();
            }
        }
    }

    if !marking_buffer.is_empty() {
        // Mark the rest
        mark_bits(
            bv.lock().unwrap().deref_mut(),
            color,
            num_colors,
            marking_buffer,
        );
    }
}

impl ColoredFinimizers {
    #[allow(clippy::type_complexity)]
    pub fn new<P: AsRef<Path> + Send + Sync>(
        filenames: &[P],
        k: usize,
        n_threads: usize,
        mem_gb: usize,
        temp_dir: &Path,
        reverse: bool,
        sparse: bool,
    ) -> Self {
        log::info!("Building for {} sequence files (colors)", filenames.len());
        let input_stream = InputStream::new(filenames, reverse, n_threads);
        let num_colors = filenames.len();
        log::info!("Building SBWT");
        let sbwt_builder = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false) // Already added in the input stream
            .k(k)
            .build_lcs(true)
            .n_threads(n_threads)
            .precalc_length(8)
            .algorithm(
                BitPackedKmerSorting::new()
                    .dedup_batches(true)
                    .temp_dir(temp_dir)
                    .mem_gb(mem_gb),
            );
        let (mut sbwt, lcs) = sbwt_builder.run(input_stream);
        let lcs = lcs.unwrap(); // Ok since used build_lcs(true) above

        if sparse {
            let (finimizer_lengths, interleaved_bitmap, sorted_lengths) =
                identify_finimizers_sparse(&sbwt, &lcs, filenames, k, n_threads, reverse);
            let concat = get_finimizer_strings(&mut sbwt, &finimizer_lengths, n_threads);
            let nonzero_lengths: Vec<u8> =
                finimizer_lengths.into_iter().filter(|x| *x > 0).collect();
            ColoredFinimizers {
                finimizer_lengths: nonzero_lengths,
                finimizers_concat: concat,
                color_sets: ColorSetStorage::Sparse(interleaved_bitmap),
                n_colors: num_colors,
                lengths_by_freq: sorted_lengths,
            }
        } else {
            let (finimizer_lengths, interleaved_bitmap, sorted_lengths) =
                identify_finimizers(&sbwt, &lcs, filenames, k, n_threads, reverse);
            let concat = get_finimizer_strings(&mut sbwt, &finimizer_lengths, n_threads);
            let nonzero_lengths: Vec<u8> =
                finimizer_lengths.into_iter().filter(|x| *x > 0).collect();
            ColoredFinimizers {
                finimizer_lengths: nonzero_lengths,
                finimizers_concat: concat,
                color_sets: ColorSetStorage::Dense(interleaved_bitmap),
                n_colors: num_colors,
                lengths_by_freq: sorted_lengths,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use bitvec::prelude::*;
    #[test]
    fn test_bit_order() {
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
