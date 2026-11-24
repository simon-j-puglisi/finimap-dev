# Finimap

**Finimap** is a tool for building and querying compressed colored finimizer indexes for fast sequence pseudoalignment.

Original paper: [**Finimap: fast and accurate single-species bacterial pseudoalignment with finimizers**](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2026.27) by J. N. Alanko, E. Biagi, S. J. Puglisi.

Finimizers paper: [**Finimizers: Variable-length bounded-frequency minimizers for k-mer sets**](https://ieeexplore.ieee.org/abstract/document/10908718) by J. N. Alanko, E. Biagi, S. J. Puglisi.

## Overview

A finimizer index enables fast querying of sequence data (e.g., DNA reads) against reference genomes or genome collections. The tool builds a compact index by:

1. **Extracting finimizers** from reference sequences
2. **Compacting the index**
3. **Querying efficiently** by extracting finimizers from query sequences and finding matches in the index

This approach provides superior speed compared to traditional k-mer indexes while maintaining query accuracy.

### Shortest Unique Finimizers

Let $G$ be the de Bruijn graph of a set of $k$-mers $R$, $t \geq 1$ be an integer, $X$ be a $k$-mer, and $Y$ be a substring of $X$. We say $Y$ is a **shortest $t$-finimizer** of $X$ with respect to the $k$-mer set $R$ if $Y$ has at most $t$ occurrences in $G$ and there does not exist a shorter substring of $X$ with at most $t$ occurrences in $G$. If $t = 1$ then we say $Y$ is a _shortest-unique finimizer_ of $X$.


## Building

### Prerequisites

- Rust compiler (for building the `finimizer_matrix` submodule)
- A C/C++ compiler (for the main project)
- Make

### Compilation Steps

First clone the repository and access it with:
```bash
git clone --recursive -b robust-query https://github.com/simon-j-puglisi/finimap-dev.git
cd finimap-dev
```

This is the development version of [ElenaBiagi/Finimap](https://github.com/ElenaBiagi/Finimap), on the branch `robust-query`.

If you forgot the flag `--recursive`, pull the submodules with:
```bash
git submodule update --init --recursive
```

Build the [finimizer_matrix](https://github.com/simon-j-puglisi/finimizer_matrix/tree/low-memory-build) submodule (a low-memory version of [jnalanko/finimizer_matrix](https://github.com/jnalanko/finimizer_matrix); see its README):
```bash
cd finimizer_matrix
cargo build -r
cd ..
```

Build the [sdsl-lite](https://github.com/iosfwd/sdsl-lite.git) dependency:
```bash
cd sdsl-lite && mkdir -p build && cd build
cmake .. -DCMAKE_CXX_FLAGS="-std=c++2a -O2 -march=native"
make -j4
cd ../..
```

Compile the main project:
```bash
make finimap
```

This will produce two main executables:
- `finimizer_matrix` (in the `finimizer_matrix/target/release/` directory)
- `finimap` (in the current directory)
## Index Construction

Building a finimizer index is a two-stage process:

1. **Stage 1**: Build a colored finimizer matrix (CFM) from reference sequences
2. **Stage 2**: Compact the CFM into a queryable finimizer index (FMIN)

### Stage 1: Building the Colored Finimizer Matrix

The first stage extracts finimizers from your reference sequences using `finimizer_matrix build`.

#### Input Format

Create a list file (`.list`) containing paths to your reference sequence files (one per line):

```bash
# example_data/coli_file_list.list
../example_data/GCA_000005845.2_ASM584v2.fna
../example_data/GCA_000006665.1_ASM666v1.fna
../example_data/GCA_000007445.1_ASM744v1.fna
```

#### Command

```bash
cd finimizer_matrix

./target/release/finimizer_matrix build [OPTIONS] --input <INPUT_LIST> --output <OUTPUT_CFM> --temp-dir <TEMP_DIR> -k <K_VALUE>
```

#### Parameters

- `-i, --input <FILE>`: Path to list file containing reference sequence paths (one per line, FASTA or FASTQ). Each file is one color
- `-o, --output <FILE>`: Output colored finimizer matrix file
- `-k <K>`: K-mer length (required; at most 31). Recommended values: 19-31
- `-t, --n-threads <THREADS>`: Number of threads (default: 4)
- `-d, --temp-dir <DIR>`: Temporary directory for intermediate files (required)
- `-m, --mem-gb <GB>`: Memory budget for building the SBWT, in GB (default: 4). Peak memory also includes the SBWT itself and the color data, so it can be well above this figure
- `--reverse`: Include reverse complement of sequences (normal for Finimap)
- `--sparse`: Store color sets as lists of integers. The default (dense) build is the one used and tested with `build-fmin`

On 3,498 *E. coli* genomes (k = 31, `-t 32 -m 64`), stage 1 takes about 12 minutes and 67 GB of memory.

#### Example

```bash
cd finimizer_matrix

./target/release/finimizer_matrix build \
  -i ../example_data/coli_file_list.list \
  -o ../example_data/coli.cfm \
  -k 31 \
  -t 4 \
  -m 8 \
  -d ./temp \
  --reverse \
  2> /tmp/build.log
```

This builds a CFM file named `coli.cfm` from the references listed in `coli_file_list.list`.

### Stage 2: Compacting the Index with `build-fmin`

The second stage compacts the CFM into an efficient queryable index format.

#### Command

```bash
cd ..

./finimap build-fmin -i <CFM_FILE> -o <INDEX_PREFIX> [OPTIONS]
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-i, --index-file` | FILE | Required | Colored finimizer matrix (CFM) file from Stage 1 |
| `-o, --out-file` | FILE | Required | Output index filename prefix (without extension) |
| `-p, --p_len` | INT | 10 | Finimizer prefix length (bits used for bucketing). Higher values → more buckets → faster but larger index |
| `-x` | INT | 16 | Short/long finimizer threshold, in characters (1 to 19). Finimizers of at most this length are stored in one structure, longer ones in another. It sets the layout of the index, not the results |
| `-k` | INT | 31 | K-mer length (must match the value used in Stage 1) |
| `-m, --meta` | flag | off | Metagenome mode: for a `.cfm` built by stage 1 with `--sparse` (color sets stored as lists). Required for such a file and refused for a dense one |
| `--index-format` | INT | 4 | Index format to write. 4 stores the color sets as meta-color sets, several times smaller than format 3; 1 to 3 are the formats of earlier versions, kept for comparison. All give the same results, and `search-fmin` reads all four |
| `-h, --help` | - | - | Display help message |

#### Parameter Tuning

- **`-p` (prefix length)**: Controls index bucket size in long finimizers. Typical values: 8-12, (p > x)
  - Smaller values: smaller index, slower queries
  - Larger values: faster queries, larger index size
  
- **`-x` (short/long threshold)**: Affects short/long finimizer representation, never the results
  - The bucket table for short finimizers takes about 4^x / 228 bytes whatever the collection:
    19 MB at 16, 302 MB at 18. Long finimizers take 8 bytes each instead of 2.
  - On 3,498 *E. coli* genomes, `-x 16` gave a 1.94 GB index against 2.94 GB at `-x 18` in
    index format 1, with identical results and about the same query time; `-x 14` was slightly
    larger and slower than `-x 16`. The paper used `-x 18`.

#### Index formats

- **Format 4** (the default) is format 3 with the color sets stored as
  meta-color sets (`include/MetaColorSets.hh`). The colors are ordered by average-linkage
  clustering of how they co-occur in the color sets (on a sample of 65,536 sets; above 16,384
  colors the input order is kept), and split into blocks of at most 256 colors that follow the
  clustering tree. Each block stores its distinct partial sets once, in the smallest of five
  forms: a list of the colors it has or of those it lacks, the runs of consecutive colors it
  has or lacks, or a bitmap. Each block has a default: its most common partial
  set. Each color set is a row: a bitmask of the blocks where its partial set is not the
  default, then the ids of those partial sets in a canonical Huffman code of each block.
  Queries add each distinct partial set of a read once, with the number of the read's
  finimizers that use it. Format 4 also stores the long finimizers in 67 - 2p bits rounded
  up to whole bytes instead of 64 (48 at the default `-p 10`). On 3,498 *E. coli* genomes
  (`-x 16`) the index takes 468 MB, against 1.55 GB in format 3; on 9,500 genomes each of
  *Salmonella*, *S. aureus* and *S. pneumoniae*, 411, 220 and 223 MB, against 2.60, 1.83 and
  1.06 GB. Query results are the same; queries take from about the same time to a third
  longer (more where the color sets vary more), and the index loads several times faster.
- **Format 3** (the default before format 4) is format 2 with three arrays packed in as many bits as their
  values need: the color set id of each finimizer (23 bits instead of 64 for 3,498 *E. coli*
  genomes), the colors of the listed color sets (12 bits instead of 16 below 4,096 colors), and
  the lengths of the short finimizers (4 bits instead of 8 at `-x 16`).
- **Format 2** stores the short finimizers in a Pino structure (*Faster Small-Alphabet Subset
  Rank Queries*, Section 4.2): the byte of low bits of each finimizer, and the start of each
  bucket of 256 values coded in about 9 bits instead of a 32-bit pointer. Files of format 2
  and later start with the signature `FINIMAP!` and the format number.
- **Format 1** is the format of earlier versions, with a 32-bit pointer per bucket (4^x / 64
  bytes). It has no signature.
- `build-fmin --index-format 1`, `2` or `3` still writes those formats, byte for byte as before.
  All four give the same results. `dev/index_sizes.py` reports the size of each part of any
  of them.

#### Examples

**Single genome index (e.g., *E. coli*):**
```bash
./finimap build-fmin \
  -i example_data/coli.cfm \
  -o example_data/coli_index \
  -k 31 \
  -p 10 \
  -x 16
```

This produces an index `example_data/coli_index.fmin`.

**Compact index with smaller memory footprint:**
```bash
./finimap build-fmin \
  -i example_data/coli.cfm \
  -o example_data/compact_index \
  -k 31 \
  -p 8 \
  -x 15
```


## Querying

The `search-fmin` command pseudoaligns your query sequences against the built index: for each
read it reports the reference genomes (colors) that share enough of the read's k-mers, judged
through the k-mers' finimizers.

### Command

```bash
./finimap search-fmin -i <INDEX_PREFIX> -q <QUERY_FILE> -o <OUTPUT_FILE> [OPTIONS]
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-i, --index-file` | FILE | Required | Index filename prefix (from `build-fmin`) |
| `-q, --query-file` | FILE | Required | Query sequences in FASTA/FASTQ format (may be gzipped). Or `.txt` file with list of query files (one per line) |
| `-o, --out-file` | FILE | stdout | Output filename. Or `.txt` file with output file list (one per line, paired with query list) |
| `-t` | FLOAT | 0 | Threshold, a fraction from 0 to 1. A color is reported for a read if its score is at least `floor(t * found)`, and at least 1, where `found` is the number of the read's k-mers that have a finimizer in the index. With 0, every color with a score of at least 1 is reported |
| `-b` | INT | 1000 | Number of reads handled together in a batch. It affects speed and memory, not the results |
| `--threads` | INT | 1 | Number of threads. Each thread handles one batch at a time; the results and their order do not depend on it |
| `--relative` | INT | off | Relative rule: report only the colors whose score is at least the read's highest score minus this many k-mers. Combines with `-t` (a color must pass both), so with the default `-t 0` it is the relative rule alone. One sequencing error, or one difference between genomes, changes up to k scores, so `--relative 31` allows about one at k = 31 |
| `-h, --help` | - | - | Display help message |

### Query files

- The format (FASTA or FASTQ) is recognised from the content, and gzipped files are read directly.
- Reads may have different lengths. A read shorter than k has no k-mers and is reported with no colors.
- Upper and lower case are both accepted. A k-mer that contains any character other than A, C, G or T is skipped.
- FASTA sequences may span several lines. FASTQ records must be four lines each.

### Output Format

The output contains one line per query sequence, in the order of the query file:
```
QUERY_ID color_1:score_1 color_2:score_2 color_3:score_3 ...
```

Where:
- **QUERY_ID**: Number of the read in the query file, starting from 0
- **color_N**: Index of a reference genome in the index
- **score_N**: Number of the read's k-mers whose finimizer occurs in that genome. It is never below the number of the read's k-mers that occur in the genome

Colors are assigned in the order references were listed in the original input file.

```color:score``` pairs are sorted by decreasing score, and by increasing color within a score.
A read with no reported color gives a line with only its number. Progress messages go to
standard error, so the output holds result lines only. At the end, standard error also gets a
summary of the time spent loading the index and in each part of the search.

### Examples

**All colors that share at least one finimizer (no minimum match requirement):**
```bash
./finimap search-fmin \
  -i example_data/coli_index \
  -q example_data/queries.fna \
  -b 4 \
  -o results.txt \
  -t 0
```

Output example:
```
0 0:50 2:50 1:32 
1 2:50 0:44 1:44 
2 2:26 0:23 1:2 
3 2:50 0:21 1:21
```

**Threshold queries (80% match required):**
```bash
./finimap search-fmin \
  -i example_data/coli_index \
  -q example_data/queries.fna \
  -b 4 \
  -o results.txt \
  -t 0.8
```

Output example:
```
0 0:50 2:50 
1 2:50 0:44 1:44 
2 
3 2:50 
```

### Found-in-index mode (experimental)

By default every k-mer of a read that contains an indexed finimizer votes, including k-mers
that are not in the index themselves (a k-mer with a sequencing error usually still contains a
finimizer). With `--found-in-index`, only the read's k-mers that occur in the index vote, and
the threshold is taken against their number: a color is reported if `score / n_pos >= t`, as
real numbers, and its score is at least 1, where `n_pos` is the number of the read's k-mers
found in the index.

The mode needs a membership test for k-mers: a filter built from the reference files and kept
in its own file beside the index. The index itself is unchanged, and the default mode does not
use the filter.

```bash
# Build the filter once: writes <INDEX_PREFIX>.fmin.filter
./finimap build-filter -i <FILE_LIST> -o <INDEX_PREFIX> --threads 8

# Query in found-in-index mode
./finimap search-fmin -i <INDEX_PREFIX> -q <QUERY_FILE> -o <OUTPUT_FILE> -t 0.9 --found-in-index
```

| `build-filter` parameter | Default | Description |
|---|---|---|
| `-i, --input` | Required | The list of reference files given to `finimizer_matrix` |
| `-o, --out-file` | Required | Index filename prefix; the filter is written to `PREFIX.fmin.filter` |
| `-k` | 31 | K-mer length (must match the index) |
| `--bits-per-kmer` | 16 | Size of the filter. 10 gives about 1% false positives, 12 about 0.5%, 16 about 0.2%, 24 about 0.07% |
| `--hashes` | 8 | Bit positions set per k-mer. Use about 6 with 10 bits, 8 with 16, 10 with 24 |
| `--distinct-kmers` | estimated | Number of distinct k-mers, if known; skips the pass that estimates it |
| `--threads` | 1 | Number of threads |

The filter is a blocked Bloom filter: a k-mer of the index always passes, and a k-mer that is
not in the index passes with the false-positive probability. Such a k-mer is counted in
`n_pos` and votes for the colors of whatever finimizer it contains, which need not include the
read's true genomes. With thresholds near 1 this matters: one such k-mer is enough to drop a
genome at `-t 1`. Use a larger filter for high thresholds.

A filter false positive usually stands alone: the k-mers on both sides of it fail. K-mers
that are truly in the index come in runs, between sequencing errors. `--min-run R` counts a
k-mer that passes only if it lies in a run of at least `R` consecutive passing k-mers of the
read (the ends of the read, and k-mers with a character other than A, C, G, T, end a run). With
`--min-run 2`, nearly all false positives are dropped, at the cost of the rare true k-mer that
stands alone.

Other options of `search-fmin` for this mode:

- `--filter FILE` uses a filter file other than `PREFIX.fmin.filter`.
- `--exact-kmers FILE_LIST` uses, instead of the filter file, the exact set of k-mers of the
  reference files, built at query time. It is for testing on small collections.
- For an exact membership test on a large collection, `build-filter --exact` writes an exact
  k-mer table instead of a Bloom filter (16 to 32 bytes per k-mer; `--filter-file` names the
  file), and `--filter` uses it like any filter.
- `--print-denominator` (in either mode) writes, after each read's number, the number its
  scores are thresholded against.

The programs in `dev/checks/` compare the output of both modes with exact k-mer counts, and the
filter with the exact set.

### Limits

- k can be at most 31.
- An index can hold at most 65,535 colors.


## Additional Information

- **Alphabet**: The code works with DNA sequences: {A, C, G, T} (case-insensitive)
- **Reverse complement**: Use `--reverse` flag in `finimizer_matrix build` to index both forward and reverse complement strands
- **Metagenomic data**: Use `--meta true` in `build-fmin` for collections of multiple genomes
- **Large databases**: For very large reference collections (>100GB), consider:
  - Using higher `-p` values for faster queries
  - Using the sparse representation
