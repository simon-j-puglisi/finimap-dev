# finimizer_matrix

Computes the colored finimizers of a collection of genomes and writes them to disk. This is stage 1 of building a [Finimap](https://github.com/simon-j-puglisi/finimap-dev) index; Finimap's `build-fmin` turns the output (a `.cfm` file) into a queryable index.

For each input file (one color), it finds the shortest unique finimizer of every k-mer with respect to the SBWT of the whole collection, and records which colors each finimizer occurs in.

This fork (branch `low-memory-build`) is based on `master_tlen_sparse` of [jnalanko/finimizer_matrix](https://github.com/jnalanko/finimizer_matrix). It writes the same `.cfm` files with much less memory and time:

- Each color keeps only a sorted, gap-encoded list of the SBWT positions of its finimizers, instead of a bit vector over all k-mers of the collection, and the finimizers × colors matrix is filled from these lists at the end. Input files are loaded one at a time instead of all being held in memory.
- The finimizer of each k-mer is picked with a sliding-window minimum.
- The input files for the SBWT build are read by background threads.

On 3,498 *E. coli* genomes (k = 31, 32 threads, `-m 64`), the build took 734 s and 66.9 GB, against about 103 minutes and 326 GB before, and wrote a byte-identical `.cfm` file.

## Building

Requires a Rust toolchain:

```bash
cargo build -r
```

The executable is `target/release/finimizer_matrix`.

## Usage

### `build`

```bash
target/release/finimizer_matrix build -i refs.list -o out.cfm -d tmpdir -k 31 [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-i, --input FILE` | required | A file listing the input FASTA/FASTQ files, one per line. Each file is one color |
| `-o, --output FILE` | required | The output `.cfm` file |
| `-d, --temp-dir DIR` | required | Directory for temporary files of the SBWT build |
| `-k K` | required | k-mer length, at most 31 for Finimap |
| `-t, --n-threads N` | 4 | Number of threads |
| `-m, --mem-gb N` | 4 | Memory budget, in GB, for building the SBWT. Peak memory also includes the SBWT itself and the per-color lists, so it can be well above this figure |
| `--reverse` | off | Also index the reverse complements of the input sequences (normal for Finimap) |
| `--sparse` | off | Store color sets as lists of integers. The memory reductions above apply to the default (dense) build, which is the one used and tested with Finimap |

Progress and the time of each phase are logged to stderr (set `RUST_LOG` to change the level).

### `verify`

```bash
target/release/finimizer_matrix verify -i out.cfm -f refs.list -s index.sbwt -l index.lcs [-t N] [--reverse]
```

Checks a `.cfm` file against an SBWT and its LCS array: that every finimizer is shortest unique, and that the colors are correct. It needs the SBWT and LCS files of the same collection, built separately.
