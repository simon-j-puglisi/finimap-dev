#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{cmp::max, fs::File, io::{BufRead, BufReader, BufWriter}, ops::Sub, path::PathBuf};
use bitvec::prelude::*;
use clap::{builder::PossibleValuesParser, Parser, Subcommand};
mod colored_finimizers;

use colored_finimizers::{verify_finimizer_colors, verify_finimizer_shortest_uniqueness, ColoredFinimizers};
use sbwt::{load_sbwt_index_variant, LcsArray, SbwtIndex, SbwtIndexVariant};

#[derive(Parser)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommands,
}

#[derive(Subcommand)]
pub enum Subcommands {
    #[command(arg_required_else_help = true)]
    Build {
        #[arg(help = "A file with one fasta/fastq filename per line", short, long, required = true)]
        input: PathBuf,

        #[arg(help = "Output filename", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Directory for temporary files", short = 'd', long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(short, required = true)]
        k: usize,

        #[arg(help = "Reverse input sequences", long = "reverse")]
        reverse: bool,

        #[arg(help = "Store color sets as sparse (= lists of integers)", long = "sparse")]
        sparse: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Memory in gigabytes", short = 'm', long = "mem-gb", default_value = "4")]
        mem_gb: usize
    },

    #[command(arg_required_else_help = true)]
    Verify {
        #[arg(help = "A colored finimizers index file", short = 'i', long, required = true)]
        input: PathBuf,

        #[arg(help = "Input sequence file-of-files, one per color", short, long, required = true)]
        fof: PathBuf,

        #[arg(help = "SBWT file", short, long, required = true)]
        sbwt: PathBuf,

        #[arg(help = "LCS file", short, long, required = true)]
        lcs: PathBuf,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Reverse input sequences", long = "reverse")]
        reverse: bool,
    },
}

fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();
    match args.command {
        Subcommands::Build { input: input_fof, output: out_path, temp_dir, k, n_threads, mem_gb, reverse, sparse} => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();

            log::info!("Building for {} colors", input_paths.len());
            let index = ColoredFinimizers::new(input_paths.as_slice(), k, n_threads, mem_gb, &temp_dir, reverse, sparse);

            log::info!("Verifying colex order");
            index.verify_colex_order();
            log::info!("Colex order ok");

            log::info!("Serializing index");
            index.serialize(&mut BufWriter::new(File::create(&out_path).unwrap()));
        },
        Subcommands::Verify { input: input_path, fof, sbwt: sbwt_path, lcs: lcs_path, n_threads, reverse} => {
            log::info!("Loading ColoredFinimizers");
            let cf = ColoredFinimizers::load(&mut BufReader::new(File::open(input_path).unwrap()));

            log::info!("Loading SBWT");
            let sbwt = load_sbwt_index_variant(&mut BufReader::new(File::open(sbwt_path).unwrap())).unwrap();

            log::info!("Loading LCS");
            let lcs = LcsArray::load(&mut BufReader::new(File::open(lcs_path).unwrap())).unwrap();

            log::info!("Parsing input file-of-files");
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();

            let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt;

            log::info!("Verifying that finimizers are shortest unique with respect the the SBWT");
            verify_finimizer_shortest_uniqueness(&cf, &sbwt, n_threads);

            log::info!("Verify that the colors of the finimizers are correct");
            verify_finimizer_colors(&cf, &sbwt, &lcs, n_threads, &input_paths, reverse);

            log::info!("Verification OK");
        }
    } 
}