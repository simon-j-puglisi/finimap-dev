#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{cmp::max, fs::File, io::{BufRead, BufReader, BufWriter}, ops::Sub, path::PathBuf};
use bitvec::prelude::*;
use clap::{builder::PossibleValuesParser, Parser, Subcommand};
mod colored_finimizers;

use colored_finimizers::ColoredFinimizers;

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

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Memory in gigabytes", short = 'm', long = "mem-gb", default_value = "4")]
        mem_gb: usize
    },
}

fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();
    match args.command {
        Subcommands::Build { input: input_fof, output: out_path, temp_dir, k, n_threads, mem_gb } => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            let index = ColoredFinimizers::new(input_paths.as_slice(), k, n_threads, mem_gb, &temp_dir);
            index.serialize(&mut BufWriter::new(File::create(&out_path).unwrap()));
        },
    } 
}