#pragma once

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <omp.h>

#include "batch_querying.hh"
#include "kmer_filter.hh"
#include "query_reader.hh"

using namespace std;

// Calls f(thread, kmer) for every canonical k-mer of every reference file, with the files
// spread over n_threads threads. Throws the first error of any thread.
template <typename callback_t>
void for_each_reference_kmer(const vector<string> &paths, const uint64_t k, const uint64_t n_threads, callback_t f)
{
    string error;
    #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 1)
    for (uint64_t i = 0; i < paths.size(); i++)
    {
        try
        {
            QueryReader reader(paths[i]);
            string seq;
            const uint64_t thread = omp_get_thread_num();
            while (reader.next(seq))
                for_each_canonical_kmer(seq, k, [&](uint64_t, const uint64_t kmer) { f(thread, kmer); });
        }
        catch (const std::exception &e)
        {
            #pragma omp critical
            if (error.empty())
                error = e.what();
        }
    }
    if (!error.empty())
        throw runtime_error(error);
}

// Builds the k-mer filter of the found-in-index query mode: PREFIX.fmin.filter.
// It needs only the reference files, in two streaming passes: one to estimate the number of
// distinct k-mers, which sets the filter's size, and one to insert them.
int build_filter(int argc, char **argv)
{
    cxxopts::Options options(argv[0], "Build the k-mer filter used by search-fmin --found-in-index.");

    options.add_options()
        ("i,input", "File with one reference file (FASTA or FASTQ, possibly gzipped) per line: the list given to finimizer_matrix.", cxxopts::value<string>())
        ("o,out-file", "Index filename prefix. The filter is written to PREFIX.fmin.filter, beside the index.", cxxopts::value<string>())
        ("k", "k-mer length (must match the index).", cxxopts::value<uint64_t>()->default_value("31"))
        ("bits-per-kmer", "Size of the filter, in bits per distinct k-mer. 10 gives about 1% false positives, 16 about 0.2%, 24 about 0.07%.", cxxopts::value<double>()->default_value("16"))
        ("hashes", "Bit positions set per k-mer. About 6 suits 10 bits per k-mer, 8 suits 16, and 10 suits 24.", cxxopts::value<uint64_t>()->default_value("8"))
        ("distinct-kmers", "Number of distinct canonical k-mers, if known. Skips the pass that estimates it.", cxxopts::value<uint64_t>())
        ("threads", "Number of threads.", cxxopts::value<uint64_t>()->default_value("1"))
        ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);
    if (old_argc == 1 || opts.count("help"))
    {
        cerr << options.help() << endl;
        exit(1);
    }
    if (!opts.count("input"))
        throw runtime_error("No reference list given (-i)");
    if (!opts.count("out-file"))
        throw runtime_error("No index prefix given (-o)");
    const uint64_t k = opts["k"].as<uint64_t>();
    if (k < 1 || k > 31)
        throw runtime_error("k must be between 1 and 31");
    const uint64_t n_threads = opts["threads"].as<uint64_t>();
    if (n_threads == 0 || n_threads > 1024)
        throw runtime_error("The number of threads (--threads) must be between 1 and 1024");
    const string out_path = opts["out-file"].as<string>() + ".fmin.filter";

    vector<string> paths;
    for (const string &line : readlines(opts["input"].as<string>()))
        if (!line.empty())
            paths.push_back(line);
    if (paths.empty())
        throw runtime_error("The reference list is empty");
    check_writable(out_path);

    const uint64_t t0 = now_ns();
    uint64_t n_distinct;
    if (opts.count("distinct-kmers"))
    {
        n_distinct = opts["distinct-kmers"].as<uint64_t>();
    }
    else
    {
        cerr << "Pass 1 of 2: estimating the number of distinct k-mers in " << paths.size() << " files" << endl;
        vector<DistinctCounter> counters(n_threads);
        for_each_reference_kmer(paths, k, n_threads, [&](const uint64_t thread, const uint64_t kmer) { counters[thread].add(kmer); });
        for (uint64_t i = 1; i < n_threads; i++)
            counters[0].merge(counters[i]);
        n_distinct = counters[0].estimate();
    }
    cerr << "Distinct canonical k-mers: about " << n_distinct << endl;
    const uint64_t t1 = now_ns();

    BlockedBloomFilter filter(k, n_distinct, opts["bits-per-kmer"].as<double>(), opts["hashes"].as<uint64_t>());
    cerr << "Inserting the k-mers into a " << filter.description() << endl;
    for_each_reference_kmer(paths, k, n_threads, [&](uint64_t, const uint64_t kmer) { filter.insert(kmer); });
    const uint64_t t2 = now_ns();

    filter.save(out_path);
    cerr << std::fixed << std::setprecision(3)
         << "False-positive rate measured on 1,000,000 random k-mers: " << 100 * filter.measure_false_positive_rate(1000000) << "%" << endl
         << std::setprecision(1)
         << "Time (s): counting " << (t1 - t0) / 1e9 << ", inserting " << (t2 - t1) / 1e9 << ", writing " << (now_ns() - t2) / 1e9 << endl
         << "Wrote " << out_path << " (" << filter.size_in_bytes() << " bytes)" << endl;
    return 0;
}
