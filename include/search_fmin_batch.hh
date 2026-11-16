#pragma once

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <cxxopts.hpp>

#include "ColoredFinimizers.hh"
#include "Colors_queries.hh"
#include "Print_output_queries.hh"
#include "query_reader.hh"

using namespace std;

// A batch is closed after this many characters even if it holds fewer reads than the batch size
constexpr uint64_t MAX_BATCH_CHARS = 16000000;

// Queries the reads of one file and writes one result line per read. Returns the number of reads.
//
// A color is reported for a read if its score is at least floor(found * t), and at least 1,
// where found is the number of the read's k-mers that have a finimizer in the index.
template <typename out_stream_t>
uint64_t run_fmin_queries_batch(const string &infile, out_stream_t &out, const CompressedColoredFinimizers &index, const float t, const uint64_t batch_size)
{
    QueryReader reader(infile);
    ResultWriter writer(out);
    BatchWorkspace ws;
    vector<string> reads;
    vector<int32_t> scores;
    vector<pair<int64_t, uint64_t>> runs;
    const uint64_t k = index.get_k();

    uint64_t n_reads_done = 0;
    uint64_t n_too_short = 0;
    string next_read;
    bool have_next = reader.next(next_read);
    while (have_next)
    {
        // Fill the batch reads[0..n), reusing the strings of the previous batch
        uint64_t n = 0, chars = 0;
        while (have_next && n < batch_size && (n == 0 || chars + next_read.size() <= MAX_BATCH_CHARS))
        {
            if (next_read.size() > UINT32_MAX)
                throw runtime_error(infile + ": read " + to_string(n_reads_done + n) + " is too long (over 4,294,967,295 characters)");
            if (n == reads.size())
                reads.emplace_back();
            std::swap(reads[n], next_read);
            chars += reads[n].size();
            n_too_short += reads[n].size() < k;
            n++;
            have_next = reader.next(next_read);
        }

        find_finimizers(reads, n, index.f8, index.pt, index.color_set_ids, k, ws,
                        [&](const uint64_t i, vector<int64_t> &color_set_ids)
                        {
                            const uint64_t found = color_set_ids.size();
                            count_colors(color_set_ids, index.CCS, index.n_colors, scores, runs);
                            const int64_t threshold = static_cast<int64_t>(found * t);
                            writer.write(n_reads_done + i, scores, max<int64_t>(threshold, 1), found);
                        });
        n_reads_done += n;
    }
    writer.flush();

    cerr << infile << ": " << n_reads_done << " reads";
    if (n_too_short > 0)
        cerr << " (" << n_too_short << " shorter than k = " << k << ", reported with no colors)";
    cerr << endl;
    return n_reads_done;
}

// Returns the number of queries executed
uint64_t run_fmin_queries(const vector<string> &infiles, const optional<vector<string>> &outfiles, const CompressedColoredFinimizers &index, const float t, const uint64_t batch_size)
{
    if (outfiles.has_value() && infiles.size() != outfiles.value().size())
    {
        throw runtime_error("Number of input and output files does not match (" + to_string(infiles.size()) + " vs " + to_string(outfiles.value().size()) + ")");
    }

    uint64_t n_queries_run = 0;
    for (uint64_t i = 0; i < infiles.size(); i++)
    {
        if (outfiles.has_value())
        {
            const string &outfile = outfiles.value()[i];
            ofstream out(outfile, ios::binary);
            if (!out)
                throw runtime_error("Could not open " + outfile + " for writing");
            n_queries_run += run_fmin_queries_batch(infiles[i], out, index, t, batch_size);
            out.close();
            if (!out)
                throw runtime_error("Could not write " + outfile);
        }
        else
        { // To stdout
            n_queries_run += run_fmin_queries_batch(infiles[i], cout, index, t, batch_size);
        }
    }
    return n_queries_run;
}

int search_fmin_batch(int argc, char **argv)
{
    cxxopts::Options options(argv[0], "Pseudoalign reads: for each read, report the colors that share enough of its k-mers' finimizers.");

    options.add_options()
        ("o,out-file", "Output filename, or stdout if not given. One line per read: the read's number (from 0), then color:score pairs, highest score first.", cxxopts::value<string>())
        ("i,index-file", "Index filename prefix (the file is PREFIX.fmin).", cxxopts::value<string>())
        ("q,query-file", "The queries in FASTA or FASTQ format, possibly gzipped. Reads may have different lengths. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("t", "Threshold, a fraction from 0 to 1. A color is reported if its score is at least floor(t * found), and at least 1, where found is the number of the read's k-mers that have a finimizer in the index.", cxxopts::value<float>()->default_value("0"))
        ("b", "Batch size, in reads", cxxopts::value<uint64_t>()->default_value("1000"))
        ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help"))
    {
        cerr << options.help() << endl;
        exit(1);
    }
    if (!opts.count("index-file"))
        throw runtime_error("No index given (-i)");
    if (!opts.count("query-file"))
        throw runtime_error("No query file given (-q)");

    const float t = opts["t"].as<float>();
    if (!(t >= 0 && t <= 1))
        throw runtime_error("The threshold (-t) must be between 0 and 1");
    const uint64_t batch_size = opts["b"].as<uint64_t>();
    if (batch_size == 0)
        throw runtime_error("The batch size (-b) must be at least 1");

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> query_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if (multi_file)
    {
        query_files = readlines(queryfile);
    }
    else
    {
        query_files = {queryfile};
    }
    for (string file : query_files)
        check_readable(file);

    // Interpret output file
    optional<vector<string>> output_files;
    if (opts.count("out-file"))
    {
        string outfile = opts["out-file"].as<string>();
        if (multi_file)
        {
            output_files = readlines(outfile);
        }
        else
        {
            output_files = {outfile};
        }
        for (string file : output_files.value())
            check_writable(file);
    }
    else
    {
        cerr << "No output file given, writing to stdout" << endl;
    }

    string index_prefix = opts["index-file"].as<string>();

    cerr << "Loading index..." << endl;
    CompressedColoredFinimizers index;
    index.load(index_prefix);
    cerr << "Index loaded" << endl;

    const uint64_t n_queries = run_fmin_queries(query_files, output_files, index, t, batch_size);
    cerr << "Total number of reads: " << n_queries << endl;

    return 0;
}
