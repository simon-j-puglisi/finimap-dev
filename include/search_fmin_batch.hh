#pragma once

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
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

// Everything one thread needs to handle a batch of reads.
struct QueryWorker
{
    vector<string> reads;     // reads[0 .. n_reads) is the batch; the strings are reused
    uint64_t n_reads = 0;
    uint64_t first_read = 0;  // Number of the batch's first read in the query file
    BatchWorkspace ws;
    ColorCounter counter;
    ResultFormatter formatter;
    string error;             // Set if handling the batch failed

    // Finds the finimizers of the batch's reads, counts colors and formats the result lines.
    //
    // A color is reported for a read if its score is at least floor(found * t), and at least
    // 1, where found is the number of the read's k-mers that have a finimizer in the index.
    void run(const CompressedColoredFinimizers &index, const float t)
    {
        try
        {
            find_finimizers(reads, n_reads, index.f8, index.pt, index.loaded_color_set_ids, index.get_k(), ws);
            for (uint64_t i = 0; i < n_reads; i++)
            {
                const uint64_t t0 = now_ns();
                int64_t *ids = ws.finimizers.data() + ws.finimizers_start[i];
                const uint64_t found = ws.finimizers_start[i + 1] - ws.finimizers_start[i];
                const int32_t *scores = counter.count(ids, found, index.CCS, index.n_colors);
                const uint64_t t1 = now_ns();
                const int64_t threshold = static_cast<int64_t>(found * t);
                formatter.write(first_read + i, scores, index.n_colors, max<int64_t>(threshold, 1), found);
                ws.times.count += t1 - t0;
                ws.times.format += now_ns() - t1;
            }
        }
        catch (const std::exception &e)
        {
            error = e.what();
        }
    }
};

struct QuerySummary
{
    uint64_t n_reads = 0;
    QueryTimes times;       // Summed over the threads
    uint64_t read_ns = 0;   // Reading the query file
    uint64_t write_ns = 0;  // Writing the output
};

// Queries the reads of one file and writes one result line per read, in the order of the file.
// The reads are taken in batches, and up to n_threads batches are handled at the same time.
template <typename out_stream_t>
QuerySummary run_fmin_queries_batch(const string &infile, out_stream_t &out, const CompressedColoredFinimizers &index, const float t, const uint64_t batch_size, const uint64_t n_threads)
{
    QueryReader reader(infile);
    vector<QueryWorker> workers(n_threads);
    const uint64_t k = index.get_k();
    QuerySummary summary;

    uint64_t n_too_short = 0;
    string next_read;
    bool have_next = reader.next(next_read);
    while (have_next)
    {
        // Fill one batch per thread, reusing the strings of the previous round
        uint64_t t0 = now_ns();
        uint64_t n_active = 0;
        for (QueryWorker &w : workers)
        {
            if (!have_next)
                break;
            w.first_read = summary.n_reads;
            w.n_reads = 0;
            uint64_t chars = 0;
            while (have_next && w.n_reads < batch_size && (w.n_reads == 0 || chars + next_read.size() <= MAX_BATCH_CHARS))
            {
                if (next_read.size() > UINT32_MAX)
                    throw runtime_error(infile + ": read " + to_string(summary.n_reads) + " is too long (over 4,294,967,295 characters)");
                if (w.n_reads == w.reads.size())
                    w.reads.emplace_back();
                std::swap(w.reads[w.n_reads], next_read);
                chars += w.reads[w.n_reads].size();
                n_too_short += w.reads[w.n_reads].size() < k;
                w.n_reads++;
                summary.n_reads++;
                have_next = reader.next(next_read);
            }
            n_active++;
        }
        uint64_t t1 = now_ns();
        summary.read_ns += t1 - t0;

        #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 1) if (n_active > 1)
        for (uint64_t i = 0; i < n_active; i++)
            workers[i].run(index, t);

        t0 = now_ns();
        for (uint64_t i = 0; i < n_active; i++)
        {
            if (!workers[i].error.empty())
                throw runtime_error(workers[i].error);
            out.write(workers[i].formatter.out.data(), workers[i].formatter.out.size());
            workers[i].formatter.out.clear();
        }
        if (!out)
            throw runtime_error("Could not write the output");
        summary.write_ns += now_ns() - t0;
    }
    uint64_t t0 = now_ns();
    out.flush();
    if (!out)
        throw runtime_error("Could not write the output");
    summary.write_ns += now_ns() - t0;
    for (const QueryWorker &w : workers)
        summary.times.add(w.ws.times);

    cerr << infile << ": " << summary.n_reads << " reads";
    if (n_too_short > 0)
        cerr << " (" << n_too_short << " shorter than k = " << k << ", reported with no colors)";
    cerr << endl;
    return summary;
}

// Runs the queries of every input file. Returns the totals.
QuerySummary run_fmin_queries(const vector<string> &infiles, const optional<vector<string>> &outfiles, const CompressedColoredFinimizers &index, const float t, const uint64_t batch_size, const uint64_t n_threads)
{
    if (outfiles.has_value() && infiles.size() != outfiles.value().size())
    {
        throw runtime_error("Number of input and output files does not match (" + to_string(infiles.size()) + " vs " + to_string(outfiles.value().size()) + ")");
    }

    QuerySummary total;
    for (uint64_t i = 0; i < infiles.size(); i++)
    {
        QuerySummary one;
        if (outfiles.has_value())
        {
            const string &outfile = outfiles.value()[i];
            ofstream out(outfile, ios::binary);
            if (!out)
                throw runtime_error("Could not open " + outfile + " for writing");
            one = run_fmin_queries_batch(infiles[i], out, index, t, batch_size, n_threads);
            out.close();
            if (!out)
                throw runtime_error("Could not write " + outfile);
        }
        else
        { // To stdout
            one = run_fmin_queries_batch(infiles[i], cout, index, t, batch_size, n_threads);
        }
        total.n_reads += one.n_reads;
        total.times.add(one.times);
        total.read_ns += one.read_ns;
        total.write_ns += one.write_ns;
    }
    return total;
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
        ("threads", "Number of threads. The results do not depend on it.", cxxopts::value<uint64_t>()->default_value("1"))
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
    const uint64_t n_threads = opts["threads"].as<uint64_t>();
    if (n_threads == 0 || n_threads > 1024)
        throw runtime_error("The number of threads (--threads) must be between 1 and 1024");

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

    const uint64_t start_ns = now_ns();
    cerr << "Loading index..." << endl;
    CompressedColoredFinimizers index;
    index.load(index_prefix);
    cerr << "Index loaded" << endl;
    const uint64_t loaded_ns = now_ns();

    const QuerySummary summary = run_fmin_queries(query_files, output_files, index, t, batch_size, n_threads);
    const uint64_t done_ns = now_ns();

    const auto seconds = [](const uint64_t ns) { return ns / 1e9; };
    cerr << "Total number of reads: " << summary.n_reads << endl;
    cerr << std::fixed << std::setprecision(2)
         << "Time (s): index load " << seconds(loaded_ns - start_ns)
         << ", search " << seconds(done_ns - loaded_ns)
         << " (wall, " << n_threads << (n_threads == 1 ? " thread)" : " threads)") << endl
         << "Search time by part (s, summed over threads): reading queries " << seconds(summary.read_ns)
         << ", preparing keys " << seconds(summary.times.prepare)
         << ", finimizer lookups " << seconds(summary.times.lookup)
         << ", picking finimizers " << seconds(summary.times.pick)
         << ", counting colors " << seconds(summary.times.count)
         << ", formatting " << seconds(summary.times.format)
         << ", writing " << seconds(summary.write_ns) << endl;
    if (summary.times.scanned_reads > 0)
        cerr << "Reads whose finimizers did not end in order (handled by the slower picker): " << summary.times.scanned_reads << endl;

    return 0;
}
