#pragma once

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
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

// How reads are scored and thresholded.
struct QueryMode
{
    // Default mode: every k-mer with a finimizer in the index votes. A color is reported if
    // its score is at least floor(found * threshold_float), and at least 1, where found is
    // the number of k-mers that voted.
    float threshold_float = 0;

    // Found-in-index mode (filter is set): only k-mers that pass the filter vote. A color is
    // reported if score / n_pos >= threshold, as real numbers, and its score is at least 1,
    // where n_pos is the number of k-mers that passed the filter.
    const KmerFilter *filter = nullptr;
    double threshold = 0;

    uint64_t min_run = 1; // Found-in-index mode: minimum run of consecutive passing k-mers

    // Relative rule (either mode), if not negative: a color must also have a score of at least
    // top - relative, where top is the read's highest score
    int64_t relative = -1;

    bool print_denominator = false; // Write found or n_pos after the read's number
};

// The smallest score s >= 1 with s / n_pos >= threshold as real numbers, for n_pos >= 1.
// It is larger than n_pos if no score can pass.
inline int64_t min_score_for_ratio(const double threshold, const uint64_t n_pos)
{
    const double n = static_cast<double>(n_pos);
    int64_t s = static_cast<int64_t>(std::ceil(threshold * n));
    s = std::max<int64_t>(s, 1);
    // The product is rounded, so settle on the comparison that defines the rule
    while (s > 1 && static_cast<double>(s - 1) / n >= threshold)
        s--;
    while (s <= static_cast<int64_t>(n_pos) && static_cast<double>(s) / n < threshold)
        s++;
    return s;
}

// Everything one thread needs to handle a batch of reads.
struct QueryWorker
{
    vector<string> reads;     // reads[0 .. n_reads) is the batch; the strings are reused
    uint64_t n_reads = 0;
    uint64_t first_read = 0;  // Number of the batch's first read in the query file
    BatchWorkspace ws;
    ColorCounter counter;
    MetaColorCounter meta_counter; // index format 4
    ResultFormatter formatter;
    string error;             // Set if handling the batch failed

    // Finds the finimizers of the batch's reads, counts colors and formats the result lines.
    void run(const CompressedColoredFinimizers &index, const QueryMode &mode)
    {
        try
        {
            find_finimizers(reads, n_reads, index.f8, index.pt, index.loaded_color_set_ids, index.get_k(), mode.filter, ws, mode.min_run);
            for (uint64_t i = 0; i < n_reads; i++)
            {
                const uint64_t t0 = now_ns();
                int64_t *ids = ws.finimizers.data() + ws.finimizers_start[i];
                const uint64_t found = ws.finimizers_start[i + 1] - ws.finimizers_start[i];
                const int32_t *scores = index.format >= 4 ? meta_counter.count(ids, found, index.MCS) : counter.count(ids, found, index.CCS, index.n_colors);
                const uint64_t t1 = now_ns();

                // No score can exceed the number of k-mers that voted
                uint64_t denominator;
                int64_t min_score;
                if (mode.filter)
                {
                    denominator = ws.n_passed[i];
                    min_score = (denominator == 0) ? 1 : min_score_for_ratio(mode.threshold, denominator);
                }
                else
                {
                    denominator = found;
                    min_score = max<int64_t>(static_cast<int64_t>(found * mode.threshold_float), 1);
                }
                if (mode.relative >= 0)
                {
                    int64_t top = 0;
                    for (uint64_t c = 0; c < index.n_colors; c++)
                        top = max<int64_t>(top, scores[c]);
                    min_score = max<int64_t>(min_score, top - mode.relative);
                }
                formatter.write(first_read + i, scores, index.n_colors, min_score, found, mode.print_denominator ? static_cast<int64_t>(denominator) : -1);
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
QuerySummary run_fmin_queries_batch(const string &infile, out_stream_t &out, const CompressedColoredFinimizers &index, const QueryMode &mode, const uint64_t batch_size, const uint64_t n_threads)
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
            workers[i].run(index, mode);

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
QuerySummary run_fmin_queries(const vector<string> &infiles, const optional<vector<string>> &outfiles, const CompressedColoredFinimizers &index, const QueryMode &mode, const uint64_t batch_size, const uint64_t n_threads)
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
            one = run_fmin_queries_batch(infiles[i], out, index, mode, batch_size, n_threads);
            out.close();
            if (!out)
                throw runtime_error("Could not write " + outfile);
        }
        else
        { // To stdout
            one = run_fmin_queries_batch(infiles[i], cout, index, mode, batch_size, n_threads);
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
        ("t", "Threshold, a fraction from 0 to 1. A color is reported if its score is at least floor(t * found), and at least 1, where found is the number of the read's k-mers that have a finimizer in the index.", cxxopts::value<string>()->default_value("0"))
        ("found-in-index", "Found-in-index mode: only the read's k-mers that occur in the index vote, and a color is reported if score / n_pos >= t, where n_pos is the number of those k-mers. Needs the k-mer filter made by build-filter (PREFIX.fmin.filter).")
        ("filter", "The k-mer filter file to use in found-in-index mode, if it is not PREFIX.fmin.filter.", cxxopts::value<string>())
        ("exact-kmers", "For testing the found-in-index mode: instead of the filter file, use the exact set of k-mers of the reference files listed in this file (one path per line, as given to finimizer_matrix). Takes 8 bytes per distinct k-mer.", cxxopts::value<string>())
        ("relative", "Relative rule: report only the colors whose score is within this many k-mers of the read's highest score (score >= top - DELTA). Combines with -t: a color must pass both, so with the default -t 0 this is the relative rule alone. Works in both modes.", cxxopts::value<int64_t>())
        ("min-run", "Found-in-index mode: a k-mer that passes the filter counts only if it lies in a run of at least this many consecutive passing k-mers of the read. 1 counts every passing k-mer.", cxxopts::value<uint64_t>()->default_value("1"))
        ("print-denominator", "Write, after each read's number, the number its scores are thresholded against: found, or n_pos in found-in-index mode.")
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

    // The threshold is read both as a float, which the default mode has always used, and as a
    // double, for the real-number comparison of the found-in-index mode
    QueryMode mode;
    {
        const string text = opts["t"].as<string>();
        char *end_f = nullptr, *end_d = nullptr;
        mode.threshold_float = std::strtof(text.c_str(), &end_f);
        mode.threshold = std::strtod(text.c_str(), &end_d);
        if (text.empty() || *end_f != '\0' || *end_d != '\0' || !(mode.threshold >= 0 && mode.threshold <= 1))
            throw runtime_error("The threshold (-t) must be a number between 0 and 1");
    }
    mode.print_denominator = opts.count("print-denominator") > 0;
    const bool found_in_index = opts.count("found-in-index") > 0;
    if (opts.count("relative"))
    {
        mode.relative = opts["relative"].as<int64_t>();
        if (mode.relative < 0)
            throw runtime_error("--relative must be at least 0");
    }
    mode.min_run = opts["min-run"].as<uint64_t>();
    if (mode.min_run == 0)
        throw runtime_error("--min-run must be at least 1");
    if (mode.min_run > 1 && !found_in_index)
        throw runtime_error("--min-run is only used with --found-in-index");
    if ((opts.count("exact-kmers") || opts.count("filter")) && !found_in_index)
        throw runtime_error("--exact-kmers and --filter are only used with --found-in-index");
    if (opts.count("exact-kmers") && opts.count("filter"))
        throw runtime_error("Give either --filter or --exact-kmers, not both");
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

    std::unique_ptr<KmerFilter> filter;
    if (found_in_index)
    {
        if (opts.count("exact-kmers"))
        {
            cerr << "Building the exact k-mer set..." << endl;
            filter = std::make_unique<ExactKmerSet>(opts["exact-kmers"].as<string>(), index.get_k());
        }
        else
        {
            const string path = opts.count("filter") ? opts["filter"].as<string>() : index_prefix + ".fmin.filter";
            if (!std::ifstream(path))
                throw runtime_error("--found-in-index needs the k-mer filter " + path + ", which does not exist. Make it with: finimap build-filter -i <reference list> -o " + index_prefix);
            uint64_t filter_k = 0;
            filter = load_kmer_filter(path, filter_k);
            if (filter_k != static_cast<uint64_t>(index.get_k()))
                throw runtime_error("The k-mer filter " + path + " was built for k = " + to_string(filter_k) + ", but the index has k = " + to_string(index.get_k()));
        }
        mode.filter = filter.get();
        cerr << "Found-in-index mode, filter: " << filter->description() << endl;
    }
    const uint64_t filter_ns = now_ns();

    const QuerySummary summary = run_fmin_queries(query_files, output_files, index, mode, batch_size, n_threads);
    const uint64_t done_ns = now_ns();

    const auto seconds = [](const uint64_t ns) { return ns / 1e9; };
    cerr << "Total number of reads: " << summary.n_reads << endl;
    cerr << std::fixed << std::setprecision(2)
         << "Time (s): index load " << seconds(loaded_ns - start_ns)
         << ", search " << seconds(done_ns - filter_ns)
         << " (wall, " << n_threads << (n_threads == 1 ? " thread)" : " threads)") << endl
         << "Search time by part (s, summed over threads): reading queries " << seconds(summary.read_ns)
         << ", preparing keys " << seconds(summary.times.prepare)
         << ", finimizer lookups " << seconds(summary.times.lookup)
         << ", k-mer filter " << seconds(summary.times.filter)
         << ", picking finimizers " << seconds(summary.times.pick)
         << ", counting colors " << seconds(summary.times.count)
         << ", formatting " << seconds(summary.times.format)
         << ", writing " << seconds(summary.write_ns) << endl;
    if (found_in_index)
    {
        cerr << "K-mers that passed the filter: " << summary.times.kmers_passed;
        if (mode.min_run > 1)
            cerr << ", of which in runs shorter than " << mode.min_run << " (not counted): " << summary.times.kmers_in_short_runs;
        cerr << endl;
    }
    if (summary.times.scanned_reads > 0)
        cerr << "Reads whose finimizers did not end in order (handled by the slower picker): " << summary.times.scanned_reads << endl;

    return 0;
}
