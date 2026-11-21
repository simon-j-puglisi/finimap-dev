#pragma once

#include <cxxopts.hpp>
#include "ColoredFinimizers.hh"

int build_fmin(int argc, char **argv)
{
    cxxopts::Options options(argv[0], "Compress colored Finimizers.");

    options.add_options()
        ("i,index-file", "ColoredFinimizers file.", cxxopts::value<string>())
        ("o,out-file", "Output index filename prefix.", cxxopts::value<string>()) // as input
        ("p,p_len", "Finimizers prefix length.", cxxopts::value<int64_t>()->default_value(std::to_string(10)))
        ("x", "Short-long finimizer threshold: finimizers of at most this many characters are short. It sets the layout of the index, not the results. The bucket table for short finimizers takes about 4^x / 228 bytes (16: 19 MB, 18: 302 MB; in index format 1, 4^x / 64 bytes).",  cxxopts::value<int64_t>()->default_value(std::to_string(16)))
	    ("k", "k-mer length.", cxxopts::value<int64_t>()->default_value(std::to_string(31)))
        ("m, meta", "Metagenome.", cxxopts::value<bool>()->default_value("false"))
        ("index-format", "Index format to write: 3 (current), or 2 or 1 (those of earlier versions of Finimap, for comparison).", cxxopts::value<int64_t>()->default_value(std::to_string(INDEX_FORMAT)))
        ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help"))
    {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    string out_prefix = opts["out-file"].as<string>();

    string indexfile = opts["index-file"].as<string>();
    int64_t prefix_len = opts["p"].as<int64_t>();
    int64_t short_long_t = opts["x"].as<int64_t>();
    int64_t format = opts["index-format"].as<int64_t>();
    if (format < 1 || format > (int64_t)INDEX_FORMAT)
        throw std::runtime_error("--index-format must be between 1 and " + std::to_string(INDEX_FORMAT));
    if (short_long_t < 1 || short_long_t > 19)
        throw std::runtime_error("-x must be between 1 and 19");

    bool meta = opts["m"].as<bool>();
    int64_t k = opts["k"].as<int64_t>();
    if (k < 1 || k > (int64_t)MAX_K)
        throw std::runtime_error("k must be between 1 and " + std::to_string(MAX_K));

    ColoredFinimizers cf;
    ifstream in(indexfile, ios::binary);
    if (!in)
        throw std::runtime_error("Could not open " + indexfile);
    cf.load(in);

    CompressedColoredFinimizers ccf(std::move(cf), prefix_len, short_long_t, k, meta, format); // Search currently works only with 31

    ccf.serialize(out_prefix);

    return 0;
}
