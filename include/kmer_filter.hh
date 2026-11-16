#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "query_reader.hh"

// Membership tests for k-mers, for the found-in-index query mode: only the k-mers of a read
// that occur in the index vote for colors, and the threshold is taken against their number.
//
// The index holds both strands of every reference, so a k-mer and its reverse complement are
// either both in it or both absent. K-mers are therefore handled in canonical form: the
// smaller of the k-mer and its reverse complement, as 2-bit codes (A=0, C=1, G=2, T=3) with
// the first character in the highest bits.

// The 2-bit code of a nucleotide, or 4 for any other character.
inline uint8_t nucleotide_code(const char c)
{
    switch (c)
    {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
    }
}

// Calls f(start, canonical_kmer) for every k-mer of seq that consists of nucleotides only,
// in order of start position. k must be at most 32.
template <typename callback_t>
inline void for_each_canonical_kmer(const std::string &seq, const uint64_t k, callback_t f)
{
    const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    uint64_t forward = 0, reverse = 0, run = 0;
    for (uint64_t i = 0; i < seq.size(); i++)
    {
        const uint64_t code = nucleotide_code(seq[i]);
        if (code > 3)
        {
            run = 0;
            continue;
        }
        forward = ((forward << 2) | code) & mask;
        reverse = (reverse >> 2) | ((3 - code) << (2 * (k - 1)));
        if (++run >= k)
            f(i + 1 - k, std::min(forward, reverse));
    }
}

class KmerFilter
{
public:
    virtual ~KmerFilter() = default;
    // Sets pass[i] to 1 if kmers[i] (canonical) may be in the index, else to 0, for i < n.
    // A k-mer of the index always passes; a filter may also pass some others.
    virtual void test(const uint64_t *kmers, const uint64_t n, uint8_t *pass) const = 0;
    virtual std::string description() const = 0;
};

// The exact set of k-mers of a list of reference files. It takes 8 bytes per distinct k-mer
// plus the duplicates while it is built, so it is meant for testing on small collections.
class ExactKmerSet : public KmerFilter
{
    std::vector<uint64_t> kmers; // Sorted, distinct

public:
    // refs_list: a file with one FASTA or FASTQ path per line (as given to finimizer_matrix)
    ExactKmerSet(const std::string &refs_list, const uint64_t k)
    {
        std::ifstream list(refs_list);
        if (!list)
            throw std::runtime_error("Could not open " + refs_list);
        std::string path, seq;
        std::vector<uint64_t> of_file, merged;
        while (std::getline(list, path))
        {
            if (!path.empty() && path.back() == '\r')
                path.pop_back();
            if (path.empty())
                continue;
            of_file.clear();
            QueryReader reader(path);
            while (reader.next(seq))
                for_each_canonical_kmer(seq, k, [&](uint64_t, const uint64_t kmer) { of_file.push_back(kmer); });
            // Merge the file's distinct k-mers into the set
            std::sort(of_file.begin(), of_file.end());
            of_file.erase(std::unique(of_file.begin(), of_file.end()), of_file.end());
            merged.clear();
            merged.reserve(kmers.size() + of_file.size());
            std::set_union(kmers.begin(), kmers.end(), of_file.begin(), of_file.end(), std::back_inserter(merged));
            kmers.swap(merged);
        }
        kmers.shrink_to_fit();
    }

    void test(const uint64_t *query, const uint64_t n, uint8_t *pass) const override
    {
        for (uint64_t i = 0; i < n; i++)
            pass[i] = std::binary_search(kmers.begin(), kmers.end(), query[i]);
    }

    std::string description() const override
    {
        return "exact set of " + std::to_string(kmers.size()) + " canonical k-mers";
    }
};
