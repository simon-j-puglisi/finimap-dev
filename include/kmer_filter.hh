#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "mapped_array.hh"
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

// A 64-bit mixing function (the finalizer of MurmurHash3).
inline uint64_t mix64(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xFF51AFD7ED558CCDULL;
    x ^= x >> 33;
    x *= 0xC4CEB9FE1A85EC53ULL;
    x ^= x >> 33;
    return x;
}

// A blocked Bloom filter over canonical k-mers.
//
// The filter is an array of 512-bit blocks, one cache line each. A k-mer is hashed to one
// block and to a fixed number of bit positions inside it, so a test costs one memory access.
// A k-mer that was inserted always passes; another k-mer passes with a small probability.
//
// File layout (PREFIX.fmin.filter): a 64-byte header of eight 64-bit integers
//   magic, k, number of blocks, bit positions per k-mer, number of k-mers the filter was
//   sized for, and three zeros
// followed by the blocks, 64 bytes each.
class BlockedBloomFilter : public KmerFilter
{
public:
    static constexpr uint64_t MAGIC = 0x31544C494B4D4946ULL; // "FIMKILT1" as bytes
    static constexpr uint64_t HEADER_WORDS = 8;
    static constexpr uint64_t BLOCK_WORDS = 8;

private:
    uint64_t k = 0;
    uint64_t n_blocks = 0;
    uint64_t n_hashes = 0;
    uint64_t sized_for = 0;
    std::vector<uint64_t> owned;                  // The blocks, when the filter is being built
    std::shared_ptr<MappedFile> mapping;          // The file, when the filter is loaded
    const uint64_t *blocks = nullptr;             // 64-byte aligned in both cases

    // The block of a k-mer, and the first position and the step of its bits inside the block
    struct Probe
    {
        uint64_t block, first, step;
    };
    Probe probe(const uint64_t kmer) const
    {
        const uint64_t h = mix64(kmer + 0x9E3779B97F4A7C15ULL);
        const uint64_t g = mix64(h ^ 0xD6E8FEB86659FD93ULL);
        const uint64_t block = static_cast<uint64_t>((static_cast<unsigned __int128>(h) * n_blocks) >> 64);
        return {block, g & 511, ((g >> 9) & 511) | 1}; // An odd step visits distinct positions
    }

public:
    BlockedBloomFilter() = default;

    // An empty filter for about n_kmers distinct k-mers.
    BlockedBloomFilter(const uint64_t k, const uint64_t n_kmers, const double bits_per_kmer, const uint64_t n_hashes)
        : k(k), n_hashes(n_hashes), sized_for(n_kmers)
    {
        if (n_hashes < 1 || n_hashes > 64)
            throw std::runtime_error("The number of bit positions per k-mer must be between 1 and 64");
        if (!(bits_per_kmer >= 1 && bits_per_kmer <= 64))
            throw std::runtime_error("The number of bits per k-mer must be between 1 and 64");
        n_blocks = static_cast<uint64_t>(std::ceil(static_cast<double>(std::max<uint64_t>(n_kmers, 1)) * bits_per_kmer / 512.0));
        n_blocks = std::max<uint64_t>(n_blocks, 1);
        owned.assign(n_blocks * BLOCK_WORDS, 0);
        blocks = owned.data();
    }

    // Loads PATH by mapping it. Throws std::runtime_error if it is not a filter file.
    BlockedBloomFilter(const std::string &path)
    {
        mapping = std::make_shared<MappedFile>(path);
        if (mapping->size() < HEADER_WORDS * 8)
            throw std::runtime_error(path + " is not a Finimap k-mer filter (too short)");
        uint64_t header[HEADER_WORDS];
        std::memcpy(header, mapping->data(), sizeof(header));
        if (header[0] != MAGIC)
            throw std::runtime_error(path + " is not a Finimap k-mer filter");
        k = header[1];
        n_blocks = header[2];
        n_hashes = header[3];
        sized_for = header[4];
        if (n_blocks == 0 || n_hashes < 1 || n_hashes > 64 || n_blocks > (mapping->size() - HEADER_WORDS * 8) / (BLOCK_WORDS * 8) ||
            mapping->size() != HEADER_WORDS * 8 + n_blocks * BLOCK_WORDS * 8)
            throw std::runtime_error(path + " is not a complete Finimap k-mer filter");
        blocks = reinterpret_cast<const uint64_t *>(mapping->data() + HEADER_WORDS * 8);
    }

    uint64_t get_k() const { return k; }
    uint64_t size_in_bytes() const { return HEADER_WORDS * 8 + n_blocks * BLOCK_WORDS * 8; }
    double bits_per_sized_kmer() const { return sized_for ? n_blocks * 512.0 / sized_for : 0; }

    // Inserts a canonical k-mer. Safe to call from several threads at once.
    void insert(const uint64_t kmer)
    {
        const Probe p = probe(kmer);
        uint64_t *block = owned.data() + p.block * BLOCK_WORDS;
        uint64_t pos = p.first;
        for (uint64_t i = 0; i < n_hashes; i++)
        {
            const uint64_t mask = 1ULL << (pos & 63);
            if ((__atomic_load_n(&block[pos >> 6], __ATOMIC_RELAXED) & mask) == 0)
                __atomic_fetch_or(&block[pos >> 6], mask, __ATOMIC_RELAXED);
            pos = (pos + p.step) & 511;
        }
    }

    bool contains(const uint64_t kmer) const
    {
        const Probe p = probe(kmer);
        const uint64_t *block = blocks + p.block * BLOCK_WORDS;
        uint64_t pos = p.first;
        for (uint64_t i = 0; i < n_hashes; i++)
        {
            if ((block[pos >> 6] & (1ULL << (pos & 63))) == 0)
                return false;
            pos = (pos + p.step) & 511;
        }
        return true;
    }

    // Tests a batch. The blocks of the next k-mers are fetched while the current ones are tested.
    void test(const uint64_t *kmers, const uint64_t n, uint8_t *pass) const override
    {
        constexpr uint64_t AHEAD = 16;
        for (uint64_t i = 0; i < n; i++)
        {
            if (i + AHEAD < n)
                __builtin_prefetch(blocks + probe(kmers[i + AHEAD]).block * BLOCK_WORDS);
            pass[i] = contains(kmers[i]);
        }
    }

    std::string description() const override
    {
        char text[160];
        snprintf(text, sizeof(text), "blocked Bloom filter, %.1f MB, %.1f bits per k-mer, %llu bit positions per k-mer",
                 size_in_bytes() / 1e6, bits_per_sized_kmer(), static_cast<unsigned long long>(n_hashes));
        return text;
    }

    // The fraction of n random k-mers that pass. Random k-mers are almost surely not in the
    // filter, so this measures its false-positive rate.
    double measure_false_positive_rate(const uint64_t n) const
    {
        const uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
        uint64_t x = 0x243F6A8885A308D3ULL, passed = 0;
        for (uint64_t i = 0; i < n; i++)
        {
            x += 0x9E3779B97F4A7C15ULL;
            passed += contains(mix64(x) & mask);
        }
        return static_cast<double>(passed) / n;
    }

    void save(const std::string &path) const
    {
        std::ofstream out(path, std::ios::binary);
        if (!out)
            throw std::runtime_error("Could not open " + path + " for writing");
        const uint64_t header[HEADER_WORDS] = {MAGIC, k, n_blocks, n_hashes, sized_for, 0, 0, 0};
        out.write(reinterpret_cast<const char *>(header), sizeof(header));
        out.write(reinterpret_cast<const char *>(blocks), n_blocks * BLOCK_WORDS * 8);
        out.close();
        if (!out)
            throw std::runtime_error("Could not write " + path);
    }
};

// Estimates the number of distinct values in a stream (HyperLogLog with 2^16 registers,
// standard error about 0.4%). Estimators of separate streams can be merged.
class DistinctCounter
{
    static constexpr uint64_t BITS = 16;
    std::vector<uint8_t> registers;

public:
    DistinctCounter() : registers(1ULL << BITS, 0) {}

    void add(const uint64_t value)
    {
        const uint64_t h = mix64(value + 0x9E3779B97F4A7C15ULL);
        const uint64_t index = h >> (64 - BITS);
        const uint64_t rest = (h << BITS) | (1ULL << (BITS - 1)); // Never zero
        const uint8_t rank = static_cast<uint8_t>(__builtin_clzll(rest) + 1);
        if (rank > registers[index])
            registers[index] = rank;
    }

    void merge(const DistinctCounter &other)
    {
        for (size_t i = 0; i < registers.size(); i++)
            registers[i] = std::max(registers[i], other.registers[i]);
    }

    uint64_t estimate() const
    {
        const double m = static_cast<double>(registers.size());
        double sum = 0;
        uint64_t zeros = 0;
        for (const uint8_t r : registers)
        {
            sum += std::ldexp(1.0, -static_cast<int>(r));
            zeros += (r == 0);
        }
        const double alpha = 0.7213 / (1.0 + 1.079 / m);
        double e = alpha * m * m / sum;
        if (e <= 2.5 * m && zeros > 0)
            e = m * std::log(m / static_cast<double>(zeros)); // Small-range correction
        return static_cast<uint64_t>(e + 0.5);
    }
};
