#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Fluke8.hh"
#include "PrefTab.hh"

// Finding the finimizers of a batch of reads.
//
// Every read is handled on its own, so the reads of a batch may have any lengths, and a read
// may be shorter than k (it then has no k-mers). A character other than A, C, G, T (in either
// case) breaks the read: no finimizer may span it, and a k-mer that contains it is skipped.

// Lookup keys hold this many characters, two bits each, starting from the top of a 64-bit word.
constexpr uint64_t FINIMIZER_KEY_CHARS = 31;

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

// Buffers that are reused from batch to batch.
struct BatchWorkspace
{
    std::vector<uint64_t> read_start;                      // Offset of each read in the batch, plus one end offset
    std::vector<uint8_t> max_len;                          // Per position: valid characters from there on, capped at k
    std::vector<uint32_t> found;                           // Per position: 1 + rank of the finimizer starting there, or 0
    std::vector<uint8_t> found_len;                        // Per position: length of that finimizer
    std::vector<std::pair<uint64_t, uint32_t>> keys;       // (lookup key, position), in read order
    std::vector<std::pair<uint64_t, uint32_t>> sorted_keys; // The same, grouped by the top 16 bits of the key
    std::vector<int64_t> finimizers;                       // Color set ids of the finimizers of one read
};

// Looks up, for every position of every read, the finimizer that starts there, and then picks
// the finimizer of every k-mer: the shortest one that lies inside the k-mer, ties broken by rank.
// Handles reads[0 .. n_reads). Calls on_read(i, ids) for each read i in order, where ids holds the color set id of the
// finimizer of each of the read's k-mers that has one (k-mers without a finimizer are left out).
template <typename callback_t>
void find_finimizers(const std::vector<std::string> &reads, const uint64_t n_reads, const Fluke8 &f8, const PrefTab &ptab,
                     const std::vector<uint64_t> &color_set_ids, const uint64_t k,
                     BatchWorkspace &ws, callback_t on_read)
{
    // Lay the reads out one after another
    ws.read_start.clear();
    uint64_t total = 0;
    for (uint64_t ri = 0; ri < n_reads; ri++)
    {
        ws.read_start.push_back(total);
        total += reads[ri].size();
    }
    ws.read_start.push_back(total);
    if (total > UINT32_MAX)
        throw std::runtime_error("A batch of reads is too long (over 4,294,967,295 characters)");

    ws.max_len.assign(total, 0);
    ws.found.assign(total, 0);
    ws.found_len.assign(total, 0);
    ws.keys.clear();

    for (uint64_t ri = 0; ri < n_reads; ri++)
    {
        const std::string &r = reads[ri];
        const uint64_t len = r.size();
        if (len < k)
            continue; // No k-mers
        const uint64_t off = ws.read_start[ri];

        // Number of valid characters starting at each position, capped at k
        uint64_t run = 0;
        for (uint64_t p = len; p-- > 0;)
        {
            run = (nucleotide_code(r[p]) < 4) ? std::min(run + 1, k) : 0;
            ws.max_len[off + p] = run;
        }

        // The key of position p holds r[p .. p+30], padded with zeros past the end of the read
        uint64_t key = 0;
        for (uint64_t q = 0; q < len + FINIMIZER_KEY_CHARS - 1; q++)
        {
            const uint64_t code = (q < len) ? (nucleotide_code(r[q]) & 3) : 0;
            key = (key << 2) | (code << 2);
            if (q + 1 >= FINIMIZER_KEY_CHARS)
            {
                const uint64_t p = q + 1 - FINIMIZER_KEY_CHARS;
                if (ws.max_len[off + p] > 0)
                    ws.keys.push_back({key, static_cast<uint32_t>(off + p)});
            }
        }
    }

    // Group the keys by their top 16 bits, so that the lookups below touch memory in order
    {
        static thread_local std::vector<uint64_t> bucket_start;
        bucket_start.assign(65536 + 1, 0);
        for (const auto &e : ws.keys)
            bucket_start[(e.first >> 48) + 1]++;
        for (uint64_t i = 0; i < 65536; i++)
            bucket_start[i + 1] += bucket_start[i];
        ws.sorted_keys.resize(ws.keys.size());
        for (const auto &e : ws.keys)
            ws.sorted_keys[bucket_start[e.first >> 48]++] = e;
    }

    // Finimizers in the prefix table have ranks after those in Fluke8
    const uint64_t n_f8 = f8.getn();
    for (const auto &[key, pos] : ws.sorted_keys)
    {
        const std::pair<int64_t, uint64_t> ret = ptab.finiLookup(key, ws.max_len[pos]);
        if (ret.first >= 0)
        {
            ws.found[pos] = static_cast<uint32_t>(ret.first + 1 + n_f8);
            ws.found_len[pos] = static_cast<uint8_t>(ret.second);
        }
    }
    for (const auto &[key, pos] : ws.sorted_keys)
    {
        if (ws.found[pos] != 0)
            continue;
        const std::pair<int64_t, uint64_t> ret = f8.finiLookup(key, ws.max_len[pos]);
        if (ret.first >= 0)
        {
            ws.found[pos] = static_cast<uint32_t>(ret.first + 1);
            ws.found_len[pos] = static_cast<uint8_t>(ret.second);
        }
    }

    // Pick the finimizer of every k-mer
    for (uint64_t ri = 0; ri < n_reads; ri++)
    {
        const uint64_t len = reads[ri].size();
        const uint64_t off = ws.read_start[ri];
        ws.finimizers.clear();
        for (uint64_t start = 0; start + k <= len; start++)
        {
            if (ws.max_len[off + start] < k)
                continue; // The k-mer contains a character that is not a nucleotide
            const uint64_t end = start + k - 1;
            // (length, rank) of the best finimizer so far
            std::pair<uint64_t, uint64_t> best = {k + 1, 0};
            for (uint64_t i = start; i <= end; i++)
            {
                if (ws.found[off + i] == 0)
                    continue;
                const uint64_t f_len = ws.found_len[off + i];
                if (i + f_len - 1 > end)
                    break; // It ends after the k-mer, and so does any finimizer that starts later
                const std::pair<uint64_t, uint64_t> candidate = {f_len, ws.found[off + i] - 1};
                if (candidate < best)
                    best = candidate;
            }
            if (best.first <= k)
            {
                if (best.second >= color_set_ids.size())
                    throw std::runtime_error("The index is inconsistent: a finimizer has no color set");
                ws.finimizers.push_back(color_set_ids[best.second]);
            }
        }
        on_read(ri, ws.finimizers);
    }
}
