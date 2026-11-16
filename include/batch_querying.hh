#pragma once

#include <algorithm>
#include <chrono>
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

// Time spent in the parts of a search, in nanoseconds.
struct QueryTimes
{
    uint64_t prepare = 0; // Building and grouping the lookup keys
    uint64_t lookup = 0;  // Finimizer lookups
    uint64_t pick = 0;    // Picking the finimizer of each k-mer
    uint64_t count = 0;   // Counting colors
    uint64_t format = 0;  // Thresholding and formatting the result lines
    uint64_t scanned_reads = 0; // Reads whose finimizers were picked by the scanning version

    void add(const QueryTimes &o)
    {
        prepare += o.prepare;
        lookup += o.lookup;
        pick += o.pick;
        count += o.count;
        format += o.format;
        scanned_reads += o.scanned_reads;
    }
};

inline uint64_t now_ns()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// A finimizer that may still be the smallest one inside some k-mer
struct PickCandidate
{
    uint64_t len, rank, pos;
};

// Buffers that are reused from batch to batch.
struct BatchWorkspace
{
    std::vector<uint64_t> read_start;                      // Offset of each read in the batch, plus one end offset
    std::vector<uint8_t> max_len;                          // Per position: valid characters from there on, capped at k
    std::vector<uint32_t> found;                           // Per position: 1 + rank of the finimizer starting there, or 0
    std::vector<uint8_t> found_len;                        // Per position: length of that finimizer
    std::vector<std::pair<uint64_t, uint32_t>> keys;       // (lookup key, position), in read order
    std::vector<std::pair<uint64_t, uint32_t>> sorted_keys; // The same, grouped by the top 16 bits of the key
    std::vector<int64_t> finimizers;                       // Color set ids of the finimizers of the k-mers, read after read
    std::vector<uint64_t> finimizers_start;                // Where each read's ids start in `finimizers`, plus one end
    std::vector<PickCandidate> candidates;                 // Sliding window of candidates, for picking
    std::vector<uint64_t> bucket_start;                    // For grouping the keys
    QueryTimes times;
};

// Picks the finimizer of every k-mer of one read: the shortest finimizer that lies inside the
// k-mer, ties broken by rank. Appends the color set id of each k-mer's finimizer to `out`
// (k-mers without a finimizer, or with a character that is not a nucleotide, are left out).
// `off` is the offset of the read in the batch arrays, `len` its length.
//
// This is the reference version: it scans the k positions of every k-mer.
inline void pick_finimizers_scan(const BatchWorkspace &ws, const uint64_t off, const uint64_t len, const uint64_t k,
                                 const std::vector<uint64_t> &color_set_ids, std::vector<int64_t> &out)
{
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
            out.push_back(color_set_ids[best.second]);
        }
    }
}

// The same result in one pass over the read, with a sliding window minimum.
//
// It relies on the finimizers of a read ending in the order in which they start, which holds
// because a finimizer is a shortest unique string: one that started later and ended earlier
// would lie inside the non-unique part of the other. Then the finimizers inside a k-mer are
// those that start in it and end by its last position, they enter the window in order as the
// k-mer moves right, and they leave it in the same order. If the ends are not in order (which
// would mean a wrong lookup), the read is handled by the scanning version instead.
inline void pick_finimizers(BatchWorkspace &ws, const uint64_t off, const uint64_t len, const uint64_t k,
                            const std::vector<uint64_t> &color_set_ids, std::vector<int64_t> &out)
{
    if (len < k)
        return;

    uint64_t last_end = 0;
    for (uint64_t p = 0; p < len; p++)
    {
        if (ws.found[off + p] == 0)
            continue;
        const uint64_t end = p + ws.found_len[off + p] - 1;
        if (end < last_end)
        {
            ws.times.scanned_reads++;
            pick_finimizers_scan(ws, off, len, k, color_set_ids, out);
            return;
        }
        last_end = end;
    }

    // candidates[head .. tail) holds the finimizers that lie inside the current k-mer or a
    // later one and are smaller than every one after them, in order of position
    if (ws.candidates.size() < len)
        ws.candidates.resize(len);
    PickCandidate *candidates = ws.candidates.data();
    uint64_t head = 0, tail = 0;
    uint64_t next = 0; // The next position to consider
    for (uint64_t end = k - 1; end < len; end++)
    {
        const uint64_t start = end + 1 - k;
        while (next <= end)
        {
            const uint64_t f = ws.found[off + next];
            if (f != 0)
            {
                const uint64_t f_len = ws.found_len[off + next];
                if (next + f_len - 1 > end)
                    break; // Not inside this k-mer yet; neither is any later one
                const uint64_t rank = f - 1;
                while (tail > head && (candidates[tail - 1].len > f_len || (candidates[tail - 1].len == f_len && candidates[tail - 1].rank >= rank)))
                    tail--;
                candidates[tail++] = {f_len, rank, next};
            }
            next++;
        }
        while (head < tail && candidates[head].pos < start)
            head++;
        if (head < tail && ws.max_len[off + start] >= k)
        {
            const uint64_t rank = candidates[head].rank;
            if (rank >= color_set_ids.size())
                throw std::runtime_error("The index is inconsistent: a finimizer has no color set");
            out.push_back(color_set_ids[rank]);
        }
    }
}

// Looks up, for every position of every read, the finimizer that starts there, and then picks
// the finimizer of every k-mer. Handles reads[0 .. n_reads). Afterwards the color set ids of
// the finimizers of read i's k-mers are ws.finimizers[ws.finimizers_start[i] .. ws.finimizers_start[i+1]).
inline void find_finimizers(const std::vector<std::string> &reads, const uint64_t n_reads, const Fluke8 &f8, const PrefTab &ptab,
                            const std::vector<uint64_t> &color_set_ids, const uint64_t k, BatchWorkspace &ws)
{
    uint64_t t0 = now_ns();

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
    ws.bucket_start.assign(65536 + 1, 0);
    for (const auto &e : ws.keys)
        ws.bucket_start[(e.first >> 48) + 1]++;
    for (uint64_t i = 0; i < 65536; i++)
        ws.bucket_start[i + 1] += ws.bucket_start[i];
    ws.sorted_keys.resize(ws.keys.size());
    for (const auto &e : ws.keys)
        ws.sorted_keys[ws.bucket_start[e.first >> 48]++] = e;

    uint64_t t1 = now_ns();
    ws.times.prepare += t1 - t0;

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

    uint64_t t2 = now_ns();
    ws.times.lookup += t2 - t1;

    ws.finimizers.clear();
    ws.finimizers_start.clear();
    for (uint64_t ri = 0; ri < n_reads; ri++)
    {
        ws.finimizers_start.push_back(ws.finimizers.size());
        pick_finimizers(ws, ws.read_start[ri], reads[ri].size(), k, color_set_ids, ws.finimizers);
    }
    ws.finimizers_start.push_back(ws.finimizers.size());

    ws.times.pick += now_ns() - t2;
}
