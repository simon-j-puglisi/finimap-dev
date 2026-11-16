#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "CompressedColorSets.hh"

using namespace std;

// BIT_EXPAND[v] has, in its byte j, bit j of v.
static constexpr std::array<uint64_t, 256> make_bit_expand_table()
{
    std::array<uint64_t, 256> table{};
    for (uint64_t v = 0; v < 256; v++)
    {
        uint64_t x = 0;
        for (uint64_t j = 0; j < 8; j++)
            x |= ((v >> j) & 1) << (8 * j);
        table[v] = x;
    }
    return table;
}
static constexpr std::array<uint64_t, 256> BIT_EXPAND = make_bit_expand_table();

// Counts, for every color, how many of a read's finimizers have that color.
//
// Color sets come in three forms (see CompressedColorSets): sparse sets as lists of colors,
// very dense sets as lists of the missing colors, and the rest as bitmaps. Lists are added to
// the counts one color at a time. Bitmaps are added eight colors at a time: each byte of a
// bitmap is spread into eight 8-bit counters held in one 64-bit word (`lanes`), and the lanes
// are emptied into the counts before any of them can exceed 255.
class ColorCounter
{
    vector<int32_t> counts;                // One count per color, padded to a multiple of 64
    vector<uint64_t> lanes;                // Eight 8-bit counters per word, for colors 8i .. 8i+7
    uint64_t pending = 0;                  // Upper bound on the value of any lane
    vector<pair<int64_t, uint64_t>> runs;  // (color set id, number of finimizers with it)

    // Adds every lane, times `scale`, to the counts, and clears the lanes.
    void flush(const int32_t scale)
    {
        if (pending == 0)
            return;
        int32_t *out = counts.data();
        for (size_t i = 0; i < lanes.size(); i++)
        {
            const uint64_t v = lanes[i];
            if (v != 0)
            {
                for (size_t j = 0; j < 8; j++)
                    out[8 * i + j] += static_cast<int32_t>((v >> (8 * j)) & 255) * scale;
                lanes[i] = 0;
            }
        }
        pending = 0;
    }

    // Adds `weight` (at most 255) to the lanes of the colors of one bitmap row.
    void add_row_to_lanes(const uint64_t *data, const uint64_t row, const uint64_t n_colors, const uint64_t weight)
    {
        const uint64_t pos = row * n_colors;
        const uint64_t *ptr = data + pos / 64;
        const uint64_t offset = pos % 64;
        const uint64_t n_words = (n_colors + 63) / 64;
        for (uint64_t w = 0; w < n_words; w++)
        {
            const uint64_t n_bits = std::min<uint64_t>(64, n_colors - w * 64);
            uint64_t word = ptr[w] >> offset;
            if (offset + n_bits > 64)
                word |= ptr[w + 1] << (64 - offset);
            if (n_bits < 64)
                word &= (1ULL << n_bits) - 1;
            if (word == 0)
                continue;
            uint64_t *lane = lanes.data() + 8 * w;
            for (uint64_t j = 0; j < 8; j++)
                lane[j] += BIT_EXPAND[(word >> (8 * j)) & 255] * weight;
        }
    }

    void add_bitmap(const uint64_t *data, const uint64_t row, const uint64_t n_colors, const uint64_t freq)
    {
        if (freq > 255)
        {
            // Too heavy for the 8-bit lanes: add the row once and scale it on the way out
            flush(1);
            add_row_to_lanes(data, row, n_colors, 1);
            pending = 1;
            flush(static_cast<int32_t>(freq));
            return;
        }
        if (pending + freq > 255)
            flush(1);
        add_row_to_lanes(data, row, n_colors, freq);
        pending += freq;
    }

public:
    // Counts the colors of the finimizers whose color set ids are ids[0 .. n_ids), which is
    // sorted in place. Returns the counts, one per color, valid until the next call.
    const int32_t *count(int64_t *ids, const size_t n_ids, const CompressedColorSets &CCS, const uint64_t n_colors)
    {
        const uint64_t n_words = (n_colors + 63) / 64;
        counts.assign(n_words * 64, 0);
        lanes.assign(n_words * 8, 0);
        pending = 0;
        if (n_ids == 0)
            return counts.data();

        // Each distinct color set once, with its number of occurrences
        std::sort(ids, ids + n_ids);
        runs.clear();
        for (size_t i = 0; i < n_ids;)
        {
            size_t j = i + 1;
            while (j < n_ids && ids[j] == ids[i])
                ++j;
            runs.emplace_back(ids[i], j - i);
            i = j;
        }

        const uint64_t *data = CCS.getBV().data();
        const vector<uint16_t> &L = CCS.getL();
        const DeltaSet &EF = CCS.getEF();

        // ids below sparse_count are sparse sets, then very dense sets up to dense_count,
        // then bitmaps; the runs are sorted, so the three kinds come one after another
        const int64_t sparse_count = CCS.sparse_count;
        const int64_t dense_count = CCS.dense_count;
        int32_t *out = counts.data();

        size_t i = 0;
        for (; i < runs.size(); i++)
        {
            const auto &[id, freq] = runs[i];
            if (!(id < sparse_count && id > 0))
                break;
            // In sparse mode, start is always > 0
            const size_t end = EF.get_start(id); // exclusive end
            size_t start = EF.get_start(id - 1); // inclusive start
            while (start < end)
                out[L[start++]] += freq;
        }

        // Very dense sets are stored as their complement: `dense` counts the finimizers with
        // such a set, and the colors missing from each are subtracted here
        uint64_t dense = 0;
        for (; i < runs.size(); i++)
        {
            const auto &[id, freq] = runs[i];
            if (!(id < dense_count))
                break;
            dense += freq;
            const size_t end = EF.get_start(id);                  // exclusive end
            size_t start = (id > 0) ? EF.get_start(id - 1) : 0;   // inclusive start
            while (start < end)
                out[L[start++]] -= freq;
        }

        for (; i < runs.size(); i++)
        {
            const auto &[id, freq] = runs[i];
            add_bitmap(data, id - dense_count, n_colors, freq);
        }
        flush(1);

        if (dense != 0)
        {
            for (uint64_t c = 0; c < n_colors; c++)
                out[c] += dense;
        }
        return out;
    }
};
