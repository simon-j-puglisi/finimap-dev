#pragma once

#include <vector>
#include <unordered_map>
//#include <limits>
#include "CompressedColorSets.hh"


using namespace std;

inline void process_word(uint64_t word, const uint64_t base, vector<int32_t> &results, const uint64_t freq)
{
    while (word)
    {
        // uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        results[base + bit] += freq;
        word &= word - 1; // clear lowest bit
    }
}

inline void read_bv(const uint64_t *data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<int32_t> &results)
{

    const uint64_t *ptr = data + (start * n_colors) / 64;
    uint64_t bit_offset = (start * n_colors) % 64;

    uint64_t color_id = 0;
    uint64_t bits_left = n_colors;

    if (bit_offset != 0)
    {
        // 1. Read the first word
        uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
        uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
        uint64_t word = (*ptr >> bit_offset) & mask;
        process_word(word, color_id, results, freq);

        ptr++;
        color_id += bits_to_read;
        bits_left -= bits_to_read;
    }

    // 2. Read aligned words in btw
    while (bits_left >= 64)
    {
        uint64_t word = *ptr++;
        process_word(word, color_id, results, freq);
        color_id += 64;
        bits_left -= 64;
    }

    // 3. Read the last word (if any)
    if (bits_left > 0)
    {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;
        process_word(word, color_id, results, freq);
    }
}


void read_verydense(const int64_t pos, const uint64_t freq, const DeltaSet &EF, const vector<uint16_t> &L, const uint64_t n_colors, vector<int32_t> &results, uint32_t &dense)
{
    dense += freq;
    const size_t end = EF.get_start(pos); // exclusive end
    size_t start = (pos > 0) ? EF.get_start(pos - 1) : 0; // inclusive start
    while (start < end)
    {
        results[L[start++]] -= freq;
    }
}

void read_colors(const CompressedColorSets &CCS, const uint64_t n_colors, vector<int32_t> &results, const vector<pair<int64_t, uint64_t>> &fmin_v, uint32_t &dense)
{
    const sdsl::bit_vector &BV = CCS.getBV();
    const uint64_t *data = BV.data();
    const vector<uint16_t> &L = CCS.getL();
    const DeltaSet &EF = CCS.getEF();

    // pos < sparse_count; [sparse]
    // sparse_count <= pos < dense_count; [very dense]
    // pos >= dense; [bitmap]

    const uint64_t sparse_count = CCS.sparse_count;
    const uint64_t dense_count = CCS.dense_count;

    // Exploit the fact that the pos are sorted
    uint64_t i;
    for (i = 0; i < fmin_v.size(); i++)
    {
        const auto &[pos, freq] = fmin_v[i];
        // sparse
        if (pos < sparse_count && pos > 0)
        {
            // read from L - in sparse mode, start is always > 0
            const size_t end = EF.get_start(pos); // exclusive end
            size_t start = EF.get_start(pos - 1); // inclusive start
            while (start < end)
            {
                results[L[start++]] += freq;
            }
        }
        else
        {
            break;
        }
    }
    // very dense
    uint64_t j;
    for (j = i; j < fmin_v.size(); j++)
    {
        const auto &[pos, freq] = fmin_v[j];
        if (pos < dense_count)
        {
            // read complementary values from L
            read_verydense(pos, freq, EF, L, n_colors, results, dense);
        }
        else
        {
            break;
        }
    }
    // read from BV
    for (auto i = j; i < fmin_v.size(); i++)
    {
        const auto &[pos, freq] = fmin_v[i];
        uint64_t start = pos - dense_count;
        read_bv(data, start, freq, n_colors, results);
    }
}

// Counts, for every color, how many of a read's finimizers have that color.
// color_set_ids holds the color set id of the finimizer of each k-mer; it is sorted in place.
// results gets one count per color.
inline void count_colors(vector<int64_t> &color_set_ids, const CompressedColorSets &CCS, const uint64_t n_colors,
                         vector<int32_t> &results, vector<pair<int64_t, uint64_t>> &runs)
{
    results.assign(n_colors, 0);
    if (color_set_ids.empty())
        return;

    // Each distinct color set once, with its number of occurrences
    std::sort(color_set_ids.begin(), color_set_ids.end());
    runs.clear();
    for (size_t i = 0; i < color_set_ids.size();)
    {
        size_t j = i + 1;
        while (j < color_set_ids.size() && color_set_ids[j] == color_set_ids[i])
            ++j;
        runs.emplace_back(color_set_ids[i], j - i);
        i = j;
    }

    // Very dense sets are stored as their complement: `dense` counts the finimizers with such
    // a set, and the colors missing from each have already been subtracted in results.
    uint32_t dense = 0;
    read_colors(CCS, n_colors, results, runs, dense);
    if (dense != 0)
    {
        for (auto &r : results)
            r += dense;
    }
}
