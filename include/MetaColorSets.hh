#pragma once

// Meta-color sets (index format 4). The colors are put in an order in which colors that occur
// in similar color sets are close (color_clustering.hh), and split into blocks of at most 256
// colors. The part of a color set that falls in a block is a partial set. Each block stores its
// distinct partial sets once, and each color set is stored as the ids of its partial sets, one
// fixed-width field per block. Most color sets share most of their partial sets with others, so
// this takes several times less space than storing each color set on its own.
//
// Counting the colors of the finimizers of a read adds, block by block, the colors of each
// distinct partial set once, multiplied by the number of the read's finimizers that use it.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "color_clustering.hh"
#include "Colors_queries.hh" // BIT_EXPAND
#include "mapped_array.hh"
#include "packed_ints.hh"

class MetaColorSets
{
public:
    static constexpr uint64_t MAX_BLOCK = 256;            // colors in a block
    static constexpr uint64_t MAX_CLUSTER_COLORS = 16384; // above this, the input order is kept
    static constexpr uint64_t SAMPLE_SETS = 65536;        // sets sampled to compare colors

    uint64_t n_colors = 0, n_sets = 0, n_blocks = 0, row_bits = 0, n_partials = 0;
    std::vector<uint16_t> order;       // order[p]: the color at position p
    std::vector<uint64_t> block_start; // n_blocks + 1 positions
    std::vector<uint64_t> width;       // bits of the field of each block (0: one partial set)
    std::vector<uint64_t> field_off;   // bit offset of each block's field in a row
    std::vector<uint64_t> pbase;       // n_blocks + 1: the first partial set of each block
    PackedInts poffset;                // n_partials + 1 byte offsets into pdata
    PackedInts pkind;                  // the kind of each partial set: LIST, ABSENT or BITMAP
    Array<uint8_t> pdata;              // the partial sets
    Array<uint8_t> rows;               // n_sets rows of row_bits bits, then 8 bytes of padding

    // Derived, not stored: for counting, each block starts at a multiple of 8 in a padded
    // space of n_padded positions
    std::vector<uint64_t> padded_start; // n_blocks + 1
    uint64_t n_padded = 0;

    void prepare()
    {
        padded_start.assign(n_blocks + 1, 0);
        for (uint64_t b = 0; b < n_blocks; b++)
            padded_start[b + 1] = padded_start[b] + (block_size(b) + 7) / 8 * 8;
        n_padded = padded_start[n_blocks];
    }

private:
    using Key = std::array<uint64_t, MAX_BLOCK / 64>;
    struct KeyHash
    {
        size_t operator()(const Key &k) const
        {
            uint64_t h = 0x9e3779b97f4a7c15ULL;
            for (uint64_t w : k)
            {
                h ^= w + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                h *= 0xbf58476d1ce4e5b9ULL;
            }
            return (size_t)(h ^ (h >> 31));
        }
    };

    static uint64_t bits_for_count(uint64_t n) // bits to tell apart n values (0 if n <= 1)
    {
        uint64_t w = 0;
        while ((uint64_t(1) << w) < n)
            w++;
        return w;
    }

    static void write_bits(uint8_t *buf, uint64_t pos, uint64_t value, uint64_t w)
    {
        for (uint64_t i = 0; i < w; i++, pos++)
            if ((value >> i) & 1)
                buf[pos >> 3] |= uint8_t(1) << (pos & 7);
    }

public:
    // A partial set is stored in the smallest of three forms: the list of its positions in the
    // block (one byte each), the list of the positions it lacks, or a bitmap of the block
    static constexpr uint64_t LIST = 0, ABSENT = 1, BITMAP = 2;

    uint64_t block_size(uint64_t b) const { return block_start[b + 1] - block_start[b]; }
    static uint64_t bitmap_bytes(uint64_t size) { return (size + 7) / 8; }

    // The id of the partial set of color set s in block b
    uint64_t field(uint64_t s, uint64_t b) const
    {
        const uint64_t w = width[b];
        if (w == 0)
            return 0;
        const uint64_t pos = s * row_bits + field_off[b];
        uint64_t v;
        std::memcpy(&v, reinterpret_cast<const uint8_t *>(rows.data()) + (pos >> 3), 8);
        return (v >> (pos & 7)) & ((uint64_t(1) << w) - 1);
    }

    // Builds the structure from n_sets color sets; get_set(i) returns the words of set i
    // (n_colors bits). Set i gets color set id i.
    template <typename GetSet>
    void build(const uint64_t sets, const uint64_t colors, GetSet get_set)
    {
        n_sets = sets;
        n_colors = colors;
        if (n_colors == 0 || n_colors > 65535)
            throw std::runtime_error("MetaColorSets: the number of colors must be 1 to 65535");
        // The order of the colors and the blocks
        std::vector<uint64_t> sizes;
        order.resize(n_colors);
        if (n_colors >= 2 && n_colors <= MAX_CLUSTER_COLORS && n_sets > 0)
        {
            uint64_t SW = 0;
            const std::vector<uint64_t> col = colorclust::columns(n_sets, n_colors, SAMPLE_SETS, get_set, SW);
            const colorclust::Tree tree = colorclust::average_linkage(col, n_colors, SW);
            for (uint64_t p = 0; p < n_colors; p++)
                order[p] = (uint16_t)tree.order[p];
            sizes = colorclust::tree_blocks(tree, MAX_BLOCK);
            std::cerr << "Meta-colors: colors ordered by average-linkage clustering" << std::endl;
        }
        else
        {
            for (uint64_t p = 0; p < n_colors; p++)
                order[p] = (uint16_t)p;
            sizes = colorclust::fixed_blocks(n_colors, MAX_BLOCK);
            std::cerr << "Meta-colors: colors kept in input order" << std::endl;
        }
        n_blocks = sizes.size();
        block_start.assign(n_blocks + 1, 0);
        for (uint64_t b = 0; b < n_blocks; b++)
            block_start[b + 1] = block_start[b] + sizes[b];

        // The distinct partial sets of each block, and each set's ids
        std::vector<std::vector<uint32_t>> ids(n_blocks);
        std::vector<std::vector<Key>> partials(n_blocks);
#pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            std::unordered_map<Key, uint32_t, KeyHash> seen;
            const Key empty{};
            seen.emplace(empty, 0); // id 0 is always the empty partial set
            partials[b].push_back(empty);
            ids[b].resize(n_sets);
            const uint64_t start = block_start[b], size = sizes[b];
            for (uint64_t i = 0; i < n_sets; i++)
            {
                const uint64_t *r = get_set(i);
                Key k{};
                for (uint64_t off = 0; off < size; off++)
                {
                    const uint64_t c = order[start + off];
                    if ((r[c / 64] >> (c % 64)) & 1)
                        k[off / 64] |= uint64_t(1) << (off % 64);
                }
                auto [it, fresh] = seen.try_emplace(k, (uint32_t)partials[b].size());
                if (fresh)
                    partials[b].push_back(k);
                ids[b][i] = it->second;
            }
        }

        // Field widths and the rows
        width.assign(n_blocks, 0);
        field_off.assign(n_blocks, 0);
        pbase.assign(n_blocks + 1, 0);
        row_bits = 0;
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            width[b] = bits_for_count(partials[b].size());
            if (width[b] > 32)
                throw std::runtime_error("MetaColorSets: too many partial sets in a block");
            field_off[b] = row_bits;
            row_bits += width[b];
            pbase[b + 1] = pbase[b] + partials[b].size();
        }
        n_partials = pbase[n_blocks];
        const uint64_t row_bytes = (n_sets * row_bits + 7) / 8 + 8;
        rows.resize(row_bytes);
        {
            std::vector<uint8_t> buf(row_bytes, 0);
#pragma omp parallel for schedule(static)
            for (uint64_t i0 = 0; i0 < n_sets; i0 += 8)
            {
                // 8 rows at a time, so that no two threads write the same byte
                for (uint64_t i = i0; i < std::min(n_sets, i0 + 8); i++)
                    for (uint64_t b = 0; b < n_blocks; b++)
                        write_bits(buf.data(), i * row_bits + field_off[b], ids[b][i], width[b]);
            }
            std::memcpy(&rows[0], buf.data(), row_bytes);
        }
        ids.clear();

        // The partial sets
        std::vector<uint64_t> offsets(n_partials + 1, 0);
        std::vector<uint8_t> kinds(n_partials, LIST);
        std::vector<uint8_t> data;
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            const uint64_t size = sizes[b], bb = bitmap_bytes(size);
            for (const Key &k : partials[b])
            {
                const uint64_t g = pbase[b] + (&k - partials[b].data());
                offsets[g] = data.size();
                uint64_t m = 0;
                for (uint64_t w : k)
                    m += __builtin_popcountll(w);
                if (m < bb || size - m < bb)
                {
                    const bool absent = m >= bb;
                    kinds[g] = absent ? ABSENT : LIST;
                    for (uint64_t off = 0; off < size; off++)
                        if ((((k[off / 64] >> (off % 64)) & 1) != 0) != absent)
                            data.push_back((uint8_t)off);
                }
                else
                {
                    kinds[g] = BITMAP;
                    for (uint64_t j = 0; j < bb; j++)
                        data.push_back((uint8_t)(k[j / 8] >> (8 * (j % 8))));
                }
            }
        }
        offsets[n_partials] = data.size();
        for (int i = 0; i < 8; i++)
            data.push_back(0); // padding for 8-byte loads
        poffset.build(offsets.data(), offsets.size(), PackedInts::bits_for(offsets[n_partials]));
        pkind.build(kinds.data(), kinds.size(), 2);
        pdata.resize(data.size());
        std::memcpy(&pdata[0], data.data(), data.size());
        prepare();
        std::cerr << "Meta-colors: " << n_blocks << " blocks of " << *std::min_element(sizes.begin(), sizes.end()) << "-"
                  << *std::max_element(sizes.begin(), sizes.end()) << " colors, " << n_partials << " distinct partial sets ("
                  << (data.size() + poffset.size_in_bytes() + pkind.size_in_bytes()) / 1e6 << " MB), " << row_bits << " bits per color set ("
                  << row_bytes / 1e6 << " MB)" << std::endl;
    }

    uint64_t size_in_bytes() const
    {
        return 6 * 8 + order.size() * 2 + (block_start.size() + width.size() + pbase.size()) * 8 + poffset.size_in_bytes() + pkind.size_in_bytes() + 48 +
               pdata.size() + 8 + rows.size() + 8;
    }

    void serialize(std::ostream &out) const
    {
        auto u64 = [&](uint64_t v)
        { out.write(reinterpret_cast<const char *>(&v), 8); };
        u64(n_colors);
        u64(n_sets);
        u64(n_blocks);
        u64(row_bits);
        u64(n_partials);
        out.write(reinterpret_cast<const char *>(order.data()), order.size() * 2);
        for (uint64_t v : block_start)
            u64(v);
        for (uint64_t v : width)
            u64(v);
        for (uint64_t v : pbase)
            u64(v);
        poffset.serialize(out);
        pkind.serialize(out);
        u64(pdata.size());
        out.write(reinterpret_cast<const char *>(pdata.data()), pdata.size());
        u64(rows.size());
        out.write(reinterpret_cast<const char *>(rows.data()), rows.size());
    }

    void load(std::istream &in)
    {
        auto u64 = [&]()
        {
            uint64_t v = 0;
            in.read(reinterpret_cast<char *>(&v), 8);
            return v;
        };
        n_colors = u64();
        n_sets = u64();
        n_blocks = u64();
        row_bits = u64();
        n_partials = u64();
        if (!in || n_colors == 0 || n_colors > 65535 || n_blocks == 0 || n_blocks > n_colors)
            throw std::runtime_error("The index is inconsistent (meta-color sets)");
        order.resize(n_colors);
        in.read(reinterpret_cast<char *>(order.data()), n_colors * 2);
        block_start.resize(n_blocks + 1);
        for (auto &v : block_start)
            v = u64();
        width.resize(n_blocks);
        for (auto &v : width)
            v = u64();
        pbase.resize(n_blocks + 1);
        for (auto &v : pbase)
            v = u64();
        field_off.assign(n_blocks, 0);
        uint64_t bits = 0;
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            if (width[b] > 32 || block_start[b + 1] <= block_start[b] || block_start[b + 1] - block_start[b] > MAX_BLOCK)
                throw std::runtime_error("The index is inconsistent (meta-color blocks)");
            field_off[b] = bits;
            bits += width[b];
        }
        if (bits != row_bits || block_start[n_blocks] != n_colors || pbase[n_blocks] != n_partials)
            throw std::runtime_error("The index is inconsistent (meta-color layout)");
        for (uint16_t c : order)
            if (c >= n_colors)
                throw std::runtime_error("The index is inconsistent (meta-color order)");
        poffset.load(in);
        pkind.load(in);
        const uint64_t pdata_size = u64();
        pdata.load(in, pdata_size);
        const uint64_t rows_size = u64();
        if (in && rows_size != (n_sets * row_bits + 7) / 8 + 8)
            throw std::runtime_error("The index is inconsistent (meta-color rows)");
        rows.load(in, rows_size);
        if (in && (poffset.size() != n_partials + 1 || pkind.size() != n_partials || poffset[n_partials] + 8 != pdata_size))
            throw std::runtime_error("The index is inconsistent (meta-color partial sets)");
        prepare();
    }
};

// Counts the colors of a read's finimizers with meta-color sets. Block by block, the read's
// color sets are mapped to their partial sets, and each distinct partial set is added once,
// with the number of finimizers that have it as its weight. Lists add to the counts one color
// at a time; bitmaps are added eight positions at a time through 8-bit lanes (as in
// ColorCounter), in the padded order of the blocks, and the lanes are emptied into the counts
// before any can exceed 255.
class MetaColorCounter
{
    std::vector<int32_t> counts;
    std::vector<uint64_t> lanes;   // zero between reads: emptying them clears them
    std::vector<uint32_t> pending; // per block: upper bound on its lanes; zero between reads
    std::vector<std::pair<int64_t, uint32_t>> runs;
    // The distinct partial sets of one block of a read and their weights: a small
    // open-addressing table (key: partial id + 1, 0 when free), cleared through `used`
    std::vector<std::pair<uint32_t, uint32_t>> table;
    std::vector<uint32_t> used;

    void flush(const MetaColorSets &M, uint64_t b, int32_t scale)
    {
        if (pending[b] == 0)
            return;
        const uint16_t *ord = M.order.data() + M.block_start[b];
        const uint64_t size = M.block_size(b);
        uint64_t *lane = lanes.data() + M.padded_start[b] / 8;
        for (uint64_t i = 0; i < (size + 7) / 8; i++)
        {
            const uint64_t v = lane[i];
            if (v == 0)
                continue;
            const uint64_t n = std::min<uint64_t>(8, size - 8 * i);
            for (uint64_t j = 0; j < n; j++)
                counts[ord[8 * i + j]] += (int32_t)((v >> (8 * j)) & 255) * scale;
            lane[i] = 0;
        }
        pending[b] = 0;
    }

    void add_bitmap_bytes(const MetaColorSets &M, uint64_t b, const uint8_t *ptr, uint64_t len, uint64_t weight)
    {
        uint64_t *lane = lanes.data() + M.padded_start[b] / 8;
        for (uint64_t j = 0; j < len; j++)
            lane[j] += BIT_EXPAND[ptr[j]] * weight;
    }

    // Adds a partial set with a weight; returns the weight to add to every color of the block
    // (for a partial set stored as the colors it lacks)
    uint64_t add_partial(const MetaColorSets &M, uint64_t b, uint64_t p, uint64_t weight)
    {
        const uint64_t g = M.pbase[b] + p;
        const uint64_t start = M.poffset[g], end = M.poffset[g + 1];
        const uint8_t *ptr = reinterpret_cast<const uint8_t *>(M.pdata.data()) + start;
        const uint64_t kind = M.pkind[g];
        if (kind != MetaColorSets::BITMAP)
        {
            const uint16_t *ord = M.order.data() + M.block_start[b];
            const int32_t w = (int32_t)weight;
            if (kind == MetaColorSets::LIST)
            {
                for (uint64_t j = 0; j < end - start; j++)
                    counts[ord[ptr[j]]] += w;
                return 0;
            }
            for (uint64_t j = 0; j < end - start; j++)
                counts[ord[ptr[j]]] -= w;
            return weight;
        }
        if (weight > 255)
        {
            flush(M, b, 1);
            add_bitmap_bytes(M, b, ptr, end - start, 1);
            pending[b] = 1;
            flush(M, b, (int32_t)weight);
            return 0;
        }
        if (pending[b] + weight > 255)
            flush(M, b, 1);
        add_bitmap_bytes(M, b, ptr, end - start, weight);
        pending[b] += (uint32_t)weight;
        return 0;
    }

public:
    // Counts the colors of the finimizers whose color set ids are ids[0 .. n_ids), which is
    // sorted in place. Returns the counts, one per color, valid until the next call.
    const int32_t *count(int64_t *ids, const size_t n_ids, const MetaColorSets &M)
    {
        counts.assign(((M.n_colors + 63) / 64) * 64, 0);
        if (n_ids == 0)
            return counts.data();
        if (lanes.size() != M.n_padded / 8)
            lanes.assign(M.n_padded / 8, 0);
        if (pending.size() != M.n_blocks)
            pending.assign(M.n_blocks, 0);
        std::sort(ids, ids + n_ids);
        runs.clear();
        for (size_t i = 0; i < n_ids;)
        {
            size_t j = i + 1;
            while (j < n_ids && ids[j] == ids[i])
                j++;
            runs.push_back({ids[i], (uint32_t)(j - i)});
            i = j;
        }
        // The table has at least twice as many slots as the read has color sets
        uint64_t tbits = 4;
        while ((uint64_t(1) << tbits) < 2 * runs.size())
            tbits++;
        if (table.size() < (uint64_t(1) << tbits))
            table.assign(uint64_t(1) << tbits, {0, 0});
        while ((uint64_t(1) << tbits) < table.size())
            tbits++;
        const uint64_t mask = table.size() - 1;
        for (uint64_t b = 0; b < M.n_blocks; b++)
        {
            if (M.width[b] == 0)
                continue; // every set has the empty partial set here
            used.clear();
            for (const auto &[s, f] : runs)
            {
                const uint64_t p = M.field((uint64_t)s, b);
                if (p == 0)
                    continue;
                const uint32_t key = (uint32_t)p + 1;
                uint64_t h = ((uint64_t)key * 0x9E3779B97F4A7C15ULL) >> (64 - tbits);
                while (true)
                {
                    auto &slot = table[h];
                    if (slot.first == 0)
                    {
                        slot = {key, f};
                        used.push_back((uint32_t)h);
                        break;
                    }
                    if (slot.first == key)
                    {
                        slot.second += f;
                        break;
                    }
                    h = (h + 1) & mask;
                }
            }
            uint64_t everywhere = 0;
            for (uint32_t h : used)
            {
                everywhere += add_partial(M, b, table[h].first - 1, table[h].second);
                table[h] = {0, 0};
            }
            flush(M, b, 1);
            if (everywhere)
            {
                const uint16_t *ord = M.order.data() + M.block_start[b];
                for (uint64_t j = 0; j < M.block_size(b); j++)
                    counts[ord[j]] += (int32_t)everywhere;
            }
        }
        return counts.data();
    }
};
