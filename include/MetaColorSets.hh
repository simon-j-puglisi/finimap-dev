#pragma once

// Meta-color sets (index format 4). The colors are put in an order in which colors that occur
// in similar color sets are close (color_clustering.hh), and split into blocks of at most 256
// colors. The part of a color set that falls in a block is a partial set. Each block stores its
// distinct partial sets once, and has a default: its most common partial set.
//
// Each color set is stored as a row: a bitmask of the blocks where its partial set is not the
// default, then the ids of those partial sets, each in a canonical Huffman code of its block
// (common partial sets get short codes). The partial sets of a block are numbered in the order
// of their codes, so that a decoded code is the partial set's id. Rows start at offsets kept for
// groups of ROW_GROUP rows plus a short offset within the group.
//
// Partial sets are stored in the smallest of five forms: the list of their positions in the
// block (one byte each), the list of the positions they lack, their runs of consecutive
// positions or the runs of positions they lack (first and last position, a byte each), or a
// bitmap. For each group of 8
// partial sets there is a record of the byte offset of the first, then one byte per partial set
// holding its length in bytes and its form.
//
// Counting the colors of the finimizers of a read adds, block by block, the colors of each
// distinct partial set once, multiplied by the number of the read's finimizers that use it.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <istream>
#include <numeric>
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
    static constexpr uint64_t MAX_CODE = 32;              // longest code
    static constexpr uint64_t TABLE_BITS = 10;            // codes up to this long decode in one lookup
    static constexpr uint64_t ROW_GROUP = 16;             // rows per full offset
    static constexpr uint64_t PART_GROUP = 8;             // partial sets per record
    static constexpr uint64_t RECORD_BYTES = 4 + PART_GROUP;
    static constexpr uint64_t NONE = ~uint64_t(0);

    // The forms of a partial set: a list of its positions in the block, one byte each; the same
    // for the positions it lacks; its runs of consecutive positions, or the runs of positions it
    // lacks, two bytes each (first and last position); a bitmap. A bitmap is the only form that
    // takes as many bytes as the bitmap of the block, so its record stores the form as LIST.
    static constexpr uint64_t LIST = 0, ABSENT = 1, RUNS = 2, ABSENT_RUNS = 3, BITMAP = 4;

    uint64_t n_colors = 0, n_sets = 0, n_blocks = 0, row_bits = 0, n_partials = 0;
    std::vector<uint16_t> order;       // order[p]: the color at position p
    std::vector<uint64_t> block_start; // n_blocks + 1 positions
    std::vector<uint64_t> pbase;       // n_blocks + 1: the first partial set of each block
    std::vector<uint64_t> defaults;    // the default partial set of each block (its last id)
    std::vector<uint64_t> empty_id;    // the id of the empty partial set of each block, or NONE
    std::vector<uint64_t> code_count;  // n_blocks * (MAX_CODE + 1): codes of each length
    Array<uint8_t> precords;           // ceil(n_partials / PART_GROUP) records, see above
    Array<uint8_t> pdata;              // the partial sets, then 8 bytes of padding
    PackedInts row_base;               // the bit offset of every ROW_GROUP-th row
    PackedInts row_rel;                // the offset of each row from its group's first
    Array<uint8_t> rows;               // the rows, then 8 bytes of padding

    // Derived, not stored
    std::vector<uint64_t> padded_start; // n_blocks + 1: blocks at multiples of 8, for counting
    uint64_t n_padded = 0;
    std::vector<uint32_t> table;        // n_blocks << TABLE_BITS: (id << 6) | length, 0 if longer
    std::vector<uint64_t> first_code, first_rank; // n_blocks * (MAX_CODE + 1), canonical codes
    std::vector<uint64_t> limit;        // n_blocks * (MAX_CODE + 1): codes of length <= l, as 32-bit
                                        // values with the first bit highest, are below limit[l]

    uint64_t block_size(uint64_t b) const { return block_start[b + 1] - block_start[b]; }
    static uint64_t bitmap_bytes(uint64_t size) { return (size + 7) / 8; }

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

    static void write_bits(uint8_t *buf, uint64_t pos, uint64_t value, uint64_t w)
    {
        for (uint64_t i = 0; i < w; i++, pos++)
            if ((value >> i) & 1)
                buf[pos >> 3] |= uint8_t(1) << (pos & 7);
    }

    // Huffman code lengths for the given frequencies (all > 0), at most max_len long
    static std::vector<uint32_t> huffman_lengths(std::vector<uint64_t> freq, uint64_t max_len)
    {
        const uint64_t n = freq.size();
        std::vector<uint32_t> len(n, 0);
        if (n == 0)
            return len;
        if (n == 1)
        {
            len[0] = 1;
            return len;
        }
        while (true)
        {
            std::vector<uint32_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::sort(idx.begin(), idx.end(), [&](uint32_t a, uint32_t b)
                      { return freq[a] < freq[b] || (freq[a] == freq[b] && a < b); });
            // Leaves 0..n-1 in increasing weight, internal nodes n..2n-2 in increasing weight
            std::vector<uint64_t> w(2 * n - 1);
            std::vector<uint32_t> parent(2 * n - 1, 0);
            for (uint64_t i = 0; i < n; i++)
                w[i] = freq[idx[i]];
            uint64_t li = 0, ii = n;
            auto take = [&](uint64_t made) -> uint64_t
            {
                if (li < n && (ii >= made || w[li] <= w[ii]))
                    return li++;
                return ii++;
            };
            for (uint64_t k = 0; k + 1 < n; k++)
            {
                const uint64_t made = n + k;
                const uint64_t a = take(made), b = take(made);
                w[made] = w[a] + w[b];
                parent[a] = parent[b] = (uint32_t)made;
            }
            std::vector<uint32_t> depth(2 * n - 1, 0);
            for (uint64_t v = 2 * n - 2; v-- > 0;)
                depth[v] = depth[parent[v]] + 1;
            uint64_t longest = 0;
            for (uint64_t i = 0; i < n; i++)
            {
                len[idx[i]] = depth[i];
                longest = std::max<uint64_t>(longest, depth[i]);
            }
            if (longest <= max_len)
                return len;
            for (auto &f : freq)
                f = (f >> 1) | 1; // flatten and try again
        }
    }

    void derive()
    {
        padded_start.assign(n_blocks + 1, 0);
        for (uint64_t b = 0; b < n_blocks; b++)
            padded_start[b + 1] = padded_start[b] + (block_size(b) + 7) / 8 * 8;
        n_padded = padded_start[n_blocks];
        // Canonical codes: within a length, codes and ids increase together
        const uint64_t L = MAX_CODE + 1;
        first_code.assign(n_blocks * L, 0);
        first_rank.assign(n_blocks * L, 0);
        limit.assign(n_blocks * L, 0);
        table.assign(n_blocks << TABLE_BITS, 0);
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            const uint64_t *count = &code_count[b * L];
            uint64_t code = 0, rank = 0;
            for (uint64_t l = 1; l <= MAX_CODE; l++)
            {
                code = (code + count[l - 1]) << 1;
                first_code[b * L + l] = code;
                first_rank[b * L + l] = rank;
                rank += count[l];
                limit[b * L + l] = (code + count[l]) << (32 - l);
            }
            for (uint64_t l = 1; l <= std::min(MAX_CODE, TABLE_BITS); l++)
                for (uint64_t j = 0; j < count[l]; j++)
                {
                    const uint64_t c = first_code[b * L + l] + j, id = first_rank[b * L + l] + j;
                    uint64_t rev = 0; // the code as it lies in the stream: first bit lowest
                    for (uint64_t t = 0; t < l; t++)
                        rev |= ((c >> (l - 1 - t)) & 1) << t;
                    for (uint64_t x = 0; x < (uint64_t(1) << (TABLE_BITS - l)); x++)
                        table[(b << TABLE_BITS) | rev | (x << l)] = (uint32_t)((id << 6) | l);
                }
        }
    }

    static uint32_t reverse32(uint32_t v)
    {
        v = ((v >> 1) & 0x55555555u) | ((v & 0x55555555u) << 1);
        v = ((v >> 2) & 0x33333333u) | ((v & 0x33333333u) << 2);
        v = ((v >> 4) & 0x0F0F0F0Fu) | ((v & 0x0F0F0F0Fu) << 4);
        return __builtin_bswap32(v);
    }

public:
    // w <= 56 bits of the rows from bit pos
    uint64_t read_bits(uint64_t pos, uint64_t w) const
    {
        uint64_t v;
        std::memcpy(&v, reinterpret_cast<const uint8_t *>(rows.data()) + (pos >> 3), 8);
        return (v >> (pos & 7)) & ((uint64_t(1) << w) - 1);
    }

    // Decodes the id of a partial set of block b at bit pos, and moves pos past it
    uint64_t decode(uint64_t b, uint64_t &pos) const
    {
        const uint32_t e = table[(b << TABLE_BITS) | read_bits(pos, TABLE_BITS)];
        if (e != 0)
        {
            pos += e & 63;
            return e >> 6;
        }
        // A code longer than TABLE_BITS: take the next 32 bits of the stream with the first bit
        // highest, and find the length whose codes include its prefix (canonical codes of a
        // length are consecutive, and below those of the next length shifted up)
        const uint64_t L = MAX_CODE + 1;
        const uint64_t x = reverse32((uint32_t)read_bits(pos, 32));
        // The length: one more than TABLE_BITS for each limit that x reaches, counted without
        // branches
        const uint64_t *lim = &limit[b * L];
        uint64_t l = TABLE_BITS + 1;
        for (uint64_t j = TABLE_BITS + 1; j < MAX_CODE; j++)
            l += x >= lim[j];
        const uint64_t code = x >> (32 - l);
        const uint64_t k = code - first_code[b * L + l];
        if (k >= code_count[b * L + l])
            throw std::runtime_error("The index is inconsistent (meta-color rows)");
        pos += l;
        return first_rank[b * L + l] + k;
    }

    // The bit position of the row of color set s, and prefetches of what finds it and the row
    uint64_t row_position(uint64_t s) const { return row_base[s / ROW_GROUP] + row_rel[s]; }
    void prefetch_position(uint64_t s) const
    {
        __builtin_prefetch(row_base.address(s / ROW_GROUP));
        __builtin_prefetch(row_rel.address(s));
    }
    void prefetch_row(uint64_t pos) const
    {
        const uint8_t *p = reinterpret_cast<const uint8_t *>(rows.data()) + (pos >> 3);
        __builtin_prefetch(p);
        __builtin_prefetch(p + 64);
    }

    // Calls f(b, p) for each block b where the color set whose row starts at bit pos has a
    // partial set p other than the block's default, in increasing order of b
    template <typename F>
    void for_each_difference_at(uint64_t pos, F f) const
    {
        uint64_t idpos = pos + n_blocks;
        for (uint64_t b0 = 0; b0 < n_blocks; b0 += 56)
            for (uint64_t m = read_bits(pos + b0, std::min<uint64_t>(56, n_blocks - b0)); m; m &= m - 1)
            {
                const uint64_t b = b0 + __builtin_ctzll(m);
                f(b, decode(b, idpos));
            }
    }

    template <typename F>
    void for_each_difference(uint64_t s, F f) const { for_each_difference_at(row_position(s), f); }

    // Partial set g (over all blocks), of block b: where its bytes are, how many, and its form
    struct Partial
    {
        const uint8_t *ptr;
        uint64_t len, form;
    };
    const void *partial_record(uint64_t g) const
    {
        return reinterpret_cast<const uint8_t *>(precords.data()) + (g / PART_GROUP) * RECORD_BYTES;
    }
    Partial partial(uint64_t g, uint64_t b) const
    {
        const uint8_t *rec = reinterpret_cast<const uint8_t *>(partial_record(g));
        uint32_t base;
        uint64_t bytes;
        std::memcpy(&base, rec, 4);
        std::memcpy(&bytes, rec + 4, 8);
        const uint64_t j = g % PART_GROUP;
        // The lengths are at most 32: summing the bytes of the record never carries over
        const uint64_t before = j ? (((bytes & 0x3F3F3F3F3F3F3F3FULL) * 0x0101010101010101ULL) >> (8 * (j - 1))) & 255 : 0;
        const uint64_t byte = (bytes >> (8 * j)) & 255;
        const uint64_t len = byte & 63;
        return {reinterpret_cast<const uint8_t *>(pdata.data()) + base + before, len,
                len == (padded_start[b + 1] - padded_start[b]) / 8 ? BITMAP : byte >> 6};
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

        // The distinct partial sets of each block, each set's ids, and the codes
        const uint64_t L = MAX_CODE + 1;
        std::vector<std::vector<uint32_t>> ids(n_blocks);
        std::vector<std::vector<Key>> partials(n_blocks);
        std::vector<std::vector<uint32_t>> code_len(n_blocks); // by new id
        defaults.assign(n_blocks, 0);
        empty_id.assign(n_blocks, NONE);
        code_count.assign(n_blocks * L, 0);
#pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            std::unordered_map<Key, uint32_t, KeyHash> seen;
            std::vector<Key> found;
            std::vector<uint64_t> uses;
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
                auto [it, fresh] = seen.try_emplace(k, (uint32_t)found.size());
                if (fresh)
                {
                    found.push_back(k);
                    uses.push_back(0);
                }
                uses[it->second]++;
                ids[b][i] = it->second;
            }
            if (found.empty())
            {
                // no sets: one empty partial set, the default
                found.push_back(Key{});
                uses.push_back(1);
                seen.emplace(Key{}, 0);
            }
            // The default, then Huffman codes for the others, numbered in code order
            const uint64_t def = std::max_element(uses.begin(), uses.end()) - uses.begin();
            std::vector<uint32_t> others;
            std::vector<uint64_t> freq;
            for (uint64_t p = 0; p < found.size(); p++)
                if (p != def)
                {
                    others.push_back((uint32_t)p);
                    freq.push_back(uses[p]);
                }
            const std::vector<uint32_t> len = huffman_lengths(freq, MAX_CODE);
            std::vector<uint32_t> rank(others.size());
            std::iota(rank.begin(), rank.end(), 0);
            std::sort(rank.begin(), rank.end(), [&](uint32_t a, uint32_t c)
                      { return len[a] < len[c] || (len[a] == len[c] && others[a] < others[c]); });
            std::vector<uint32_t> new_id(found.size());
            partials[b].resize(found.size());
            code_len[b].assign(found.size(), 0);
            for (uint64_t r = 0; r < rank.size(); r++)
            {
                const uint32_t old = others[rank[r]];
                new_id[old] = (uint32_t)r;
                partials[b][r] = found[old];
                code_len[b][r] = len[rank[r]];
                code_count[b * L + len[rank[r]]]++;
            }
            new_id[def] = (uint32_t)others.size();
            partials[b][others.size()] = found[def];
            defaults[b] = others.size();
            for (auto &id : ids[b])
                id = new_id[id];
            auto e = seen.find(Key{});
            if (e != seen.end())
                empty_id[b] = new_id[e->second];
        }
        pbase.assign(n_blocks + 1, 0);
        for (uint64_t b = 0; b < n_blocks; b++)
            pbase[b + 1] = pbase[b] + partials[b].size();
        n_partials = pbase[n_blocks];
        derive();

        // The codes of each block's ids, as they are written: first bit lowest
        std::vector<std::vector<uint32_t>> code_rev(n_blocks);
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            code_rev[b].assign(partials[b].size(), 0);
            for (uint64_t id = 0; id < defaults[b]; id++)
            {
                const uint64_t l = code_len[b][id];
                const uint64_t c = first_code[b * L + l] + (id - first_rank[b * L + l]);
                uint64_t rev = 0;
                for (uint64_t t = 0; t < l; t++)
                    rev |= ((c >> (l - 1 - t)) & 1) << t;
                code_rev[b][id] = (uint32_t)rev;
            }
        }

        // The rows
        std::vector<uint64_t> offsets(n_sets + 1, 0);
        for (uint64_t i = 0; i < n_sets; i++)
        {
            uint64_t bits = n_blocks;
            for (uint64_t b = 0; b < n_blocks; b++)
                if (ids[b][i] != defaults[b])
                    bits += code_len[b][ids[b][i]];
            offsets[i + 1] = offsets[i] + bits;
        }
        row_bits = offsets[n_sets];
        const uint64_t row_bytes = (row_bits + 7) / 8 + 8;
        {
            std::vector<uint8_t> buf(row_bytes, 0);
            for (uint64_t i = 0; i < n_sets; i++)
            {
                const uint64_t pos = offsets[i];
                uint64_t idpos = pos + n_blocks;
                for (uint64_t b = 0; b < n_blocks; b++)
                {
                    const uint32_t id = ids[b][i];
                    if (id == defaults[b])
                        continue;
                    buf[(pos + b) >> 3] |= uint8_t(1) << ((pos + b) & 7);
                    write_bits(buf.data(), idpos, code_rev[b][id], code_len[b][id]);
                    idpos += code_len[b][id];
                }
            }
            rows.resize(row_bytes);
            std::memcpy(&rows[0], buf.data(), row_bytes);
        }
        {
            const uint64_t n_groups = (n_sets + ROW_GROUP - 1) / ROW_GROUP;
            std::vector<uint64_t> base(n_groups), rel(n_sets);
            uint64_t max_rel = 0;
            for (uint64_t i = 0; i < n_sets; i++)
            {
                if (i % ROW_GROUP == 0)
                    base[i / ROW_GROUP] = offsets[i];
                rel[i] = offsets[i] - base[i / ROW_GROUP];
                max_rel = std::max(max_rel, rel[i]);
            }
            row_base.build(base.data(), base.size(), PackedInts::bits_for(row_bits));
            row_rel.build(rel.data(), rel.size(), PackedInts::bits_for(max_rel));
        }
        ids.clear();

        // The partial sets and their records
        const uint64_t n_records = (n_partials + PART_GROUP - 1) / PART_GROUP;
        std::vector<uint8_t> data, records(n_records * RECORD_BYTES, 0);
        std::vector<uint64_t> forms(5, 0);
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            const uint64_t size = sizes[b], bb = bitmap_bytes(size);
            for (uint64_t id = 0; id < partials[b].size(); id++)
            {
                const Key &k = partials[b][id];
                const uint64_t g = pbase[b] + id;
                uint8_t *rec = &records[(g / PART_GROUP) * RECORD_BYTES];
                if (g % PART_GROUP == 0)
                {
                    if (data.size() >= (uint64_t(1) << 32))
                        throw std::runtime_error("MetaColorSets: the partial sets take more than 4 GB");
                    const uint32_t base = (uint32_t)data.size();
                    std::memcpy(rec, &base, 4);
                }
                uint64_t m = 0;
                for (uint64_t w : k)
                    m += __builtin_popcountll(w);
                const uint64_t before = data.size();
                // The smallest of the five forms (a bitmap takes bb bytes, at most 32)
                std::vector<uint32_t> present, absent;
                for (uint64_t off = 0; off < size; off++)
                    (((k[off / 64] >> (off % 64)) & 1) ? present : absent).push_back((uint32_t)off);
                auto runs = [](const std::vector<uint32_t> &pos)
                {
                    std::vector<uint8_t> r;
                    for (uint64_t i = 0; i < pos.size(); i++)
                        if (i == 0 || pos[i] != pos[i - 1] + 1)
                            r.push_back((uint8_t)pos[i]), r.push_back((uint8_t)pos[i]);
                        else
                            r.back() = (uint8_t)pos[i];
                    return r;
                };
                const std::vector<uint8_t> present_runs = runs(present), absent_runs = runs(absent);
                uint64_t form = BITMAP, best = bb;
                if (present.size() < best)
                    form = LIST, best = present.size();
                if (absent.size() < best)
                    form = ABSENT, best = absent.size();
                if (present_runs.size() < best)
                    form = RUNS, best = present_runs.size();
                if (absent_runs.size() < best)
                    form = ABSENT_RUNS, best = absent_runs.size();
                if (form == LIST || form == ABSENT)
                    for (uint32_t off : form == LIST ? present : absent)
                        data.push_back((uint8_t)off);
                else if (form == RUNS || form == ABSENT_RUNS)
                {
                    const std::vector<uint8_t> &r = form == RUNS ? present_runs : absent_runs;
                    data.insert(data.end(), r.begin(), r.end());
                }
                else
                    for (uint64_t j = 0; j < bb; j++)
                        data.push_back((uint8_t)(k[j / 8] >> (8 * (j % 8))));
                forms[form]++;
                rec[4 + g % PART_GROUP] = (uint8_t)((data.size() - before) | ((form == BITMAP ? LIST : form) << 6));
            }
        }
        for (int i = 0; i < 8; i++)
            data.push_back(0); // padding for 8-byte loads
        precords.resize(records.size());
        if (!records.empty())
            std::memcpy(&precords[0], records.data(), records.size());
        pdata.resize(data.size());
        std::memcpy(&pdata[0], data.data(), data.size());
        std::cerr << "Meta-colors: " << n_blocks << " blocks of " << *std::min_element(sizes.begin(), sizes.end()) << "-"
                  << *std::max_element(sizes.begin(), sizes.end()) << " colors, " << n_partials << " distinct partial sets ("
                  << (data.size() + records.size()) / 1e6 << " MB: " << forms[LIST] << " lists, " << forms[ABSENT] << " lists of absent colors, " << forms[RUNS]
                  << " as runs, " << forms[ABSENT_RUNS] << " as runs of absent colors, " << forms[BITMAP] << " bitmaps), " << (n_sets ? (double)row_bits / (double)n_sets : 0.0) << " bits per color set ("
                  << (row_bytes + row_base.size_in_bytes() + row_rel.size_in_bytes()) / 1e6 << " MB with offsets)" << std::endl;
    }

    uint64_t size_in_bytes() const
    {
        return 5 * 8 + order.size() * 2 + (block_start.size() + pbase.size() + defaults.size() + empty_id.size() + code_count.size()) * 8 +
               precords.size() + pdata.size() + row_base.size_in_bytes() + row_rel.size_in_bytes() + rows.size() + 7 * 8;
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
        for (const auto *v : {&block_start, &pbase, &defaults, &empty_id, &code_count})
            for (uint64_t x : *v)
                u64(x);
        u64(precords.size());
        out.write(reinterpret_cast<const char *>(precords.data()), precords.size());
        u64(pdata.size());
        out.write(reinterpret_cast<const char *>(pdata.data()), pdata.size());
        row_base.serialize(out);
        row_rel.serialize(out);
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
        auto fail = [](const char *what)
        { throw std::runtime_error(std::string("The index is inconsistent (meta-color ") + what + ")"); };
        n_colors = u64();
        n_sets = u64();
        n_blocks = u64();
        row_bits = u64();
        n_partials = u64();
        if (!in || n_colors == 0 || n_colors > 65535 || n_blocks == 0 || n_blocks > n_colors)
            fail("sets");
        order.resize(n_colors);
        in.read(reinterpret_cast<char *>(order.data()), n_colors * 2);
        const uint64_t L = MAX_CODE + 1;
        block_start.resize(n_blocks + 1);
        pbase.resize(n_blocks + 1);
        defaults.resize(n_blocks);
        empty_id.resize(n_blocks);
        code_count.resize(n_blocks * L);
        for (auto *v : {&block_start, &pbase, &defaults, &empty_id, &code_count})
            for (auto &x : *v)
                x = u64();
        if (!in)
            fail("layout");
        for (uint16_t c : order)
            if (c >= n_colors)
                fail("order");
        if (block_start[0] != 0 || block_start[n_blocks] != n_colors || pbase[0] != 0 || pbase[n_blocks] != n_partials)
            fail("layout");
        for (uint64_t b = 0; b < n_blocks; b++)
        {
            if (block_start[b + 1] <= block_start[b] || block_size(b) > MAX_BLOCK || pbase[b + 1] <= pbase[b])
                fail("blocks");
            // The default is the last id, every other id has a code, and the codes are prefix-free
            uint64_t codes = 0;
            double kraft = 0;
            if (code_count[b * L] != 0)
                fail("codes");
            for (uint64_t l = 1; l <= MAX_CODE; l++)
            {
                codes += code_count[b * L + l];
                kraft += (double)code_count[b * L + l] / (double)(uint64_t(1) << l);
            }
            if (defaults[b] != pbase[b + 1] - pbase[b] - 1 || codes != defaults[b] || kraft > 1.0 + 1e-9 ||
                (empty_id[b] != NONE && empty_id[b] > defaults[b]))
                fail("codes");
        }
        const uint64_t rec_size = u64();
        if (in && rec_size != (n_partials + PART_GROUP - 1) / PART_GROUP * RECORD_BYTES)
            fail("partial set records");
        precords.load(in, rec_size);
        const uint64_t pdata_size = u64();
        pdata.load(in, pdata_size);
        row_base.load(in);
        row_rel.load(in);
        const uint64_t rows_size = u64();
        if (in && (rows_size != (row_bits + 7) / 8 + 8 || row_rel.size() != n_sets ||
                   row_base.size() != (n_sets + ROW_GROUP - 1) / ROW_GROUP))
            fail("rows");
        rows.load(in, rows_size);
        if (!in)
            return; // the caller reports the short file
        derive();
    }
};

// Counts the colors of a read's finimizers with meta-color sets. Block by block, the read's
// color sets are mapped to their partial sets, and each distinct partial set is added once,
// with the number of finimizers that have it as its weight. Rows and partial sets are fetched
// in stages, with prefetches, so that their memory latencies overlap. Lists add to the counts
// one color at a time; runs and bitmaps are added eight positions at a time through 8-bit lanes (as in
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
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> acc; // per block: (partial id, weight)
    std::vector<uint64_t> diff_weight; // per block: the weight of the partial sets in acc
    std::vector<uint64_t> positions;   // the row of each run
    struct Todo
    {
        uint64_t g, weight;
        uint32_t block;
        MetaColorSets::Partial part;
    };
    std::vector<Todo> todo;           // the distinct partial sets to add
    std::vector<uint64_t> everywhere; // per block: weight to add to every color (absent lists)

    static uint64_t clock_ns()
    {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
    }

public:
    // Time spent, in ns: finding the rows (with prefetching), reading them, and adding the
    // partial sets. The caller may read and reset these.
    uint64_t t_rows = 0, t_read = 0, t_add = 0;

private:
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

    // The bitmap (in words) of the runs of len bytes at ptr, or of the positions between them
    // if absent, in a block of size positions: each run toggles a bit at its start and one
    // past its end, and a prefix xor of the toggles fills the runs in, without branches
    static void runs_bitmap(const uint8_t *ptr, uint64_t len, bool absent, uint64_t size, uint64_t *bits)
    {
        uint64_t t[5] = {0, 0, 0, 0, 0};
        for (uint64_t j = 0; j < len; j += 2)
        {
            const uint64_t from = ptr[j], to = uint64_t(ptr[j + 1]) + 1;
            t[from >> 6] ^= uint64_t(1) << (from & 63);
            t[to >> 6] ^= uint64_t(1) << (to & 63);
        }
        uint64_t carry = absent ? ~uint64_t(0) : 0;
        for (uint64_t w = 0; w < 4; w++)
        {
            uint64_t x = t[w];
            x ^= x << 1;
            x ^= x << 2;
            x ^= x << 4;
            x ^= x << 8;
            x ^= x << 16;
            x ^= x << 32;
            x ^= carry;
            carry = uint64_t(0) - (x >> 63);
            const uint64_t keep = size >= 64 * (w + 1) ? ~uint64_t(0) : size <= 64 * w ? 0 : (uint64_t(1) << (size - 64 * w)) - 1;
            bits[w] = x & keep;
        }
    }

    void add_partial(const MetaColorSets &M, uint64_t b, const MetaColorSets::Partial &part, uint64_t weight)
    {
        if (part.form == MetaColorSets::RUNS || part.form == MetaColorSets::ABSENT_RUNS)
        {
            // Runs that cover few positions are added one color at a time, like lists
            const uint64_t bb = (M.block_size(b) + 7) / 8;
            uint64_t covered = 0;
            for (uint64_t j = 0; j < part.len; j += 2)
                covered += uint64_t(part.ptr[j + 1]) - part.ptr[j] + 1;
            if (covered <= bb)
            {
                const uint16_t *ord = M.order.data() + M.block_start[b];
                int32_t *c = counts.data();
                const bool absent = part.form == MetaColorSets::ABSENT_RUNS;
                const int32_t w = (int32_t)(absent ? -weight : weight);
                for (uint64_t j = 0; j < part.len; j += 2)
                    for (uint64_t p = part.ptr[j]; p <= part.ptr[j + 1]; p++)
                        c[ord[p]] += w;
                if (absent)
                    everywhere[b] += weight;
                return;
            }
        }
        if (part.form == MetaColorSets::LIST || part.form == MetaColorSets::ABSENT)
        {
            const uint16_t *ord = M.order.data() + M.block_start[b];
            int32_t *c = counts.data();
            const int32_t w = (int32_t)(part.form == MetaColorSets::LIST ? weight : -weight);
            for (uint64_t j = 0; j < part.len; j++)
                c[ord[part.ptr[j]]] += w;
            if (part.form == MetaColorSets::ABSENT)
                everywhere[b] += weight;
            return;
        }
        // Runs and bitmaps go to the lanes, runs as the bitmap they make
        auto add = [&](uint64_t wt)
        {
            if (part.form == MetaColorSets::BITMAP)
                add_bitmap_bytes(M, b, part.ptr, part.len, wt);
            else
            {
                uint64_t bits[4];
                const uint64_t size = M.block_size(b);
                runs_bitmap(part.ptr, part.len, part.form == MetaColorSets::ABSENT_RUNS, size, bits);
                add_bitmap_bytes(M, b, reinterpret_cast<const uint8_t *>(bits), (size + 7) / 8, wt);
            }
        };
        if (weight > 255)
        {
            flush(M, b, 1);
            add(1);
            pending[b] = 1;
            flush(M, b, (int32_t)weight);
            return;
        }
        if (pending[b] + weight > 255)
            flush(M, b, 1);
        add(weight);
        pending[b] += (uint32_t)weight;
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
        const uint64_t t0 = clock_ns();
        // Find the rows, prefetching what each step needs before it is read
        for (const auto &r : runs)
            M.prefetch_position((uint64_t)r.first);
        positions.resize(runs.size());
        for (size_t i = 0; i < runs.size(); i++)
        {
            positions[i] = M.row_position((uint64_t)runs[i].first);
            M.prefetch_row(positions[i]);
        }
        const uint64_t t1 = clock_ns();
        // The partial sets other than the defaults, block by block
        acc.resize(M.n_blocks);
        diff_weight.assign(M.n_blocks, 0);
        uint64_t total = 0;
        for (size_t i = 0; i < runs.size(); i++)
        {
            const uint32_t f = runs[i].second;
            total += f;
            M.for_each_difference_at(positions[i], [&](uint64_t b, uint64_t p)
                                     { acc[b].push_back({(uint32_t)p, f});
                                       diff_weight[b] += f; });
        }
        const uint64_t t2 = clock_ns();
        // The distinct partial sets with their weights, the defaults included
        todo.clear();
        everywhere.assign(M.n_blocks, 0);
        for (uint64_t b = 0; b < M.n_blocks; b++)
        {
            auto &v = acc[b];
            const uint64_t default_weight = total - diff_weight[b];
            if (default_weight && M.defaults[b] != M.empty_id[b])
                todo.push_back({M.pbase[b] + M.defaults[b], default_weight, (uint32_t)b, {}});
            if (v.empty())
                continue;
            used.clear();
            for (const auto &[p, f] : v)
            {
                if (p == M.empty_id[b])
                    continue;
                const uint32_t key = p + 1;
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
            for (uint32_t h : used)
            {
                todo.push_back({M.pbase[b] + table[h].first - 1, table[h].second, (uint32_t)b, {}});
                table[h] = {0, 0};
            }
            v.clear();
        }
        // Fetch the partial sets in stages, then add them
        for (const Todo &t : todo)
            __builtin_prefetch(M.partial_record(t.g));
        for (Todo &t : todo)
        {
            t.part = M.partial(t.g, t.block);
            __builtin_prefetch(t.part.ptr);
        }
        for (const Todo &t : todo)
            add_partial(M, t.block, t.part, t.weight);
        for (uint64_t b = 0; b < M.n_blocks; b++)
        {
            flush(M, b, 1);
            if (everywhere[b])
            {
                const uint16_t *ord = M.order.data() + M.block_start[b];
                for (uint64_t j = 0; j < M.block_size(b); j++)
                    counts[ord[j]] += (int32_t)everywhere[b];
            }
        }
        const uint64_t t3 = clock_ns();
        t_rows += t1 - t0;
        t_read += t2 - t1;
        t_add += t3 - t2;
        return counts.data();
    }
};
