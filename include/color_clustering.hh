#pragma once

// Orders the colors so that colors that occur in similar color sets are close together, and
// splits them into blocks: the basis of the meta-color sets of index format 4.
//
// The colors are compared on a sample of the color sets: the column of a color tells which of
// the sampled sets contain it, and the distance of two colors is the Hamming distance of their
// columns. They are clustered by average linkage (nearest-neighbor chain algorithm), and their
// order is the left-to-right order of the leaves of the tree. The blocks are the largest
// subtrees with at most max_block leaves, with neighboring small ones merged while they fit.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

namespace colorclust
{
    inline uint64_t hamming(const uint64_t *a, const uint64_t *b, const uint64_t W)
    {
        uint64_t d = 0;
        for (uint64_t w = 0; w < W; w++)
            d += __builtin_popcountll(a[w] ^ b[w]);
        return d;
    }

    // The columns: for each color, a bitmap over the sampled sets. get_set(i) must return the
    // words of set i (n_colors bits). Returns C * SW words.
    template <typename GetSet>
    std::vector<uint64_t> columns(const uint64_t n_sets, const uint64_t n_colors, const uint64_t sample_sets, GetSet get_set, uint64_t &SW)
    {
        std::mt19937_64 rng(7);
        std::vector<uint64_t> pick;
        const uint64_t n_pick = std::min(n_sets, sample_sets);
        for (uint64_t i = 0; i < n_pick; i++)
            pick.push_back(n_sets <= sample_sets ? i : rng() % n_sets);
        SW = (pick.size() + 63) / 64;
        if (SW == 0)
            SW = 1;
        const uint64_t W = (n_colors + 63) / 64;
        std::vector<uint64_t> col(n_colors * SW, 0);
        for (uint64_t j = 0; j < pick.size(); j++)
        {
            const uint64_t *r = get_set(pick[j]);
            for (uint64_t w = 0; w < W; w++)
            {
                uint64_t word = r[w];
                if (w == W - 1 && n_colors % 64)
                    word &= (uint64_t(1) << (n_colors % 64)) - 1;
                for (; word; word &= word - 1)
                {
                    const uint64_t c = w * 64 + __builtin_ctzll(word);
                    col[c * SW + j / 64] |= uint64_t(1) << (j % 64);
                }
            }
        }
        return col;
    }

    struct Tree
    {
        uint64_t n_leaves = 0;
        std::vector<uint32_t> order;                    // the leaves, left to right
        std::vector<std::array<uint32_t, 3>> node;      // internal node v (C + i): {left, right, leaves}

        uint32_t root() const { return node.empty() ? 0 : (uint32_t)(n_leaves + node.size() - 1); }
        uint64_t leaves(uint32_t v) const { return v < n_leaves ? 1 : node[v - n_leaves][2]; }
    };

    // Average-linkage clustering of the C colors. Needs C * C floats of memory.
    inline Tree average_linkage(const std::vector<uint64_t> &col, const uint64_t C, const uint64_t SW)
    {
        Tree t;
        t.n_leaves = C;
        if (C == 0)
            return t;
        std::vector<float> D(C * C);
#pragma omp parallel for schedule(dynamic, 16)
        for (uint64_t a = 0; a < C; a++)
        {
            D[a * C + a] = 0;
            for (uint64_t b = a + 1; b < C; b++)
            {
                const float d = (float)hamming(&col[a * SW], &col[b * SW], SW);
                D[a * C + b] = d;
                D[b * C + a] = d;
            }
        }
        std::vector<uint8_t> active(C, 1);
        std::vector<uint32_t> size(C, 1), id(C);
        std::iota(id.begin(), id.end(), 0);
        std::vector<uint64_t> chain;
        uint64_t n_active = C;
        while (n_active > 1)
        {
            if (chain.empty())
                for (uint64_t a = 0; a < C; a++)
                    if (active[a])
                    {
                        chain.push_back(a);
                        break;
                    }
            const uint64_t a = chain.back();
            const uint64_t prev = chain.size() >= 2 ? chain[chain.size() - 2] : UINT64_MAX;
            float bd = INFINITY;
            uint64_t bb = UINT64_MAX;
            for (uint64_t b = 0; b < C; b++)
            {
                if (!active[b] || b == a)
                    continue;
                const float d = D[a * C + b];
                if (d < bd || (d == bd && b == prev))
                {
                    bd = d;
                    bb = b;
                }
            }
            if (bb == prev)
            {
                chain.pop_back();
                chain.pop_back();
                const uint64_t x = std::min(a, prev), y = std::max(a, prev);
                const float nx = (float)size[x], ny = (float)size[y];
                for (uint64_t k = 0; k < C; k++)
                {
                    if (!active[k] || k == x || k == y)
                        continue;
                    const float d = (nx * D[x * C + k] + ny * D[y * C + k]) / (nx + ny);
                    D[x * C + k] = d;
                    D[k * C + x] = d;
                }
                t.node.push_back({id[x], id[y], size[x] + size[y]});
                id[x] = (uint32_t)(C + t.node.size() - 1);
                size[x] += size[y];
                active[y] = 0;
                n_active--;
            }
            else
                chain.push_back(bb);
        }
        std::vector<uint32_t> stack = {t.root()};
        while (!stack.empty())
        {
            const uint32_t v = stack.back();
            stack.pop_back();
            if (v < C)
                t.order.push_back(v);
            else
            {
                stack.push_back(t.node[v - C][1]);
                stack.push_back(t.node[v - C][0]);
            }
        }
        return t;
    }

    // Block sizes, in leaf order: the largest subtrees of at most max_block leaves, neighbors
    // merged while they fit in max_block together
    inline std::vector<uint64_t> tree_blocks(const Tree &t, const uint64_t max_block)
    {
        std::vector<uint64_t> sizes;
        if (t.n_leaves == 0)
            return sizes;
        std::vector<uint32_t> stack = {t.root()};
        while (!stack.empty())
        {
            const uint32_t v = stack.back();
            stack.pop_back();
            const uint64_t sz = t.leaves(v);
            if (sz <= max_block)
                sizes.push_back(sz);
            else
            {
                stack.push_back(t.node[v - t.n_leaves][1]);
                stack.push_back(t.node[v - t.n_leaves][0]);
            }
        }
        std::vector<uint64_t> merged;
        for (uint64_t sz : sizes)
        {
            if (!merged.empty() && merged.back() + sz <= max_block)
                merged.back() += sz;
            else
                merged.push_back(sz);
        }
        return merged;
    }

    inline std::vector<uint64_t> fixed_blocks(const uint64_t C, const uint64_t block)
    {
        std::vector<uint64_t> sizes;
        for (uint64_t p = 0; p < C; p += block)
            sizes.push_back(std::min(block, C - p));
        return sizes;
    }
}
