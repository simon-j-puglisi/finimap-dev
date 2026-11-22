// Estimates how much space the color sets of a Finimap index would take under cross-set
// encodings, without changing the index. It reads the color sets of an index and reports:
//
//   1. the current size of the color sets in the index;
//   2. per-set coding: each set alone in the best of a bitmap, an Elias-Fano list of its
//      colors, or one of the colors it lacks, and the entropy bound log2(C choose |S|);
//   3. reference coding: sets grouped around a representative (one level, no chains), each
//      stored as the list of colors where it differs from its representative;
//   4. meta-colors: the colors split into blocks, the distinct partial sets of each block
//      stored once, and each set stored as the ids of its partial sets.
//
// Every estimate uses the same costs, so that they can be compared with each other: a list of
// m values in a universe of size u takes m * (2 + floor(log2(u / m))) bits (Elias-Fano), a
// bitmap u bits. Offsets of variable-length sets, which every scheme needs, are left out.
// Optionally the sets are also written to a file, to compress with a general-purpose
// compressor.
//
// Build (from the repository root):
//   g++ -std=c++2a -O3 -march=native -fopenmp -I include -I SDSL/include -o dev/colorsets/cs_analyze \
//       dev/colorsets/cs_analyze.cpp src/utilities.cpp src/seqio.cpp -L SDSL/lib -lsdsl -lz
// Run: cs_analyze INDEX_PREFIX [--threads N] [--dump FILE] [--quick]
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <omp.h>

#include "ColoredFinimizers.hh"
#include "Colors_queries.hh"

using std::vector;

static double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static uint64_t hmix(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Bits of an Elias-Fano list of m values below u
static double ef_bits(uint64_t m, uint64_t u)
{
    if (m == 0)
        return 0;
    if (m >= u)
        return (double)u; // then a bitmap is no larger
    const double low = std::floor(std::log2((double)u / (double)m));
    return (double)m * (2.0 + std::max(0.0, low));
}

// The cheapest of a bitmap, a list of the members and a list of the non-members
static double set_bits(uint64_t m, uint64_t u)
{
    return std::min({(double)u, ef_bits(m, u), ef_bits(u - m, u)});
}

static double log2_binomial(uint64_t n, uint64_t k)
{
    return (std::lgamma((double)n + 1) - std::lgamma((double)k + 1) - std::lgamma((double)(n - k) + 1)) / std::log(2.0);
}

static uint64_t bits_for(uint64_t n) // bits to tell apart n values
{
    uint64_t b = 0;
    while (b < 64 && (uint64_t(1) << b) < n)
        b++;
    return b;
}

static std::string mb(double bits)
{
    char buf[64];
    snprintf(buf, sizeof(buf), "%.1f MB", bits / 8e6);
    return buf;
}

struct Sets
{
    uint64_t C = 0;      // number of colors
    uint64_t W = 0;      // 64-bit words per set
    uint64_t n = 0;      // number of distinct sets
    vector<uint64_t> bits; // n * W words
    vector<uint32_t> size; // number of colors in each set
    vector<uint64_t> weight; // number of finimizers with each set
    vector<uint8_t> kind;  // 0: listed, 1: listed as complement, 2: bitmap (as in the index)

    const uint64_t *row(uint64_t i) const { return bits.data() + i * W; }
    uint64_t *row(uint64_t i) { return bits.data() + i * W; }
};

static uint64_t hamming(const uint64_t *a, const uint64_t *b, uint64_t W)
{
    uint64_t d = 0;
    for (uint64_t w = 0; w < W; w++)
        d += __builtin_popcountll(a[w] ^ b[w]);
    return d;
}

// Reads every distinct color set of the index into bitmaps
static Sets read_sets(const CompressedColoredFinimizers &index)
{
    Sets s;
    s.C = index.n_colors;
    s.W = (s.C + 63) / 64;
    const CompressedColorSets &CCS = index.CCS;
    const uint64_t sparse_count = CCS.sparse_count, dense_count = CCS.dense_count;
    uint64_t max_id = 0;
    const PackedInts &ids = index.loaded_color_set_ids;
    for (uint64_t i = 0; i < ids.size(); i++)
        max_id = std::max<uint64_t>(max_id, ids[i]);
    const uint64_t n_lists = sparse_count > 0 ? sparse_count - 1 : 0;
    const uint64_t n_complements = dense_count - sparse_count;
    const uint64_t n_bitmaps = max_id >= dense_count ? max_id - dense_count + 1 : 0;
    // id 0 is not a set: sparse ids start at 1, then complements, then bitmaps
    s.n = n_lists + n_complements + n_bitmaps;
    std::cerr << "Sets: " << n_lists << " listed, " << n_complements << " listed as complements, " << n_bitmaps << " bitmaps, " << s.C << " colors" << std::endl;
    s.bits.assign(s.n * s.W, 0);
    s.size.assign(s.n, 0);
    s.weight.assign(s.n, 0);
    s.kind.assign(s.n, 0);
    const DeltaSet &EF = CCS.getEF();
    const PackedInts &L = CCS.list_colors();
    const Array<uint64_t>::element *data = CCS.bitmap_words();
    auto index_of = [&](uint64_t id) -> int64_t
    {
        if (id >= 1 && id < sparse_count)
            return (int64_t)id - 1;
        if (id >= sparse_count && id < dense_count)
            return (int64_t)(n_lists + (id - sparse_count));
        if (id >= dense_count)
            return (int64_t)(n_lists + n_complements + (id - dense_count));
        return -1;
    };
#pragma omp parallel for schedule(dynamic, 4096)
    for (uint64_t id = 1; id <= max_id; id++)
    {
        const int64_t i = index_of(id);
        if (i < 0 || (uint64_t)i >= s.n)
            continue;
        uint64_t *r = s.row(i);
        if (id < dense_count)
        {
            const uint64_t start = id > 0 ? EF.get_start(id - 1) : 0, end = EF.get_start(id);
            for (uint64_t p = start; p < end; p++)
            {
                const uint64_t c = L[p];
                r[c / 64] |= uint64_t(1) << (c % 64);
            }
            if (id >= sparse_count)
            {
                // the list holds the colors the set lacks
                for (uint64_t w = 0; w < s.W; w++)
                    r[w] = ~r[w];
                if (s.C % 64)
                    r[s.W - 1] &= (uint64_t(1) << (s.C % 64)) - 1;
                s.kind[i] = 1;
            }
        }
        else
        {
            const uint64_t pos = (id - dense_count) * s.C;
            for (uint64_t w = 0; w < s.W; w++)
            {
                const uint64_t bit = pos + 64 * w;
                const uint64_t n_bits = std::min<uint64_t>(64, s.C - 64 * w);
                uint64_t word = (uint64_t)data[bit / 64] >> (bit % 64);
                if (bit % 64 && (bit % 64) + n_bits > 64)
                    word |= (uint64_t)data[bit / 64 + 1] << (64 - bit % 64);
                if (n_bits < 64)
                    word &= (uint64_t(1) << n_bits) - 1;
                r[w] = word;
            }
            s.kind[i] = 2;
        }
        uint64_t m = 0;
        for (uint64_t w = 0; w < s.W; w++)
            m += __builtin_popcountll(r[w]);
        s.size[i] = (uint32_t)m;
    }
    for (uint64_t f = 0; f < ids.size(); f++)
    {
        const int64_t i = index_of(ids[f]);
        if (i >= 0 && (uint64_t)i < s.n)
            s.weight[i]++;
    }
    return s;
}

// ---- reference coding ----------------------------------------------------------------------

struct RefResult
{
    uint64_t n_reps = 0;
    double rep_bits = 0, diff_bits = 0, id_bits = 0;
    uint64_t diff_total = 0, members = 0;
};

// Locality-sensitive keys for finding similar sets. Sets with few or many colors get MinHash
// keys of the colors they have (or lack); the others get keys of sampled bit positions, as
// their Hamming distance is what the difference lists cost.
struct LshKeys
{
    uint64_t bands, rows;
    vector<uint32_t> perm;           // a random order of the colors (MinHash)
    vector<vector<uint32_t>> sample; // positions sampled for each band (bit sampling)

    LshKeys(uint64_t C, uint64_t bands, uint64_t rows, uint64_t seed) : bands(bands), rows(rows)
    {
        std::mt19937_64 rng(seed);
        perm.resize(C);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        sample.resize(bands);
        for (auto &v : sample)
            for (uint64_t r = 0; r < std::min<uint64_t>(3 * rows, 64); r++)
                v.push_back((uint32_t)(rng() % C));
    }

    void keys(const Sets &s, uint64_t i, uint64_t *out) const
    {
        const uint64_t *r = s.row(i);
        const uint64_t m = s.size[i];
        const bool lacking = m > s.C / 2;
        const bool listy = m < s.C / 16 || m > s.C - s.C / 16;
        if (listy)
        {
            // One-permutation MinHash of the colors it has (or lacks): the colors are put in a
            // random order and split into bands * rows bins, and each bin keeps its first color
            const uint64_t K = bands * rows;
            vector<uint32_t> mins(K, UINT32_MAX);
            for (uint64_t w = 0; w < s.W; w++)
            {
                uint64_t word = lacking ? ~r[w] : r[w];
                if (w == s.W - 1 && s.C % 64)
                    word &= (uint64_t(1) << (s.C % 64)) - 1;
                for (; word; word &= word - 1)
                {
                    const uint32_t rank = perm[w * 64 + __builtin_ctzll(word)];
                    uint32_t &v = mins[rank % K];
                    v = std::min(v, rank);
                }
            }
            for (uint64_t b = 0; b < bands; b++)
            {
                uint64_t h = hmix(1000 + b + (lacking ? 77 : 0));
                for (uint64_t j = b * rows; j < (b + 1) * rows; j++)
                    h = hmix(h ^ mins[j]);
                out[b] = h;
            }
        }
        else
        {
            // Bit sampling: the values of 3 * rows sampled colors
            for (uint64_t b = 0; b < bands; b++)
            {
                const auto &v = sample[b];
                uint64_t key = 0;
                for (uint64_t j = 0; j < v.size(); j++)
                    key |= ((r[v[j] / 64] >> (v[j] % 64)) & 1) << j;
                out[b] = hmix(key ^ hmix(5000 + b));
            }
        }
    }
};

// Groups the sets around representatives. The sets are taken in decreasing order of weight;
// each joins the existing representative that minimizes its cost, if that is cheaper than
// storing it alone, and becomes a representative otherwise. Representatives are existing sets
// (leaders); with majority, each group's representative is then replaced by the colorwise
// majority of its members where that is cheaper.
static RefResult reference_coding(const Sets &s, uint64_t bands, uint64_t rows, uint64_t max_candidates, bool majority, vector<int64_t> *assignment_out = nullptr)
{
    const uint64_t n = s.n, W = s.W, C = s.C;
    LshKeys lsh(C, bands, rows, 42);
    vector<uint64_t> keys(n * bands);
#pragma omp parallel for schedule(dynamic, 4096)
    for (uint64_t i = 0; i < n; i++)
        lsh.keys(s, i, &keys[i * bands]);

    vector<uint64_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](uint64_t a, uint64_t b)
                     { return s.weight[a] > s.weight[b]; });

    vector<std::unordered_map<uint64_t, vector<uint32_t>>> buckets(bands);
    vector<int64_t> leader_of(n, -1); // -1: is a representative
    vector<uint32_t> reps;
    // A rough id width for the decisions; the final cost uses the real number of
    // representatives
    const double id_guess = (double)bits_for(n / 4 + 1);

    const uint64_t CHUNK = 8192;
    vector<int64_t> best(CHUNK);
    for (uint64_t c0 = 0; c0 < n; c0 += CHUNK)
    {
        const uint64_t c1 = std::min(n, c0 + CHUNK);
#pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t t = c0; t < c1; t++)
        {
            const uint64_t i = order[t];
            int64_t bi = -1;
            uint64_t bd = UINT64_MAX;
            uint64_t seen = 0;
            vector<uint32_t> tried;
            for (uint64_t b = 0; b < bands && seen < max_candidates; b++)
            {
                auto it = buckets[b].find(keys[i * bands + b]);
                if (it == buckets[b].end())
                    continue;
                const auto &v = it->second;
                // the most recent leaders first; they are the lighter ones, but any is fine
                for (uint64_t j = 0; j < v.size() && seen < max_candidates; j++)
                {
                    const uint32_t cand = v[v.size() - 1 - j];
                    if (std::find(tried.begin(), tried.end(), cand) != tried.end())
                        continue;
                    tried.push_back(cand);
                    seen++;
                    const uint64_t d = hamming(s.row(i), s.row(cand), W);
                    if (d < bd)
                    {
                        bd = d;
                        bi = cand;
                    }
                }
            }
            if (bi >= 0 && ef_bits(bd, C) + id_guess < set_bits(s.size[i], C))
                best[t - c0] = bi;
            else
                best[t - c0] = -1;
        }
        for (uint64_t t = c0; t < c1; t++)
        {
            const uint64_t i = order[t];
            if (best[t - c0] >= 0)
                leader_of[i] = best[t - c0];
            else
            {
                leader_of[i] = -1;
                reps.push_back((uint32_t)i);
                for (uint64_t b = 0; b < bands; b++)
                {
                    auto &v = buckets[b][keys[i * bands + b]];
                    if (v.size() < 4 * max_candidates)
                        v.push_back((uint32_t)i);
                }
            }
        }
    }

    // Groups
    vector<int64_t> rep_index(n, -1);
    for (uint64_t r = 0; r < reps.size(); r++)
        rep_index[reps[r]] = (int64_t)r;
    vector<vector<uint32_t>> members(reps.size());
    for (uint64_t i = 0; i < n; i++)
        if (leader_of[i] >= 0)
            members[rep_index[leader_of[i]]].push_back((uint32_t)i);

    RefResult res;
    res.n_reps = reps.size();
    const double id_bits = (double)bits_for(reps.size());
    std::atomic<uint64_t> diff_total{0}, n_members{0};
    double rep_bits = 0, diff_bits = 0;
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : rep_bits, diff_bits)
    for (uint64_t r = 0; r < reps.size(); r++)
    {
        const uint64_t leader = reps[r];
        const auto &mem = members[r];
        // Cost with the leader as representative
        double leader_rep = set_bits(s.size[leader], C), leader_diff = 0;
        uint64_t leader_dsum = 0;
        for (uint32_t i : mem)
        {
            const uint64_t d = hamming(s.row(i), s.row(leader), W);
            leader_dsum += d;
            leader_diff += ef_bits(d, C);
        }
        double best_rep = leader_rep, best_diff = leader_diff;
        uint64_t best_dsum = leader_dsum;
        if (majority && mem.size() >= 2)
        {
            // The majority of the group, the leader included; the leader is then a member too
            vector<uint32_t> count(C, 0);
            const uint64_t g = mem.size() + 1;
            auto add = [&](uint64_t i)
            {
                const uint64_t *row = s.row(i);
                for (uint64_t w = 0; w < W; w++)
                    for (uint64_t word = row[w]; word; word &= word - 1)
                        count[w * 64 + __builtin_ctzll(word)]++;
            };
            add(leader);
            for (uint32_t i : mem)
                add(i);
            vector<uint64_t> maj(W, 0);
            uint64_t msize = 0;
            for (uint64_t c = 0; c < C; c++)
                if (2 * count[c] > g)
                {
                    maj[c / 64] |= uint64_t(1) << (c % 64);
                    msize++;
                }
            double mrep = set_bits(msize, C), mdiff = 0;
            uint64_t mdsum = 0;
            auto one = [&](uint64_t i)
            {
                const uint64_t d = hamming(s.row(i), maj.data(), W);
                mdsum += d;
                mdiff += ef_bits(d, C);
            };
            one(leader);
            for (uint32_t i : mem)
                one(i);
            // the leader now needs its own id and difference
            if (mrep + mdiff + id_bits < best_rep + best_diff)
            {
                best_rep = mrep;
                best_diff = mdiff + id_bits;
                best_dsum = mdsum;
            }
        }
        rep_bits += best_rep;
        diff_bits += best_diff;
        diff_total += best_dsum;
        n_members += mem.size();
    }
    res.rep_bits = rep_bits;
    res.diff_bits = diff_bits;
    res.members = n_members;
    res.diff_total = diff_total;
    res.id_bits = id_bits * (double)n_members;
    if (assignment_out)
        *assignment_out = leader_of;
    return res;
}

// ---- meta-colors ---------------------------------------------------------------------------

// An order of the colors in which colors that occur in similar sets are close: a greedy
// nearest-neighbor chain over the colors, compared on a sample of the sets.
static vector<uint32_t> color_order(const Sets &s, uint64_t sample_sets)
{
    const uint64_t C = s.C, n = s.n;
    std::mt19937_64 rng(7);
    vector<uint64_t> pick;
    for (uint64_t i = 0; i < std::min(n, sample_sets); i++)
        pick.push_back(n <= sample_sets ? i : rng() % n);
    const uint64_t SW = (pick.size() + 63) / 64;
    vector<uint64_t> col(C * SW, 0); // column c: which sampled sets contain color c
    for (uint64_t j = 0; j < pick.size(); j++)
    {
        const uint64_t *r = s.row(pick[j]);
        for (uint64_t w = 0; w < s.W; w++)
            for (uint64_t word = r[w]; word; word &= word - 1)
            {
                const uint64_t c = w * 64 + __builtin_ctzll(word);
                col[c * SW + j / 64] |= uint64_t(1) << (j % 64);
            }
    }
    vector<uint32_t> order;
    vector<uint8_t> used(C, 0);
    uint64_t cur = 0;
    used[0] = 1;
    order.push_back(0);
    for (uint64_t step = 1; step < C; step++)
    {
        uint64_t best = UINT64_MAX, bc = 0;
#pragma omp parallel
        {
            uint64_t lb = UINT64_MAX, lc = 0;
#pragma omp for nowait
            for (uint64_t c = 0; c < C; c++)
            {
                if (used[c])
                    continue;
                const uint64_t d = hamming(&col[cur * SW], &col[c * SW], SW);
                if (d < lb || (d == lb && c < lc))
                {
                    lb = d;
                    lc = c;
                }
            }
#pragma omp critical
            if (lb < best || (lb == best && lc < bc))
            {
                best = lb;
                bc = lc;
            }
        }
        used[bc] = 1;
        order.push_back((uint32_t)bc);
        cur = bc;
    }
    return order;
}

struct MetaResult
{
    uint64_t block = 0, n_blocks = 0, distinct = 0;
    double partial_bits = 0, ref_bits = 0;
};

// Sizes with the colors, in the given order, split into blocks of `block` colors
static MetaResult meta_colors(const Sets &s, const vector<uint32_t> &order, uint64_t block)
{
    const uint64_t C = s.C, n = s.n;
    const uint64_t n_blocks = (C + block - 1) / block;
    const uint64_t BW = (block + 63) / 64; // words of a partial set
    vector<uint32_t> pos(C);
    for (uint64_t p = 0; p < C; p++)
        pos[order[p]] = (uint32_t)p;
    // Each set's partial sets, hashed, block by block
    vector<uint64_t> hashes(n * n_blocks);
    vector<uint16_t> psize(n * n_blocks);
#pragma omp parallel
    {
        vector<uint64_t> tmp(n_blocks * BW);
#pragma omp for schedule(dynamic, 4096)
        for (uint64_t i = 0; i < n; i++)
        {
            std::fill(tmp.begin(), tmp.end(), 0);
            const uint64_t *r = s.row(i);
            for (uint64_t w = 0; w < s.W; w++)
                for (uint64_t word = r[w]; word; word &= word - 1)
                {
                    const uint64_t p = pos[w * 64 + __builtin_ctzll(word)];
                    const uint64_t b = p / block, off = p % block;
                    tmp[b * BW + off / 64] |= uint64_t(1) << (off % 64);
                }
            for (uint64_t b = 0; b < n_blocks; b++)
            {
                uint64_t h = hmix(b + 1), m = 0;
                for (uint64_t w = 0; w < BW; w++)
                {
                    h = hmix(h ^ tmp[b * BW + w]);
                    m += __builtin_popcountll(tmp[b * BW + w]);
                }
                hashes[i * n_blocks + b] = h;
                psize[i * n_blocks + b] = (uint16_t)m;
            }
        }
    }
    MetaResult res;
    res.block = block;
    res.n_blocks = n_blocks;
    vector<uint64_t> distinct(n_blocks);
    double partial_bits = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : partial_bits)
    for (uint64_t b = 0; b < n_blocks; b++)
    {
        const uint64_t bsize = std::min(block, C - b * block);
        std::unordered_map<uint64_t, uint16_t> seen;
        seen.reserve(1024);
        for (uint64_t i = 0; i < n; i++)
            seen.emplace(hashes[i * n_blocks + b], psize[i * n_blocks + b]);
        distinct[b] = seen.size();
        for (auto &[h, m] : seen)
            partial_bits += set_bits(m, bsize);
    }
    res.partial_bits = partial_bits;
    for (uint64_t b = 0; b < n_blocks; b++)
        res.distinct += distinct[b];
    // Each set: the id of its partial set in every block, or the list of its nonempty blocks
    // with their ids, whichever is smaller
    vector<double> idw(n_blocks);
    double fixed = 0;
    for (uint64_t b = 0; b < n_blocks; b++)
    {
        idw[b] = (double)bits_for(distinct[b]);
        fixed += idw[b];
    }
    const double block_id = (double)bits_for(n_blocks);
    double ref_bits = 0;
#pragma omp parallel for schedule(static) reduction(+ : ref_bits)
    for (uint64_t i = 0; i < n; i++)
    {
        double listed = 0;
        for (uint64_t b = 0; b < n_blocks; b++)
            if (psize[i * n_blocks + b] != 0)
                listed += block_id + idw[b];
        ref_bits += std::min(fixed, listed) + 1; // one bit to tell which
    }
    res.ref_bits = ref_bits;
    return res;
}

// ---- main ----------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: cs_analyze INDEX_PREFIX [--threads N] [--dump FILE] [--quick]" << std::endl;
        return 1;
    }
    std::string prefix = argv[1], dump;
    int threads = omp_get_max_threads();
    bool quick = false;
    for (int a = 2; a < argc; a++)
    {
        std::string arg = argv[a];
        if (arg == "--threads" && a + 1 < argc)
            threads = std::atoi(argv[++a]);
        else if (arg == "--dump" && a + 1 < argc)
            dump = argv[++a];
        else if (arg == "--quick")
            quick = true;
        else
        {
            std::cerr << "Unknown argument " << arg << std::endl;
            return 1;
        }
    }
    omp_set_num_threads(threads);
    double t0 = now();
    CompressedColoredFinimizers index;
    index.load(prefix);
    const Sets s = read_sets(index);
    const uint64_t C = s.C, n = s.n;
    std::cerr << "Read " << n << " sets in " << now() - t0 << " s" << std::endl;

    // Check the sets read against the query's own counting, for a sample of finimizers
    {
        ColorCounter counter;
        const PackedInts &ids = index.loaded_color_set_ids;
        const uint64_t sparse_count = index.CCS.sparse_count, dense_count = index.CCS.dense_count;
        const uint64_t n_lists = sparse_count > 0 ? sparse_count - 1 : 0, n_complements = dense_count - sparse_count;
        std::mt19937_64 rng(3);
        uint64_t bad = 0, checked = 0;
        for (uint64_t t = 0; t < 2000 && ids.size() > 0; t++)
        {
            int64_t id = (int64_t)ids[rng() % ids.size()];
            const uint64_t uid = (uint64_t)id;
            const uint64_t i = uid < sparse_count ? uid - 1 : uid < dense_count ? n_lists + (uid - sparse_count) : n_lists + n_complements + (uid - dense_count);
            const int32_t *c = counter.count(&id, 1, index.CCS, C);
            for (uint64_t col = 0; col < C; col++)
                bad += (uint64_t)(c[col] != (int32_t)((s.row(i)[col / 64] >> (col % 64)) & 1));
            checked++;
        }
        printf("CS_CHECK %lu finimizers' sets against the query counter: %s\n", checked, bad ? "MISMATCH" : "ok");
        if (bad)
            return 2;
    }

    // 1. Current size, from the index's own arrays
    const CompressedColorSets &CCS = index.CCS;
    const double cur_lists = (double)CCS.list_colors().size_in_bytes() * 8;
    uint64_t n_bitmaps = 0, n_listed = 0, n_complements = 0, listed_colors = 0;
    for (uint64_t i = 0; i < n; i++)
    {
        if (s.kind[i] == 2)
            n_bitmaps++;
        else if (s.kind[i] == 1)
        {
            n_complements++;
            listed_colors += C - s.size[i];
        }
        else
        {
            n_listed++;
            listed_colors += s.size[i];
        }
    }
    const double cur_bitmaps = (double)n_bitmaps * C;
    uint64_t total_weight = 0;
    for (uint64_t w : s.weight)
        total_weight += w;

    printf("CS_INDEX %s: %lu colors, %lu distinct color sets (%lu listed, %lu as complements, %lu bitmaps), %lu finimizers\n",
           prefix.c_str(), C, n, n_listed, n_complements, n_bitmaps, total_weight);
    printf("CS_CURRENT bitmaps %s + lists %s (%lu colors) = %s\n", mb(cur_bitmaps).c_str(), mb(cur_lists).c_str(), listed_colors, mb(cur_bitmaps + cur_lists).c_str());

    // 2. Per-set coding
    double per_set = 0, entropy = 0, per_set_bitmaps = 0, per_set_lists = 0;
    for (uint64_t i = 0; i < n; i++)
    {
        const double b = set_bits(s.size[i], C);
        per_set += b;
        (s.kind[i] == 2 ? per_set_bitmaps : per_set_lists) += b;
        entropy += log2_binomial(C, s.size[i]);
    }
    printf("CS_PER_SET best of bitmap / EF list / EF complement: %s (the sets now in bitmaps %s, the listed ones %s); entropy bound %s\n",
           mb(per_set).c_str(), mb(per_set_bitmaps).c_str(), mb(per_set_lists).c_str(), mb(entropy).c_str());
    // How much the dense sets differ from each other: sizes of the sets now in bitmaps
    {
        vector<uint32_t> sz;
        for (uint64_t i = 0; i < n; i++)
            if (s.kind[i] == 2)
                sz.push_back(s.size[i]);
        std::sort(sz.begin(), sz.end());
        if (!sz.empty())
            printf("CS_BITMAP_SIZES colors in the bitmap sets: 10%% %u, median %u, 90%% %u (of %lu)\n", sz[sz.size() / 10], sz[sz.size() / 2], sz[sz.size() * 9 / 10], C);
    }
    fflush(stdout);

    // Dump the sets, as bitmaps in the index's order, for a general-purpose compressor
    if (!dump.empty())
    {
        FILE *f = fopen(dump.c_str(), "wb");
        if (!f)
        {
            std::cerr << "Cannot write " << dump << std::endl;
            return 1;
        }
        const uint64_t bytes = (C + 7) / 8;
        for (uint64_t i = 0; i < n; i++)
            fwrite(s.row(i), 1, bytes, f);
        fclose(f);
        printf("CS_DUMP %s: %lu sets of %lu bytes\n", dump.c_str(), n, bytes);
        fflush(stdout);
    }

    // 3. Reference coding
    struct Setting
    {
        uint64_t bands, rows, cand;
    };
    vector<Setting> settings = quick ? vector<Setting>{{8, 8, 64}} : vector<Setting>{{8, 8, 64}, {16, 6, 256}};
    for (const Setting &st : settings)
    {
        for (bool majority : {false, true})
        {
            const double t = now();
            const RefResult r = reference_coding(s, st.bands, st.rows, st.cand, majority);
            const double total = r.rep_bits + r.diff_bits + r.id_bits;
            printf("CS_REF bands=%lu rows=%lu candidates=%lu %s: %s (representatives %lu: %s; %lu members: differences %s, representative ids %s); mean difference %.1f colors; %.0f s\n",
                   st.bands, st.rows, st.cand, majority ? "majority" : "leaders", mb(total).c_str(), r.n_reps, mb(r.rep_bits).c_str(), r.members,
                   mb(r.diff_bits).c_str(), mb(r.id_bits).c_str(), r.members ? (double)r.diff_total / (double)r.members : 0.0, now() - t);
            fflush(stdout);
        }
    }

    // 4. Meta-colors
    {
        double t = now();
        vector<uint32_t> identity(C);
        std::iota(identity.begin(), identity.end(), 0);
        const vector<uint32_t> ordered = color_order(s, 32768);
        printf("CS_ORDER colors ordered by similarity in %.0f s\n", now() - t);
        vector<uint64_t> blocks = quick ? vector<uint64_t>{128} : vector<uint64_t>{64, 128, 256, 512};
        for (int o = 0; o < 2; o++)
            for (uint64_t b : blocks)
            {
                t = now();
                const MetaResult m = meta_colors(s, o ? ordered : identity, b);
                printf("CS_META order=%s block=%lu: %s (%lu blocks, %lu distinct partial sets: %s; partial set ids of each set %s); %.0f s\n",
                       o ? "similarity" : "input", b, mb(m.partial_bits + m.ref_bits).c_str(), m.n_blocks, m.distinct, mb(m.partial_bits).c_str(),
                       mb(m.ref_bits).c_str(), now() - t);
                fflush(stdout);
            }
    }
    printf("CS_DONE in %.0f s\n", now() - t0);
    return 0;
}
