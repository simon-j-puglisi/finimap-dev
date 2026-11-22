// Checks MetaColorSets and MetaColorCounter against direct counting, on random color sets of
// many shapes: sparse, dense, nearly full, clustered, with few and many colors, and with
// finimizer multisets whose weights pass 255. Also checks that the structure survives
// serialization.
// Build from the repository root:
//   g++ -std=c++2a -O2 -march=native -fopenmp -I include -I SDSL/include -o dev/checks/meta_check dev/checks/meta_check.cpp -L SDSL/lib -lsdsl
// Run: dev/checks/meta_check [ROUNDS] 2>/dev/null
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "MetaColorSets.hh"

static std::mt19937_64 rng(99);
static uint64_t rnd(uint64_t n) { return n ? rng() % n : 0; }

int main(int argc, char **argv)
{
    const int rounds = argc > 1 ? std::atoi(argv[1]) : 40;
    uint64_t failures = 0, kinds[3] = {0, 0, 0};
    for (int r = 0; r < rounds; r++)
    {
        const uint64_t C = 1 + rnd(r % 5 == 0 ? 2000 : 700);
        const uint64_t W = (C + 63) / 64;
        const uint64_t n_sets = 1 + rnd(3000);
        // A few "genome groups": sets take whole groups, with noise
        std::vector<uint64_t> group(C);
        const uint64_t n_groups = 1 + rnd(12);
        for (auto &g : group)
            g = rnd(n_groups);
        std::vector<uint64_t> words(n_sets * W, 0);
        for (uint64_t i = 0; i < n_sets; i++)
        {
            const int shape = (int)rnd(5);
            const uint64_t take = rng();
            for (uint64_t c = 0; c < C; c++)
            {
                bool in;
                if (shape == 0)
                    in = rnd(100) < 3; // sparse
                else if (shape == 1)
                    in = rnd(100) < 97; // nearly full
                else if (shape == 2)
                    in = rnd(2); // random half
                else
                    in = ((take >> (group[c] % 64)) & 1) != (rnd(100) < 2); // groups, with noise
                if (in)
                    words[i * W + c / 64] |= uint64_t(1) << (c % 64);
            }
        }
        MetaColorSets built;
        built.build(n_sets, C, [&](uint64_t i)
                    { return (const uint64_t *)words.data() + i * W; });
        std::stringstream ss;
        built.serialize(ss);
        MetaColorSets M;
        M.load(ss);
        for (uint64_t g = 0; g < M.n_partials; g++)
            kinds[M.pkind[g]]++;
        MetaColorCounter counter;
        for (int q = 0; q < 30; q++)
        {
            // A multiset of ids: a few sets with large weights, many with small
            std::vector<int64_t> ids;
            const uint64_t distinct = 1 + rnd(q % 3 == 0 ? 3 : 200);
            for (uint64_t d = 0; d < distinct; d++)
            {
                const int64_t id = (int64_t)rnd(n_sets);
                const uint64_t w = rnd(4) == 0 ? 1 + rnd(1000) : 1 + rnd(5);
                for (uint64_t k = 0; k < w; k++)
                    ids.push_back(id);
            }
            if (q == 0)
                ids.clear(); // no finimizers
            std::vector<int64_t> want(C, 0);
            for (int64_t id : ids)
                for (uint64_t c = 0; c < C; c++)
                    want[c] += (words[id * W + c / 64] >> (c % 64)) & 1;
            std::vector<int64_t> copy = ids;
            const int32_t *got = counter.count(copy.data(), copy.size(), M);
            for (uint64_t c = 0; c < C; c++)
                if (got[c] != want[c])
                {
                    if (failures++ < 10)
                        std::cout << "FAIL: round " << r << " (C = " << C << ", " << n_sets << " sets) query " << q << " color " << c << ": " << got[c] << ", expected " << want[c] << "\n";
                    break;
                }
        }
    }
    std::cout << "MetaColorSets counting against direct counting: " << (failures ? "FAILED" : "ok") << " (partial sets stored as lists " << kinds[0]
              << ", as lists of absent colors " << kinds[1] << ", as bitmaps " << kinds[2] << ")\n";
    return failures ? 1 : 0;
}
