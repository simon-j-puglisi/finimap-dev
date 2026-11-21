// Times predecessor queries on Pino against Pred8v2, on a random set with the density of the
// short finimizers of an E. coli index at -x 16 (about 41.7 million keys in 4^16 values).
// Build from the repository root:
//   g++ -std=c++2a -O2 -march=native -I include -o dev/checks/pino_bench dev/checks/pino_bench.cpp
// Run: dev/checks/pino_bench [N_KEYS [N_QUERIES]] 2>/dev/null
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "Pino.hh"
#include "Pred8v2.hh"

int main(int argc, char **argv)
{
    const uint64_t n = argc > 1 ? std::stoull(argv[1]) : 41700000;
    const uint64_t nq = argc > 2 ? std::stoull(argv[2]) : 20000000;
    const uint64_t u = uint64_t(1) << 32;
    std::mt19937_64 rng(7);
    std::vector<uint64_t> v(n);
    for (auto &x : v)
        x = rng() % u;
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    if (v[0] >= 256)
        v.insert(v.begin(), 0); // Pred8v2 expects a non-empty first bucket
    Pred8v2 old_p(v);
    Pino pino(v);
    std::cout << v.size() << " keys; Pred8v2 " << old_p.sizeInBytes() / 1e6 << " MB, Pino " << pino.sizeInBytes() / 1e6 << " MB\n";

    std::vector<uint64_t> q(nq);
    for (auto &x : q)
        x = rng() % u;
    for (int sorted = 0; sorted < 2; sorted++)
    {
        if (sorted)
        {
            // search-fmin groups its keys by their top 16 bits; sorting them is close to that
            std::sort(q.begin(), q.end());
        }
        for (int round = 0; round < 2; round++)
        {
            auto t0 = std::chrono::steady_clock::now();
            uint64_t sum = 0;
            for (uint64_t x : q)
            {
                auto r = old_p.getPredPrefix(x);
                sum += r.first + r.second;
            }
            auto t1 = std::chrono::steady_clock::now();
            uint64_t sum2 = 0;
            for (uint64_t x : q)
            {
                uint64_t value;
                bool known;
                sum2 += pino.predecessor(x, value, known) + known;
            }
            auto t2 = std::chrono::steady_clock::now();
            const double a = std::chrono::duration<double>(t1 - t0).count() * 1e9 / nq;
            const double b = std::chrono::duration<double>(t2 - t1).count() * 1e9 / nq;
            std::cout << (sorted ? "sorted" : "random") << " queries: Pred8v2 " << a << " ns, Pino " << b << " ns (" << sum % 10 << sum2 % 10 << ")\n";
        }
    }
}
