// Checks the Pino structure and the short-finimizer lookups built on it.
//  1. Pino::predecessor and Pino::in_bucket against a sorted array, on random sets of many
//     densities (equal elements included).
//  2. PackedInts at every width, against the values packed.
//  3. Fluke8::finiLookup with Pino against Fluke8::finiLookup with Pred8v2 (index format 1),
//     on random prefix-free sets of strings, as the short finimizers are.
// Build from the repository root:
//   g++ -std=c++2a -O2 -march=native -I include -o dev/checks/pino_check dev/checks/pino_check.cpp
// Run: dev/checks/pino_check [ROUNDS]
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include "Fluke8.hh"

static std::mt19937_64 rng(12345);

static uint64_t rnd(uint64_t n) { return n ? rng() % n : 0; }

static int failures = 0;
static void fail(const std::string &what)
{
    if (failures++ < 20)
        std::cout << "FAIL: " << what << "\n";
}

// The reference: the first element equal to key if there is one, else the last element
// smaller than key, or -1
static int64_t ref_pred(const std::vector<uint64_t> &v, uint64_t key)
{
    auto it = std::lower_bound(v.begin(), v.end(), key);
    if (it != v.end() && *it == key)
        return it - v.begin();
    if (it == v.begin())
        return -1;
    return (it - v.begin()) - 1;
}

static void check_set(const std::vector<uint64_t> &v, uint64_t n_queries)
{
    Pino p(v);
    std::vector<uint64_t> queries;
    for (uint64_t x : v)
    {
        queries.push_back(x);
        if (x > 0)
            queries.push_back(x - 1);
        queries.push_back(x + 1);
        queries.push_back((x | 255));
        queries.push_back((x & ~uint64_t(255)));
    }
    const uint64_t top = v.back() + 600;
    for (uint64_t i = 0; i < n_queries; i++)
        queries.push_back(rnd(top));
    queries.push_back(0);
    queries.push_back(top + 100000);
    for (uint64_t q : queries)
    {
        uint64_t value = 0;
        bool known = false;
        int64_t got = p.predecessor(q, value, known);
        int64_t want = ref_pred(v, q);
        if (got != want)
        {
            fail("predecessor(" + std::to_string(q) + ") = " + std::to_string(got) + ", expected " + std::to_string(want) + " (n = " + std::to_string(v.size()) + ")");
            continue;
        }
        if (got < 0)
            continue;
        if (known && value != v[got])
            fail("predecessor value of " + std::to_string(q));
        if (!known && (v[got] >> 8) >= (q >> 8))
            fail("predecessor reported unknown in the same bucket, key " + std::to_string(q));
        if (p.low(got) != (v[got] & 255))
            fail("low byte");
        const uint64_t b = v[got] >> 8;
        if (!p.in_bucket(got, b))
            fail("in_bucket of its own bucket");
        if (b > 0 && p.in_bucket(got, b - 1))
            fail("in_bucket of the bucket before");
        if (b + 1 < p.buckets() && p.in_bucket(got, b + 1))
            fail("in_bucket of the bucket after");
    }
    for (uint64_t j = 0; j <= p.buckets(); j++)
    {
        uint64_t want = std::lower_bound(v.begin(), v.end(), j << 8) - v.begin();
        if (p.bucket_start(j) != want)
        {
            fail("bucket_start(" + std::to_string(j) + ")");
            break;
        }
    }
}

static void check_pino(int rounds)
{
    // Densities from sparse to dense, plus clustered sets and runs of equal values
    const double densities[] = {0.0005, 0.003, 0.01, 0.04, 0.2, 1.0, 3.0};
    for (int r = 0; r < rounds; r++)
    {
        for (double d : densities)
        {
            const uint64_t u = 256 * (1 + rnd(4000));
            const uint64_t n = std::max<uint64_t>(1, (uint64_t)(u * d / (d > 1 ? 4 : 1)));
            std::vector<uint64_t> v;
            const int kind = r % 3;
            for (uint64_t i = 0; i < n; i++)
            {
                if (kind == 0)
                    v.push_back(rnd(u));
                else if (kind == 1)
                    v.push_back((rnd(8) * u / 8) + rnd(u / 64 + 1)); // clusters with long gaps
                else
                    v.push_back(rnd(u / 3) * 3); // multiples, with equal values when dense
            }
            std::sort(v.begin(), v.end());
            if (r % 2)
                v.erase(std::unique(v.begin(), v.end()), v.end());
            check_set(v, 2000);
        }
    }
    // Edge cases: one element, all in one bucket, a full bucket of 256 and more
    check_set({0}, 100);
    check_set({255}, 100);
    check_set({256}, 100);
    check_set({1000000}, 100);
    {
        std::vector<uint64_t> v;
        for (uint64_t i = 0; i < 256; i++)
            v.push_back(512 + i);
        check_set(v, 1000);
        for (uint64_t i = 0; i < 300; i++)
            v.push_back(1024 + 3); // many equal values
        std::sort(v.begin(), v.end());
        check_set(v, 1000);
    }
    {
        std::vector<uint64_t> v;
        for (uint64_t i = 0; i < 100000; i++)
            v.push_back(i * 256 * 7); // one element every 7 buckets: long runs of empty buckets
        check_set(v, 20000);
    }
}

// A random prefix-free set of strings over ACGT with lengths in [lo, hi], sorted
static std::vector<std::string> prefix_free_set(uint64_t n, uint64_t lo, uint64_t hi, bool skewed)
{
    static const char ACGT[] = "ACGT";
    std::vector<std::string> all;
    for (uint64_t i = 0; i < n; i++)
    {
        uint64_t len = lo + rnd(hi - lo + 1);
        std::string s;
        for (uint64_t c = 0; c < len; c++)
            s += ACGT[skewed ? (rnd(10) < 7 ? 0 : 1 + rnd(3)) : rnd(4)];
        all.push_back(s);
    }
    std::sort(all.begin(), all.end());
    std::vector<std::string> out;
    for (const std::string &s : all)
    {
        if (!out.empty() && s.compare(0, out.back().size(), out.back()) == 0)
            continue; // out.back() is a prefix of s (or equal to it)
        out.push_back(s);
    }
    return out;
}

static uint64_t encode_key(const std::string &s)
{
    // 2 bits a character from the top bit, as search-fmin builds its keys
    uint64_t key = 0;
    for (size_t i = 0; i < s.size() && i < 32; i++)
        key |= (uint64_t)char2bits(s[i]) << (62 - 2 * i);
    return key;
}

static uint64_t n_cross = 0, n_found_all = 0;

static void check_fluke(int rounds)
{
    static const char ACGT[] = "ACGT";
    for (int r = 0; r < rounds; r++)
    {
        const uint64_t x = 4 + rnd(9); // 4..12, so that the old layout's table stays small
        const uint64_t lo = 1 + rnd(x);
        const uint64_t n = 1 + rnd(r % 4 == 0 ? 50 : 30000);
        std::vector<std::string> fins = prefix_free_set(n, lo, x, r % 3 == 1);
        std::vector<std::string_view> views(fins.begin(), fins.end());
        Fluke8 old_f(views, 0, x, true);
        Fluke8 new_f(views, 0, x, false);
        uint64_t n_found = 0;
        auto compare = [&](const std::string &q, uint64_t len)
        {
            const uint64_t key = encode_key(q);
            auto a = old_f.finiLookup(key, len);
            auto b = new_f.finiLookup(key, len);
            if (a != b)
                fail("x = " + std::to_string(x) + " query " + q + " len " + std::to_string(len) + ": Pred8v2 " + std::to_string(a.first) + "/" + std::to_string(a.second) + ", Pino " + std::to_string(b.first) + "/" + std::to_string(b.second));
            n_found += a.first >= 0;
            if (a.first >= 0 && (encode_key(fins[a.first]) >> (64 - 2 * x)) >> 8 != (key >> (64 - 2 * x)) >> 8)
                n_cross++; // found in an earlier bucket than the query's
        };
        for (const std::string &f : fins)
        {
            std::string q = f;
            while (q.size() < 32)
                q += ACGT[rnd(4)];
            compare(q, 31);
            compare(q, f.size());
            compare(q, f.size() - 1);
            std::string q2 = q;
            q2[rnd(f.size())] = ACGT[rnd(4)];
            compare(q2, 31);
        }
        for (uint64_t i = 0; i < 20000; i++)
        {
            std::string q;
            for (int c = 0; c < 32; c++)
                q += ACGT[r % 3 == 1 ? (rnd(10) < 7 ? 0 : 1 + rnd(3)) : rnd(4)];
            compare(q, 1 + rnd(31));
        }
        compare(std::string(32, 'A'), 31);
        compare(std::string(32, 'T'), 31);
        n_found_all += n_found;
        if (r == 0)
            std::cout << "  example round: x = " << x << ", " << fins.size() << " finimizers, " << n_found << " lookups found one\n";
    }
}

static void check_packed()
{
    for (uint64_t w = 1; w <= 64; w++)
    {
        if (w > 57 && w < 64)
            continue; // not a width PackedInts uses
        for (uint64_t n : {0ULL, 1ULL, 7ULL, 8ULL, 9ULL, 1000ULL})
        {
            std::vector<uint64_t> v(n);
            const uint64_t mask = w == 64 ? ~0ULL : (1ULL << w) - 1;
            for (auto &x : v)
                x = rng() & mask;
            if (n > 0)
                v[n - 1] = mask; // the largest value, at the end
            PackedInts p;
            p.build(v.data(), n, w);
            for (uint64_t i = 0; i < n; i++)
                if (p[i] != v[i])
                {
                    fail("PackedInts width " + std::to_string(w) + " index " + std::to_string(i));
                    break;
                }
        }
    }
    if (PackedInts::bits_for(0) != 1 || PackedInts::bits_for(1) != 1 || PackedInts::bits_for(2) != 2 || PackedInts::bits_for(4725954) != 23 || PackedInts::bits_for(~0ULL) != 64 || PackedInts::bits_for(1ULL << 57) != 64 || PackedInts::bits_for((1ULL << 57) - 1) != 57)
        fail("PackedInts::bits_for");
}

int main(int argc, char **argv)
{
    int rounds = argc > 1 ? std::atoi(argv[1]) : 30;
    // The constructors write messages to stderr: run with 2>/dev/null
    check_pino(rounds);
    std::cout << "Pino against a sorted array: " << (failures ? "FAILED" : "ok") << "\n";
    int before = failures;
    check_packed();
    std::cout << "PackedInts: " << (failures > before ? "FAILED" : "ok") << "\n";
    before = failures;
    check_fluke(rounds * 4);
    std::cout << "Lookups with Pino against Pred8v2: " << (failures > before ? "FAILED" : "ok") << " (" << n_found_all << " lookups found a finimizer, " << n_cross << " of them in an earlier bucket than the key's)\n";
    return failures ? 1 : 0;
}
