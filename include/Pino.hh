#ifndef _PINO_H_
#define _PINO_H_

// Pino: a predecessor structure for a sorted set of integers, from "Faster Small-Alphabet
// Subset Rank Queries" (Section 4.2). It keeps the byte of low bits of each element, as Pred8
// does, but replaces Pred8's array of 32-bit bucket pointers with a much smaller structure.
//
// The universe is split into buckets of 256 values. Bucket j holds c_j elements, and start(j)
// is the number of elements in the buckets before j. Conceptually there is a bitvector
// H' = 1 0^c_0 1 0^c_1 ... 1 0^c_{nb-1} 1, in which the jth 1 is at position P_j = j + start(j).
// The positions P_j are stored in a Pred8 structure that also supports select:
//   - _L2 is divided into blocks of 36 bytes, one for every 32 consecutive ones of H'. The
//     first 4 bytes of block t hold P_{32t} >> 8, and the other 32 the low bytes of P_{32t},
//     ..., P_{32t+31}.
//   - _H2[b] is the number of ones of H' at positions below 256 b.
// So select1(j) = P_j reads one block of _L2 and scans _H2 forward from the bucket that the
// block's first 1 lies in. In a Finimap index H' is dense (there are about 2.5 short
// finimizers per bucket at -x 16), so the scan is short and _H2 is small.

#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>

#include "mapped_array.hh"

class Pino
{
    uint64_t _n = 0;   // number of elements
    uint64_t _min = 0; // smallest element
    uint64_t _max = 0; // largest element
    uint64_t _nb = 0;  // number of buckets of the elements: (_max >> 8) + 1
    uint64_t _nb2 = 0; // number of buckets of H': (P_{_nb} >> 8) + 1
    Array<uint8_t> _L;   // the low byte of each element
    Array<uint8_t> _L2;  // the blocks of 36 bytes described above, for the _nb + 1 ones of H'
    Array<uint32_t> _H2; // _nb2 + 1 entries; _H2[_nb2] = _nb + 1

    static constexpr uint64_t BLOCK_BYTES = 36;
    static constexpr uint64_t BLOCK_ONES = 32;

    static uint32_t read32(const uint8_t *p)
    {
        uint32_t v;
        std::memcpy(&v, p, sizeof(v));
        return v;
    }

    const uint8_t *L() const { return reinterpret_cast<const uint8_t *>(_L.data()); }
    const uint8_t *L2() const { return reinterpret_cast<const uint8_t *>(_L2.data()); }

    // P_j, the position of the jth 1 of H', for j <= _nb. b is set to P_j >> 8.
    uint64_t select1(const uint64_t j, uint64_t &b) const
    {
        const uint8_t *block = L2() + (j / BLOCK_ONES) * BLOCK_BYTES;
        b = read32(block);
        const uint32_t *H2 = reinterpret_cast<const uint32_t *>(_H2.data());
        while (H2[b + 1] <= j)
            b++;
        return (b << 8) | block[4 + j % BLOCK_ONES];
    }

    // P_{j+1}, given b = P_j >> 8: the scan continues from there.
    uint64_t select1_next(const uint64_t j1, uint64_t &b) const
    {
        const uint32_t *H2 = reinterpret_cast<const uint32_t *>(_H2.data());
        while (H2[b + 1] <= j1)
            b++;
        return (b << 8) | L2()[(j1 / BLOCK_ONES) * BLOCK_BYTES + 4 + j1 % BLOCK_ONES];
    }

public:
    Pino() = default;

    // data must be sorted (equal values are allowed)
    explicit Pino(const std::vector<uint64_t> &data)
    {
        _n = data.size();
        if (_n == 0)
            throw std::runtime_error("Pino: the set is empty");
        if (_n >= (uint64_t(1) << 40))
            throw std::runtime_error("Pino: too many elements");
        _min = data[0];
        _max = data[_n - 1];
        _nb = (_max >> 8) + 1;
        if (_nb + 1 >= (uint64_t(1) << 32))
            throw std::runtime_error("Pino: the universe is too large (at most 2^40 values)");

        _L.resize(_n);
        for (uint64_t i = 0; i < _n; i++)
        {
            if (i > 0 && data[i] < data[i - 1])
                throw std::runtime_error("Pino: the elements are not sorted");
            _L[i] = static_cast<uint8_t>(data[i] & 255);
        }

        // The positions of the ones of H': P_j = j + start(j), for j = 0..._nb
        const uint64_t n_ones = _nb + 1;
        const uint64_t last_pos = _nb + _n;
        _nb2 = (last_pos >> 8) + 1;
        _H2.resize(_nb2 + 1);
        for (uint64_t b = 0; b <= _nb2; b++)
            _H2[b] = 0;
        const uint64_t n_blocks = (n_ones + BLOCK_ONES - 1) / BLOCK_ONES;
        _L2.resize(n_blocks * BLOCK_BYTES);
        for (uint64_t i = 0; i < n_blocks * BLOCK_BYTES; i++)
            _L2[i] = 0;

        uint64_t i = 0; // index into data
        for (uint64_t j = 0; j < n_ones; j++)
        {
            while (i < _n && (data[i] >> 8) < j)
                i++;
            const uint64_t pos = j + i; // i = start(j)
            _H2[(pos >> 8) + 1]++;
            uint8_t *block = &_L2[(j / BLOCK_ONES) * BLOCK_BYTES];
            if (j % BLOCK_ONES == 0)
            {
                const uint32_t high = static_cast<uint32_t>(pos >> 8);
                std::memcpy(block, &high, sizeof(high));
            }
            block[4 + j % BLOCK_ONES] = static_cast<uint8_t>(pos & 255);
        }
        for (uint64_t b = 0; b < _nb2; b++)
            _H2[b + 1] += _H2[b];
    }

    uint64_t size() const { return _n; }
    uint64_t min() const { return _min; }
    uint64_t max() const { return _max; }

    // start(j): the number of elements in the buckets before bucket j, for j <= _nb
    uint64_t bucket_start(const uint64_t j) const
    {
        uint64_t b;
        return select1(j, b) - j;
    }

    // The predecessor of key: the index of the largest element <= key, or -1 if there is none.
    // If there are equal elements, this is the first one equal to key, or else the last one
    // smaller than key, as with Pred8v2.
    // value is set to that element when it is known without more work: always, unless the
    // element lies in an earlier bucket than key. known tells which.
    int64_t predecessor(const uint64_t key, uint64_t &value, bool &known) const
    {
        known = false;
        if (key < _min)
            return -1;
        const uint64_t j = key >> 8;
        if (j >= _nb)
        {
            // key is larger than every element
            value = _max;
            known = true;
            return static_cast<int64_t>(_n - 1);
        }
        uint64_t b;
        const uint64_t s = select1(j, b) - j;
        const uint64_t e = select1_next(j + 1, b) - (j + 1);
        const uint64_t k = key & 255;
        const uint8_t *low = L();
        uint64_t p = s;
        while (p < e && low[p] < k)
            p++;
        if (p < e && low[p] == k)
        {
            value = key;
            known = true;
            return static_cast<int64_t>(p);
        }
        if (p > s)
        {
            value = (j << 8) | low[p - 1];
            known = true;
            return static_cast<int64_t>(p - 1);
        }
        // Every element of bucket j is larger than key, or the bucket is empty: the predecessor
        // is the last element of an earlier bucket. s > 0, as _min <= key.
        return static_cast<int64_t>(s - 1);
    }

    // The low byte of element idx
    uint8_t low(const uint64_t idx) const { return L()[idx]; }

    // Whether element idx lies in bucket j (j < _nb)
    bool in_bucket(const uint64_t idx, const uint64_t j) const
    {
        uint64_t b;
        const uint64_t s = select1(j, b) - j;
        if (idx < s)
            return false;
        const uint64_t e = select1_next(j + 1, b) - (j + 1);
        return idx < e;
    }

    uint64_t buckets() const { return _nb; }

    uint64_t sizeInBytes() const
    {
        return 5 * sizeof(uint64_t) + 3 * sizeof(uint64_t) + _L2.size() + _H2.size() * sizeof(uint32_t) + _L.size();
    }

    void serialize(std::ostream &os) const
    {
        const uint64_t header[5] = {_n, _min, _max, _nb, _nb2};
        os.write(reinterpret_cast<const char *>(header), sizeof(header));
        const uint64_t l2_size = _L2.size(), h2_size = _H2.size(), l_size = _L.size();
        os.write(reinterpret_cast<const char *>(&l2_size), sizeof(l2_size));
        os.write(reinterpret_cast<const char *>(_L2.data()), l2_size);
        os.write(reinterpret_cast<const char *>(&h2_size), sizeof(h2_size));
        os.write(reinterpret_cast<const char *>(_H2.data()), h2_size * sizeof(uint32_t));
        os.write(reinterpret_cast<const char *>(&l_size), sizeof(l_size));
        os.write(reinterpret_cast<const char *>(_L.data()), l_size);
    }

    void load(std::istream &is)
    {
        uint64_t header[5];
        is.read(reinterpret_cast<char *>(header), sizeof(header));
        _n = header[0];
        _min = header[1];
        _max = header[2];
        _nb = header[3];
        _nb2 = header[4];
        uint64_t l2_size, h2_size, l_size;
        is.read(reinterpret_cast<char *>(&l2_size), sizeof(l2_size));
        if (!is || l2_size != (_nb + BLOCK_ONES) / BLOCK_ONES * BLOCK_BYTES)
            throw std::runtime_error("The index is inconsistent (short finimizers: blocks)");
        _L2.load(is, l2_size);
        is.read(reinterpret_cast<char *>(&h2_size), sizeof(h2_size));
        if (!is || h2_size != _nb2 + 1)
            throw std::runtime_error("The index is inconsistent (short finimizers: bucket counts)");
        _H2.load(is, h2_size);
        is.read(reinterpret_cast<char *>(&l_size), sizeof(l_size));
        if (!is || l_size != _n)
            throw std::runtime_error("The index is inconsistent (short finimizers: low bytes)");
        _L.load(is, l_size);
        if (!is)
            return; // the caller reports the short file
        const Array<uint32_t> &H2 = _H2; // read-only: the arrays may lie in a mapped file
        if (_n == 0 || _nb != (_max >> 8) + 1 || _nb2 != ((_nb + _n) >> 8) + 1 || H2[_nb2] != _nb + 1 || _min > _max)
            throw std::runtime_error("The index is inconsistent (short finimizers)");
    }
};

#endif
