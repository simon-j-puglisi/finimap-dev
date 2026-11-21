#pragma once

#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>

#include "mapped_array.hh"

// An array of unsigned integers of a fixed number of bits each (1 to 57, or 64), packed one after
// another from bit 0 of the first byte. Like Array, it either owns its bytes (while an index is
// built) or points into a mapped file.
//
// It also reads the plain arrays of index formats 1 and 2 (integers of 8, 16 or 64 bits, one
// after another), which are the same layout with a width of 8, 16 or 64.
class PackedInts
{
    Array<uint8_t> bytes;
    uint64_t n = 0;
    uint64_t width = 64;
    uint64_t mask = ~uint64_t(0);

    static uint64_t width_mask(const uint64_t w) { return w >= 64 ? ~uint64_t(0) : (uint64_t(1) << w) - 1; }

    // A value is read with one 8-byte load from the byte where it starts, so it must fit in the
    // 64 bits from there: widths up to 57, or 64 (where every value starts at a byte)
    void set_width(const uint64_t w)
    {
        if (w < 1 || (w > 57 && w != 64))
            throw std::runtime_error("The index is inconsistent (width of a packed array: " + std::to_string(w) + ")");
        width = w;
        mask = width_mask(w);
    }

public:
    // The bytes that n values of width bits take, plus 8 so that every value can be read with
    // one 8-byte load
    static uint64_t packed_bytes(const uint64_t n, const uint64_t w) { return (n * w + 7) / 8 + 8; }

    // The width for values up to max_value: the number of bits they need (at least 1), or 64
    // if that is more than 57
    static uint64_t bits_for(const uint64_t max_value)
    {
        uint64_t w = 1;
        while (w < 64 && (max_value >> w) != 0)
            w++;
        return w > 57 ? 64 : w;
    }

    size_t size() const { return n; }
    uint64_t bits() const { return width; }
    uint64_t size_in_bytes() const { return bytes.size(); }

    uint64_t operator[](const size_t i) const
    {
        const uint64_t pos = i * width;
        const uint64_t byte = pos >> 3;
        const uint8_t *base = reinterpret_cast<const uint8_t *>(bytes.data());
        uint64_t v;
        if (byte + 8 <= bytes.size())
            std::memcpy(&v, base + byte, 8);
        else
        {
            // Near the end of a plain array of format 1 or 2, which has no padding
            v = 0;
            std::memcpy(&v, base + byte, bytes.size() - byte);
        }
        return (v >> (pos & 7)) & mask;
    }

    // Packs values, each of which must fit in w bits
    template <typename T>
    void build(const T *values, const size_t count, const uint64_t w)
    {
        set_width(w);
        n = count;
        const uint64_t nb = packed_bytes(n, width);
        bytes.resize(nb);
        std::vector<uint8_t> buf(nb, 0);
        for (size_t i = 0; i < n; i++)
        {
            const uint64_t v = static_cast<uint64_t>(values[i]);
            if ((v & mask) != v)
                throw std::logic_error("PackedInts: a value does not fit in its width");
            uint64_t pos = i * width;
            uint64_t left = width, val = v;
            while (left > 0)
            {
                const uint64_t bit = pos & 7;
                const uint64_t take = std::min<uint64_t>(8 - bit, left);
                buf[pos >> 3] |= static_cast<uint8_t>((val & width_mask(take)) << bit);
                val >>= take;
                pos += take;
                left -= take;
            }
        }
        for (uint64_t i = 0; i < nb; i++)
            bytes[i] = buf[i];
    }

    // Format 3 and later: the number of values, the width, the number of bytes, the bytes
    void serialize(std::ostream &out) const
    {
        const uint64_t header[3] = {n, width, bytes.size()};
        out.write(reinterpret_cast<const char *>(header), sizeof(header));
        out.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    }

    void load(std::istream &in)
    {
        uint64_t header[3] = {0, 0, 0};
        in.read(reinterpret_cast<char *>(header), sizeof(header));
        if (!in)
            return;
        n = header[0];
        set_width(header[1]);
        if (n > (UINT64_MAX - 64) / width || header[2] != packed_bytes(n, width))
            throw std::runtime_error("The index is inconsistent (size of a packed array)");
        bytes.load(in, header[2]);
    }

    // Formats 1 and 2: count values of element_bytes bytes each (1, 2 or 8)
    void load_plain(std::istream &in, const uint64_t count, const uint64_t element_bytes)
    {
        if (count > UINT64_MAX / 8)
            throw std::length_error("plain array length");
        set_width(8 * element_bytes);
        n = count;
        bytes.load(in, count * element_bytes);
    }
};
