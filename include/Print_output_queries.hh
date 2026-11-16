#pragma once

#include <charconv>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// The decimal digits of every number below 10,000, so that most numbers of a result line
// are copied instead of converted.
struct SmallNumberText
{
    static constexpr uint32_t LIMIT = 10000;
    char digits[LIMIT][4];
    uint8_t length[LIMIT];

    SmallNumberText()
    {
        for (uint32_t x = 0; x < LIMIT; x++)
        {
            char buffer[8];
            const auto end = std::to_chars(buffer, buffer + sizeof(buffer), x).ptr;
            length[x] = static_cast<uint8_t>(end - buffer);
            for (uint8_t i = 0; i < 4; i++)
                digits[x][i] = (i < length[x]) ? buffer[i] : ' ';
        }
    }
};
inline const SmallNumberText SMALL_NUMBER_TEXT;

// Writes x in decimal at p and returns the position after it. Needs room for 20 characters.
inline char *write_number(char *p, const uint64_t x)
{
    if (x < SmallNumberText::LIMIT)
    {
        std::memcpy(p, SMALL_NUMBER_TEXT.digits[x], 4);
        return p + SMALL_NUMBER_TEXT.length[x];
    }
    return std::to_chars(p, p + 20, x).ptr;
}

// Formats result lines into a string.
//
// One line per read: the read's number (from 0, in input order), then "color:score" for each
// reported color, by decreasing score and then increasing color. Every field is followed by
// a space. With a denominator, it is written after the read's number.
class ResultFormatter
{
    std::vector<uint64_t> bucket_start;
    std::vector<std::pair<uint32_t, uint32_t>> ordered; // (color, score)

public:
    std::string out; // The formatted lines; the caller writes and clears it

    // Appends the line of one read. Colors with a score of at least min_score are reported;
    // min_score must be at least 1, and max_score at least the largest score.
    // If denominator is not negative, it is written after the read's number.
    void write(const uint64_t read_number, const int32_t *scores, const uint64_t n_colors, const int64_t min_score, const int64_t max_score, const int64_t denominator = -1)
    {
        uint64_t n_reported = 0;
        if (min_score <= max_score)
        {
            // Counting sort by decreasing score; colors stay in increasing order within a score
            bucket_start.assign(max_score - min_score + 2, 0);
            for (uint64_t color = 0; color < n_colors; color++)
            {
                const int32_t score = scores[color];
                if (score >= min_score)
                {
                    if (score > max_score)
                        throw std::runtime_error("A score exceeds the number of finimizers of its read");
                    bucket_start[max_score - score + 1]++;
                    n_reported++;
                }
            }
            if (n_reported > 0)
            {
                for (size_t i = 1; i < bucket_start.size(); i++)
                    bucket_start[i] += bucket_start[i - 1];
                ordered.resize(n_reported);
                for (uint32_t color = 0; color < n_colors; color++)
                {
                    const int32_t score = scores[color];
                    if (score >= min_score)
                        ordered[bucket_start[max_score - score]++] = {color, static_cast<uint32_t>(score)};
                }
            }
        }

        // Room for the longest possible line, then cut back to what was written
        const size_t old_size = out.size();
        out.resize(old_size + 48 + n_reported * 32);
        char *p = out.data() + old_size;
        p = write_number(p, read_number);
        *p++ = ' ';
        if (denominator >= 0)
        {
            p = write_number(p, denominator);
            *p++ = ' ';
        }
        for (uint64_t i = 0; i < n_reported; i++)
        {
            p = write_number(p, ordered[i].first);
            *p++ = ':';
            p = write_number(p, ordered[i].second);
            *p++ = ' ';
        }
        *p++ = '\n';
        out.resize(p - out.data());
    }
};
