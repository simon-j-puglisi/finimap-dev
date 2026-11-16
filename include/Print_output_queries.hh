#pragma once

#include <charconv>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Formats result lines into a string.
//
// One line per read: the read's number (from 0, in input order), then "color:score" for each
// reported color, by decreasing score and then increasing color. Every field is followed by
// a space.
class ResultFormatter
{
    std::vector<uint64_t> bucket_start;
    std::vector<std::pair<uint32_t, uint32_t>> ordered; // (color, score)

    void append(const uint64_t x)
    {
        char digits[24];
        const auto end = std::to_chars(digits, digits + sizeof(digits), x).ptr;
        out.append(digits, end - digits);
    }

public:
    std::string out; // The formatted lines; the caller writes and clears it

    // Appends the line of one read. Colors with a score of at least min_score are reported;
    // min_score must be at least 1, and max_score at least the largest score.
    void write(const uint64_t read_number, const int32_t *scores, const uint64_t n_colors, const int64_t min_score, const int64_t max_score)
    {
        append(read_number);
        out.push_back(' ');

        if (min_score <= max_score)
        {
            // Counting sort by decreasing score; colors stay in increasing order within a score
            bucket_start.assign(max_score - min_score + 2, 0);
            uint64_t n_reported = 0;
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
                for (const auto &[color, score] : ordered)
                {
                    append(color);
                    out.push_back(':');
                    append(score);
                    out.push_back(' ');
                }
            }
        }
        out.push_back('\n');
    }
};
