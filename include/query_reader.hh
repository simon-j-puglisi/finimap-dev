#pragma once

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

// Reads a file line by line through zlib, which passes plain files through unchanged and
// decompresses gzipped ones.
class LineSource
{
    gzFile file = nullptr;
    std::string name;
    std::vector<char> buffer;
    size_t pos = 0, end = 0;

    // Refills the buffer. False at the end of the file.
    bool fill()
    {
        const int n = gzread(file, buffer.data(), buffer.size());
        if (n < 0)
        {
            int code = 0;
            const char *message = gzerror(file, &code);
            throw std::runtime_error("Could not read " + name + ": " + (message ? message : "read error"));
        }
        pos = 0;
        end = n;
        return n > 0;
    }

public:
    LineSource(const std::string &filename) : name(filename), buffer(1 << 18)
    {
        file = gzopen(filename.c_str(), "rb");
        if (!file)
            throw std::runtime_error("Could not open " + filename);
    }
    ~LineSource()
    {
        if (file)
            gzclose(file);
    }
    LineSource(const LineSource &) = delete;
    LineSource &operator=(const LineSource &) = delete;

    // Reads the next line into `line`, without its line feed. False at the end of the file.
    bool getline(std::string &line)
    {
        line.clear();
        bool any = false;
        while (true)
        {
            if (pos == end && !fill())
                return any;
            any = true;
            const char *start = buffer.data() + pos;
            const char *newline = static_cast<const char *>(memchr(start, '\n', end - pos));
            if (newline)
            {
                line.append(start, newline - start);
                pos += (newline - start) + 1;
                return true;
            }
            line.append(start, end - pos);
            pos = end;
        }
    }
};

// Reads query sequences from a FASTA or FASTQ file, possibly gzipped, one record at a time.
//
// The format is taken from the first character of the input ('>' or '@'), not from the
// file name. FASTA sequences may span several lines. FASTQ records must be four lines.
// Windows line endings and blank lines between records are accepted. An empty input holds
// no reads, and a record may have an empty sequence.
class QueryReader
{
    LineSource in;
    std::string name;      // For error messages
    std::string line;      // The current line, when have_line is set
    bool have_line = false;
    bool started = false;
    bool fastq = false;
    uint64_t line_number = 0;

    // Reads the next line into `line`, without its line ending. False at the end of input.
    bool next_line()
    {
        if (!in.getline(line))
            return false;
        line_number++;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return true;
    }

    [[noreturn]] void fail(const std::string &what) const
    {
        throw std::runtime_error(name + ", line " + std::to_string(line_number) + ": " + what);
    }

public:
    QueryReader(const std::string &filename) : in(filename), name(filename) {}

    // Reads the next record's sequence into seq. Returns false when there are no more records.
    bool next(std::string &seq)
    {
        seq.clear();
        if (!started)
        {
            started = true;
            // Skip blank lines before the first record
            while ((have_line = next_line()) && line.empty())
            {
            }
            if (!have_line)
                return false; // Empty input
            if (line[0] == '@')
                fastq = true;
            else if (line[0] != '>')
                fail("not a FASTA or FASTQ file (the first record starts with neither '>' nor '@')");
        }
        if (!have_line)
            return false;

        if (fastq)
        {
            if (line[0] != '@')
                fail("expected a FASTQ header starting with '@'");
            if (!next_line())
                fail("FASTQ record ends after its header");
            seq = line;
            if (!next_line() || line.empty() || line[0] != '+')
                fail("expected the '+' line of a FASTQ record (multi-line FASTQ is not supported)");
            if (!next_line())
                fail("FASTQ record has no quality line");
            if (line.size() != seq.size())
                fail("FASTQ quality line and sequence differ in length (multi-line FASTQ is not supported)");
            // Move to the next header, skipping blank lines
            while ((have_line = next_line()) && line.empty())
            {
            }
        }
        else
        {
            // `line` holds the header. The sequence is every line up to the next header.
            while ((have_line = next_line()) && (line.empty() || line[0] != '>'))
                seq += line;
        }
        return true;
    }
};
