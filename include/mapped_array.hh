#pragma once

#include <cstdint>
#include <istream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// Reading an index without copying it: the file is mapped into memory, and the large arrays
// of the index point into the mapping instead of holding their own copy.

// A whole file mapped read-only. The pages are read in when the file is mapped, so the time
// to load an index is spent here and not at the first lookups.
class MappedFile
{
    const char *base = nullptr;
    size_t length = 0;

public:
    // Throws std::runtime_error if the file cannot be opened or mapped.
    MappedFile(const std::string &filename)
    {
        const int fd = open(filename.c_str(), O_RDONLY);
        if (fd < 0)
            throw std::runtime_error("Could not open " + filename);
        struct stat st;
        if (fstat(fd, &st) != 0 || st.st_size <= 0)
        {
            close(fd);
            throw std::runtime_error("Could not map " + filename);
        }
        void *p = mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
        close(fd);
        if (p == MAP_FAILED)
            throw std::runtime_error("Could not map " + filename);
        base = static_cast<const char *>(p);
        length = st.st_size;
    }
    ~MappedFile()
    {
        if (base)
            munmap(const_cast<char *>(base), length);
    }
    MappedFile(const MappedFile &) = delete;
    MappedFile &operator=(const MappedFile &) = delete;

    const char *data() const { return base; }
    size_t size() const { return length; }
};

// A stream buffer over a block of memory, so that the usual load(std::istream&) functions
// can read from a mapped file. Array::load recognises it and takes a pointer instead of a copy.
class MemoryStreamBuf : public std::streambuf
{
public:
    MemoryStreamBuf(const char *data, const size_t size)
    {
        char *p = const_cast<char *>(data);
        setg(p, p, p + size);
    }
    const char *current() const { return gptr(); }
    size_t remaining() const { return egptr() - gptr(); }
    // Moves forward by n bytes. Returns false, at the end of the data, if fewer are left.
    bool skip(const size_t n)
    {
        if (n > remaining())
        {
            setg(eback(), egptr(), egptr());
            return false;
        }
        setg(eback(), gptr() + n, egptr());
        return true;
    }
};

// An array that either owns its elements (while an index is built) or points into a mapped
// file (when one is loaded). Mapped elements may lie at any address, so they are read through
// a type without alignment requirements. Only an owned array can be changed.
template <typename T>
class Array
{
public:
    typedef T element __attribute__((aligned(1), may_alias));

private:
    std::vector<T> owned;
    const element *view = nullptr; // Set when the elements lie in a mapped file
    size_t n = 0;

public:
    size_t size() const { return n; }
    const element *data() const { return view ? view : reinterpret_cast<const element *>(owned.data()); }
    T operator[](const size_t i) const { return data()[i]; }

    void resize(const size_t count)
    {
        if (view)
            throw std::logic_error("A mapped array cannot be resized");
        owned.resize(count);
        n = count;
    }
    T &operator[](const size_t i)
    {
        if (view)
            throw std::logic_error("A mapped array cannot be changed");
        return owned[i];
    }

    // Reads count elements from the stream. From a MemoryStreamBuf nothing is copied.
    void load(std::istream &in, const size_t count)
    {
        if (MemoryStreamBuf *memory = dynamic_cast<MemoryStreamBuf *>(in.rdbuf()))
        {
            const char *start = memory->current();
            if (count > memory->remaining() / sizeof(T) || !memory->skip(count * sizeof(T)))
            {
                memory->skip(memory->remaining());
                in.setstate(std::ios::failbit | std::ios::eofbit);
                return;
            }
            owned.clear();
            view = reinterpret_cast<const element *>(start);
            n = count;
        }
        else
        {
            view = nullptr;
            owned.resize(count);
            n = count;
            in.read(reinterpret_cast<char *>(owned.data()), count * sizeof(T));
        }
    }
};
