#pragma once

#include <memory>
#include <vector>
#include <optional>
#include <unordered_map>
#include <map>
#include <utility>

#include "sdsl/bit_vectors.hpp"
#include "CompressedColorSets.hh"

#include "Fluke8.hh"
#include "PrefTab.hh"
#include "mapped_array.hh"

#include "batch_querying.hh"

#include <chrono>

using namespace std;

static std::chrono::nanoseconds time_rarest_fmin(0);
static std::chrono::nanoseconds time_rarest_fmin_rc(0);
static std::chrono::nanoseconds time_combine(0);
static std::chrono::nanoseconds time_index_loading(0);
static std::chrono::nanoseconds time_output(0);

void print_search_timing_stats()
{
    using namespace std::chrono;
    std::cerr << "Time to load the index: "
              << duration_cast<milliseconds>(time_index_loading).count() << " ms\n";
    std::cerr << "Time in rarest_fmin_streaming_search (fwd): "
              << duration_cast<milliseconds>(time_rarest_fmin).count() << " ms\n";
    std::cerr << "Time in combine_f_rc: "
              << duration_cast<milliseconds>(time_combine).count() << " ms\n";
    std::cerr << "Time to output results: "
              << duration_cast<milliseconds>(time_output).count() << " ms\n";
}

// Colors are stored as 16-bit numbers in the compressed color sets
constexpr uint64_t MAX_COLORS = 65535;
// Finimizers are looked up through 62-bit keys, two bits per character
constexpr uint64_t MAX_K = 31;

// Colored finimizers without much compression
class ColoredFinimizers
{
public:
    vector<char> concat; // Finimizers concatenated in lexicographic order (ASCII characters).
    vector<uint8_t> lengths;
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors.
    vector<uint8_t> lengths_by_freq;

    struct VecConcat {
        std::vector<size_t> concat;
        std::vector<size_t> ends;
    };
    
    VecConcat sparse_colors;  // Sparse color representation (when meta=true)

    VecConcat load_sparse_colors(std::istream& file) {
        // uint8_t is_sparse;
        // file.read(reinterpret_cast<char*>(&is_sparse), 1);
        
        // if (is_sparse != 1) {
        //     throw std::runtime_error("Expected sparse format");
        // }
        
        uint64_t n_elements, n_sets;
        file.read(reinterpret_cast<char*>(&n_elements), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&n_sets), sizeof(uint64_t));
        cerr << n_sets << endl;
        
        std::vector<size_t> concat(n_elements);
        file.read(reinterpret_cast<char*>(concat.data()), n_elements * sizeof(size_t));
        
        std::vector<size_t> ends(n_sets);
        file.read(reinterpret_cast<char*>(ends.data()), n_sets * sizeof(size_t));
        
        return {concat, ends};
    }
    // Loads from the format output by the Rust CLI command `finimizer_matrix` with option --reverse.
    // That format has colexicographically sorted reverse finimizers. We reverse them to
    // get lex-sorted finimizers.
    void load(std::istream &in)
    {
        cerr << "Loading uncompressed finimizers..." << endl;
        uint64_t n_finimizers;
        in.read(reinterpret_cast<char *>(&n_finimizers), sizeof(n_finimizers));

        uint64_t finimizer_total_length;
        in.read(reinterpret_cast<char *>(&finimizer_total_length), sizeof(finimizer_total_length));

        uint64_t n_colors;
        in.read(reinterpret_cast<char *>(&n_colors), sizeof(n_colors));
        if (!in)
            throw std::runtime_error("Could not read the finimizer matrix file");
        if (n_colors > MAX_COLORS)
            throw std::runtime_error("The input has " + to_string(n_colors) + " colors, but at most " + to_string(MAX_COLORS) + " are supported");

        lengths.resize(n_finimizers);
        in.read(reinterpret_cast<char *>(lengths.data()), n_finimizers * sizeof(uint8_t));

        concat.resize(finimizer_total_length);
        in.read(reinterpret_cast<char *>(concat.data()), finimizer_total_length * sizeof(char));


        uint8_t is_sparse;
        in.read(reinterpret_cast<char *>(&is_sparse), 1);
        //cerr << (int)is_sparse << endl;
        if (is_sparse == 1) {
            cerr << "meta" << endl;
            sparse_colors = load_sparse_colors(in);
        }
        else{
            int64_t n_bits = n_finimizers * n_colors;
            // The bits are in u64 Lsb format, which is also the layout of sdsl::bit_vector,
            // so the words are read straight into it.
            int64_t n_words = (n_bits + 63) / 64; // Ceil div by 64
            color_sets_concat.resize(n_bits);
            in.read(reinterpret_cast<char *>(color_sets_concat.data()), n_words * sizeof(uint64_t));
        }
        

        uint64_t lengths_by_freq_size;
        in.read(reinterpret_cast<char *>(&lengths_by_freq_size), sizeof(lengths_by_freq_size));
        lengths_by_freq.resize(lengths_by_freq_size);
        in.read(reinterpret_cast<char *>(lengths_by_freq.data()), lengths_by_freq_size * sizeof(uint8_t));

        cerr << "Reversing finimizer strings..." << endl;
        int64_t start_in_concat = 0;
        for (uint64_t f_idx = 0; f_idx < lengths.size(); f_idx++)
        {
            int64_t s = start_in_concat;
            int64_t e = start_in_concat + lengths[f_idx];
            std::reverse(concat.begin() + s, concat.begin() + e);
            start_in_concat = e;
        }
        cerr << "Uncompressed index loaded" << endl;
    }

    vector<pair<int64_t, uint64_t>> fmin_stats()
    {
        cerr << "Printing finimizers stats:" << endl;
        std::sort(lengths.begin(), lengths.end());
        vector<pair<int64_t, uint64_t>> lengths_and_freq;
        lengths_and_freq.reserve(lengths_by_freq.size());

        for (size_t i = 0; i < lengths.size();)
        {
            size_t j = i + 1;
            while (j < lengths.size() && lengths[j] == lengths[i])
                ++j;
            lengths_and_freq.emplace_back(lengths[i], j - i);
            i = j;
        }

        return lengths_and_freq;
    }
};

void true_or_crash(bool b, const char *error_message)
{
    if (!b)
    {
        cerr << error_message << endl;
        exit(1);
    }
}

//inline int16_t pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results, const float t);

//inline void pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results);

class CompressedColoredFinimizers
{

public:
    CompressedColorSets CCS; // L, EF, BV

    uint64_t n_colors;
    uint64_t n_finimizers;
    uint64_t plen;
    uint64_t short_long_t;
    uint64_t k;
    bool meta;

    Fluke8 f8;
    PrefTab pt;

    vector<uint64_t> color_set_ids;       // The color set of each finimizer, while the index is built
    Array<uint64_t> loaded_color_set_ids; // The same, in a loaded index
    std::shared_ptr<MappedFile> mapping;  // The index file, if it is loaded by mapping it

    int get_k() const { return k; }

    CompressedColoredFinimizers() = default;

    CompressedColoredFinimizers(ColoredFinimizers &&cf, int64_t prefix_len, int64_t short_long_t, uint64_t kmer_size, bool meta)
    {
        cerr << "Let's compress it!" << endl;
        this->plen = prefix_len;
        cerr << "Prefix length: " << plen << endl;
        this->short_long_t = short_long_t;
        cerr << "Threshold : " << short_long_t << endl;
        this->k = kmer_size;
        cerr << "k-mer size: " << k << endl;
        this->meta = meta;

        n_finimizers = cf.lengths.size();
        cerr << "total finimizers: " << (int)n_finimizers << endl;
        true_or_crash(n_finimizers > 0, "ERROR: 0 finimizers");

        if (meta) {
            // Sparse deduplication
            cerr << "Deduplicate sparse color sets..." << endl;
            n_colors = cf.sparse_colors.ends.size();
            map<vector<size_t>, vector<size_t>> color_set_to_finimizers;
            
            size_t start = 0;
            for (size_t i = 0; i < n_finimizers; ++i)
            {
                size_t end = cf.sparse_colors.ends[i];
                
                //colors for this finimizer
                vector<size_t> color_set(cf.sparse_colors.concat.begin() + start, cf.sparse_colors.concat.begin() + end);
                
                color_set_to_finimizers[color_set].push_back(i);
                start = end;
            }
            
            cerr << "Unique color sets: " << color_set_to_finimizers.size() << endl;
            
            this->color_set_ids.resize(n_finimizers);
            vector<pair<vector<size_t>, vector<size_t>>> deduplicated_cs_sparse;
            uint64_t color_set_id = 0;
            for (auto &[color_set, finimizer_indices] : color_set_to_finimizers)
            {
                for (size_t fm_idx : finimizer_indices)
                {
                    this->color_set_ids[fm_idx] = color_set_id;
                }
                color_set_id++;
                
                deduplicated_cs_sparse.emplace_back(color_set, std::move(finimizer_indices));

            }
            CompressedColorSets CCS(deduplicated_cs_sparse, n_colors, this->color_set_ids);
            this->CCS = std::move(CCS);
        } else {
            true_or_crash(cf.color_sets_concat.size() % n_finimizers == 0, "ERROR: color set bitmap length not divisible by finimizer count");
            n_colors = cf.color_sets_concat.size() / n_finimizers;
            cerr << "n_colors: " << (int)n_colors << endl;
            //cerr << sdsl::util::cnt_one_bits(cf.color_sets_concat) << endl;

            cerr << "Deduplicate color sets" << endl;
            const uint64_t *data = cf.color_sets_concat.data();
            vector<pair<sdsl::bit_vector, vector<size_t>>> deduplicated_cs = deduplicate_color_sets(data, n_colors, n_finimizers);

            cerr << "Unique color sets: " << deduplicated_cs.size() << endl;

            this->color_set_ids.resize(n_finimizers);
            CompressedColorSets CCS(deduplicated_cs, n_colors, this->color_set_ids);
            this->CCS = std::move(CCS);
        }

        // this->color_set_ids is in the right order for finimizers
        vector<uint64_t> short_color_set_ids;
        vector<uint64_t> long_color_set_ids;

        vector<std::string_view> short_fmins;
        vector<std::string_view> long_fmins;

        // TODO get rid of short_color_set_ids by overwriting color_set_ids
        //NB: Finimizers are assumed to be in ascending lex order

        int64_t f_start = 0;
        for (int64_t i = 0; i < n_finimizers; i++)
        {
            if (cf.lengths[i] <= short_long_t)
            {
                // fluke8
                short_color_set_ids.push_back(this->color_set_ids[i]);
                short_fmins.push_back(std::string_view(cf.concat.data() + f_start, cf.lengths[i]));
            }
            else
            {
                // prefix tab
                long_color_set_ids.push_back(this->color_set_ids[i]);
                long_fmins.push_back(std::string_view(cf.concat.data() + f_start, cf.lengths[i]));
            }
            f_start += cf.lengths[i];
        }
        uint64_t n_short_fmin = short_color_set_ids.size();
        cerr << "n_short_fmin: "<< n_short_fmin << endl;
        
        short_color_set_ids.insert(short_color_set_ids.end(), long_color_set_ids.begin(), long_color_set_ids.end());
        this->color_set_ids = short_color_set_ids;
        

        // build fluke
        this->f8 = std::move(Fluke8(short_fmins, 0, short_long_t));

        // build prefix table using plen
        this->pt = std::move(PrefTab(long_fmins, plen));

    }

    // Returns the word_idx-th 64-bit word of the color set of the given finimizer.
    // Color c of the set is bit c % 64 of word c / 64. Bits past the last color are zero.
    static inline uint64_t color_set_word(const uint64_t *data, const uint64_t n_colors, const uint64_t finimizer, const uint64_t word_idx)
    {
        const uint64_t n_bits = std::min<uint64_t>(64, n_colors - word_idx * 64);
        const uint64_t pos = finimizer * n_colors + word_idx * 64;
        const uint64_t offset = pos % 64;
        uint64_t word = data[pos / 64] >> offset;
        if (offset + n_bits > 64)
            word |= data[pos / 64 + 1] << (64 - offset);
        if (n_bits < 64)
            word &= (1ULL << n_bits) - 1;
        return word;
    }

    // True if the color set of finimizer a is smaller than that of finimizer b, in the order of
    // sdsl::bit_vector::operator<: the sets are compared color by color, and absent < present.
    static inline bool color_set_less(const uint64_t *data, const uint64_t n_colors, const uint64_t a, const uint64_t b)
    {
        const uint64_t n_words = (n_colors + 63) / 64;
        for (uint64_t w = 0; w < n_words; w++)
        {
            const uint64_t x = color_set_word(data, n_colors, a, w);
            const uint64_t y = color_set_word(data, n_colors, b, w);
            if (x != y)
                return ((x >> __builtin_ctzll(x ^ y)) & 1) == 0;
        }
        return false;
    }

    static inline bool color_set_equal(const uint64_t *data, const uint64_t n_colors, const uint64_t a, const uint64_t b)
    {
        const uint64_t n_words = (n_colors + 63) / 64;
        for (uint64_t w = 0; w < n_words; w++)
        {
            if (color_set_word(data, n_colors, a, w) != color_set_word(data, n_colors, b, w))
                return false;
        }
        return true;
    }

    // Returns the distinct color sets in increasing order (the order of sdsl::bit_vector::operator<),
    // each with the increasing list of the finimizers that have it.
    // Finimizers are grouped by a hash of their color set, so that whole color sets are compared
    // only to confirm that the sets of a group are equal, and to sort the distinct sets.
    vector<pair<sdsl::bit_vector, vector<size_t>>> deduplicate_color_sets(const uint64_t *data, const uint64_t n_colors, const uint64_t n_finimizers)
    {
        const uint64_t n_words = (n_colors + 63) / 64;

        // (hash of color set, finimizer), sorted
        vector<pair<uint64_t, uint64_t>> hashed(n_finimizers);
        for (uint64_t i = 0; i < n_finimizers; i++)
        {
            uint64_t h = 0x9E3779B97F4A7C15ULL;
            for (uint64_t w = 0; w < n_words; w++)
            {
                h = (h ^ color_set_word(data, n_colors, i, w)) * 0xFF51AFD7ED558CCDULL;
                h ^= h >> 32;
            }
            hashed[i] = {h, i};
        }
        std::sort(hashed.begin(), hashed.end());

        // Split into groups of equal color sets: (first finimizer of the group, range in `hashed`)
        struct Group { uint64_t first, begin, end; };
        vector<Group> groups;
        for (uint64_t run_begin = 0; run_begin < n_finimizers;)
        {
            uint64_t run_end = run_begin + 1;
            bool all_equal = true;
            while (run_end < n_finimizers && hashed[run_end].first == hashed[run_begin].first)
            {
                all_equal &= color_set_equal(data, n_colors, hashed[run_begin].second, hashed[run_end].second);
                run_end++;
            }
            if (!all_equal)
            {
                // Different color sets with the same hash: order this run by color set, then finimizer
                std::sort(hashed.begin() + run_begin, hashed.begin() + run_end,
                          [&](const pair<uint64_t, uint64_t> &x, const pair<uint64_t, uint64_t> &y) {
                              if (color_set_less(data, n_colors, x.second, y.second)) return true;
                              if (color_set_less(data, n_colors, y.second, x.second)) return false;
                              return x.second < y.second;
                          });
            }
            for (uint64_t g_begin = run_begin; g_begin < run_end;)
            {
                uint64_t g_end = g_begin + 1;
                while (g_end < run_end && (all_equal || color_set_equal(data, n_colors, hashed[g_begin].second, hashed[g_end].second)))
                    g_end++;
                groups.push_back({hashed[g_begin].second, g_begin, g_end});
                g_begin = g_end;
            }
            run_begin = run_end;
        }

        std::sort(groups.begin(), groups.end(), [&](const Group &x, const Group &y) {
            return color_set_less(data, n_colors, x.first, y.first);
        });

        vector<pair<sdsl::bit_vector, vector<size_t>>> deduplicated_cs;
        deduplicated_cs.reserve(groups.size());
        for (const Group &g : groups)
        {
            sdsl::bit_vector bv(n_colors, 0);
            for (uint64_t w = 0; w < n_words; w++)
                bv.set_int(w * 64, color_set_word(data, n_colors, g.first, w), std::min<uint64_t>(64, n_colors - w * 64));
            vector<size_t> finimizer_indices(g.end - g.begin);
            for (uint64_t j = g.begin; j < g.end; j++)
                finimizer_indices[j - g.begin] = hashed[j].second;
            deduplicated_cs.emplace_back(std::move(bv), std::move(finimizer_indices));
        }
        return deduplicated_cs;
    }

    void read_colors_to_bv(const uint64_t *data, const uint64_t n_colors, const uint64_t start, sdsl::bit_vector &bv)
    {
        sdsl::util::set_to_value(bv, 0);
        const uint64_t *ptr = data + (start * n_colors) / 64;
        uint64_t bit_offset = (start * n_colors) % 64;

        uint64_t color_id = 0;
        uint64_t bits_left = n_colors;

        if (bit_offset != 0)
        {
            // 1. Read the first word
            uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
            uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
            uint64_t word = (*ptr >> bit_offset) & mask;

            while (word != 0)
            {
                uint64_t bit = __builtin_ctzll(word);
                bv[bit] = 1;
                word &= word - 1;
            }

            ++ptr;
            color_id += bits_to_read;
            bits_left -= bits_to_read;
        }

        // 2. Read aligned words in btw
        while (bits_left >= 64)
        {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;)
            {
                uint64_t bit = __builtin_ctzll(w);
                bv[color_id + bit] = 1;
                w &= w - 1;
            }
            color_id += 64;
            bits_left -= 64;
        }

        // 3. Read the last word (if any)
        if (bits_left > 0)
        {
            uint64_t mask = ((1ULL << bits_left) - 1);
            uint64_t word = *ptr & mask;

            while (word != 0)
            {
                uint64_t bit = __builtin_ctzll(word);
                bv[color_id + bit] = 1;
                word &= word - 1;
            }
        }
    }




    void serialize(const std::string &index_prefix) const
    {
        string filename = index_prefix + ".fmin";
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            cerr << "Error: Could not open file for writing: " << filename << endl;
            return;
        }

        cerr << "Save the index to " << filename << endl;

        CCS.serialize(out);

        uint64_t color_set_ids_size = color_set_ids.size();
        out.write(reinterpret_cast<const char *>(&color_set_ids_size), sizeof(uint64_t));
        out.write(reinterpret_cast<const char *>(color_set_ids.data()), color_set_ids_size * sizeof(uint64_t));

        f8.serialize(out);
        pt.serialize(out);

        // metadata
        out.write(reinterpret_cast<const char *>(&n_colors), sizeof(n_colors));
        out.write(reinterpret_cast<const char *>(&n_finimizers), sizeof(n_finimizers)); // this is not needed
        out.write(reinterpret_cast<const char *>(&k), sizeof(k));

        out.close();
        cerr << "DONE" << endl;
    }
    // Loads the index from PREFIX.fmin. The file is mapped into memory and its large arrays
    // are used in place; if it cannot be mapped, it is read in the ordinary way.
    void load(const std::string &index_prefix)
    {
        string filename = index_prefix + ".fmin";
        std::ifstream file(filename, std::ios::binary);
        if (!file)
        {
            throw std::runtime_error("Could not open the index file " + filename);
        }
        std::unique_ptr<MemoryStreamBuf> memory;
        try
        {
            mapping = std::make_shared<MappedFile>(filename);
            memory = std::make_unique<MemoryStreamBuf>(mapping->data(), mapping->size());
        }
        catch (const std::runtime_error &)
        {
            mapping.reset(); // Read from `file` instead
        }
        std::istream in(memory ? static_cast<std::streambuf *>(memory.get()) : file.rdbuf());

        // The file has no signature, so check what can be checked. A read past the end of
        // the file raises an exception at once, before a half-read structure can be used.
        try
        {
            in.exceptions(std::ios::failbit | std::ios::badbit | std::ios::eofbit);
            cerr << "Loading index from " << filename << endl;

            CCS.load(in);

            uint64_t color_set_ids_size;
            in.read(reinterpret_cast<char *>(&color_set_ids_size), sizeof(uint64_t));
            cerr << "color_set_ids_size: " << color_set_ids_size << endl; // Do we need this?

            loaded_color_set_ids.load(in, color_set_ids_size);

            f8.load(in);
            pt.load(in);

            // metadata
            in.read(reinterpret_cast<char *>(&n_colors), sizeof(n_colors));
            //this->n_colors = 10000;
            cerr << "n_colors: " << n_colors << endl;
            in.read(reinterpret_cast<char *>(&n_finimizers), sizeof(n_finimizers));
            cerr << "n_finimizers: " << n_finimizers << endl; // this is not needed
            //in.read(reinterpret_cast<char *>(&plen), sizeof(plen));
            in.read(reinterpret_cast<char *>(&k), sizeof(k));

            in.exceptions(std::ios::goodbit);
        }
        catch (const std::ios_base::failure &)
        {
            throw std::runtime_error(filename + " is not a complete Finimap index (the file ends too early)");
        }
        catch (const std::bad_alloc &)
        {
            throw std::runtime_error(filename + " is not a valid Finimap index");
        }
        catch (const std::length_error &)
        {
            throw std::runtime_error(filename + " is not a valid Finimap index");
        }

        if (in.peek() != std::ifstream::traits_type::eof())
            throw std::runtime_error(filename + " is not a Finimap index (unexpected data after the end)");
        if (n_colors == 0 || n_colors > MAX_COLORS)
            throw std::runtime_error(filename + " is not a valid Finimap index (number of colors: " + to_string(n_colors) + ")");
        if (k == 0 || k > MAX_K)
            throw std::runtime_error(filename + " is not a valid Finimap index (k = " + to_string(k) + ")");

        cerr << "DONE" << endl;
    }
};
