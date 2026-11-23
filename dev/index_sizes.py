#!/usr/bin/env python3
"""Report the size of each part of a Finimap index file (PREFIX.fmin), and of the k-mer filter
beside it if there is one. Reads only the headers of the parts.

usage: index_sizes.py PREFIX [PREFIX ...]
"""
import os, struct, sys

def u64(f):
    return struct.unpack("<Q", f.read(8))[0]

def report(prefix):
    path = prefix + ".fmin"
    total = os.path.getsize(path)
    parts = []
    with open(path, "rb") as f:
        def skip(name, nbytes, note=""):
            parts.append((name, nbytes, note))
            f.seek(nbytes, 1)

        fmt = 1
        if f.read(8) == b"FINIMAP!":
            fmt = u64(f)
            if fmt not in (2, 3, 4):
                sys.exit("%s: index format %d is not known to this script" % (path, fmt))
        else:
            f.seek(0)
        sparse_count = dense_count = bits = 0
        meta = None
        if fmt >= 4:
            mc_colors, mc_sets, mc_blocks, mc_row_bits, mc_partials = [u64(f) for _ in range(5)]
            skip("color sets: meta-color layout (color order, blocks, codes)", mc_colors * 2 + (37 * mc_blocks + 2) * 8)
            rec = u64(f)
            skip("color sets: partial set records (offsets, lengths, forms)", rec)
            pd = u64(f)
            skip("color sets: partial sets", pd, "%d distinct partial sets in %d blocks" % (mc_partials, mc_blocks))
            rb_n, rb_bits, rb_bytes = u64(f), u64(f), u64(f)
            skip("color sets: row offsets, every 16th row", rb_bytes, "%d bits each" % rb_bits)
            rr_n, rr_bits, rr_bytes = u64(f), u64(f), u64(f)
            skip("color sets: row offsets within groups of 16", rr_bytes, "%d bits each" % rr_bits)
            rows = u64(f)
            skip("color sets: rows (blocks off the default, Huffman-coded partial set ids)", rows, "%d color sets, %.1f bits each on average" % (mc_sets, mc_row_bits / max(1, mc_sets)))
            meta = (mc_sets, mc_blocks, mc_partials)
        else:
            sparse_count, dense_count = u64(f), u64(f)
            bits = u64(f)
            skip("color sets: bitmaps (dense sets)", (bits + 63) // 64 * 8, "%d bits" % bits)
            ef_n, ef_period = u64(f), u64(f)
            skip("color sets: list offsets", (ef_n // ef_period + 1) * 8 + ef_n * 2, "%d lists" % ef_n)
        if fmt >= 4:
            n_ids, id_bits, id_bytes = u64(f), u64(f), u64(f)
            skip("color set id of each finimizer", id_bytes, "%d finimizers, %d bits each" % (n_ids, id_bits))
        elif fmt >= 3:
            l_size, l_bits, l_bytes = u64(f), u64(f), u64(f)
            skip("color sets: lists (sparse sets, complements of very dense ones)", l_bytes, "%d colors listed, %d bits each" % (l_size, l_bits))
            n_ids, id_bits, id_bytes = u64(f), u64(f), u64(f)
            skip("color set id of each finimizer", id_bytes, "%d finimizers, %d bits each" % (n_ids, id_bits))
        else:
            l_size = u64(f)
            skip("color sets: lists (sparse sets, complements of very dense ones)", l_size * 2, "%d colors listed" % l_size)
            n_ids = u64(f)
            skip("color set id of each finimizer", n_ids * 8, "%d finimizers" % n_ids)
        f8_u, f8_n, f8_min, upperlen, lowerlen = [u64(f) for _ in range(5)]
        if fmt == 1:
            p_u, p_n, p_min, nblocks = [u64(f) for _ in range(4)]
            skip("short finimizers: bucket table (size set by -x)", (nblocks + 1) * 4, "%d buckets" % nblocks)
            u64(f)
            skip("short finimizers: low bytes", p_n, "%d short finimizers (length <= %d)" % (p_n, upperlen))
        else:
            p_n, p_min, p_max, nblocks, nb2 = [u64(f) for _ in range(5)]
            l2 = u64(f)
            skip("short finimizers: bucket table, Pino blocks (size set by -x)", l2, "%d buckets" % nblocks)
            h2 = u64(f)
            skip("short finimizers: bucket table, Pino counts", h2 * 4)
            l = u64(f)
            skip("short finimizers: low bytes", l, "%d short finimizers (length <= %d)" % (p_n, upperlen))
        if fmt >= 3:
            len_n, len_bits, len_bytes = u64(f), u64(f), u64(f)
            skip("short finimizers: lengths", len_bytes, "%d bits each" % len_bits)
        else:
            skip("short finimizers: lengths", f8_n)
        preflen, nprefs, nfins = u64(f), u64(f), u64(f)
        skip("long finimizers: prefix table (size set by -p)", (nprefs + 1) * 4, "%d prefixes of length %d" % (nprefs, preflen))
        if fmt >= 4 and 5 <= preflen <= 30:
            k_n, k_bits, k_bytes = u64(f), u64(f), u64(f)
            skip("long finimizers: keys", k_bytes, "%d long finimizers, %d bits each" % (nfins, k_bits))
        else:
            skip("long finimizers: keys", nfins * 8, "%d long finimizers" % nfins)
        n_colors, n_finimizers, k = u64(f), u64(f), u64(f)
        end = f.tell()
    accounted = sum(p[1] for p in parts)
    n_bitmaps = bits // n_colors if n_colors else 0
    if meta:
        print("%s: %.3f GB, format %d, %d colors, %d finimizers, k = %d | distinct color sets: %d, as meta-color sets in %d blocks"
              % (path, total / 1e9, fmt, n_colors, n_finimizers, k, meta[0], meta[1]))
    else:
        print("%s: %.3f GB, format %d, %d colors, %d finimizers, k = %d | distinct color sets: %d listed (sparse), %d listed as complements (very dense), %d as bitmaps"
              % (path, total / 1e9, fmt, n_colors, n_finimizers, k, sparse_count, dense_count - sparse_count, n_bitmaps))
    for name, nbytes, note in parts:
        print("  %8.1f MB  %5.1f%%  %s%s" % (nbytes / 1e6, 100.0 * nbytes / total, name, ("  (" + note + ")") if note else ""))
    print("  %8.1f MB  %5.1f%%  headers" % ((total - accounted) / 1e6, 100.0 * (total - accounted) / total))
    if end != total:
        print("  WARNING: the parts end at byte %d, but the file has %d bytes" % (end, total))
    if os.path.exists(path + ".filter"):
        print("  filter beside it: %s, %.1f MB" % (path + ".filter", os.path.getsize(path + ".filter") / 1e6))

for prefix in sys.argv[1:]:
    report(prefix[:-5] if prefix.endswith(".fmin") else prefix)
