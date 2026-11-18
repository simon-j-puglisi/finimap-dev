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

        sparse_count, dense_count = u64(f), u64(f)
        bits = u64(f)
        skip("color sets: bitmaps (dense sets)", (bits + 63) // 64 * 8, "%d bits" % bits)
        ef_n, ef_period = u64(f), u64(f)
        skip("color sets: list offsets", (ef_n // ef_period + 1) * 8 + ef_n * 2, "%d lists" % ef_n)
        l_size = u64(f)
        skip("color sets: lists (sparse sets, complements of very dense ones)", l_size * 2, "%d colors listed" % l_size)
        n_ids = u64(f)
        skip("color set id of each finimizer", n_ids * 8, "%d finimizers" % n_ids)
        f8_u, f8_n, f8_min, upperlen, lowerlen = [u64(f) for _ in range(5)]
        p_u, p_n, p_min, nblocks = [u64(f) for _ in range(4)]
        skip("short finimizers: bucket table (size set by -x)", (nblocks + 1) * 4, "%d buckets" % nblocks)
        u64(f)
        skip("short finimizers: low bytes", p_n, "%d short finimizers (length <= %d)" % (p_n, upperlen))
        skip("short finimizers: lengths", f8_n)
        preflen, nprefs, nfins = u64(f), u64(f), u64(f)
        skip("long finimizers: prefix table (size set by -p)", (nprefs + 1) * 4, "%d prefixes of length %d" % (nprefs, preflen))
        skip("long finimizers: keys", nfins * 8, "%d long finimizers" % nfins)
        n_colors, n_finimizers, k = u64(f), u64(f), u64(f)
        end = f.tell()
    accounted = sum(p[1] for p in parts)
    n_bitmaps = bits // n_colors if n_colors else 0
    print("%s: %.3f GB, %d colors, %d finimizers, k = %d | distinct color sets: %d listed (sparse), %d listed as complements (very dense), %d as bitmaps"
          % (path, total / 1e9, n_colors, n_finimizers, k, sparse_count, dense_count - sparse_count, n_bitmaps))
    for name, nbytes, note in parts:
        print("  %8.1f MB  %5.1f%%  %s%s" % (nbytes / 1e6, 100.0 * nbytes / total, name, ("  (" + note + ")") if note else ""))
    print("  %8.1f MB  %5.1f%%  headers" % ((total - accounted) / 1e6, 100.0 * (total - accounted) / total))
    if end != total:
        print("  WARNING: the parts end at byte %d, but the file has %d bytes" % (end, total))
    if os.path.exists(path + ".filter"):
        print("  filter beside it: %s, %.1f MB" % (path + ".filter", os.path.getsize(path + ".filter") / 1e6))

for prefix in sys.argv[1:]:
    report(prefix[:-5] if prefix.endswith(".fmin") else prefix)
