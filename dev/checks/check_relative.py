#!/usr/bin/env python3
"""Check search-fmin --relative against the rule derived from a -t 0 run.

usage: check_relative.py BINARY INDEX_PREFIX QUERIES [extra search-fmin options ...]

Runs search-fmin once at -t 0 with --print-denominator, which lists every color with a score
of at least 1, and derives from it, for several thresholds t and deltas d, the colors that
must be reported: score >= top - d, where top is the read's highest score, and the threshold
rule of the mode (default: score >= max(floor(found * t), 1) in single precision; with
--found-in-index: score / n_pos >= t as real numbers). Then it runs search-fmin with -t t
--relative d for each pair and compares. Extra options (such as --found-in-index and its
filter) are passed to every run.
"""
import os, struct, subprocess, sys, tempfile

binary, index, queries = sys.argv[1:4]
extra = sys.argv[4:]
found_mode = "--found-in-index" in extra
tmp = tempfile.mkdtemp()

def f32(x):
    return struct.unpack("f", struct.pack("f", x))[0]

def min_score(t, denom):
    if not found_mode:
        return max(int(f32(f32(float(denom)) * f32(t))), 1)
    if denom <= 0:
        return 1
    n = float(denom)
    s = max(int(-(-t * n // 1)), 1)
    while s > 1 and (s - 1) / n >= t:
        s -= 1
    while s <= denom and s / n < t:
        s += 1
    return s

def run(args, name):
    out = os.path.join(tmp, name)
    subprocess.check_call([binary, "search-fmin", "-i", index, "-q", queries, "-o", out] + args + extra, stderr=subprocess.DEVNULL)
    return [l.rstrip("\n") for l in open(out)]

base = []
for line in run(["-t", "0", "--print-denominator"], "base.txt"):
    f = line.split()
    base.append((int(f[0]), int(f[1]), [(int(c), int(s)) for c, s in (x.split(":") for x in f[2:])]))

bad = 0
for t in ("0", "0.8", "0.96"):
    for d in (0, 1, 31, 62, 155):
        want = []
        for rid, denom, pairs in base:
            top = max((s for _, s in pairs), default=0)
            m = max(min_score(float(t), denom), top - d)
            want.append("%d %s" % (rid, "".join("%d:%d " % (c, s) for c, s in pairs if s >= m)))
        got = run(["-t", t, "--relative", str(d)], "rel.txt")
        diff = sum(1 for a, b in zip(want, got) if a != b) + abs(len(want) - len(got))
        pairs = sum(w.count(":") for w in want)
        print("t=%s delta=%d: reads=%d reported pairs=%d differing lines=%d" % (t, d, len(want), pairs, diff))
        bad += diff
sys.exit(1 if bad else 0)
