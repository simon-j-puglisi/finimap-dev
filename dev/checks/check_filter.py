#!/usr/bin/env python3
"""Check the found-in-index mode with the real k-mer filter against the exact k-mer set.

usage: check_filter.py BINARY INDEX_PREFIX REFS_LIST QUERIES K

Runs search-fmin --found-in-index twice at -t 0, once with the filter PREFIX.fmin.filter and
once with --exact-kmers, and checks that the filter never loses anything: for every read,
n_pos and every color's score are at least those of the exact set. Then it reports how many
k-mers that are not in the index passed the filter, which should match the filter's
false-positive rate.
"""
import os, re, subprocess, sys, tempfile
binary, index, refs, queries, k = sys.argv[1:6]
k = int(k)
tmp = tempfile.mkdtemp()

def run(extra, name):
    out = os.path.join(tmp, name)
    subprocess.check_call([binary, "search-fmin", "-i", index, "-q", queries, "-t", "0", "--found-in-index",
                           "--print-denominator", "-o", out] + extra, stderr=subprocess.DEVNULL)
    rows = []
    for line in open(out):
        f = line.split()
        rows.append((int(f[1]), dict((int(c), int(s)) for c, s in (x.split(":") for x in f[2:]))))
    return rows

reads, cur = [], None
for line in open(queries):
    line = line.strip()
    if line.startswith(">"):
        reads.append("")
    elif reads:
        reads[-1] += line
valid = re.compile("[ACGTacgt]{%d}" % k)
n_valid = [sum(1 for i in range(len(r) - k + 1) if valid.fullmatch(r, i, i + k)) for r in reads]

exact = run(["--exact-kmers", refs], "exact.txt")
bloom = run([], "filter.txt")
assert len(exact) == len(bloom) == len(reads)
lost_npos = lost_score = extra_pairs = 0
absent = passed_absent = 0
score_exact = score_bloom = 0
for (ne, se), (nb, sb), nv in zip(exact, bloom, n_valid):
    lost_npos += nb < ne
    lost_score += any(sb.get(c, 0) < s for c, s in se.items())
    extra_pairs += sum(1 for c in sb if c not in se)
    absent += nv - ne
    passed_absent += nb - ne
    score_exact += sum(se.values())
    score_bloom += sum(sb.values())
print("reads %d | k-mers: %d valid, %d in the index, %d not" % (len(reads), sum(n_valid), sum(n for n, _ in exact), absent))
print("reads with n_pos below the exact set's: %d | reads with a score below the exact set's: %d" % (lost_npos, lost_score))
print("k-mers not in the index that passed the filter: %d of %d (%.3f%%)" % (passed_absent, absent, 100.0 * passed_absent / max(absent, 1)))
print("total score: exact set %d, filter %d (+%.3f%%) | read-colour pairs reported only with the filter: %d"
      % (score_exact, score_bloom, 100.0 * (score_bloom - score_exact) / max(score_exact, 1), extra_pairs))
sys.exit(1 if (lost_npos or lost_score) else 0)
