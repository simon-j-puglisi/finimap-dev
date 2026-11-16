#!/usr/bin/env python3
"""Check the found-in-index mode against a reference built from the default mode.

usage: check_found.py NEW_BINARY DEFAULT_BINARY INDEX_PREFIX REFS_LIST QUERIES K TAU [TAU ...]

Reference, per read: P = the read's k-mers that are in the references; each k-mer of P is
queried on its own with DEFAULT_BINARY at -t 1 and the votes are added up per colour; colours
with score / |P| >= tau (and score >= 1) are reported, by decreasing score then colour.
The new mode, with the exact k-mer set as its filter, must print the same lines.

Needs the foundref program (foundref.cpp) built next to this script.
"""
import os, subprocess, sys, tempfile
new, default, index, refs, queries, k = sys.argv[1:7]
taus = sys.argv[7:]
here = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp()
kmers, npos_file = os.path.join(tmp, "kmers.fna"), os.path.join(tmp, "npos.txt")
subprocess.check_call([os.path.join(here, "foundref"), refs, queries, k, kmers, npos_file])
npos = [int(l) for l in open(npos_file)]
owner = [int(l[1:]) for l in open(kmers) if l.startswith(">")]
votes_out = os.path.join(tmp, "votes.txt")
scores = [dict() for _ in npos]
if owner:
    p = subprocess.run([default, "search-fmin", "-i", index, "-q", kmers, "-t", "1"], stdout=open(votes_out, "w"), stderr=subprocess.DEVNULL)
    if p.returncode != 0:
        sys.exit("default binary failed on the k-mer file (exit %d)" % p.returncode)
    lines = [l.split() for l in open(votes_out) if l[:1].isdigit()]
    assert len(lines) == len(owner), (len(lines), len(owner))
    for j, f in enumerate(lines):
        assert int(f[0]) == j % 65536 or int(f[0]) == j
        for tok in f[1:]:
            c, s = tok.split(":")
            assert s == "1", tok
            d = scores[owner[j]]
            d[int(c)] = d.get(int(c), 0) + 1
bad = 0
for tau in taus:
    t = float(tau)
    want = []
    for r, n in enumerate(npos):
        items = [(c, s) for c, s in scores[r].items() if n > 0 and s >= 1 and s / n >= t]
        items.sort(key=lambda x: (-x[1], x[0]))
        want.append("%d %d %s" % (r, n, "".join("%d:%d " % x for x in items)))
    out = os.path.join(tmp, "new.txt")
    subprocess.check_call([new, "search-fmin", "-i", index, "-q", queries, "-t", tau, "--found-in-index", "--exact-kmers", refs, "--print-denominator", "-o", out], stderr=subprocess.DEVNULL)
    got = [l.rstrip("\n") for l in open(out)]
    diff = [i for i in range(max(len(got), len(want))) if i >= len(got) or i >= len(want) or got[i] != want[i]]
    reported = sum(l.count(":") for l in want)
    print("tau=%s reads=%d k-mers in index=%d of all, reported pairs=%d, differing lines=%d" % (tau, len(want), sum(npos), reported, len(diff)))
    for i in diff[:3]:
        print("  read %d\n    want: %s\n    got:  %s" % (i, want[i][:160] if i < len(want) else None, got[i][:160] if i < len(got) else None))
    bad += len(diff)
sys.exit(1 if bad else 0)
