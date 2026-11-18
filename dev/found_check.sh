#!/bin/bash -l
# The found-in-index mode at scale: build k-mer filters (timed), then run every query file in
# the default mode and in the found-in-index mode with each filter, minimum run length and
# threshold (timed), and report times, n_pos and how many colours are reported.
# The index, the reference files and the query files are only read. Safe to resubmit.
#
# Settings (environment variables):
#   FM_REPO      path of this repository checkout (required)
#   FM_INDEX     index prefix, i.e. the path without .fmin (required)
#   FM_REFS      the list of reference files the index was built from (required)
#   FM_QUERIES   directory of query files *.fna (required; *.templates.fna are ignored)
#   FM_REF_OUT   directory with reference results of the default mode, NAME.txt (optional)
#   FM_OUT       output directory (default: the current directory)
#   FM_FILTERS   filters to build: "bits:hashes" pairs for Bloom filters, and "exact" for the
#                exact k-mer table (default "10:6 16:8 24:10")
#   FM_MIN_RUNS  values of --min-run (default "1")
#   FM_TAUS      thresholds for the found-in-index mode (default "0.8 0.9 0.96 1")
#   FM_TAU       threshold for the default mode (default 0.8)
#   FM_K         k of the index (default 31); the filters are built with the same k
#   FM_CMAKE_MODULE, FM_GCC_MODULE   as in stage2.sh
#
# Outputs: default/NAME.txt, and found/NAME.F.rR.tTAU.txt for filter F (b10, b16, ..., or
# exact), minimum run R and threshold TAU, written with --print-denominator: the second field
# of each line is n_pos.
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
INDEX=${FM_INDEX:?set FM_INDEX}
REFS=$(readlink -f "${FM_REFS:?set FM_REFS}")
QUERIES=$(readlink -f "${FM_QUERIES:?set FM_QUERIES}")
REF_OUT=${FM_REF_OUT:-}
OUT=$(readlink -f "${FM_OUT:-$PWD}")
FILTERS=${FM_FILTERS:-10:6 16:8 24:10}
MIN_RUNS=${FM_MIN_RUNS:-1}
TAUS=${FM_TAUS:-0.8 0.9 0.96 1}
TAU=${FM_TAU:-0.8}
K=${FM_K:-31}
JOBS=${SLURM_CPUS_PER_TASK:-4}
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)

echo "FM_FOUND_CHECK start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$COMMIT index=$INDEX refs=$REFS ($(wc -l < "$REFS") files) queries=$QUERIES out=$OUT k=$K filters='$FILTERS' min_runs='$MIN_RUNS' taus='$TAUS' threads=$JOBS"
[ -e "$REPO/dev/timeit.py" ] || { echo "FM_ERROR repo not found at $REPO"; exit 1; }
[ -s "$INDEX.fmin" ] || { echo "FM_ERROR index not found: $INDEX.fmin"; exit 1; }
mkdir -p "$OUT/default" "$OUT/found"

. "$REPO/dev/build_finimap.sh"
fm_load_modules
fm_build_finimap
if [ "$FM_REBUILT" = 1 ]; then rm -f "$OUT"/*.time.txt "$OUT"/default/*.time.txt "$OUT"/found/*.time.txt "$OUT"/filter_*.fmin.filter; fi
BIN=$OUT/finimap

# Reference paths may be relative to the directory of the list
if [ ! -e "$(head -n 1 "$REFS")" ]; then cd "$(dirname "$REFS")"; fi
[ -e "$(head -n 1 "$REFS")" ] || { echo "FM_ERROR first reference not found: $(head -n 1 "$REFS")"; exit 1; }

name_of() { if [ "$1" = exact ]; then echo exact; else echo "b${1%:*}"; fi; }

echo "--- building the filters"
distinct=""
for f in $FILTERS; do
    F=$(name_of "$f")
    P=$OUT/filter_$F
    if [ -s "$P.time.txt" ] && grep -q '^exit_code 0$' "$P.time.txt" && [ -s "$P.fmin.filter" ]; then
        echo "filter $F: already built"
    else
        # The first build estimates the number of distinct k-mers; the others reuse it
        extra=""; [ -z "$distinct" ] || extra="--distinct-kmers $distinct"
        if [ "$f" = exact ]; then kind="--exact"; else kind="--bits-per-kmer ${f%:*} --hashes ${f#*:}"; fi
        python3 "$REPO/dev/timeit.py" "$P.time.part" "$BIN" build-filter -i "$REFS" -o "$P" -k "$K" \
            $kind --threads "$JOBS" $extra 2> "$P.log" \
            || { echo "FM_ERROR build-filter failed:"; tail -n 5 "$P.log"; exit 1; }
        mv "$P.time.part" "$P.time.txt"
    fi
    [ -n "$distinct" ] || distinct=$(grep -oE 'about [0-9]+' "$P.log" | grep -oE '[0-9]+' | head -n 1)
    echo "$F | $(grep -E 'Distinct|False-positive|Time \(s\)|Wrote' "$P.log" | tr '\n' '|') $(tr '\n' ' ' < "$P.time.txt")"
done
cd "$OUT"

run() {  # name of the run, then the search-fmin arguments
    local name=$1; shift
    if [ -s "$name.time.txt" ] && grep -q '^exit_code 0$' "$name.time.txt"; then return 0; fi
    python3 "$REPO/dev/timeit.py" "$name.time.part" "$BIN" search-fmin -i "$INDEX" -o "$name.txt.part" "$@" > /dev/null 2> "$name.stderr" \
        || { echo "FM_ERROR search-fmin failed for $name:"; tail -n 5 "$name.stderr"; exit 1; }
    mv "$name.txt.part" "$name.txt"
    mv "$name.time.part" "$name.time.txt"
}

echo "--- running the queries"
n_files=0
for q in "$QUERIES"/*.fna; do
    case "$q" in *.templates.fna) continue ;; esac
    [ -e "$q" ] || continue
    name=$(basename "$q" .fna)
    n_files=$((n_files + 1))
    run "default/$name" -q "$q" -t "$TAU"
    for f in $FILTERS; do
        F=$(name_of "$f")
        for R in $MIN_RUNS; do
            for t in $TAUS; do
                run "found/$name.$F.r$R.t$t" -q "$q" -t "$t" --found-in-index --filter "$OUT/filter_$F.fmin.filter" --min-run "$R" --print-denominator
            done
        done
    done
done
[ "$n_files" -gt 0 ] || { echo "FM_ERROR no query files in $QUERIES"; exit 1; }
echo "ran $n_files query files"

echo "--- report"
python3 - "$QUERIES" "$OUT" "$REF_OUT" "$FILTERS" "$MIN_RUNS" "$TAUS" <<'PY'
import os, re, sys
queries, out, ref_out, filters, min_runs, taus = sys.argv[1:7]
names_f = ["exact" if f == "exact" else "b" + f.split(":")[0] for f in filters.split()]
min_runs = min_runs.split()
taus = taus.split()

def field(path, key):
    for l in open(path):
        if l.startswith(key):
            return l.split()[-1]
    return "n/a"

def part(path, label):
    m = re.search(label + r" ([0-9.]+)", open(path).read())
    return m.group(1) if m else "n/a"

def summary(path):
    n_pos = pairs = hit_reads = 0
    for l in open(path):
        p = l.split()
        n_pos += int(p[1])
        pairs += len(p) - 2
        hit_reads += len(p) > 2
    return n_pos, pairs, hit_reads

names = sorted(n[:-4] for n in os.listdir(queries) if n.endswith(".fna") and not n.endswith(".templates.fna"))
unexpected = 0
for name in names:
    d = os.path.join(out, "default", name)
    default_lines = [" ".join(l.split()) for l in open(d + ".txt")]
    line = "FM_FOUND %s reads=%d | default: wall_s=%s search_s=%s" % (name, len(default_lines), field(d + ".time.txt", "wall_s"), part(d + ".stderr", "search"))
    if ref_out and os.path.exists(os.path.join(ref_out, name + ".txt")):
        ref = [" ".join(l.split()) for l in open(os.path.join(ref_out, name + ".txt"))]
        differ = sum(1 for a, b in zip(default_lines, ref) if a != b) + abs(len(default_lines) - len(ref))
        unexpected += differ
        line += " lines_differing_from_reference=%d" % differ
    print(line)
    for R in min_runs:
        exact_npos = None
        if "exact" in names_f:
            exact_npos = summary(os.path.join(out, "found", "%s.exact.r%s.t%s.txt" % (name, R, taus[0])))[0]
        for F in names_f:
            cells = []
            for t in taus:
                n_pos, pairs, hit_reads = summary(os.path.join(out, "found", "%s.%s.r%s.t%s.txt" % (name, F, R, t)))
                cells.append("t=%s reads_with_a_colour=%d pairs=%d" % (t, hit_reads, pairs))
            f0 = os.path.join(out, "found", "%s.%s.r%s.t%s" % (name, F, R, taus[0]))
            m = re.search(r"passed the filter: ([0-9]+)(?:, of which in runs shorter than [0-9]+ \(not counted\): ([0-9]+))?", open(f0 + ".stderr").read())
            passed, dropped = (int(m.group(1)), int(m.group(2) or 0)) if m else (0, 0)
            excess = "" if exact_npos is None else " excess_over_exact=%d" % (n_pos - exact_npos)
            print("  %s r=%s n_pos_total=%d%s passed=%d dropped_by_min_run=%d wall_s=%s search_s=%s filter_s=%s | %s"
                  % (F, R, n_pos, excess, passed, dropped, field(f0 + ".time.txt", "wall_s"), part(f0 + ".stderr", "search"), part(f0 + ".stderr", "k-mer filter"), " | ".join(cells)))
print("FM_FOUND_TOTAL files=%d default_lines_differing_from_reference=%d" % (len(names), unexpected))
PY
ls -l "$OUT"/filter_*.fmin.filter "$INDEX.fmin" | awk '{print $5, $9}'
echo "FM_FOUND_CHECK_DONE"
