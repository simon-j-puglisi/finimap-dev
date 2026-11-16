#!/bin/bash -l
# Run search-fmin on every query file of a directory, timed, and compare each result with a
# reference output of an earlier version. Differences are split into reads that hold only
# A, C, G, T (none are expected) and reads with other characters (expected: k-mers that
# contain such a character are now skipped). The index and the reference files are only read.
# Safe to resubmit: finished query files are skipped.
#
# Settings (environment variables):
#   FM_REPO      path of this repository checkout (required)
#   FM_INDEX     index prefix, i.e. the path without .fmin (required)
#   FM_QUERIES   directory of query files *.fna (required; *.templates.fna are ignored)
#   FM_REF_OUT   directory with the reference result of each query file, NAME.txt (required)
#   FM_REF_TIME  directory with the reference timings, finimap_query_NAME.txt (optional)
#   FM_OUT       output directory (default: the current directory)
#   FM_TAU       threshold (default 0.8)
#   FM_THREADS   threads for search-fmin (default 1)
#   FM_CMAKE_MODULE, FM_GCC_MODULE   as in stage2.sh
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
INDEX=${FM_INDEX:?set FM_INDEX}
QUERIES=$(readlink -f "${FM_QUERIES:?set FM_QUERIES}")
REF_OUT=$(readlink -f "${FM_REF_OUT:?set FM_REF_OUT}")
REF_TIME=${FM_REF_TIME:-}
OUT=$(readlink -f "${FM_OUT:-$PWD}")
TAU=${FM_TAU:-0.8}
THREADS=${FM_THREADS:-1}
JOBS=${SLURM_CPUS_PER_TASK:-4}
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)

echo "FM_QUERY_CHECK start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$COMMIT index=$INDEX queries=$QUERIES ref_out=$REF_OUT out=$OUT tau=$TAU threads=$THREADS"
[ -e "$REPO/dev/timeit.py" ] || { echo "FM_ERROR repo not found at $REPO"; exit 1; }
[ -s "$INDEX.fmin" ] || { echo "FM_ERROR index not found: $INDEX.fmin"; exit 1; }
mkdir -p "$OUT"

. "$REPO/dev/build_finimap.sh"
fm_load_modules
fm_build_finimap
if [ "$FM_REBUILT" = 1 ]; then rm -f "$OUT"/*.time.txt; fi

n_files=0
for q in "$QUERIES"/*.fna; do
    case "$q" in *.templates.fna) continue ;; esac
    [ -e "$q" ] || continue
    name=$(basename "$q" .fna)
    n_files=$((n_files + 1))
    if [ -s "$OUT/$name.time.txt" ] && grep -q '^exit_code 0$' "$OUT/$name.time.txt"; then
        echo "$name: already done"
    else
        python3 "$REPO/dev/timeit.py" "$OUT/$name.time.part" \
            "$OUT/finimap" search-fmin -i "$INDEX" -q "$q" -o "$OUT/$name.txt.part" -t "$TAU" --threads "$THREADS" \
            > "$OUT/$name.stdout" 2> "$OUT/$name.stderr" \
            || { echo "FM_ERROR search-fmin failed on $name; end of its messages:"; tail -n 5 "$OUT/$name.stderr"; exit 1; }
        mv "$OUT/$name.txt.part" "$OUT/$name.txt"
        mv "$OUT/$name.time.part" "$OUT/$name.time.txt"
    fi
done
[ "$n_files" -gt 0 ] || { echo "FM_ERROR no query files in $QUERIES"; exit 1; }

echo "--- comparison with the reference results"
python3 - "$QUERIES" "$REF_OUT" "$OUT" "$REF_TIME" <<'PY'
import os, re, sys
queries, ref_out, out, ref_time = sys.argv[1:5]
acgt = re.compile("^[ACGTacgt]*$")

def reads_of(path):
    reads, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if cur is not None:
                    reads.append("".join(cur))
                cur = []
            elif cur is not None:
                cur.append(line)
    if cur is not None:
        reads.append("".join(cur))
    return reads

def lines_of(path):
    with open(path) as f:
        return [" ".join(l.split()) for l in f]

def field(path, key):
    if not path or not os.path.exists(path):
        return "n/a"
    for l in open(path):
        if key in l:
            return l.strip().split()[-1]
    return "n/a"

total_bad = 0
names = sorted(n[:-4] for n in os.listdir(queries) if n.endswith(".fna") and not n.endswith(".templates.fna"))
for name in names:
    reads = reads_of(os.path.join(queries, name + ".fna"))
    new = lines_of(os.path.join(out, name + ".txt"))
    ref_path = os.path.join(ref_out, name + ".txt")
    t_new = field(os.path.join(out, name + ".time.txt"), "wall_s")
    m_new = field(os.path.join(out, name + ".time.txt"), "peak_rss_gb")
    t_ref = field(os.path.join(ref_time, "finimap_query_%s.txt" % name) if ref_time else "", "Elapsed")
    lens = sorted(set(len(r) for r in reads))
    head = "FM_QCHECK %s reads=%d lengths=%s new_lines=%d" % (name, len(reads), lens[0] if len(lens) == 1 else "%d..%d" % (lens[0], lens[-1]), len(new))
    stdout_bytes = os.path.getsize(os.path.join(out, name + ".stdout"))
    zero = sum(1 for l in new if re.search(r":0( |$)", l))
    if not os.path.exists(ref_path):
        print(head + " NO REFERENCE wall_s=%s peak_gb=%s" % (t_new, m_new))
        continue
    ref = lines_of(ref_path)
    n = min(len(new), len(ref), len(reads))
    same = diff_other = diff_acgt = 0
    examples = []
    for i in range(n):
        if new[i] == ref[i]:
            same += 1
        elif acgt.match(reads[i]):
            diff_acgt += 1
            if len(examples) < 2:
                examples.append((i, ref[i][:150], new[i][:150]))
        else:
            diff_other += 1
    with_other = sum(1 for r in reads if not acgt.match(r))
    bad = diff_acgt + abs(len(new) - len(reads)) + (1 if stdout_bytes else 0) + zero
    total_bad += bad
    print(head + " ref_lines=%d identical=%d differ_reads_with_other_chars=%d (of %d such reads) differ_acgt_only=%d zero_scores=%d stdout_bytes=%d wall_s=%s ref_wall_s=%s peak_gb=%s"
          % (len(ref), same, diff_other, with_other, diff_acgt, zero, stdout_bytes, t_new, t_ref, m_new))
    for i, r, w in examples:
        print("  read %d\n    ref: %s\n    new: %s" % (i, r, w))
print("FM_QCHECK_TOTAL files=%d unexpected=%d" % (len(names), total_bad))
PY
echo "--- where the time went, per query file"
for f in "$OUT"/*.stderr; do
    echo "$(basename "$f" .stderr):"
    grep -E '^(Time|Search time|Reads whose)' "$f" | sed 's/^/  /' || true
done
echo "FM_QUERY_CHECK_DONE"
