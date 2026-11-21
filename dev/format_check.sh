#!/bin/bash -l
# Index format 2 stores the short finimizers in a Pino structure instead of format 1's array
# of 32-bit bucket pointers. This job builds both formats from one .cfm, with the same -x, and
# runs every query file against each index twice, in the order 1, 2, 1, 2, all on one node.
# Results are compared with reference outputs; sizes, build costs and query times are
# reported side by side. The .cfm, the query files and the reference outputs are only read.
#
# Settings (environment variables):
#   FM_REPO       path of this repository checkout (required)
#   FM_CFM        the stage 1 output to build from (required)
#   FM_QUERIES    directory of query files *.fna (required)
#   FM_REF_OUT    directory with reference results, NAME.txt (required)
#   FM_X          -x for both indexes (default 16)
#   FM_P          prefix length (default 10)     FM_K   k (default 31)
#   FM_OUT        output directory (default: the current directory)
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
CFM=$(readlink -f "${FM_CFM:?set FM_CFM}")
OUT=$(readlink -f "${FM_OUT:-$PWD}")
X=${FM_X:-16}
export FM_REPO FM_QUERIES FM_REF_OUT FM_REF_TIME=""
: "${FM_QUERIES:?set FM_QUERIES}" "${FM_REF_OUT:?set FM_REF_OUT}"

echo "FM_FORMAT_CHECK start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$(git -C "$REPO" rev-parse --short HEAD) cfm=$CFM x=$X p=${FM_P:-10} k=${FM_K:-31} out=$OUT"
mkdir -p "$OUT"

for F in 1 2; do
    echo "=== building format $F with -x $X"
    mkdir -p "$OUT/f$F"
    (cd "$OUT/f$F" && FM_CFM=$CFM FM_X=$X FM_FORMAT=$F FM_OUT=$OUT/f$F bash "$REPO/dev/stage2.sh") > "$OUT/f$F.build.out" 2>&1 \
        || { echo "FM_ERROR build failed for format $F:"; tail -n 20 "$OUT/f$F.build.out"; exit 1; }
    grep -E '^(wall_s|peak_rss_gb)|n_short_fmin' "$OUT/f$F.build.out" | tr '\n' ' '; echo
done
python3 "$REPO/dev/index_sizes.py" "$OUT/f1/idx" "$OUT/f2/idx"

for RUN in f1.1 f2.1 f1.2 f2.2; do
    F=${RUN%%.*}
    echo "=== queries, run $RUN"
    FM_INDEX=$OUT/$F/idx FM_OUT=$OUT/q$RUN bash "$REPO/dev/query_check.sh" > "$OUT/q$RUN.out" 2>&1 \
        || { echo "FM_ERROR query check failed for $RUN:"; tail -n 20 "$OUT/q$RUN.out"; exit 1; }
    grep -E '^FM_QCHECK_TOTAL' "$OUT/q$RUN.out"
done

echo "--- summary (times in seconds; lookups = finimizer lookups, the part that changed)"
python3 - "$OUT" <<'PY'
import os, re, sys
out = sys.argv[1]
runs = ["f1.1", "f2.1", "f1.2", "f2.2"]
for f in ("f1", "f2"):
    print("index %s: %d bytes" % (f, os.path.getsize(os.path.join(out, f, "idx.fmin"))))
rows = {}
for run in runs:
    for line in open(os.path.join(out, "q" + run + ".out")):
        if line.startswith("FM_QCHECK "):
            name = line.split()[1]
            fields = dict(kv.split("=", 1) for kv in line.split()[2:] if "=" in kv)
            se = open(os.path.join(out, "q" + run, name + ".stderr")).read()
            m = re.search(r"index load ([0-9.]+), search ([0-9.]+)", se)
            lk = re.search(r"finimizer lookups ([0-9.]+)", se)
            rows.setdefault(name, {})[run] = "identical=%s wall=%s load=%s search=%s lookups=%s" % (
                fields.get("identical"), fields.get("wall_s"), m.group(1) if m else "n/a", m.group(2) if m else "n/a", lk.group(1) if lk else "n/a")
for name in sorted(rows):
    print("FM_FORMAT %s | " % name + " | ".join("%s: %s" % (run, rows[name][run]) for run in runs if run in rows[name]))
PY
echo "FM_FORMAT_CHECK_DONE"
