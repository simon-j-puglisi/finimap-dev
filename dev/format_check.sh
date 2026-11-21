#!/bin/bash -l
# Compares index formats (--index-format of build-fmin). This job builds each format given from
# one .cfm, with the same -x, and runs every query file against each index twice, the formats
# taking turns (for formats "2 3": 2, 3, 2, 3), all on one node. Results are compared with
# reference outputs; sizes, build costs and query times are reported side by side. The .cfm,
# the query files and the reference outputs are only read.
#
# Settings (environment variables):
#   FM_REPO       path of this repository checkout (required)
#   FM_CFM        the stage 1 output to build from (required)
#   FM_QUERIES    directory of query files *.fna (required)
#   FM_REF_OUT    directory with reference results, NAME.txt (required)
#   FM_FORMATS    formats to compare (default "1 2")
#   FM_X          -x for every index (default 16)
#   FM_P          prefix length (default 10)     FM_K   k (default 31)
#   FM_OUT        output directory (default: the current directory)
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
CFM=$(readlink -f "${FM_CFM:?set FM_CFM}")
OUT=$(readlink -f "${FM_OUT:-$PWD}")
X=${FM_X:-16}
FORMATS=${FM_FORMATS:-1 2}
export FM_REPO FM_QUERIES FM_REF_OUT FM_REF_TIME=""
: "${FM_QUERIES:?set FM_QUERIES}" "${FM_REF_OUT:?set FM_REF_OUT}"

echo "FM_FORMAT_CHECK start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$(git -C "$REPO" rev-parse --short HEAD) cfm=$CFM formats='$FORMATS' x=$X p=${FM_P:-10} k=${FM_K:-31} out=$OUT"
mkdir -p "$OUT"

INDEXES=""
RUNS=""
for F in $FORMATS; do
    INDEXES="$INDEXES $OUT/f$F/idx"
done
for R in 1 2; do
    for F in $FORMATS; do
        RUNS="$RUNS f$F.$R"
    done
done

for F in $FORMATS; do
    echo "=== building format $F with -x $X"
    mkdir -p "$OUT/f$F"
    (cd "$OUT/f$F" && FM_CFM=$CFM FM_X=$X FM_FORMAT=$F FM_OUT=$OUT/f$F bash "$REPO/dev/stage2.sh") > "$OUT/f$F.build.out" 2>&1 \
        || { echo "FM_ERROR build failed for format $F:"; tail -n 20 "$OUT/f$F.build.out"; exit 1; }
    grep -E '^(wall_s|peak_rss_gb)|n_short_fmin' "$OUT/f$F.build.out" | tr '\n' ' '; echo
done
python3 "$REPO/dev/index_sizes.py" $INDEXES

for RUN in $RUNS; do
    F=${RUN%%.*}
    echo "=== queries, run $RUN"
    FM_INDEX=$OUT/$F/idx FM_OUT=$OUT/q$RUN bash "$REPO/dev/query_check.sh" > "$OUT/q$RUN.out" 2>&1 \
        || { echo "FM_ERROR query check failed for $RUN:"; tail -n 20 "$OUT/q$RUN.out"; exit 1; }
    grep -E '^FM_QCHECK_TOTAL' "$OUT/q$RUN.out"
done

echo "--- summary (times in seconds)"
python3 - "$OUT" "$FORMATS" <<'PY'
import os, re, sys
out = sys.argv[1]
formats = sys.argv[2].split()
runs = ["f%s.%d" % (f, r) for r in (1, 2) for f in formats]
for f in ["f" + x for x in formats]:
    print("index %s: %d bytes" % (f, os.path.getsize(os.path.join(out, f, "idx.fmin"))))
rows = {}
for run in runs:
    for line in open(os.path.join(out, "q" + run + ".out")):
        if line.startswith("FM_QCHECK "):
            name = line.split()[1]
            fields = dict(kv.split("=", 1) for kv in line.split()[2:] if "=" in kv)
            se = open(os.path.join(out, "q" + run, name + ".stderr")).read()
            m = re.search(r"index load ([0-9.]+), search ([0-9.]+)", se)
            part = lambda p: (re.search(p + r" ([0-9.]+)", se) or [None, "n/a"])[1]
            rows.setdefault(name, {})[run] = "identical=%s wall=%s load=%s search=%s lookups=%s pick=%s count=%s" % (
                fields.get("identical"), fields.get("wall_s"), m.group(1) if m else "n/a", m.group(2) if m else "n/a",
                part("finimizer lookups"), part("picking finimizers"), part("counting colors"))
for name in sorted(rows):
    print("FM_FORMAT %s | " % name + " | ".join("%s: %s" % (run, rows[name][run]) for run in runs if run in rows[name]))
PY
echo "FM_FORMAT_CHECK_DONE"
