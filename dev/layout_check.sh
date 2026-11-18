#!/bin/bash -l
# The short/long finimizer threshold (-x of build-fmin) only decides how finimizers are laid
# out in the lookup tables. This job builds the index from one .cfm at each value given, and
# runs every query file against each index and against a reference index, all on one node,
# comparing results with reference outputs. Index sizes, build costs and query times are
# reported side by side. The .cfm, the reference index and the query files are only read.
#
# Settings (environment variables):
#   FM_REPO       path of this repository checkout (required)
#   FM_CFM        the stage 1 output to build from (required)
#   FM_INDEX_REF  reference index prefix, built with the usual settings (required)
#   FM_QUERIES    directory of query files *.fna (required)
#   FM_REF_OUT    directory with reference results, NAME.txt (required)
#   FM_REF_TIME   directory with reference timings (optional)
#   FM_XS         values of -x to try (default "16 14")
#   FM_P          prefix length (default 10)     FM_K   k (default 31)
#   FM_OUT        output directory (default: the current directory)
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
CFM=$(readlink -f "${FM_CFM:?set FM_CFM}")
INDEX_REF=${FM_INDEX_REF:?set FM_INDEX_REF}
OUT=$(readlink -f "${FM_OUT:-$PWD}")
XS=${FM_XS:-16 14}
export FM_REPO FM_QUERIES FM_REF_OUT FM_REF_TIME="${FM_REF_TIME:-}"
: "${FM_QUERIES:?set FM_QUERIES}" "${FM_REF_OUT:?set FM_REF_OUT}"

echo "FM_LAYOUT_CHECK start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$(git -C "$REPO" rev-parse --short HEAD) cfm=$CFM index_ref=$INDEX_REF xs='$XS' p=${FM_P:-10} k=${FM_K:-31} out=$OUT"
mkdir -p "$OUT"

for X in $XS; do
    echo "=== building with -x $X"
    mkdir -p "$OUT/x$X"
    (cd "$OUT/x$X" && FM_CFM=$CFM FM_X=$X FM_OUT=$OUT/x$X bash "$REPO/dev/stage2.sh") > "$OUT/x$X.build.out" 2>&1 \
        || { echo "FM_ERROR build failed for -x $X:"; tail -n 20 "$OUT/x$X.build.out"; exit 1; }
    grep -E '^(wall_s|peak_rss_gb)|n_short_fmin' "$OUT/x$X.build.out" | tr '\n' ' '; echo
done

for X in ref $XS; do
    if [ "$X" = ref ]; then I=$INDEX_REF; else I=$OUT/x$X/idx; fi
    echo "=== queries against the index $X"
    FM_INDEX=$I FM_OUT=$OUT/q$X bash "$REPO/dev/query_check.sh" > "$OUT/q$X.out" 2>&1 \
        || { echo "FM_ERROR query check failed for $X:"; tail -n 20 "$OUT/q$X.out"; exit 1; }
    grep -E '^FM_QCHECK_TOTAL' "$OUT/q$X.out"
done

echo "--- summary"
python3 - "$OUT" "$INDEX_REF" "$XS" <<'PY'
import os, re, sys
out, index_ref, xs = sys.argv[1], sys.argv[2], sys.argv[3].split()
labels = ["ref"] + xs
size = {"ref": os.path.getsize(index_ref + ".fmin")}
for x in xs:
    size[x] = os.path.getsize(os.path.join(out, "x" + x, "idx.fmin"))
rows = {}
for lab in labels:
    for line in open(os.path.join(out, "q" + lab + ".out")):
        if line.startswith("FM_QCHECK ") :
            f = dict(kv.split("=", 1) for kv in line.split()[2:] if "=" in kv)
            name = line.split()[1]
            se = open(os.path.join(out, "q" + lab, name + ".stderr")).read()
            m = re.search(r"index load ([0-9.]+), search ([0-9.]+)", se)
            lk = re.search(r"finimizer lookups ([0-9.]+)", se)
            rows.setdefault(name, {})[lab] = (f.get("identical"), f.get("differ_acgt_only"), f.get("wall_s"), m.group(1) if m else "n/a", m.group(2) if m else "n/a", lk.group(1) if lk else "n/a")
print("index size (bytes): " + "  ".join("%s=%d" % (lab, size[lab]) for lab in labels))
for name in sorted(rows):
    print("FM_LAYOUT %s | " % name + " | ".join(
        "%s: identical=%s differ=%s wall_s=%s load_s=%s search_s=%s lookups_s=%s" % ((lab,) + rows[name][lab]) for lab in labels if lab in rows[name]))
PY
echo "FM_LAYOUT_CHECK_DONE"
