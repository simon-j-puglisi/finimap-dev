#!/bin/bash -l
# Stage 1 of the Finimap build (finimizer_matrix build), timed, with an optional
# byte comparison against a reference .cfm file. Safe to resubmit: a finished
# stage is skipped.
#
# Settings (environment variables):
#   FM_REPO     path of this repository checkout (required under sbatch)
#   FM_REFS     file of genome paths, one per line (required)
#   FM_OUT      output directory (default: the current directory)
#   FM_REF_CFM  reference .cfm to compare with (optional)
#   FM_T        threads (default 32)      FM_K     k (default 31)
#   FM_SBWT_GB  SBWT build memory budget in GB (default 64)
#   FM_KEEP     set to 1 to keep the .cfm even when it equals the reference
set -euo pipefail
export LC_ALL=C

REPO=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
# Under sbatch the script is copied, so take the repo from FM_REPO if given
REPO=${FM_REPO:-$REPO}
BIN=$REPO/finimizer_matrix/target/release/finimizer_matrix
REFS=$(readlink -f "${FM_REFS:?set FM_REFS}")
OUT=$(readlink -f "${FM_OUT:-$PWD}")
T=${FM_T:-32}; K=${FM_K:-31}; SBWT_GB=${FM_SBWT_GB:-64}
REF_CFM=${FM_REF_CFM:-}

echo "FM_STAGE1 start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$(git -C "$REPO" rev-parse --short HEAD) submodule=$(git -C "$REPO/finimizer_matrix" rev-parse --short HEAD)"
echo "refs=$REFS genomes=$(wc -l < "$REFS") out=$OUT threads=$T k=$K sbwt_gb=$SBWT_GB ref_cfm=${REF_CFM:-none}"
[ -e "$REPO/dev/timeit.py" ] || { echo "FM_ERROR repo not found at $REPO: set FM_REPO"; exit 1; }
[ -x "$BIN" ] || { echo "FM_ERROR binary not built: $BIN"; exit 1; }

# Genome paths may be relative to the directory of the list
if [ ! -e "$(head -n 1 "$REFS")" ]; then cd "$(dirname "$REFS")"; fi
[ -e "$(head -n 1 "$REFS")" ] || { echo "FM_ERROR first genome not found: $(head -n 1 "$REFS")"; exit 1; }

mkdir -p "$OUT/tmp"
if [ -s "$OUT/stage1.time.txt" ] && grep -q '^exit_code 0$' "$OUT/stage1.time.txt"; then
    echo "stage 1 already done, skipping the build"
else
    rm -f "$OUT/idx.cfm" "$OUT/idx.cfm.part"
    python3 "$REPO/dev/timeit.py" "$OUT/stage1.time.part" \
        "$BIN" build -i "$REFS" -o "$OUT/idx.cfm.part" -k "$K" -t "$T" -d "$OUT/tmp" -m "$SBWT_GB" --reverse \
        2> "$OUT/stage1.log"
    mv "$OUT/idx.cfm.part" "$OUT/idx.cfm"
    mv "$OUT/stage1.time.part" "$OUT/stage1.time.txt"
fi
rm -rf "$OUT/tmp"

echo "--- time and memory"
cat "$OUT/stage1.time.txt"
echo "--- build log (without per-file lines)"
tr '\r' '\n' < "$OUT/stage1.log" | grep ' INFO ' | grep -v 'Processing file' || true

if [ -e "$OUT/idx.cfm" ]; then
    echo "--- output"
    ls -l "$OUT/idx.cfm"
    if [ -n "$REF_CFM" ] && [ ! -e "$REF_CFM" ]; then
        echo "FM_COMPARE skipped: reference not found: $REF_CFM"
    elif [ -n "$REF_CFM" ]; then
        ls -l "$REF_CFM"
        if cmp -s "$OUT/idx.cfm" "$REF_CFM"; then
            echo "FM_COMPARE IDENTICAL to reference"
            if [ "${FM_KEEP:-0}" != 1 ]; then rm -f "$OUT/idx.cfm"; echo "removed the identical copy"; fi
        else
            echo "FM_COMPARE DIFFERENT from reference"
            cmp "$OUT/idx.cfm" "$REF_CFM" | head -n 1 || true
        fi
    fi
else
    echo "output already removed after an earlier identical comparison"
fi
echo "FM_STAGE1_DONE"
