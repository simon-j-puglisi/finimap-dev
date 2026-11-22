#!/bin/bash -l
# Runs dev/format_check.sh for several species of the benchmark, one after another: for each,
# both formats are built from the benchmark's stage 1 output and queried with its query files,
# and the results are compared with its reference outputs. Species whose files are missing are
# skipped. The benchmark's files are only read.
#
# Settings (environment variables):
#   FM_REPO      path of this repository checkout (required)
#   FM_BENCH     the benchmark's run directory (required), holding for each species SP:
#                  SP/index/finimap/build/idx.cfm, SP/queries/*.fna, SP/out/finimap/*.txt
#   FM_SPECIES   species to run, separated by spaces (required)
#   FM_FORMATS   formats to compare (default "3 4")
#   FM_KEEP      1 to keep the indexes built (default: they are deleted after each species)
#   FM_OUT       output directory (default: the current directory)
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
BENCH=${FM_BENCH:?set FM_BENCH}
SPECIES=${FM_SPECIES:?set FM_SPECIES}
OUT=$(readlink -f "${FM_OUT:-$PWD}")
export FM_FORMATS=${FM_FORMATS:-3 4}
KEEP=${FM_KEEP:-0}

echo "FM_SPECIES_CHECK start host=$(hostname) date=$(date '+%F %T') commit=$(git -C "$REPO" rev-parse --short HEAD) species='$SPECIES' formats='$FM_FORMATS'"
for SP in $SPECIES; do
    CFM=$BENCH/$SP/index/finimap/build/idx.cfm
    QUERIES=$BENCH/$SP/queries
    REF=$BENCH/$SP/out/finimap
    echo "##### $SP"
    if [ ! -f "$CFM" ] || [ ! -d "$QUERIES" ] || [ ! -d "$REF" ]; then
        echo "FM_SKIP $SP: missing $( [ -f "$CFM" ] || echo "$CFM " )$( [ -d "$QUERIES" ] || echo "$QUERIES " )$( [ -d "$REF" ] || echo "$REF" )"
        continue
    fi
    mkdir -p "$OUT/$SP"
    if FM_REPO=$REPO FM_CFM=$CFM FM_QUERIES=$QUERIES FM_REF_OUT=$REF FM_OUT=$OUT/$SP \
        bash "$REPO/dev/format_check.sh" > "$OUT/$SP/format_check.out" 2>&1; then
        sed -n '/^=== building/,$p' "$OUT/$SP/format_check.out" | grep -v '^=== queries'
    else
        echo "FM_ERROR format check failed for $SP:"
        tail -n 20 "$OUT/$SP/format_check.out"
    fi
    if [ "$KEEP" != 1 ]; then
        rm -f "$OUT/$SP"/f*/idx.fmin
    fi
done
echo "FM_SPECIES_CHECK_DONE"
