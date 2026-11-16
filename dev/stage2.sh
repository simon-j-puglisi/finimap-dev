#!/bin/bash -l
# Stage 2 of the Finimap build (finimap build-fmin), timed, with a time stamp on every
# log line and an optional byte comparison against a reference .fmin file.
# The C++ program is compiled inside the job, on the node that runs it.
# Safe to resubmit: finished steps are skipped.
#
# Settings (environment variables):
#   FM_REPO      path of this repository checkout (required)
#   FM_CFM       the .cfm file from stage 1 (required; only read)
#   FM_OUT       output directory (default: the current directory)
#   FM_REF_FMIN  reference .fmin to compare with (optional)
#   FM_K, FM_P, FM_X   k, prefix length, short/long threshold (defaults 31, 10, 18)
#   FM_CMAKE_MODULE    module providing cmake (default CMake/4.2.1; "none" to skip)
#   FM_GCC_MODULE      compiler module to load first (default none)
#   FM_KEEP      set to 1 to keep the .fmin even when it equals the reference
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
CFM=$(readlink -f "${FM_CFM:?set FM_CFM}")
OUT=$(readlink -f "${FM_OUT:-$PWD}")
K=${FM_K:-31}; P=${FM_P:-10}; X=${FM_X:-18}
REF_FMIN=${FM_REF_FMIN:-}
JOBS=${SLURM_CPUS_PER_TASK:-4}
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)

echo "FM_STAGE2 start host=$(hostname) date=$(date '+%F %T')"
echo "repo=$REPO commit=$COMMIT cfm=$CFM out=$OUT k=$K p=$P x=$X ref_fmin=${REF_FMIN:-none}"
[ -e "$REPO/dev/timeit.py" ] || { echo "FM_ERROR repo not found at $REPO"; exit 1; }
[ -s "$CFM" ] || { echo "FM_ERROR cfm not found: $CFM"; exit 1; }
mkdir -p "$OUT"

. "$REPO/dev/build_finimap.sh"
fm_load_modules
fm_build_finimap
if [ "$FM_REBUILT" = 1 ]; then rm -f "$OUT/stage2.time.txt"; fi

if [ -s "$OUT/stage2.time.txt" ] && grep -q '^exit_code 0$' "$OUT/stage2.time.txt"; then
    echo "stage 2 already done, skipping"
else
    rm -f "$OUT/idx.fmin"
    python3 "$REPO/dev/timeit.py" "$OUT/stage2.time.part" bash -c \
        'set -o pipefail; "$1" build-fmin -i "$2" -o "$3" -k "$4" -p "$5" -x "$6" 2>&1 | python3 "$7" > "$8"' \
        _ "$OUT/finimap" "$CFM" "$OUT/idx" "$K" "$P" "$X" "$REPO/dev/stamp.py" "$OUT/stage2.log" \
        || { echo "FM_ERROR build-fmin failed; end of its log:"; tail -n 20 "$OUT/stage2.log"; exit 1; }
    mv "$OUT/stage2.time.part" "$OUT/stage2.time.txt"
fi

echo "--- time and memory"
cat "$OUT/stage2.time.txt"
echo "--- build-fmin log (seconds elapsed, then the program's line)"
grep -v '^ *[0-9.]* *$' "$OUT/stage2.log" | cut -c 1-160 || true

if [ -e "$OUT/idx.fmin" ]; then
    echo "--- output"
    ls -l "$OUT/idx.fmin"
    if [ -n "$REF_FMIN" ] && [ ! -e "$REF_FMIN" ]; then
        echo "FM_COMPARE skipped: reference not found: $REF_FMIN"
    elif [ -n "$REF_FMIN" ]; then
        ls -l "$REF_FMIN"
        if cmp -s "$OUT/idx.fmin" "$REF_FMIN"; then
            echo "FM_COMPARE IDENTICAL to reference"
            if [ "${FM_KEEP:-0}" != 1 ]; then rm -f "$OUT/idx.fmin"; echo "removed the identical copy"; fi
        else
            echo "FM_COMPARE DIFFERENT from reference"
            cmp "$OUT/idx.fmin" "$REF_FMIN" | head -n 1 || true
        fi
    fi
else
    echo "output already removed after an earlier identical comparison"
fi
echo "FM_STAGE2_DONE"
