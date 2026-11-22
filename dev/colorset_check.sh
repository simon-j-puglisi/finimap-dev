#!/bin/bash -l
# Estimates how much smaller the color sets of existing Finimap indexes would be under
# cross-set encodings (dev/colorsets/cs_analyze.cpp), and how far a general-purpose compressor
# gets on them. The indexes are only read.
#
# Settings (environment variables):
#   FM_REPO      path of this repository checkout (required)
#   FM_INDEXES   index prefixes to analyse, separated by spaces (required; any index format)
#   FM_OUT       output directory (default: the current directory)
#   FM_COMPRESS  1 to also compress each index's sets with zstd or xz (default 1)
#   FM_CS_ARGS   more options for cs_analyze, such as --no-ref (default none)
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
OUT=$(readlink -f "${FM_OUT:-$PWD}")
INDEXES=${FM_INDEXES:?set FM_INDEXES}
COMPRESS=${FM_COMPRESS:-1}
CS_ARGS=${FM_CS_ARGS:-}
JOBS=${SLURM_CPUS_PER_TASK:-4}
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
. "$REPO/dev/build_finimap.sh"

echo "FM_CS_CHECK start host=$(hostname) date=$(date '+%F %T') commit=$COMMIT threads=$JOBS args='$CS_ARGS' compress=$COMPRESS"
mkdir -p "$OUT"
fm_load_modules

# Build sdsl-lite and the analysis program in the output directory, from a copy of the sources
t0=$(date +%s)
rm -rf "$OUT/sdsl-build" "$OUT/src"
mkdir -p "$OUT/sdsl-build" "$OUT/src"
(cd "$OUT/sdsl-build" && cmake "$REPO/sdsl-lite" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_CXX_FLAGS="-std=c++2a -O2 -march=native -Wno-template-body" && make -j "$JOBS") > "$OUT/build.log" 2>&1 \
    || { echo "FM_ERROR sdsl-lite build failed:"; tail -n 20 "$OUT/build.log"; exit 1; }
cp -r "$REPO/src" "$REPO/include" "$REPO/dev/colorsets" "$OUT/src/"
(cd "$OUT/src" && g++ -std=c++2a -O3 -march=native -fopenmp -DNDEBUG -Wno-deprecated-declarations -Wno-template-body \
    -I include -I "$REPO/sdsl-lite/include" -o "$OUT/cs_analyze" colorsets/cs_analyze.cpp src/utilities.cpp src/seqio.cpp \
    -L "$OUT/sdsl-build/lib" -lsdsl -lz) >> "$OUT/build.log" 2>&1 \
    || { echo "FM_ERROR cs_analyze build failed:"; grep -E 'error' -A 3 "$OUT/build.log" | head -n 30; exit 1; }
rm -rf "$OUT/src" "$OUT/sdsl-build"
echo "built cs_analyze in $(( $(date +%s) - t0 )) s"

COMPRESSOR=""
if [ "$COMPRESS" = 1 ]; then
    if command -v zstd > /dev/null; then COMPRESSOR=zstd; elif command -v xz > /dev/null; then COMPRESSOR=xz; fi
    echo "compressor: ${COMPRESSOR:-none found}"
fi

for IDX in $INDEXES; do
    NAME=$(basename "$(dirname "$IDX")")_$(basename "$IDX")
    echo "=== $IDX ($NAME)"
    if [ ! -f "$IDX.fmin" ]; then
        echo "FM_ERROR $IDX.fmin not found; skipped"
        continue
    fi
    DUMP_ARGS=()
    [ -n "$COMPRESSOR" ] && DUMP_ARGS=(--dump "$OUT/$NAME.sets.bin")
    python3 "$REPO/dev/timeit.py" "$OUT/$NAME.time" "$OUT/cs_analyze" "$IDX" --threads "$JOBS" ${DUMP_ARGS[@]+"${DUMP_ARGS[@]}"} $CS_ARGS \
        > "$OUT/$NAME.out" 2> "$OUT/$NAME.err" \
        || { echo "FM_ERROR cs_analyze failed on $IDX:"; tail -n 20 "$OUT/$NAME.err"; continue; }
    cat "$OUT/$NAME.out"
    grep -E '^(wall_s|peak_rss_gb)' "$OUT/$NAME.time" | tr '\n' ' '; echo
    if [ -n "$COMPRESSOR" ] && [ -f "$OUT/$NAME.sets.bin" ]; then
        t1=$(date +%s)
        RAW=$(stat -c %s "$OUT/$NAME.sets.bin")
        if [ "$COMPRESSOR" = zstd ]; then
            zstd -q -19 --long=31 -T"$JOBS" -c "$OUT/$NAME.sets.bin" > "$OUT/$NAME.sets.zst"
            PACKED=$(stat -c %s "$OUT/$NAME.sets.zst")
            rm -f "$OUT/$NAME.sets.zst"
        else
            xz -q -9 -T"$JOBS" -c "$OUT/$NAME.sets.bin" > "$OUT/$NAME.sets.xz"
            PACKED=$(stat -c %s "$OUT/$NAME.sets.xz")
            rm -f "$OUT/$NAME.sets.xz"
        fi
        rm -f "$OUT/$NAME.sets.bin"
        python3 -c "import sys; r, p = int(sys.argv[1]), int(sys.argv[2]); print('CS_GENERIC %s: all sets as bitmaps %.1f MB -> %.1f MB (%.1fx), %s s' % (sys.argv[3], r / 1e6, p / 1e6, r / p, sys.argv[4]))" \
            "$RAW" "$PACKED" "$COMPRESSOR" "$(( $(date +%s) - t1 ))"
    fi
done
echo "FM_CS_CHECK_DONE"
