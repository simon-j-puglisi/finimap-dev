# Shared by the job scripts: load the build modules and compile finimap into the output
# directory, on the node that will run it. Source this file, then call fm_build_finimap.
#
# Uses: REPO (repository checkout), OUT (output directory), COMMIT (short commit of REPO),
#       JOBS (parallel compile jobs)
# Settings: FM_CMAKE_MODULE (default CMake/4.2.1; "none" to skip), FM_GCC_MODULE (default none)
# Result: $OUT/finimap, rebuilt only when COMMIT differs from $OUT/finimap.commit.
#         FM_REBUILT is 1 if it was rebuilt, else 0.

fm_load_modules() {
    if ! type module > /dev/null 2>&1; then
        for f in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
            if [ -r "$f" ]; then set +u; . "$f"; set -u; break; fi
        done
    fi
    if type module > /dev/null 2>&1; then
        set +u
        [ "${FM_GCC_MODULE:-none}" = none ] || module load "$FM_GCC_MODULE" || echo "WARN: module load $FM_GCC_MODULE failed"
        [ "${FM_CMAKE_MODULE:-CMake/4.2.1}" = none ] || module load "${FM_CMAKE_MODULE:-CMake/4.2.1}" || echo "WARN: module load failed"
        set -u
    else
        echo "WARN: no module command; using whatever is on PATH"
    fi
    echo "compiler: $(g++ --version | head -n 1)"
}

fm_build_finimap() {
    FM_REBUILT=0
    if [ "$(cat "$OUT/finimap.commit" 2>/dev/null)" = "$COMMIT" ] && [ -x "$OUT/finimap" ]; then
        echo "finimap already built at $COMMIT"
        return 0
    fi
    local t0
    t0=$(date +%s)
    rm -rf "$OUT/sdsl-build" "$OUT/finimap" "$OUT/src"
    mkdir -p "$OUT/sdsl-build"
    (cd "$OUT/sdsl-build" && cmake "$REPO/sdsl-lite" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_CXX_FLAGS="-std=c++2a -O2 -march=native -Wno-template-body" && make -j "$JOBS") > "$OUT/build.log" 2>&1 \
        || { echo "FM_ERROR sdsl-lite build failed:"; tail -n 20 "$OUT/build.log"; exit 1; }
    # Compile from a copy of the sources, so that the repository checkout is not written to
    mkdir -p "$OUT/src"
    cp -r "$REPO/Makefile" "$REPO/src" "$REPO/include" "$OUT/src/"
    make -C "$OUT/src" finimap CXX="g++ -Wno-template-body" SDSL_INCLUDE="$REPO/sdsl-lite/include/" SDSL_LIB="$OUT/sdsl-build/lib/" >> "$OUT/build.log" 2>&1 \
        || { echo "FM_ERROR finimap build failed:"; grep -E 'error' -A 3 "$OUT/build.log" | head -n 30; exit 1; }
    mv "$OUT/src/finimap" "$OUT/finimap"
    rm -rf "$OUT/src" "$OUT/sdsl-build"
    echo "$COMMIT" > "$OUT/finimap.commit"
    FM_REBUILT=1
    echo "built finimap at $COMMIT in $(( $(date +%s) - t0 )) s"
}
