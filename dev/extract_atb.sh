#!/bin/bash -l
# Extract a seeded random sample of genomes from AllTheBacteria assembly archives
# (tar.xz files with one SAMPLE.fa per genome) and write a list of their paths.
# The archives are only read. Safe to resubmit: a finished extraction is skipped.
#
# Settings (environment variables):
#   FM_REPO      path of this repository checkout (required)
#   FM_ARCHIVES  file listing the archives, one path per line (required)
#   FM_SAMPLES   file of sample names to choose from, one per line (required)
#   FM_OUT       output directory (default: the current directory); gets fa/ and refs.list
#   FM_N         number of genomes wanted (default 10000)
#   FM_SEED      random seed (default 1)
set -euo pipefail
export LC_ALL=C

REPO=${FM_REPO:?set FM_REPO}
ARCHIVES=$(readlink -f "${FM_ARCHIVES:?set FM_ARCHIVES}")
SAMPLES=$(readlink -f "${FM_SAMPLES:?set FM_SAMPLES}")
OUT=$(readlink -f "${FM_OUT:-$PWD}")
N=${FM_N:-10000}; SEED=${FM_SEED:-1}
JOBS=${SLURM_CPUS_PER_TASK:-4}

echo "FM_EXTRACT start host=$(hostname) date=$(date '+%F %T')"
echo "archives=$ARCHIVES ($(wc -l < "$ARCHIVES")) samples=$SAMPLES ($(wc -l < "$SAMPLES")) out=$OUT n=$N seed=$SEED jobs=$JOBS"
mkdir -p "$OUT/fa"
cd "$OUT"

if [ -s refs.list ] && [ "$(wc -l < refs.list)" -eq "$N" ]; then
    echo "refs.list already has $N genomes, skipping"
else
    t0=$(date +%s)
    # Candidates in random order: a few more than needed, in case some are not in the archives
    python3 - "$SAMPLES" "$N" "$SEED" > candidates.txt <<'PY'
import random, sys
names = sorted(set(l.strip() for l in open(sys.argv[1]) if l.strip()))
n = int(sys.argv[2])
random.Random(int(sys.argv[3])).shuffle(names)
print("\n".join(names[:n + n // 50 + 10]))
PY
    echo "candidates: $(wc -l < candidates.txt)"

    # One streaming pass over each archive, several archives at a time
    cat > extract_one.py <<'PY'
import os, sys, tarfile
archive, wanted_file, outdir = sys.argv[1:4]
wanted = set(l.strip() for l in open(wanted_file))
n = 0
with tarfile.open(archive, "r|xz") as tar:
    for member in tar:
        name = os.path.basename(member.name)
        if not member.isfile() or not name.endswith(".fa") or name[:-3] not in wanted:
            continue
        dest = os.path.join(outdir, name)
        if not os.path.exists(dest):
            src = tar.extractfile(member)
            with open(dest + ".part", "wb") as f:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            os.rename(dest + ".part", dest)
        n += 1
print("%s: %d genomes" % (os.path.basename(archive), n))
PY
    xargs -a "$ARCHIVES" -P "$JOBS" -I{} python3 extract_one.py {} candidates.txt fa
    rm -f fa/*.part

    # Keep the first N candidates that were found; remove the others
    python3 - candidates.txt fa "$N" > refs.list.part <<'PY'
import os, sys
cands = [l.strip() for l in open(sys.argv[1])]
fa, n = os.path.abspath(sys.argv[2]), int(sys.argv[3])
found = [c for c in cands if os.path.exists(os.path.join(fa, c + ".fa"))]
for c in found[n:]:
    os.remove(os.path.join(fa, c + ".fa"))
sys.stderr.write("candidates found in the archives: %d of %d; kept %d\n" % (len(found), len(cands), min(n, len(found))))
print("\n".join(os.path.join(fa, c + ".fa") for c in sorted(found[:n])))
PY
    mv refs.list.part refs.list
    rm -f extract_one.py
    echo "extraction took $(( $(date +%s) - t0 )) s"
fi

echo "--- result"
echo "genomes in refs.list: $(wc -l < refs.list)"
echo "files in fa/: $(ls -f fa | grep -c '\.fa$')"
echo "size of fa/: $(du -sh fa | cut -f 1)"
echo "first: $(head -n 1 refs.list)"
echo "sequences and bases in the first genome: $(grep -c '>' "$(head -n 1 refs.list)") $(grep -v '>' "$(head -n 1 refs.list)" | tr -d '\n' | wc -c)"
if [ "$(wc -l < refs.list)" -eq "$N" ]; then echo "FM_EXTRACT_DONE"; else echo "FM_ERROR fewer than $N genomes found"; exit 1; fi
