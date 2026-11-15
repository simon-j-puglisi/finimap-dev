#!/bin/bash -l
# Make a collection of N genomes of one species from the AllTheBacteria files on disk:
# pick a seeded random sample of the species' genomes, extract them from the assembly
# archives with tar, and write refs.list. The ATB files are only read.
# Safe to resubmit: a finished collection is skipped.
#
# The species' genomes are those whose species call (column sylph_species_pre_202505 of
# species_calls.tsv.gz) equals FM_SPECIES; the archive of each comes from the ATB file list.
#
# Settings (environment variables):
#   FM_SPECIES  species name as in the species calls, e.g. "Salmonella enterica" (required)
#   FM_OUT      output directory (default: the current directory); gets fa/ and refs.list
#   FM_N        number of genomes wanted (default 10000)
#   FM_SEED     random seed (default 1)
#   FM_ATB      downloaded ATB data (default ~/proj/akc_data/atb)
set -euo pipefail
export LC_ALL=C

SPECIES=${FM_SPECIES:?set FM_SPECIES}
ATB=${FM_ATB:-$HOME/proj/akc_data/atb}
OUT=$(readlink -f "${FM_OUT:-$PWD}")
N=${FM_N:-10000}; SEED=${FM_SEED:-1}
JOBS=${SLURM_CPUS_PER_TASK:-4}
CALLS=$ATB/meta/Aggregated/Latest_2025-05/species_calls.tsv.gz
FLIST=$ATB/lists/File_Lists/file_list.all.latest.tsv.gz

echo "FM_COLLECT start host=$(hostname) date=$(date '+%F %T')"
echo "species=$SPECIES atb=$ATB out=$OUT n=$N seed=$SEED jobs=$JOBS"
for f in "$CALLS" "$FLIST"; do [ -s "$f" ] || { echo "FM_ERROR missing $f"; exit 1; }; done
mkdir -p "$OUT/fa" "$OUT/members"
cd "$OUT"

if [ -s refs.list ] && [ "$(wc -l < refs.list)" -eq "$N" ]; then
    echo "refs.list already has $N genomes, skipping"
else
    t0=$(date +%s)
    zcat "$CALLS" | awk -F'\t' -v sp="$SPECIES" '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
                  if (!("sample_accession" in c) || !("sylph_species_pre_202505" in c)) { print "FM_ERROR unexpected columns in species calls" > "/dev/stderr"; exit 1 }
                  next }
        $c["sylph_species_pre_202505"] == sp { print $c["sample_accession"] }' > samples.txt
    echo "genomes of the species: $(wc -l < samples.txt)"
    [ -s samples.txt ] || { echo "FM_ERROR no samples for '$SPECIES'"; exit 1; }

    # Candidates in random order: a few more than needed, in case some cannot be extracted
    python3 - samples.txt "$N" "$SEED" > candidates.txt <<'PY'
import random, sys
names = sorted(set(l.strip() for l in open(sys.argv[1]) if l.strip()))
n = int(sys.argv[2])
random.Random(int(sys.argv[3])).shuffle(names)
print("\n".join(names[:n + n // 50 + 10]))
PY
    echo "candidates: $(wc -l < candidates.txt)"

    # The archive of each candidate, from the ATB file list: lines "archive name <TAB> sample"
    zcat "$FLIST" | awk -F'\t' '
        FILENAME == ARGV[1] { want[$1] = 1; next }
        FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
                   if (!("sample" in c) || !("tar_xz" in c)) { print "FM_ERROR unexpected columns in the file list" > "/dev/stderr"; exit 1 }
                   next }
        ($c["sample"] in want) { print $c["tar_xz"] "\t" $c["sample"] }' candidates.txt - | sort -u > candidate_archives.tsv
    echo "candidates in the file list: $(cut -f 2 candidate_archives.tsv | sort -u | wc -l), in $(cut -f 1 candidate_archives.tsv | sort -u | wc -l) archives"

    # One list of members per archive; a member is <archive name without .tar.xz>/<sample>.fa
    rm -f members/*
    find "$ATB/asm" -name '*.tar.xz' > all_archives.txt
    awk -F'\t' '
        FILENAME == ARGV[1] { n = split($0, p, "/"); path[p[n]] = $0; next }
        { if (!($1 in path)) { missing[$1] = 1; next }
          base = $1; sub(/\.tar\.xz$/, "", base)
          print base "/" $2 ".fa" > ("members/" $1 ".txt")
          if (!($1 in seen)) { seen[$1] = 1; print path[$1] "\t" "members/" $1 ".txt" > "jobs.tsv" } }
        END { for (a in missing) print "WARN archive not on disk: " a }' all_archives.txt candidate_archives.tsv
    echo "archives to read: $(wc -l < jobs.tsv)"

    # tar stops reading an archive once every listed member has been extracted
    extract_one() {
        tar -xJf "$1" -C fa --strip-components=1 --occurrence=1 -T "$2" 2> "$2.err" \
            || echo "WARN tar reported a problem for $(basename "$1"): $(head -n 1 "$2.err")"
    }
    export -f extract_one
    tr '\t' '\n' < jobs.tsv | xargs -d '\n' -n 2 -P "$JOBS" bash -c 'extract_one "$1" "$2"' _

    # Keep the first N candidates that were extracted; remove the others
    python3 - candidates.txt fa "$N" > refs.list.part <<'PY'
import os, sys
cands = [l.strip() for l in open(sys.argv[1])]
fa, n = os.path.abspath(sys.argv[2]), int(sys.argv[3])
found = [c for c in cands if os.path.exists(os.path.join(fa, c + ".fa"))]
for c in found[n:]:
    os.remove(os.path.join(fa, c + ".fa"))
sys.stderr.write("candidates extracted: %d of %d; kept %d\n" % (len(found), len(cands), min(n, len(found))))
print("\n".join(os.path.join(fa, c + ".fa") for c in sorted(found[:n])))
PY
    mv refs.list.part refs.list
    rm -rf members jobs.tsv all_archives.txt
    echo "collection took $(( $(date +%s) - t0 )) s"
fi

echo "--- result"
echo "genomes in refs.list: $(wc -l < refs.list)"
echo "files in fa/: $(ls -f fa | grep -c '\.fa$')"
echo "size of fa/: $(du -sh fa | cut -f 1)"
echo "first: $(head -n 1 refs.list)"
echo "sequences and bases in the first genome: $(grep -c '>' "$(head -n 1 refs.list)") $(grep -v '>' "$(head -n 1 refs.list)" | tr -d '\n' | wc -c)"
if [ "$(wc -l < refs.list)" -eq "$N" ]; then echo "FM_COLLECT_DONE"; else echo "FM_ERROR fewer than $N genomes found"; exit 1; fi
