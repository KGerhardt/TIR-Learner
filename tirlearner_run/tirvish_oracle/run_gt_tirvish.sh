#!/usr/bin/env bash
# Run gt suffixerator + gt tirvish on each oracle chunk with the EXACT flags
# TIR-Learner's tirvish_new.one_tirvish uses. Captures raw GFF3 stdout (the gold
# set) and per-chunk wall time. Idempotent: skips a chunk whose .gff already exists.
set -euo pipefail

GT="${GT:-/mnt/c/Users/kenji/Desktop/conda_envs/tirlearner/bin/gt}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p gold idx

for fa in chunks/chunk*.fa; do
  base="$(basename "$fa" .fa)"
  gff="gold/${base}.tirvish.gff"
  if [[ -f "$gff" ]]; then echo "skip $base (gff exists)"; continue; fi
  idx="idx/${base}"

  # suffixerator: build the mirrored enhanced suffix array (exact pipeline flags)
  "$GT" suffixerator -db "$fa" -indexname "$idx" \
        -tis -suf -lcp -des -ssp -sds -dna -mirrored >/dev/null 2>&1

  # tirvish: exact pipeline params; -seqids yes so GFF carries the `name;;0` ids
  t0=$(date +%s.%N)
  "$GT" tirvish -index "$idx" -seed 20 -mintirlen 10 -maxtirlen 1000 \
        -mintirdist 10 -maxtirdist 5000 -similar 80 -mintsd 2 -maxtsd 11 \
        -vic 13 -seqids "yes" > "$gff" 2> "gold/${base}.tirvish.err"
  t1=$(date +%s.%N)

  # clean the suffix-array files (as one_tirvish does)
  rm -f "${idx}".{des,esq,lcp,llv,md5,prj,sds,suf,ssp,al1,ois}
  n=$(grep -vc '^#' "$gff" || true)
  printf '%s\twall_s=%.1f\tgff_feature_lines=%s\n' "$base" "$(echo "$t1 - $t0" | bc)" "$n"
done
