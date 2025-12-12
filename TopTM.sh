#!/usr/bin/env bash

# Exit on error or unset variable
set -euo pipefail

# Check input
if [ "$#" -ne 1 ]; then
  echo "Usage: bash TopTM.sh inputfile.m8"
  exit 1
fi

INPUT="$1"
OUTPUT="${INPUT%.m8}_bestHit.m8"

awk '{
  q=$1; t=$2; tms=$3; qc=$4; tc=$5; e=$6; b=$7; ql=$8; tl=$9
  key=q
  if (!(key in best) || tms > bestTMS[key]) {
    best[key]    = $0
    bestTMS[key] = tms
  }
}
END {
  for (k in best) print best[k]
}' "$INPUT" \
| sort > "$OUTPUT"
