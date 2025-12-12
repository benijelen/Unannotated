#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <interproscan.tsv>"
  exit 1
fi

IPRSCAN="$1"

# Strip directory and extension safely
BASENAME=$(basename "$IPRSCAN")
PREFIX="${BASENAME%.tsv}"

PAIR_OUT="${PREFIX}.seq_ipr_pairs.tsv"
IDS_OUT="${PREFIX}.seq_ipr_ids.tsv"

echo "Input:  $IPRSCAN"
echo "Pairs:  $PAIR_OUT"
echo "IDs:    $IDS_OUT"

# Step 1: extract (query_id, IPRxxxxxx) pairs
cut -f1,12 "$IPRSCAN" \
  | awk -F'\t' '$2 ~ /^IPR/ {print $1"\t"$2}' \
  | sort -u \
  > "$PAIR_OUT"

# Step 2: collapse to query_id <tab> IPRa;IPRb;...
awk -F'\t' '
BEGIN { OFS="\t" }
{
  if (NR == 1) {
    prev = $1
    list = $2
    next
  }
  if ($1 != prev) {
    print prev, list
    prev = $1
    list = $2
  } else {
    list = list ";" $2
  }
}
END {
  if (NR > 0) print prev, list
}
' "$PAIR_OUT" > "$IDS_OUT"

echo "Done."
