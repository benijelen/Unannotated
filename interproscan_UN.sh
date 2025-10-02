#!/bin/bash

indir="/home/yggshare/current_projects/GB4_AT50-22/bins/interproscan_5.75-106.0_2025July"
faadir="/home/yggshare/ALL_DATA/02.bin_data/04.AT50-22/2025-08-15_all_quality_bins/faa"
outdir="/home/bj22232/InterproUnannotated"
outfile="$outdir/interpro_coverage.csv"

echo "bin,total_proteins,with_ipr,with_ipr_pct,no_ipr,no_ipr_pct" > "$outfile"

for tsv in "$indir"/*.tsv; do
    base=$(basename "$tsv" .tsv)
    faa="${faadir}/${base}"

    if [ ! -f "$faa" ]; then
        echo "Warning: no FASTA found for $base at $faa" >&2
        continue
    fi

    total=$(grep -c '^>' "$faa")
    with_ipr=$(awk -F'\t' '$12 ~ /^IPR/' "$tsv" | cut -f1 | sort -u | wc -l)
    no_ipr=$(( total - with_ipr ))

    with_pct=$(awk -v w=$with_ipr -v t=$total 'BEGIN{if(t>0) printf "%.2f", 100*w/t; else print "0.00"}')
    no_pct=$(awk -v n=$no_ipr  -v t=$total 'BEGIN{if(t>0) printf "%.2f", 100*n/t; else print "0.00"}')

    echo "$base,$total,$with_ipr,$with_pct,$no_ipr,$no_pct" >> "$outfile"
done

echo "Done. Results written to $outfile"
