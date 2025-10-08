#!/bin/bash
# ----------------------------------------------------------
# Extract one ESMFold-ready sequence from a single InterProScan/FASTA pair
# ----------------------------------------------------------

# --- user inputs ---
bin="GB4-5258-C15-3_Bin578.1L.faa"   # 👈 change this for another bin
indir="/home/yggshare/current_projects/GB4_AT50-22/foundation/bins/interproscan_5.75-106.0_2025July"
faadir="/home/yggshare/ALL_DATA/02.bin_data/04.AT50-22/2025-08-15_all_quality_bins/faa"
outdir="/home/bj22232/InterproUnannotated"

tsv="$indir/$bin.tsv"
faa="$faadir/$bin"

mkdir -p "$outdir"

if [ ! -f "$tsv" ] || [ ! -f "$faa" ]; then
    echo "❌ Missing input files for $bin"
    exit 1
fi

echo "🔹 Extracting no-IPR sequences for $bin ..."

# 1. IDs with IPRs
awk -F'\t' '$12 ~ /^IPR/' "$tsv" | cut -f1 | sort -u > "$outdir/with_ipr.ids"

# 2. All IDs
grep '^>' "$faa" | sed 's/^>//' > "$outdir/all.ids"

# 3. IDs without IPRs
grep -F -v -f "$outdir/with_ipr.ids" "$outdir/all.ids" > "$outdir/no_ipr.ids"

# 4. Extract no-IPR sequences
awk 'NR==FNR {keep[$1]; next}
     /^>/ {header=$0; seq=""; next}
     {seq=seq$0}
     /^>/ || ENDFILE {
         if (header) {
             id=substr(header,2);
             gsub(/\r/,"",id);
             if (id in keep) print header"\n"seq
         }
     }' "$outdir/no_ipr.ids" "$faa" > "$outdir/${bin%.faa}_noIPR.faa"

# 5. Filter 120–450 aa
awk 'BEGIN{RS=">"; ORS=""}
NR>1{
    split($0,a,"\n"); header=a[1]; seq="";
    for(i=2;i<=length(a);i++) seq=seq a[i];
    gsub(/\r/,"",seq);
    len=length(seq);
    if(len>=120 && len<=450)
        print ">"header"\n"seq"\n"
}' "$outdir/${bin%.faa}_noIPR.faa" > "$outdir/${bin%.faa}_noIPR_120to450.faa"

# 6. Pick one random sequence
grep -n "^>" "$outdir/${bin%.faa}_noIPR_120to450.faa" | shuf -n1 | cut -d: -f1 |
while read -r start; do
    end=$(grep -n "^>" "$outdir/${bin%.faa}_noIPR_120to450.faa" | awk -F: -v s="$start" '$1>s {print $1; exit}')
    if [ -z "$end" ]; then
        tail -n +"$start" "$outdir/${bin%.faa}_noIPR_120to450.faa"
    else
        sed -n "${start},$(($end-1))p" "$outdir/${bin%.faa}_noIPR_120to450.faa"
    fi
done > "$outdir/ESM_candidate_1.faa"

# 7. Quick stats
num_total=$(grep -c "^>" "$faa")
num_noipr=$(grep -c "^>" "$outdir/${bin%.faa}_noIPR.faa")
num_filt=$(grep -c "^>" "$outdir/${bin%.faa}_noIPR_120to450.faa")

echo "✅ $num_noipr / $num_total proteins lack IPRs"
echo "✅ $num_filt are in 120–450 aa range"
echo "✅ Candidate saved to: $outdir/ESM_candidate_1.faa"
