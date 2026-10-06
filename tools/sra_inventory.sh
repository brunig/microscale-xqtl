#!/bin/bash
# Run ON HOFFMAN2 as a batch job:   qsub sra_inventory.sh
# Read-only inventory of every FASTQ to be deposited at SRA: name, size and
# MD5 (SRA checks uploads against the MD5s in the metadata sheet).
# Writes ~/xqtl_sra/fastq_inventory.tsv and ~/xqtl_sra/summary.txt.
#$ -cwd
#$ -o sra_inventory.log
#$ -j y
#$ -l h_rt=6:00:00,h_data=2G
#$ -pe shared 8
set -uo pipefail
OUT=~/xqtl_sra; mkdir -p "$OUT"
POOLS=/u/project/kruglyak/bruni/projects/xQTLpaper/fastq
PARENTS=/u/project/kruglyak/bruni/projects/A1ref/fastq

{ printf 'set\tpath\tbytes\n'
  { find "$POOLS"   -name '*.fastq.gz' -printf 'pools\t%p\t%s\n'
    find "$PARENTS" -name '*.fastq.gz' -printf 'parents\t%p\t%s\n'; } | sort
} > "$OUT/fastq_sizes.tsv"

## MD5s, 8 at a time
tail -n +2 "$OUT/fastq_sizes.tsv" | cut -f2 | tr '\n' '\0' | xargs -0 -P 8 -n 1 md5sum > "$OUT/fastq_md5.txt"

{
  echo "== files per set"; tail -n +2 "$OUT/fastq_sizes.tsv" | cut -f1 | sort | uniq -c
  echo "== total size per set"
  tail -n +2 "$OUT/fastq_sizes.tsv" | awk -F'\t' '{s[$1]+=$3} END {for (k in s) printf "%s %.1f GB\n", k, s[k]/1e9}'
  echo "== lanes and reads (pools)"
  grep -oE '_L00[0-9]_R[12]_' "$OUT/fastq_sizes.tsv" | sort | uniq -c
  echo "== pool files not matching SampleNNN_SNNN_L00N_RN_001"
  awk -F'\t' '$1=="pools"' "$OUT/fastq_sizes.tsv" | cut -f2 | xargs -n1 basename \
    | grep -vE '^Sample[0-9]{3}_S[0-9]+_L00[0-9]_R[12]_001\.fastq\.gz$'
  echo "== samples missing a lane or read (pools)"
  awk -F'\t' '$1=="pools"' "$OUT/fastq_sizes.tsv" | cut -f2 | xargs -n1 basename \
    | grep -oE '^Sample[0-9]{3}' | sort | uniq -c | awk '$1!=4'
  echo "== smallest 5 pool files"; awk -F'\t' '$1=="pools"' "$OUT/fastq_sizes.tsv" | sort -t$'\t' -k3,3n | head -5
  echo "== parent dir layout (first 6)"; awk -F'\t' '$1=="parents"' "$OUT/fastq_sizes.tsv" | head -6
  echo "== reads per file, first two pool samples, both lanes"
  for f in $(awk -F'\t' '$1=="pools"' "$OUT/fastq_sizes.tsv" | cut -f2 | grep -E 'Sample00[12]_' | grep R1); do
    echo "$(basename $f) $(( $(zcat "$f" | wc -l) / 4 ))"; done
  echo "== first header line, lane 1 and lane 2"
  for f in $(awk -F'\t' '$1=="pools"' "$OUT/fastq_sizes.tsv" | cut -f2 | grep Sample001_ | grep R1); do
    zcat "$f" | head -1; done
} > "$OUT/summary.txt" 2>&1
echo done
