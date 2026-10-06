#!/bin/bash
# Run ON HOFFMAN2. Collects the code files the repo is missing, and sizes the
# data for the Zenodo / SRA deposit. Writes everything to ~/xqtl_v7_collect/.
set -uo pipefail
P=/u/project/kruglyak/bruni/projects/xQTLpaper
OUT=~/xqtl_v7_collect; mkdir -p "$OUT"; cd "$P"

## 1. code: the exact files v7 used, plus the Nextflow front end
{
  printf '%s\n' files/qc_bam_generate.sh files/qc_pipeline.R files/hamming_export.py \
    files/hamming_integrate.R files/round3_apply.R files/verify_parents.R \
    files/flatkey_A1.csv runs/v7/scripts/xQTL_analysis7.R
  find . -maxdepth 3 \( -path ./work -o -path ./runs -o -path ./fastq -o -name .ipynb_checkpoints \) -prune -o \
    -type f \( -name '*.nf' -o -name 'nextflow.config' -o -name '*.config' \) -print
  find . -maxdepth 3 -path ./work -prune -o -type d \( -name bin -o -name modules \) -print
} | sed 's|^\./||' | sort -u > "$OUT/code_list.txt"
tar czf "$OUT/xqtl_v7_code.tgz" -T "$OUT/code_list.txt"
tar tzf "$OUT/xqtl_v7_code.tgz" | grep -v '/$' | xargs md5sum > "$OUT/code_md5.txt"

## 2. data sizes, for planning the deposit
{
  echo "== top-level directories";  du -sh -- */ 2>/dev/null | sort -h
  echo "== runs/v7";                du -sh runs/v7/* 2>/dev/null | sort -h
  echo "== runs/v6 (min_med_depth = 5 comparator)"; du -sh runs/v6/* 2>/dev/null | sort -h
  echo "== fastq: files, total";    ls fastq | wc -l; du -sh fastq
  echo "== fastq: first names";     ls fastq | head -8
  echo "== per-sample count tables"; find . -maxdepth 3 -path ./work -prune -o -type d -iname '*count*' -print \
                                       | while read d; do du -sh "$d"; ls "$d" | wc -l; done
} > "$OUT/data_sizes.txt" 2>&1


## 3. which driver version produced which stage (see RATIONALE.md section 2)
ls -la --time-style=full-iso runs/v7/scripts/ runs/v7/*.log \
  runs/v7/landscape/ runs/v7/candidates/*/ runs/v7/go/ runs/v7/contrasts/*.csv \
  > "$OUT/v7_timestamps.txt" 2>&1
echo "done: $OUT"; ls -la "$OUT"
