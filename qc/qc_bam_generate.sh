#!/bin/bash
# ══════════════════════════════════════════════════════════════════
# qc_bam_generate.sh - run samtools flagstat + stats on all BAMs
#
# Outputs per BAM (in qc/bam_stats/):
#   Sample001.flagstat    - short summary (mapping %, dup %, etc.)
#   Sample001.stats.txt   - full samtools stats (SN, IS, GC, etc.)
#
# Idempotent: skips samples whose outputs already exist.
#
# Usage:
#   ./qc_bam_generate.sh                     # use defaults
#   ./qc_bam_generate.sh <bam_dir> <out_dir> <ncores>
#
# On Hoffman2 (interactive):
#   qrsh -l h_data=8G,h_rt=1:00:00 -pe shared 8
#   ./qc_bam_generate.sh
#
# Or submit as SGE job (uncomment header below).
# ══════════════════════════════════════════════════════════════════

#$ -N xqtl_bam_qc
#$ -cwd
#$ -o qc_bam_qc.$JOB_ID.log
#$ -j y
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 8

set -euo pipefail

BAM_DIR="${1:-/u/project/kruglyak/bruni/projects/xQTLpaper/bam_files}"
OUT_DIR="${2:-/u/project/kruglyak/bruni/projects/xQTLpaper/qc/bam_stats}"
NCORES="${3:-8}"

mkdir -p "$OUT_DIR"

# Activate conda env if samtools isn't on PATH
if ! command -v samtools >/dev/null 2>&1; then
    source /u/local/Modules/default/init/bash
    module load mamba
    mamba activate genomicsV2
fi

echo "BAM_DIR : $BAM_DIR"
echo "OUT_DIR : $OUT_DIR"
echo "NCORES  : $NCORES"
echo "samtools: $(which samtools)  ($(samtools --version | head -1))"
echo

# Worker: runs flagstat + stats on one BAM, skips if already done
process_one() {
    local bam="$1"
    local out_dir="$2"
    local sid
    sid=$(basename "$bam" .bam)
    local fs="$out_dir/${sid}.flagstat"
    local st="$out_dir/${sid}.stats.txt"

    if [[ -s "$fs" && -s "$st" ]]; then
        echo "[skip] $sid"
        return 0
    fi

    samtools flagstat "$bam" > "${fs}.tmp" && mv "${fs}.tmp" "$fs"
    samtools stats     "$bam" > "${st}.tmp" && mv "${st}.tmp" "$st"
    echo "[done] $sid"
}
export -f process_one

# Parallelize with xargs
find "$BAM_DIR" -maxdepth 1 -name "Sample*.bam" -print0 \
    | xargs -0 -n1 -P "$NCORES" -I {} bash -c 'process_one "$@"' _ {} "$OUT_DIR"

echo
echo "==== summary ===="
n_flagstat=$(find "$OUT_DIR" -name '*.flagstat' | wc -l)
n_stats=$(find "$OUT_DIR" -name '*.stats.txt' | wc -l)
n_bams=$(find "$BAM_DIR" -maxdepth 1 -name 'Sample*.bam' | wc -l)
echo "BAMs found         : $n_bams"
echo "flagstat outputs   : $n_flagstat"
echo "stats.txt outputs  : $n_stats"
du -sh "$OUT_DIR"
