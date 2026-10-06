#!/bin/bash
# Run ON HOFFMAN2 (login or dtn node is fine; nothing heavy):  bash hoffman_versions.sh
# Records the software versions of the conda envs the pipeline ran in.
# Writes ~/xqtl_versions.txt.
P=/u/project/kruglyak/bruni/projects/xQTLpaper
OUT=~/xqtl_versions.txt
source /u/local/Modules/default/init/bash
module load mamba
{
for ENV in genomicsV2 genomics; do
  echo "######## conda env: $ENV"
  mamba activate "$ENV" 2>/dev/null || { echo "(env not found)"; continue; }
  echo "== tools"
  for t in nextflow bwa samtools gatk bcftools sambamba java python R; do
    printf '%-10s %s\n' "$t" "$(command -v $t || echo -)"
  done
  nextflow -version 2>&1 | grep -i version
  bwa 2>&1 | grep -i '^Version'
  samtools --version 2>/dev/null | head -1
  gatk --version 2>/dev/null | grep -iE 'gatk|htsjdk'
  bcftools --version 2>/dev/null | head -1
  java -version 2>&1 | head -1
  python --version 2>&1
  echo "== conda list (key packages)"
  mamba list 2>/dev/null | grep -E '^(nextflow|bwa|samtools|gatk4|bcftools|sambamba|htslib|openjdk|python|r-base) '
  echo "== R packages"
  Rscript -e '
    cat(R.version.string, "\n")
    pk <- c("data.table","dplyr","tidyr","readr","purrr","ggplot2","ggrepel","patchwork","scales",
            "qs2","pbmcapply","IRanges","GenomicRanges","BiocGenerics","ComplexHeatmap","circlize",
            "xQTLStats","jsonlite","ragg")
    for (p in pk) {
      v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
      d <- tryCatch(packageDescription(p), error = function(e) NULL)
      src <- if (!is.null(d$RemoteSha)) paste0(" (", d$RemoteUsername, "/", d$RemoteRepo, "@", substr(d$RemoteSha, 1, 10), ")") else ""
      cat(sprintf("%-15s %s%s\n", p, v, src))
    }' 2>&1
  mamba deactivate
done
echo "######## Nextflow run reports in pipeline_info/"
ls -la "$P/pipeline_info" | head -20
grep -hoE 'nextflow[^<"]{0,40}version[^<"]{0,30}|Nextflow version[^<]{0,30}' "$P"/pipeline_info/*.html 2>/dev/null | sort -u | head
echo "######## Driver frozen in runs/v7: head of manifest-relevant lines"
grep -nE '^\s*library\(|xQTLStats' "$P/runs/v7/scripts/xQTL_analysis7.R" | head -20
} > "$OUT" 2>&1
echo "wrote $OUT"
