#!/usr/bin/env Rscript
## Run ON HOFFMAN2 (genomicsV2 active), from the directory holding acl_chrVIII_dead_sites.tsv:
##   qsub -V -cwd -o aclcheck.log -j y -l h_rt=1:00:00,h_data=16G -N aclcheck -b y \
##        "env OPENBLAS_NUM_THREADS=1 Rscript acl_chrVIII_a1_check.R"
##
## Question: at the 222 chrVIII sites where no ACL x BBQ pool ever shows the ACL
## allele (>= 20 reads, depth 0.90 of neighbours, so not read loss), what does
## OUR OWN sequencing of the ACL parent (A1 panel) carry?
##   - BBQ's allele  -> our ACL stock differs from its 1011 entry there; the
##                      genotype-dead flag is the right description.
##   - ACL's 1011 allele -> our stock does carry it and the pools lost it;
##                      the flag is mislabelled and something biological is going on.
## Reads only. Nothing is combined or patched.
suppressMessages(library(data.table))
P  <- "/u/project/kruglyak/bruni/projects/xQTLpaper"
A1 <- "/u/project/kruglyak/bruni/projects/A1ref"
sites <- fread("acl_chrVIII_dead_sites.tsv")
Zo <- qs2::qs_read(file.path(P, "files", "Z_final.qs")); Z <- Zo$Z; vid <- Zo$variant_ids
m  <- regmatches(vid, regexec("^([^_:]+)[_:]([0-9]+)", vid))
zs <- data.table(j = seq_along(vid), chrom = vapply(m, `[`, "", 2), pos = as.integer(vapply(m, `[`, "", 3)))
zs[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
zs <- zs[chrom == "chrVIII"]
vals <- sort(unique(as.numeric(Z["BBQ", zs$j]))); LO <- min(vals); HI <- max(vals)
zs[, `:=`(acl = as.numeric(Z["ACL", j]), bbq = as.numeric(Z["BBQ", j]))]
ref <- fread(file.path(A1, "20251119_refA1key.csv"))
lib <- ref[Standardized_name == "ACL" & STRAIN != "EMPTY", Sample_ID]
cat("A1 libraries for ACL:", paste(lib, collapse = ", "), "\n")
a1 <- rbindlist(lapply(lib, function(s) { f <- file.path(A1, "count_files", paste0(s, ".txt")); if (!file.exists(f)) return(NULL)
  x <- fread(f, select = c("contig", "position", "refCount", "altCount")); setnames(x, c("chrom", "pos", "ref", "alt"))
  x[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]; x[chrom == "chrVIII"][, lib := s] }))
a1 <- a1[, .(ref = sum(ref), alt = sum(alt)), by = .(chrom, pos)]
judge <- function(s) {
  x <- merge(merge(s, zs, by = c("chrom", "pos")), a1, by = c("chrom", "pos"), all.x = TRUE)
  x[, acl_allele := fifelse(acl == HI, "alt", fifelse(acl == LO, "ref", "het/NA"))]
  x[, a1_reads := fifelse(is.na(ref), 0L, as.integer(ref + alt))]
  x[, a1_call := fifelse(a1_reads < 3, "too few reads",
                 fifelse((acl_allele == "alt" & alt / a1_reads >= 0.9) | (acl_allele == "ref" & ref / a1_reads >= 0.9), "carries ACL-1011 allele",
                 fifelse((acl_allele == "alt" & alt / a1_reads <= 0.1) | (acl_allele == "ref" & ref / a1_reads <= 0.1), "carries BBQ allele", "mixed")))]
  x }
d <- judge(sites[, .(chrom, pos)])
cat("\n== the dead sites: what our ACL parent carries\n"); print(d[, .N, by = a1_call][order(-N)])
## control: live chrVIII sites where ACL and BBQ differ in 1011 and the pools DO see both alleles
live <- zs[acl %in% c(LO, HI) & bbq %in% c(LO, HI) & acl != bbq & !pos %in% sites$pos, .(chrom, pos)]
cat("\n== control, live differing chrVIII sites\n"); print(judge(live)[, .N, by = a1_call][order(-N)])
fwrite(d, "acl_chrVIII_a1_check.csv"); cat("\nwrote acl_chrVIII_a1_check.csv\n")
