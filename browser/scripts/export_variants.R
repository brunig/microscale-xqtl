suppressMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
RUN  <- if (length(args) >= 1) args[1] else "hoffmanDrop/notAF/v7"   # a run umbrella, e.g. runs/v7
SITE <- normalizePath(if (length(args) >= 2) args[2] else "browser/docs", mustWork = FALSE)
GFF  <- if (length(args) >= 3) args[3] else Sys.getenv("XQTL_GFF", "")
setwd(RUN)
OUT <- file.path(SITE, "v")
dir.create(OUT, showWarnings=FALSE)
CONS <- c(synonymous=1L, nonsynonymous=2L, frameshift=3L, stop_gained=4L)
gtof <- function(g, who) { m <- regmatches(g, regexpr(paste0("(^|;)", who, "= ?[0-9]"), g)); x <- rep(NA_integer_, length(g)); ok <- grepl(paste0("(^|;)", who, "= ?[0-9]"), g); x[ok] <- as.integer(sub(".*= ?", "", m)); x }
tot <- 0
for (set in c("fitness","stress")) {
  v <- fread(paste0("contrasts/qtl_variants_", set, ".csv"))
  v[, wild := sub("xBBQ$", "", comparison)]
  v[, gB := gtof(genotypes, "BBQ")]
  v[, gW := mapply(function(g, w) { m <- regmatches(g, regexpr(paste0("(^|;)", w, "= ?[0-9]"), g)); if (length(m)) as.integer(sub(".*= ?", "", m)) else NA_integer_ }, genotypes, wild)]
  v[, ra := sub("^[^:]+:", "", idx)]
  v[, aa := ifelse(grepl("character\\(0\\)|numeric\\(0\\)", AAeffect), "", sub("^[^:]+:([^:]*):([^:]*):(.*)$", "\\1\\2\\3", AAeffect))]
  v[, cons := fifelse(is.na(CONSEQUENCE), 0L, CONS[CONSEQUENCE])]
  v[, hiW := as.integer(high.parent != "BBQ")]
  cn <- fread(paste0("contrasts/qtl_cnvs_", set, ".csv"))
  for (ch in unique(v$chr)) {
    x <- v[chr == ch][order(pos)]
    keys <- unique(x[, .(comparison, peak, hiW, peak.pos)]); keys[, k := .I - 1L]
    x <- merge(x, keys[, .(comparison, peak, k)], by = c("comparison", "peak"), sort = FALSE)[order(pos)]
    y <- cn[chr == ch]
    out <- list(
      peaks = list(cross = keys$comparison, key = keys$peak, hiW = keys$hiW, summit = keys$peak.pos),
      v = list(k = x$k, pos = x$pos, ra = x$ra, g = fifelse(is.na(x$GENEID), "", x$GENEID), c = x$cons, aa = x$aa,
               pv = round(x$provean, 2), ph = round(x$phyloP, 2), gB = x$gB, gW = x$gW, ind = as.integer(x$type == "indel"),
               dub = as.integer(x$dubious %in% TRUE), ess = as.integer(x$essential %in% TRUE)),
      cnv = list(cross = y$comparison, key = y$peak, a = y$start, e = y$end, strain = y$Strain, type = y$CNV_type, rd = round(y$normalized_RD, 2), q = signif(y$q0, 2)))
    f <- file.path(OUT, paste0(ifelse(set == "fitness", "fit", "stress"), "_", ch, ".json"))
    writeLines(toJSON(out, auto_unbox = FALSE, digits = NA, na = "null"), f); tot <- tot + file.size(f)
  }
  cat(set, "rows", nrow(v), " het-wild rows", sum(v$gW == 1, na.rm=T), " NA gW", sum(is.na(v$gW)), "\n")
}
cand <- rbind(fread("candidates/fitness/qtl_candidates_fitness.csv")[, set := "fit"], fread("candidates/stress/qtl_candidates_stress.csv")[, set := "stress"], fill = TRUE)
keep <- c("set","locus","condition","chrom","summit","cons_kb","locus_crosses","locus_pools","gene","name","orf","dist_kb","crosses_cov","n_seg","alleles_wild","crosses_wild","dam_frac","q_allelic","n_lof","n_reg","min_provean","max_phylop","score","rank_in_locus","published","wild_strains")
cand <- cand[, intersect(keep, names(cand)), with = FALSE]
f <- file.path(OUT, "candidates.json"); writeLines(toJSON(as.list(cand), digits = NA, na = "null"), f); tot <- tot + file.size(f)
cat("candidate rows", nrow(cand), "by set", table(cand$set), " conditions", unique(cand$condition), "\n")
cat("total MB", round(tot/1e6,1), " files", length(list.files(OUT)), "\n")
