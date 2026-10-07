#!/usr/bin/env Rscript
## Run ON HOFFMAN2 (qsub, genomicsV2 active):
##   qsub -V -cwd -o parentcheck.log -j y -l h_rt=4:00:00,h_data=16G -N parentcheck -b y \
##        "env OPENBLAS_NUM_THREADS=1 Rscript parent_genotype_check.R"
##
## Question: when a wild parent's 1011-genomes genotype (the Z row the driver
## phases against) disagrees with OUR OWN sequencing of that parent (the A1
## haploid panel), which one do the pools agree with?
##
## Measure: "dead" sites. At a site where the two parents differ, every pool
## of the cross must show the wild allele somewhere. A site with >= 20 reads
## summed over the cross's pools and ZERO wild-allele reads is a site where the
## assumed wild genotype is wrong (reads are present, at full depth, so this is
## not read loss). Fewer dead sites = better genotype source.
##
## Reads only. Writes ~/xqtl_parentcheck/ (small CSVs) and prints a summary.
suppressMessages(library(data.table))
P    <- "/u/project/kruglyak/bruni/projects/xQTLpaper"
A1   <- "/u/project/kruglyak/bruni/projects/A1ref"
OUT  <- path.expand("~/xqtl_parentcheck"); dir.create(OUT, showWarnings = FALSE)
POOL_COUNTS <- file.path(P, "v8", "count_files")   # both lanes
A1_COUNTS   <- file.path(A1, "count_files")
MIN_DEPTH_A1 <- 3      # reads needed to call a haploid A1 genotype
MIN_CROSS    <- 20     # reads summed over a cross's pools to judge a site
say <- function(...) cat(sprintf(...), "\n")

## ---- 1. Z: what it is and how it codes genotypes -----------------------------
Zo <- qs2::qs_read(file.path(P, "files", "Z_final.qs"))
Z  <- Zo$Z; vid <- Zo$variant_ids
say("Z class %s, %d strains x %d sites; variant_ids e.g. %s", class(Z)[1], nrow(Z), ncol(Z),
    paste(head(vid, 3), collapse = " | "))
smp <- as.vector(as.matrix(Z[, sample(ncol(Z), min(ncol(Z), 20000))]))
say("Z value table (20k-site sample):"); print(table(smp, useNA = "ifany"))
vals <- sort(unique(smp[!is.na(smp)]))
LO <- min(vals); HI <- max(vals)
say("coding assumed: %s = hom ref, %s = hom alt, anything between = het, NA = missing", LO, HI)
m <- regmatches(vid, regexec("^([^_:]+)[_:]([0-9]+)", vid))
zsite <- data.table(j = seq_along(vid), chrom = vapply(m, `[`, "", 2), pos = as.integer(vapply(m, `[`, "", 3)))
if (anyNA(zsite$pos)) stop("could not parse positions from variant_ids; see the examples above")
zsite[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
zcall <- function(code) { g <- as.numeric(Z[code, ]); fifelse(is.na(g), NA_real_, fifelse(g == LO, 0, fifelse(g == HI, 1, 0.5))) }

## ---- 2. keys ------------------------------------------------------------------
fk  <- fread(file.path(P, "files", "flatkey_A1.csv"))
ref <- fread(file.path(A1, "20251119_refA1key.csv"))[Project != "" & STRAIN != "EMPTY"]
setnames(ref, "Standardized_name", "code", skip_absent = TRUE)
ref <- ref[code != ""]
say("A1 panel: %d libraries, %d distinct 1011 codes", nrow(ref), uniqueN(ref$code))

readc <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- fread(f, select = c("contig", "position", "refCount", "altCount"), showProgress = FALSE)
  setnames(x, c("chrom", "pos", "ref", "alt")); x[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]; x
}

## ---- 3. A1 haploid genotypes, per code (replicate wells summed) ---------------
a1 <- rbindlist(lapply(split(ref, ref$code), function(r) {
  x <- rbindlist(lapply(file.path(A1_COUNTS, paste0(r$Sample_ID, ".txt")), readc))
  if (!nrow(x)) return(NULL)
  x <- x[, .(ref = sum(ref), alt = sum(alt)), by = .(chrom, pos)][ref + alt >= MIN_DEPTH_A1]
  x[, g := fifelse(alt / (ref + alt) >= 0.9, 1, fifelse(alt / (ref + alt) <= 0.1, 0, 0.5))]
  x[, code := r$code[1]]; x
}))
say("A1 genotype calls: %d site calls over %d codes", nrow(a1), uniqueN(a1$code))

## ---- 4. per parent: is the 1011 entry heterozygous where our haploid is not? --
codes <- intersect(unique(c(fk$parent_final, fk$parent_nominal)), intersect(rownames(Z), unique(a1$code)))
conc <- rbindlist(lapply(codes, function(cd) {
  z <- zsite[, .(chrom, pos, gz = zcall(cd))][!is.na(gz)]
  x <- merge(a1[code == cd & g != 0.5, .(chrom, pos, ga = g)], z, by = c("chrom", "pos"))
  data.table(code = cd, sites_compared = nrow(x),
             z_het_pct = round(100 * mean(z$gz == 0.5), 2),
             disagree_at_z_hom_pct = round(100 * mean(x[gz != 0.5, ga != gz]), 2),
             a1_hom_at_z_het = x[gz == 0.5, .N])
}))
setorder(conc, -disagree_at_z_hom_pct)
fwrite(conc, file.path(OUT, "parent_concordance_A1_vs_1011.csv"))
say("\n== our haploid parent (A1) vs its 1011 entry, per code"); print(conc)

## ---- 5. per cross: dead-site share phasing against 1011 vs against A1 ---------
## BBQ from Z (it is the common parent everywhere). Wild genotype from each source.
bbq <- zsite[, .(chrom, pos, gb = zcall("BBQ"))][gb %in% c(0, 1)]
well_counts <- function(well) {
  pools <- fk[Wells == well, Sample_ID]
  x <- rbindlist(lapply(file.path(POOL_COUNTS, paste0(pools, ".txt")), readc))
  if (!nrow(x)) return(NULL)
  x[, .(ref = sum(ref), alt = sum(alt)), by = .(chrom, pos)]
}
dead_share <- function(x, wild_g) {
  ## x: summed pool counts for one plate well (all 4 conditions); wild_g: (chrom,pos,gw), gw in {0,1}
  if (is.null(x)) return(list(n = 0L, dead = NA_real_))
  s <- merge(wild_g, bbq, by = c("chrom", "pos"))[gw != gb]
  s <- merge(s, x, by = c("chrom", "pos"))[ref + alt >= MIN_CROSS]
  s[, wild_reads := fifelse(gw == 1, alt, ref)]
  list(n = nrow(s), dead = round(100 * mean(s$wild_reads == 0), 2))
}
src_z  <- function(cd) if (cd %in% rownames(Z)) zsite[, .(chrom, pos, gw = zcall(cd))][gw %in% c(0, 1)] else NULL
src_a1 <- function(cd) if (cd %in% a1$code) a1[code == cd & g != 0.5, .(chrom, pos, gw = g)] else NULL
wells <- unique(fk[, .(Wells, parent_final, parent_nominal)])
if (anyDuplicated(wells$Wells)) say('!! some wells carry more than one parent call across plates; each row is judged on all 4 pools of the well')
res <- rbindlist(lapply(seq_len(nrow(wells)), function(i) {
  w <- wells[i]; out <- list(well = w$Wells, final = w$parent_final, nominal = w$parent_nominal)
  x <- well_counts(w$Wells)
  for (cand in unique(c(w$parent_final, w$parent_nominal))) for (src in c("1011", "A1")) {
    g <- if (src == "1011") src_z(cand) else src_a1(cand)
    r <- if (is.null(g)) list(n = NA, dead = NA) else dead_share(x, g)
    out[[paste0(if (cand == w$parent_final) "final" else "nominal", "_", src, "_dead_pct")]] <- r$dead
    out[[paste0(if (cand == w$parent_final) "final" else "nominal", "_", src, "_sites")]] <- r$n
  }
  as.data.table(out)
}), fill = TRUE)
fwrite(res, file.path(OUT, "dead_share_by_well_and_source.csv"))
say("\n== dead-site share (%%) per plate well, by which parent and which genotype source")
print(res)
say("\nA10 and A11 (the two AAL wells):"); print(res[well %in% c("A10", "A11")])
say("\nwrote %s", OUT)
