suppressMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
RUN  <- if (length(args) >= 1) args[1] else "hoffmanDrop/notAF/v7"   # a run umbrella, e.g. runs/v7
SITE <- normalizePath(if (length(args) >= 2) args[2] else "browser/docs", mustWork = FALSE)
GFF  <- if (length(args) >= 3) args[3] else Sys.getenv("XQTL_GFF", "")
TAG  <- if (length(args) >= 4) args[4] else "v7"   # run tag: cache file prefix and the build label
setwd(RUN)
OUT <- file.path(SITE, "data.json")
BW <- 2000
CHRLEN <- c(chrI=230218, chrII=813184, chrIII=316620, chrIV=1531933, chrV=576874, chrVI=270161, chrVII=1090940,
  chrVIII=562643, chrIX=439888, chrX=745751, chrXI=666816, chrXII=1078177, chrXIII=924431, chrXIV=784333, chrXV=1091291, chrXVI=948066)
nb <- ceiling(CHRLEN/BW); off <- c(0, cumsum(nb))[1:16]; names(off) <- names(CHRLEN); NB <- sum(nb)
binidx <- function(ch, pos) off[as.character(ch)] + floor(pos/BW) + 1
q01 <- function(v) { r <- round(pmin(pmax(v,0),1)*250); r[is.na(r)] <- 255L; as.integer(r) }
qpm <- function(v) { r <- round((pmin(pmax(v,-1),1)+1)*125); r[is.na(r)] <- 255L; as.integer(r) }
b64 <- function(i) base64_enc(as.raw(i))
binmean <- function(idx, v) { s <- rep(NA_real_, NB); ok <- !is.na(v); a <- tapply(v[ok], idx[ok], mean); s[as.integer(names(a))] <- a; s }
binsum <- function(idx, v) { s <- numeric(NB); a <- tapply(v, idx, sum); s[as.integer(names(a))] <- a; s }

k <- read.csv("qc/results_key.csv", stringsAsFactors=FALSE)
fit <- qs2::qs_read(paste0("cache/", TAG, "_fit_results.qs")); drg <- qs2::qs_read(paste0("cache/", TAG, "_drug_results.qs"))
pools <- list(); tr <- list(); raw <- list()
for (s in names(fit)) { d <- fit[[s]]; idx <- binidx(d$chrom, d$physical.position)
  af <- d[[paste0("afd_",s)]]; p1 <- d[[paste0("p1_",s)]]; p2 <- d[[paste0("p2_",s)]]
  tr[[s]] <- b64(q01(binmean(idx, af)))
  n1 <- binsum(idx,p1); nt <- binsum(idx,p1+p2); rf <- ifelse(nt>=4, n1/nt, NA); raw[[s]] <- b64(q01(rf))
  h <- !is.na(d$LOD) & d$LOD>=5 & !is.na(d$contrast.beta); pools[[s]] <- list(se=signif(median(d[[paste0("afd.se_",s)]],na.rm=T),4), th=if(any(h)) signif(min(abs(d$contrast.beta[h])),4) else NULL) }
con <- list()
for (nm in names(drg)) { e <- drg[[nm]]; d <- e$result; idx <- binidx(d$chrom, d$physical.position)
  h <- !is.na(d$LOD) & d$LOD>=5 & !is.na(d$contrast.beta); con[[e$meta$L]] <- list(ctrl=e$meta$R, se=signif(median(d$contrast.beta.se,na.rm=T),4), th=if(any(h)) signif(min(abs(d$contrast.beta[h])),4) else NULL, b=b64(qpm(binmean(idx, d$contrast.beta)))) }
cn <- c(Glu="Glu",Ace="Ace",Caf="Caf",NaCl="NaCl")
meta <- lapply(seq_len(nrow(k)), function(i) { r <- k[i,]; s <- r$Sample_ID
  list(id=s, plate=r$Plate, well=r$Wells, cross=r$cross, alt=r$parent_a, nominal=r$oldLabel_parent_a, cond=r$condition, bin=r$databin,
       flags=ifelse(is.na(r$flag_reasons),"",r$flag_reasons), med=r$med_depth, mean=round(r$mean_depth_seg,1),
       reassign=ifelse(is.na(r$reassign_class),"",r$reassign_class), decision=ifelse(is.na(r$decision),"",r$decision),
       nfit=r$n_peaks_fit, nstress=r$n_peaks_stress,
       se=if (!is.null(pools[[s]])) pools[[s]]$se else NULL, th=if (!is.null(pools[[s]])) pools[[s]]$th else NULL,
       tr=tr[[s]], raw=raw[[s]], con=con[[s]]) })
MR <- read.csv("landscape/masked_regions.csv"); MR <- MR[MR$what %in% c("MAT","HO"),]
maskit <- function(ch, st, en, pos) { m <- rep(FALSE, length(ch)); for (j in seq_len(nrow(MR))) { lo <- MR$lo[j]; hi <- MR$hi[j]
  ov <- pmax(0, pmin(en,hi) - pmax(st,lo)); m <- m | (ch==MR$chrom[j] & ((pos>=lo & pos<=hi) | ov > 0.5*(hi-lo))) }; m }
fp <- read.csv("contrasts/fit_peaks.csv"); sp <- read.csv("contrasts/stress_peaks.csv")
fpo <- data.frame(s=fp$sample, c=fp$seqnames, a=fp$start, e=fp$end, p=fp$peak.pos, b=round(fp$peak.beta,4), m=as.integer(maskit(fp$seqnames,fp$start,fp$end,fp$peak.pos)))
spo <- data.frame(s=sp$L, c=sp$seqnames, a=sp$start, e=sp$end, p=sp$peak.pos, b=round(sp$peak.beta,4), m=as.integer(maskit(sp$seqnames,sp$start,sp$end,sp$peak.pos)))
cat("fit main peaks after mask:", sum(fp$databin=="main" & fpo$m==0), "\n")
pub <- read.csv("landscape/published_loci_used.csv")
pubo <- data.frame(q=pub$qtl, cond=pub$cond, c=pub$chrom, a=pub$lo, e=pub$hi, p=pub$maxpos, n=pub$ncross, lod=pub$lod, st=pub$study)
g <- read.delim(gzfile(GFF), header=FALSE, comment.char="#", quote="", stringsAsFactors=FALSE, fill=TRUE)
g <- g[g$V3=="gene" & g$V1 %in% names(CHRLEN),]
att <- function(a, key) { m <- regmatches(a, regexpr(paste0(key,"=[^;]+"), a)); out <- rep("", length(a)); out[grepl(paste0(key,"="),a)] <- sub(".*=","",m); utils::URLdecode(out) }
gid <- att(g$V9,"ID"); gnm <- att(g$V9,"gene"); cls <- att(g$V9,"orf_classification")
cls <- vapply(cls, utils::URLdecode, "")
genes <- data.frame(c=g$V1, a=g$V4, e=g$V5, s=g$V7, id=gid, n=gnm, k=substr(cls,1,1))
cat("genes:", nrow(genes), table(genes$k), "\n")
out <- list(build=TAG, bw=BW, chroms=as.list(CHRLEN), off=as.list(off), nbins=NB,
  masks=read.csv("landscape/masked_regions.csv"), pools=meta, fitpeaks=fpo, stresspeaks=spo, published=pubo, genes=genes)
writeLines(toJSON(out, auto_unbox=TRUE, dataframe="columns", digits=NA, null="null", na="null"), OUT)
cat("bytes:", file.size(OUT), "\n")
