## Effective segregant number per pool, from the smoothed allele-frequency trace.
##
## LOGIC. In an unselected pool of N independent haploid segregants, the allele
## frequency at a locus is a mean of N Bernoulli(0.5) draws, so it deviates from
## 0.5 with SD = 0.5/sqrt(N). Reading that pool at depth d adds a further
## sqrt(0.25/d) of sampling noise, but the loess fit averages over many sites and
## removes most of it, so the SMOOTHED trace at a locus under no selection is
## close to the pool's true allele frequency. Take its SD across UNLINKED,
## PEAK-FREE loci and invert:  N_eff = 0.25 / sigma^2.
##
## The diagnostic that decides whether this is measuring segregants at all is
## whether sigma falls with sequencing depth. If it does, we are still reading
## noise. If it is flat, we have hit the floor set by the pool itself.
suppressPackageStartupMessages({library(data.table); library(qs2)})
RUN <- "hoffmanDrop/notAF/v7"
key <- fread(file.path(RUN,"qc/results_key.csv"))[databin=="main" & Drug3Let=="Glu"]
pk  <- fread(file.path(RUN,"contrasts/fit_peaks.csv"))[databin=="main"]
cf  <- Sys.glob(file.path(RUN,"cache","*_fit_results.qs"))
cf  <- cf[order(file.info(cf)$size, decreasing=TRUE)][1]
cat("cache:", basename(cf), "\n"); fitr <- qs_read(cf)

EXCL <- data.table(chrom=c("chrIII","chrIV","chrVII"),
                   lo=c(125000,15000,440000), hi=c(215000,60000,600000))
CHRLEN <- c(chrI=230218,chrII=813184,chrIII=316620,chrIV=1531933,chrV=576874,
  chrVI=270161,chrVII=1090940,chrVIII=562643,chrIX=439888,chrX=745751,chrXI=666816,
  chrXII=1078177,chrXIII=924431,chrXIV=784333,chrXV=1091291,chrXVI=948066)

est <- rbindlist(lapply(key$Sample_ID, function(sid){
  r <- fitr[[sid]]; if (is.null(r)) return(NULL)
  af <- r[[grep("^afd_", names(r), value=TRUE)[1]]]
  d  <- data.table(chrom=as.character(r$chrom), pos=r$physical.position, af=af)
  d  <- d[!is.na(af) & af>0 & af<1]
  ## drop the design-selected regions and every called peak interval of THIS pool
  for(i in 1:nrow(EXCL)) d <- d[!(chrom==EXCL$chrom[i] & pos>=EXCL$lo[i] & pos<=EXCL$hi[i])]
  pp <- pk[sample==sid]
  if (nrow(pp)) for(i in 1:nrow(pp))
    d <- d[!(chrom==as.character(pp$seqnames[i]) & pos>=pp$start[i] & pos<=pp$end[i])]
  if (nrow(d) < 2000) return(NULL)
  ## one value per 200 kb block: far enough apart to be near-unlinked in a
  ## segregant pool (yeast is ~2-4 cM/kb-scale; 200 kb is many crossovers)
  d[, blk := paste(chrom, pos %/% 200000)]
  v <- d[, .(af=median(af)), by=blk]$af
  if (length(v) < 12) return(NULL)
  ## centre on the pool OWN median, not 0.5: a pool whose whole trace sits off
  ## 0.5 has a skew problem, which is a different artefact and is flagged
  ## separately. Centring on 0.5 would score that offset as segregant spread.
  sigma <- mad(v)
  sigma0 <- mad(v, center=0.5)
  data.table(sample=sid, cross=key$cross[match(sid,key$Sample_ID)],
             depth=key$med_depth[match(sid,key$Sample_ID)],
             blocks=length(v), sigma=sigma, sigma0=sigma0,
             N_eff=0.25/sigma^2, N_eff0=0.25/sigma0^2)
}))
setorder(est, -depth)
fwrite(est, "draft/derived/pool_size_estimates.csv")
cat("wrote draft/derived/pool_size_estimates.csv\n")
cat("\n== per-pool effective segregant number, glucose pools ==\n")
print(est[, .(sample, cross, depth, blocks, sigma=round(sigma,4), N_eff=round(N_eff))], nrows=25)
cat(sprintf("\nmedian N_eff = %.0f   (IQR %.0f - %.0f, n = %d pools)\n",
    median(est$N_eff), quantile(est$N_eff,.25), quantile(est$N_eff,.75), nrow(est)))
cat("\n== does sigma fall with depth?  (if yes we are still reading noise) ==\n")
print(est[, .(depth_band=cut(depth, c(0,5,10,15,25,Inf)),
              sigma, N_eff)][, .(pools=.N, median_sigma=round(median(sigma),4),
                                 median_N=round(median(N_eff))), by=depth_band][order(depth_band)])
ct <- cor.test(est$depth, est$sigma, method="spearman")
cat(sprintf("\nSpearman(depth, sigma) = %.3f, P = %.3g\n", ct$estimate, ct$p.value))

cat("\n== headline, pools at 10x or better (where sigma has plateaued) ==\n")
hi <- est[depth >= 10]
cat(sprintf("  n = %d pools\n  sigma  median %.4f  ->  N_eff median %.0f  (IQR %.0f - %.0f)\n",
    nrow(hi), median(hi$sigma), median(hi$N_eff),
    quantile(hi$N_eff,.25), quantile(hi$N_eff,.75)))
cat(sprintf("  centred on 0.5 instead: N_eff median %.0f  (the gap is per-pool skew, not segregants)\n",
    median(hi$N_eff0)))
cat("\n== what the pipeline ASSUMES for the standard error ==\n")
cat("  afd_sample_size 1e4 x sel_strength 0.9 = 9000, capped by total depth / 600.\n")
cat(sprintf("  At a typical 20x over ~24,400 informative sites that cap is ~%.0f.\n", 20*24434/600))
cat(sprintf("  Measured N_eff is ~%.0f, so the SE is too small by about %.1f-fold,\n",
    median(hi$N_eff), sqrt((20*24434/600)/median(hi$N_eff))))
cat("  which is a quantified form of the missing overdispersion term.\n")
