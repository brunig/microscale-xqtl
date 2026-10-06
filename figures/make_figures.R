#!/usr/bin/env Rscript
## make_figures.R — every manuscript figure from run v6b, in one go.
##
##   Rscript draft/figures/make_figures.R
##
## Figure 1 is hand-drawn and is not built here.  Figures 4D and 5E are GO-term
## panels that cannot be built: nothing in v6b emits GO terms.  They are drawn as
## explicitly empty placeholders rather than filled with a substitute result.
##
## Conventions enforced here, from CLAUDE.md:
##   - condition order is glucose, acetate, caffeine, sodium chloride
##   - `main` pools only; MASK_REPORT (MAT + HO, shoulder 0.5) for every figure
##   - percentages quoted against the RETAINED set, and labelled as such
##   - effect size (beta) is the plotted quantity; LOD is only ever a cut

suppressPackageStartupMessages({
  library(IRanges); library(GenomicRanges)      # loaded first: dplyr must mask
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(readr)
})

## ----------------------------------------------------------------- config --
ROOT <- "/Users/kruglyakbws/Desktop/AFxQTL"
RUN  <- file.path(ROOT, "hoffmanDrop/notAF/v7")   # v6b + the V7 additions
OUT  <- file.path(ROOT, "draft/figures/generated")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LOD_CUT <- 5                      # the peak caller's cut; thresholds derive from it
CONDS   <- c("Glu", "Ace", "Caf", "NaCl")
LONG    <- c(Glu = "glucose", Ace = "acetate",
             Caf = "caffeine", NaCl = "NaCl")
CCOL    <- c("glucose" = "#666666", "acetate" = "#E69F00",
             "caffeine" = "#009E73", "NaCl" = "#0072B2")

CHRLEN <- c(chrI=230218, chrII=813184, chrIII=316620, chrIV=1531933, chrV=576874,
            chrVI=270161, chrVII=1090940, chrVIII=562643, chrIX=439888, chrX=745751,
            chrXI=666816, chrXII=1078177, chrXIII=924431, chrXIV=784333,
            chrXV=1091291, chrXVI=948066)
CHROMS <- names(CHRLEN); GAP <- 120000
OFFSET <- cumsum(c(0, head(CHRLEN + GAP, -1))); names(OFFSET) <- CHROMS
GEND   <- sum(CHRLEN + GAP) - GAP

## Regions the experiment selects on.  chrVII is the shared-parent sweep: it is a
## finding, not an artefact, so it is censored only where asked.
MASK <- tibble::tribble(~name,  ~chrom,   ~lo,     ~hi,
                        "MAT",  "chrIII", 125000,  215000,
                        "HO",   "chrIV",   15000,   60000)
CHR7 <- tibble::tibble(name = "chrVII sweep", chrom = "chrVII", lo = 440000, hi = 600000)
SHOULDER_FRAC <- 0.5

theme_pub <- function(base = 9)
  theme_classic(base_size = base) +
  theme(plot.tag = element_text(face = "bold", size = base + 3),
        plot.title = element_text(size = base, face = "bold"),
        plot.subtitle = element_text(size = base - 1, colour = "grey30"),
        axis.text = element_text(colour = "black"),
        ## the two-line contrast title was sitting a long way off its ticks
        axis.title.y = element_text(margin = margin(r = 1)),
        axis.text.x = element_text(size = base - 3),
        axis.title.x = element_text(margin = margin(t = 2)),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key.size = unit(3.2, "mm"))

## Crosses are written BBQ x WILD with a true multiplication sign, never the
## internal pool id: "Sample161" tells a reader nothing and cannot be looked up
## outside this project.
cross_label <- function(cross) sub("^(.*)xBBQ$", "BBQ \u00d7 \\1", cross)

say <- function(...) cat(sprintf(...), "\n", sep = "")
## Each figure ships its caption beside it, generated from the same values that
## made the figure, so the two cannot drift.
write_legend <- function(file, ...) {
  writeLines(c(...), file.path(OUT, file))
  say("  wrote %-42s (legend)", file)
}
wrap_para <- function(x, width = 96) paste(strwrap(x, width), collapse = "\n")
save_fig <- function(p, file, w, h, dpi = 300) {
  ggsave(file.path(OUT, file), p, width = w, height = h, units = "in",
         dpi = dpi, device = ragg::agg_png, bg = "white", limitsize = FALSE)
  say("  wrote %-42s %.1f x %.1f in", file, w, h)
}

gx <- function(chrom, pos) unname(OFFSET[as.character(chrom)]) + pos
chrom_axis <- function() {
  mids <- OFFSET + CHRLEN / 2
  list(scale_x_continuous(breaks = mids, labels = sub("^chr", "", CHROMS),
                          limits = c(0, GEND), expand = expansion(mult = 0.004)),
       theme(panel.grid = element_blank()))
}
chrom_bands <- function(ymin, ymax) {
  b <- tibble::tibble(chrom = CHROMS, xmin = OFFSET, xmax = OFFSET + CHRLEN,
                      i = seq_along(CHROMS)) |> filter(i %% 2 == 0)
  geom_rect(data = b, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey96", colour = NA, inherit.aes = FALSE)
}

## ------------------------------------------------------------------- load --
say("reading run v7")
key <- read_csv(file.path(RUN, "qc/results_key.csv"), show_col_types = FALSE)
main <- key |> filter(databin == "main")
fitp <- read_csv(file.path(RUN, "contrasts/fit_peaks.csv"), show_col_types = FALSE)
pub  <- read_csv(file.path(RUN, "landscape/published_loci_used.csv"), show_col_types = FALSE)
enr  <- read_csv(file.path(RUN, "landscape/published_overlap_enrichment.csv"), show_col_types = FALSE)
rec  <- read_csv(file.path(RUN, "landscape/published_recovery.csv"), show_col_types = FALSE)


## MASK_REPORT: summit inside a masked region, OR interval covering more than
## SHOULDER_FRAC of it.  Same rule as qtl_landscape.R.
apply_mask <- function(pk, regions) {
  drop <- rep(FALSE, nrow(pk))
  for (i in seq_len(nrow(regions))) {
    r <- regions[i, ]
    same <- pk$seqnames == r$chrom
    inside <- same & pk$peak.pos >= r$lo & pk$peak.pos <= r$hi
    ov <- pmax(0, pmin(pk$end, r$hi) - pmax(pk$start, r$lo))
    drop <- drop | inside | (same & ov > SHOULDER_FRAC * (r$hi - r$lo))
  }
  pk[!drop, ]
}
fit_main <- fitp |> filter(databin == "main")
fit_keep <- apply_mask(fit_main, MASK)
say("  %d fitness peaks in main, %d retained after MASK_REPORT, in %d pools / %d crosses",
    nrow(fit_main), nrow(fit_keep), n_distinct(fit_keep$sample), n_distinct(fit_keep$cross))
## V7. The published comparisons now exist on both peak sets and both units of
## detection. Treusch et al. mapped with MULTIPOOL in contrast mode, so the
## comparable set is the paired stress contrast, and the unit is the CROSS --
## replicate wells of one cross are not independent evidence. Figure 4 uses
## peak_set == "contrast" and unit == "cross" throughout; the by-pool numbers
## are carried in the legend only, because the gap between them is the argument
## for the cross being the unit.
enr_s <- read_csv(file.path(RUN, "landscape/published_overlap_enrichment_bysets.csv"),
                  show_col_types = FALSE)
rec_s <- read_csv(file.path(RUN, "landscape/published_recovery_bysets.csv"),
                  show_col_types = FALSE)
strp  <- read_csv(file.path(RUN, "contrasts/stress_peaks.csv"), show_col_types = FALSE)
str_main <- strp |> filter(databin == "main")
str_keep <- apply_mask(str_main, MASK)
## The cache filename prefix is NOT the run tag -- it is a filename convention
## inherited from V4/V5, and it travels with a hardlinked cache. runs/v7/cache
## ended up holding BOTH v6_* (hardlinked from v6b) and v7_* (written by the
## from-scratch run), and only the v7_ set is complete. Resolve it rather than
## hardcoding either, or the script breaks depending on which run last wrote.
cache_file <- function(stem) {
  hits <- Sys.glob(file.path(RUN, "cache", paste0("*_", stem, ".qs")))
  if (!length(hits)) stop("no cache file for '", stem, "' under ",
                          file.path(RUN, "cache"), call. = FALSE)
  if (length(hits) > 1) {
    hits <- hits[order(file.info(hits)$size, decreasing = TRUE)]
    say("  note: %d cache files match '%s'; using %s",
        length(hits), stem, basename(hits[1]))
  }
  hits[1]
}

CELL <- function(d, ps = "contrast", un = "cross")
  d |> filter(peak_set == ps, unit == un)
say("  %d contrast peaks in main, %d retained, in %d pools / %d crosses",
    nrow(str_main), nrow(str_keep), n_distinct(paste(str_keep$cross, str_keep$well, str_keep$condition)),
    n_distinct(str_keep$cross))

## ------------------------------------------------- per-site data on demand --
## the fit_results cache is 437 MB and 1.7 GB in memory, so it is read once and
## only the pools the figures need are kept.
FIT_CACHE <- new.env(parent = emptyenv())
load_sites <- function(samples) {
  need <- setdiff(samples, ls(FIT_CACHE))
  if (!length(need)) return(invisible(NULL))
  say("  loading per-site data for %s", paste(need, collapse = ", "))
  x <- qs2::qs_read(cache_file("fit_results"))
  for (s in need) {
    if (!s %in% names(x)) { warning("no per-site data for ", s); next }
    d <- as.data.frame(x[[s]])
    p1 <- d[[paste0("p1_", s)]]; p2 <- d[[paste0("p2_", s)]]
    afd <- d[[paste0("afd_", s)]]
    tot <- p1 + p2
    a <- p1 / tot; b <- p2 / tot            # phase orientation is not documented:
    ca <- suppressWarnings(cor(a, afd, use = "complete.obs"))
    cb <- suppressWarnings(cor(b, afd, use = "complete.obs"))
    raw <- if (is.na(ca) || (!is.na(cb) && cb > ca)) b else a
    assign(s, tibble::tibble(chrom = d$chrom, pos = d$physical.position,
                             raw = raw, afd = afd,
                             afd.se = d[[paste0("afd.se_", s)]],
                             beta = d$contrast.beta, beta.se = d$contrast.beta.se,
                             LOD = d$LOD) |> filter(!is.na(afd)), envir = FIT_CACHE)
  }
  rm(x); gc(verbose = FALSE); invisible(NULL)
}

## The significance threshold is derived from the peak caller itself: the smallest
## effect size that reaches the LOD cut in that pool.  Expressed in allele-frequency
## units so it can be drawn on both panels.
beta_threshold <- function(d) {
  hit <- d |> filter(!is.na(LOD), LOD >= LOD_CUT, !is.na(beta))
  if (nrow(hit)) min(abs(hit$beta)) else 4.3 * median(d$beta.se, na.rm = TRUE)
}

## Panel builders.  The AF panel is the raw per-site frequency with the smoothed
## fit over it; the delta panel is a contrast trace with no points.  Both carry the
## same threshold and the same shading of the regions the design selects on.
sel_bands <- function(regions) {
  if (is.null(regions) || !nrow(regions)) return(NULL)
  r <- regions |> mutate(g0 = gx(chrom, lo), g1 = gx(chrom, hi))
  ## labelled, not just shaded: MAT sits left of its band and HO right of its
  ## own, because the two bands are only ~250 kb apart on the genome axis
  r <- r |> mutate(hj = ifelse(name == "MAT", 1.1, -0.1),
                   xl = ifelse(name == "MAT", g0, g1))
  list(geom_rect(data = r, aes(xmin = g0, xmax = g1, ymin = -Inf, ymax = Inf),
                 fill = "#d9c9a3", alpha = .55, colour = NA, inherit.aes = FALSE),
       geom_text(data = r, aes(x = xl, y = Inf, label = name, hjust = hj),
                 vjust = 1.3, size = 2.5, fontface = "italic", colour = "#8a7442",
                 inherit.aes = FALSE))
}
mark_layers <- function(mark) {
  if (is.null(mark)) return(list(NULL, NULL))
  list(geom_vline(data = mark, aes(xintercept = gx(chrom, pos)), colour = "#B2182B",
                  linewidth = .4, linetype = "22", inherit.aes = FALSE),
       geom_text(data = mark, aes(x = gx(chrom, pos), y = Inf, label = label),
                 vjust = 1.4, hjust = -0.06, size = 2.6, colour = "#B2182B",
                 fontface = "italic", inherit.aes = FALSE))
}
## Two-line contrast axis title: "D allele frequency" over "<arm> - <arm>".
delta_lab <- function(minuend, subtrahend = "glucose")
  bquote(atop(Delta*" allele frequency", .(sprintf("%s - %s", minuend, subtrahend))))

COL_FIT <- "#2a78d6"     # allele frequency
COL_DEL <- "#D55E00"     # a contrast is a different quantity, so a different hue

panel_af <- function(d, th, mark = NULL, regions = MASK, parents = NULL) {
  ml <- mark_layers(mark)
  .fmt <- function(x) sprintf("%.2f", x)
  ysc <- if (is.null(parents))
    scale_y_continuous(name = "allele frequency", labels = .fmt)
  else
    scale_y_continuous(name = "allele frequency", labels = .fmt,
      sec.axis = sec_axis(~ ., breaks = c(0.04, 0.96),
        labels = c(sprintf("%s\nallele", parents[2]), sprintf("%s\nallele", parents[1]))))
  ggplot(d, aes(g)) +
    chrom_bands(min(d$raw, na.rm = TRUE), max(d$raw, na.rm = TRUE)) +
    sel_bands(regions) +
    geom_point(aes(y = raw), colour = "grey58", size = .18, alpha = .30, stroke = 0) +
    geom_hline(yintercept = .5, colour = "grey35", linewidth = .3) +
    geom_hline(yintercept = c(.5 - th, .5 + th), colour = "#B2182B",
               linewidth = .3, linetype = "42") +
    ml[[1]] +
    geom_line(aes(y = afd, group = chrom), colour = COL_FIT, linewidth = .45) +
    chrom_axis() + ysc + labs(x = NULL) +
    theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.title.y.right = element_blank(),
                        axis.text.y.right = element_text(size = 6.2, colour = "grey25",
                                                         lineheight = .85),
                        axis.ticks.y.right = element_blank(),
                        axis.line.y.right = element_blank())
}
panel_delta <- function(d, th, mark = NULL, regions = MASK,
                        ylab = expression(Delta*" allele frequency"),
                        parents = NULL, colour = COL_DEL, up_in = "treated") {
  ml <- mark_layers(mark)
  lo <- min(d$beta, na.rm = TRUE); hi <- max(d$beta, na.rm = TRUE)
  .fmt <- function(x) sprintf("%.2f", x)
  ysc <- if (is.null(parents))
    scale_y_continuous(name = ylab, labels = .fmt)
  else
    scale_y_continuous(name = ylab, labels = .fmt,
      sec.axis = sec_axis(~ ., breaks = c(lo + .06 * (hi - lo), hi - .06 * (hi - lo)),
        labels = c(sprintf("%s up\nin %s", parents[2], up_in),
                   sprintf("%s up\nin %s", parents[1], up_in))))
  ggplot(d, aes(g)) +
    chrom_bands(min(d$beta, na.rm = TRUE), max(d$beta, na.rm = TRUE)) +
    sel_bands(regions) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = .3) +
    geom_hline(yintercept = c(-th, th), colour = "#B2182B", linewidth = .3, linetype = "42") +
    ml[[1]] +
    geom_line(aes(y = beta, group = chrom), colour = colour, linewidth = .45) +
    chrom_axis() + ysc + labs(x = "chromosome") +
    theme_pub() + theme(axis.title.y.right = element_blank(),
                        axis.text.y.right = element_text(size = 6.2, colour = "grey25",
                                                         lineheight = .85),
                        axis.ticks.y.right = element_blank(),
                        axis.line.y.right = element_blank())
}

## Fitness pair of panels, used by Figures 4 and 5.
af_panel <- function(sample, title = NULL, mark = NULL,
                     ylab_delta = expression(Delta*" allele frequency"),
                     delta_only = FALSE, which = "both") {
  load_sites(sample)
  d  <- get(sample, envir = FIT_CACHE) |> mutate(g = gx(chrom, pos))
  th <- beta_threshold(d)
  .cr  <- main$cross[match(sample, main$Sample_ID)]
  par2 <- if (length(.cr) && !is.na(.cr)) c(sub("xBBQ$", "", .cr), "BBQ") else NULL
  p <- if (identical(which, "af"))
    panel_af(d, th, mark, parents = par2) +
      labs(x = "chromosome") +
      theme(axis.text.x = element_text(size = 6, colour = "black"),
            axis.ticks.x = element_line(colour = "black"))
  else if (identical(which, "delta") || isTRUE(delta_only))
    panel_delta(d, th, mark, colour = COL_FIT, ylab = ylab_delta,
                parents = par2, up_in = "glucose")
  else
    panel_af(d, th, mark) / panel_delta(d, th, mark, colour = COL_FIT,
                                        ylab = ylab_delta) +
      plot_layout(heights = c(1, .85))
  if (!is.null(title))
    p <- p + plot_annotation(title = title,
      subtitle = sprintf("Grey points, per-site allele frequency; blue, the smoothed fit. Dashed red, the smallest effect reaching LOD %g in this pool, %+.3f in allele-frequency units;\nthe standard error carries no overdispersion term, so that threshold is a lower bound. Tan bands, MAT and HO, the two regions the experiment selects on by design. The chromosome VII region is NOT banded: it is a finding of this experiment, not a masked artefact.",
                         LOD_CUT, th),
      theme = theme_pub())
  attr(p, "threshold") <- th
  p
}

## ---------------------------------------------- the stress contrast on demand --
STR_CACHE <- new.env(parent = emptyenv())
PAIR_INDEX <- NULL
load_pair <- function(keys) {
  need <- setdiff(keys, ls(STR_CACHE))
  if (!length(need) && !is.null(PAIR_INDEX)) return(invisible(NULL))
  say("  loading stress contrasts: %s", paste(need, collapse = ", "))
  x <- qs2::qs_read(cache_file("drug_results"))
  if (is.null(PAIR_INDEX))
    PAIR_INDEX <<- bind_rows(lapply(names(x), function(k)
      tibble::tibble(key = k, L = x[[k]]$meta$L, R = x[[k]]$meta$R,
                     condition = as.character(x[[k]]$meta$condition),
                     cross = x[[k]]$meta$cross, well = x[[k]]$meta$well)))
  for (k in need) {
    e <- x[[k]]
    assign(k, list(meta = e$meta,
                   d = tibble::tibble(chrom = e$result$chrom, pos = e$result$physical.position,
                                      beta = e$result$contrast.beta,
                                      beta.se = e$result$contrast.beta.se,
                                      LOD = e$result$LOD) |> filter(!is.na(beta))),
           envir = STR_CACHE)
  }
  rm(x); gc(verbose = FALSE); invisible(NULL)
}

## A treated pool over its own stress contrast — the Figure 2 shape, reused by
## Figure 4 so that a condition-specific locus can be seen surviving the contrast
## while the design-selected regions cancel.
treated_over_contrast <- function(treated, pair, mark, delta_only = FALSE) {
  load_sites(treated); load_pair(pair)
  dt <- get(treated, envir = FIT_CACHE) |> mutate(g = gx(chrom, pos))
  sp <- get(pair, envir = STR_CACHE)
  ds <- sp$d |> mutate(g = gx(chrom, pos))
  tht <- beta_threshold(dt); ths <- beta_threshold(ds)
  kk <- main |> filter(Sample_ID == treated)
  par2 <- c(sub("xBBQ$", "", kk$cross[1]), "BBQ")
  top <- panel_af(dt, tht, mark = mark, regions = MASK, parents = par2) +
    labs(subtitle = sprintf("%s in %s: allele frequency against 0.5",
                            cross_label(kk$cross[1]), unname(LONG[[kk$condition[1]]])))
  bot <- panel_delta(ds, ths, mark = mark, regions = MASK,
                     ylab = delta_lab(unname(LONG[[kk$condition[1]]])),
                     parents = par2, up_in = unname(LONG[[kk$condition[1]]])) +
    labs(subtitle = NULL)
  ## delta_only: Figures 4 draws just the contrast panel. The allele-frequency
  ## panel above it repeats what Figure 2 already established about a single
  ## pool, and the point of Figure 4 is the contrast.
  if (isTRUE(delta_only))
    return(list(plot = bot, th = c(tht, ths)))
  list(plot = (top / bot) + plot_layout(heights = c(1, .85)), th = c(tht, ths))
}

## =================================================================== FIG 2 ==
## A, allele frequency in a glucose pool.  B, the paired stress contrast from the
## SAME well: treated minus its own glucose control.  MAT, HO and the chromosome
## VII sweep are carried by both members of the pair and cancel in B, which is the
## internal control; the script checks that and prints what it finds.
say("Figure 2")
FIG2 <- list(list(glu = "Sample302", pair = "Sample014_vs_Sample302"),
             list(glu = "Sample303", pair = "Sample207_vs_Sample303"),
             list(glu = "Sample302", pair = "Sample206_vs_Sample302"))
## CANCEL is the set of regions the paired contrast should cancel, and it is
## used for the ARITHMETIC CHECK below -- not for shading. chrVII is a finding,
## not an artefact: the text says only MAT and HO are masked, so banding chrVII
## tan tells the reader it was masked and contradicts the Results. Bands come
## from MASK alone, everywhere.
CANCEL <- bind_rows(MASK, CHR7)

fig2_build <- function(glu, pair, stem) {
  load_sites(glu); load_pair(pair)
  dg <- get(glu, envir = FIT_CACHE) |> mutate(g = gx(chrom, pos))
  sp <- get(pair, envir = STR_CACHE)
  ds <- sp$d |> mutate(g = gx(chrom, pos))
  thg <- beta_threshold(dg); ths <- beta_threshold(ds)

  ## does the contrast actually cancel the design-selected regions?
  chk <- CANCEL |> rowwise() |> mutate(
    fit_max = max(abs(dg$beta[dg$chrom == chrom & dg$pos >= lo & dg$pos <= hi]), na.rm = TRUE),
    str_max = max(abs(ds$beta[ds$chrom == chrom & ds$pos >= lo & ds$pos <= hi]), na.rm = TRUE)) |>
    ungroup()
  for (i in seq_len(nrow(chk)))
    say("    %-12s fitness |dAF| up to %.3f, stress contrast up to %.3f  (%.0f-fold lower%s)",
        chk$name[i], chk$fit_max[i], chk$str_max[i], chk$fit_max[i] / chk$str_max[i],
        if (chk$str_max[i] < ths) ", below this pair's threshold" else ", small residual above it")

  kg   <- main |> filter(Sample_ID == glu)
  wild <- sub("xBBQ$", "", kg$cross[1])
  par2 <- c(wild, "BBQ")

  ## 2026-10-05: Figure 2 is the allele-frequency panel alone. The paired
  ## contrast that used to sit under it as B is gone; the cancellation check
  ## above still prints, so the fact stays available to the text.
  fig <- panel_af(dg, thg, regions = MASK, parents = par2) +
    labs(x = "chromosome") +
    theme(axis.text.x = element_text(size = 6, colour = "black"),
          axis.ticks.x = element_line(colour = "black"))

  cond <- unname(LONG[[as.character(sp$meta$condition)]])
  write_legend(sprintf("%s_legend.txt", stem),
    "FIGURE 2.  What one pool yields.",
    "",
    wrap_para(sprintf(
      "Allele frequency across the genome in pool %s (cross %s, glucose, well %s). Grey points are the per-site estimates and the blue line the smoothed fit, estimated separately for each chromosome and drawn with a break at every chromosome boundary. The allele frequency is that of the wild parent %s, so the trace runs to 1 where the %s allele is fixed and to 0 where the BBQ allele is fixed. The dashed red lines mark the smallest effect that reaches LOD %g in this pool, %+.3f in allele-frequency units; because the standard error carries no overdispersion term, that threshold is a lower bound. Tan bands mark the mating-type locus (MAT, chromosome III) and the HO locus (chromosome IV), which the experiment selects on by design, so both run to the BBQ allele.",
      glu, kg$cross[1], kg$Wells[1], wild, wild, LOD_CUT, thg)),
    "",
    sprintf("Source: run v7.  Built by draft/figures/make_figures.R on %s.", Sys.Date()))
  fig
}
## Only the first pairing is built: with B gone the alternates differed only in
## a panel that no longer exists.
p <- fig2_build(FIG2[[1]]$glu, FIG2[[1]]$pair, "Figure2")
save_fig(p, "Figure2.png", 7.5, 2.7)

## =================================================================== FIG 3 ==
## protoC reworked: colour is the effect size of the strongest peak covering the
## bin, not the share of pools.  LOD is deliberately not the colour scale — it is
## nominal in this pipeline and would not be comparable between pools.
## chromosome VII is RETAINED in every panel (the user's call, 2026-09-18).
## Every bin in the genome gets a tile, so a bin with no peak reads as its own
## colour rather than as blank paper.
say("Figure 3")
NOPEAK <- "grey93"

genome_bins <- function(bin) {
  bind_rows(lapply(CHROMS, function(c) {
    st <- seq(0, CHRLEN[[c]] - 1, by = bin)
    tibble::tibble(seqnames = c, bin_id = st / bin, gmid = gx(c, st + bin / 2))
  }))
}
bin_beta <- function(pk, bin) {
  pk |>
    rowwise() |>
    mutate(b0 = floor(start / bin), b1 = floor(end / bin)) |>
    ungroup() |>
    mutate(bin_id = purrr::map2(b0, b1, seq)) |>
    tidyr::unnest(bin_id) |>
    group_by(condition, seqnames, bin_id) |>
    ## V7: n_distinct(cross), not max|beta|. Colouring by the largest effect in
    ## a bin makes every chromosome light up, because the union of 1,472
    ## intervals covers the genome 5.15 times over -- the panel was true by
    ## construction and said nothing. Counting the CROSSES with a peak in the
    ## bin asks the question the Results actually ask, and makes the chromosome
    ## VII shared-parent block legible as what it is.
    summarise(n_cross = n_distinct(cross), max_beta = max(abs(peak.beta)),
              .groups = "drop")
}
fig3_variant <- function(bin, lod_min, subtitle, palette = "mako", bands = TRUE,
                         encoding = c("crosses", "beta")) {
  encoding <- match.arg(encoding)
  d <- apply_mask(fit_main |> filter(peak.LOD >= lod_min), MASK) |>
    mutate(condition = factor(unname(LONG[condition]), unname(LONG[CONDS])))
  grid <- tidyr::expand_grid(condition = factor(unname(LONG[CONDS]), unname(LONG[CONDS])),
                             genome_bins(bin))
  b <- grid |> left_join(bin_beta(d, bin), by = c("condition", "seqnames", "bin_id"))
  LEGEND_NAME <- if (encoding == "crosses") "crosses with a peak in the bin"
                 else expression("max |AF - 0.5| in bin")
  scale_fn <- if (palette == "BuGn")
    scale_fill_distiller(palette = "BuGn", direction = 1, na.value = NOPEAK,
                         name = LEGEND_NAME,
                         guide = guide_colourbar(barwidth = 9, barheight = .4))
  else
    scale_fill_viridis_c(option = "G", direction = -1, na.value = NOPEAK,
                         name = LEGEND_NAME,
                         guide = guide_colourbar(barwidth = 9, barheight = .4))
  mk <- MASK |> mutate(g0 = gx(chrom, lo), g1 = gx(chrom, hi))
  band_layer <- if (bands)
    geom_rect(data = mk, aes(xmin = g0, xmax = g1, ymin = -Inf, ymax = Inf),
              fill = "#d9c9a3", colour = NA, inherit.aes = FALSE) else NULL
  fill_var <- if (encoding == "crosses") "n_cross" else "max_beta"
  ggplot(b) +
    geom_tile(aes(gmid, condition, fill = .data[[fill_var]]), height = .82, width = bin) +
    band_layer +
    scale_fn + chrom_axis() +
    scale_y_discrete(limits = rev(unname(LONG[CONDS]))) +
    labs(x = "chromosome", y = NULL, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom",
                        axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
## Chosen 2026-09-18: mako, 50 kb bins, peaks at LOD >= 5, no censor bands.
## Set MAKE_VARIANTS to TRUE to re-emit the bin-size / threshold comparison sheet.
MAKE_VARIANTS <- FALSE
VARIANTS <- list(list(20000,  5, "20 kb bins, peaks at LOD >= 5"),
                 list(50000,  5, "50 kb bins, peaks at LOD >= 5"),
                 list(20000, 20, "20 kb bins, peaks at LOD >= 20"),
                 list(50000, 20, "50 kb bins, peaks at LOD >= 20"))

## 2026-10-05: Figure 3 is CONDITION-SPECIFIC. It is drawn from the paired
## stress-contrast peaks (stress pool minus its own glucose control), so a bin
## lights up only where a stress shifts allele frequency relative to glucose.
## Glucose is the reference of every contrast and so has no row. The fitness
## version above (fig3_variant) is kept as a function but no longer emitted.
## Colour encodings on offer, pending the user's choice:
##   mako  -- crosses with a contrast peak in the bin, one sequential scale
##   hue   -- the same count, each row in its own condition colour
##   dir   -- signed: crosses where the wild allele rises minus crosses where
##            BBQ rises; a cross counts once, in the direction of its largest
##            peak over the bin
FIG3_BIN <- 50000; FIG3_ENC <- "mako"     # the default until one is chosen
SCOND <- c("Ace", "Caf", "NaCl")
str_bins <- function(bin) {
  d <- str_keep |> filter(condition %in% SCOND) |>
    mutate(condition = factor(unname(LONG[condition]), unname(LONG[SCOND])),
           b0 = floor(start / bin), b1 = floor(end / bin))
  d <- d |> mutate(bin_id = purrr::map2(b0, b1, seq)) |> tidyr::unnest(bin_id) |>
    group_by(condition, seqnames, bin_id, cross) |>
    slice_max(abs(peak.beta), n = 1, with_ties = FALSE) |>
    group_by(condition, seqnames, bin_id) |>
    summarise(n_cross = n_distinct(cross),
              net = sum(peak.beta > 0) - sum(peak.beta < 0), .groups = "drop")
  tidyr::expand_grid(condition = factor(unname(LONG[SCOND]), unname(LONG[SCOND])),
                     genome_bins(bin)) |>
    left_join(d, by = c("condition", "seqnames", "bin_id"))
}
fig3_contrast <- function(bin, enc = "mako", subtitle = NULL) {
  b <- str_bins(bin)
  base <- ggplot(b) + chrom_axis() +
    scale_y_discrete(limits = rev(unname(LONG[SCOND]))) +
    labs(x = "chromosome", y = NULL, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom", axis.line.y = element_blank(),
                        axis.ticks.y = element_blank())
  gb <- guide_colourbar(barwidth = 9, barheight = .4)
  if (enc == "mako")
    base + geom_tile(aes(gmid, condition, fill = n_cross), height = .82, width = bin) +
      scale_fill_viridis_c(option = "G", direction = -1, na.value = NOPEAK,
                           name = "crosses with a condition-specific peak", guide = gb)
  else if (enc == "hue")
    base + geom_tile(aes(gmid, condition), fill = NOPEAK, height = .82, width = bin) +
      geom_tile(data = b |> filter(!is.na(n_cross)),
                aes(gmid, condition, fill = condition, alpha = n_cross),
                height = .82, width = bin) +
      scale_fill_manual(values = CCOL, guide = "none") +
      scale_alpha_continuous(range = c(.25, 1), name = "crosses with a condition-specific peak")
  else {
    lim <- max(abs(b$net), na.rm = TRUE)
    base + geom_tile(aes(gmid, condition, fill = net), height = .82, width = bin) +
      scale_fill_distiller(palette = "PuOr", direction = 1, limits = c(-lim, lim),
                           na.value = NOPEAK, guide = gb,
                           name = "crosses, wild allele up minus BBQ allele up")
  }
}
fig3 <- fig3_contrast(FIG3_BIN, FIG3_ENC)
save_fig(fig3, "Figure3.png", 7, 2.6)
## The option sheet, 2026-10-06: bins fixed at 50 kb; colour schemes only.
## The first sheet's light-green low end was hard to tell from the no-peak grey,
## so every scheme here either starts its ramp well away from the light end, or
## changes the no-peak colour instead. Each is also written alone at 600 dpi.
fig3_scheme <- function(fill_scale, nopeak = NOPEAK, subtitle, stepped = FALSE) {
  b <- str_bins(50000)
  if (stepped) b <- b |> mutate(n_cross = factor(ifelse(n_cross >= 5, "5+", n_cross),
                                                 c("1", "2", "3", "4", "5+")))
  ggplot(b) +
    geom_tile(aes(gmid, condition), fill = nopeak, height = .82, width = 50000) +
    geom_tile(data = b |> filter(!is.na(n_cross)),
              aes(gmid, condition, fill = n_cross), height = .82, width = 50000) +
    fill_scale + chrom_axis() +
    scale_y_discrete(limits = rev(unname(LONG[SCOND]))) +
    labs(x = "chromosome", y = NULL, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom", axis.line.y = element_blank(),
                        axis.ticks.y = element_blank())
}
LEG <- "crosses with a condition-specific peak"
gb  <- guide_colourbar(barwidth = 9, barheight = .4)
vir <- function(opt, begin, end) scale_fill_viridis_c(option = opt, direction = -1,
         begin = begin, end = end, name = LEG, guide = gb)
step5 <- function(cols) scale_fill_manual(values = setNames(cols, c("1","2","3","4","5+")),
         name = LEG, drop = FALSE)
SCHEMES <- list(
  list("mako_trimmed",  vir("G", .05, .70), NOPEAK,  "A. mako, light end trimmed off; no peak grey"),
  list("mako_white",    vir("G", 0, .90),   "white", "B. mako, full ramp; no peak WHITE instead of grey"),
  list("viridis",       vir("D", 0, .85),   NOPEAK,  "C. viridis, yellow (1) to purple (many)"),
  list("magma",         vir("A", .10, .80), NOPEAK,  "D. magma, orange (1) to deep purple (many)"),
  list("rocket",        vir("F", .10, .75), NOPEAK,  "E. rocket, salmon (1) to dark red (many)"),
  list("cividis",       vir("E", 0, .85),   NOPEAK,  "F. cividis, khaki (1) to navy (many); safest for colour-blind readers"),
  list("YlOrRd_steps",  step5(c("#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#800026")), NOPEAK,
       "G. YlOrRd, five discrete steps (1, 2, 3, 4, 5+)", TRUE),
  list("Blues_steps",   step5(c("#9ECAE1","#6BAED6","#3182BD","#08519C","#08306B")), "grey97",
       "H. Blues, five discrete steps; no peak near-white", TRUE))
opts <- lapply(SCHEMES, function(sc) {
  p <- fig3_scheme(sc[[2]], sc[[3]], sc[[4]], stepped = length(sc) > 4)
  save_fig(p, sprintf("Figure3_option_%s.png", sc[[1]]), 7, 2.6, dpi = 600)
  p
})
save_fig(wrap_plots(opts, ncol = 1), "Figure3_OPTIONS.png", 7, 2.6 * length(opts), dpi = 400)
say("  contrast peaks binned: %d in %d crosses",
    nrow(str_keep |> filter(condition %in% SCOND)), n_distinct(str_keep$cross))
write_legend("Figure3_legend.txt",
  "FIGURE 3.  Where each stress moves allele frequency away from glucose.",
  "",
  wrap_para(sprintf("The genome is divided into %d kb bins and each bin is coloured by the number of DISTINCT CROSSES with a condition-specific peak whose interval covers it. A condition-specific peak is one called on the paired stress contrast: a stress pool's allele frequency minus that of its own glucose control from the same well, so loci that act equally in glucose and in the stress cancel and do not appear. Glucose is the reference of every contrast and therefore has no row. Bins covered by no peak are light grey, which is distinct from the low end of the scale. Peaks inside the mating-type region and the HO locus are removed before binning; both cancel in the contrast in any case. The cross is the unit because replicate wells of one cross share parents and segregant material.", FIG3_BIN / 1000)),
  "",
  wrap_para(sprintf("Conditions are ordered acetate, caffeine, sodium chloride. Drawn from %d condition-specific peaks after the mating-type and HO mask: acetate %d, caffeine %d, sodium chloride %d.",
    nrow(str_keep |> filter(condition %in% SCOND)),
    sum(str_keep$condition == "Ace"), sum(str_keep$condition == "Caf"), sum(str_keep$condition == "NaCl"))),
  "",
  sprintf("Source: run v7.  Built by draft/figures/make_figures.R on %s.", Sys.Date()))

## =================================================================== FIG 4 ==
say("Figure 4")
pub_at <- function(pattern) pub |> filter(grepl(pattern, qtl)) |> slice(1)
ena  <- pub_at("^NaCl\\.chr04\\.1$")
whi2 <- pub_at("^Caff\\.chr15\\.3$")
tor1 <- pub_at("^Caff\\.chr10\\.1$")   # the caffeine hallmark; reported even though it fails

## The pool drawn must cover the published position AND have a matched glucose
## control, so the contrast panel exists.  Deepest such pool, not the most extreme:
## shallow pools inflate |beta|.
load_pair(character(0))                       # builds PAIR_INDEX
depth_of <- setNames(key$mean_depth_seg, key$Sample_ID)
pool_over <- function(cond, chrom, pos) {
  cand <- fit_keep |>
    filter(condition == cond, seqnames == chrom, start <= pos, end >= pos) |>
    distinct(sample) |>
    inner_join(PAIR_INDEX |> filter(condition == cond), by = c("sample" = "L")) |>
    mutate(dep = unname(depth_of[sample])) |>
    arrange(desc(dep))
  if (!nrow(cand)) return(NULL)
  list(sample = cand$sample[1], pair = cand$key[1], dep = cand$dep[1])
}
pool_over_dir <- function(cond, chrom, pos, want) {
  cand <- fit_keep |>
    filter(condition == cond, seqnames == chrom, start <= pos, end >= pos,
           if (want == "BBQ") peak.beta < 0 else peak.beta > 0) |>
    group_by(sample) |> slice_max(abs(peak.beta), n = 1, with_ties = FALSE) |> ungroup() |>
    ## keep only what is needed before the join: PAIR_INDEX carries its own
    ## `cross` and `condition`, which would otherwise be suffixed away
    select(sample, peak.beta, our_cross = cross) |>
    inner_join(PAIR_INDEX |> filter(condition == cond), by = c("sample" = "L")) |>
    mutate(dep = unname(depth_of[sample])) |>
    ## DEEPEST, not most extreme. Shallow pools inflate |beta|, so ranking on
    ## effect size systematically draws the worst-supported pool -- the first
    ## pass here picked a 5x pool over a 26x one. Require the locus to be
    ## visible at all, then take the best-sequenced pool that shows it.
    filter(!is.na(dep), abs(peak.beta) >= 0.10) |> arrange(desc(dep))
  if (!nrow(cand)) return(NULL)
  list(sample = cand$sample[1], pair = cand$key[1], dep = cand$dep[1],
       beta = cand$peak.beta[1], cross = cand$our_cross[1])
}
e_sel <- pool_over_dir("NaCl", ena$chrom, ena$maxpos, "BBQ")
w_sel <- pool_over_dir("NaCl", ena$chrom, ena$maxpos, "wild")
if (is.null(w_sel)) { w_sel <- pool_over("Caf", whi2$chrom, whi2$maxpos)
                      say("  !! no wild-direction ENA pool with a control; falling back to WHI2") }
say("  ENA toward BBQ : %s (%s) beta %+.3f at %.0fx", e_sel$sample, e_sel$cross, e_sel$beta, e_sel$dep)
say("  ENA toward wild: %s (%s) beta %+.3f at %.0fx", w_sel$sample, w_sel$cross, w_sel$beta, w_sel$dep)
say("  ENA:  %s (%.0fx) vs its control, pair %s", e_sel$sample, e_sel$dep, e_sel$pair)
say("  WHI2: %s (%.0fx) vs its control, pair %s", w_sel$sample, w_sel$dep, w_sel$pair)

a4 <- treated_over_contrast(e_sel$sample, e_sel$pair, delta_only = TRUE,
        tibble::tibble(chrom = ena$chrom, pos = ena$maxpos, label = "ENA"))
## B is ENA again, in a cross where the WILD allele wins. Same locus, opposite
## direction: that is the panel, not an accident of pool choice.
b4 <- treated_over_contrast(w_sel$sample, w_sel$pair, delta_only = TRUE,
        tibble::tibble(chrom = ena$chrom, pos = ena$maxpos, label = "ENA"))
p4a <- a4$plot; p4b <- b4$plot

## 2026-10-05: panel C (the condition-matched test and the recurrence bars) is
## deleted. Those statistics live in the text and Table form now.
fig4 <- (wrap_elements(p4a) / wrap_elements(p4b)) +
  plot_annotation(tag_levels = "A", theme = theme_pub())
save_fig(fig4, "Figure4.png", 7, 5.4)
## The legend MUST read from the same cell the panel draws, or the figure
## contradicts its own caption. Both are contrast peaks, unit = cross; the
## by-pool figure is quoted alongside because the gap is the argument.
enr_row  <- function(q) CELL(enr_s)               |> filter(pool_set == "matched", qtl == q)
enr_pool <- function(q) CELL(enr_s, un = "pool")  |> filter(pool_set == "matched", qtl == q)
write_legend("Figure4_legend.txt",
  "FIGURE 4.  ENA on the paired contrast, in both directions.",
  "",
  wrap_para(sprintf("(A) The ENA locus on chromosome IV at %s, in pool %s (%s, sodium chloride, %.0fx), shown as that pool's own stress contrast against its matched glucose control %s from the same well. (B) The SAME LOCUS in a cross where the wild allele wins instead: pool %s (%s, sodium chloride, %.0fx) against %s, where ENA shifts %+.2f rather than %+.2f. Treusch et al. report ENA as multi-allelic, and the two panels are that: one locus pulling in opposite directions in different crosses. Each panel is the pool's paired stress contrast; a condition-specific locus survives the subtraction while the mating-type region and HO, banded in tan, cancel, as does the chromosome VII region, which is not banded because it is reported as a finding rather than masked. Dashed red lines mark the smallest effect reaching LOD %g in that pool.",
    format(round(ena$maxpos), big.mark = ","), e_sel$sample, main$cross[match(e_sel$sample, main$Sample_ID)], e_sel$dep, PAIR_INDEX$R[match(e_sel$pair, PAIR_INDEX$key)],
    w_sel$sample, w_sel$cross, w_sel$dep, PAIR_INDEX$R[match(w_sel$pair, PAIR_INDEX$key)],
    w_sel$beta, e_sel$beta, LOD_CUT)),
  "",
  sprintf("Source: run v7.  Built by draft/figures/make_figures.R on %s.", Sys.Date()))

## =================================================================== FIG 5 ==
## 2026-10-05: Figure 5 is DROPPED. A search of the contrast peaks for a large,
## multi-cross, condition-specific locus (acetate especially) found none worth a
## figure. The block is kept, switched off, in case that changes.
MAKE_FIG5 <- FALSE
if (MAKE_FIG5) {
## Largest-effect locus per condition with no published counterpart.  Treusch et al.
## assayed only caffeine and NaCl, so glucose and acetate have NO COMPARATOR —
## which is not the same as having no counterpart, and the labels say so.
say("Figure 5")
loci_for <- function(cond) {
  ## The panel for a stress condition shows the PAIRED CONTRAST, so the locus
  ## must be called in the contrast too. Selecting from the fitness peaks and
  ## then drawing the contrast put the acetate locus on a panel where it had
  ## zero contrast peaks and so could not reach the threshold by construction.
  ## Glucose is the reference arm and has no contrast, so it keeps fitness.
  d <- if (cond == "Glu") fit_keep |> filter(condition == cond) else
       str_keep |> filter(condition == cond) |>
         mutate(sample = paste(cross, well, condition, sep = "_"))
  gr <- reduce(GRanges(d$seqnames, IRanges(pmax(d$start, 1), d$end)))
  hits <- findOverlaps(GRanges(d$seqnames, IRanges(d$peak.pos, d$peak.pos)), gr)
  d$locus <- subjectHits(hits)[match(seq_len(nrow(d)), queryHits(hits))]
  ## "no published locus explains this" has to be tested against the published
  ## INTERVALS and against ALL of them, not against the max-LOD position of the
  ## matching condition only. Two failures of the narrower test, both caught
  ## 2026-09-22: the largest unmatched acetate locus was WHI2 -- Treusch never
  ## assayed acetate, so nothing in that condition could ever match -- and the
  ## NaCl one sat 11 kb outside a published interval while containing no
  ## published summit. Glucose and acetate have no published set of their own,
  ## so under the narrow test every locus in them was "unmatched" by
  ## construction, which is not a finding.
  covered <- unique(subjectHits(findOverlaps(
    GRanges(pub$chrom, IRanges(pub$lo, pub$hi)), gr)))
  ## Shallow pools inflate |beta| and call the most peaks, so the locus is ranked
  ## on the median effect across the pools that detect it, and the pool drawn is
  ## the deepest of them rather than the most extreme.
  depth <- setNames(key$mean_depth_seg, key$Sample_ID)
  ## depth belongs to the treated member of a contrast pair
  d <- d |> mutate(.id = if (cond == "Glu") sample else L)
  d |> mutate(dep = unname(depth[.id])) |>
    group_by(locus) |>
    summarise(chrom = dplyr::first(seqnames), pos = peak.pos[which.max(abs(peak.beta))],
              med_beta = median(abs(peak.beta)), max_beta = max(abs(peak.beta)),
              pools = n_distinct(sample), crosses = n_distinct(cross),
              sample = sample[which.max(ifelse(is.na(dep), -1, dep))],
              ## rep_depth must describe the pool that will actually be DRAWN,
              ## which is one covering `pos`. Taking the maximum over the whole
              ## merged block instead let a locus through on the strength of a
              ## deep pool that does not cover the marked position at all: the
              ## acetate panel was filtered on a 28x pool and then drawn from
              ## the only pool covering its locus, at 8x.
              ## NB peak.pos, not pos: `pos` is created earlier in this same
              ## summarise() and is already a scalar, so indexing it again
              ## returned NA and made rep_depth -Inf for nearly every locus.
              rep_depth = suppressWarnings(max(
                dep[start <= peak.pos[which.max(abs(peak.beta))] &
                    end   >= peak.pos[which.max(abs(peak.beta))]], na.rm = TRUE)),
              pools_at_pos = sum(start <= peak.pos[which.max(abs(peak.beta))] &
                                 end   >= peak.pos[which.max(abs(peak.beta))]),
              .groups = "drop") |>
    mutate(condition = cond, matched = locus %in% covered,
           in_sweep = chrom == CHR7$chrom & pos >= CHR7$lo & pos <= CHR7$hi)
}
all_loci <- bind_rows(lapply(CONDS, loci_for))
## SELECTION, tightened. A panel exists to be looked at, so a locus qualifies
## only if some pool can actually show it:
##   rep_depth >= 15  the drawn pool is deep enough that the trace is a trace.
##                    NO cross threshold: a locus seen in one cross may be a
##                    variant private to that parent, which is exactly what is
##                    worth following up, so requiring three crosses was
##                    filtering out the interesting case.
##   not subtelomeric within TELO_BP of a chromosome end, where marker density
##                    falls away and the loess runs out of data on one side --
##                    three of the four previously drawn loci were there
TELO_BP <- 75000   # 35 kb passed a locus 44 kb from a chromosome end
all_loci <- all_loci |>
  mutate(telomeric = pos < TELO_BP | pos > (unname(CHRLEN[as.character(chrom)]) - TELO_BP))
## Spread the panels across the genome. The four top-ranked loci happened to be
## three on chromosome IV and one elsewhere, which reads as a statement about
## chromosome IV rather than about four conditions. Each condition takes its
## best locus on a chromosome no earlier panel has used; if every eligible
## locus for a condition is on an already-used chromosome it takes its best
## anyway rather than going empty, and the legend names the rank so the cost of
## the constraint is visible.
.elig <- all_loci |> filter(!matched, !in_sweep, rep_depth >= 15, !telomeric) |>
  mutate(condition = factor(condition, CONDS)) |>
  arrange(condition, desc(med_beta)) |>
  group_by(condition) |> mutate(rank_in_cond = row_number()) |> ungroup()
.used <- character(0); .picks <- list()
for (.cd in CONDS) {
  .c  <- .elig |> filter(as.character(condition) == .cd)
  if (!nrow(.c)) next
  .free <- .c |> filter(!(as.character(chrom) %in% .used))
  .row <- if (nrow(.free)) .free[1, ] else .c[1, ]
  .used <- c(.used, as.character(.row$chrom))
  .picks[[.cd]] <- .row
}
top_un <- bind_rows(.picks) |> mutate(condition = factor(condition, CONDS)) |> arrange(condition)
for (i in seq_len(nrow(top_un)))
  say("  %-9s %s %.0f kb  (rank %d of %d for this condition)",
      unname(LONG[[as.character(top_un$condition[i])]]), top_un$chrom[i],
      top_un$pos[i]/1000, top_un$rank_in_cond[i],
      sum(as.character(.elig$condition) == as.character(top_un$condition[i])))
load_pair(character(0))
panels <- lapply(seq_len(nrow(top_un)), function(i) {
  r    <- top_un[i, ]
  cond <- as.character(r$condition)
  comparator <- if (cond %in% c("Caf", "NaCl")) "no published counterpart" else
                "no published comparator for this condition"
  mk <- tibble::tibble(chrom = r$chrom, pos = r$pos,
                       label = sprintf("%s %.0f kb — %s", r$chrom, r$pos/1000, comparator))
  ## pools detecting this locus, deepest first; prefer one with a matched control
  cand <- fit_keep |>
    filter(condition == cond, seqnames == r$chrom, start <= r$pos, end >= r$pos) |>
    distinct(sample) |> mutate(dep = unname(depth_of[sample])) |> arrange(desc(dep))
  paired <- cand |> inner_join(PAIR_INDEX |> filter(condition == cond),
                               by = c("sample" = "L")) |> arrange(desc(dep))
  if (cond != "Glu" && nrow(paired)) {
    pick <- paired[1, ]
    body <- treated_over_contrast(pick$sample, pick$key, mk, delta_only = TRUE)$plot
    note <- sprintf("lower panel, %s minus %s, its own glucose control", pick$sample, pick$R)
    drawn <- pick$sample; dep <- pick$dep
  } else {
    pick <- cand[1, ]
    ## Glucose: ONE panel. Its allele-frequency trace and its fitness contrast
    ## are the same series shifted by 0.5, so showing both said one thing twice.
    body <- af_panel(pick$sample, mark = mk, which = "af")
    note <- if (cond == "Glu") "lower panel, the fitness contrast against 0.5 — glucose is the reference arm and has no stress contrast"
            else "lower panel, the fitness contrast against 0.5 — no matched glucose control for any pool detecting this locus"
    drawn <- pick$sample; dep <- pick$dep
  }
  ttl <- sprintf("%s — %s %.0f kb", unname(LONG[[cond]]), r$chrom, r$pos/1000)
  body + labs(subtitle = ttl) +
    theme(plot.subtitle = element_text(face = "bold", size = 8, colour = "grey15",
                                       margin = margin(b = 2)))
})
## Panel E, rebuilt. Three things were wrong with it.
##  1. The axis said |dAF|, which in every other panel means treated minus
##     glucose. These are FITNESS peaks, so the quantity is |AF - 0.5|. The two
##     were being given one name.
##  2. It plotted EVERY locus while A-D are filtered to >= 3 crosses, adequate
##     depth and non-telomeric, so the largest point in each row was a locus
##     A-D had deliberately excluded and the panel contradicted the ones above
##     it. It now draws the same filtered set, with the excluded loci behind in
##     pale grey so nothing is hidden.
##  3. Glucose and acetate showed as all-unmatched, which reads as "we found
##     nothing published there" when the truth is that nothing published exists
##     to match. They are now a third category rather than a colour absence.
sm <- all_loci |>
  mutate(condition = factor(unname(LONG[condition]), unname(LONG[CONDS])),
         eligible = !in_sweep & rep_depth >= 15 & !telomeric,
         status = dplyr::case_when(
           !as.character(condition) %in% unname(LONG[c("Caf", "NaCl")]) ~ "no published set for this condition",
           matched ~ "matches a published locus",
           TRUE    ~ "no published counterpart"))
E_LV <- c("matches a published locus", "no published counterpart",
          "no published set for this condition",
          "below the threshold used for A-D")
sm <- sm |> mutate(key = factor(ifelse(eligible, status, E_LV[4]), E_LV))
p5e <- ggplot() +
  geom_point(data = sm |> filter(!eligible),
             aes(med_beta, condition, colour = key),
             position = position_jitter(height = .18), size = .6, alpha = .30) +
  geom_point(data = sm |> filter(eligible),
             aes(med_beta, condition, colour = key),
             position = position_jitter(height = .18), size = 1.3, alpha = .85) +
  scale_colour_manual(values = setNames(c("#CC79A7", "#0072B2", "#009E73", "grey78"), E_LV),
                      drop = FALSE, name = NULL,
                      guide = guide_legend(override.aes = list(size = 2, alpha = 1),
                                           ncol = 1)) +
  scale_y_discrete(limits = rev(unname(LONG[CONDS]))) +
  labs(x = "median |AF - 0.5| across the pools detecting the locus", y = NULL) +
  theme_pub() + theme(legend.position = "inside",
                      legend.position.inside = c(.99, .55),
                      legend.justification = c(1, .5),
                      legend.background = element_rect(fill = alpha("white", .8), colour = NA),
                      legend.key = element_blank(),
                      legend.text = element_text(size = 6.5),
                      legend.key.size = unit(2.8, "mm"))
## V7: GO panel dropped, as in Figure 4 -- the test is null everywhere.
MISSING <- setdiff(CONDS, as.character(top_un$condition))
for (cd in MISSING)
  say("  !! %s: no locus meets the criteria; drawing a labelled blank", unname(LONG[[cd]]))
blank_panel <- function(cd)
  ggplot() + theme_void() +
    annotate("text", 0, 0, size = 2.9, colour = "grey35", lineheight = 1.2,
      label = sprintf("%s: no qualifying locus.\nEvery candidate is covered only by pools too shallow to draw.",
                      unname(LONG[[cd]]))) +
    theme(panel.border = element_rect(colour = "grey85", fill = NA, linetype = "22"))
ordered <- lapply(CONDS, function(cd) {
  i <- match(cd, as.character(top_un$condition))
  if (is.na(i)) blank_panel(cd) else panels[[i]]
})
H <- vapply(CONDS, function(cd) if (cd %in% MISSING) 0.30 else 1, 0)
fig5 <- (ordered[[1]] / ordered[[2]] / ordered[[3]] / ordered[[4]] / p5e) +
  plot_layout(heights = c(H, 1.25)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E")),
    title = "Figure 5.  The largest-effect locus of each condition that no published locus explains",
    theme = theme_pub())
save_fig(fig5, "Figure5.png", 7, 9.5)
write_legend("Figure5_legend.txt",
  "FIGURE 5.  The largest-effect locus of each condition that no published locus explains.",
  "",
  wrap_para("WHY THE DASHED LINES SIT AT DIFFERENT HEIGHTS IN EACH PANEL. They are not a fixed significance level. Each pair marks the smallest allele-frequency shift that reaches LOD 5 IN THAT POOL, and the LOD a given shift earns depends on that pool's own standard error, which in turn depends on its sequencing depth. A deeply sequenced pool therefore has a tighter threshold than a shallow one, and the same excursion can clear the line in one panel and not in another. The thresholds are comparable within a panel, never between them. They are also a lower bound: the standard error is a flat per-pool scalar with no overdispersion term, and the effective pool size measured from the traces is about a quarter of the value the error model assumes."),
  "",
  if (length(MISSING))
    wrap_para(sprintf("NO PANEL FOR %s. Every locus of that condition meeting the criteria below is covered only by pools too shallow to draw, and a panel drawn from such a trace would be showing noise and labelling it a locus, so the slot is left explicitly empty. That is a statement about the depth of that arm, not about whether it has loci.",
                      paste(LONG[MISSING], collapse = " and ")))
  else
    wrap_para("EVERY CONDITION HAS A PANEL. Acetate received one selective passage rather than two, because its 96 h cultures did not survive, and it is the shallowest-covered arm after caffeine; its panel rests on a single pool and should be read with that in mind."),  "",
  "",
  wrap_para("(A) glucose, (B) acetate, (C) caffeine, (D) sodium chloride. Peaks were grouped into loci by merging overlapping intervals within a condition. A locus qualifies here if it overlaps NO published interval from either mapped condition, lies outside the masked regions and the chromosome VII shared-parent block, sits more than 75 kb from a chromosome end, and is covered by at least one pool sequenced deeply enough to draw. There is deliberately no threshold on the number of crosses: a locus seen in a single cross may rest on a variant private to that parent, which is the case most worth following up, and requiring several crosses removed exactly those. Loci are ranked on the median absolute allele-frequency shift across the pools that detect them rather than the maximum, and the pool drawn is the deepest of those pools rather than the most extreme, because shallow pools inflate both the effect size and the peak count. Each stress panel shows one pool's paired contrast against its own glucose control from the same well. Glucose is the reference arm of that contrast and so has none; its single panel is the fitness contrast against 0.5, which is its allele-frequency trace shifted by a half, and is labelled as such. Loci for the stress conditions are called on the contrast peaks rather than the fitness peaks, so the locus drawn is one the displayed panel actually contains; selecting on fitness and drawing the contrast previously put a locus on a panel that held no contrast peak for it at all."),
  "",
  wrap_para("Each panel draws its condition's highest-ranked qualifying locus on a chromosome no earlier panel has used. Without that constraint three of the four fell on chromosome IV, which reads as a claim about chromosome IV rather than about four conditions. The cost is small and is stated per locus below: glucose and acetate take their top-ranked locus, caffeine and sodium chloride their second, at 0.183 against 0.193 and 0.163 against 0.181 respectively."),
  "",
  wrap_para("CAVEAT ON THE GLUCOSE PANEL. The chromosome VIII locus carries the largest shift in the figure, 0.50, but the trace immediately around it is noisy rather than a smooth selected haplotype, and the locus rests on five pools from a single cross. Read it as a candidate to check rather than as an established locus; a cleaner alternative at a smaller effect is available if it does not survive scrutiny."),
  "",
  wrap_para("The dashed vertical line in each panel marks the locus named in that panel's heading. No published locus overlaps any of them: for caffeine and sodium chloride that means no interval from Treusch et al. covers the position, and for glucose and acetate it means no comparable published set exists at all, since neither condition was assayed there."),
  "",
  "Loci drawn:",
  paste0("  ", apply(top_un, 1, function(r) sprintf("%-16s %s %.0f kb, %s pools in %s crosses, median |dAF| %.2f",
         unname(LONG[[as.character(r[["condition"]])]]), r[["chrom"]],
         as.numeric(r[["pos"]])/1000, r[["pools"]], r[["crosses"]], as.numeric(r[["med_beta"]])))),
  "",
  wrap_para("(E) Every locus, by median effect size. Loci meeting the criteria used for A-D are drawn solid in one of three colours: those containing a published maximum-LOD position for their condition, those containing none, and those in a condition for which no published set exists at all -- Treusch et al. assayed only caffeine and sodium chloride, so glucose and acetate cannot be matched, which is not the same as having no counterpart. Loci that fail those criteria are drawn in pale grey and labelled in the key; they are shown so the filtered set can be seen against everything that was excluded, and they are the reason the largest point in a row may not be one of the loci in A-D."),
  "",
  wrap_para("These are reported as detections, not discoveries. With 1,472 intervals across 179 pools totalling 5.15 times the length of the genome, so the absence of a published locus under a peak is weak evidence that the position is novel. No locus here was tested by allele replacement."),
  "",
  sprintf("Source: run v7.  Built by draft/figures/make_figures.R on %s.", Sys.Date()))

write_csv(all_loci, file.path(OUT, "figure5_loci_table.csv"))
}
say("done — %s", OUT)
