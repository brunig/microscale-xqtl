library(xQTLStats)
library(qs2)
library(parallel)
library(pbmcapply)
library(patchwork)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
# ════════════════════════════════════════════════════════════════════
# ══ PACKAGE OVERRIDES (V6.1) — read this before trusting any LOD ══
# ════════════════════════════════════════════════════════════════════
# xQTLStats is installed as a conda package on Hoffman, inside the env's
# own library, so its R sources cannot be edited in place and a fix means
# rebuilding the package. These six functions are therefore carried here
# instead, exactly as plot_overrides.R used to carry plotSummarySignedShade.
#
# WHY THIS WORKS: R resolves a name from the calling frame outward, and
# the global environment is searched before any attached package. Every
# one of these six is called DIRECTLY from this file, so these definitions
# win. It would NOT work for a function the package calls internally --
# those resolve inside the package namespace and would still get the
# installed version. If a future change needs one of those, it has to go
# through assignInNamespace() or a real reinstall; masking silently would
# not take effect, which is the dangerous failure.
#
# WHAT THEY CHANGE (all four are one bug or its consequence):
#   LODfromZ    NEW. The LOD in the installed package is
#               PvalToLOD(2*pnorm(-|z|)), and p underflows to exactly 0 at
#               |z| = 38.4674 -- an ordinary strong peak, not an extreme
#               one -- so qchisq() returns Inf. 594 of 2411 main fitness
#               peaks in v5. The damage is not the Inf itself but the TIE:
#               peak detection needs lod[i] > lod[i-1], Inf > Inf is FALSE,
#               so the summit pins to the leftmost marker of the plateau
#               and the "3-LOD interval" becomes the width of a
#               floating-point artefact. Evaluated in log space the same
#               function stays finite past |z| = 1e154. The log(4)
#               reproduces the package's existing double-two-sided
#               convention ON PURPOSE, so every currently-finite LOD is
#               preserved against LOD_threshold = 5; fixing that is a
#               separate change that moves every peak call.
#   afdSE       NEW. The second Inf source: afd.se = sqrt(afd*(1-afd)/n)
#               is exactly 0 where the loess allele frequency saturates at
#               0 or 1, so z = beta/0. Clamps afd to 1/(2n) .. 1-1/(2n).
#               ONLY BITES WHEN THE AFDs ARE RECOMPUTED -- a run off a
#               cached AFD still carries the zeros. See the note at
#               `force` in §0.
#   calcAFD     identical to the installed version except afd.se=afdSE(afd,n)
#   calcContrastStats  identical except LOD=LODfromZ(z). `p` is untouched
#               and still underflows; nothing downstream reads it.
#   findQTLPeaks  fixed LOD drop replaced by LOD.drop.frac of the peak's
#               own height, floored and capped (§0). Guards a non-finite
#               peak, which would otherwise make drop_level NaN and abort
#               the walk on NA.
#   plotSummarySignedShade  what plot_overrides.R carried: signed statistic
#               capped at +/-50 so one underflowed p cannot flatten the
#               panel, off-scale runs marked with a dashed line.
#
# SOURCE OF TRUTH: ~/Documents/Git/xQTLStats/R/xQTL_stats.R. These are
# verbatim copies. When that version is what conda installs, DELETE this
# whole section -- two live copies of a statistic, one of which silently
# wins, is how a figure stops matching its own caption.
# -Claude generated

LODfromZ = function(z){ l = qchisq(log(4) + pnorm(-abs(z), log.p=TRUE), df=1,
                                   lower.tail=FALSE, log.p=TRUE)/(2*log(10))
                        l[is.na(l)]=0
                        return(l) }

afdSE = function(afd, n){ a = pmin(pmax(afd, 1/(2*n)), 1-1/(2*n))
                          sqrt((a*(1-a))/n) }

calcAFD=function(phasedCounts, experiment.name='', 
                           sample.size=250000, 
                           sel.strength=.5,
                           bin.width=500,
                           eff.length=600,
                           gmap=NULL,
                           vcf.cross=NULL,
                           uchr = paste0('chr', as.roman(1:16))
                          ){
  
    chrom=data.table::tstrsplit(phasedCounts$ID, '_')[[1]]
    #uchr=paste0('chr', as.roman(1:16))
    chrom=factor(chrom, levels=uchr)
    physical.position=data.table::tstrsplit(phasedCounts$ID, '_', type.convert=T)[[2]]
    coreCols=data.frame(chrom=chrom, 
                        physical.position=physical.position)
    
    if(!is.null(gmap) & !is.null(vcf.cross)){
        genetic.position=stack(jitterGmapVector(getGmapPositions(vcf.cross, gmap, uchr)))$values
        coreCols$genetic.position= genetic.position
    }

    pC = dplyr::bind_cols(coreCols, phasedCounts) %>% 
            dplyr::mutate(nindv=round(sample.size*sel.strength)) %>%
            dplyr::mutate(ndepth=round(sum(p1+p2)/eff.length) ) %>%
            dplyr::mutate(n=min(nindv,ndepth)) %>%
            dplyr::group_by(chrom) %>%
            #dplyr::group_map(~dplyr::mutate(., 
            #   afd=doBlockLoess(.x$p1, .x$p2,.x$physical.position)), .keep=T) %>%
            #dplyr::bind_rows() %>% 
            dplyr::mutate(afd=doBlockLoess(p1, p2,physical.position, bin.width=bin.width)) %>% 
            dplyr::bind_rows() %>%
            dplyr::mutate(afd.se=afdSE(afd, n)) 
       if(!is.null(gmap) & !is.null(vcf.cross)){
            pC=pC %>% dplyr::rename_with(.,~gsub('$', paste0('_', experiment.name), .x), .cols=c(-ID,-chrom, -physical.position,-genetic.position))
       } else{
            pC=pC %>% dplyr::rename_with(.,~gsub('$', paste0('_', experiment.name), .x), .cols=c(-ID,-chrom, -physical.position))

       }
    attr(pC, 'experiment.name')=experiment.name
    return(pC)
   
}

calcContrastStats=function(results, L='meta_high', R='meta_low'){

    combined.results = results %>% purrr::reduce(dplyr::left_join) %>% dplyr::ungroup() #, by=c('ID', 'chrom', 'physical.position') #%>% dplyr::ungroup()
    #contrast.label=paste0(L,'-',R)
    if(!is.null(R)){
        if(is.numeric(R)){
            # R is a constant (e.g. 0.5) with no SE -- test L against constant expectation
            if(grepl('^_', L)){
                combined.results$contrast.beta=as.numeric(unlist(combined.results[,paste0('afd',L)]))-R
                combined.results$contrast.beta.se = as.numeric(unlist(combined.results[,paste0('afd.se',L)]))
            }else{
                combined.results$contrast.beta=as.numeric(unlist(combined.results[,L]))-R
                combined.results$contrast.beta.se = as.numeric(unlist(combined.results[,paste0(L, '.se')]))
            }
        }else if(is.numeric(L)){
            # L is a constant (e.g. 0.5) with no SE -- test R against constant expectation
            if(grepl('^_', R)){
                combined.results$contrast.beta=L-as.numeric(unlist(combined.results[,paste0('afd',R)]))
                combined.results$contrast.beta.se = as.numeric(unlist(combined.results[,paste0('afd.se',R)]))
            }else{
                combined.results$contrast.beta=L-as.numeric(unlist(combined.results[,R]))
                combined.results$contrast.beta.se = as.numeric(unlist(combined.results[,paste0(R, '.se')]))
            }
        }else if(grepl('^_', L)){
            #calc contrasts without meta analysis
            combined.results$contrast.beta=as.numeric(unlist(combined.results[,paste0('afd',L)]-combined.results[,paste0('afd',R)]))
            combined.results$contrast.beta.se = as.numeric(unlist(sqrt(combined.results[,paste0('afd.se',L)]^2+combined.results[,paste0('afd.se',R)]^2)))
        }else{
            combined.results$contrast.beta=as.numeric(unlist(combined.results[,L]-combined.results[,R]))
            combined.results$contrast.beta.se = as.numeric(unlist(sqrt(combined.results[,paste0(L, '.se')]^2+combined.results[,paste0(R,'.se')]^2)))
        }
       combined.results = combined.results %>% dplyr::mutate(z=contrast.beta/contrast.beta.se)
    } else {
        combined.results$z = as.numeric(unlist(combined.results[,L]/combined.results[,paste0(L, '.se')]))
    }
    combined.results=combined.results %>%
        dplyr::mutate(p= 2*pnorm(abs(z), lower.tail=F) ) %>%
        dplyr::mutate(LOD=LODfromZ(z)) %>%suppressWarnings()

    return(combined.results)
}

plotSummarySignedShade <- function(results, simulatedQTL = NULL, effective.n.tests = 600,
                                   p1.name = NULL, p2.name = NULL, use.LOD = FALSE,
                                   cap = 50, mark.offscale = TRUE,
                                   offscale.color = 'grey80') {

    if (use.LOD) {
        raw_y     = sign(results$contrast.beta) * results$LOD
        y_label   = "signed LOD"
        threshold = PvalToLOD(0.05 / effective.n.tests)
    } else {
        raw_y     = sign(results$contrast.beta) * -log10(results$p)
        y_label   = "signed -log10(p)"
        threshold = -log10(0.05 / effective.n.tests)
    }

    if (is.null(cap) || !is.finite(cap)) cap = Inf
    results$signed_y = pmax(pmin(raw_y, cap), -cap)
    if (is.finite(cap)) y_label = sprintf("%s (capped +/-%g)", y_label, cap)

    # An off-scale marker is one the cap actually moved, Inf included: p == 0
    # underflows to Inf, and testing is.finite() would silently drop exactly
    # the peaks most worth flagging. One vline per contiguous run, placed at
    # the run's argmax, so a broad clipped plateau gets one line not hundreds.
    offscale = !is.na(raw_y) & abs(raw_y) > cap
    vl = NULL
    if (mark.offscale && any(offscale)) {
        d   = data.frame(chrom = results$chrom,
                         physical.position = results$physical.position,
                         mag = abs(raw_y), off = offscale,
                         stringsAsFactors = FALSE)
        d   = d[order(d$chrom, d$physical.position), , drop = FALSE]
        n   = nrow(d)
        run = cumsum(c(TRUE, d$chrom[-1] != d$chrom[-n] | d$off[-1] != d$off[-n]))
        vl  = do.call(rbind, lapply(split(d, run), function(r) {
            if (!r$off[1]) return(NULL)
            r[which.max(r$mag), c("chrom", "physical.position"), drop = FALSE]
        }))
    }

    # Parent labels sit at the plotted range, which the cap makes finite.
    y_range = range(results$signed_y, na.rm = TRUE)
    has_parents = !is.null(p1.name) || !is.null(p2.name)
    if (has_parents) {
        sec_breaks = y_range
        sec_labels = c(
            if (!is.null(p2.name)) p2.name else "",
            if (!is.null(p1.name)) p1.name else ""
        )
        y_sec = sec_axis(~., breaks = sec_breaks, labels = sec_labels)
    }

    ggplot(results, aes(x = physical.position, y = signed_y)) +
        annotate("rect",
                 xmin = -Inf, xmax = Inf,
                 ymin = 0, ymax = Inf,
                 fill = "mistyrose", alpha = 0.3) +
        annotate("rect",
                 xmin = -Inf, xmax = Inf,
                 ymin = -Inf, ymax = 0,
                 fill = "lightblue", alpha = 0.3) +
        {
            if (!is.null(vl) && nrow(vl))
                geom_vline(data = vl,
                           aes(xintercept = physical.position),
                           color = offscale.color, linetype = "dashed")
        } +
        geom_line(size = 1, color = "black") +
        {
            if (!is.null(simulatedQTL))
                geom_vline(data = simulatedQTL,
                           aes(xintercept = physical.position),
                           color = "blue")
        } +
        facet_grid(~chrom, scales = "free_x", space = "free_x") +
        geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
        geom_hline(yintercept =  threshold, color = "red") +
        geom_hline(yintercept = -threshold, color = "red") +
        {if(has_parents) scale_y_continuous(name = y_label, sec.axis = y_sec) else scale_y_continuous(name = y_label)}+
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            title = paste(y_label, "with Shaded Backgrounds"),
            x = "Position"
        )
}

findQTLPeaks <- function(results, LOD.threshold = 4, LOD.drop = 2,
                          LOD.drop.frac = 0.1,
                          min.interval = 20000, max.interval = 200000,
                          p1.name = "p1", p2.name = "p2") {

    if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
        !requireNamespace("IRanges", quietly = TRUE)) {
        stop("GenomicRanges and IRanges packages are required. ",
             "Install with: BiocManager::install(c('GenomicRanges', 'IRanges'))")
    }

    if (is.null(min.interval)) min.interval = 0
    if (is.null(max.interval)) max.interval = Inf

    # Confidence interval around one summit.
    #
    # A fixed LOD drop does not survive the dynamic range this statistic
    # produces. Measured over 6,544 fitness summits: below LOD 10 a 3-LOD drop
    # gives a median interval of 159 kb, while above LOD 1,000 it gives 0.4 kb
    # and returns a zero-width interval for 27% of peaks -- a 3-LOD descent off
    # a peak of 10,000 is 0.03% of its height and does not reach the next
    # marker. Descending a fixed fraction of the peak's own height instead
    # holds the median between 51 and 20 kb across the same range. The floor
    # and cap bound what is left: neither a single marker nor half a chromosome
    # is a reportable interval.
    dropInterval = function(pi, lod, pos, n) {
        peak = lod[pi]
        drop = if (is.null(LOD.drop.frac)) LOD.drop else LOD.drop.frac * peak
        drop_level = peak - drop

        # A peak at Inf -- the smoothed allele frequency saturated at exactly 0
        # or 1, so afd.se is 0 and z is +/-Inf -- makes drop and drop_level NaN,
        # and the walk below would then compare against NA and abort. Fall back
        # to the fixed drop there: Inf - LOD.drop is still Inf, so the walk
        # crosses the whole Inf plateau exactly as the fixed-drop code did, and
        # the cap below bounds the result. Clamping afd.se is the actual fix
        # (see afdSE), but it only takes effect when the AFDs are recomputed.
        if (!is.finite(drop_level)) drop_level = peak - LOD.drop

        li = pi
        while (li > 1 && lod[li - 1] >= drop_level) li = li - 1
        ri = pi
        while (ri < n && lod[ri + 1] >= drop_level) ri = ri + 1

        l = pos[li]
        r = pos[ri]
        pp = pos[pi]
        if (min.interval > 0 && (r - l) < min.interval) {
            l = min(l, pp - min.interval/2)
            r = max(r, pp + min.interval/2)
        }
        if (is.finite(max.interval) && (r - l) > max.interval) {
            l = max(l, pp - max.interval/2)
            r = min(r, pp + max.interval/2)
        }
        c(max(l, pos[1]), min(r, pos[n]))
    }

    chroms = unique(results$chrom)
    peak_list = list()

    for (chr in chroms) {
        idx = which(results$chrom == chr)
        pos = results$physical.position[idx]
        lod = results$LOD[idx]
        beta = results$contrast.beta[idx]
        n = length(lod)
        if (n < 3) next

        # Order by position within chromosome
        ord = order(pos)
        pos = pos[ord]
        lod = lod[ord]
        beta = beta[ord]

        # Replace NA/NaN LOD with 0 so comparisons work cleanly
        lod[is.na(lod)] = 0

        # Find local maxima above threshold
        # A local maximum: LOD[i] > LOD[i-1] and LOD[i] >= LOD[i+1]
        is_peak = rep(FALSE, n)
        for (i in 2:(n - 1)) {
            if (lod[i] >= LOD.threshold &&
                lod[i] > lod[i - 1] &&
                lod[i] >= lod[i + 1]) {
                is_peak[i] = TRUE
            }
        }
        # Check endpoints
        if (n >= 2 && lod[1] >= LOD.threshold && lod[1] >= lod[2]) is_peak[1] = TRUE
        if (n >= 2 && lod[n] >= LOD.threshold && lod[n] > lod[n - 1]) is_peak[n] = TRUE

        peak_idx = which(is_peak)
        if (length(peak_idx) == 0) next

        # Merge peaks that share overlapping LOD drop intervals
        # First pass: compute raw intervals
        raw_intervals = data.frame(
            peak_i = peak_idx,
            peak_pos = pos[peak_idx],
            peak_lod = lod[peak_idx],
            peak_beta = beta[peak_idx],
            left = NA_real_,
            right = NA_real_
        )

        for (j in seq_len(nrow(raw_intervals))) {
            iv = dropInterval(raw_intervals$peak_i[j], lod, pos, n)
            raw_intervals$left[j]  = iv[1]
            raw_intervals$right[j] = iv[2]
        }

        # Merge overlapping intervals, keeping the highest peak
        raw_intervals = raw_intervals[order(raw_intervals$left), ]
        merged = list()
        current = raw_intervals[1, ]
        if (nrow(raw_intervals) > 1) {
            for (j in 2:nrow(raw_intervals)) {
                if (raw_intervals$left[j] <= current$right) {
                    # Overlapping — extend right boundary, keep higher peak
                    current$right = max(current$right, raw_intervals$right[j])
                    if (raw_intervals$peak_lod[j] > current$peak_lod) {
                        current$peak_i = raw_intervals$peak_i[j]
                        current$peak_pos = raw_intervals$peak_pos[j]
                        current$peak_lod = raw_intervals$peak_lod[j]
                        current$peak_beta = raw_intervals$peak_beta[j]
                    }
                } else {
                    merged[[length(merged) + 1]] = current
                    current = raw_intervals[j, ]
                }
            }
        }
        merged[[length(merged) + 1]] = current
        merged = do.call(rbind, merged)

        # Recompute interval boundaries from the winning peak after merge
        for (j in seq_len(nrow(merged))) {
            iv = dropInterval(merged$peak_i[j], lod, pos, n)
            merged$left[j]  = iv[1]
            merged$right[j] = iv[2]
        }

        merged$chrom = chr
        peak_list[[chr]] = merged
    }

    if (length(peak_list) == 0) {
        gr = GenomicRanges::GRanges()
        GenomicRanges::mcols(gr) = data.frame(
            peak.pos = integer(0),
            peak.LOD = numeric(0),
            peak.beta = numeric(0),
            direction = integer(0),
            high.parent = character(0),
            low.parent = character(0),
            stringsAsFactors = FALSE
        )
        return(gr)
    }

    all_peaks = do.call(rbind, peak_list)
    direction = sign(all_peaks$peak_beta)
    high.parent = ifelse(direction > 0, p1.name, p2.name)
    low.parent = ifelse(direction > 0, p2.name, p1.name)

    gr = GenomicRanges::GRanges(
        seqnames = all_peaks$chrom,
        ranges = IRanges::IRanges(start = all_peaks$left, end = all_peaks$right),
        peak.pos = all_peaks$peak_pos,
        peak.LOD = all_peaks$peak_lod,
        peak.beta = all_peaks$peak_beta,
        direction = direction,
        high.parent = high.parent,
        low.parent = low.parent
    )

    gr = sort(gr)
    return(gr)
}

cat('package overrides active (V6.1): LODfromZ, afdSE, calcAFD,',
    'calcContrastStats,\n  findQTLPeaks, plotSummarySignedShade',
    '— masking xQTLStats', as.character(utils::packageVersion('xQTLStats')), '\n')
# ════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════
# xQTLpaper — stress-vs-glucose xQTL analysis (V7)
# ══════════════════════════════════════════════════════════════════════
# V7 is V6.1 with five changes. None of them touches §1-15, so a V7 run
# at v6b's thresholds reproduces v6b's tables; that equality is the
# acceptance test and should be checked before v6b stops being the run
# of record. Verified on the Mac 2026-09-21: every landscape CSV byte-
# identical except last-digit float noise in p_pois / q_bh / p_binom.
#
#   1. The `which_cond` protoF fix, ported in from the standalone
#      qtl_landscape.R. It was made there on 2026-09-18, one day AFTER
#      V6.1 inlined that file, so it had never reached the driver: the
#      published overlay was being drawn on glucose and acetate panels,
#      which Treusch et al. never assayed. This is why v6b's
#      protoF_surface_*.png cannot be used as rendered.
#
#   2. §18, GO-term enrichment (qtl_go_run). The one downstream child
#      V6.1 did not absorb. Reads the peak tables and the GFF only.
#
#   3. The published-concordance tests run on BOTH peak sets and BOTH
#      units of detection, written as *_bysets.csv beside the originals.
#      Treusch et al. mapped with MULTIPOOL in `-m contrast` mode, so
#      the comparable set here is the stress contrast, not the fitness
#      peaks against 0.5. Nothing above that block changed.
#
#   4. Two stale things fixed: protoD's hardcoded "LOD is Inf for
#      816/3308 peaks" subtitle (both numbers from a superseded run),
#      and the uniform-rank expectations for the candidate-gene
#      self-check, which the draft quoted but no code emitted.
#
#   5. run_inline() now handles the message condition directly. Tested
#      2026-09-21: the sink DOES capture message() under Rscript, so the
#      old note claiming run_inline "swallows" it was a misdiagnosis --
#      but v6b's logs really did lose those lines, so the run happened
#      somewhere the sink did not apply (IRkernel is the likely path).
#      Handling the condition makes the log right either way.
#
# ---- V6's changes, unchanged -----------------------------------------
# V6 is V5 with two changes and no change to the analysis itself. The
# thresholds, the flag-and-fork architecture and every number this
# script computes are V5's; a V6 run of an unchanged tag reproduces V5
# exactly. What changed is WHERE the outputs go and WHAT runs after.
#
#   1. Run umbrella (§0). V5 wrote every run into the same contrasts/,
#      cache/ and qc/ directories. Nothing on disk recorded which run
#      made it, so a half-finished re-run left fit_peaks.csv from today
#      next to qtl_variants_fitness.csv from last week and the two were
#      indistinguishable. V6 writes everything under runs/<run_tag>/ and
#      writes nothing outside it. A table's directory is now its
#      provenance.
#
#   2. Downstream scripts (§16-17). qtl_landscape.R and qtl_candidates.R
#      previously found their input by SEARCHING the filesystem two
#      directory levels deep for fit_peaks.csv — which is exactly how a
#      figure ends up describing a peak table nobody can identify. They
#      now run at the end of this script, in their own R processes, with
#      QTL_PEAK_DIR pinned to this run's contrasts/.
#
# One consequence worth stating: the cache lives inside the run tree, so
# a NEW run_tag recomputes counts and AFDs from scratch. That is the
# point — a tag means "these tables were all made by one invocation".
# To warm-start anyway, see the hardlink note at the end of §0.
#
# ---- V5's architecture, unchanged ------------------------------------
# V5 is V4 with one architectural change: NOTHING IS REMOVED.
#
# V4 gated. Each QC step subset `samples`, so a sample that failed a
# threshold vanished and left only a line in an attrition ledger. That
# is wrong for a published pipeline for two reasons. The samples nearest
# a threshold are exactly the ones whose classification is least
# certain, and a gate deletes precisely those. And a reader cannot
# audit a cut whose casualties were never computed — "what did the
# excluded ones look like?" has no answer.
#
# V5 runs everything and LABELS it. Every sample that can produce a
# result does produce one; each check writes a flag instead of removing
# a row; and one wide results key (qc/results_key.csv, one row per
# keyfile row, keyfile columns carried verbatim so it merges straight
# back) records where each sample landed and why. Plots for flagged
# samples are written to figures/flagged/ mirroring the main tree, so
# they can be looked at rather than imagined.
#
#   databin == 'main'     ran, no flags. The analysis dataset.
#   databin == 'flagged'  ran, >=1 flag. Peaks ARE called and tagged,
#                         plots go to figures/flagged/. Held separate
#                         from main; see flag_reasons for why.
#   databin == 'no_data'  no analysable result exists (no parseable
#                         count table, parent absent from Z, or AFD did
#                         not converge). Row still present in the key.
#
# The distinction that matters: 'flagged' is a judgement — the sample
# ran and you may disagree with the call. 'no_data' is structural —
# there is no number to disagree with. Do not collapse them.
#
# Peaks are called ONCE over everything and carry a `databin` column;
# the fork is a filter at the point of use, not a second pipeline.
# Downstream consumers (reproducibility, variants, CNV, enrichment)
# restrict to main; the peak CSVs contain both, labelled.
#
# Threshold changes in V5 (see §0): `min_depth` no longer removes
# anything, it only raises flag_low_depth; a sample too thin to fit
# falls out on its own as afd_failed, which is the honest reason.
#
# Pipeline stages, each cached to disk so a re-run resumes rather than
# recomputing (see `stages` / `force` / `cache_prefix` in §0):
#    counts   phase allele counts per sample against its parents
#    afd      smooth allele frequency per chromosome (loess)
#    contrast fitness (vs 0.5) and stress-vs-control contrasts
#    peaks    LOD peaks + intervals, tagged with databin
#    figures  per-pair / per-sample panels and summary grids
#    repro    cross-level reproducibility heatmaps and dot-plots (main)
#    variants segregating variants under each interval (main)
#    cnv      parental CNVs under each interval (main)
#    enrich   gene-set enrichment in intervals (main, optional)
#
# Original single-cross workflow and the xQTLStats calls are Josh
# Bloom's; the batch driver, flags and caching are -Claude generated.
# ══════════════════════════════════════════════════════════════════════


# ── 0. Configuration ──────────────────────────────────────────────────
# Where the QC/analysis boundary sits: pre-flight QC (BAM metrics, raw
# coverage, plate diagnostics, parent identity) is parent-agnostic, runs
# once per sequencing run, and in V4 *gated* which samples entered here.
# In V5 it no longer gates — qc_pipeline.R's verdict arrives as a flag
# column like any other. It still has to run first, because identity
# errors change parent_a and therefore change §3 and §4.
# -Claude generated
env = 'hoffman2'   # 'hoffman2' | 'local'

paths_by_env = list(
    hoffman2 = list(
        Z.qs.file    = '/u/project/kruglyak/bruni/projects/xQTLpaper/files/Z_final.qs',
        vlookup      = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/variant_lookup.qs',
        count_dir    = '/u/project/kruglyak/bruni/projects/xQTLpaper/v8/count_files/',   # v8
        key.file     = '/u/project/kruglyak/bruni/projects/xQTLpaper/files_v8/flatkey_A1_v8b.csv',   # v8b
        qc_dir       = '/u/project/kruglyak/bruni/projects/xQTLpaper/v8/qc_v8b',         # v8b: v8 QC + corrected unresolved list
        script_dir   = '/u/project/kruglyak/bruni/projects/xQTLpaper/files_v8',   # v8
        gff_file     = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz',
        cnv.file     = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/full3039Matrix.CNVnatorResults.tsv',
        fungamr_file = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/fungAMR/41564_2025_2084_MOESM3_ESM.xlsx',
        # V7 §18. Neither file is on Hoffman yet -- scp both beside the driver:
        #   scp ~/Desktop/yeast/S288C_reference_genome_R64-5-1_20240529/gene_association_R64-5-1_20240529.sgd \
        #       ~/Desktop/yeast/S288C_reference_genome_R64-5-1_20240529/go-basic.obo \
        #       hoffman:/u/project/kruglyak/bruni/projects/xQTLpaper/files/
        # Flat files rather than org.Sc.sgd.db + GO.db deliberately: those two
        # are installed on neither machine, and the GAF is R64-5-1, the same
        # build as the GFF above, so annotation and gene models agree.
        gaf_file     = '/u/project/kruglyak/bruni/projects/xQTLpaper/files/gene_association_R64-5-1_20240529.sgd',
        obo_file     = '/u/project/kruglyak/bruni/projects/xQTLpaper/files/go-basic.obo'
    ),
    local = list(
        Z.qs.file    = '/media/hoffman2/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/Z_1011_filtered.qs',
        vlookup      = '/data1/yeast/reference/3000/VCFs_SNPs/variant_lookup.qs',
        count_dir    = NA_character_,
        key.file     = NA_character_,
        qc_dir       = NA_character_,
        script_dir   = NA_character_,
        gff_file     = '/data1/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz',
        cnv.file     = '/media/hoffman2/PUBLIC_SHARED/yeast/reference/3000/full3039Matrix.CNVnatorResults.tsv',
        fungamr_file = '/data1/yeast/reference/fungAMR/41564_2025_2084_MOESM3_ESM.xlsx',
        gaf_file     = '~/Desktop/yeast/S288C_reference_genome_R64-5-1_20240529/gene_association_R64-5-1_20240529.sgd',
        obo_file     = '~/Desktop/yeast/S288C_reference_genome_R64-5-1_20240529/go-basic.obo'
    )
)
stopifnot(env %in% names(paths_by_env))
for (.k in names(paths_by_env[[env]]))
    assign(.k, paths_by_env[[env]][[.k]])

# ---- Run umbrella ----------------------------------------------------
# run_tag is the only thing to change between runs. Reusing a tag
# RESUMES that run from its own cache; a new tag is a clean run. For one
# tree per invocation:  run_tag = format(Sys.time(), '%Y%m%d_%H%M')
# -Claude generated
# V7: 'v7', NOT 'v6'. This file was copied from the working V6 driver,
# which was sitting in comparator configuration (run_tag 'v6',
# min_med_depth 5). Running it unchanged would have written into
# runs/v6 and destroyed the min_med_depth = 5 comparator that the
# sensitivity analysis depends on.
run_tag  = 'v8b'  # v8b: v8 + parent fixes A11->AAM, G10->ADD, A10 unresolved (tools/parent_genotype_check.R)

proj_dir = dirname(count_dir)                  # NA_character_ in env 'local'
run_dir  = file.path(proj_dir, 'runs', run_tag)

# qc_pipeline.R's outputs are INPUTS to this script (sample_allowlist.csv,
# unresolved_wells_A1.txt); they are parent-agnostic, run once per
# sequencing run and are shared by every run_tag. What THIS script writes
# about QC (results_key.csv, reproducibility_main.csv) is a per-run
# output. V5 pointed both at one `qc_dir`, so re-running with different
# thresholds rewrote the key that the previous run's tables were built
# against — and left no trace that it had. The two are now separate.
qc_in_dir  = qc_dir                            # read-only, shared
qc_out_dir = file.path(run_dir, 'qc')          # written by this run

out_dir     = file.path(run_dir, 'contrasts')
cache_dir   = file.path(run_dir, 'cache')
fig_dir     = file.path(out_dir, 'figures')
flag_dir    = file.path(fig_dir, 'flagged')
land_dir    = file.path(run_dir, 'landscape')   # §16, qtl_landscape.R
cand_dir    = file.path(run_dir, 'candidates')  # §17, qtl_candidates.R
go_dir      = file.path(run_dir, 'go')          # §18, qtl_go_run output (V7)
script_snap = file.path(run_dir, 'scripts')     # copies of what ran

# V6.1: the downstream analysis is inlined below (§16-17), so the code
# this run executed is this file. script_dir is otherwise unused now.
# V7: this MUST name the file you are actually running. It was still
# 'xQTL_analysis6.R' in the first V7 draft, which would have archived V6
# into runs/v7/scripts/ as "the code that produced this run" -- the exact
# provenance failure the run umbrella exists to prevent.
driver_file = if (is.na(script_dir)) NA_character_ else
              file.path(script_dir, 'xQTL_analysis8b.R')

for (.d in Filter(function(x) !is.na(x),
                  c(run_dir, out_dir, cache_dir, fig_dir, flag_dir,
                    qc_out_dir, land_dir, cand_dir, go_dir, script_snap)))
    dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# Recomputing counts and AFDs for a new tag costs hours. To start a new
# tag from a finished run's cache without copying the data, hardlink it
# first — outputs still land in the new tree:
#     cp -al <proj>/runs/v6/cache <proj>/runs/v6_try2/cache
# Do this only when the inputs and thresholds that produced those caches
# are the ones you are about to run with. If in doubt, don't.

# ---- Design ----------------------------------------------------------
# 96 wells x 4 plates. A well identifies the source pool and carries the
# same parent across all four plates, so the well is the matching key
# for a stress-vs-control pair AND the replicate unit for a cross.
condition_col = 'Drug3Let'
control_cond  = 'Glu'                       # HCDM, no stressor
match_keys    = c('parent_a', 'parent_alpha', 'Wells')

# ---- Stages ----------------------------------------------------------
# Data stages (counts, afd, contrast, peaks) are always cache-backed.
# Output stages below are optional; drop a name to skip that work.
# 'landscape' and 'candidates' are V6's §16-17; they need 'variants'
# to have run (qtl_candidates.R consumes qtl_variants_<tag>.csv).
# 'go' is V7's §18 and needs only the peak tables -- no variants, no cache.
stages = c('figures', 'repro', 'variants', 'cnv', 'enrich',
           'landscape', 'candidates', 'go')
force  = character(0)   # e.g. c('afds', 'fit_results', 'drug_results')

# For the V6.1 rerun off the v6 cache:
#     force = c('fit_results', 'drug_results')
# recomputes the contrasts (new LOD) and the peaks (new interval rule)
# while reusing the cached counts and AFDs. That fixes the pnorm-underflow
# Inf, which is the one that affects 594 of 2411 peaks.
#
# It does NOT fix the second Inf source: afd.se is already 0 in the cached
# AFDs wherever the loess saturated, and afdSE() only clamps when calcAFD
# actually runs. Those summits stay at LOD = Inf and findQTLPeaks falls
# back to the fixed drop for them. Add 'afds' to force as well to clear
# them -- that re-smooths every sample but reuses the phased counts, which
# is the expensive half. It affected 2 peak summits and ~100 markers in v5,
# so it is a judgement call, not a blocker.

# The prefix is vestigial in V6 — each run_tag has its own cache
# directory, so there is nothing to collide with. It is kept so the V5
# regression check below can still read a V4/V5 cache by path.
#
# V7: LEAVE THIS AT 'v6_'. It looks wrong next to run_tag = 'v7' and it is
# not: the prefix is NOT the run tag, it is the cache filename convention,
# and it is what makes a cache interchangeable between tags. A hardlinked
# cache carries files named v6_countdfs.qs / v6_afds.qs / v6_drug_results.qs;
# setting this to 'v7_' would make V7 look for v7_*.qs, miss every one of
# them, and silently recompute the whole pipeline from scratch — which on
# 2026-09-21 is exactly what killed 142 forked workers in §8. The only cost
# of leaving it is that runs/v7/cache/ holds files whose names start v6_.
# Clear a whole run with one line:   rm -r <proj>/runs/<run_tag>
cache_prefix = 'v8b_'

# ---- Thresholds ------------------------------------------------------
# NOTE the changed meaning of min_depth. In V4 it REMOVED samples at or
# below it. Here it only raises flag_low_depth: the sample still gets an
# AFD, a contrast, peaks and a plot in figures/flagged/. Samples too
# thin to fit at all fall out as afd_failed, which is a fact rather than
# a threshold. Raising this number now moves plots between folders; it
# no longer deletes anything.
min_depth      = 1      # mean per-site depth over segregating sites
LOD_threshold  = 5
LOD_drop       = 3      # fixed drop; used only if LOD_drop_frac is NULL

# Interval rule. The fixed 3-LOD drop returns a ZERO-WIDTH interval for
# 27% of peaks above LOD 1000, because the LOD scale here is inflated
# (nominal binomial SE, no overdispersion term) and a 3-LOD drop off a
# LOD-10000 summit is nothing. A drop proportional to the peak is
# scale-free; the floor and cap keep it interpretable at both ends.
# Evidence: ~/Desktop/AFxQTL/misc/interval_test/. Set LOD_drop_frac=NULL,
# min_interval=0, max_interval=Inf to get the old fixed-drop behaviour.
LOD_drop_frac  = 0.1
min_interval   = 20000
max_interval   = 200000
afd_sample_size = 1e4   # raise to 5e4 to recover marginal-depth samples
afd_sel_strength = 0.9

# Genome-wide AF baseline. A sample whose smoothed AF sits away from 0.5
# is either mis-parented or (working hypothesis, 2026-08-28) carries a
# TRIPLOID parent, where 2:1 dosage puts the expectation at 1/3 or 2/3.
# Both break the 1:1 segregation calcAFD assumes, so the affected
# samples are flagged and forked rather than having their null rescaled
# — rescaling would make a wrong parent look clean.
#
# Leave NA to flag nothing. Set it from evidence: run detect_skew.R,
# read the sorted per-well table, and put the cut in the GAP between the
# clean bulk and the 1/3 cluster. Do not pick a round number.
# Set 2026-08-28 from 91 well-level medians. The distribution is
# unimodal at ~0.49 with a continuous left tail and NO mode at 1/3, so
# the triploid hypothesis is not supported and there is no natural gap
# to cut in. This number is therefore empirical, not mechanistic: from
# the cross-level model in qc/skew_vs_depth.csv the skew buys ~281
# peaks per unit, so 0.05 is the point past which the artifact
# contributes ~14 spurious peaks. Flags 12 wells; AAL sits at 0.042
# to 0.045 and stays in main.
skew_cut = 0.05

# Loess undershoot at chromosome ends puts healthy samples slightly
# outside [0,1], so divergence is a FRACTION, not a min/max test —
# af_min < 0 would mislabel ~60 healthy samples.
diverged_frac = 0.01

# Analysis-quality thresholds, applied on the AFD's own median depth
# rather than the count table's mean. These reproduce the ladder
# assess_fit.R used, so that databin 'main' means what flag == 'ok'
# meant there; without them main is ~60 samples more permissive.
min_med_depth = 2      # median (p1 + p2) over sites
                       # V7: 2, the run-of-record setting. The V6 working
                       # copy this was cloned from had 5 (the comparator).
                       # The V7 acceptance test -- v7/landscape must equal
                       # v6b/landscape -- is meaningless at any other value.
                       # Set to 5 and retag to rebuild the comparator.
max_excess    = 0.12   # mean|AF - 0.5| above the binomial floor

# mc.cores. Per the 2026-05-06 debugging, 8 needs at least -m 16 in
# h2jupynb; forks duplicate cross.Z slices and get OOM-reaped below
# that. Drop to 4 on -m 8, and to 1 when bisecting a worker crash.
ncores = 8

# The count-table step (§4) is the heaviest: every fork touches a
# different cross.Z slice, so copy-on-write does not help and RSS scales
# with cores. A killed worker is NOT an R error — it surfaces as
# mclapply's "N function calls resulted in an error". Two causes, told
# apart by how the count responds to this number:
#   scales with cores  -> memory. Lower it, or ask for more with -m.
#   barely moves       -> the pbmclapply SIGALRM/ZMQ collision described
#                         at par_lapply() in §1; set to 1 to confirm.
# V5 builds count tables for every sample, not just the gated survivors,
# so this step does ~11% more work than V4.
count_cores = 1

# ---- Downstream scripts (§16-17) -------------------------------------
# qtl_candidates.R ranks genes under loci. chrVII 440-600 kb is the BBQ
# growth QTL: a finding, not an artefact, so KEEPING it is defensible —
# but it is detected in most pools and therefore sits at the top of any
# recurrence-weighted list, and the two settings give materially
# different allele-series numbers. V5 left this decided inside that
# script, where the run had no record of it. It is now a driver setting
# and goes into the run manifest. Set it deliberately.
keep_chr7 = TRUE

# Which variant dumps to rank genes in. Both need `variants` to have run.
cand_tags = c('fitness', 'stress')

# ---- GO-term enrichment (V7 §18) -------------------------------------
# The first thing in this project that emits GO terms. Needs only the peak
# tables and the GFF, so it can run anywhere the run tree can be read.
#
# go_source '' picks orgdb when org.Sc.sgd.db is installed and gaf
# otherwise. Set it explicitly to pin which annotation a run used --
# 'gaf' is the decided path (see gaf_file / obo_file in §0), because the
# Bioconductor packages are on neither machine and the GAF matches the
# R64-5-1 build of gff_file.
#
# go_peak_sets: BOTH, deliberately. Figures 4D and 5D ask condition-
# specific questions, and only the matched contrast supports a per-locus
# condition-specific claim; fitness loci carry everything that also acts
# in glucose. Quote the contrast rows for those panels.
#
# Runtime is go_perms x conditions x peak sets locus-shift passes. 1000 x
# 4 x 2 is the expensive part of this section; drop go_perms to 50 for a
# smoke test and read the console summary before trusting a full run.
go_source    = 'gaf'          # 'gaf' | 'orgdb' | '' to auto-detect
go_ont       = 'BP'           # BP | MF | CC
go_perms     = 1000
go_min_term  = 5              # ignore tiny terms
go_max_term  = 500            # and unspecific ones
go_peak_sets = c('fitness', 'contrast')

# ---- Gene-set enrichment (optional) ----------------------------------
# Named list of character vectors of systematic gene names. Left NULL
# the enrich stage is skipped. -Claude generated
gene_sets = NULL
n_perm    = 10000


# ── 1. Helpers ────────────────────────────────────────────────────────
# One failure predicate. Hard worker crashes come back as try-error
# objects (atomic), soft failures as a list with .failed.
ok = function(x) {
    if (inherits(x, 'try-error')) return(FALSE)
    if (is.list(x) && '.failed' %in% names(x)) return(!isTRUE(x$.failed))
    TRUE
}

# Lazy disk cache. `expr` is a promise, so a cache hit never evaluates
# the body. `force` names are the UNPREFIXED stage names.
cache = function(name, expr) {
    f = file.path(cache_dir, paste0(cache_prefix, name, '.qs'))
    if (!(name %in% force) && file.exists(f)) {
        cat('  [cache] ', name, '\n', sep = '')
        return(qs2::qs_read(f))
    }
    v = expr
    qs2::qs_save(v, f)
    cat('  [saved] ', name, '\n', sep = '')
    v
}

running = function(stage) stage %in% stages

# Parallel map that is safe inside IRkernel. pbmclapply drives its
# progress bar from a SIGALRM timer; forked children inherit the
# notebook's ZMQ socket, the timer interrupts their writes (EINTR), and
# the children die — which looks exactly like an OOM without being one.
# The tell is that the failure count barely moves with mc.cores.
in_jupyter = isTRUE('IRkernel' %in% loadedNamespaces())
#
# V7: KILLED WORKERS ARE RETRIED SERIALLY. Two kinds of failure reach this
# function and they are not the same thing:
#
#   a caught error   every call site wraps its body in tryCatch and returns
#                    list(.failed = TRUE, msg = ...). That is a DATA failure --
#                    this sample has no parseable count file, this AFD did not
#                    converge -- and it is a legitimate result that ok() routes
#                    to the no_data bin. Keep it.
#
#   a killed worker  the child died by signal, so the tryCatch inside it never
#                    ran and mclapply hands back a try-error. That is an
#                    ENVIRONMENT failure -- the OOM killer, or a fork that
#                    inherited a notebook's ZMQ socket -- and it says nothing
#                    about the data.
#
# Before V7 the second was indistinguishable from the first downstream: ok()
# filtered both, the run carried on, and §9 called peaks over the survivors.
# On 2026-09-21 a from-scratch run lost 142 of the stress pairs that way and
# would have written a stress_peaks.csv missing them, with one line in a log
# as the only trace. A published table must never be quietly partial.
#
# Retrying in the parent is the fix rather than a workaround: there is one
# process, so whatever memory the forks could not have is available, and the
# answer is identical because FUN is deterministic. It costs wall-clock on a
# run that was going to be wrong otherwise. It is also agnostic to WHICH of
# the two environment causes bit, which matters because one run cannot tell
# them apart. If everything dies, the stop() below is the honest outcome --
# better than a table nobody can reproduce.
# -Claude generated
par_lapply = function(X, FUN, cores) {
    if (cores <= 1) return(lapply(X, FUN))
    out = if (in_jupyter)
        parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
    else
        pbmclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)

    # A killed child comes back EITHER as a try-error OR as a plain NULL,
    # depending on how far it got before the signal -- mclapply warns
    # "did not deliver a result" for the NULL case. Testing only for
    # try-error misses it, and a NULL then flows downstream as a legitimate
    # empty result. Verified 2026-09-21 against a worker killed mid-call.
    .lost = function(x) is.null(x) || inherits(x, 'try-error')
    dead = which(vapply(out, .lost, logical(1)))
    if (length(dead)) {
        cat('!! ', length(dead), ' of ', length(X), ' workers were killed ',
            '(not an R error — the child died, so its tryCatch never ran).\n',
            '   Cause is memory or the notebook fork, never the data.\n',
            '   Retrying those ', length(dead), ' serially in the parent.\n', sep = '')
        for (i in dead)
            out[[i]] = tryCatch(FUN(X[[i]]),
                                error = function(e) list(.failed = TRUE,
                                                         msg = conditionMessage(e)))
        still = which(vapply(out, .lost, logical(1)))
        if (length(still))
            stop(sprintf(paste0('%d of %d elements could not be computed even ',
                 'serially. Refusing to continue: the tables this run would ',
                 'write would be silently incomplete. Lower ncores, ask for ',
                 'more memory, or run outside the notebook.'),
                 length(still), length(X)), call. = FALSE)
        cat('   recovered all ', length(dead), ' on retry.\n', sep = '')
    }
    out
}
if (in_jupyter)
    cat('running under IRkernel — using mclapply without a progress bar\n')

.short = function(s) sub('^Sample', 'S', s)

# ---- The results key -------------------------------------------------
# One row per keyfile row, keyfile columns carried verbatim so it merges
# straight back onto the flat key. Every check appends a logical
# flag_<name> column; nothing is ever removed. -Claude generated
key = NULL
flag_add = function(ids, name) {
    col = paste0('flag_', name)
    k = key
    if (!col %in% names(k)) k[[col]] = FALSE
    hit = k$Sample_ID %in% ids
    k[[col]][hit] = TRUE
    key <<- k
    cat('  flag ', name, ': ', sum(hit), ' samples\n', sep = '')
    invisible(sum(hit))
}
key_set = function(ids, col, values) {
    k = key
    if (!col %in% names(k)) k[[col]] = NA
    k[[col]][match(ids, k$Sample_ID)] = values
    key <<- k
    invisible(NULL)
}

# Figure path resolver. Every ggsave/png in this script goes through
# here, so the main/flagged split is one function rather than a
# condition at 20 call sites. -Claude generated
fig_path = function(bin, ...) {
    p = file.path(if (identical(bin, 'main')) fig_dir else flag_dir, ...)
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    p
}


# ── 2. Load reference and keyfile ─────────────────────────────────────
gmaps = readRDS(system.file('reference', 'yeast_gmaps.RDS',
                            package = 'xQTLStats'))
gmap  = gmaps[['A']]
uchr  = paste0('chr', as.character(as.roman(1:16)))

samples = as.data.frame(readr::read_csv(key.file, show_col_types = FALSE))
stopifnot(all(c('Sample_ID', 'Wells', 'parent_a', 'parent_alpha',
                condition_col) %in% names(samples)))
samples$parent_a     = trimws(samples$parent_a)
samples$parent_alpha = trimws(samples$parent_alpha)
samples$cross        = paste(samples$parent_a, samples$parent_alpha, sep = 'x')
samples$condition    = samples[[condition_col]]

# The key starts as the whole keyfile, verbatim. It only ever gains
# columns.
key = samples

cat('keyfile: ', nrow(samples), ' samples, ',
    length(unique(samples$Wells)), ' wells, ',
    length(unique(samples$cross)), ' crosses, conditions: ',
    paste(sort(unique(samples$condition)), collapse = ','), '\n', sep = '')

# ---- Flag: pre-flight QC allowlist -----------------------------------
# qc_pipeline.R writes qc/sample_allowlist.csv (Sample_ID + qc_pass).
# In V4 this removed samples. Here it flags them: they run, and their
# plots land in figures/flagged/ so the call can be checked.
allow.file = file.path(qc_in_dir, 'sample_allowlist.csv')
if (file.exists(allow.file)) {
    allow  = readr::read_csv(allow.file, show_col_types = FALSE)
    passed = allow$Sample_ID[allow$qc_pass %in% TRUE]
    flag_add(setdiff(samples$Sample_ID, passed), 'qc_allowlist')
} else {
    cat('!! no ', allow.file, ' — pre-flight QC not applied.\n', sep = '')
    key$flag_qc_allowlist = FALSE
}

# ---- Flag: wells with no defensible parent ---------------------------
# These phase against the best available parent, which is not the right
# one. Their peaks are misphased BY CONSTRUCTION — positioned by the
# wrong variant set, not merely noisy. They are computed and plotted so
# the failure is visible, and flag_reasons must be loud enough that
# nobody merges them into a main-databin table by accident.
unres.file = file.path(qc_in_dir, 'unresolved_wells_A1.txt')
if (file.exists(unres.file)) {
    unres = trimws(readLines(unres.file, warn = FALSE)); unres = unres[nzchar(unres)]
    flag_add(samples$Sample_ID[samples$Wells %in% unres], 'parent_unresolved')
} else {
    key$flag_parent_unresolved = FALSE
}

# ---- Structural: parents present in Z --------------------------------
# Not a judgement call. Without a Z row there is no cross.Z slice, so
# makeCountTablesGATK2 has nothing to phase against and no result can
# exist. These become databin 'no_data'.
Z.obj = qs2::qs_read(Z.qs.file)
zrows = rownames(Z.obj$Z)
notZ  = !(samples$parent_a %in% zrows) | !(samples$parent_alpha %in% zrows)
if (any(notZ)) {
    cat('parents absent from Z: ',
        paste(unique(c(samples$parent_a[notZ], samples$parent_alpha[notZ])),
              collapse = ','), '\n', sep = '')
    flag_add(samples$Sample_ID[notZ], 'parent_not_in_Z')
} else {
    key$flag_parent_not_in_Z = FALSE
}
runnable = samples[!notZ, , drop = FALSE]
stopifnot(nrow(runnable) > 0)


# ── 3. Per-cross Z slices ─────────────────────────────────────────────
# The full Z is too large to fork safely, so reduce it to the two parent
# rows per cross before any mclapply and drop the original. -JBloom
unique.crosses = unique(paste(runnable$parent_a, runnable$parent_alpha,
                              sep = '|'))
cross.Z = lapply(unique.crosses, function(cr) {
    p.names = strsplit(cr, '\\|')[[1]]
    list(Z = Z.obj$Z[p.names, ], variant_ids = Z.obj$variant_ids)
})
names(cross.Z) = unique.crosses
rm(Z.obj); invisible(gc())
cat('cross.Z slices:', length(cross.Z), '\n')


# ── 4. Phased count tables ────────────────────────────────────────────
# Reads each sample's ASEReadCounter output, finds sites segregating
# between its two parents, and phases counts so p1 = parent_a and
# p2 = parent_alpha. -JBloom
#
# Pre-fork memory check. Read the cross.Z figure with suspicion:
# object.size() counts a shared character vector once per reference, and
# all cross.Z entries point at the SAME variant_ids, so it reports ~12
# GB where resident size is a few hundred MB. Useful for spotting a
# stray Z left in the session; not a memory budget. -Claude generated
.big = sort(vapply(ls(envir = globalenv()),
    function(nm) as.numeric(object.size(get(nm, envir = globalenv()))), 0),
    decreasing = TRUE)
.big = .big[.big > 50e6]
if (length(.big)) {
    cat('\nlarge objects in the session before forking:\n')
    for (nm in names(.big))
        cat(sprintf('  %-20s %6.0f MB\n', nm, .big[[nm]] / 1e6))
    cat(sprintf('  total %.1f GB x %d workers.\n',
                sum(.big) / 1e9, count_cores))
    if (sum(.big) > 1.5e9)
        cat('  !! consider a fresh session — this step OOMs easily.\n')
}

countdfs = cache('countdfs', {
    cd = par_lapply(seq_len(nrow(runnable)), function(n) {
        p.names = c(runnable$parent_a[n], runnable$parent_alpha[n])
        tryCatch(
            makeCountTablesGATK2(count_dir,
                paste0(runnable$Sample_ID[n], '.txt'),
                p.names, cross.Z[[paste(p.names, collapse = '|')]], gmap),
            error = function(e) list(.failed = TRUE,
                sample = runnable$Sample_ID[n], msg = conditionMessage(e)))
    }, count_cores)
    names(cd) = runnable$Sample_ID
    cd
})

# ---- Structural: count tables that did not parse ---------------------
# An empty or truncated count file yields character ref/alt columns
# rather than numeric, which then blows up on `x$ref + x$alt`. No data,
# not a threshold. -Claude generated
good_cd = vapply(countdfs, function(x)
    ok(x) && is.data.frame(x) && is.numeric(x$ref) && is.numeric(x$alt),
    logical(1))
flag_add(names(countdfs)[!good_cd], 'countfile_malformed')
countdfs = countdfs[good_cd]


# ── 5. Depth ──────────────────────────────────────────────────────────
# Mean per-site depth over segregating sites. Note this is the
# parent-dependent depth — it is not the raw coverage qc_pipeline.R
# reports, because the site set depends on which parents were assigned.
#
# V5 change: this FLAGS, it does not filter. Everything still gets an
# AFD attempt. -JBloom (metric) / -Claude generated (flag)
tdepth = vapply(countdfs, function(x) sum(x$ref + x$alt) / nrow(x), 0)
key_set(names(tdepth), 'mean_depth_seg', unname(tdepth))
flag_add(names(countdfs)[tdepth <= min_depth], 'low_depth')
cat('depth: ', sum(tdepth > min_depth), '/', length(countdfs),
    ' above min_depth (', min_depth, '); the rest are flagged, not dropped\n',
    sep = '')


# ── 6. Smoothed allele frequency (AFD) ────────────────────────────────
# Per sample: bin counts, loess per chromosome with AICc-chosen span,
# return smoothed AF and its SE. Run on every sample with a count table,
# including thin ones — a sample too thin to fit fails here, and
# 'afd_failed' is a more honest label than a depth cut. -JBloom
afds = cache('afds', {
    a = par_lapply(names(countdfs), function(snn) {
        tryCatch(
            calcAFD(countdfs[[snn]], experiment.name = snn,
                    sample.size = afd_sample_size,
                    sel.strength = afd_sel_strength),
            error = function(e) list(.failed = TRUE, sample = snn,
                                     msg = conditionMessage(e)))
    }, ncores)
    names(a) = names(countdfs)
    a
})

# ---- Structural: AFD that did not converge ---------------------------
# doBlockLoess fails when a chromosome has too little signal to fit.
# Usually thin coverage — but not always, and a good-depth sample that
# still fails is worth looking at, which is why depth is a separate flag.
good_afd = vapply(afds, function(x)
    ok(x) && is.data.frame(x) && 'chrom' %in% colnames(x), logical(1))
flag_add(names(afds)[!good_afd], 'afd_failed')
afds = afds[good_afd]
keep = names(afds)
samples_kept = samples[samples$Sample_ID %in% keep, , drop = FALSE]
cat('samples with a usable AFD:', length(keep), '/', nrow(samples), '\n')


# ── 6b. AF baseline diagnostics ───────────────────────────────────────
# Computed here rather than in a downstream script so the results key is
# complete in one pass, and so skew is a config threshold in §0 like any
# other rather than a second file to keep in sync. -Claude generated
.afcol = function(a) {
    en = attr(a, 'experiment.name')
    cn = paste0('afd_', en)
    if (!is.null(en) && cn %in% names(a)) return(cn)
    grep('^afd_', names(a), value = TRUE)[1]
}
.depthcols = function(a) {
    en = attr(a, 'experiment.name')
    p1c = paste0('p1_', en); p2c = paste0('p2_', en)
    if (!all(c(p1c, p2c) %in% names(a))) {
        p1c = grep('^p1_', names(a), value = TRUE)[1]
        p2c = grep('^p2_', names(a), value = TRUE)[1]
    }
    c(p1c, p2c)
}

stopifnot(length(keep) > 0)
afstat = do.call(rbind, lapply(keep, function(sid) {
    a  = afds[[sid]]
    af = a[[.afcol(a)]]
    dc = .depthcols(a)
    dep = if (!any(is.na(dc))) a[[dc[1]]] + a[[dc[2]]] else NA_real_
    data.frame(Sample_ID = sid,
               n_sites   = nrow(a),
               med_depth = median(dep, na.rm = TRUE),
               af_median = median(af, na.rm = TRUE),
               af_min    = min(af, na.rm = TRUE),
               af_max    = max(af, na.rm = TRUE),
               frac_out  = mean(af < 0 | af > 1, na.rm = TRUE),
               d_half    = mean(abs(af - 0.5), na.rm = TRUE),
               stringsAsFactors = FALSE)
}))
# excess over the binomial noise floor (0.399/sqrt(n)); raw d_half is
# depth-confounded at Spearman -0.75 and must never be used bare.
afstat$excess = afstat$d_half - 0.399 / sqrt(pmax(afstat$med_depth, 1))
afstat$skew   = abs(afstat$af_median - 0.5)
# 1/3 and 2/3 are the 2:1-dosage expectations; closer to one of those
# than to 1/2 is the triploid signature.
afstat$d_third = pmin(abs(afstat$af_median - 1/3), abs(afstat$af_median - 2/3))
afstat$triploid_like = afstat$d_third < afstat$skew

for (cn in setdiff(names(afstat), 'Sample_ID'))
    key_set(afstat$Sample_ID, cn, afstat[[cn]])

flag_add(afstat$Sample_ID[afstat$frac_out > diverged_frac], 'afd_diverged')
flag_add(afstat$Sample_ID[afstat$med_depth < min_med_depth], 'shallow')
flag_add(afstat$Sample_ID[afstat$excess   > max_excess],    'poor_fit')

# Skew is a property of the WELL, not the sample: a well's conditions
# are independent post-split cultures and show the same offset in each.
# Flagging per sample would split a well across databins and break the
# replicate structure the whole design rests on. Diverged samples are
# excluded from the well median because their af_median is meaningless.
if (!is.na(skew_cut)) {
    aw   = key$Wells[match(afstat$Sample_ID, key$Sample_ID)]
    usable = afstat$frac_out <= diverged_frac & afstat$med_depth >= min_med_depth
    wsk  = tapply(afstat$skew[usable], aw[usable], median, na.rm = TRUE)
    bad  = names(wsk)[!is.na(wsk) & wsk >= skew_cut]
    cat('  af_skew: ', length(bad), ' of ', length(wsk), ' wells over ',
        skew_cut, ': ', paste(sort(bad), collapse = ', '), '\n', sep = '')
    flag_add(key$Sample_ID[key$Wells %in% bad], 'af_skew')
} else {
    key$flag_af_skew = FALSE
    cat('  skew_cut is NA — no skew flag applied.\n')
}


# ── 7. databin ────────────────────────────────────────────────────────
# Resolved here because everything from this point on is either labelled
# by it or filtered on it. Order matters: no_data is checked first, so a
# thin sample whose AFD also failed reads as no_data with
# flag_reasons 'low_depth;afd_failed' rather than as a flagged result
# that does not exist. -Claude generated
flag_cols = grep('^flag_', names(key), value = TRUE)
for (cn in flag_cols) key[[cn]][is.na(key[[cn]])] = FALSE
fm = as.matrix(as.data.frame(key[flag_cols]))
key$n_flags      = rowSums(fm)
key$flag_reasons = apply(fm, 1, function(r)
    paste(sub('^flag_', '', flag_cols)[which(r)], collapse = ';'))
key$has_result   = key$Sample_ID %in% keep
key$databin      = ifelse(!key$has_result, 'no_data',
                   ifelse(key$n_flags > 0, 'flagged', 'main'))

bin_of  = function(sid) key$databin[match(sid, key$Sample_ID)]
# A contrast is only as clean as its worse half.
bin_pair = function(L, R) {
    b = c(bin_of(L), bin_of(R))
    if (any(b == 'no_data')) 'no_data' else if (any(b == 'flagged')) 'flagged' else 'main'
}

cat('\n══ databin ══\n'); print(table(key$databin))
cat('\n══ flags (samples may carry several) ══\n')
print(sort(colSums(fm), decreasing = TRUE))

main_ids = key$Sample_ID[key$databin == 'main']
cat('\nmain dataset: ', length(main_ids), ' samples, ',
    length(unique(key$Wells[key$databin == 'main'])), ' wells, ',
    length(unique(key$cross[key$databin == 'main'])), ' crosses\n', sep = '')

p1_of = function(sid) samples$parent_a[match(sid, samples$Sample_ID)]
p2_of = function(sid) samples$parent_alpha[match(sid, samples$Sample_ID)]
meta_of = function(sid) samples[match(sid, samples$Sample_ID), , drop = FALSE]


# ── 8. Contrasts ──────────────────────────────────────────────────────
# (i)  fitness: every sample with an AFD against the constant R = 0.5,
#      the one-sided analogue of a no-selection null. -JBloom
#
#      R stays a hardcoded 0.5 DELIBERATELY. Passing each sample's own
#      genome-wide median instead would silence the skew — but it would
#      do so identically whether the cause is triploidy, a residual
#      wrong parent, or contamination, and a wrong parent would come out
#      looking clean. The skew is a flag to record, not a nuisance to
#      normalise away. See flag_af_skew in §6b. -Claude generated
# (ii) stress vs control: every non-control sample against the control
#      sample sharing its (parent_a, parent_alpha, Wells). There is no
#      T0 in this design, so the glucose plate is the reference.
fit_results = cache('fit_results', {
    r = par_lapply(keep, function(sid) {
        tryCatch(
            calcContrastStats(results = list(afds[[sid]]),
                              L = paste0('_', sid), R = 0.5),
            error = function(e) list(.failed = TRUE, sample = sid,
                                     msg = conditionMessage(e)))
    }, ncores)
    names(r) = keep
    r
})

mk = function(df) do.call(paste, c(unname(as.list(df[match_keys])), sep = '|'))
treated  = samples_kept[samples_kept$condition != control_cond, , drop = FALSE]
controls = samples_kept[samples_kept$condition == control_cond, , drop = FALSE]
treated$.k = mk(treated); controls$.k = mk(controls)

pairs = lapply(seq_len(nrow(treated)), function(i) {
    ctl = controls$Sample_ID[controls$.k == treated$.k[i]]
    if (!length(ctl)) return(NULL)
    list(L = treated$Sample_ID[i], R = ctl[1],
         condition = treated$condition[i],
         cross     = treated$cross[i],
         well      = treated$Wells[i])
})
unmatched = treated$Sample_ID[vapply(pairs, is.null, logical(1))]
pairs = pairs[!vapply(pairs, is.null, logical(1))]
# Recorded as a column, NOT a flag: having no control blocks the stress
# contrast only. The sample's fitness result is unaffected, so it has no
# business being pushed out of the main databin over it.
key$stress_status = NA_character_
key$stress_status[key$has_result] = 'na_control_sample'   # i.e. a Glu sample
key_set(treated$Sample_ID, 'stress_status', 'paired')
key_set(unmatched, 'stress_status', 'no_matched_control')
cat('stress-vs-', control_cond, ' pairs: ', length(pairs),
    '  (', length(unmatched), ' unmatched)\n', sep = '')

drug_results = cache('drug_results', {
    r = par_lapply(pairs, function(pair) {
        tryCatch({
            res = calcContrastStats(
                results = list(afds[[pair$L]], afds[[pair$R]]),
                L = paste0('_', pair$L), R = paste0('_', pair$R))
            list(meta = pair, result = res)
        }, error = function(e) list(.failed = TRUE, meta = pair,
                                   msg = conditionMessage(e)))
    }, ncores)
    names(r) = vapply(pairs, function(p) paste0(p$L, '_vs_', p$R),
                      character(1))
    r
})

cat('fitness contrasts: ', length(fit_results), ' (',
    sum(!vapply(fit_results, ok, logical(1))), ' failed)\n', sep = '')
cat('stress contrasts:  ', length(drug_results), ' (',
    sum(!vapply(drug_results, ok, logical(1))), ' failed)\n', sep = '')


# ── 9. QTL peaks ──────────────────────────────────────────────────────
# Called ONCE over everything and tagged with databin. The fork is a
# filter at the point of use, not a second pipeline — two copies of the
# peak-calling code would drift. -JBloom (findQTLPeaks) / -Claude (tag)
.tag = function(gr, fields) {
    if (is.null(gr) || length(gr) == 0) return(GRanges())
    for (k in names(fields)) mcols(gr)[[k]] = fields[[k]]
    gr
}

peaks_all = cache('peaks', {
    fp = lapply(names(fit_results), function(sid) {
        if (!ok(fit_results[[sid]])) return(GRanges())
        sm = meta_of(sid)
        .tag(findQTLPeaks(fit_results[[sid]],
                LOD.threshold = LOD_threshold, LOD.drop = LOD_drop,
                LOD.drop.frac = LOD_drop_frac,
                min.interval = min_interval, max.interval = max_interval,
                p1.name = p1_of(sid), p2.name = p2_of(sid)),
            list(sample = sid, contrast = 'fitness',
                 condition = sm$condition, cross = sm$cross,
                 well = sm$Wells, databin = bin_of(sid)))
    })
    dp = lapply(names(drug_results), function(k) {
        it = drug_results[[k]]
        if (!ok(it)) return(GRanges())
        .tag(findQTLPeaks(it$result,
                LOD.threshold = LOD_threshold, LOD.drop = LOD_drop,
                LOD.drop.frac = LOD_drop_frac,
                min.interval = min_interval, max.interval = max_interval,
                p1.name = p1_of(it$meta$L), p2.name = p2_of(it$meta$L)),
            list(L = it$meta$L, R = it$meta$R,
                 contrast = paste0('stress_vs_', control_cond),
                 condition = it$meta$condition, cross = it$meta$cross,
                 well = it$meta$well,
                 databin = bin_pair(it$meta$L, it$meta$R)))
    })
    list(fit = unlist(GRangesList(fp)), drug = unlist(GRangesList(dp)))
})
fit_peaks  = peaks_all$fit
drug_peaks = peaks_all$drug

# The main-databin views. Everything downstream of here uses these;
# the CSVs written below carry both, labelled.
fit_main  = if (length(fit_peaks))  fit_peaks[mcols(fit_peaks)$databin == 'main']   else GRanges()
drug_main = if (length(drug_peaks)) drug_peaks[mcols(drug_peaks)$databin == 'main'] else GRanges()

cat('fitness peaks: ', length(fit_peaks), ' (', length(fit_main), ' main)',
    ' | stress peaks: ', length(drug_peaks), ' (', length(drug_main),
    ' main)\n', sep = '')

peak_table = function(gr) {
    if (length(gr) == 0) return(data.frame())
    df = as.data.frame(gr)
    if ('peak.LOD' %in% names(df)) df = df[order(-df$peak.LOD), ]
    df
}
write.csv(peak_table(fit_peaks),  file.path(out_dir, 'fit_peaks.csv'),
          row.names = FALSE)
write.csv(peak_table(drug_peaks), file.path(out_dir, 'stress_peaks.csv'),
          row.names = FALSE)

# Peak counts back into the key, split by databin so a reader can see
# that flagged samples call more peaks rather than having to infer it.
key$n_peaks_fit    = NA_integer_
key$n_peaks_stress = NA_integer_
if (length(fit_peaks)) {
    n = table(mcols(fit_peaks)$sample)
    key_set(names(n), 'n_peaks_fit', as.integer(n))
}
if (length(drug_peaks)) {
    n = table(mcols(drug_peaks)$L)
    key_set(names(n), 'n_peaks_stress', as.integer(n))
}
# A sample that ran and called nothing is a 0, not an NA. A sample that
# never ran stays NA -- the two must not be averaged together.
key$n_peaks_fit[key$has_result & is.na(key$n_peaks_fit)] = 0L
key$n_peaks_stress[key$stress_status %in% 'paired' &
                   is.na(key$n_peaks_stress)] = 0L

# Recurrent intervals: positions called in several independent pairs.
# Main only — a recurrent interval driven by flagged samples is not
# evidence of anything.
if (length(drug_main)) {
    hits = reduce(drug_main, with.revmap = TRUE)
    mcols(hits)$n_pairs = lengths(mcols(hits)$revmap)
    recurrent = hits[mcols(hits)$n_pairs >= 3]
    cat('recurrent intervals (>=3 pairs, main):', length(recurrent), '\n')
    write.csv(data.frame(chrom   = as.character(seqnames(recurrent)),
                         start   = start(recurrent),
                         end     = end(recurrent),
                         n_pairs = mcols(recurrent)$n_pairs),
              file.path(out_dir, 'recurrent_intervals.csv'), row.names = FALSE)
}
peaks_by_cross     = if (length(drug_main)) split(drug_main, mcols(drug_main)$cross) else list()
fit_peaks_by_cross = if (length(fit_main))  split(fit_main,  mcols(fit_main)$cross)  else list()


# ── 10. Figures ───────────────────────────────────────────────────────
# Four-panel per stress/control pair, two-panel per fitness sample, and
# one summary grid per (condition, cross). Every write goes through
# fig_path(), so flagged samples land in figures/flagged/ with the same
# subtree shape and nothing is silently absent.
# -JBloom (panels) / -Claude generated (batch driver, routing)
plotPairFourPanel = function(k) {
    it = drug_results[[k]]
    if (!ok(it)) return(NULL)
    L = it$meta$L; R = it$meta$R
    p1 = p1_of(L); p2 = p2_of(L)
    pA = plotIndividualExperiment(afds[[L]], L, p1.name = p1, p2.name = p2)
    pB = plotIndividualExperiment(afds[[R]], R, p1.name = p1, p2.name = p2)
    pC = plotContrast(it$result, L, R, p1.name = p1, p2.name = p2)
    bin = bin_pair(L, R)
    pD = plotSummarySignedShade(it$result, p1.name = p1, p2.name = p2) +
         ggtitle(sprintf('%s vs %s — %s (well %s)  [%s vs %s]%s',
             it$meta$condition, control_cond, it$meta$cross,
             it$meta$well, L, R,
             if (bin == 'main') '' else paste0('  FLAGGED: ',
                 key$flag_reasons[match(L, key$Sample_ID)])))
    (pA | pB) / (pC | pD)
}

plotFitnessTwoPanel = function(sid) {
    res = fit_results[[sid]]
    if (!ok(res)) return(NULL)
    sm = meta_of(sid)
    fr = key$flag_reasons[match(sid, key$Sample_ID)]
    pA = plotIndividualExperiment(afds[[sid]], sid,
             p1.name = p1_of(sid), p2.name = p2_of(sid))
    pB = plotSummarySignedShade(res, p1.name = p1_of(sid),
             p2.name = p2_of(sid)) +
         ggtitle(sprintf('%s — fitness vs 0.5 (%s, %s, well %s)%s',
                         sid, sm$condition, sm$cross, sm$Wells,
                         if (nzchar(fr)) paste0('  FLAGGED: ', fr) else ''))
    pA / pB
}

if (running('figures')) {
    for (k in names(drug_results)) {
        it = drug_results[[k]]; if (!ok(it)) next
        p = plotPairFourPanel(k); if (is.null(p)) next
        ggsave(fig_path(bin_pair(it$meta$L, it$meta$R),
                   'stress_vs_control_fourpanel',
                   sprintf('%s-%s_w%s-%sm%s.png',
                   it$meta$cross, it$meta$condition, it$meta$well,
                   .short(it$meta$L), .short(it$meta$R))),
               plot = p, width = 60, height = 30, units = 'cm')
    }
    for (sid in names(fit_results)) {
        if (!ok(fit_results[[sid]])) next
        p = plotFitnessTwoPanel(sid); if (is.null(p)) next
        sm = meta_of(sid)
        ggsave(fig_path(bin_of(sid), 'fitness_vs_0.5_twopanel',
                   sprintf('%s-%s-%s_w%s.png',
                   .short(sid), sm$cross, sm$condition, sm$Wells)),
               plot = p, width = 60, height = 30, units = 'cm')
    }

    .save_grid = function(panels, path, ncol_grid = 2) {
        if (!length(panels)) return(invisible(NULL))
        ggsave(path, plot = wrap_plots(panels, ncol = ncol_grid),
               width = 60, height = 15 * ceiling(length(panels) / ncol_grid),
               units = 'cm', limitsize = FALSE)
    }
    # Summary grids are per (condition, cross) AND per databin, so a
    # grid never mixes main and flagged wells in one figure — the whole
    # point of the grid is comparing replicate wells that are alike.
    dm = do.call(rbind, lapply(names(drug_results), function(k) {
        it = drug_results[[k]]; if (!ok(it)) return(NULL)
        data.frame(key = k, condition = it$meta$condition,
                   cross = it$meta$cross, well = it$meta$well,
                   bin = bin_pair(it$meta$L, it$meta$R),
                   stringsAsFactors = FALSE)
    }))
    if (!is.null(dm)) {
        combos = unique(dm[, c('condition', 'cross', 'bin')])
        for (i in seq_len(nrow(combos))) {
            sub = dm[dm$condition == combos$condition[i] &
                     dm$cross == combos$cross[i] & dm$bin == combos$bin[i], ]
            sub = sub[order(sub$well), , drop = FALSE]
            .save_grid(lapply(sub$key, function(k) {
                it = drug_results[[k]]
                plotSummarySignedShade(it$result, p1.name = p1_of(it$meta$L),
                    p2.name = p2_of(it$meta$L)) +
                    ggtitle(sprintf('%s — %s (well %s)', it$meta$condition,
                                    it$meta$cross, it$meta$well))
            }), fig_path(combos$bin[i], 'summaries',
                 sprintf('%s-%s_stress_summaries.png',
                         combos$cross[i], combos$condition[i])))
        }
    }
    fmeta = samples_kept[samples_kept$Sample_ID %in%
                names(fit_results)[vapply(fit_results, ok, logical(1))], ]
    fmeta$bin = bin_of(fmeta$Sample_ID)
    combos = unique(fmeta[, c('condition', 'cross', 'bin')])
    for (i in seq_len(nrow(combos))) {
        sub = fmeta[fmeta$condition == combos$condition[i] &
                    fmeta$cross == combos$cross[i] & fmeta$bin == combos$bin[i], ]
        sub = sub[order(sub$Wells), , drop = FALSE]
        .save_grid(Map(function(sid, well)
            plotSummarySignedShade(fit_results[[sid]], p1.name = p1_of(sid),
                p2.name = p2_of(sid)) +
                ggtitle(sprintf('%s — %s, %s (well %s)', sid,
                        combos$condition[i], combos$cross[i], well)),
            sub$Sample_ID, sub$Wells),
            fig_path(combos$bin[i], 'summaries',
                sprintf('%s-%s_fitness_summaries.png',
                        combos$cross[i], combos$condition[i])))
    }
}


# ── 11. Reproducibility within a cross ────────────────────────────────
# Main databin only, by design. Replicate wells are the unit, and a
# flagged well dragged into the comparison would depress a cross's
# agreement for a reason that has nothing to do with reproducibility.
have_ch = requireNamespace('ComplexHeatmap', quietly = TRUE) &&
          requireNamespace('circlize', quietly = TRUE)

.sample_depth = function(sid) {
    a = afds[[sid]]
    if (is.null(a) || !ok(a)) return(NA_real_)
    dc = .depthcols(a)
    if (any(is.na(dc))) return(NA_real_)
    median(a[[dc[1]]] + a[[dc[2]]], na.rm = TRUE)
}

# Bin signed -log10(p) onto a common genomic grid; within a bin keep the
# largest-magnitude value, sign preserved. Infinite LODs (AFD pinned at
# a boundary so SE -> 0) are clamped before binning so they cannot
# swallow a bin.
.signalMatrix = function(rows, get_result, bin_size = 5000, inf_cap = 50) {
    if (!nrow(rows)) return(NULL)
    ex = get_result(rows$id[1])
    grid = ex %>% mutate(bin = floor(physical.position / bin_size)) %>%
        distinct(chrom, bin) %>% arrange(chrom, bin)
    col_key = paste(grid$chrom, grid$bin, sep = ':')
    M = matrix(NA_real_, nrow(rows), length(col_key),
               dimnames = list(rows$row_label, col_key))
    for (i in seq_len(nrow(rows))) {
        r  = get_result(rows$id[i])
        sl = sign(r$contrast.beta) * (-log10(r$p))
        im = is.infinite(sl); sl[im] = sign(r$contrast.beta[im]) * inf_cap
        b = data.frame(chrom = as.character(r$chrom),
                       bin = floor(r$physical.position / bin_size),
                       sl = sl) %>%
            group_by(chrom, bin) %>%
            summarise(sl = { v = sl[is.finite(sl)]
                             if (!length(v)) NA_real_ else v[which.max(abs(v))] },
                      .groups = 'drop')
        M[i, paste(b$chrom, b$bin, sep = ':')] = b$sl
    }
    attr(M, 'chrom') = grid$chrom
    M
}

# `name` must stay ASCII: a newline in it becomes part of the viewport
# path and downViewport() then fails on every heatmap. The two-line
# legend title belongs in heatmap_legend_param. clamp stays at 20 —
# raising it to 50 washes out typical peaks (3-15) against the ramp.
.drawHeatmap = function(M, depth, title, path, threshold = 3, clamp = 20) {
    if (is.null(M) || !have_ch) return(invisible(NULL))
    hi = which(abs(M) > clamp); M[hi] = sign(M[hi]) * clamp
    M[which(abs(M) < threshold)] = NA_real_      # sub-threshold -> grey
    dr = range(depth, na.rm = TRUE)
    if (any(!is.finite(dr)) || diff(dr) == 0)
        dr = c(0, max(c(depth, 1), na.rm = TRUE) + 1)
    ht = ComplexHeatmap::Heatmap(M, name = 'signed_log10p',
        heatmap_legend_param = list(title = 'signed\n-log10(p)'),
        col = circlize::colorRamp2(c(-clamp, 0, clamp),
                                   c('royalblue', 'white', 'firebrick')),
        na_col = 'grey90',
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE,
        column_split = factor(attr(M, 'chrom'), levels = uchr),
        column_gap = grid::unit(0.5, 'mm'),
        column_title_gp = grid::gpar(fontsize = 7),
        row_names_gp = grid::gpar(fontsize = 6), row_names_side = 'left',
        right_annotation = ComplexHeatmap::rowAnnotation(
            depth = depth,
            col = list(depth = circlize::colorRamp2(dr, c('grey95','black'))),
            annotation_name_gp = grid::gpar(fontsize = 8),
            width = grid::unit(6, 'mm')),
        column_title = title)
    png(path, width = 20, height = max(5, 0.18 * nrow(M) + 2),
        units = 'in', res = 150, type = 'cairo')
    ComplexHeatmap::draw(ht, merge_legend = TRUE)
    dev.off()
    invisible(path)
}

.dotplot = function(gr, id_col, title, path, lod_cap = 50) {
    if (!length(gr)) return(invisible(NULL))
    df = as.data.frame(gr)
    df$row_label = sprintf('%s w%s [%s]', df$condition, df$well,
                           .short(df[[id_col]]))
    df = df[order(df$condition, df$well), , drop = FALSE]
    df$row_label = factor(df$row_label, levels = unique(df$row_label))
    df$seqnames  = factor(as.character(df$seqnames), levels = uchr)
    df$dir = ifelse(df$peak.beta > 0, '+ (parent_a)', '- (parent_alpha)')
    df$lodc = pmin(df$peak.LOD, lod_cap)
    df$lodc[!is.finite(df$lodc)] = lod_cap

    p_dot = ggplot(df, aes(peak.pos, row_label, size = lodc, color = dir)) +
        geom_point(alpha = 0.8) +
        facet_grid(. ~ seqnames, scales = 'free_x', space = 'free_x',
                   switch = 'x') +
        scale_color_manual(values = c('+ (parent_a)' = 'firebrick',
                                      '- (parent_alpha)' = 'royalblue')) +
        scale_size_continuous(range = c(0.5, 5),
                              name = sprintf('LOD (cap %d)', lod_cap)) +
        theme_bw(base_size = 8) +
        theme(panel.spacing.x = grid::unit(0.05, 'lines'),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              strip.text = element_text(size = 6),
              strip.placement = 'outside',
              axis.text.y = element_text(size = 6)) +
        labs(x = NULL, y = NULL, title = title)

    cov = unique(df[, c('row_label', id_col)])
    cov$depth = vapply(cov[[id_col]], .sample_depth, 0)
    p_cov = ggplot(cov, aes(1, row_label, fill = depth)) + geom_tile() +
        scale_fill_gradient(low = 'grey95', high = 'black',
                            name = 'med depth') +
        theme_void(base_size = 8) +
        theme(legend.position = 'right', axis.text = element_blank())

    ggsave(path, plot = p_dot + p_cov + plot_layout(widths = c(20, 1)),
           width = 60, height = max(15, 0.5 * nlevels(df$row_label) + 4),
           units = 'cm', limitsize = FALSE)
    invisible(path)
}

if (running('repro')) {
    if (!have_ch)
        cat('!! ComplexHeatmap/circlize unavailable — dot-plots only\n')
    smain = samples_kept[samples_kept$Sample_ID %in% main_ids, , drop = FALSE]

    for (cr in unique(smain$cross)) {
        ks = names(drug_results)[vapply(drug_results, function(x)
                 ok(x) && x$meta$cross == cr &&
                 bin_pair(x$meta$L, x$meta$R) == 'main', logical(1))]
        if (length(ks)) {
            rows = data.frame(id = ks,
                condition = vapply(ks, function(k) drug_results[[k]]$meta$condition, ''),
                well      = vapply(ks, function(k) drug_results[[k]]$meta$well, ''),
                L         = vapply(ks, function(k) drug_results[[k]]$meta$L, ''),
                stringsAsFactors = FALSE)
            rows = rows[order(rows$condition, rows$well), , drop = FALSE]
            rows$row_label = sprintf('%s w%s [%s]', rows$condition,
                                     rows$well, .short(rows$L))
            .drawHeatmap(.signalMatrix(rows, function(k) drug_results[[k]]$result),
                vapply(rows$L, .sample_depth, 0),
                sprintf('%s — stress vs %s (|signed -log10 p| >= 3, 5 kb bins)',
                        cr, control_cond),
                fig_path('main', 'reproducibility',
                         sprintf('%s_stress_heatmap.png', cr)))
            .dotplot(drug_main[mcols(drug_main)$cross == cr], 'L',
                sprintf('QTL peaks — %s stress vs %s', cr, control_cond),
                fig_path('main', 'reproducibility',
                         sprintf('%s_stress_dotplot.png', cr)))
        }
        sf = smain[smain$cross == cr & smain$Sample_ID %in%
                 names(fit_results)[vapply(fit_results, ok, logical(1))], ]
        if (nrow(sf)) {
            sf = sf[order(sf$condition, sf$Wells), , drop = FALSE]
            rows = data.frame(id = sf$Sample_ID,
                row_label = sprintf('%s w%s [%s]', sf$condition, sf$Wells,
                                    .short(sf$Sample_ID)),
                stringsAsFactors = FALSE)
            .drawHeatmap(.signalMatrix(rows, function(s) fit_results[[s]]),
                vapply(sf$Sample_ID, .sample_depth, 0),
                sprintf('%s — fitness vs 0.5 (|signed -log10 p| >= 3, 5 kb bins)', cr),
                fig_path('main', 'reproducibility',
                         sprintf('%s_fitness_heatmap.png', cr)))
            .dotplot(fit_main[mcols(fit_main)$cross == cr], 'sample',
                sprintf('QTL peaks — %s fitness vs 0.5', cr),
                fig_path('main', 'reproducibility',
                         sprintf('%s_fitness_dotplot.png', cr)))
        }
    }

    # Cross x condition replicate agreement, main only. Reported as
    # ABSOLUTE n_shared as well as the fraction: frac_shared is
    # mechanically depressed by peak count (Spearman -0.588 vs
    # n_intervals) and has no permutation null yet, so it must not be
    # compared across groups that differ in peak count. -Claude generated
    rep = list()
    fpm = as.data.frame(fit_main)
    for (cc in unique(fpm$cross)) for (cd in unique(fpm$condition[fpm$cross == cc])) {
        sub = fpm[fpm$cross == cc & fpm$condition == cd, ]
        nw = length(unique(sub$well)); if (nw < 2) next
        gr = GRanges(sub$seqnames, IRanges(sub$start, sub$end))
        rd = reduce(gr, with.revmap = TRUE)
        nwell = vapply(mcols(rd)$revmap, function(i) length(unique(sub$well[i])), 0L)
        rep[[length(rep)+1]] = data.frame(cross = cc, condition = cd,
            n_wells = nw, n_intervals = length(rd), n_shared = sum(nwell > 1),
            frac_shared = round(mean(nwell > 1), 3))
    }
    if (length(rep)) {
        R = do.call(rbind, rep)
        R = R[order(-R$n_shared), ]
        write.csv(R, file.path(qc_out_dir, 'reproducibility_main.csv'),
                  row.names = FALSE)
        cat('\n── replicate reproducibility, main databin ──\n')
        print(R, row.names = FALSE)
    }
}


# ── 12. Segregating variants under each interval ──────────────────────
# Main databin only. extractQTLVariants queries strain genotypes, so
# intervals must be split by cross. suppressMessages matters: it emits
# three message() calls per peak, and on thousands of peaks the stderr
# flushes dominate runtime on an SGE shell. -Claude generated
if (running('variants')) {
    lookup = qs2::qs_read(vlookup)
    jobs = c(
        setNames(lapply(names(peaks_by_cross), function(cr)
            list(type = 'stress', pk = peaks_by_cross[[cr]])),
            paste0('stress|', names(peaks_by_cross))),
        setNames(lapply(names(fit_peaks_by_cross), function(cr)
            list(type = 'fitness', pk = fit_peaks_by_cross[[cr]])),
            paste0('fitness|', names(fit_peaks_by_cross))))
    jobs = jobs[vapply(jobs, function(j) length(j$pk) > 0, logical(1))]

    vout = cache('qtl_variants', par_lapply(jobs, function(j)
        tryCatch(suppressMessages(extractQTLVariants(j$pk, lookup)),
                 error = function(e) list(.failed = TRUE,
                                          msg = conditionMessage(e))),
        ncores))

    nf = sum(!vapply(vout, ok, logical(1)))
    if (nf) cat('extractQTLVariants failures:', nf, '/', length(vout), '\n')

    .by_type = function(tag) {
        k = grep(paste0('^', tag, '\\|'), names(vout), value = TRUE)
        o = vout[k]; names(o) = sub(paste0('^', tag, '\\|'), '', k); o
    }
    for (tag in c('stress', 'fitness')) {
        v = Filter(ok, .by_type(tag))
        if (!length(v)) next
        write.csv(flattenQTLVariants(v, lookup),
            file.path(out_dir, sprintf('qtl_variants_%s.csv', tag)),
            row.names = FALSE)
    }
    rm(lookup); invisible(gc())
}


# ── 13. Parental CNVs under each interval ─────────────────────────────
# Experimental — CNVnator calls for the parent strains of each cross.
# Main databin only. -JBloom
if (running('cnv')) {
    cnv_gr = loadCNVnator(cnv.file)
    for (tag in c('stress', 'fitness')) {
        by_cross = if (tag == 'stress') peaks_by_cross else fit_peaks_by_cross
        res = list()
        for (cr in names(by_cross)) {
            if (!length(by_cross[[cr]])) next
            res[[cr]] = tryCatch(extractQTLCNVs(by_cross[[cr]], cnv_gr),
                                 error = function(e) NULL)
        }
        res = Filter(Negate(is.null), res)
        if (!length(res)) next
        write.csv(flattenQTLCNVs(res),
            file.path(out_dir, sprintf('qtl_cnvs_%s.csv', tag)),
            row.names = FALSE)
    }
    rm(cnv_gr); invisible(gc())
}


# ── 14. Gene-set enrichment in QTL intervals ──────────────────────────
# Circular-genome permutation: shift all intervals per chromosome by a
# random offset, preserving interval count and sizes, and count gene-set
# overlaps to build a null. Main databin only. NOTE: `subset` here is
# the stratification of the TEST (overall / per condition / per cross),
# unrelated to databin. -JBloom (test) / -Claude generated
if (running('enrich') && !is.null(gene_sets) && length(drug_main)) {
    run_one = function(gr, gs) {
        if (!length(gr)) return(list(observed = NA_real_,
            n_gene_set = length(gs), p_value = NA_real_,
            observed_genes = character(0), null_distribution = numeric(0)))
        testQTLGeneEnrichment(gr, gs, n_perm = n_perm, gff_path = gff_file)
    }
    tests = list()
    add = function(scope, subset, gr, gs) {
        tests[[paste(scope, subset, sep = '|')]] <<-
            list(scope = scope, subset = subset, n_peaks = length(gr),
                 result = run_one(gr, gs))
    }
    for (nm in names(gene_sets)) {
        add('overall_stress',  nm, drug_main, gene_sets[[nm]])
        add('overall_fitness', nm, fit_main,  gene_sets[[nm]])
        for (cn in unique(mcols(drug_main)$condition))
            add('per_condition', paste(cn, nm, sep = ':'),
                drug_main[mcols(drug_main)$condition == cn], gene_sets[[nm]])
        for (cr in unique(mcols(drug_main)$cross))
            add('per_cross', paste(cr, nm, sep = ':'),
                drug_main[mcols(drug_main)$cross == cr], gene_sets[[nm]])
    }
    qs2::qs_save(tests, file.path(out_dir, 'enrichment_tests.qs'))

    edf = do.call(rbind, lapply(tests, function(t) with(t, data.frame(
        scope = scope, subset = subset, n_peaks = n_peaks,
        n_genes = result$n_gene_set, observed = result$observed,
        p_value = result$p_value,
        genes = paste(result$observed_genes, collapse = ',')))))
    edf = edf[order(edf$p_value), ]
    write.csv(edf, file.path(out_dir, 'enrichment_summary.csv'),
              row.names = FALSE)
    print(utils::head(edf, 20))

    fdf = do.call(rbind, lapply(tests, function(t) {
        nd = t$result$null_distribution
        if (!length(nd)) return(NULL)
        q = quantile(nd, c(0.025, 0.975), na.rm = TRUE)
        data.frame(scope = t$scope, subset = t$subset,
                   n_peaks = t$n_peaks, observed = t$result$observed,
                   null_mean = mean(nd, na.rm = TRUE),
                   null_lo = q[[1]], null_hi = q[[2]],
                   p_value = t$result$p_value)
    }))
    if (!is.null(fdf)) {
        fdf$label = factor(sprintf('%s - %s (n=%d)', fdf$scope, fdf$subset,
                                   fdf$n_peaks))
        p_forest = ggplot(fdf, aes(y = label)) +
            geom_segment(aes(x = null_lo, xend = null_hi, yend = label),
                         color = 'grey60', linewidth = 1.2) +
            geom_point(aes(x = null_mean), color = 'grey40', shape = 18,
                       size = 2.5) +
            geom_point(aes(x = observed, color = scope), size = 3) +
            geom_text(aes(x = observed, label = sprintf(' p=%.3g', p_value)),
                      hjust = 0, size = 3) +
            labs(x = 'Gene-set members overlapping QTL intervals', y = NULL,
                 title = 'QTL x gene-set enrichment (main databin)',
                 subtitle = 'grey bar = 95% permutation null   diamond = null mean   point = observed',
                 color = NULL) +
            theme_minimal(base_size = 11) +
            theme(panel.grid.major.y = element_blank(),
                  legend.position = 'bottom')
        ggsave(fig_path('main', 'enrichment_forest.png'), p_forest,
               width = 24, height = max(8, 0.55 * nrow(fdf) + 4),
               units = 'cm', dpi = 150)
    }
} else if (running('enrich')) {
    cat('enrich: skipped (set `gene_sets` in §0 to enable)\n')
}


# ── 15. Results key and run manifest ──────────────────────────────────
# The key is the handshake with qc_pipeline.R and the audit trail for
# every judgement this script made. One row per keyfile row, keyfile
# columns verbatim, so merging it back onto the flat key is a no-op.
# -Claude generated
key$flag_reasons[is.na(key$flag_reasons)] = ''
readr::write_csv(key, file.path(qc_out_dir, 'results_key.csv'))

cat('\n══ results key ══\n')
print(as.data.frame(key %>% count(databin, name = 'n_samples') %>%
                    arrange(desc(n_samples))), row.names = FALSE)
cat('\nreasons, most common first:\n')
print(as.data.frame(key %>% filter(nzchar(flag_reasons)) %>%
                    count(flag_reasons, name = 'n') %>% arrange(desc(n))),
      row.names = FALSE)
cat('\npeaks per sample by databin (fitness):\n')
print(key %>% filter(has_result) %>% group_by(databin) %>%
      summarise(n = n(), median_peaks = median(n_peaks_fit, na.rm = TRUE),
                mean_peaks = round(mean(n_peaks_fit, na.rm = TRUE), 1),
                .groups = 'drop'))

manifest = list(
    run_time     = format(Sys.time()),
    version      = 'V6 (V5 analysis; run umbrella + downstream scripts)',
    run_tag      = run_tag,
    run_dir      = run_dir,
    env          = env,
    key.file     = key.file,
    cache_prefix = cache_prefix,
    n_keyfile    = nrow(key),
    n_main       = sum(key$databin == 'main'),
    n_flagged    = sum(key$databin == 'flagged'),
    n_no_data    = sum(key$databin == 'no_data'),
    n_pairs      = length(pairs),
    n_fit_peaks  = length(fit_peaks),
    n_fit_main   = length(fit_main),
    n_drug_peaks = length(drug_peaks),
    n_drug_main  = length(drug_main),
    thresholds   = list(min_depth = min_depth, LOD_threshold = LOD_threshold,
                        LOD_drop = LOD_drop, LOD_drop_frac = LOD_drop_frac,
                        min_interval = min_interval,
                        max_interval = max_interval, skew_cut = skew_cut,
                        diverged_frac = diverged_frac,
                        afd_sample_size = afd_sample_size,
                        afd_sel_strength = afd_sel_strength),
    stages_run   = stages, forced = force,
    keep_chr7    = keep_chr7, cand_tags = cand_tags,
    sessionInfo  = utils::capture.output(utils::sessionInfo()))
qs2::qs_save(manifest, file.path(out_dir, 'run_manifest.qs'))

cat('\nrun:      ', run_dir, '\noutputs:  ', out_dir, '\nfigures:  ', fig_dir,
    '\nflagged:  ', flag_dir, '\ncache:    ', cache_dir,
    '\nkey:      ', file.path(qc_out_dir, 'results_key.csv'), '\n', sep = '')


# ════════════════════════════════════════════════════════════════════
# ══ 16-17. Downstream analysis, INLINED (V6.1) ══
# ════════════════════════════════════════════════════════════════════
# Until V6.0 these two sections shelled out to qtl_landscape.R and
# qtl_candidates.R in <script_dir> via system2(). That is how the v6
# run produced no landscape figure and no candidate tables: neither
# file had been copied to Hoffman, and a missing child script is a
# one-line `!!` in the log, not an error — the run still exited
# `complete`. The code now lives here, so the driver is the unit that
# gets copied and there is nothing left to forget.
#
# The bodies below are VERBATIM copies of those two scripts from their
# first line of real code onward, at their original indentation, so a
# diff against the originals shows only the configuration headers.
# What was dropped in each case:
#   qtl_landscape.R  lines 1-38  — QTL_PEAK_DIR/QTL_OUT_DIR env vars
#                                  and find_peaks(), the two-levels-
#                                  deep search for fit_peaks.csv that
#                                  gave this project figures and tables
#                                  from different runs. The driver
#                                  knows the path, so the search is
#                                  deleted rather than pinned.
#   qtl_candidates.R lines 1-101 — QTL_* env vars and its own
#                                  paths_by_env; peak_dir, gff and
#                                  keep_chr7 now arrive as arguments
#                                  from §0, so the manifest and the
#                                  tables cannot disagree about them.
#
# Isolation used to come from a separate R process. Function scope is
# enough: neither body assigns with <<-, and none calls quit(),
# setwd() or rm(list=), so their objects die with the call. What does
# change is that an error in either one is now an error in THIS
# process — hence run_inline()'s tryCatch, which keeps the old
# contract that a failed section returns non-zero into the manifest
# and the run carries on to save it.
#
# Memory: these bodies now run alongside the driver's own objects
# instead of in a fresh process. The peak tables and the variant dump
# are the large inputs and both were already being read by the child,
# so the ceiling is roughly unchanged, but a job sized tightly with
# -l h_data has less headroom here than it did.
# -Claude generated

# Runs one inlined section with the child-process behaviour preserved:
# all of its output to <run_dir>/<label>.log, an exit code rather than
# an abort, and the last 15 log lines echoed on failure. Sinks and any
# graphics device the body leaves open are unwound on exit, so one
# section failing mid-plot cannot swallow the rest of the run's output.
run_inline = function(label, fn) {
    logf = file.path(run_dir, paste0(label, '.log'))
    cat('\n── ', label, ' ──\n', sep = '')
    t0   = Sys.time()
    dev0 = grDevices::dev.cur()
    con  = file(logf, open = 'wt')
    # Restore the sink depth we were CALLED at, not zero. Under Rscript
    # that is 0; inside an IRkernel cell evaluate() already holds an
    # output sink, and popping it sends the rest of the run's stdout to
    # the kernel's terminal instead of the notebook. -Claude generated
    sink0 = sink.number()
    msg0  = sink.number(type = 'message')
    unwind = function() {
        while (sink.number() > sink0) sink()
        if (sink.number(type = 'message') != msg0) sink(type = 'message')
        while (grDevices::dev.cur() > dev0 && grDevices::dev.cur() > 1L)
            try(grDevices::dev.off(), silent = TRUE)
    }
    # close() invalidates the connection, so a later isOpen(con) ERRORS
    # rather than returning FALSE. The success path below closes con, so
    # the guard has to tolerate an already-closed connection.
    on.exit({ unwind(); try(close(con), silent = TRUE) }, add = TRUE)
    sink(con); sink(con, type = 'message')
    # V7: catch message() explicitly rather than trusting the message sink.
    # TESTED 2026-09-21: sink(con, type='message') does capture message() under
    # Rscript, with or without an outer sink already in place -- so run_inline
    # was never broken, and the earlier note in CLAUDE.md saying it "swallows
    # message()" misread the symptom. But v6b's candidates_*.log and
    # landscape.log really do contain only cat()/print() output, which means
    # that run was invoked somewhere the message stream did not reach the sink
    # (the IRkernel path this file detects at `in_jupyter` is the likely one).
    # Those logs lost every "peaks in / masked / remaining / N loci /
    # consensus IQR" line and the counts had to be rebuilt from the CSVs.
    # Handling the condition directly makes the log correct however the driver
    # is invoked; muffleMessage stops the sink from writing it a second time.
    st = tryCatch(
             withCallingHandlers({ fn(); 0L },
                 message = function(m) {
                     cat(conditionMessage(m), file = con, sep = '')
                     invokeRestart('muffleMessage')
                 }),
             error = function(e) {
                 cat('ERROR: ', conditionMessage(e), '\n', file = con, sep = ''); 1L })
    unwind()
    close(con)
    cat('  exit ', st, '  (', round(as.numeric(difftime(Sys.time(), t0,
        units = 'secs'))), 's)  log: ', logf, '\n', sep = '')
    if (!identical(st, 0L)) {
        cat('!! ', label, ' FAILED — its outputs are absent or stale.\n',
            '   Last lines of the log:\n', sep = '')
        cat(paste0('   ', tail(readLines(logf, warn = FALSE), 15)), sep = '\n')
        cat('\n')
    }
    st
}

# The run tree still records the code that produced it; the unit is now
# the driver itself rather than the child scripts it used to copy in.
if (!is.na(driver_file) && file.exists(driver_file))
    file.copy(driver_file, script_snap, overwrite = TRUE)


# ════════════════════════════════════════════════════════════════════
# qtl_landscape.R — genome-wide, per-condition views of the peak tables
# ════════════════════════════════════════════════════════════════════
qtl_landscape_run = function(PEAK_DIR, OUT_DIR) {

## ---- configuration, replaces qtl_landscape.R lines 1-38 ------------------
FIT_CSV    <- file.path(PEAK_DIR, "fit_peaks.csv")
STRESS_CSV <- file.path(PEAK_DIR, "stress_peaks.csv")
.miss <- c(FIT_CSV, STRESS_CSV)[!file.exists(c(FIT_CSV, STRESS_CSV))]
if (length(.miss))
  stop("peak table(s) absent — run the peaks stage first:\n  ",
       paste(.miss, collapse = "\n  "), call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("peak tables : ", PEAK_DIR)
message("figures  -> : ", OUT_DIR)

## ---------------------------------------------------------------------------
## ---- shared setup for QTL landscape prototypes -----------------------------
suppressMessages({
  library(data.table); library(ggplot2); library(GenomicRanges)
  library(ggrepel); library(scales); library(patchwork)
})

CHRLEN <- c(chrI=230218, chrII=813184, chrIII=316620, chrIV=1531933, chrV=576874,
            chrVI=270161, chrVII=1090940, chrVIII=562643, chrIX=439888, chrX=745751,
            chrXI=666816, chrXII=1078177, chrXIII=924431, chrXIV=784333,
            chrXV=1091291, chrXVI=948066)
CHROMS <- names(CHRLEN)
GAP    <- 120000                      # visual gap between chromosomes, bp
OFFSET <- cumsum(c(0, head(CHRLEN + GAP, -1))); names(OFFSET) <- CHROMS
CHRMID <- OFFSET + CHRLEN/2
GENOME_END <- unname(OFFSET[length(OFFSET)] + CHRLEN[length(CHRLEN)])

## diverging poles: BBQ-favoured (blue) vs wild-parent-favoured (red), gray midpoint
COL_BBQ  <- "#2a78d6"; COL_WILD <- "#e34948"; COL_MID <- "#f0efec"
COL_INK  <- "#0b0b0b"; COL_INK2 <- "#52514e"; COL_GRID <- "#e6e5e1"

gx <- function(chrom, pos) unname(OFFSET[chrom]) + pos   # genome coordinate

read_peaks <- function(path) {
  p <- fread(path)
  p <- p[seqnames %in% CHROMS]
  if ("databin" %in% names(p)) p <- p[databin == "main"]   # V5 output; absent pre-V5
  p[, seqnames := factor(seqnames, levels = CHROMS)]
  p[, `:=`(gstart = gx(as.character(seqnames), start),
           gend   = gx(as.character(seqnames), end),
           gpeak  = gx(as.character(seqnames), peak.pos))]
  p[, favoured := fifelse(direction < 0, "BBQ", "wild parent")]
  p[]
}

## binned interval depth, split by direction; returns fraction of samples in condition
bin_depth <- function(p, binsize = 2000) {
  out <- list()
  for (cd in unique(p$condition)) for (dirn in c(-1, 1)) {
    x <- p[condition == cd & direction == dirn]
    if (!nrow(x)) next
    n_samp <- uniqueN(p[condition == cd]$sample)
    gr <- GRanges(factor(as.character(x$seqnames), levels = CHROMS),
                  IRanges(pmax(x$start, 1), x$end),
                  seqinfo = Seqinfo(CHROMS, unname(CHRLEN)))
    cv <- coverage(gr)
    for (ch in names(cv)) {
      r  <- cv[[ch]]; L <- CHRLEN[[ch]]
      bs <- seq(1, L, by = binsize); be <- pmin(bs + binsize - 1, L)
      v  <- viewMeans(Views(r, start = bs, end = pmin(be, length(r))))
      out[[length(out) + 1]] <- data.table(condition = cd, direction = dirn,
        chrom = ch, bin_mid = (bs + be)/2, depth = as.numeric(v),
        frac = as.numeric(v)/n_samp)
    }
  }
  d <- rbindlist(out)
  d[, gmid := gx(chrom, bin_mid)][]
}

## chromosome backdrop: alternating bands + labels
chrom_bands <- function(ymin, ymax) {
  bands <- data.table(chrom = CHROMS, x1 = OFFSET, x2 = OFFSET + CHRLEN,
                      i = seq_along(CHROMS))[i %% 2 == 0]
  geom_rect(data = bands, aes(xmin = x1, xmax = x2, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#000000", alpha = 0.035)
}
chrom_axis <- function() {
  scale_x_continuous(breaks = CHRMID, labels = sub("^chr", "", CHROMS),
                     expand = expansion(mult = 0.004), limits = c(0, GENOME_END))
}
theme_gw <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_line(colour = COL_GRID, linewidth = 0.3),
          axis.text.x  = element_text(size = base - 2, colour = COL_INK2),
          axis.title   = element_text(colour = COL_INK2),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", colour = COL_INK),
          strip.text.x = element_text(hjust = 0, face = "bold", colour = COL_INK),
          plot.title    = element_text(face = "bold", colour = COL_INK),
          plot.subtitle = element_text(colour = COL_INK2, size = base - 1),
          plot.caption  = element_text(colour = COL_INK2, size = base - 3, hjust = 0),
          legend.position = "top", legend.justification = "left",
          legend.title = element_text(colour = COL_INK2, size = base - 2),
          legend.text  = element_text(size = base - 2))
}

## beta is NEGATIVE when the BBQ allele is favoured, so low = BBQ pole
scale_beta <- function(aes_type = c("colour", "fill")) {
  f <- if (match.arg(aes_type) == "colour") scale_colour_gradient2 else scale_fill_gradient2
  f(low = COL_BBQ, mid = COL_MID, high = COL_WILD, midpoint = 0,
    limits = c(-0.5, 0.5), name = "allele-frequency shift at peak",
    breaks = c(-0.5, 0, 0.5), labels = c("BBQ 0.5", "0", "wild 0.5"),
    guide = guide_colourbar(barwidth = 9, barheight = 0.4, title.vjust = 1))
}

## ---------------------------------------------------------------------------
## ---- previously published xQTL ---------------------------------------------
## Treusch et al. 2015, PLoS Genet 11:e1004913, Supporting Table S4 "Grouped RR QTL"
## (pgen.1004913.s008.xlsx): the 2-LOD confidence intervals for all 60 grouped QTL
## from their round-robin X-QTL across 12 crosses - 23 caffeine, 37 NaCl.  These are
## real published intervals, not gene spans: median width 45.7 kb (range 7.4-99),
## against our own median LOD-drop width of 29.3 kb in v6b, so the two are directly
## comparable.  Coordinates are S288C R64; every named candidate gene falls inside
## its own interval except the ENA locus (their bound stops 3.3 kb short of ENA1).
## `maxpos` is their max-LOD position and is what the enrichment test below uses.
COL_PUB <- "#7a6a55"        # warm brown: outside the blue/red/purple data palette

TREUSCH <- fread(sep = "|", text = "
qtl|cond|chrN|lo|maxpos|hi|lod|ncross|crosses|qtg
Caff.chr03.1|Caf|3|267550|288775|291500|8.22|1|8|
Caff.chr04.1|Caf|4|442650|463300|502800|6.05|1|1|
Caff.chr04.2|Caf|4|866200|877925|896475|38.12|2|2 3|MSS4
Caff.chr07.1|Caf|7|468500|513400|567450|8.04|2|2 3|MSB2,MTL1,GSC2
Caff.chr08.1|Caf|8|18175|25675|41725|14.76|1|5|GUT1,VMR1
Caff.chr08.2|Caf|8|133700|178850|207650|8.53|1|10|SLT2
Caff.chr09.1|Caf|9|397400|412900|413000|5.58|1|5|
Caff.chr10.1|Caf|10|547725|561025|569800|29.08|8|1 2 3 5 6 8 9 10|TOR1 (confirmed)
Caff.chr10.2|Caf|10|431300|459000|477200|28.77|2|5 6|CYR1
Caff.chr11.1|Caf|11|46400|59650|73900|10.80|5|1 2 3 6 11|TOR2
Caff.chr11.2|Caf|11|613500|633650|634000|7.23|2|4 5|
Caff.chr12.1|Caf|12|614300|653850|682850|9.44|1|9|HAP1,GSY2
Caff.chr12.2|Caf|12|446550|460250|484850|6.89|1|10|
Caff.chr13.1|Caf|13|14900|26975|51250|10.07|1|10|RSC9,ATR1,GTR1
Caff.chr14.1|Caf|14|207250|227800|276150|5.64|1|1|
Caff.chr14.2|Caf|14|453125|468650|482700|19.62|1|10|
Caff.chr15.1|Caf|15|604100|625275|648050|8.61|4|2 4 9 10|PDR5
Caff.chr15.2|Caf|15|40650|45200|49150|8.41|1|3|
Caff.chr15.3|Caf|15|395000|406200|421500|38.55|3|5 7 8|WHI2 (confirmed)
Caff.chr15.4|Caf|15|115125|136725|163375|8.94|1|6|WSC3,PKH2
Caff.chr15.5|Caf|15|978600|1005250|1033900|7.43|1|10|
Caff.chr16.1|Caf|16|389900|416875|435400|21.14|2|6 8|BRO1,RGL1
Caff.chr16.2|Caf|16|187200|208000|236800|6.84|1|8|TCO89
NaCl.chr02.1|NaCl|2|510450|538250|567050|11.84|3|10 11 12|IRA1
NaCl.chr02.2|NaCl|2|189950|238950|269050|6.46|1|11|
NaCl.chr03.1|NaCl|3|14900|33600|54750|5.62|1|7|
NaCl.chr03.2|NaCl|3|263700|289800|291500|6.55|1|8|
NaCl.chr04.1|NaCl|4|491650|509600|531850|41.24|11|1 2 3 5 6 7 8 9 10 11 12|ENA locus
NaCl.chr04.2|NaCl|4|411600|421750|429050|18.21|1|3|GPD1
NaCl.chr04.3|NaCl|4|142400|192600|230050|6.05|1|7|
NaCl.chr04.4|NaCl|4|797100|834700|861150|8.99|1|7|
NaCl.chr04.5|NaCl|4|1368600|1397800|1412700|6.55|1|7|
NaCl.chr05.1|NaCl|5|81750|107825|133825|15.35|1|4|
NaCl.chr05.2|NaCl|5|493975|514325|546450|8.21|2|4 7|BCK2
NaCl.chr05.3|NaCl|5|387600|418800|441000|9.87|1|12|SHO1
NaCl.chr07.1|NaCl|7|898400|932050|958400|10.66|3|4 5 7|
NaCl.chr07.2|NaCl|7|511300|536050|571500|8.20|3|4 6 11|MSB2,MTL1,GSC2
NaCl.chr07.3|NaCl|7|14900|15000|22300|6.20|1|7|
NaCl.chr08.1|NaCl|8|27250|51250|65125|13.02|3|4 5 10|GUT1,WSC4,VMR1
NaCl.chr10.1|NaCl|10|530050|548300|561300|8.93|2|3 4|OSM1,STK2
NaCl.chr10.2|NaCl|10|420250|434650|451625|28.50|2|5 6|CYR1
NaCl.chr10.3|NaCl|10|182800|204725|234375|21.11|2|6 10|YJL107c & YJL108c
NaCl.chr10.4|NaCl|10|78050|104300|150900|6.08|1|9|
NaCl.chr10.5|NaCl|10|14900|37400|53550|6.51|1|11|
NaCl.chr11.1|NaCl|11|49425|66850|78250|12.84|3|8 10 11|TOR2
NaCl.chr11.2|NaCl|11|542750|598050|639900|5.46|1|9|
NaCl.chr11.3|NaCl|11|206100|237750|262650|7.24|1|10|
NaCl.chr12.1|NaCl|12|854850|872550|891450|10.71|1|7|
NaCl.chr12.2|NaCl|12|428400|457350|475850|13.79|1|12|
NaCl.chr13.1|NaCl|13|255700|288300|338200|8.09|3|3 4 12|
NaCl.chr13.2|NaCl|13|532600|571100|615800|6.35|1|4|
NaCl.chr13.3|NaCl|13|671575|706900|720300|15.45|2|7 11|EFR3
NaCl.chr14.1|NaCl|14|550350|582200|605000|8.01|3|4 8 12|CRZ1
NaCl.chr14.2|NaCl|14|132300|139400|145900|6.28|1|11|
NaCl.chr14.3|NaCl|14|225900|263050|299900|10.61|1|12|
NaCl.chr15.1|NaCl|15|57650|73800|87150|21.98|3|6 11 12|YGK3
NaCl.chr15.2|NaCl|15|398850|408100|412700|57.10|3|7 8 12|WHI2 (confirmed)
NaCl.chr15.3|NaCl|15|992200|1020100|1042025|9.66|2|8 10|PDE2
NaCl.chr16.1|NaCl|16|399150|434250|481050|6.80|3|1 7 8|
NaCl.chr16.2|NaCl|16|200625|206950|212175|43.89|2|7 8|TCO89,RTT10
")
TREUSCH[, `:=`(study = "Treusch 2015", coord = "2-LOD interval",
               chrom = CHROMS[chrN], crosses = gsub(" ", ",", crosses))]

## Ehrenreich et al. 2010, Nature 464:1039 (BYxRM X-QTL, 17 compounds).  That paper
## tabulates no QTL coordinates, so these are R64 spans of the four genes it names -
## a weaker kind of coordinate, kept deliberately: NEITHER caffeine NOR NaCl was
## assayed there, so these four are the negative control that calibrates the null.
EHRENREICH <- data.table(
  qtl    = c("Ehr.RAD5", "Ehr.MKT1", "Ehr.CUP1", "Ehr.HAP1"),
  cond   = "",                       # none of our conditions
  chrom  = c("chrXII", "chrXIV", "chrVIII", "chrXII"),
  lo     = c(204991, 467131, 212535, 646415),
  hi     = c(208500, 469623, 214718, 650923),
  lod    = NA_real_, ncross = NA_integer_, crosses = "",
  qtg    = c("RAD5", "MKT1", "CUP1", "HAP1"),
  study  = "Ehrenreich 2010", coord = "gene span")
EHRENREICH[, maxpos := (lo + hi)/2]

PUB <- rbind(TREUSCH[, .(qtl, cond, chrom, lo, maxpos, hi, lod, ncross, crosses, qtg, study, coord)],
             EHRENREICH[, .(qtl, cond, chrom, lo, maxpos, hi, lod, ncross, crosses, qtg, study, coord)])
PUB[, `:=`(glo = gx(chrom, lo), ghi = gx(chrom, hi), gmax = gx(chrom, maxpos),
           lane = fifelse(cond == "Caf", 2, fifelse(cond == "NaCl", 1, 0)))]

PUB_NEAR <- 0   # bp of slack allowed around the published max-LOD position

## --- our peaks against theirs ----------------------------------------------
## Primary call: does our LOD-drop interval cover their max-LOD position?  That keeps
## the null exact (see pub_enrich) and does not reward us for their interval being
## wide.  Interval-vs-interval overlap is carried alongside as a descriptive.
pub_hits <- function(p, near = PUB_NEAR) {
  a <- copy(p)[, chrom := as.character(seqnames)]
  b <- PUB[, .(chrom, qtl, pub_cond = cond, qtg, study, coord, lod_pub = lod,
               ncross_pub = ncross, pub_lo = lo, pub_hi = hi, pub_max = maxpos)]
  h <- merge(a, b, by = "chrom", allow.cartesian = TRUE)
  h <- h[start - near <= pub_max & end + near >= pub_max]
  if (!nrow(h)) return(h)
  h[, `:=`(ov_bp      = pmax(pmin(end, pub_hi) - pmax(start, pub_lo) + 1, 0),
           peak_off_kb = round((peak.pos - pub_max)/1000, 1),
           cond_match  = nzchar(pub_cond) & as.character(condition) == pub_cond)]
  h[, `:=`(ov_frac_pub = round(ov_bp/(pub_hi - pub_lo + 1), 3),
           ov_frac_ours = round(ov_bp/width, 3))]
  h[order(qtl, condition, cross, well)][]
}

## collapse onto the peak table as flag columns
pub_flag <- function(p, near = PUB_NEAR) {
  h <- pub_hits(p, near)
  p[, `:=`(pub_qtl = NA_character_, pub_qtg = NA_character_, pub_study = NA_character_,
           pub_cond_match = FALSE, pub_n = 0L)]
  if (!nrow(h)) return(invisible(p))
  key <- c("seqnames", "start", "end", "condition", "cross", "well")
  g <- h[, .(pub_qtl        = paste(unique(qtl), collapse = ","),
             pub_qtg        = paste(unique(qtg[nzchar(qtg)]), collapse = ","),
             pub_study      = paste(sort(unique(study)), collapse = ","),
             pub_cond_match = any(cond_match),
             pub_n          = .N), by = key]
  p[g, on = key, `:=`(pub_qtl = i.pub_qtl, pub_qtg = i.pub_qtg, pub_study = i.pub_study,
                      pub_cond_match = i.pub_cond_match, pub_n = i.pub_n)]
  invisible(p)
}

## Is a hit more than our intervals being wide?  For each pool, the chance that a
## randomly placed set of intervals of the same total length covers a given point is
## that pool's covered fraction of the genome.  Expected pools = sum of those
## fractions; observed is compared against it as a Poisson (the exact null is a
## Poisson-binomial, which this approximates well at these rates).  Pools of one
## cross are not independent, so read the q values as a ranking.
pub_cov <- function(p, samp_col = "sample") {
  GLEN <- sum(CHRLEN)
  cv <- p[, .(cov = {
    gr <- GRanges(factor(as.character(seqnames), levels = CHROMS),
                  IRanges(pmax(start, 1), end), seqinfo = Seqinfo(CHROMS, unname(CHRLEN)))
    sum(width(reduce(gr)))/GLEN
  }), by = c(samp_col, "condition")]
  setnames(cv, samp_col, "samp")[]
}

pub_enrich <- function(p, samp_col = "sample") {
  cv <- pub_cov(p, samp_col)
  pp <- copy(p)[, samp := get(samp_col)][, chrom := as.character(seqnames)]
  res <- rbindlist(lapply(seq_len(nrow(PUB)), function(i) {
    L <- PUB[i]
    covering <- unique(pp[chrom == L$chrom & start <= L$maxpos & end >= L$maxpos, samp])
    sets <- list(all = cv)
    if (nzchar(L$cond) && L$cond %in% as.character(cv$condition))
      sets[["matched"]] <- cv[condition == L$cond]
    rbindlist(lapply(names(sets), function(nm) {
      s <- sets[[nm]]; obs <- sum(covering %in% s$samp); exp <- sum(s$cov)
      data.table(qtl = L$qtl, qtg = L$qtg, study = L$study, coord = L$coord,
                 pub_cond = L$cond, lod_pub = L$lod, ncross_pub = L$ncross,
                 pool_set = nm, pools_tested = nrow(s), obs_pools = obs,
                 exp_pools = round(exp, 2), fold = round(obs/exp, 2),
                 p_pois = ppois(obs - 1, exp, lower.tail = FALSE))
    }))
  }))
  res[, q_bh := p.adjust(p_pois, "BH"), by = pool_set]
  res[order(pool_set, p_pois)][]
}

## How much of their catalogue do we recover, and how much would we recover by chance?
pub_recovery <- function(p, samp_col = "sample", kmin = 1:3) {
  cv <- pub_cov(p, samp_col)
  pp <- copy(p)[, samp := get(samp_col)][, chrom := as.character(seqnames)]
  T4 <- PUB[study == "Treusch 2015"]
  rbindlist(lapply(unique(T4$cond), function(cd) {
    s  <- cv[condition == cd]; qs <- T4[cond == cd]
    obs <- vapply(seq_len(nrow(qs)), function(i)
      uniqueN(pp[samp %in% s$samp & chrom == qs$chrom[i] &
                 start <= qs$maxpos[i] & end >= qs$maxpos[i], samp]), integer(1))
    lam <- sum(s$cov)            # expected pools covering ANY given point
    rbindlist(lapply(kmin, function(k) {
      pk <- ppois(k - 1, lam, lower.tail = FALSE)   # P(a random point is hit by >= k pools)
      data.table(condition = cd, our_pools = nrow(s), their_qtl = nrow(qs), min_pools = k,
                 recovered = sum(obs >= k), exp_recovered = round(nrow(qs) * pk, 1),
                 fold = round(sum(obs >= k)/(nrow(qs) * pk), 2),
                 p_binom = binom.test(sum(obs >= k), nrow(qs), pk, "greater")$p.value)
    }))
  }))
}

## ---- regions whose recurrence is already explained -------------------------
## Bounds are the 50%-of-maximum plateau of the pooled pileup.  Two of these come
## out of every analysis; the third comes out of the concordance test only.
##   chrIII 125-215 kb  MAT (200.4 kb, 51% of pools) - mating-type selection,
##                      reinforced by the drugs; sweeps the whole right arm.
##                      Artefact of the experiment.  Always masked.
##   chrIV   15-60  kb  HO (47.1 kb, 41% of pools) - selected by the drugs.
##                      Artefact of the experiment.  Always masked.
##   chrVII 440-600 kb  the BBQ large-effect growth QTL: BBQ-favoured in 184 of
##                      185 peaks, 34 of 35 crosses, every condition, median
##                      |beta| 0.37.  This is a FINDING, not an artefact, so it is
##                      reported - MASK_REPORT keeps it, and it appears in the peak
##                      counts, in our_loci() and on the figures.
##                      It is still withheld from the published-overlap test
##                      (MASK_PUB drops it), because it sits on top of two of
##                      Treusch's loci - Caff.chr07.1 (468,500-567,450) and
##                      NaCl.chr07.2 (511,300-571,500) - and we cannot tell the two
##                      apart there.  Scoring them as recovered would credit our own
##                      shared-parent sweep as replication of theirs.
MASK_ALL <- data.table(
  what  = c("MAT", "HO", "BBQ growth QTL"),
  note  = c("mating-type selection", "drug selection",
            "our locus - reported, withheld from concordance test"),
  chrom = c("chrIII", "chrIV", "chrVII"),
  lo    = c(125000, 15000, 440000), hi = c(215000, 60000, 600000),
  scope = c("always", "always", "concordance test only"))
MASK_ALL[, `:=`(glo = gx(chrom, lo), ghi = gx(chrom, hi), gmid = gx(chrom, (lo + hi)/2))]

MASK_REPORT <- MASK_ALL[scope == "always"]   # peak counts, our_loci(), figures
MASK_PUB    <- MASK_ALL                      # published-overlap test only

## The summit test alone leaks: a peak called just outside a mask can carry a
## 3-LOD interval that covers the whole masked region, and the sweep comes through
## intact.  chrIII 90-110 kb is the worst case - 36 peaks, BBQ-favoured in 36 of 36,
## median interval 99.9-215.3 kb against a mask edge at 215 kb.  So a peak is also
## masked when its interval covers more than SHOULDER_FRAC of a masked region.
## Costs 81 peaks genome-wide (62 MAT, 19 chrVII), 80 of them BBQ-favoured.
SHOULDER_FRAC <- 0.5

in_mask <- function(chrom, pos, mask = MASK_REPORT, start = NULL, end = NULL) {
  g   <- gx(as.character(chrom), pos)
  hit <- Reduce(`|`, lapply(seq_len(nrow(mask)), function(i) g >= mask$glo[i] & g <= mask$ghi[i]))
  if (!is.null(start)) {                       # interval given: also drop shoulders
    gs <- gx(as.character(chrom), start); ge <- gx(as.character(chrom), end)
    hit <- hit | Reduce(`|`, lapply(seq_len(nrow(mask)), function(i)
      pmax(0, pmin(ge, mask$ghi[i]) - pmax(gs, mask$glo[i])) >
        SHOULDER_FRAC * (mask$ghi[i] - mask$glo[i])))
  }
  hit
}

## ---- cross-level replication test ------------------------------------------
## A QTL only segregates in a cross whose parents differ at it, and we share no wild
## parent with Treusch, so "how many of OUR POOLS hit it" is the wrong unit - it asks
## the locus to be polymorphic and large-effect in many of our crosses.  The unit is
## the CROSS: did any pool of cross C, in the matching condition, put an interval over
## their max-LOD position?  33 of their 60 grouped QTL were themselves seen in exactly
## one of their 12 crosses, so cross-private is the normal case on both sides.
##
## Null: circular-shift each pool's intervals around the genome as a block, preserving
## interval count, width and spacing.  Poisson placement would understate the null
## because our intervals are autocorrelated.
pub_crosslevel <- function(p, nperm = 1000, seed = 1) {
  GL  <- sum(CHRLEN)
  gof <- cumsum(c(0, head(CHRLEN, -1))); names(gof) <- CHROMS   # gapless genome
  pk  <- p[!in_mask(seqnames, peak.pos, MASK_PUB, start, end)]   # chrVII withheld here only
  pk[, `:=`(a = gof[as.character(seqnames)] + start, b = gof[as.character(seqnames)] + end)]
  qs  <- PUB[study == "Treusch 2015"][!in_mask(chrom, maxpos, MASK_PUB)]
  qs[, m := gof[chrom] + maxpos]

  count <- function(a, b, cond, cross) vapply(seq_len(nrow(qs)), function(i) {
    k <- cond == qs$cond[i] & ((a <= qs$m[i] & b >= qs$m[i]) |
                               (b > GL & a <= qs$m[i] + GL & b >= qs$m[i] + GL))
    uniqueN(cross[k])
  }, integer(1))

  obs <- count(pk$a, pk$b, as.character(pk$condition), pk$cross)
  det <- data.table(qtl = qs$qtl, qtg = qs$qtg, cond = qs$cond, their_lod = qs$lod,
                    their_crosses = qs$ncross, our_crosses = obs)
  det[, our_pools := vapply(seq_len(.N), function(i)
    uniqueN(pk[condition == cond[i] & a <= qs$m[i] & b >= qs$m[i], sample]), integer(1))]

  set.seed(seed)
  idx <- split(seq_len(nrow(pk)), pk$sample)
  perm <- vapply(seq_len(nperm), function(z) {
    a <- pk$a; b <- pk$b
    for (k in idx) { sh <- runif(1, 0, GL); w <- b[k] - a[k]
                     a[k] <- (a[k] + sh) %% GL; b[k] <- a[k] + w }
    cc <- count(a, b, as.character(pk$condition), pk$cross)
    c(sum(cc), sum(cc > 0))
  }, integer(2))

  tot <- sum(obs); anyk <- sum(obs > 0)
  pt  <- perm[1, ]; pa <- perm[2, ]
  ct  <- suppressWarnings(cor.test(det$their_crosses, det$our_crosses, method = "spearman"))
  cl  <- suppressWarnings(cor.test(det$their_lod,     det$our_crosses, method = "spearman"))
  ## permute the pairing: confirms the dose-response is not an artefact of the marginals
  set.seed(seed + 1)
  pr  <- replicate(5000, suppressWarnings(
           cor(sample(det$their_crosses), det$our_crosses, method = "spearman")))
  list(detections = det[order(-our_crosses)],
       summary = data.table(
         loci_tested = nrow(qs), peaks_used = nrow(pk), peaks_masked = nrow(p) - nrow(pk),
         obs_detections = tot, null_mean = round(mean(pt), 1),
         null_lo = quantile(pt, .025), null_hi = quantile(pt, .975),
         fold = round(tot/mean(pt), 2), p_perm = (1 + sum(pt >= tot))/(nperm + 1),
         obs_any = anyk, null_any_mean = round(mean(pa), 1),
         p_any = (1 + sum(pa >= anyk))/(nperm + 1),
         rho_their_crosses = round(ct$estimate, 3), p_rho = signif(ct$p.value, 3),
         p_rho_labelperm = round((1 + sum(pr >= ct$estimate))/5001, 4),
         rho_their_lod = round(cl$estimate, 3), p_rho_lod = signif(cl$p.value, 3)),
       dose = det[, .(loci = .N, our_crosses_mean = round(mean(our_crosses), 2),
                      our_pools_mean = round(mean(our_pools), 2)),
                  by = .(their = fifelse(their_crosses == 1, "1 cross",
                         fifelse(their_crosses == 2, "2 crosses", "3+ crosses")))][order(their)],
       null = pt)
}

## a guide rule only for the QTL they saw in 3 or more crosses - 60 rules is noise
## `which_cond` restricts the overlay to the loci published for that condition.
## Treusch et al. assayed only caffeine and NaCl, so a glucose or acetic acid
## panel gets no published overlay at all rather than a mismatched one.
## V7: ported from the standalone qtl_landscape.R, where this fix was made on
## 2026-09-18 -- one day AFTER the V6.1 inline, so it never reached the driver
## and v6b's protoF_surface_*.png carry the old mismatched overlays.
pub_rules <- function(which_cond = NULL) {
  d <- PUB[!is.na(ncross) & ncross >= 3]
  if (!is.null(which_cond)) d <- d[cond == which_cond]
  geom_vline(data = d, aes(xintercept = gmax), inherit.aes = FALSE,
             colour = COL_PUB, linewidth = 0.25, linetype = "22", alpha = 0.7)
}

## --- annotation track -------------------------------------------------------
## Their intervals as two lanes of bars on the same x scale, so ours can be read
## against theirs directly.  Filled = we cover their max-LOD position with at least
## two pools of the matching condition; hollow = we do not.
PUB_LANES <- c("NaCl", "Caf")
pub_track <- function(p = NULL, samp_col = "sample", size = 2.4, min_pools = 2,
                      which_cond = NULL) {
  d <- PUB[study == "Treusch 2015"]
  if (!is.null(which_cond)) d <- d[cond == which_cond]
  d[, hit := FALSE]
  if (!is.null(p)) {
    pp <- copy(p)[, samp := get(samp_col)][, chrom := as.character(seqnames)]
    d[, hit := vapply(seq_len(.N), function(i)
      uniqueN(pp[as.character(condition) == cond[i] & chrom == chrom[i] &
                 start <= maxpos[i] & end >= maxpos[i], samp]) >= min_pools, logical(1))]
  }
  lanes <- if (is.null(which_cond)) PUB_LANES else intersect(PUB_LANES, which_cond)
  d[, y := match(cond, lanes)]
  mk <- copy(MASK_REPORT)[, y := 0.15]
  ## label only their multi-cross loci, first gene name, one label per gene
  lab <- d[nzchar(qtg) & ncross >= 3][
    , .(gmax, y, lod, g = sub("[, ].*$", "", sub(" [(]confirmed[)]", "*", qtg)))][
    order(-lod)][, .SD[1], by = g]
  ggplot(d) +
    chrom_bands(-0.25, 2.55) +
    geom_rect(data = mk, aes(xmin = glo, xmax = ghi, ymin = -0.2, ymax = 0.1),
              fill = "#c9c3b6", colour = NA, inherit.aes = FALSE) +
    geom_text(data = mk, aes(gmid, -0.05, label = what), inherit.aes = FALSE,
              size = 2.1, colour = "#5d564a", vjust = 0.5, hjust = -0.08) +
    geom_rect(aes(xmin = glo, xmax = ghi, ymin = y - 0.3, ymax = y + 0.3,
                  fill = hit), colour = COL_PUB, linewidth = 0.25) +
    geom_segment(aes(gmax, y - 0.3, xend = gmax, yend = y + 0.3),
                 colour = COL_PUB, linewidth = 0.4) +
    ggrepel::geom_text_repel(data = lab, aes(gmax, 2.35, label = g), size = size,
      colour = COL_PUB, fontface = "italic", direction = "x", vjust = 0,
      ylim = c(2.62, 3.1), segment.size = 0.2, segment.colour = COL_PUB,
      min.segment.length = 0, box.padding = 0.15, max.overlaps = Inf, seed = 3) +
    scale_fill_manual(values = c(`TRUE` = COL_PUB, `FALSE` = "#ffffff"), guide = "none") +
    scale_y_continuous(breaks = c(-0.05, seq_along(PUB_LANES)), labels = c("explained", PUB_LANES),
                       limits = c(-0.28, 3.15), expand = expansion(mult = 0)) +
    chrom_axis() + labs(y = NULL) +
    theme_gw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          panel.grid.major.y = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 7, colour = COL_PUB),
          plot.margin = margin(0, 0, 2, 0))
}
PUB_CAPTION <- paste(
  "Track: the 60 grouped QTL of Treusch 2015 (S4 Table) as their published 2-LOD intervals, caffeine and NaCl lanes.",
  "\nTick = their max-LOD position; filled = we cover it with 2+ pools of the matching condition.",
  "Labelled: their candidate gene where they mapped the QTL in 3 or more crosses (* = allele-replacement confirmed).")

OUT <- OUT_DIR
fit <- read_peaks(FIT_CSV)
fit[, condition := factor(condition, levels = c("Glu", "Ace", "Caf", "NaCl"))]
nsamp <- fit[, .(n = uniqueN(sample)), by = condition]
labf  <- setNames(sprintf("%s\n(n=%d)", nsamp$condition, nsamp$n), nsamp$condition)

## ---- flag our peaks against the published catalogue ------------------------
pub_flag(fit)
hits <- pub_hits(fit)
fwrite(PUB, file.path(OUT, "published_loci_used.csv"))
fwrite(hits[, .(qtl, qtg, study, coord, pub_cond, lod_pub, ncross_pub,
                chrom, pub_lo, pub_max, pub_hi,
                condition, cross, well, sample, start, end, width, peak.pos,
                peak.LOD, peak.beta, direction, favoured,
                cond_match, peak_off_kb, ov_bp, ov_frac_pub, ov_frac_ours)],
       file.path(OUT, "published_overlap_peaks.csv"))
fwrite(fit, file.path(OUT, "fit_peaks_pubflagged.csv"))

pub_en  <- pub_enrich(fit);   fwrite(pub_en,  file.path(OUT, "published_overlap_enrichment.csv"))
pub_rec <- pub_recovery(fit); fwrite(pub_rec, file.path(OUT, "published_recovery.csv"))

message(sprintf("published QTL touched: %d of %d; %d of %d of our peaks flagged (%d in the matching condition)",
                uniqueN(hits$qtl), nrow(PUB), sum(!is.na(fit$pub_qtl)), nrow(fit),
                sum(fit$pub_cond_match)))
cat("\n-- how much of Treusch S4 we recover, vs the same-coverage null --\n"); print(pub_rec)
cat("\n-- per-QTL enrichment, pools of the matching condition --\n")
print(pub_en[pool_set == "matched"][order(p_pois)][1:15])
cat("\n-- Ehrenreich negative control (conditions never assayed here) --\n")
print(pub_en[study == "Ehrenreich 2010"])

## the test that matters: cross-level, with the three explained regions removed
cl <- pub_crosslevel(fit)
fwrite(cl$detections, file.path(OUT, "published_crosslevel_detections.csv"))
fwrite(cl$summary,    file.path(OUT, "published_crosslevel_summary.csv"))
fwrite(MASK_ALL,      file.path(OUT, "masked_regions.csv"))   # `scope` says which test drops which
cat("\n== cross-level replication, MAT / HO / chrVII removed (chrVII for this test only) ==\n")
print(cl$summary)
cat("\ndose-response - their cross count vs ours (random intervals cannot produce this):\n")
print(cl$dose)
cat("\ntop loci:\n"); print(cl$detections[1:12])

## ---- how many loci we find, counted the way Treusch counted theirs ---------
## Their rule: within a condition, merge overlapping 2-LOD intervals across crosses;
## each merged block is one grouped QTL.  Same rule applied to our intervals, with
## the two artefact regions removed (MASK_REPORT); the chrVII locus is kept.
our_loci <- function(p) {
  q <- p[!in_mask(seqnames, peak.pos, start = start, end = end)]
  tq <- PUB[study == "Treusch 2015"]
  rbindlist(lapply(levels(droplevels(q$condition)), function(cd) {
    d  <- q[condition == cd]
    gr <- reduce(GRanges(as.character(d$seqnames),
                         IRanges(d$start, pmin(d$end, CHRLEN[as.character(d$seqnames)]))))
    xl <- d[, .(n = length(reduce(GRanges(as.character(seqnames),
              IRanges(start, pmin(end, CHRLEN[as.character(seqnames)])))))), by = cross]
    tc <- if (cd %in% tq$cond) tq[cond == cd] else NULL
    ov <- if (is.null(tc)) rep(0L, length(gr)) else
            countOverlaps(gr, GRanges(tc$chrom, IRanges(tc$lo, tc$hi)))
    data.table(condition = cd, pools = uniqueN(d$sample), crosses = uniqueN(d$cross),
               peaks = nrow(d), grouped_loci = length(gr),
               median_kb = round(median(width(gr))/1000, 1),
               cross_level_detections = sum(xl$n),
               loci_overlapping_Treusch = sum(ov > 0),
               loci_not_in_Treusch = if (is.null(tc)) NA_integer_ else sum(ov == 0),
               their_grouped_loci = if (is.null(tc)) NA_integer_ else nrow(tc),
               their_cross_level = if (is.null(tc)) NA_integer_ else sum(tc$ncross))
  }))
}
ol <- our_loci(fit); fwrite(ol, file.path(OUT, "our_loci_summary.csv"))
cat("\n== our loci, grouped by Treusch's own rule (overlap-merge within a condition) ==\n")
print(ol)

## ═══════════════════════════════════════════════════════════════════════════
## V7: the same published tests on BOTH peak sets and BOTH units of detection
## ═══════════════════════════════════════════════════════════════════════════
## WHY. Treusch et al. mapped with MULTIPOOL in `-m contrast` mode -- each
## selected pool against its matched YPD control (their Methods; our peaks
## against 0.5 have no equivalent in that paper).  So the set of ours that is
## comparable to theirs is the STRESS contrast, treated minus its paired glucose
## pool, not the fitness peaks.  A fitness peak carries every locus that also
## acts in glucose, which is precisely what their control subtraction removes.
##
## Everything above this block is UNCHANGED and still runs on `fit`, so a v7 run
## reproduces v6b's tables exactly and the two answers can be read side by side.
## Nothing here overwrites anything there -- these files all end in `_bysets`.
##
## UNIT.  `sample` counts pools, `cross` counts distinct crosses.  Pools of one
## cross share parents and segregant material, so three replicate wells of one
## cross agreeing is not the evidence that three crosses agreeing is.  Both
## pub_enrich() and pub_recovery() already take samp_col, and pub_cov() then
## builds the null from each cross's own union coverage, which is the right
## denominator for a by-cross test.
##
## TWO CONSEQUENCES, to carry into the prose as facts rather than as bugs:
##   - glucose leaves the published comparison entirely.  It is the control on
##     both sides, so there is no contrast to test and no Treusch set to test it
##     against.  Their only YPD QTL are the six MAT-dependent ones in S8 Table.
##   - acetate has no Treusch counterpart at all; they never assayed it.
## Only caffeine and NaCl are testable on contrast peaks.  That is a narrowing
## of what the section can claim, and it is the correct narrowing.

PSETS <- list(
  fitness  = fit,
  contrast = {
    .s <- read_peaks(STRESS_CSV)
    ## stress_peaks.csv carries L/R, not a single pool id; same construction the
    ## protoE block below uses, so the two agree on what a contrast pool is.
    .s[, sample := paste(cross, well, condition, sep = "_")]
    .s[, condition := factor(condition, levels = c("Ace", "Caf", "NaCl"))]
    ## RENAME L/R. pub_enrich() does `L <- PUB[i]` and then filters with
    ## `pp[chrom == L$chrom & ...]`; inside [.data.table the table's OWN column
    ## wins over the local variable, so a column called `L` turns L$chrom into
    ## "$ operator is invalid for atomic vectors" and the whole contrast peak
    ## set silently produced zero enrichment rows on 2026-09-21. The fitness
    ## table has no such column, which is why only this set broke. Renaming
    ## here rather than editing pub_enrich keeps the fitness path byte-identical
    ## to v6b, which is the acceptance test.
    setnames(.s, c("L", "R"), c("pool_L", "pool_R"), skip_absent = TRUE)
    .s
  })

## THE UNIT COLUMN MUST BE CONDITION-UNIQUE.  pub_cov() groups by
## (samp_col, condition) and pub_enrich()/pub_recovery() then select the
## condition-matched units with `samp %in% s$samp`.  That filter enforces
## condition matching only because a `sample` id is unique to one
## (cross, well, condition).  A bare `cross` is NOT: the same cross appears in
## all four conditions, so passing samp_col = "cross" lets a locus detected in
## GLUCOSE count as a caffeine detection.  Caught 2026-09-21 by the fact that it
## recovered 21 of 23 caffeine loci at k >= 1 where the pool unit recovered 13 --
## impossible, since a locus covered by one pool is covered by one cross.
## `cross_cond` keeps the intent (replicate wells of one cross are one unit,
## not three) while staying inside the condition.  At k >= 1 it must agree with
## the pool unit exactly; that equality is the regression test for this bug.
for (.nm in names(PSETS)) PSETS[[.nm]][, cross_cond := paste(cross, condition, sep = "|")]
UNITS <- c(pool = "sample", cross = "cross_cond")

cat("\n\n══════ V7: published tests by peak set and unit ══════\n")
for (ps in names(PSETS))
  cat(sprintf("  %-8s %5d peaks  %3d pools  %2d crosses  (%s)\n", ps,
              nrow(PSETS[[ps]]), uniqueN(PSETS[[ps]]$sample),
              uniqueN(PSETS[[ps]]$cross),
              paste(levels(droplevels(PSETS[[ps]]$condition)), collapse = "/")))

## ---- per-locus enrichment, and recovery against a same-coverage null -------
grid_over <- function(f, rename = NULL) rbindlist(lapply(names(PSETS), function(ps)
  rbindlist(lapply(names(UNITS), function(un) {
    r <- tryCatch(f(PSETS[[ps]], samp_col = UNITS[[un]]), error = function(e) {
      message("  !! ", ps, "/", un, ": ", conditionMessage(e)); NULL })
    if (is.null(r) || !nrow(r)) return(NULL)
    if (!is.null(rename)) setnames(r, names(rename), unname(rename), skip_absent = TRUE)
    r[, `:=`(peak_set = ps, unit = un)][]
  }))), fill = TRUE)

pub_en_sets  <- grid_over(pub_enrich)
fwrite(pub_en_sets, file.path(OUT, "published_overlap_enrichment_bysets.csv"))

pub_rec_sets <- grid_over(pub_recovery,
                          rename = c(our_pools = "our_units", min_pools = "min_units"))
fwrite(pub_rec_sets, file.path(OUT, "published_recovery_bysets.csv"))

cat("\n-- recovery vs a same-coverage null, every peak set x unit --\n")
print(pub_rec_sets[order(peak_set, unit, condition, min_units)])

## Regression guard for the condition-leak described at UNITS above: at a
## threshold of one unit, "detected by >= 1 pool" and "detected by >= 1 cross"
## are the same set of loci.  If these ever disagree the unit column has stopped
## being condition-unique and every by-cross number in this block is inflated.
.k1 <- dcast(pub_rec_sets[min_units == 1], peak_set + condition ~ unit,
             value.var = "recovered")
if (nrow(.k1) && !isTRUE(all.equal(.k1$pool, .k1$cross))) {
  print(.k1)
  stop("by-pool and by-cross recovery disagree at k = 1: the unit column is ",
       "leaking across conditions.  See the UNITS comment above.", call. = FALSE)
}
cat("  [ok] pool and cross units agree at k = 1 (", paste(.k1$pool, collapse = ", "), ")\n")

cat("\n-- per-locus enrichment, condition-matched, q < 0.1 in any cell --\n")
print(pub_en_sets[pool_set == "matched" & q_bh < 0.1][order(peak_set, unit, p_pois)],
      nrows = 40)

## TOR1 and ENA by name.  TOR1 is the caffeine hallmark in Treusch (8 of 12
## crosses) and the draft currently carries a bare "[result]" where its number
## should be; ENA is the one locus the paper claims individually.
cat("\n-- the two named loci, every cell --\n")
print(pub_en_sets[grepl("TOR1|ENA", qtg) & pool_set == "matched",
                  .(peak_set, unit, qtl, qtg, pub_cond, units_tested = pools_tested,
                    obs = obs_pools, exp = exp_pools, fold, p_pois, q_bh)][order(qtg, peak_set, unit)])

## ---- cross-level replication, both peak sets ------------------------------
## The unit is already the cross inside pub_crosslevel(), so this varies the
## peak set only.  The permutation blocks by pool, which is why the synthesized
## `sample` column above matters.
cl_sets <- rbindlist(lapply(names(PSETS), function(ps) {
  z <- tryCatch(pub_crosslevel(PSETS[[ps]]), error = function(e) {
    message("  !! crosslevel/", ps, ": ", conditionMessage(e)); NULL })
  if (is.null(z)) return(NULL)
  fwrite(z$detections, file.path(OUT, sprintf("published_crosslevel_detections_%s.csv", ps)))
  z$summary[, peak_set := ps][]
}), fill = TRUE)
fwrite(cl_sets, file.path(OUT, "published_crosslevel_bysets.csv"))
cat("\n-- cross-level replication, fitness vs contrast peaks --\n")
print(cl_sets)
cat("\nNOTE: null_lo/null_hi are the 2.5%/97.5% QUANTILES of the permuted totals,",
    "\nnot the range.  p_perm is the add-one tail and is the quotable number.\n")

## ---- loci counted Treusch's way, both peak sets ---------------------------
ol_sets <- rbindlist(lapply(names(PSETS), function(ps)
  our_loci(PSETS[[ps]])[, peak_set := ps][]), fill = TRUE)
fwrite(ol_sets, file.path(OUT, "our_loci_summary_bysets.csv"))
cat("\n-- grouped loci, both peak sets --\n"); print(ol_sets)

## ---- how much genome our intervals actually cover -------------------------
## The draft quotes "the genome about 3.6 times over", computed by hand as
## 1,472 x the 29.3 kb median / 12.1 Mb.  A median times a count is not a total.
## This emits the real summed width, and the union, per peak set.
GLEN <- sum(CHRLEN)
cov_sets <- rbindlist(lapply(names(PSETS), function(ps) {
  p <- PSETS[[ps]][!in_mask(seqnames, peak.pos, start = start, end = end)]
  gr <- GRanges(factor(as.character(p$seqnames), levels = CHROMS),
                IRanges(pmax(p$start, 1), pmin(p$end, CHRLEN[as.character(p$seqnames)])),
                seqinfo = Seqinfo(CHROMS, unname(CHRLEN)))
  data.table(peak_set = ps, peaks = nrow(p),
             total_interval_Mb = round(sum(width(gr))/1e6, 1),
             union_Mb          = round(sum(width(reduce(gr)))/1e6, 2),
             genome_Mb         = round(GLEN/1e6, 1),
             times_over        = round(sum(width(gr))/GLEN, 2),
             union_frac        = round(sum(width(reduce(gr)))/GLEN, 3),
             median_kb         = round(median(width(gr))/1000, 1))
}))
fwrite(cov_sets, file.path(OUT, "interval_coverage_bysets.csv"))
cat("\n-- interval coverage (replaces the hand-computed 3.6x) --\n"); print(cov_sets)

## ---- does the shared parent explain the cross-level excess? ---------------
## Every cross shares BBQ, so repeated detection of a published locus may be one
## BBQ allele seen many times rather than independent replication.  The draft
## carries an unfilled "[N] of the 71 detections favored the BBQ allele".
bbq_share <- rbindlist(lapply(names(PSETS), function(ps) {
  p <- PSETS[[ps]][!in_mask(seqnames, peak.pos, MASK_PUB, start, end)]
  h <- pub_hits(p)
  if (!nrow(h)) return(NULL)
  h <- h[cond_match == TRUE & study == "Treusch 2015"]
  if (!nrow(h)) return(NULL)
  ## One row per (locus, cross) -- the unit the cross-level test counts.  A cross
  ## can put peaks of BOTH directions over one published position, so deduping on
  ## (qtl, cross, favoured) returns more rows than there are detections and the
  ## total will not reconcile with obs_detections.  Break the tie on the largest
  ## |shift|, which is the peak that would be reported for that locus anyway.
  setorder(h, qtl, cross, -ov_bp)
  d <- h[order(qtl, cross, -abs(peak.beta))][, .SD[1], by = .(qtl, cross)]
  data.table(peak_set = ps, detections = nrow(d),
             bbq_favoured = sum(d$favoured == "BBQ"),
             pct_bbq = round(100 * mean(d$favoured == "BBQ"), 1))
}))
fwrite(bbq_share, file.path(OUT, "published_crosslevel_direction.csv"))
cat("\n-- direction of the cross-level detections (fills the [N] in the draft) --\n")
print(bbq_share)
cat("\n══════ end V7 published-test grid ══════\n\n")

## ---------- A. signed recurrence profile -----------------------------------
d <- bin_depth(fit, binsize = 2000)
d[, condition := factor(condition, levels = levels(fit$condition))]
d[, signed := fifelse(direction < 0, 100 * frac, -100 * frac)]
ylo <- min(d$signed); yhi <- max(d$signed)
pad <- 0.04 * (yhi - ylo); ylo <- ylo - pad; yhi <- yhi + pad + 10  # headroom for labels
pA <- ggplot(d, aes(gmid, signed)) +
  chrom_bands(ylo, yhi) +
  pub_rules() +
  geom_area(aes(fill = factor(direction, c(-1, 1))), position = "identity") +
  geom_hline(yintercept = 0, colour = COL_INK2, linewidth = 0.3) +
  facet_grid(condition ~ ., labeller = labeller(condition = labf)) +
  scale_fill_manual(values = c(`-1` = COL_BBQ, `1` = COL_WILD),
                    labels = c("BBQ allele favoured", "wild-parent allele favoured"),
                    name = NULL) +
  scale_y_continuous(labels = function(x) paste0(abs(x), "%"),
                     limits = c(ylo, yhi), expand = expansion(mult = 0)) +
  chrom_axis() +
  labs(
       x = "chromosome", y = "% of pools") +
  theme_gw()
## local maxima of the BBQ arm, >=25% of pools, >=50 kb apart
tops <- d[direction == -1 & frac >= 0.25][order(condition, -frac)]
keep <- tops[, {
  sel <- integer(0)
  for (i in seq_len(.N)) if (!any(abs(gmid[i] - gmid[sel]) < 50000)) sel <- c(sel, i)
  .SD[head(sel, 6)]
}, by = condition]
keep[, lab := sprintf("%s:%.2f Mb", sub("^chr", "", chrom), bin_mid/1e6)]
pA <- pA + ggrepel::geom_text_repel(data = keep, aes(gmid, signed, label = lab),
        size = 2.5, colour = COL_INK2, segment.colour = COL_INK2, segment.size = 0.2,
        min.segment.length = 0, box.padding = 0.2, point.padding = 0.15,
        direction = "y", nudge_y = 6, max.overlaps = Inf, seed = 1)
pA <- (pub_track(fit) / pA) + plot_layout(heights = c(1.3, 10)) +
  plot_annotation(title = "A. Signed QTL recurrence across segregant pools",
    subtitle = "Height = % of pools in that condition whose LOD-drop interval covers the position",
    caption = PUB_CAPTION, theme = theme_gw())
ggsave(file.path(OUT, "protoA_recurrence.png"), pA, width = 13, height = 7.5, dpi = 140, bg = "white")

## ---------- B. per-pool interval stack (karyogram) --------------------------
mk_stack <- function(cond) {
  x <- fit[condition == cond]
  ord <- unique(x[order(cross, sample), .(sample, cross)])
  ord[, y := .I]
  x <- merge(x, ord, by = c("sample", "cross"))
  crosslab  <- ord[, .(ymin = min(y) - 0.5, ymax = max(y) + 0.5, mid = mean(y)), by = cross]
  crossband <- crosslab[seq(1, .N, 2)]
  ggplot(x) +
    geom_rect(data = crossband, aes(xmin = 0, xmax = GENOME_END, ymin = ymin, ymax = ymax),
              fill = "#000000", alpha = 0.04, inherit.aes = FALSE) +
    chrom_bands(0.5, nrow(ord) + 0.5) +
    geom_segment(aes(x = gstart, xend = gend, y = y, yend = y, colour = peak.beta),
                 linewidth = 1.6, lineend = "butt") +
    geom_point(aes(gpeak, y), colour = COL_INK, size = 0.25) +
    scale_beta("colour") +
    scale_y_continuous(breaks = crosslab$mid, labels = crosslab$cross,
                       expand = expansion(add = 0.6)) +
    chrom_axis() +
    labs(title = paste0("B. One row per pool - ", cond),
         subtitle = paste0("Bar spans the LOD-drop interval, tick marks the peak position; ",
                           nrow(ord), " pools grouped by cross"),
         x = "chromosome", y = "cross") +
    theme_gw() + theme(axis.text.y = element_text(size = 6.5, colour = COL_INK2),
                       panel.grid.major.y = element_blank())
}
pB <- mk_stack("Glu")
ggsave(file.path(OUT, "protoB_stack_Glu.png"), pB, width = 13, height = 9, dpi = 140, bg = "white")

## ---------- C. compact condition x genome heat strip ------------------------
CBIN <- 20000
d10 <- bin_depth(fit, binsize = CBIN)
d10[, condition := factor(condition, levels = levels(fit$condition))]
net <- dcast(d10, condition + chrom + gmid ~ direction, value.var = "frac", fill = 0)
setnames(net, c("-1", "1"), c("bbq", "wild"))
net[, netfrac := 100 * (bbq - wild)]
lim <- as.numeric(quantile(abs(net$netfrac), 0.995))
net[, ynum := as.integer(factor(condition, levels = rev(levels(fit$condition))))]
pC <- ggplot(net) +
  geom_rect(aes(xmin = gmid - CBIN*0.55, xmax = gmid + CBIN*0.55,
                ymin = ynum - 0.42, ymax = ynum + 0.42, fill = netfrac)) +
  scale_fill_gradient2(low = COL_WILD, mid = COL_MID, high = COL_BBQ, midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish,
                       name = paste0("net % of pools (BBQ - wild), capped at ", round(lim)),
                       guide = guide_colourbar(barwidth = 10, barheight = 0.4, title.vjust = 1)) +
  scale_y_continuous(breaks = seq_along(levels(fit$condition)),
                     labels = labf[rev(levels(fit$condition))], expand = expansion(add = 0.15)) +
  chrom_axis() +
  labs(title = "C. Condition x genome, one strip each",
       subtitle = paste0(CBIN/1000, " kb bins; colour = net share of pools with an interval covering the bin"),
       x = "chromosome", y = NULL) +
  theme_gw() + theme(panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "protoC_strip.png"), pC, width = 13, height = 3.6, dpi = 140, bg = "white")

## ---------- D. Manhattan + interval whiskers (the proposed form) ------------
pD <- ggplot(fit) +
  chrom_bands(0, 0.55) +
  pub_rules() +
  geom_segment(aes(x = gstart, xend = gend, y = abs(peak.beta), yend = abs(peak.beta),
                   colour = peak.beta), linewidth = 0.3, alpha = 0.5) +
  geom_point(aes(gpeak, abs(peak.beta), colour = peak.beta), size = 0.7) +
  facet_grid(condition ~ ., labeller = labeller(condition = labf)) +
  scale_beta("colour") +
  chrom_axis() +
  labs(
       x = "chromosome", y = "|AF shift|") +
  theme_gw()
pD <- (pub_track(fit) / pD) + plot_layout(heights = c(1.3, 10)) +
  plot_annotation(title = "D. Manhattan-style: peak point + interval whisker",
    ## V7: was a hardcoded "LOD is Inf for 816/3308 peaks" -- both numbers from a
    ## superseded run and false for v6b, which has 0 non-finite LODs.  Computed
    ## now, and the clause disappears entirely when there is nothing to report.
    subtitle = paste0("y = |allele-frequency shift| at the peak",
                      { .k <- sum(!is.finite(fit$peak.LOD))
                        if (.k) sprintf("  (LOD non-finite for %d of %d peaks)",
                                        .k, nrow(fit)) else "" }),
    caption = PUB_CAPTION, theme = theme_gw())
ggsave(file.path(OUT, "protoD_manhattan.png"), pD, width = 13, height = 7.5, dpi = 140, bg = "white")


## ---------- E. condition-specific contrasts (stress vs Glu) -----------------
st <- read_peaks(STRESS_CSV)
st[, sample := paste(cross, well, condition, sep = "_")]   # no single 'sample' col in stress output
st[, condition := factor(condition, levels = c("Ace", "Caf", "NaCl"))]
ns <- st[, .(n = uniqueN(sample)), by = condition]
labs <- setNames(sprintf("%s vs Glu\n(n=%d)", ns$condition, ns$n), ns$condition)
ds <- bin_depth(st, binsize = 2000)
ds[, condition := factor(condition, levels = levels(st$condition))]
ds[, signed := fifelse(direction < 0, 100 * frac, -100 * frac)]
slo <- min(ds$signed) * 1.04; shi <- max(ds$signed) * 1.04
pE <- ggplot(ds, aes(gmid, signed)) +
  chrom_bands(slo, shi) +
  geom_area(aes(fill = factor(direction, c(-1, 1))), position = "identity") +
  geom_hline(yintercept = 0, colour = COL_INK2, linewidth = 0.3) +
  facet_grid(condition ~ ., labeller = labeller(condition = labs)) +
  scale_fill_manual(values = c(`-1` = COL_BBQ, `1` = COL_WILD),
                    labels = c("BBQ allele favoured in stress", "wild-parent allele favoured in stress"),
                    name = NULL) +
  scale_y_continuous(labels = function(x) paste0(abs(x), "%"), limits = c(slo, shi),
                     expand = expansion(mult = 0)) +
  chrom_axis() +
  labs(title = "E. Condition-SPECIFIC QTL: stress vs glucose contrast",
       subtitle = "Same encoding as A, but on the paired stress-vs-Glu contrast - direction is balanced here",
       x = "chromosome", y = "% of pools") +
  theme_gw()
ggsave(file.path(OUT, "protoE_stress.png"), pE, width = 13, height = 6, dpi = 140, bg = "white")


## ---- the chrVII BBQ growth QTL, now reported as a finding ------------------
## Unmasked in MASK_REPORT.  Two numbers carry the argument that this is a real
## growth locus of BBQ's and not a reference or alignment artefact:
##   (i)  in the fitness contrast it is BBQ-favoured in essentially every pool
##        that has it, across all four conditions including plain glucose;
##   (ii) in the paired stress-vs-glucose contrast it almost vanishes - which is
##        exactly what a locus shared by every condition must do when one
##        condition is contrasted against another.  A reference artefact would
##        not cancel that way.
## Delete this block if the locus gets its own script.
bbq    <- MASK_ALL[what == "BBQ growth QTL"]
in_bbq <- function(d) as.character(d$seqnames) == bbq$chrom &
                      d$peak.pos >= bbq$lo & d$peak.pos <= bbq$hi

bbq_rows <- function(d, lab) {
  agg <- function(x, cd) { h <- x[in_bbq(x)]; data.table(
    contrast = lab, condition = cd,
    pools_with_locus = uniqueN(h$sample),   pools_total = uniqueN(x$sample),
    pct_pools        = round(100 * uniqueN(h$sample) / uniqueN(x$sample), 1),
    crosses_with_locus = uniqueN(h$cross),  crosses_total = uniqueN(x$cross),
    peaks_at_locus   = nrow(h),             peaks_total = nrow(x),
    pct_peaks        = round(100 * nrow(h) / nrow(x), 1),
    pct_BBQ_favoured = if (nrow(h)) round(100 * mean(h$favoured == "BBQ"), 1) else NA_real_,
    median_abs_beta  = if (nrow(h)) round(median(abs(h$peak.beta)), 3) else NA_real_,
    median_peak_pos  = if (nrow(h)) as.integer(median(h$peak.pos)) else NA_integer_) }
  rbindlist(c(lapply(levels(droplevels(d$condition)),
                     function(cd) agg(d[condition == cd], cd)),
              list(agg(d, "ALL"))))
}
bbq_tab <- rbind(bbq_rows(fit, "fitness vs 0.5"), bbq_rows(st, "stress vs Glu"))
fwrite(bbq_tab, file.path(OUT, "bbq_chr07_locus.csv"))
cat("\n== chrVII", format(bbq$lo, big.mark = ","), "-", format(bbq$hi, big.mark = ","),
    "BBQ growth QTL - reported, not masked ==\n")
print(bbq_tab)

## peak.LOD saturates to Inf in a non-trivial share of peaks; beta does not.
## Quote beta, not LOD, until this is traced back through the caller.
ninf <- sum(!is.finite(fit$peak.LOD))
if (ninf) message(sprintf(
  "NOTE: %d of %d fitness peaks have non-finite peak.LOD (%d of them at the chrVII locus). Do not quote a LOD summary.",
  ninf, nrow(fit), sum(!is.finite(fit[in_bbq(fit)]$peak.LOD))))

## ---- F. cross-resolved recurrence surface ---------------------------------
## x = genome position, y/z = cross, height = fraction of that cross's pools
## whose LOD-drop interval covers the position.  Marginals: top = pooled over
## crosses (what panel A shows), right = pooled over the genome.
BIN <- 10000
AMP <- 2.1                                   # ridge height in row units
COL_MIX <- "#ad6ee0"   # diverging midpoint here means "peaks of BOTH directions
                       # stack here" - a real state, so it gets a real colour,
                       # unlike the gray midpoint of the beta scales in B/D
                       # where the midpoint just means "no shift".

## ---- binned interval depth, per cross, split by direction ------------------
bin_cross <- function(p, binsize = BIN) {
  out <- list()
  for (cr in unique(p$cross)) {
    n_c <- uniqueN(p[cross == cr]$sample)
    for (dirn in c(-1, 1)) {
      x <- p[cross == cr & direction == dirn]
      if (!nrow(x)) next
      gr <- GRanges(factor(as.character(x$seqnames), levels = CHROMS),
                    IRanges(pmax(x$start, 1), x$end),
                    seqinfo = Seqinfo(CHROMS, unname(CHRLEN)))
      cv <- coverage(gr)
      for (ch in CHROMS) {
        r  <- cv[[ch]]; L <- CHRLEN[[ch]]
        bs <- seq(1, L, by = binsize); be <- pmin(bs + binsize - 1, L)
        v  <- viewMeans(Views(r, start = bs, end = pmin(be, length(r))))
        out[[length(out) + 1]] <- data.table(cross = cr, n_pool = n_c, direction = dirn,
          chrom = ch, bin_mid = (bs + be)/2, frac = as.numeric(v)/n_c)
      }
    }
  }
  d <- dcast(rbindlist(out), cross + n_pool + chrom + bin_mid ~ direction,
             value.var = "frac", fill = 0)
  for (nm in c("-1", "1")) if (!nm %in% names(d)) d[, (nm) := 0]
  setnames(d, c("-1", "1"), c("f_bbq", "f_wild"))
  d[, tot := f_bbq + f_wild]
  d[, net := fifelse(tot > 0, (f_bbq - f_wild)/tot, 0)]   # +1 all BBQ, -1 all wild
  d[, gmid := gx(chrom, bin_mid)][]
}

## same, pooled over every cross -> the top marginal
bin_all <- function(p, binsize = BIN) {
  n_tot <- uniqueN(p$sample); out <- list()
  for (dirn in c(-1, 1)) {
    x  <- p[direction == dirn]
    if (!nrow(x)) next
    gr <- GRanges(factor(as.character(x$seqnames), levels = CHROMS),
                  IRanges(pmax(x$start, 1), x$end),
                  seqinfo = Seqinfo(CHROMS, unname(CHRLEN)))
    cv <- coverage(gr)
    for (ch in CHROMS) {
      r  <- cv[[ch]]; L <- CHRLEN[[ch]]
      bs <- seq(1, L, by = binsize); be <- pmin(bs + binsize - 1, L)
      v  <- viewMeans(Views(r, start = bs, end = pmin(be, length(r))))
      out[[length(out) + 1]] <- data.table(direction = dirn, chrom = ch,
        bin_mid = (bs + be)/2, frac = as.numeric(v)/n_tot)
    }
  }
  d <- rbindlist(out)
  d[, gmid := gx(chrom, bin_mid)]
  d[, signed := fifelse(direction < 0, frac, -frac)][]   # BBQ up, wild down
}

scale_net <- function() scale_fill_gradient2(
  low = COL_WILD, mid = COL_MIX, high = COL_BBQ, midpoint = 0, limits = c(-1, 1),
  breaks = c(-1, 0, 1), labels = c("all wild", "mixed", "all BBQ"),
  name = "direction of the peaks stacked here",
  guide = guide_colourbar(barwidth = 9, barheight = 0.4, title.vjust = 1))

## one surface: ridges + both marginals, for whatever slice of the table is passed
make_surface <- function(p, title, sub, fade_below = 4, which_cond = NULL) {
  ## no published comparator for this condition -> no overlay, no track
  has_pub <- is.null(which_cond) || which_cond %in% PUB_LANES
  dd  <- bin_cross(p)
  ord <- unique(dd[, .(cross, n_pool)])[order(-n_pool, cross)][, row := .I][]
  dd  <- ord[dd, on = .(cross, n_pool)]
  NR  <- nrow(ord); YTOP <- NR + AMP + 1.2

  ## ridges, opaque, drawn back to front.  The white silhouette is built from
  ## EVERY bin, zeros included, so it drops to the baseline between peaks
  ## instead of flying straight across the empty stretches.
  ridge_layers <- unlist(lapply(rev(seq_len(NR)), function(i) {
    s <- dd[row == i][order(gmid)]
    a <- if (s$n_pool[1] < fade_below) 0.45 else 1
    sil <- rbind(s[, .(gmid, y = row + tot*AMP)],
                 data.table(gmid = c(0, GENOME_END), y = i))[order(gmid)]
    list(geom_line(data = sil, aes(gmid, y), inherit.aes = FALSE,
                   colour = "white", linewidth = 0.4, lineend = "round"),
         geom_rect(data = s[tot > 0], aes(xmin = gmid - BIN/2, xmax = gmid + BIN/2,
                                          ymin = row, ymax = row + tot*AMP, fill = net),
                   alpha = a, inherit.aes = FALSE))
  }), recursive = FALSE)

  pRidge <- ggplot() +
    chrom_bands(0.4, YTOP) +
    pub_rules(which_cond) +
    geom_hline(yintercept = seq_len(NR), colour = COL_GRID, linewidth = 0.25) +
    ridge_layers + scale_net() + chrom_axis() +
    scale_y_continuous(breaks = ord$row,
                       labels = sprintf("%s  (%d)", sub("xBBQ$", "", ord$cross), ord$n_pool),
                       expand = expansion(mult = c(0.005, 0.01)), limits = c(0.4, YTOP)) +
    labs(x = "chromosome", y = "cross  (pools)") +
    theme_gw() + theme(panel.grid.major.y = element_blank(),
                       axis.text.y = element_text(size = 7.5, colour = COL_INK2))

  ad <- bin_all(p)
  pTop <- ggplot(ad) +
    chrom_bands(-max(ad$frac), max(ad$frac)) +
    pub_rules(which_cond) +
    geom_rect(aes(xmin = gmid - BIN/2, xmax = gmid + BIN/2, ymin = 0, ymax = signed,
                  fill = fifelse(direction < 0, "BBQ", "wild")), show.legend = FALSE) +
    scale_fill_manual(values = c(BBQ = COL_BBQ, wild = COL_WILD)) +
    geom_hline(yintercept = 0, colour = COL_INK2, linewidth = 0.3) +
    chrom_axis() +
    scale_y_continuous(labels = function(v) percent(abs(v), accuracy = 1),
                       expand = expansion(mult = 0.02)) +
    labs(y = sprintf("all %d pools", uniqueN(p$sample))) +
    theme_gw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                       plot.margin = margin(2, 2, 0, 2))

  GLEN <- sum(CHRLEN)
  cov_frac <- p[, {
    gr <- GRanges(factor(as.character(seqnames), levels = CHROMS),
                  IRanges(pmax(start, 1), end), seqinfo = Seqinfo(CHROMS, unname(CHRLEN)))
    .(covered = sum(width(reduce(gr)))/GLEN,
      net = (sum(direction < 0) - sum(direction > 0))/.N)
  }, by = cross]
  marg <- ord[cov_frac, on = "cross"]

  pRight <- ggplot(marg) +
    geom_col(aes(covered, row, fill = net), width = 0.72, orientation = "y",
             show.legend = FALSE) +
    scale_net() +
    scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1),
                       limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(0.4, YTOP), expand = expansion(mult = c(0.005, 0.01))) +
    labs(x = "genome\nunder a peak") +
    theme_gw() + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                       panel.grid.major.y = element_blank(),
                       plot.margin = margin(2, 2, 0, 0))

  ptrack <- if (has_pub) pub_track(p, size = 2.2, which_cond = which_cond) else plot_spacer()
  list(plot = (ptrack + plot_spacer() + pTop + plot_spacer() + pRidge + pRight +
               plot_layout(widths = c(6, 1), heights = c(1.3, 1.4, 9), guides = "collect")) +
         plot_annotation(title = title, subtitle = sub,
           caption = paste("crosses ordered by pool count, best-powered in front.  Top: the same profile pooled over every pool in the panel, BBQ up / wild down.",
                           "\nRight: share of the genome each cross puts under at least one interval.  Faded rows are crosses with fewer than",
                           fade_below, "pools here - at one pool every interval it calls is a full-height ridge.",
                           "\n", PUB_CAPTION),
           theme = theme_gw()),
       nrow = NR)
}

## all conditions together
sf <- make_surface(fit,
  "QTL recurrence surface: every cross, every condition",
  sprintf("ridge height = fraction of that cross's pools with a LOD-drop interval covering the position (%d kb bins); colour = which parent the stacked peaks favour", BIN/1000),
  fade_below = 4)
ggsave(file.path(OUT, "protoF_surface.png"), sf$plot,
       width = 14, height = 2.9 + 0.24*sf$nrow, dpi = 140, bg = "white")

## and one per condition
for (cd in levels(fit$condition)) {
  s <- make_surface(fit[condition == cd],
    sprintf("QTL recurrence surface: %s", cd),
    sprintf("%s only - %d pools across %d crosses; ridge height = fraction of that cross's %s pools with an interval covering the position (%d kb bins)",
            cd, uniqueN(fit[condition == cd]$sample), uniqueN(fit[condition == cd]$cross), cd, BIN/1000),
    fade_below = 2, which_cond = cd)
  ggsave(file.path(OUT, sprintf("protoF_surface_%s.png", cd)), s$plot,
         width = 14, height = 2.9 + 0.24*s$nrow, dpi = 140, bg = "white")
}

cat("done\n")


invisible(TRUE)
}


# ════════════════════════════════════════════════════════════════════
# qtl_candidates.R — rank candidate causal genes under each locus
# ════════════════════════════════════════════════════════════════════
qtl_candidates_run = function(TAG, PEAK_DIR, OUT, GFF, KEEP_CHR7) {

## ---- configuration, replaces qtl_candidates.R lines 1-101 ----------------
suppressMessages({library(data.table); library(GenomicRanges)})
.t0 <- Sys.time(); tick <- function(what)
  message(sprintf("  [%5.1fs] %s", as.numeric(difftime(Sys.time(), .t0, units = "secs")), what))
PROVEAN_CUT <- -2.5     # standard PROVEAN deleterious threshold
PHYLOP_MIN  <- 1.0      # a PROVEAN-deleterious call must also sit at a conserved site
DROP_DUBIOUS<- TRUE     # dubious ORFs are annotation noise, not candidate genes
PROM_BP     <- 1000     # upstream window for assigning intergenic variants
LOCUS_W     <- 25000    # summits within this distance are one locus
TOP_N       <- 5        # candidates listed per locus in the *_top file
MIN_CROSS   <- 2        # a locus needs this many crosses before it is ranked
CONS_FRAC   <- 0.5      # consensus region = where interval depth >= this x the locus max
MIN_COV_FRAC<- 0.5      # a candidate must lie inside this fraction of the locus's intervals

PEAKS <- file.path(PEAK_DIR,
                   if (TAG == "fitness") "fit_peaks.csv" else "stress_peaks.csv")
VARS  <- file.path(PEAK_DIR, sprintf("qtl_variants_%s.csv", TAG))
.miss <- c(PEAKS, VARS)[!file.exists(c(PEAKS, VARS))]
if (length(.miss))
  stop("not in ", PEAK_DIR, ":\n  ", paste(basename(.miss), collapse = "\n  "),
       call. = FALSE)
if (!file.exists(GFF)) stop("no S288C GFF at ", GFF, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
message("peaks    : ", PEAKS)
message("gff      : ", GFF)
message("tables  -> ", OUT)

CHRLEN <- c(chrI=230218, chrII=813184, chrIII=316620, chrIV=1531933, chrV=576874,
            chrVI=270161, chrVII=1090940, chrVIII=562643, chrIX=439888, chrX=745751,
            chrXI=666816, chrXII=1078177, chrXIII=924431, chrXIV=784333,
            chrXV=1091291, chrXVI=948066)

## regions whose recurrence is already explained - see qtl_landscape.R.
## chrVII is the BBQ growth QTL: a finding, not an artefact, so it is KEPT here and
## its genes are ranked like any other locus.  It is present in 185 of 217 pools, so
## it will sit at the top of any recurrence-weighted list - set KEEP_CHR7 to FALSE to
## take it back out if it swamps the ranking.
## V6: settable from the driver (xQTL_analysis6.R §0 keep_chr7) so the run
## manifest records which way it ran. Default is unchanged.
## KEEP_CHR7 is a function argument in V6.1 (driver §0 keep_chr7).
MASK <- data.table(what  = c("MAT", "HO", "BBQ growth QTL"),
                   chrom = c("chrIII", "chrIV", "chrVII"),
                   lo    = c(125000, 15000, 440000),
                   hi    = c(215000, 60000, 600000))
if (KEEP_CHR7) MASK <- MASK[what != "BBQ growth QTL"]

## Same shoulder rule as qtl_landscape.R: the summit test alone lets a sweep through
## when the peak is called just outside the mask but its interval spans it.
SHOULDER_FRAC <- 0.5

## published candidates, for annotation only - never used in the ranking
PUB_GENES <- c("TOR1","CYR1","WHI2","YGK3","PDE2","ENA1","ENA2","ENA5","OSM1","STK2",
               "MSB2","MTL1","GSC2","TCO89","RTT10","RAD5","MKT1","CUP1","HAP1")

## ---- 1. gene models --------------------------------------------------------
message("reading gene models ...")
ZCAT <- if (grepl("\\.gz$", GFF)) {
  cand <- c("gzcat", "zcat", "gunzip -c")          # macOS, Linux, fallback
  ok   <- cand[nzchar(Sys.which(sub(" .*", "", cand)))]
  if (!length(ok)) stop("no gzip decompressor on PATH (tried gzcat, zcat, gunzip)")
  ok[1]
} else "cat"
gff <- fread(cmd = sprintf("%s %s | awk -F'\\t' '$3==\"gene\"'", ZCAT, shQuote(GFF)),
             sep = "\t", header = FALSE, quote = "", showProgress = FALSE)
if (!nrow(gff))
  stop("no 'gene' records parsed from ", GFF, " via '", ZCAT, "'")
setnames(gff, c("chr","src","type","start","end","score","strand","phase","attr"))
gff <- gff[chr %in% names(CHRLEN)]
att <- function(k) ifelse(grepl(paste0("(^|;)", k, "="), gff$attr),
                          sub(paste0(".*(^|;)", k, "=([^;]+).*"), "\\2", gff$attr), NA_character_)
GENES <- gff[, .(chrom = chr, gstart = start, gend = end, strand,
                 gene = att("ID"), name = att("gene"), orf = att("orf_classification"))]
GENES[is.na(name) | !nzchar(name), name := gene]
GENES[, `:=`(len = gend - gstart + 1, tss = fifelse(strand == "+", gstart, gend),
             gmid = (gstart + gend)/2)]
GENES <- GENES[!is.na(gene)]
message(sprintf("  %d genes (%d named, %d dubious)", nrow(GENES),
                sum(GENES$name != GENES$gene), sum(GENES$orf %in% "Dubious")))

## ---- 2. peaks -> loci ------------------------------------------------------
## A locus is a cluster of PEAK SUMMITS within LOCUS_W, per condition.  Chaining
## overlapping intervals instead would fuse whole chromosome arms: in glucose the
## intervals are wide enough to merge chrXVI into one 838 kb "locus".
tick("gene models"); message("reading peaks ...")
pk <- fread(PEAKS,
            showProgress = FALSE)
if ("databin" %in% names(pk)) pk <- pk[databin == "main"]
## stress_peaks.csv has no `sample` column - build a pool id so the two tables
## can share the rest of the code (and so uniqueN never picks up base::sample)
if (!"sample" %in% names(pk)) pk[, sample := paste(condition, cross, well, sep = "|")]
pk[, masked := FALSE]
for (i in seq_len(nrow(MASK))) {
  pk[seqnames == MASK$chrom[i] & peak.pos >= MASK$lo[i] & peak.pos <= MASK$hi[i], masked := TRUE]
  ## shoulder: interval covers more than half the masked region
  pk[seqnames == MASK$chrom[i] &
     pmax(0, pmin(end, MASK$hi[i]) - pmax(start, MASK$lo[i])) >
       SHOULDER_FRAC * (MASK$hi[i] - MASK$lo[i]), masked := TRUE]
}
message(sprintf("  masked %d of %d peaks (summit or >%.0f%% shoulder over MAT/HO%s)",
                sum(pk$masked), nrow(pk), 100 * SHOULDER_FRAC,
                if (KEEP_CHR7) "" else "/chrVII"))
pk <- pk[!(masked)]
pk[, end := pmin(end, CHRLEN[as.character(seqnames)])]

setorder(pk, condition, seqnames, peak.pos)
pk[, locus_id := cumsum(c(TRUE, diff(peak.pos) > LOCUS_W)), by = .(condition, seqnames)]
pk[, locus := sprintf("%s|%s:%.0fkb", condition, seqnames, median(peak.pos)/1e3),
   by = .(condition, seqnames, locus_id)]
LOCI <- pk[, .(crosses = uniqueN(cross), pools = uniqueN(sample),
               summit = as.integer(median(peak.pos)),
               lo = min(start), hi = max(end),
               bbq_low = mean(low.parent == "BBQ"), max_lod = max(peak.LOD[is.finite(peak.LOD)])),
           by = .(locus, condition, chrom = seqnames)]
message(sprintf("  %d peaks -> %d loci (%d with >= %d crosses)",
                nrow(pk), nrow(LOCI), sum(LOCI$crosses >= MIN_CROSS), MIN_CROSS))

## ---- 3. consensus region and gene coverage --------------------------------
## The union of a locus's intervals is far too generous - one 400 kb interval from
## a single cross drags in half a chromosome.  Take instead the CONSENSUS REGION:
## where the interval pileup is at least CONS_FRAC of its maximum, the same
## 50%-of-max plateau rule used to place the masks.  Genes outside it are not
## candidates for that locus however wide somebody's interval was.
## coverage/slice/start/end are namespace-qualified: dplyr masks slice() and an
## attached tidyverse turns the Rle into an unhandled class at this line.
CONS <- rbindlist(lapply(split(pk, pk$locus), function(d) {
  ch  <- as.character(d$seqnames[1])
  cv  <- GenomicRanges::coverage(GRanges(ch, IRanges(d$start, pmin(d$end, CHRLEN[[ch]])),
                                         seqinfo = Seqinfo(ch, CHRLEN[[ch]])))[[ch]]
  thr <- max(cv) * CONS_FRAC
  sl  <- IRanges::slice(cv, lower = thr, rangesOnly = TRUE)
  data.table(locus = d$locus[1], chrom = ch,
             cstart = BiocGenerics::start(sl), cend = BiocGenerics::end(sl))
}))
CONSW <- CONS[, .(cons_kb = round(sum(cend - cstart + 1)/1000, 1)), by = locus]
message(sprintf("  consensus regions: median %.1f kb (IQR %.1f-%.1f)",
                median(CONSW$cons_kb), quantile(CONSW$cons_kb, .25), quantile(CONSW$cons_kb, .75)))

## genes inside the consensus region of their locus
hcg <- findOverlaps(GRanges(CONS$chrom, IRanges(CONS$cstart, CONS$cend)),
                    GRanges(GENES$chrom, IRanges(GENES$gstart, GENES$gend)))
INCONS <- unique(data.table(locus = CONS$locus[queryHits(hcg)],
                            gene  = GENES$gene[subjectHits(hcg)]))

## coverage denominator: how many of the locus's crosses have that gene in their own
## interval.  From the INTERVALS, not the variant table - a covered gene with no
## segregating variant still belongs in the denominator.
cover <- unique(pk[, .(locus, condition, cross, seqnames, start, end)])
hc <- findOverlaps(GRanges(as.character(cover$seqnames), IRanges(cover$start, cover$end)),
                   GRanges(GENES$chrom, IRanges(GENES$gstart, GENES$gend)))
COVER <- unique(data.table(locus = cover$locus[queryHits(hc)],
                           condition = cover$condition[queryHits(hc)],
                           cross = cover$cross[queryHits(hc)],
                           gene  = GENES$gene[subjectHits(hc)]))
COVER <- merge(COVER, INCONS, by = c("locus","gene"))
message(sprintf("  %d genes in consensus regions, %d (locus,cross,gene) coverage rows",
                uniqueN(INCONS$gene), nrow(COVER)))

## ---- 4. variants, classified and sided ------------------------------------
tick("loci + coverage"); message("reading variants ...")
v <- fread(VARS, showProgress = FALSE)
setnames(v, "comparison", "cross")
v <- v[chr %in% names(CHRLEN)]

## the dump carries no condition column - recover it from the peak table
v[, c("pchr","prange") := tstrsplit(peak, ":", fixed = TRUE)]
v[, c("pstart","pend")  := tstrsplit(prange, "-", fixed = TRUE)]
v[, `:=`(pstart = as.integer(pstart), pend = as.integer(pend))]
key <- unique(pk[, .(cross, pchr = as.character(seqnames), pstart = start, pend = end,
                     condition, locus)])
v <- merge(v, key, by = c("cross","pchr","pstart","pend"), allow.cartesian = TRUE)
v <- v[!is.na(locus)]
message(sprintf("  %d variant rows after masking and the condition join", nrow(v)))

## intergenic variants -> nearest TSS within PROM_BP
v[, prom_dist := NA_integer_]
noid <- is.na(v$GENEID) | !nzchar(v$GENEID)
if (any(noid)) {
  iv <- v[noid, .(chr, pos)]
  nr <- distanceToNearest(GRanges(iv$chr, IRanges(iv$pos, iv$pos)),
                          GRanges(GENES$chrom, IRanges(GENES$tss, GENES$tss)))
  gs <- rep(NA_character_, nrow(iv)); ds <- rep(NA_integer_, nrow(iv))
  gs[queryHits(nr)] <- GENES$gene[subjectHits(nr)]
  ds[queryHits(nr)] <- mcols(nr)$distance
  v[noid, `:=`(GENEID = fifelse(!is.na(ds) & ds <= PROM_BP, gs, NA_character_),
               prom_dist = ds)]
}
v <- v[!is.na(GENEID) & nzchar(GENEID)]

## Every damaging call must sit at a CONSERVED site, LOF included.  Per distinct
## allele the stop/frameshift rate is 0.4% in verified genes, 2.5% in
## uncharacterized ORFs and 7.2% in dubious ones, and the median phyloP at those
## calls is 1.12 / 0.13 / 0.20 - a "truncated" ORF that nothing else in the clade
## conserves is an ORF that is not a gene.  Ungated, YCR006C (uncharacterized,
## 1.1 kb, 11 LOF alleles = 9.9 per kb, median phyloP 0.13) outranks every real
## gene in its locus.  The same phyloP floor already gates deleterious missense;
## applying it to LOF as well costs about half the verified-gene LOF calls.
v[, class := fifelse(CONSEQUENCE %in% c("stop_gained","frameshift") &
                     !is.na(phyloP) & phyloP >= PHYLOP_MIN, "lof",
              fifelse(CONSEQUENCE %in% c("stop_gained","frameshift"), "lof_unconserved",
              fifelse(CONSEQUENCE == "nonsynonymous" & !is.na(provean) & provean <= PROVEAN_CUT &
                      !is.na(phyloP) & phyloP >= PHYLOP_MIN, "deleterious",
              fifelse(CONSEQUENCE == "nonsynonymous", "missense",
              fifelse(is.na(CONSEQUENCE) & !is.na(prom_dist) & prom_dist <= PROM_BP,
                      "regulatory", "silent")))))]
v[, damaging := class %in% c("lof","deleterious")]
## The phyloP gate is not cosmetic.  PROVEAN calls 39% of variants in dubious ORFs
## deleterious against 3.4% in verified genes, and the median phyloP at those calls
## is 0.17 against 6.68 - it is scoring alignments that do not exist.  Ungated, the
## whole candidate list is dubious and uncharacterised ORFs.

## which parent carries the ALT - i.e. the allele the annotation describes
gt <- tstrsplit(v$genotypes, ";", fixed = TRUE)
p1 <- sub("=.*", "", gt[[2]]); d1 <- suppressWarnings(as.integer(sub(".*=", "", gt[[2]])))
p2 <- sub("=.*", "", gt[[3]]); d2 <- suppressWarnings(as.integer(sub(".*=", "", gt[[3]])))
v[, alt_parent := fifelse(d1 > d2, p1, fifelse(d2 > d1, p2, NA_character_))]
v[, side := fifelse(alt_parent == "BBQ", "bbq", "wild")]
v[, coherent := !is.na(alt_parent) & alt_parent == low.parent]   # damaging allele on the low-fitness parent
message("  classes: ", paste(sprintf("%s=%d", names(table(v$class)), table(v$class)), collapse = "  "),
        "\n  damaging by side: ",
        paste(sprintf("%s=%d", names(table(v[(damaging)]$side)), table(v[(damaging)]$side)), collapse = "  "))

tick("variants classified")
## ---- 5. per (locus, gene) evidence ----------------------------------------
## Distinct ALLELES, not variant rows and not crosses: a wild haplotype shared by
## twenty related strains is one allele, and a BBQ allele is one allele in all 32.
dam <- v[(damaging)]
E <- dam[, .(
      alleles_wild = uniqueN(gEdit[side == "wild"]),
      crosses_wild = uniqueN(cross[side == "wild"]),
      alleles_bbq  = uniqueN(gEdit[side == "bbq"]),
      n_lof        = uniqueN(gEdit[class == "lof"]),
      n_lof_wild   = uniqueN(gEdit[class == "lof" & side == "wild"]),
      coh_wild     = uniqueN(cross[side == "wild" & coherent]),
      min_provean  = suppressWarnings(min(provean, na.rm = TRUE)),
      max_phylop   = suppressWarnings(max(phyloP, na.rm = TRUE)),
      wild_strains = paste(sort(unique(sub("xBBQ$", "", cross[side == "wild"]))), collapse = ",")),
    by = .(locus, condition, gene = GENEID)]
R <- v[class == "regulatory", .(n_reg = uniqueN(gEdit)), by = .(locus, condition, gene = GENEID)]
S <- v[, .(n_seg = uniqueN(gEdit)), by = .(locus, condition, gene = GENEID)]   # all segregating variants

D <- merge(COVER[, .(locus, condition, gene, crosses_cov = .N),
                 by = .(locus, condition, gene)][, .(locus, condition, gene, crosses_cov)],
           E, by = c("locus","condition","gene"), all.x = TRUE)
D <- merge(D, R, by = c("locus","condition","gene"), all.x = TRUE)
D <- merge(D, S, by = c("locus","condition","gene"), all.x = TRUE)
for (j in c("alleles_wild","crosses_wild","alleles_bbq","n_lof","n_lof_wild","coh_wild","n_reg","n_seg"))
  set(D, which(is.na(D[[j]])), j, 0L)
D[is.na(wild_strains), wild_strains := ""]
D[!is.finite(min_provean), min_provean := NA_real_]
D[!is.finite(max_phylop),  max_phylop  := NA_real_]
D <- merge(D, GENES[, .(gene, name, chrom, gstart, gend, len, orf)], by = "gene")
D <- merge(D, LOCI[, .(locus, summit, locus_crosses = crosses, locus_pools = pools,
                       bbq_low, locus_lod = max_lod)], by = "locus")
D <- merge(D, CONSW, by = "locus", all.x = TRUE)
D[, dist_kb := round(abs((gstart + gend)/2 - summit)/1000, 1)]

tick("evidence assembled")
## ---- 6. null for the allelic series ---------------------------------------
## Poisson on the count of distinct damaging wild alleles, with the rate set by how
## polymorphic the gene actually is at this locus (n_seg) times the genome-wide
## damaging fraction.  Conditioning on n_seg is what stops long, diverged genes
## from topping the list on divergence alone.
frac_dam <- uniqueN(dam[side == "wild"]$gEdit) / uniqueN(v$gEdit)
D[, exp_wild := round(n_seg * frac_dam, 3)]
D[, p_allelic := ppois(alleles_wild - 1, pmax(exp_wild, 1e-6), lower.tail = FALSE)]
D[alleles_wild == 0, p_allelic := 1]
message(sprintf("  genome-wide damaging fraction of segregating alleles: %.4f", frac_dam))

if (DROP_DUBIOUS) {
  nd <- sum(D$orf %in% "Dubious")
  D  <- D[!(orf %in% "Dubious")]
  message(sprintf("  dropped %d (locus,gene) rows on dubious ORFs", nd))
}
D[, dam_frac := round(alleles_wild/pmax(n_seg, 1), 3)]
CAND <- D[locus_crosses >= MIN_CROSS]
## the causal gene has to be inside most of the detecting crosses' intervals
CAND[, cov_frac := round(crosses_cov/locus_crosses, 2)]
CAND <- CAND[cov_frac >= MIN_COV_FRAC]
## BH within a locus: the question is "which gene in THIS interval", not a
## genome-wide scan, and correcting across every gene of every locus buries it
CAND[, q_allelic := p.adjust(p_allelic, "BH"), by = locus]
## coherence: is the damaging allele on the low-fitness parent more often than chance
CAND[, p_coh := fifelse(crosses_wild > 0,
        pbinom(coh_wild - 1, crosses_wild, 0.5, lower.tail = FALSE), NA_real_)]

## ---- 7. composite, used to order within a locus ---------------------------
z <- function(x, lo, hi) pmin(pmax((x - lo)/(hi - lo), 0), 1)
## BUG FIX 2026-09-22.  min_provean and max_phylop are NA whenever a gene has
## no wild MISSENSE variant to score -- PROVEAN scores missense only.  Those two
## NAs propagated through the sum, and the `is.na(score) := 0` below then turned
## the whole composite into a hard zero.  4,173 of 5,160 rows (81%) were scored
## 0 for that reason alone, INCLUDING 36 of the 65 rows carrying a wild
## loss-of-function allele -- so the term weighted 1.0 for "outright loss of
## function" could never fire on its own, and the most damaging allele class was
## the one most reliably discarded.  WHI2, whose confirmed allele in Treusch
## et al. is a frameshift, was zeroed in both instances where one of our wild
## parents (ABP) carries a LOF allele.
## A missing term means "no evidence of this kind", which on the z() scale is 0,
## not "no evidence of any kind".  fcoalesce gives each term its own neutral
## value instead of letting one NA erase the other six.
CAND[, score := round(
    3.0 * z(-log10(pmax(p_allelic, 1e-12)), 0, 6) +   # allelic series
    1.2 * z(alleles_wild, 1, 4) +                     # how many independent alleles
    1.0 * z(n_lof_wild, 0, 2) +                       # outright loss of function
    0.7 * fcoalesce(z(-min_provean, 2.5, 8), 0) +     # severity (NA = no missense)
    0.5 * fcoalesce(z(max_phylop, 1, 4), 0) +         # conservation (NA = none scored)
    0.5 * fifelse(!is.na(p_coh) & p_coh < 0.1, 1, 0) +# damaging allele on the low parent
    0.5 * z(10 - dist_kb, 0, 10) -                    # proximity to the summit
    1.5 * fifelse(orf %in% "Dubious", 1, 0), 2)]
## retained as a backstop; with the fcoalesce above it should now never fire
if (anyNA(CAND$score)) {
  warning(sprintf("score still NA for %d rows -- a new NA input has appeared",
                  sum(is.na(CAND$score))))
  CAND[is.na(score), score := 0]
}
CAND[, published := name %in% PUB_GENES]
setorder(CAND, locus, -score, q_allelic)
CAND[, rank_in_locus := seq_len(.N), by = locus]

tick("scored")
## ---- 8. output -------------------------------------------------------------
cols <- c("locus","condition","chrom","summit","cons_kb","locus_crosses","locus_pools","locus_lod",
          "gene","name","orf","len","dist_kb","crosses_cov","cov_frac","n_seg",
          "alleles_wild","crosses_wild","exp_wild","dam_frac","p_allelic","q_allelic",
          "alleles_bbq","bbq_low","n_lof","n_lof_wild","n_reg",
          "coh_wild","p_coh","min_provean","max_phylop",
          "score","rank_in_locus","published","wild_strains")
fwrite(CAND[, ..cols], file.path(OUT, sprintf("qtl_candidates_%s.csv", TAG)))
fwrite(CAND[rank_in_locus <= TOP_N, ..cols],
       file.path(OUT, sprintf("qtl_candidates_%s_top%d.csv", TAG, TOP_N)))

message(sprintf("\n%d (locus,gene) pairs across %d loci; %d alleles-series hits at q < 0.05",
                nrow(CAND), uniqueN(CAND$locus), sum(CAND$q_allelic < 0.05)))
fmt <- function(d) d[, .(locus, gene = name, orf, rk = rank_in_locus,
                         cross = locus_crosses, cov = crosses_cov,
                         seg = n_seg, alW = alleles_wild, exp = exp_wild,
                         q = signif(q_allelic, 2), lofW = n_lof_wild, bbq = alleles_bbq,
                         coh = sprintf("%d/%d", coh_wild, crosses_wild),
                         prov = round(min_provean, 1), kb = dist_kb, score, published)]
cat("\n== top of the list, all conditions ==\n")
print(head(fmt(CAND[order(-score)]), 25))
cat("\n== drug conditions only ==\n")
print(head(fmt(CAND[condition %in% c("Caf","NaCl")][order(-score)]), 20))
cat("\n== where published candidates land in their own locus ==\n")
print(fmt(CAND[published == TRUE][order(rank_in_locus, -score)])[1:20])

## ---- 9. self-check: does the ranking recover genes other people confirmed? --
## The only honest way to know whether the ranking works.  Treusch's named
## candidates were never used to build it, so where they land in their own locus
## is an out-of-sample test.  RAD5/MKT1/CUP1/HAP1 come from conditions this
## experiment never assayed and should sit at chance - they are the control.
CAND[, n_genes := .N, by = locus]
CAND[, pct := (rank_in_locus - 0.5)/n_genes]      # 0 = top of its locus, 1 = bottom
pubs <- CAND[published == TRUE]
if (nrow(pubs) >= 10) {
  set.seed(7)
  gl   <- split(CAND$pct, CAND$locus)[pubs$locus]
  perm <- colMeans(vapply(seq_len(10000), function(i)
            vapply(gl, function(x) x[sample.int(length(x), 1L)], numeric(1)),
            numeric(length(gl))))
  w <- suppressWarnings(wilcox.test(pubs$pct, mu = 0.5, alternative = "less"))
  ## V7: the uniform-rank expectations for rank1 and top5.  These were quoted in
  ## the draft (0.89 and 4.44 for v6b fitness) but emitted by nothing -- they had
  ## to be reconstructed by hand, and an earlier prose value of 4.47 came from a
  ## different run.  Under a uniform rank within each locus, the chance an
  ## instance lands first is 1/n_i and the chance it lands in the top five is
  ## min(5, n_i)/n_i, with n_i the gene count of that instance's locus; the
  ## expected counts are the sums over instances.
  exp_rank1 <- sum(1 / pubs$n_genes)
  exp_top5  <- sum(pmin(5, pubs$n_genes) / pubs$n_genes)
  VAL <- data.table(
    instances = nrow(pubs), genes = uniqueN(pubs$name), loci = uniqueN(pubs$locus),
    median_genes_per_locus = median(pubs$n_genes),
    mean_pct = round(mean(pubs$pct), 3), rank1 = sum(pubs$rank_in_locus == 1),
    exp_rank1 = round(exp_rank1, 2),
    top5 = sum(pubs$rank_in_locus <= 5),
    exp_top5 = round(exp_top5, 2),
    p_wilcox = signif(w$p.value, 3),
    p_perm = signif((1 + sum(perm <= mean(pubs$pct)))/10001, 3))
  fwrite(VAL, file.path(OUT, sprintf("qtl_candidates_%s_validation.csv", TAG)))
  cat("\n== self-check: where published candidates land in their own locus ==\n")
  print(VAL)
  cat("(mean_pct 0.5 = no better than chance; lower is better)\n")
  print(pubs[order(pct), .(locus, gene = name, rank = rank_in_locus, of = n_genes,
                           pct = round(pct, 3), alW = alleles_wild, score)][1:15])
}

message("\nwrote ", file.path(OUT, sprintf("qtl_candidates_%s.csv", TAG)))


invisible(TRUE)
}




# ════════════════════════════════════════════════════════════════════
# qtl_go.R — GO-term enrichment of the genes under our loci  (V7, §18)
# ════════════════════════════════════════════════════════════════════
# Inlined the way §16-17 were, so the driver is again the unit that gets
# copied. Body is qtl_go.R from its first line of real code onward, with
# four changes, all marked "V7:" below:
#
#   1. config by argument instead of QTL_* environment variables, so the
#      run manifest records which annotation ran;
#   2. it runs over BOTH peak sets. Figures 4D and 5D ask CONDITION-
#      SPECIFIC questions -- cell-wall genes in caffeine, salt genes in
#      NaCl, acetate stress -- and by this project's own rule only the
#      matched contrast supports a per-locus condition-specific claim.
#      Fitness loci carry every locus that also acts in glucose, which is
#      the opposite of what those panels are asking. Both are written;
#      quote the contrast one for 4D/5D.
#   3. a dead `expanded <- ... # placeholder` line removed;
#   4. the ancestor expansion rewritten. It built one data.table per
#      (orf, GO) pair and rbindlist-ed tens of thousands of them, which
#      is minutes of work for a one-line result.
#
# TWO STATISTICS, AND THEY ARE NOT INTERCHANGEABLE.
#   p_hyper  the usual hypergeometric. It assumes genes are exchangeable.
#            They are not: our loci are contiguous intervals, so the genes
#            under one are physical neighbours and a single interval over
#            a gene family gives a small p on its own. Ranking only.
#   p_perm   the positional null -- the whole locus set shifted rigidly on
#            a gapless concatenation of the sixteen chromosomes, wrapping
#            at the end. Count, width, spacing and local gene density are
#            preserved and only position is randomised. THIS is the
#            quotable one, and it is the same null the landscape section
#            uses for published-locus recovery.
# -Claude generated
qtl_go_run = function(PEAK_DIR, OUT, GFF, KEEP_CHR7, GAF = "", OBO = "",
                      GO_SOURCE = "", ONT = "BP", NPERM = 1000,
                      MIN_TERM = 5, MAX_TERM = 500, LOD_CUT = 5,
                      PEAK_SETS = c("fitness", "contrast")) {

suppressPackageStartupMessages({library(data.table)})
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CHRLEN <- c(chrI=230218, chrII=813184, chrIII=316620, chrIV=1531933, chrV=576874,
            chrVI=270161, chrVII=1090940, chrVIII=562643, chrIX=439888, chrX=745751,
            chrXI=666816, chrXII=1078177, chrXIII=924431, chrXIV=784333,
            chrXV=1091291, chrXVI=948066)
CHROMS <- names(CHRLEN)
OFFSET <- cumsum(c(0, head(CHRLEN, -1))); names(OFFSET) <- CHROMS   # gapless
GLEN   <- sum(CHRLEN)

## glucose first, then alphabetical.  House order, everywhere.
CONDS <- c("Glu", "Ace", "Caf", "NaCl")
LONG  <- c(Glu = "glucose", Ace = "acetic acid", Caf = "caffeine", NaCl = "sodium chloride")

MASK <- data.table(name = c("MAT", "HO"), chrom = c("chrIII", "chrIV"),
                   lo = c(125000, 15000), hi = c(215000, 60000))
if (!KEEP_CHR7)
  MASK <- rbind(MASK, data.table(name = "chrVII sweep", chrom = "chrVII",
                                 lo = 440000, hi = 600000))
SHOULDER_FRAC <- 0.5

## The three biology questions the figures ask.  Keyword filters over term
## names, not a curated list of GO IDs: naming a wrong ID silently would be
## worse than casting a slightly wide net and showing what came back.
## V7 FIX: these are matched with \\b word boundaries. Without them the
## patterns match INSIDE longer words and the filter quietly stops meaning
## anything: "cation" matches repli-CATION, communi-CATION, appli-CATION, so
## the 2026-09-21 run classified "DNA replication checkpoint signaling" as a
## salt/ion term and swept 129 terms into that panel. "pH" was worse -- it is
## case-insensitive here, so it matched every "phosphate", "morphogenesis",
## "sphingolipid" and "phenotype" in the ontology. A keyword filter that
## silently over-matches is more dangerous than a narrow one, because the
## panel still looks populated.
FOCUS <- list(
  `cell wall`   = "\\bcell wall\\b|\\bglucan\\b|\\bchitin\\b|\\bmannoprotein\\b|\\bcell surface\\b",
  `salt / ion`  = "\\bsalt\\b|\\bosmotic\\b|\\bosmolarity\\b|\\bion homeostasis\\b|\\bcation\\b|\\bsodium\\b|\\bpotassium\\b|\\bion transport\\b|\\bmetal ion\\b",
  `acetate`     = "\\bacetate\\b|\\bacetic\\b|\\borganic acid\\b|\\bcarboxylic acid\\b|\\bmonocarboxylic\\b|\\bweak acid\\b|\\bacidity\\b|\\bpH\\b")

say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
say("qtl_go (V7 §18)")
say("  peaks    : %s", PEAK_DIR)
say("  out      : %s", OUT)
say("  keep_chr7: %s   ontology: %s   permutations: %d", KEEP_CHR7, ONT, NPERM)
say("  peak sets: %s", paste(PEAK_SETS, collapse = ", "))

## --------------------------------------------------------------------- GFF --
if (!file.exists(GFF)) stop("no S288C GFF at ", GFF, call. = FALSE)
say("  gff      : %s", GFF)
ZCAT <- if (grepl("\\.gz$", GFF)) "gzcat" else "cat"
if (ZCAT == "gzcat" && system("command -v gzcat > /dev/null 2>&1") != 0) ZCAT <- "zcat"
gff <- fread(cmd = sprintf("%s %s | awk -F'\t' '$3==\"gene\"'", ZCAT, shQuote(GFF)),
             header = FALSE, sep = "\t", quote = "", showProgress = FALSE)
setnames(gff, c("chr","src","type","start","end","score","strand","phase","attr"))
gff <- gff[chr %in% CHROMS]
att <- function(k) {
  v <- sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", gff$attr)
  ifelse(grepl(paste0("(^|;)", k, "="), gff$attr), v, NA_character_)
}
genes <- data.table(orf = att("ID"), name = att("gene"), chrom = gff$chr,
                    start = gff$start, end = gff$end, note = att("orf_classification"))
genes[is.na(name) | !nzchar(name), name := orf]
genes <- genes[is.na(note) | note != "Dubious"]     # as in qtl_candidates_run
genes <- genes[!is.na(orf)]
genes[, gmid := OFFSET[chrom] + (start + end) / 2]
say("  genes    : %d after dropping dubious ORFs", nrow(genes))

## --------------------------------------------------------- GO annotation ----
src <- GO_SOURCE
if (!nzchar(src)) src <- if (requireNamespace("org.Sc.sgd.db", quietly = TRUE)) "orgdb" else "gaf"
say("  GO source: %s", src)

go2genes <- NULL; term_name <- NULL
if (src == "orgdb") {
  if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE))
    stop("GO_SOURCE=orgdb needs org.Sc.sgd.db and GO.db installed.\n",
         "  install with BiocManager::install(c('org.Sc.sgd.db','GO.db'))\n",
         "  or set GO_SOURCE='gaf' with gaf_file and obo_file in §0.", call. = FALSE)
  suppressPackageStartupMessages({library(org.Sc.sgd.db); library(GO.db)})
  all2 <- as.list(org.Sc.sgdGO2ALLORFS)
  ont_of <- AnnotationDbi::Ontology(GO.db::GOTERM)
  keep <- names(all2)[!is.na(ont_of[names(all2)]) & ont_of[names(all2)] == ONT]
  go2genes <- lapply(all2[keep], unique)
  term_name <- AnnotationDbi::Term(GO.db::GOTERM)[keep]
} else {
  if (!nzchar(GAF) || !nzchar(OBO))
    stop("GO_SOURCE='gaf' needs gaf_file (SGD gene_association) and obo_file (go-basic.obo).\n",
         "  GAF: http://sgd-archive.yeastgenome.org/curation/literature/gene_association.sgd.gz\n",
         "  OBO: http://purl.obolibrary.org/obo/go/go-basic.obo", call. = FALSE)
  .m <- c(GAF, OBO)[!file.exists(c(GAF, OBO))]
  if (length(.m)) stop("GO annotation file(s) absent:\n  ", paste(.m, collapse = "\n  "), call. = FALSE)
  gz <- function(f) if (grepl("\\.gz$", f)) sprintf("%s %s", ZCAT, shQuote(f)) else sprintf("cat %s", shQuote(f))
  gaf <- fread(cmd = sprintf("%s | grep -v '^!'", gz(GAF)), header = FALSE, sep = "\t",
               quote = "", fill = TRUE, showProgress = FALSE)
  ## GAF 2.x: 3 symbol, 5 GO ID, 7 evidence, 9 aspect, 11 synonyms
  gaf <- gaf[V9 == substr(ONT, 1, 1) | (ONT == "BP" & V9 == "P") |
             (ONT == "MF" & V9 == "F") | (ONT == "CC" & V9 == "C")]
  gaf <- gaf[!V7 %in% c("ND")]
  ## map to systematic names: the ORF is usually in the synonym field
  syn <- strsplit(gaf$V11, "\\|")
  orf <- vapply(seq_along(syn), function(i) {
    s <- syn[[i]]; h <- s[grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", s)]
    if (length(h)) h[1] else gaf$V3[i] }, "")
  ann <- data.table(orf = orf, go = gaf$V5)[orf %in% genes$orf]
  say("  GAF      : %d annotations on %d of our genes", nrow(ann), uniqueN(ann$orf))

  ## propagate up the ontology: a gene annotated to a child counts for the parent
  say("  parsing %s for is_a / part_of ...", basename(OBO))
  ol <- readLines(OBO, warn = FALSE)
  starts <- which(ol == "[Term]"); ids <- character(); nm <- character()
  parents <- list()
  want <- unname(c(BP = "biological_process", MF = "molecular_function",
                   CC = "cellular_component")[ONT])
  for (b in seq_along(starts)) {
    from <- starts[b] + 1
    to   <- if (b < length(starts)) starts[b + 1] - 1 else length(ol)
    blk  <- ol[from:to]
    id   <- sub("^id: ", "", grep("^id: GO:", blk, value = TRUE))[1]
    if (is.na(id)) next
    if (any(grepl("^is_obsolete: true", blk))) next
    ns <- sub("^namespace: ", "", grep("^namespace: ", blk, value = TRUE))[1]
    if (!identical(ns, want)) next
    ids <- c(ids, id)
    nm  <- c(nm, sub("^name: ", "", grep("^name: ", blk, value = TRUE))[1])
    p1 <- sub("^is_a: (GO:[0-9]+).*", "\\1", grep("^is_a: GO:", blk, value = TRUE))
    p2 <- sub("^relationship: part_of (GO:[0-9]+).*", "\\1",
              grep("^relationship: part_of GO:", blk, value = TRUE))
    parents[[id]] <- unique(c(p1, p2))
  }
  names(nm) <- ids; term_name <- nm
  say("  ontology : %d %s terms, %d with a parent", length(ids), ONT, length(parents))

  anc_cache <- new.env(parent = emptyenv())
  ancestors <- function(id, seen = character()) {
    if (!is.null(a <- anc_cache[[id]])) return(a)
    p <- parents[[id]]; p <- p[!p %in% seen]
    a <- unique(c(p, unlist(lapply(p, ancestors, seen = c(seen, id, p)))))
    assign(id, a, envir = anc_cache); a
  }
  ## V7: was one data.table per (orf, GO) pair fed to rbindlist -- tens of
  ## thousands of allocations for a two-column result.  Expand each DISTINCT
  ## term once, then index, which is the same answer in a fraction of the time.
  pairs  <- unique(ann)
  uterms <- unique(pairs$go)
  expand <- lapply(uterms, function(g) unique(c(g, ancestors(g))))
  names(expand) <- uterms
  reps  <- lengths(expand[pairs$go])
  full  <- unique(data.table(orf = rep(pairs$orf, reps),
                             go  = unlist(expand[pairs$go], use.names = FALSE)))
  full  <- full[go %in% ids]
  go2genes <- split(full$orf, full$go)
  say("  propagated: %d (gene, term) pairs from %d direct", nrow(full), nrow(pairs))
}
universe <- sort(unique(genes$orf))
go2genes <- lapply(go2genes, function(g) intersect(g, universe))
sz <- lengths(go2genes)
go2genes <- go2genes[sz >= MIN_TERM & sz <= MAX_TERM]
term_name <- term_name[names(go2genes)]
say("  GO terms : %d in %s with %d-%d genes, universe %d genes",
    length(go2genes), ONT, MIN_TERM, MAX_TERM, length(universe))
if (!length(go2genes)) stop("no usable GO terms; check the annotation source", call. = FALSE)

## --------------------------------------------------- genes under the loci ----
genes_in <- function(iv) {
  if (!nrow(iv)) return(character())
  unique(unlist(lapply(seq_len(nrow(iv)), function(i)
    genes[chrom == iv$chrom[i] & end >= iv$lo[i] & start <= iv$hi[i], orf])))
}
## Circular shift of the whole locus set on the gapless genome, as in the
## landscape section: count, width and spacing preserved, position randomised.
## An interval running past a chromosome end continues into the next one and is
## split there rather than clipped -- clipping loses about 8% of the total width
## and biases the null low, which would make the observed counts look enriched.
CHR_END <- OFFSET + CHRLEN[CHROMS]
split_genome_iv <- function(a, b) {          # [a, b) in genome coordinates
  out <- vector("list", 0L)
  segs <- if (b <= GLEN) list(c(a, b)) else list(c(a, GLEN), c(0, b - GLEN))
  for (sg in segs) {
    s0 <- sg[1]; s1 <- sg[2]
    if (s1 <= s0) next
    i0 <- findInterval(s0, OFFSET); i1 <- findInterval(s1 - 1e-9, OFFSET)
    for (i in i0:i1) {
      lo <- max(s0, OFFSET[i]); hi <- min(s1, CHR_END[i])
      if (hi > lo) out[[length(out) + 1L]] <-
        list(chrom = CHROMS[i], lo = lo - OFFSET[i] + 1, hi = hi - OFFSET[i])
    }
  }
  out
}
shift_iv <- function(iv, off) {
  a <- (OFFSET[iv$chrom] + iv$lo - 1 + off) %% GLEN
  w <- iv$hi - iv$lo
  rbindlist(lapply(seq_along(a), function(i) rbindlist(split_genome_iv(a[i], a[i] + w[i]))))
}
shift_genes <- function(iv, off) genes_in(shift_iv(iv, off))

## Loci: overlapping intervals merged within a condition, the rule of
## Treusch et al. (2015), as in the landscape section.
merge_iv <- function(dt) {
  if (!nrow(dt)) return(data.table(chrom = character(), lo = numeric(), hi = numeric()))
  out <- dt[order(seqnames, start)][, {
    lo <- start[1]; hi <- end[1]; los <- c(); his <- c()
    for (k in seq_along(start)) {
      if (start[k] <= hi) { hi <- max(hi, end[k]) } else { los <- c(los, lo); his <- c(his, hi); lo <- start[k]; hi <- end[k] }
    }
    .(lo = c(los, lo), hi = c(his, hi))
  }, by = seqnames]
  setnames(out, "seqnames", "chrom")
  out[, `:=`(lo = pmax(lo, 1), hi = pmin(hi, CHRLEN[chrom]))]
  out[]
}
in_mask <- function(dt) {
  drop <- rep(FALSE, nrow(dt))
  for (i in seq_len(nrow(MASK))) {
    r <- MASK[i]
    same <- dt$seqnames == r$chrom
    inside <- same & dt$peak.pos >= r$lo & dt$peak.pos <= r$hi
    ov <- pmax(0, pmin(dt$end, r$hi) - pmax(dt$start, r$lo))
    drop <- drop | inside | (same & ov > SHOULDER_FRAC * (r$hi - r$lo))
  }
  drop
}

## ------------------------------------------------------- one peak set --------
run_one_set <- function(pset) {
  f <- file.path(PEAK_DIR, if (pset == "fitness") "fit_peaks.csv" else "stress_peaks.csv")
  if (!file.exists(f)) { say("  %s: %s absent, skipped", pset, basename(f)); return(NULL) }
  pk <- fread(f, showProgress = FALSE)
  pk <- pk[databin == "main" & peak.LOD >= LOD_CUT]
  if (!"sample" %in% names(pk)) pk[, sample := paste(cross, well, condition, sep = "_")]
  pk <- pk[!in_mask(pk)]
  say("")
  say("── peak set: %s ──", pset)
  say("  peaks    : %d retained in %d pools, %d crosses",
      nrow(pk), uniqueN(pk$sample), uniqueN(pk$cross))
  conds <- intersect(CONDS, unique(as.character(pk$condition)))
  loci <- rbindlist(lapply(conds, function(c) {
    m <- merge_iv(pk[condition == c]); if (nrow(m)) m[, condition := c]; m
  }), fill = TRUE)
  say("  loci     : %s", paste(sprintf("%s %d", LONG[conds],
      vapply(conds, function(c) sum(loci$condition == c), 0L)), collapse = ", "))

  results <- list()
  for (cd in conds) {
    iv <- loci[condition == cd]
    obs_genes <- genes_in(iv)
    n <- length(obs_genes); N <- length(universe)
    if (!n) { say("  %s: no loci, skipped", LONG[[cd]]); next }
    say("  %s: %d loci, %d genes under them", LONG[[cd]], nrow(iv), n)

    k <- vapply(go2genes, function(g) sum(obs_genes %in% g), 0L)
    K <- lengths(go2genes)
    p_hyper <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    fold <- (k / n) / (K / N)

    set.seed(1234)                       # seeded: a re-run reproduces
    hits <- matrix(0L, nrow = length(go2genes), ncol = NPERM)
    perm_n <- integer(NPERM)
    for (b in seq_len(NPERM)) {
      pg <- shift_genes(iv, sample.int(GLEN, 1))
      perm_n[b] <- length(pg)
      hits[, b] <- vapply(go2genes, function(g) sum(pg %in% g), 0L)
      if (b %% 200 == 0) say("      permutation %d / %d", b, NPERM)
    }
    p_perm <- (1 + rowSums(hits >= k)) / (NPERM + 1)
    exp_perm <- rowMeans(hits)

    res <- data.table(peak_set = pset, condition = cd, go = names(go2genes),
                      term = unname(term_name[names(go2genes)]),
                      genes_in_locus = k, term_size = K, genes_tested = n, universe = N,
                      fold = round(fold, 3), exp_perm = round(exp_perm, 2),
                      p_hyper = p_hyper, q_hyper = p.adjust(p_hyper, "BH"),
                      p_perm = p_perm, q_perm = p.adjust(p_perm, "BH"))
    res[, focus := ""]
    for (fo in names(FOCUS)) res[grepl(FOCUS[[fo]], term, ignore.case = TRUE), focus := fo]
    setorder(res, p_perm, p_hyper)
    fwrite(res, file.path(OUT, sprintf("go_%s_%s_%s.csv", ONT, pset, cd)))
    results[[cd]] <- res
    say("      median permuted gene count %d against %d observed",
        as.integer(median(perm_n)), n)
    top <- res[q_perm < 0.1][1:min(5, .N)]
    if (nrow(top) && !is.na(top$go[1]))
      for (i in seq_len(nrow(top)))
        say("      %-52s fold %.2f  p_perm %.3g", substr(top$term[i], 1, 52),
            top$fold[i], top$p_perm[i])
    else say("      nothing at q_perm < 0.1")
  }
  rbindlist(results)
}

all_res <- rbindlist(lapply(PEAK_SETS, run_one_set), fill = TRUE)
if (!nrow(all_res)) stop("no GO results produced", call. = FALSE)
fwrite(all_res, file.path(OUT, sprintf("go_%s_all.csv", ONT)))
focus <- all_res[focus != ""]
setorder(focus, peak_set, condition, p_perm)
fwrite(focus, file.path(OUT, sprintf("go_%s_focus_terms.csv", ONT)))

say("")
say("FIGURE PANELS — the question each was asked to answer.")
say("These are CONDITION-SPECIFIC claims, so read the `contrast` rows:")
for (q in list(c("Caf",  "cell wall",  "Figure 4D, caffeine"),
               c("NaCl", "salt / ion", "Figure 4D, sodium chloride"),
               c("Ace",  "acetate",    "Figure 5D, acetic acid"))) {
  for (ps in PEAK_SETS) {
    d <- focus[peak_set == ps & condition == q[1] & focus == q[2]]
    if (!nrow(d)) { say("  %-28s %-8s no terms matched the filter", q[3], ps); next }
    best <- d[which.min(p_perm)]
    say("  %-28s %-8s %d terms; best: %s, fold %.2f, p_perm %.3g, q_perm %.3g",
        q[3], ps, nrow(d), substr(best$term, 1, 40), best$fold, best$p_perm, best$q_perm)
  }
}
say("")
say("Wrote %s", OUT)
say("Quote p_perm, not p_hyper.  p_hyper assumes genes are exchangeable and they")
say("are not: the genes under one locus are physical neighbours.")

invisible(TRUE)
}

# ── 16. Downstream: genome-wide landscape ────────────────────────
st_land = NA_integer_
if (running('landscape'))
    st_land = run_inline('landscape',
                         function() qtl_landscape_run(PEAK_DIR = out_dir,
                                                      OUT_DIR  = land_dir))


# ── 17. Downstream: candidate genes under each locus ────────────────
# One run per variant dump. keep_chr7 comes straight from §0, so the
# manifest and the tables cannot disagree about it.
st_cand = list()
if (running('candidates')) {
    for (.tag in cand_tags) {
        vfile = file.path(out_dir, sprintf('qtl_variants_%s.csv', .tag))
        if (!file.exists(vfile)) {
            cat('!! ', vfile, " absent — run the 'variants' stage first;",
                ' skipping candidates_', .tag, '\n', sep = '')
            st_cand[[.tag]] = NA_integer_
            next
        }
        .tg = .tag
        st_cand[[.tg]] = run_inline(paste0('candidates_', .tg),
            function() qtl_candidates_run(TAG       = .tg,
                                          PEAK_DIR  = out_dir,
                                          OUT       = file.path(cand_dir, .tg),
                                          GFF       = gff_file,
                                          KEEP_CHR7 = keep_chr7))
    }
}

# ── 18. Downstream: GO-term enrichment ──────────────────────────────
# V7. Needs only the peak tables and the GFF, so unlike §17 it has no
# dependency on the variants stage. Which annotation ran goes into the
# manifest beside the exit code -- a GO result whose source nobody can
# identify is not quotable, and this is the first section in the project
# whose inputs live outside the run tree.
st_go = NA_integer_
if (running('go')) {
    .gomiss = if (identical(go_source, 'orgdb')) character(0) else
              c(gaf_file, obo_file)[!file.exists(path.expand(
                  c(gaf_file, obo_file)))]
    if (length(.gomiss)) {
        cat('!! GO annotation absent, skipping §18:\n   ',
            paste(.gomiss, collapse = '\n   '), '\n',
            '   GAF: http://sgd-archive.yeastgenome.org/curation/literature/gene_association.sgd.gz\n',
            '   OBO: http://purl.obolibrary.org/obo/go/go-basic.obo\n', sep = '')
    } else {
        st_go = run_inline('go',
            function() qtl_go_run(PEAK_DIR  = out_dir,
                                  OUT       = go_dir,
                                  GFF       = gff_file,
                                  KEEP_CHR7 = keep_chr7,
                                  GAF       = path.expand(gaf_file),
                                  OBO       = path.expand(obo_file),
                                  GO_SOURCE = go_source,
                                  ONT       = go_ont,
                                  NPERM     = go_perms,
                                  MIN_TERM  = go_min_term,
                                  MAX_TERM  = go_max_term,
                                  PEAK_SETS = go_peak_sets))
    }
}

# The manifest is re-saved with the downstream exit codes. §15 already
# wrote one, so a crash in §16-18 still leaves a manifest behind; this
# overwrites it with the complete record.
manifest$go = list(source = go_source, ont = go_ont, perms = go_perms,
                   peak_sets = go_peak_sets, gaf = gaf_file, obo = obo_file,
                   min_term = go_min_term, max_term = go_max_term)
manifest$downstream = list(landscape = st_land, candidates = st_cand, go = st_go,
                           land_dir = land_dir, cand_dir = cand_dir,
                           go_dir = go_dir)
manifest$run_end    = format(Sys.time())
qs2::qs_save(manifest, file.path(out_dir, 'run_manifest.qs'))

cat('\n══ run ', run_tag, ' complete ══\n', sep = '')
cat('everything this run produced is under:\n  ', run_dir, '\n', sep = '')
.bad = c(st_land, unlist(st_cand), st_go)
.bad = .bad[!is.na(.bad) & .bad != 0L]
if (length(.bad))
    cat('!! ', length(.bad), ' downstream section(s) exited non-zero — ',
        'the tables in landscape/, candidates/ and go/ are NOT all from this run.\n',
        sep = '')


# ── A. Optional interactive spot-checks ───────────────────────────────
if (FALSE) {
    sid = main_ids[1]
    plotFitnessTwoPanel(sid)

    k = names(drug_results)[vapply(drug_results, ok, logical(1))][1]
    plotPairFourPanel(k)

    # What is in the flagged bin, and why
    subset(key, databin == 'flagged',
           c(Sample_ID, Wells, cross, condition, af_median, med_depth,
             n_peaks_fit, flag_reasons))

    # Regression check against V4: main-databin fitness peaks should
    # match V4's for the same samples. A difference means the
    # restructure changed something it should not have.
    v4 = qs2::qs_read(file.path(cache_dir, 'peaks.qs'))$fit
    a = as.data.frame(fit_main)[, c('sample','seqnames','start','end')]
    b = as.data.frame(v4[mcols(v4)$sample %in% main_ids])[, c('sample','seqnames','start','end')]
    cat('V5 main:', nrow(a), ' V4 same samples:', nrow(b), ' identical:',
        isTRUE(all.equal(a[order(a$sample,a$seqnames,a$start),],
                         b[order(b$sample,b$seqnames,b$start),],
                         check.attributes = FALSE)), '\n')
}

# Discard a run entirely — tables, figures, cache and logs together,
# which is the only safe way to throw one away:
#   rm -r /u/project/kruglyak/bruni/projects/xQTLpaper/runs/v6
