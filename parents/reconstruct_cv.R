# ══════════════════════════════════════════════════════════════════════
# reconstruct_cv.R — use the reconstructed parent DIRECTLY, honestly
# ══════════════════════════════════════════════════════════════════════
# reconstruct_parents.R derived each parent from its pool and then NAMED
# it against the 1011 panel. That only helps when the true parent is in
# the panel. For E03 and G10 the whole panel is flat -- 763 candidates
# within 0.005 of each other -- which means no named strain is right.
#
# The obvious next move is to skip naming and phase the well against its
# own reconstructed genotype. The trap is circularity: derive a parent
# from a pool, score that same pool against it, and d_half goes to ~0 BY
# CONSTRUCTION. Every well would look perfect, including a dead one.
#
# So: 4-fold cross-validation on the four conditions. Reconstruct the
# genotype from three samples, score the held-out fourth, rotate. The
# held-out sample never contributed to the genotype, so its d_half and
# od are honest.
#
#   real cross, unnamed parent -> reconstruction PREDICTS the held-out
#                                 sample; excess and od come out low
#   dead pool                  -> nothing to predict; held-out scores
#                                 stay bad no matter how it was built
#
# Known-good wells are included as positive controls: whatever CV score
# B01/H07/G11 produce is what "passing" looks like on this plate. A
# target is only recoverable if it reaches that range.
#
# d_half is depth-corrected throughout (excess = d_half - 0.399/sqrt(n)),
# since the binomial floor alone gives 0.16 at 6x and 0.07 at 30x, and
# these wells differ several-fold in depth.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
base      <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir <- file.path(base, 'notAF/count_files')
Z  <- qs2::qs_read(file.path(base, 'A1ref/files/Z_final.qs'))
k  <- read.csv(file.path(base, 'A1ref/files/flatkey_A1.csv'), check.names = FALSE)
out.file <- file.path(base, 'A1ref/qc/reconstruct_cv.csv')

TARGETS  <- c('A05','C02','E03','E12','G10','H01','H09')
CONTROLS <- c('B01','H07','G11','A02')
MIN_DEPTH <- 6; LO <- 0.15; HI <- 0.35; WIN <- 501

tw <- function(x) trimws(as.character(x))
sid_v <- sub(':.*$', '', Z$variant_ids); m <- length(sid_v)
chr_v <- sub('_[0-9]+$', '', sid_v); pos_v <- as.integer(sub('^.*_', '', sid_v))
od_of <- function(af, n, ch, po) {
    num <- c(); den <- c()
    for (cc in unique(ch)) {
        i <- which(ch == cc); i <- i[order(po[i])]
        if (length(i) < 2*WIN) next
        loc <- as.numeric(stats::filter(af[i], rep(1/WIN, WIN), sides = 2))
        ok <- !is.na(loc) & loc > 0 & loc < 1; if (!any(ok)) next
        r <- af[i][ok] - loc[ok]
        num <- c(num, r^2); den <- c(den, loc[ok]*(1-loc[ok])/n[i][ok])
    }
    if (!length(num)) NA_real_ else mean(num)/mean(den)
}
read1 <- function(sid) {
    f <- file.path(count_dir, paste0(sid, '.txt')); if (!file.exists(f)) return(NULL)
    d <- try(fread(f, select = c('contig','position','refCount','altCount'),
                   showProgress = FALSE), silent = TRUE)
    if (inherits(d, 'try-error') || !nrow(d)) return(NULL)
    j <- match(paste0(d$contig, '_', d$position), sid_v); ok <- !is.na(j)
    list(j = j[ok], ref = d$refCount[ok], alt = d$altCount[ok])
}

res <- list()
for (w in c(TARGETS, CONTROLS)) {
    sub <- k[tw(k$Wells) == w, , drop = FALSE]; if (!nrow(sub)) next
    alpha <- tw(sub$parent_alpha)[1]; assigned <- tw(sub$parent_a)[1]
    if (!alpha %in% rownames(Z$Z)) next
    b <- as.numeric(Z$Z[alpha, ])
    cn <- lapply(tw(sub$Sample_ID), read1); names(cn) <- tw(sub$Sample_ID)
    cn <- cn[!vapply(cn, is.null, logical(1))]
    if (length(cn) < 2) { cat('  ', w, ': <2 readable samples\n', sep=''); next }

    for (h in names(cn)) {
        # ---- reconstruct from every sample EXCEPT the held-out one ----
        ref <- numeric(m); alt <- numeric(m)
        for (s in setdiff(names(cn), h)) {
            ref[cn[[s]]$j] <- ref[cn[[s]]$j] + cn[[s]]$ref
            alt[cn[[s]]$j] <- alt[cn[[s]]$j] + cn[[s]]$alt
        }
        tot <- ref + alt; use <- tot >= MIN_DEPTH & b != 1
        if (!any(use)) next
        pa <- ifelse(b == 0, alt/pmax(tot,1), 1 - alt/pmax(tot,1))
        g <- rep(1L, m)
        d1 <- use & pa >= HI; s1 <- use & pa <= LO
        g[d1] <- ifelse(b[d1] == 0, 2L, 0L); g[s1] <- as.integer(b[s1])

        # ---- score the held-out sample against it ----
        seg <- g != b & g != 1 & b != 1
        x <- cn[[h]]; sel <- seg[x$j]
        if (!any(sel)) next
        jj <- x$j[sel]; r <- x$ref[sel]; a2 <- x$alt[sel]
        tt <- r + a2; keep <- tt >= 3
        if (sum(keep) < 5*WIN) next
        jj <- jj[keep]; tt <- tt[keep]
        p1 <- ifelse(g[jj] == 0, r[keep], a2[keep]); af <- p1/tt
        o <- order(chr_v[jj], pos_v[jj]); dp <- mean(tt)
        res[[length(res)+1]] <- data.frame(Wells = w, assigned = assigned,
            heldout = h, role = ifelse(w %in% CONTROLS, 'control', 'target'),
            n_seg = sum(sel), d_half = mean(abs(af-0.5)),
            excess = mean(abs(af-0.5)) - 0.399/sqrt(dp),
            od = od_of(af[o], tt[o], chr_v[jj][o], pos_v[jj][o]), depth = dp)
    }
    r <- do.call(rbind, res[vapply(res, function(z) z$Wells[1] == w, logical(1))])
    if (is.null(r)) { cat('  ', w, ': no fold scored\n', sep=''); next }
    cat(sprintf('  %-4s %-7s folds %d  n_seg %6.0f  excess %6.3f  od %5.2f\n',
        w, r$role[1], nrow(r), mean(r$n_seg), mean(r$excess, na.rm=TRUE),
        mean(r$od, na.rm=TRUE)))
}
R <- do.call(rbind, res); write.csv(R, out.file, row.names = FALSE)
cat('\n── mean over folds, by well ──\n')
agg <- aggregate(cbind(n_seg, excess, od) ~ Wells + role + assigned, R, mean)
print(agg[order(agg$role, agg$excess), ], row.names = FALSE, digits = 3)
cat('\ncontrol excess range:',
    paste(round(range(agg$excess[agg$role=='control']),3), collapse=' - '),
    '| control od range:',
    paste(round(range(agg$od[agg$role=='control']),2), collapse=' - '), '\n')
cat('\nwrote ', out.file, '\n', sep = '')
