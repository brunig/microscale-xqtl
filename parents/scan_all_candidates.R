# ══════════════════════════════════════════════════════════════════════
# scan_all_candidates.R — score every parent against the unresolved wells
# ══════════════════════════════════════════════════════════════════════
# test_parents_A1.R only ever compared TWO candidates per well: the
# nominal and whatever round 3 chose. For the 15 unresolved wells that
# established the assigned parent is wrong; it did not establish that
# no parent fits. This scores all 47 rows of Z_hybrid_maj against each
# well's deepest sample, so a well is only declared unrecoverable after
# every available strain has actually been tried.
#
# F09/G08 matter most here. Their nominal is YBB, which has no genotype
# anywhere, so they have never been scored against ANY candidate -- not
# once in the whole override chain. If some strain fits them, those 8
# samples come back from nothing.
#
# Same two statistics as before, and the same asymmetry: many
# segregating sites plus AF far off 0.5 means wrong. A candidate is
# only a real answer if it has BOTH a competitive d_half AND enough
# sites to mean something -- a thin site set makes any candidate look
# tidy, which is the trap that made ACQ appear to beat ADF at 4,950
# sites vs 53,538.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages(library(Matrix))

base      <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir <- file.path(base, 'notAF/count_files')
z.file    <- file.path(base, 'A1ref/files/Z_hybrid_maj.qs')
key.file  <- file.path(base, 'A1ref/files/flatkey_A1.csv')
unr.file  <- file.path(base, 'A1ref/qc/unresolved_wells_A1.txt')
out.file  <- file.path(base, 'A1ref/qc/candidate_scan.csv')
MIN_DEPTH <- 3; WIN <- 501

Z <- qs2::qs_read(z.file); zr <- rownames(Z$Z)
stopifnot(length(zr) > 40)
sid_v <- sub(':.*$', '', Z$variant_ids)
chr_v <- sub('_[0-9]+$', '', sid_v)
pos_v <- as.integer(sub('^.*_', '', sid_v))

tw <- function(x) trimws(as.character(x))
k  <- read.csv(key.file, check.names = FALSE)
wells <- tw(readLines(unr.file)); wells <- wells[nzchar(wells)]
cat('scanning ', length(wells), ' wells x ', length(zr), ' candidates\n', sep = '')

gcache <- new.env(parent = emptyenv())
geno <- function(s) { if (is.null(gcache[[s]])) gcache[[s]] <- as.numeric(Z$Z[s, ]); gcache[[s]] }

od_of <- function(af, n, ch, po) {
    num <- c(); den <- c()
    for (cc in unique(ch)) {
        i <- which(ch == cc); i <- i[order(po[i])]
        if (length(i) < 2 * WIN) next
        loc <- as.numeric(stats::filter(af[i], rep(1/WIN, WIN), sides = 2))
        ok <- !is.na(loc) & loc > 0 & loc < 1
        if (!any(ok)) next
        r <- af[i][ok] - loc[ok]
        num <- c(num, r^2); den <- c(den, loc[ok]*(1-loc[ok])/n[i][ok])
    }
    if (!length(num)) NA_real_ else mean(num)/mean(den)
}

res <- list()
for (w in wells) {
    sub <- k[tw(k$Wells) == w, , drop = FALSE]
    alpha <- tw(sub$parent_alpha)[1]
    best <- NULL; bn <- -1
    for (sid in tw(sub$Sample_ID)) {
        f <- file.path(count_dir, paste0(sid, '.txt'))
        if (!file.exists(f)) next
        d <- try(read.delim(f, stringsAsFactors = FALSE), silent = TRUE)
        if (inherits(d, 'try-error') || !nrow(d)) next
        j <- match(paste0(d$contig, '_', d$position), sid_v); ok <- !is.na(j)
        tt <- sum(d$refCount[ok] + d$altCount[ok])
        if (tt > bn) { bn <- tt; best <- list(sid = sid, j = j[ok],
                        ref = d$refCount[ok], alt = d$altCount[ok]) }
    }
    if (is.null(best)) { cat('  ', w, ': no readable counts\n', sep=''); next }
    b <- geno(alpha)
    for (cand in zr) {
        if (cand == alpha) next
        a <- geno(cand)
        seg <- a != b & a != 1 & b != 1
        s <- seg[best$j]
        if (!any(s)) next
        j <- best$j[s]; r <- best$ref[s]; al <- best$alt[s]
        tot <- r + al; keep <- tot >= MIN_DEPTH
        n_seg <- sum(s)
        if (sum(keep) < 5*WIN) {
            res[[length(res)+1]] <- data.frame(Wells=w, sample=best$sid,
                cand=cand, n_seg=n_seg, n_used=sum(keep), d_half=NA, od=NA, depth=NA)
            next
        }
        j <- j[keep]; tot <- tot[keep]
        p1 <- ifelse(a[j] == 0, r[keep], al[keep])
        af <- p1/tot; o <- order(chr_v[j], pos_v[j])
        res[[length(res)+1]] <- data.frame(Wells=w, sample=best$sid, cand=cand,
            n_seg=n_seg, n_used=sum(keep), d_half=mean(abs(af-0.5)),
            od=od_of(af[o], tot[o], chr_v[j][o], pos_v[j][o]), depth=mean(tot))
    }
    R <- do.call(rbind, res[sapply(res, function(x) x$Wells[1] == w)])
    R <- R[!is.na(R$d_half) & R$n_seg >= 10000, ]
    R <- R[order(R$d_half), ]
    cur <- tw(sub$parent_a)[1]
    cat(sprintf('  %-4s (assigned %-4s)  best: %s\n', w, cur,
        if (nrow(R)) paste(sprintf('%s %.3f/od %.2f (%d sites)', R$cand[1],
            R$d_half[1], R$od[1], R$n_seg[1])) else 'nothing scoreable'))
}
A <- do.call(rbind, res)
write.csv(A, out.file, row.names = FALSE)
cat('\nwrote ', out.file, '\n', sep = '')
