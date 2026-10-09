# ══════════════════════════════════════════════════════════════════════
# scan_919.R — widen the candidate set to the whole 1011 panel
# ══════════════════════════════════════════════════════════════════════
# scan_all_candidates.R only tried the 47 strains used in this
# experiment. If a well received a tube from elsewhere in the
# collection, no scan so far could have found it. This tries all 919.
#
# Candidate genotypes are NOT cached: 919 rows x 1.24M doubles would be
# ~9 GB. Each row is materialised, used, and dropped.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages(library(Matrix))
base      <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir <- file.path(base, 'notAF/count_files')
panel     <- '/Users/kruglyakbws/Desktop/hoffmanDrop/Z_1011_filtered.qs'
key.file  <- file.path(base, 'A1ref/files/flatkey_A1.csv')
unr.file  <- file.path(base, 'A1ref/qc/unresolved_wells_A1.txt')
out.file  <- file.path(base, 'A1ref/qc/scan_919.csv')
MIN_DEPTH <- 3; WIN <- 501

P <- qs2::qs_read(panel); zr <- rownames(P$Z)
sid_v <- sub(':.*$', '', P$variant_ids)
chr_v <- sub('_[0-9]+$', '', sid_v); pos_v <- as.integer(sub('^.*_', '', sid_v))
tw <- function(x) trimws(as.character(x))
k <- read.csv(key.file, check.names = FALSE)
wells <- tw(readLines(unr.file)); wells <- wells[nzchar(wells)]
cat('scanning ', length(wells), ' wells x ', length(zr), ' panel strains\n', sep='')

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

res <- list()
for (w in wells) {
    sub <- k[tw(k$Wells) == w, , drop = FALSE]
    alpha <- tw(sub$parent_alpha)[1]
    if (!alpha %in% zr) { cat('  ', w, ': alpha not in panel\n', sep=''); next }
    best <- NULL; bn <- -1
    for (sid in tw(sub$Sample_ID)) {
        f <- file.path(count_dir, paste0(sid, '.txt')); if (!file.exists(f)) next
        d <- try(read.delim(f, stringsAsFactors = FALSE), silent = TRUE)
        if (inherits(d, 'try-error') || !nrow(d)) next
        j <- match(paste0(d$contig,'_',d$position), sid_v); ok <- !is.na(j)
        tt <- sum(d$refCount[ok] + d$altCount[ok])
        if (tt > bn) { bn <- tt; best <- list(sid=sid, j=j[ok],
                        ref=d$refCount[ok], alt=d$altCount[ok]) }
    }
    if (is.null(best)) { cat('  ', w, ': no counts\n', sep=''); next }
    b <- as.numeric(P$Z[alpha, ])
    rw <- list()
    for (cand in zr) {
        if (cand == alpha) next
        a <- as.numeric(P$Z[cand, ])
        seg <- a != b & a != 1 & b != 1
        s <- seg[best$j]; if (!any(s)) next
        j <- best$j[s]; r <- best$ref[s]; al <- best$alt[s]
        tot <- r + al; keep <- tot >= MIN_DEPTH
        if (sum(keep) < 5*WIN) next
        jj <- j[keep]; tt <- tot[keep]
        p1 <- ifelse(a[jj] == 0, r[keep], al[keep])
        af <- p1/tt; o <- order(chr_v[jj], pos_v[jj])
        rw[[length(rw)+1]] <- data.frame(Wells=w, sample=best$sid, cand=cand,
            n_seg=sum(s), d_half=mean(abs(af-0.5)),
            od=od_of(af[o], tt[o], chr_v[jj][o], pos_v[jj][o]))
    }
    Rw <- do.call(rbind, rw)
    if (is.null(Rw) || !nrow(Rw)) {
        # No candidate in 919 yields a scoreable site set: not a mislabel,
        # the well has nothing to phase against anything.
        cat('  ', w, ': NO scoreable candidate in the whole panel\n', sep = '')
        next
    }
    res[[w]] <- Rw
    q <- Rw[!is.na(Rw$od) & Rw$n_seg >= 15000, ]; q <- q[order(q$d_half), ]
    cat(sprintf('  %-4s  top: %s\n', w, if (nrow(q))
        paste(sprintf('%s %.3f/od %.2f (%d)', q$cand[1], q$d_half[1], q$od[1], q$n_seg[1]))
        else 'none with >=15k sites'))
}
A <- do.call(rbind, res)
if (is.null(A)) stop('no results at all')
write.csv(A, out.file, row.names = FALSE)
cat('\nwrote ', out.file, '\n', sep='')
