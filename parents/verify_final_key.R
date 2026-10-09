# ══════════════════════════════════════════════════════════════════════
# verify_final_key.R — does every usable well still hold under Z_final?
# ══════════════════════════════════════════════════════════════════════
# The key's parent assignments were decided under Z_hybrid_maj. Z_final
# changes the genotype source for 9 disputed strains plus the 11 that
# agreed with the panel, which moves the segregating-site set for ~49 of
# the 90 usable wells. An assignment made on one site set is not
# automatically valid on another, so re-score every usable well against
# the parent the key now assigns it.
#
# This is a check, not a search: it does not propose alternatives, it
# asks whether what we already decided still stands.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
base      <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir <- file.path(base, 'notAF/count_files')
Z  <- qs2::qs_read(file.path(base, 'A1ref/files/Z_final.qs'))
k  <- read.csv(file.path(base, 'A1ref/files/flatkey_A1.csv'), check.names = FALSE)
out.file <- file.path(base, 'A1ref/qc/final_key_verification.csv')
tw <- function(x) trimws(as.character(x))
sid_v <- sub(':.*$', '', Z$variant_ids)
chr_v <- sub('_[0-9]+$', '', sid_v); pos_v <- as.integer(sub('^.*_', '', sid_v))
WIN <- 501
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
u <- unique(data.frame(W = tw(k$Wells), p = tw(k$parent_a),
                       a = tw(k$parent_alpha), un = k$unresolved,
                       stringsAsFactors = FALSE))
u <- u[!u$un, ]
cat('verifying ', nrow(u), ' usable wells under Z_final\n', sep = '')
res <- list()
for (i in seq_len(nrow(u))) {
    w <- u$W[i]; st <- u$p[i]; al <- u$a[i]
    if (!(st %in% rownames(Z$Z)) || !(al %in% rownames(Z$Z))) {
        cat('  ', w, ': parent missing from Z_final\n', sep = ''); next }
    sub <- k[tw(k$Wells) == w, , drop = FALSE]
    best <- NULL; bn <- -1
    for (sid in tw(sub$Sample_ID)) {
        f <- file.path(count_dir, paste0(sid, '.txt')); if (!file.exists(f)) next
        d <- try(fread(f, select = c('contig','position','refCount','altCount'),
                       showProgress = FALSE), silent = TRUE)
        if (inherits(d, 'try-error') || !nrow(d)) next
        j <- match(paste0(d$contig, '_', d$position), sid_v); ok <- !is.na(j)
        tt <- sum(d$refCount[ok] + d$altCount[ok])
        if (tt > bn) { bn <- tt; best <- list(sid=sid, j=j[ok],
                        ref=d$refCount[ok], alt=d$altCount[ok]) }
    }
    if (is.null(best)) next
    a <- as.numeric(Z$Z[st, ]); b <- as.numeric(Z$Z[al, ])
    seg <- a != b & a != 1 & b != 1; s <- seg[best$j]
    if (!any(s)) next
    j <- best$j[s]; r <- best$ref[s]; alt <- best$alt[s]
    tot <- r + alt; keep <- tot >= 3
    if (sum(keep) < 5*WIN) next
    jj <- j[keep]; tt <- tot[keep]
    p1 <- ifelse(a[jj] == 0, r[keep], alt[keep]); af <- p1/tt
    o <- order(chr_v[jj], pos_v[jj])
    dh <- mean(abs(af - 0.5)); od <- od_of(af[o], tt[o], chr_v[jj][o], pos_v[jj][o])
    res[[length(res)+1]] <- data.frame(Wells=w, parent=st, sample=best$sid,
        n_seg=sum(s), d_half=dh, od=od, depth=mean(tt))
    cat(sprintf('  %-4s %-4s  n_seg %6d  d_half %.3f  od %6.2f\n', w, st, sum(s), dh, od))
}
R <- do.call(rbind, res)
R$flag <- ifelse(is.na(R$od), 'no_od',
          ifelse(R$d_half > 0.20 | R$od > 2.5, 'FAIL',
          ifelse(R$d_half > 0.16 | R$od > 1.75, 'watch', 'clean')))
write.csv(R, out.file, row.names = FALSE)
cat('\n── verdict ──\n'); print(table(R$flag))
cat('\nmedian d_half:', round(median(R$d_half, na.rm=TRUE),3),
    ' median od:', round(median(R$od, na.rm=TRUE),2), '\n')
if (any(R$flag == 'FAIL')) { cat('\nFAIL wells:\n')
    print(R[R$flag=='FAIL', c('Wells','parent','n_seg','d_half','od')], row.names=FALSE, digits=3) }
cat('\nwrote ', out.file, '\n', sep = '')
