# ══════════════════════════════════════════════════════════════════════
# test_row_source.R — A1-sequenced row or 1011 panel row, per strain
# ══════════════════════════════════════════════════════════════════════
# The hybrid currently picks a genotype source by "did the A1 library
# work". That is the wrong criterion, and we have counterexamples in
# both directions:
#   AAM  98.52% concordant -- the A1 row was RIGHT. The panel row is
#        what broke A01/A02/A12; the A1 row rescued all three.
#   ACV  94.81% concordant -- the A1 row was WRONG. E08 recovered onto
#        panel-derived ACT while A1-ACV was rejected at d_half 0.411.
# So a low concordance says "these are different strains" without
# saying which one is this well's parent. Only the reads can.
#
# Strains at >=99.9% need no test: same genotype, and the panel row
# simply has more callable sites because het-suppression ate the A1
# ones (ABS 8,672 -> 15,600; ABH 20,769 -> 35,980). Those are free.
#
# This scores the 9 DISPUTED strains both ways on every well that uses
# them, then reports a per-strain verdict. Wells of one strain should
# agree with each other -- if they split, the strain is not the issue,
# that well is.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
base      <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir <- file.path(base, 'notAF/count_files')
Zh <- qs2::qs_read(file.path(base, 'A1ref/files/Z_hybrid_maj.qs'))
P  <- qs2::qs_read('/Users/kruglyakbws/Desktop/hoffmanDrop/Z_1011_filtered.qs')
g  <- read.csv(file.path(base, 'A1ref/qc/parent_genotyping_maj.csv'), stringsAsFactors = FALSE)
pv <- read.csv(file.path(base, 'A1ref/qc/Z_hybrid_provenance.csv'), stringsAsFactors = FALSE)
k  <- read.csv(file.path(base, 'A1ref/files/flatkey_A1.csv'), check.names = FALSE)
out.file <- file.path(base, 'A1ref/qc/row_source_test.csv')

tw <- function(x) trimws(as.character(x))
sid_v <- sub(':.*$', '', Zh$variant_ids)
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

seqs <- pv$strain[pv$source == 'A1_sequenced']
disp <- intersect(seqs[ g$panel_concord[match(seqs, g$strain)] < 99.9 ], rownames(P$Z))
disp <- disp[!is.na(disp)]
use <- unique(data.frame(W = tw(k$Wells), p = tw(k$parent_a), u = k$unresolved,
                         stringsAsFactors = FALSE))
use <- use[!use$u & use$p %in% disp, ]
cat('disputed strains: ', paste(disp, collapse = ', '), '\n',
    'wells to test  : ', nrow(use), '\n\n', sep = '')

bq <- as.numeric(Zh$Z['BBQ', ])
res <- list()
for (i in seq_len(nrow(use))) {
    w <- use$W[i]; st <- use$p[i]
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
    row <- data.frame(Wells=w, strain=st, sample=best$sid,
                      concord=g$panel_concord[match(st, g$strain)])
    for (src in c('A1','panel')) {
        a <- if (src == 'A1') as.numeric(Zh$Z[st, ]) else as.numeric(P$Z[st, ])
        seg <- a != bq & a != 1 & bq != 1; s <- seg[best$j]
        if (!any(s)) { row[[paste0('n_',src)]] <- 0; row[[paste0('d_',src)]] <- NA
                       row[[paste0('od_',src)]] <- NA; next }
        j <- best$j[s]; r <- best$ref[s]; al <- best$alt[s]
        tot <- r + al; keep <- tot >= 3
        if (sum(keep) < 5*WIN) { row[[paste0('n_',src)]] <- sum(s)
            row[[paste0('d_',src)]] <- NA; row[[paste0('od_',src)]] <- NA; next }
        jj <- j[keep]; tt <- tot[keep]
        p1 <- ifelse(a[jj] == 0, r[keep], al[keep]); af <- p1/tt
        o <- order(chr_v[jj], pos_v[jj])
        row[[paste0('n_',src)]]  <- sum(s)
        row[[paste0('d_',src)]]  <- mean(abs(af-0.5))
        row[[paste0('od_',src)]] <- od_of(af[o], tt[o], chr_v[jj][o], pos_v[jj][o])
    }
    res[[length(res)+1]] <- row
    cat(sprintf('  %-4s %-4s  A1: %6d sites d %.3f od %5.2f | panel: %6d sites d %.3f od %5.2f\n',
        w, st, row$n_A1, row$d_A1, row$od_A1, row$n_panel, row$d_panel, row$od_panel))
}
R <- do.call(rbind, res)
# Winner needs BOTH statistics, not just the tidier mean: a thin site
# set flatters d_half regardless of correctness.
R$winner <- ifelse(is.na(R$d_panel), 'A1',
            ifelse(is.na(R$d_A1), 'panel',
            ifelse(R$d_panel < R$d_A1 - 0.01, 'panel',
            ifelse(R$d_A1 < R$d_panel - 0.01, 'A1', 'tie'))))
write.csv(R, out.file, row.names = FALSE)
cat('\n── per-well winner ──\n'); print(table(R$winner))
cat('\n── per-strain ──\n')
print(as.data.frame.matrix(table(R$strain, R$winner)))
cat('\nwrote ', out.file, '\n', sep = '')
