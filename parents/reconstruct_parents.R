# ══════════════════════════════════════════════════════════════════════
# reconstruct_parents.R — derive parent_a's genotype from the pool itself
# ══════════════════════════════════════════════════════════════════════
# Every method so far asked "which reference strain best explains these
# reads". That fails whenever the true parent is in no reference — YBB
# is absent from the 1011 panel AND its A1 library is empty, so F09/G08
# have never been scored against anything that could be right.
#
# But the pool IS the parent, half the time. An xQTL pool is
# parent_a x BBQ, so at a site where BBQ carries reference, every alt
# read came from parent_a. With BBQ's genotype known, parent_a's
# genotype falls out of the pool's allele frequencies — no reference
# strain needed:
#
#   BBQ = ref  and pool alt fraction ~0.5  ->  parent_a is ALT   (2)
#   BBQ = ref  and pool alt fraction ~0    ->  parent_a is REF   (0)
#   BBQ = alt  -> uninformative here (mirrored logic, handled too)
#
# Four samples per well are summed first: the conditions share a parent,
# so this is four independent looks at one genotype.
#
# ── What the shape of the histogram tells you ────────────────────────
# The alt-fraction distribution at BBQ-ref sites is the diagnostic that
# separates the two failure modes we currently cannot tell apart:
#   BIMODAL near 0 and 0.5  -> a real cross. The 0.5 mode IS parent_a's
#                              alt sites, so the parent is recoverable
#                              even if it matches no known strain.
#   UNIMODAL near 0         -> the pool collapsed to BBQ (or the well
#                              never crossed). Nothing to recover; this
#                              is a dead pool, not a mislabel.
# `frac_het_like` below is that mode's mass. High = recoverable.
#
# ── The honest caveat ────────────────────────────────────────────────
# Selection distorts pool frequencies, so a region under strong
# selection reconstructs badly and the genome-wide call is conservative
# there. A CALL_RATE well under the cohort is itself informative rather
# than a failure.
#
# Also emits, for every well including the ones that already work, the
# agreement between the reconstructed parent and the ASSIGNED one. That
# is the only available check on the 14 strains whose genotype still
# comes from the 1011 panel because their A1 library failed — ABA and
# ABB serve 5 wells each on nothing but panel rows. If five ABA wells
# reconstruct to the same thing and it disagrees with panel-ABA, that is
# the next AAM.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(Matrix); library(data.table) })

base      <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir <- file.path(base, 'notAF/count_files')
z.file    <- file.path(base, 'A1ref/files/Z_hybrid_maj.qs')
panel     <- '/Users/kruglyakbws/Desktop/hoffmanDrop/Z_1011_filtered.qs'
key.file  <- file.path(base, 'A1ref/files/flatkey_A1.csv')
out.qs    <- file.path(base, 'A1ref/files/Z_reconstructed.qs')
rep.file  <- file.path(base, 'A1ref/qc/reconstruction.csv')
match.file<- file.path(base, 'A1ref/qc/reconstruction_matches.csv')

MIN_DEPTH  <- 8      # summed over the well's samples, per site
LO         <- 0.15   # alt frac below this -> parent is ref
HI         <- 0.35   # alt frac above this -> parent carries alt
ONLY_WELLS <- c('F09','G08','A05','G10','C02','C12','E12','H01','H09','B01','B11','C11','G09','C01','B02','B12','E07','H07','H08')  # 9 unresolved + ABA/ABB wells to validate the panel rows

Z <- qs2::qs_read(z.file); zr <- rownames(Z$Z)
sid_v <- sub(':.*$', '', Z$variant_ids); m <- length(sid_v)
tw <- function(x) trimws(as.character(x))
k <- read.csv(key.file, check.names = FALSE)
wells <- unique(tw(k$Wells)); if (!is.null(ONLY_WELLS)) wells <- ONLY_WELLS
cat('reconstructing ', length(wells), ' wells\n', sep = '')

i_all <- integer(0); j_all <- integer(0); x_all <- numeric(0)
rep <- list(); rn <- character(0)

for (wi in seq_along(wells)) {
    w <- wells[wi]
    sub <- k[tw(k$Wells) == w, , drop = FALSE]
    alpha <- tw(sub$parent_alpha)[1]
    if (!alpha %in% zr) next
    b <- as.numeric(Z$Z[alpha, ])

    # Sum the four conditions: same parent, four independent looks.
    ref <- numeric(m); alt <- numeric(m)
    for (sid in tw(sub$Sample_ID)) {
        f <- file.path(count_dir, paste0(sid, '.txt')); if (!file.exists(f)) next
        d <- try(fread(f, select = c('contig','position','refCount','altCount'),
                       showProgress = FALSE), silent = TRUE)
        if (inherits(d, 'try-error') || !nrow(d)) next
        j <- match(paste0(d$contig, '_', d$position), sid_v); ok <- !is.na(j)
        ref[j[ok]] <- ref[j[ok]] + d$refCount[ok]
        alt[j[ok]] <- alt[j[ok]] + d$altCount[ok]
    }
    tot <- ref + alt
    use <- tot >= MIN_DEPTH & b != 1          # BBQ must be confidently called
    if (!any(use)) next
    af <- alt / pmax(tot, 1)

    # At BBQ-ref sites the alt reads are parent_a's. At BBQ-alt sites the
    # logic mirrors: it is the REF reads that are parent_a's.
    pa <- ifelse(b == 0, af, 1 - af)   # frequency of the NON-BBQ allele
    differs <- use & pa >= HI          # parent_a carries the other allele
    matches <- use & pa <= LO          # parent_a carries what BBQ carries
    g <- rep(1L, m)                    # 1 = no-call, dropped downstream
    g[differs] <- ifelse(b[differs] == 0, 2L, 0L)
    g[matches] <- as.integer(b[matches])

    n_use <- sum(use)
    het_like <- sum(use & pa > LO & pa < HI)
    nm <- paste0('POOL_', w)
    nz <- g != 0L
    i_all <- c(i_all, rep_len(wi, sum(nz))); j_all <- c(j_all, which(nz))
    x_all <- c(x_all, g[nz]); rn <- c(rn, nm)

    rep[[length(rep)+1]] <- data.frame(Wells = w, assigned = tw(sub$parent_a)[1],
        n_sites_used = n_use, call_rate = round(sum(use & g != 1L)/max(1,n_use), 3),
        n_differs = sum(differs),
        frac_het_like = round(het_like / max(1, n_use), 3),
        median_depth = median(tot[use]),
        verdict = ifelse(het_like / max(1, n_use) < 0.02,
                         'DEAD POOL (no 0.5 mode)', 'recoverable'))
    cat(sprintf('  %-4s assigned %-4s  sites %7d  het-like %.3f  %s\n',
        w, tw(sub$parent_a)[1], n_use, het_like/max(1,n_use),
        rep[[length(rep)]]$verdict))
}

R <- do.call(rbind, rep)
Zr <- sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(length(rn), m),
                   dimnames = list(rn, NULL))
qs2::qs_save(list(Z = Zr, variant_ids = Z$variant_ids), out.qs)

# ── Name the reconstruction: nearest panel strain over confident sites ──
P <- qs2::qs_read(panel)
best <- list()
for (i in seq_along(rn)) {
    a <- as.numeric(Zr[i, ]); ok <- which(a != 1)
    if (length(ok) < 5000) next
    if (length(ok) > 20000) ok <- sort(sample(ok, 20000))
    S <- as.matrix(P$Z[, ok, drop = FALSE])
    av <- a[ok]
    d <- rowSums(S != rep(av, each = nrow(S)) & S != 1)
    o <- order(d)
    best[[length(best)+1]] <- data.frame(Wells = R$Wells[i],
        assigned = R$assigned[i], n_cmp = length(ok),
        top1 = rownames(P$Z)[o[1]], d1 = d[o[1]],
        top2 = rownames(P$Z)[o[2]], d2 = d[o[2]],
        d_assigned = if (R$assigned[i] %in% rownames(P$Z))
                        d[R$assigned[i]] else NA_real_)
}
B <- do.call(rbind, best)
R$top_match <- B$top1[match(R$Wells, B$Wells)]
R$d_top     <- B$d1[match(R$Wells, B$Wells)]
R$d_assigned<- B$d_assigned[match(R$Wells, B$Wells)]
write.csv(R, rep.file, row.names = FALSE)
write.csv(B, match.file, row.names = FALSE)
cat('\nwrote ', out.qs, '\n      ', rep.file, '\n      ', match.file, '\n', sep='')
