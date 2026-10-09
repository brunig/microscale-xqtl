# ══════════════════════════════════════════════════════════════════════
# build_final_key.R — the keyfile, decided on sequenced parents
# ══════════════════════════════════════════════════════════════════════
# Supersedes the override1/2/3 chain. Those rounds picked each well's
# parent_a by nearest match against the 1011 PANEL; this one judges the
# candidates against the physically sequenced A1 haploids, using the
# well's own reads (test_parents_A1.R).
#
# The evidence is ASYMMETRIC and the rules below respect that. A
# candidate yielding MANY segregating sites while the pool sits far off
# 0.5 is definitively wrong -- a wrong parent invents segregating sites.
# The converse does not hold: few sites plus a tidy d_half means
# uninformative, never "confirmed". So the test REJECTS parents; it
# does not confirm them, and a well with nothing rejected keeps
# whatever the od rounds already chose rather than being churned.
#
# Precedence, highest first:
#   1. A1 read test rejects one candidate  -> use the other
#   2. A1 read test rejects both           -> unresolved, exclude
#   3. A1 test undecided                   -> keep the round-3 parent
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages(library(Matrix))   # rownames() on the
# sparse Z needs Matrix ATTACHED, not merely loaded: without it the call
# returns NULL, `%in% NULL` is all-FALSE, and every well silently reads
# as "no genotype". Guarded below so it can never fail quietly again.
base   <- '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
o2.f   <- file.path(base, 'notAF/files/flatkey_override2.csv')
nom.f  <- file.path(base, 'notAF/files/flatkey.csv')
od.f   <- file.path(base, 'notAF/partb.tsv')            # od table, round 3 input
a1.f   <- file.path(base, 'A1ref/qc/parent_test_A1_maj.csv')
out.f  <- file.path(base, 'A1ref/files/flatkey_A1.csv')
log.f  <- file.path(base, 'A1ref/qc/key_changes.csv')
unr.f  <- file.path(base, 'A1ref/qc/unresolved_wells_A1.txt')

OD_APPLY <- 1.45; OD_EXCLUDE <- 2.50; MIN_MARGIN <- 0.02
SEG_MANY <- 20000; DHALF_BAD <- 0.20    # rejection thresholds

tw <- function(x) trimws(as.character(x))
k  <- read.csv(o2.f, check.names = FALSE)
nk <- read.csv(nom.f, check.names = FALSE)
k$Wells <- tw(k$Wells)
k$parent_nominal <- tw(nk$parent_a)[match(tw(k$Sample_ID), tw(nk$Sample_ID))]

# ── Reconstruct round 3 (verified against the 2026-08-26 Hoffman run:
# 10 wells changed, same wells) so the diff has a real baseline. ───────
R <- read.table(od.f, header = TRUE, stringsAsFactors = FALSE, na.strings = 'NA')
pick <- function(nom, bst, odn, odb) {
    if (is.na(odb) && is.na(odn)) return(c(nom, 'hold_no_information', NA))
    if (is.na(odb))               return(c(nom, 'keep_nominal', odn))
    better <- is.na(odn) || odb < odn * (1 - MIN_MARGIN)
    if (!better) return(c(nom, 'keep_nominal', if (is.na(odn)) odb else odn))
    if (identical(tw(nom), tw(bst))) return(c(nom, 'confirmed', odb))
    if (odb <= OD_APPLY) return(c(bst, 'apply_best', odb))
    c(bst, 'keep_best_flagged', odb)
}
d <- t(mapply(pick, R$nominal, R$best, R$od_nominal, R$od_best))
R$r3_parent <- d[, 1]; R$r3_decision <- d[, 2]
R$od_final <- suppressWarnings(as.numeric(d[, 3]))
i <- match(k$Wells, R$Wells)
k$parent_r3 <- ifelse(is.na(i), tw(k$parent_a), R$r3_parent[i])
k$od_final  <- ifelse(is.na(i), NA, R$od_final[i])

# ── Apply the A1 read test ────────────────────────────────────────────
A <- read.csv(a1.f, stringsAsFactors = FALSE)
rej <- function(n, dh) !is.na(n) & !is.na(dh) & n > SEG_MANY & dh > DHALF_BAD
A$nom_rej <- rej(A$n_seg_nom, A$dhalf_nom)
A$cur_rej <- rej(A$n_seg_cur, A$dhalf_cur)
# Catastrophic d_half counts as rejection even on a thin site set: G10
# sits at 0.47 on both candidates with only 7,993 sites on the nominal.
A$nom_rej <- A$nom_rej | (!is.na(A$dhalf_nom) & A$dhalf_nom > 0.40)
A$cur_rej <- A$cur_rej | (!is.na(A$dhalf_cur) & A$dhalf_cur > 0.40)

j <- match(k$Wells, tw(A$Wells))
k$a1_nom <- ifelse(is.na(j), NA, A$nominal[j])
k$a1_cur <- ifelse(is.na(j), NA, A$current[j])
k$parent_final <- k$parent_r3
k$decision <- 'r3_kept (A1 undecided)'
k$unresolved <- FALSE

for (r in seq_len(nrow(k))) {
    jj <- j[r]; if (is.na(jj)) { k$decision[r] <- 'no_A1_data'; k$unresolved[r] <- TRUE; next }
    nr <- A$nom_rej[jj]; cr <- A$cur_rej[jj]
    if (nr && cr) { k$decision[r] <- 'A1: both rejected'; k$unresolved[r] <- TRUE }
    else if (cr)  { k$parent_final[r] <- A$nominal[jj]; k$decision[r] <- 'A1: current rejected -> nominal' }
    else if (nr)  { k$parent_final[r] <- A$current[jj]; k$decision[r] <- 'A1: nominal rejected -> current' }
}
# Wells the od rounds could not resolve stay excluded unless A1 settled them.
k$unresolved[!is.na(k$od_final) & k$od_final > OD_EXCLUDE &
             !grepl('^A1: (current|nominal) rejected', k$decision)] <- TRUE

# Grounds carried forward from the od rounds that the A1 test cannot
# see, and which would otherwise be silently dropped:
#  - the chosen parent has no genotype anywhere (YBB: absent from the
#    1011 panel AND its A1 library is empty, so F09/G08 can never be
#    phased). Without this they read as "A1 undecided" and slip through.
#  - Part A's NOT_CLONAL verdict on A05: the candidate cluster is not
#    near-clonal, so the WELL is suspect (mixed or failed pool), which
#    is a stronger finding than a mislabel and is invisible to a test
#    that only compares two parents.
zc <- rownames(qs2::qs_read(file.path(base, 'A1ref/files/Z_hybrid_maj.qs'))$Z)
stopifnot(length(zc) > 40, 'BBQ' %in% zc)
no_geno <- !(tw(k$parent_final) %in% zc) | !(tw(k$parent_alpha) %in% zc)
if (any(no_geno)) {
    k$unresolved[no_geno] <- TRUE
    k$decision[no_geno] <- 'no genotype for chosen parent'
}
NOT_CLONAL <- c('A05')
k$unresolved[k$Wells %in% NOT_CLONAL] <- TRUE
k$decision[k$Wells %in% NOT_CLONAL] <- 'PartA: well not clonal (mixed/failed pool)'

k$parent_a <- k$parent_final
k$cross <- paste(tw(k$parent_a), tw(k$parent_alpha), sep = 'x')
write.csv(k, out.f, row.names = FALSE)

w <- unique(k[, c('Wells','parent_nominal','parent_r3','parent_final','decision','od_final','unresolved')])
w <- w[order(w$Wells), ]
write.csv(w, log.f, row.names = FALSE)
ex <- sort(unique(k$Wells[k$unresolved]))
writeLines(ex, unr.f)

cat('── vs ORIGINAL flatkey.csv ──\n')
cat(sum(tw(w$parent_final) != tw(w$parent_nominal)), 'of', nrow(w), 'wells differ\n')
cat('── vs ROUND 3 ──\n')
ch <- w[tw(w$parent_final) != tw(w$parent_r3), ]
cat(nrow(ch), 'wells differ\n'); if (nrow(ch)) print(ch, row.names = FALSE)
cat('\n── decisions ──\n'); print(table(w$decision))
cat('\nunresolved:', length(ex), 'wells /', 4 * length(ex), 'samples\n  ',
    paste(ex, collapse = ', '), '\n')
cat('\nwrote', out.f, '\n     ', log.f, '\n     ', unr.f, '\n')
