# ══════════════════════════════════════════════════════════════════════
# verify_parents.R — d_Bref screen + sibling check on the current key
# ══════════════════════════════════════════════════════════════════════
# Needs nothing but fit_results.qs and the keyfile. No Z, no scan.
#
# d_Bref / d_Balt: depth-weighted parent_a frequency minus 0.5, split by
#   whether the fixed parent (BBQ) carries ref or alt at each site. The
#   phasing encodes this directly: p2 == ref means BBQ is ref there.
#   Clean-cohort floor measured 2026-08-25: d_Bref ~ -0.022, d_Balt ~ 0.
#
# bref_frac is NOT a misassignment flag — it is a property of the strain
#   pair and is near-constant within a parent group. Only bref_frac ~ 0
#   is diagnostic: it means a reference-like strain got assigned.
# ──────────────────────────────────────────────────────────────────────
library(qs2)
base    <- '/u/project/kruglyak/bruni/projects/xQTLpaper'
keyfile <- file.path(base, 'files/flatkey_override2.csv')   # <- current key
COND    <- 'Glu'
FLOOR   <- -0.05    # at/near the clean floor
BAD     <- -0.15    # clearly wrong
DIVERGE <-  0.05    # how far from clean siblings counts as divergent

f <- qs_read(file.path(base, 'cache/fit_results.qs'))
k <- read.csv(keyfile, check.names = FALSE)
g <- k[k$Drug3Let == COND, ]

r <- do.call(rbind, lapply(intersect(g$Sample_ID, names(f)), function(sid) {
    d <- f[[sid]]; if (!is.data.frame(d)) return(NULL)
    q <- function(p) d[[paste0(p, '_', sid)]]
    ref <- q('ref'); alt <- q('alt'); p1 <- q('p1'); p2 <- q('p2')
    if (is.null(ref) || is.null(p1)) return(NULL)
    tot <- p1 + p2
    ok <- is.finite(tot) & tot > 0 & ref != alt
    B  <- ok & p2 == ref & p1 == alt      # BBQ carries reference
    A  <- ok & p2 == alt & p1 == ref      # BBQ carries alt
    j  <- match(sid, g$Sample_ID)
    data.frame(sid = sid, well = g$Wells[j], par = g$parent_a[j],
               n_sites = sum(ok), bref_frac = sum(B) / sum(ok),
               d_Bref = sum(p1[B]) / sum(tot[B]) - 0.5,
               d_Balt = sum(p1[A]) / sum(tot[A]) - 0.5,
               stringsAsFactors = FALSE)
}))
stopifnot(nrow(r) > 0)

r$class <- ifelse(r$bref_frac < 0.05, 'REF_LIKE_PARENT',
           ifelse(r$d_Bref > FLOOR, 'clean',
           ifelse(r$d_Bref > BAD,   'mild', 'BAD')))

# ── sibling check: wells sharing a parent should behave alike ─────────
r$sib_base <- NA_real_; r$sib_note <- ''
for (p in unique(r$par)) {
    i  <- which(r$par == p)
    cl <- i[r$d_Bref[i] > FLOOR]
    if (!length(cl)) {
        r$sib_note[i] <- if (length(i) > 1) 'no_clean_member_whole_strain'
                         else 'single_well_no_sibling'
        next
    }
    b <- median(r$d_Bref[cl])
    r$sib_base[i] <- b
    r$sib_note[i] <- ifelse(r$d_Bref[i] < b - DIVERGE, 'DIVERGES_from_siblings', 'agrees')
}

o <- order(r$d_Bref)
print(cbind(r[o, c('sid','well','par','n_sites')],
            round(r[o, c('bref_frac','d_Bref','d_Balt')], 4),
            r[o, c('class','sib_note')]), row.names = FALSE)

cat('\n── class tally ──\n'); print(table(r$class))
cat('\n── wells diverging from their siblings ──\n')
dv <- r[r$sib_note == 'DIVERGES_from_siblings', ]
if (nrow(dv)) print(cbind(dv[c('well','par')],
    round(dv[c('d_Bref','sib_base')], 4)), row.names = FALSE) else cat('  none\n')

cat('\n── parents holding BOTH clean and ruined wells (sink signature) ──\n')
hit <- FALSE
for (p in unique(r$par)) {
    v <- r[r$par == p, ]
    if (nrow(v) > 1 && any(v$d_Bref > FLOOR) && any(v$d_Bref < BAD)) {
        hit <- TRUE
        cat(sprintf('  %-4s clean %-20s ruined %s\n', p,
            paste(v$well[v$d_Bref > FLOOR], collapse = ','),
            paste(v$well[v$d_Bref < BAD],  collapse = ',')))
    }
}
if (!hit) cat('  none\n')

cat('\n── whole strains with no clean well (wrong stock / bad reference) ──\n')
ws <- unique(r$par[r$sib_note == 'no_clean_member_whole_strain'])
if (length(ws)) for (p in ws) {
    v <- r[r$par == p, ]
    cat(sprintf('  %-4s wells %-20s d_Bref %s\n', p,
        paste(v$well, collapse = ','),
        paste(sprintf('%.3f', v$d_Bref), collapse = ' ')))
} else cat('  none\n')

# ── did round 2 land? ────────────────────────────────────────────────
r2 <- c('E06','E07','D06','G07','E05','D10','F12','B02','C05','F02','H08','B12',
        'H01','H09','A11','D09')
cat('\n── the 16 wells changed in round 2 ──\n')
v <- r[r$well %in% r2, ]
print(cbind(v[c('well','par')], round(v[c('d_Bref','d_Balt')], 4), v['class']),
      row.names = FALSE)
cat(sprintf('\n  clean: %d / %d\n', sum(v$class == 'clean'), nrow(v)))

out <- file.path(base, 'qc', 'parent_verify.csv')
write.csv(r, out, row.names = FALSE)
cat('\nwrote', out, '\n')
