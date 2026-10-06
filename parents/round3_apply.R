# ══════════════════════════════════════════════════════════════════════
# round3_apply.R — build the key from the od table, one rule, no hand-edits
# ══════════════════════════════════════════════════════════════════════
# Run AFTER the extended Part B (recheck_held_wells.R with
#   lk_rows = wd %>% filter(!is.na(best), nzchar(best))
# so every well gets an od_nominal/od_best pair, not just the 33 held).
#
# ── The assignment rule ───────────────────────────────────────────────
# Calibrated 2026-08-25 against d_Bref on 9 wells: every well with
# od_best <= 1.45 came back at the clean floor; every well at >= 1.54
# kept residual skew.
#
#   od_best <= 1.45 and better than nominal  -> APPLY best
#   od_best  > 1.45 but better than nominal  -> keep best, flag
#   nominal better than best                 -> KEEP nominal
#
# The middle case is the one I got wrong by hand: "not the right parent"
# is not the same as "worse than what we had". Reverting H01/H09/A11 to
# nominal sent them back to od 6.30/2.36/3.96 from 2.21/1.54/1.98.
#
# ── The exclusion rule is NOT the same threshold ──────────────────────
# OD_APPLY answers "do I trust this reassignment"; it does not answer
# "is this parent wrong". The 2026-08-26 full-plate run gives the null
# directly: among the 34 wells the scan CONFIRMED (nominal == best), od
# runs median 1.46, p90 2.12, p95 2.36, max 3.96 — half of all
# known-good wells sit above 1.45. Excluding there would drop 48 of 92
# wells, 18 of them scan-confirmed. The floor is per-parent, not
# per-well (ABG 2.34/2.40, ABP 1.76/2.23, AVB 1.50/1.86, ABS 1.54/1.86
# across their own separate wells), i.e. reference quality, not a
# mislabel. OD_EXCLUDE is ~p97 of that null.
#
# ── no_call is not evidence of anything ──────────────────────────────
# Part B declines to call two ways, and they mean different things:
#   no_call_low_depth (D06, G02, H06): pmax(depth) < 8. The statistic
#     is unpowered below ~12x, so an od here is a coin flip. Change
#     NOTHING — keep whatever the current key says — and watch.
#   no_call_insufficient_sites: fires when EITHER candidate has < 5*win
#     informative sites. If only the NOMINAL is NA (E03) that is itself
#     informative — the nominal cannot be fit at all — so the normal
#     rule applies and best wins. If BOTH are NA (A05, B06) there is no
#     information; keep current and watch.
# A05 is excluded regardless, on Part A's independent NOT_CLONAL verdict.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
base <- '/u/project/kruglyak/bruni/projects/xQTLpaper'
cur  <- file.path(base, 'files/flatkey_override2.csv')
rec  <- file.path(base, 'qc/held_well_recheck.csv')
out  <- file.path(base, 'files/flatkey_override3.csv')
OD_APPLY   <- 1.45   # trust a REASSIGNMENT at or below this
OD_EXCLUDE <- 2.50   # ~p97 of the confirmed-well null; above = drop well
# Minimum relative improvement before a label is allowed to flip. The
# 2026-08-26 run flipped D03 ACG->ACH on od 1.239678 vs 1.239172 — a
# 0.04% margin between two strains that differ at ONE site of 665,567.
# That is noise choosing a label. Real corrections clear this by miles
# (1.08 vs 4.53 and the like), so the guard only suppresses churn.
MIN_MARGIN <- 0.02

stopifnot(file.exists(cur), file.exists(rec))
k   <- read.csv(cur, check.names = FALSE)
all <- read.csv(rec, check.names = FALSE)
tw  <- function(x) trimws(as.character(x))

R <- all[all$part == 'linkage' & !is.na(all$Wells), ]
R <- R[!duplicated(R$Wells), ]
if (!'verdict' %in% names(R)) R$verdict <- NA_character_

# Part A's own verdict, independent of od: a well whose candidate
# cluster is NOT near-clonal has a well problem, not a label problem.
A <- all[all$part == 'cluster' & !is.na(all$Wells), ]
notclonal <- tw(A$Wells[A$verdict %in% 'NOT_CLONAL_investigate_well'])

cat(sprintf('recheck covers %d wells (of %d on the plate)\n',
            nrow(R), length(unique(k$Wells))))
if (length(notclonal))
    cat('Part A flags NOT_CLONAL (well problem, not label): ',
        paste(sort(notclonal), collapse = ', '), '\n', sep = '')

# `now` = what the CURRENT key says for this well, which is what a
# no_call must fall back to. Falling back to `nominal` instead would
# silently undo round-2 corrections on the strength of no evidence.
now <- setNames(tw(k$parent_a)[match(tw(R$Wells), tw(k$Wells))], tw(R$Wells))

pick <- function(w, nom, bst, odn, odb, vd) {
    nom <- tw(nom); bst <- tw(bst); cu <- unname(now[tw(w)])
    if (isTRUE(grepl('^no_call_low_depth', vd)))
        return(c(cu, 'hold_low_depth', NA))
    if (is.na(odb) && is.na(odn))  return(c(cu, 'hold_no_information', NA))
    if (is.na(odb))                return(c(nom, 'keep_nominal', odn))
    better <- is.na(odn) || odb < odn * (1 - MIN_MARGIN)
    if (!better)                   return(c(nom, 'keep_nominal',
                                            if (is.na(odn)) odb else odn))
    if (identical(nom, bst))       return(c(nom, 'confirmed', odb))
    if (odb <= OD_APPLY)           return(c(bst, 'apply_best', odb))
    c(bst, 'keep_best_flagged', odb)
}
d <- t(mapply(pick, R$Wells, R$nominal, R$best,
              R$od_nominal, R$od_best, R$verdict))
R$parent_new <- d[, 1]
R$decision   <- d[, 2]
R$od_final   <- suppressWarnings(as.numeric(d[, 3]))

k$parent_a_prev <- k$parent_a
i   <- match(tw(k$Wells), tw(R$Wells))
hit <- !is.na(i)
k$parent_a[hit]  <- R$parent_new[i[hit]]
k$od_decision    <- NA_character_
k$od_decision[hit] <- R$decision[i[hit]]
k$od_final       <- NA_real_
k$od_final[hit]  <- R$od_final[i[hit]]
k$cross <- paste(trimws(k$parent_a), trimws(k$parent_alpha), sep = 'x')
write.csv(k, out, row.names = FALSE)

cat('\n── decisions (per well) ──\n'); print(table(R$decision))
ch <- unique(k[tw(k$parent_a) != tw(k$parent_a_prev),
               c('Wells','parent_a_prev','parent_a','od_decision','od_final')])
cat(sprintf('\n%d wells change vs override2 (%d samples)\n',
            nrow(ch), 4 * nrow(ch)))
print(ch[order(ch$Wells), ], row.names = FALSE)

# ── Exclusion ─────────────────────────────────────────────────────────
# Three independent grounds, none of them "od is merely elevated":
#   1. the parent we ended up using does not phase cleanly
#   2. Part A says the well itself is suspect
#   3. the well never got a recheck row at all (parents absent from Z)
ex_od   <- tw(k$Wells[!is.na(k$od_final) & k$od_final > OD_EXCLUDE])
ex_well <- intersect(notclonal, tw(k$Wells))
ex_none <- setdiff(tw(k$Wells), tw(R$Wells))
ex <- sort(unique(c(ex_od, ex_well, ex_none)))

cat('\n── wells to EXCLUDE from QTL analysis ──\n')
bad <- R[tw(R$Wells) %in% ex,
         c('Wells','nominal','best','od_nominal','od_best',
           'parent_new','decision','od_final','verdict')]
print(bad[order(bad$od_final, decreasing = TRUE), ], row.names = FALSE)
cat(sprintf('\n  od_final > %.2f      : %s\n', OD_EXCLUDE,
            paste(sort(ex_od), collapse = ', ')))
cat(sprintf('  Part A well problem : %s\n',
            paste(sort(ex_well), collapse = ', ')))
cat(sprintf('  no recheck row      : %s\n',
            paste(sort(ex_none), collapse = ', ')))
cat(sprintf('\n%d wells / %d samples excluded\n', length(ex), 4 * length(ex)))
writeLines(ex, file.path(base, 'qc/unresolved_wells.txt'))

# Kept but not pristine: best available parent, inside the confirmed-well
# null, but above the reassignment threshold — or held for no evidence.
watch <- sort(unique(tw(k$Wells[
    !(tw(k$Wells) %in% ex) &
    ((!is.na(k$od_final) & k$od_final > OD_APPLY) |
     k$od_decision %in% c('hold_low_depth','hold_no_information'))])))
cat(sprintf('\n%d wells kept on watch: %s\n',
            length(watch), paste(watch, collapse = ', ')))
writeLines(watch, file.path(base, 'qc/watch_wells.txt'))
cat('\nwrote', out, '\n  qc/unresolved_wells.txt, qc/watch_wells.txt\n')
