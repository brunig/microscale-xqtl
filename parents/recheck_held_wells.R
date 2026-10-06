# ══════════════════════════════════════════════════════════════════════
# recheck_held_wells.R — resolve the wells reassign_parents.R held back
# ══════════════════════════════════════════════════════════════════════
# The 2026-08-24 run applied 25 wells and held 33: 19 as
# `hold_cluster` (near-clonal group, unresolved), 13 as `hold_linkage`
# (scan proposed a change, the linkage check contradicted it), 1 as
# `hold_numeric` (A01 — see diag_well_A01.R). This script attacks the
# first two groups, each for a different reason.
#
# ── Part A: are the clusters real? ───────────────────────────────────
# The scan reports a cluster whenever the top candidates are within
# ll_cluster nats of each other, but nothing so far has checked whether
# those candidates are ACTUALLY near-identical strains. That distinction
# decides everything:
#   near-identical  -> no statistic can separate them, and none needs
#                      to: any member phases the sample the same way, so
#                      the well is resolved, not questionable.
#   far apart       -> a tiny d_2nd means something is wrong with the
#                      WELL (mixed pool, failed cross, contamination),
#                      which is a more serious finding than a mislabel.
# Five unrelated nominals resolving into ABB/ACK/ADH is what makes this
# worth checking before anything else.
#
# ── Part B: a linkage statistic with actual power ────────────────────
# The check that held 13 wells was the lag-1 autocorrelation of binned
# allele frequency. It returned 0.08-0.31 for BOTH candidates on every
# well — in a real cross adjacent 4 kb bins sit inside one linkage block
# and should correlate far higher than that, so the statistic was
# reading mostly noise while acting as sole arbiter.
#
# Replaced with a statistic aimed at what actually distinguishes a wrong
# parent. A wrong-but-related candidate still shows the real linkage
# structure (it phases most sites correctly), so "is there structure" is
# the wrong question. What it does is misassign alleles at the sites
# where it differs from the true parent, injecting site-level noise
# around the local trend. So: measure the residual around a local smooth
# and compare it to what binomial sampling alone would give.
#   overdispersion ~ 1  -> residuals are pure sampling noise: correct
#   overdispersion >> 1 -> excess scatter: misphased sites
#
# Simulated on a chromosome with realistic linkage structure, varying
# the fraction of misphased sites and the depth (od for the correct
# parent, then for a wrong one):
#
#   misphased   5x           12x          25x
#   5%          0.99 / 1.04  0.99 / 1.13  0.99 / 1.34
#   15%         1.01 / 1.14  0.99 / 1.37  1.01 / 1.90
#   35%         0.99 / 1.22  1.02 / 1.67  1.00 / 2.53
#   50%         1.01 / 1.26  0.98 / 1.72  1.00 / 2.67
#
# Two things follow. The null is calibrated — the correct parent gives
# 1.00 in every cell — so an ABSOLUTE threshold is meaningful and the
# statistic does not depend on having a good comparator. And the power
# is real at 12x and above but thin at 5x, so low-depth wells get no
# call rather than a coin flip.
#
# A per-site outlier count is reported alongside but deliberately NOT
# used for the verdict: it stops being monotonic once misphasing is
# heavy, because the local smooth is then itself corrupted.
#
# Self-contained: rebuilds the scan's objects from
# scan_parent_identity.R if they are not in the session. Writes
# qc/held_well_recheck.csv and changes nothing else.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════

scan_src = '/u/project/kruglyak/bruni/projects/xQTLpaper/files/scan_parent_identity.R'

if (!exists('score_sample') || !exists('Z') || !exists('sel_idx')) {
    cat('scan objects not in session — rebuilding setup from\n  ',
        scan_src, '\n', sep = '')
    if (!file.exists(scan_src)) stop('cannot find ', scan_src)
    .ln = readLines(scan_src, warn = FALSE)
    .at = grep('^# ── S6\\.', .ln)
    if (!length(.at))
        stop('could not find the "# ── S6." marker in ', scan_src)
    eval(parse(text = paste(.ln[seq_len(.at[1] - 1L)], collapse = '\n')),
         envir = globalenv())
    cat('setup rebuilt.\n\n')
}
stopifnot(exists('Z'), exists('Zu'), exists('sel_idx'), exists('read_counts'))

suppressPackageStartupMessages({ library(readr); library(dplyr) })

well.file = file.path(qc_dir, 'parent_identity_scan_well.csv')
out.file  = file.path(qc_dir, 'held_well_recheck.csv')
log.file  = file.path(qc_dir, 'parent_reassignment_log.csv')
stopifnot(file.exists(well.file))
wells_df = read_csv(well.file, show_col_types = FALSE, progress = FALSE)

# Tier labels come from reassign_parents.R; fall back to the scan's own
# verdicts if that log has not been written.
tiers = NULL
if (file.exists(log.file)) {
    lg = read_csv(log.file, show_col_types = FALSE, progress = FALSE)
    if ('tier_label' %in% names(lg))
        tiers = lg %>% distinct(Wells, tier_label)
}
wd = if (is.null(tiers)) wells_df else left_join(wells_df, tiers, by = 'Wells')
if (!'tier_label' %in% names(wd))
    wd$tier_label = ifelse(wd$verdict == 'ambiguous_cluster',
                           'hold_cluster', NA_character_)

# Site set: the same one the scan scored on.
scored = logical(m_all); scored[sel_idx] = TRUE


# ── A1. Genotype distance between two strains ─────────────────────────
# Counted only where BOTH strains are homozygous, on the scored sites,
# so it is the same comparison the likelihood makes.
geno_of = local({
    cache = new.env(parent = emptyenv())
    function(s) {
        if (!is.null(cache[[s]])) return(cache[[s]])
        gu = as.numeric(Zu[s, ])
        gh = if (has_het) (as.numeric(Z[s, ]) - 2 * gu) > 0
             else rep(FALSE, m_all)
        v = list(gu = gu, hom = !gh)
        cache[[s]] = v
        v
    }
})

geno_dist = function(s1, s2) {
    if (!(s1 %in% zrows) || !(s2 %in% zrows))
        return(c(d = NA_real_, n = NA_real_))
    g1 = geno_of(s1); g2 = geno_of(s2)
    both = scored & g1$hom & g2$hom
    c(d = sum(both & (g1$gu != g2$gu)), n = sum(both))
}


# ── A2. Calibration: what does "far apart" look like here? ────────────
# Distances mean nothing without a scale, so measure a sample of random
# strain pairs from Z on the same sites.
set.seed(1)
cal_strains = sample(zrows, min(24, length(zrows)))
cal = c()
for (i in seq_along(cal_strains))
    for (j in seq_len(i - 1L))
        cal = c(cal, geno_dist(cal_strains[i], cal_strains[j])[['d']])
cal = cal[is.finite(cal)]
cat('── calibration: genotype distance between random 1011 strains ──\n')
print(round(quantile(cal, c(0, .05, .25, .5, .75, 1), na.rm = TRUE)))
cal_med = median(cal, na.rm = TRUE)
cat('median unrelated distance:', format(round(cal_med), big.mark = ','),
    'sites of', format(sum(scored), big.mark = ','), 'scored\n\n')

# Empirical conversion from the 2026-08-24 simulations: a candidate
# differing at D sites separates by roughly 1.2 nats per read at those
# sites, i.e. Delta ~ D * depth * 1.2. Used only to sanity-check whether
# an observed d_2nd is the size the genotype distance predicts.
predict_delta = function(D, depth) D * depth * 1.2


# ── A3. Cluster reality check ─────────────────────────────────────────
clust_rows = wd %>% filter(!is.na(cluster), nzchar(cluster),
                           n_cluster > 1)
cat('── Part A: cluster reality (', nrow(clust_rows), ' wells ) ──\n',
    sep = '')
A = list()
for (i in seq_len(nrow(clust_rows))) {
    w  = clust_rows$Wells[i]
    mem = trimws(strsplit(clust_rows$cluster[i], '/')[[1]])
    mem = mem[mem %in% zrows]
    if (length(mem) < 2) next
    dd = c(); for (p in seq_along(mem)) for (q in seq_len(p - 1L))
        dd = c(dd, geno_dist(mem[p], mem[q])[['d']])
    depth = if ('n_reads' %in% names(clust_rows))
        clust_rows$n_reads[i] / max(1, sum(scored)) else NA_real_
    A[[length(A) + 1]] = tibble(
        Wells = w, nominal = clust_rows$nominal[i],
        best = clust_rows$best[i],
        cluster = clust_rows$cluster[i],
        max_member_dist = max(dd, na.rm = TRUE),
        med_unrelated   = cal_med,
        rel_to_unrelated = max(dd, na.rm = TRUE) / cal_med,
        observed_d_2nd  = clust_rows$d_2nd[i],
        predicted_d_2nd = predict_delta(max(dd, na.rm = TRUE), depth),
        verdict = NA_character_)
}
A = bind_rows(A)
if (nrow(A)) {
    A = A %>% mutate(verdict = case_when(
        rel_to_unrelated <= 0.02 ~ 'clonal_group_pick_any',
        rel_to_unrelated <= 0.10 ~ 'close_relatives_low_impact',
        TRUE                     ~ 'NOT_CLONAL_investigate_well'))
    print(as.data.frame(A %>% select(Wells, nominal, best, cluster,
        max_member_dist, rel_to_unrelated, observed_d_2nd,
        predicted_d_2nd, verdict) %>% arrange(desc(rel_to_unrelated))),
        row.names = FALSE)
    cat('\nclonal_group_pick_any      : members differ at <2% of the\n',
        '  unrelated baseline. Nothing can separate them and nothing\n',
        '  needs to — apply `best` and record the group.\n', sep = '')
    cat('NOT_CLONAL_investigate_well: members are genuinely different\n',
        '  strains, so a d_2nd of a few nats is not a resolution limit.\n',
        '  Suspect a mixed or failed pool in that well.\n', sep = '')
}


# ── B1. Local-smooth overdispersion ───────────────────────────────────
# f = frequency of the FIXED parent's allele at sites segregating
# between the fixed parent and `cand`. Compare the scatter of f around a
# local mean to the binomial expectation. Correct parent -> ~1.
linkage_fit = function(a, b, r, cand, win = 501, min_n = 3) {
    g = geno_of(cand)
    seg = scored & g$hom & ((r == 0 & g$gu == 1) | (r == 1 & g$gu == 0))
    idx = which(seg & (a + b) >= min_n)
    if (length(idx) < 5 * win)
        return(c(od = NA_real_, frac_out = NA_real_, n = length(idx),
                 depth = NA_real_))

    n  = a[idx] + b[idx]
    f  = ifelse(r[idx] == 0, b[idx], a[idx]) / n
    ch = site_chr[idx]; po = site_pos[idx]

    od_num = c(); od_den = c(); outl = c()
    for (cc in unique(ch)) {
        j = which(ch == cc); j = j[order(po[j])]
        if (length(j) < 2 * win) next
        fj = f[j]; nj = n[j]
        # local mean over a centred window of `win` sites
        loc = as.numeric(stats::filter(fj, rep(1 / win, win), sides = 2))
        keep = !is.na(loc) & loc > 0 & loc < 1
        if (!any(keep)) next
        res = fj[keep] - loc[keep]
        exp_var = loc[keep] * (1 - loc[keep]) / nj[keep]
        od_num = c(od_num, res^2)
        od_den = c(od_den, exp_var)
        outl   = c(outl, abs(res) / sqrt(exp_var) > 4)
    }
    if (!length(od_num)) return(c(od = NA_real_, frac_out = NA_real_,
                                  n = length(idx), depth = NA_real_))
    c(od = mean(od_num) / mean(od_den),
      frac_out = mean(outl),
      n = length(od_num),
      depth = mean(n))
}

# Below this the statistic has too little power to arbitrate (see the
# 5x column of the table in the header).
min_od_depth = 8

# Deepest sample of a well, so the statistic has something to work with.
deepest_sample = function(w) {
    sub = samples[samples$Wells == w, , drop = FALSE]
    best = NULL; bestn = -1
    for (sid in sub$Sample_ID) {
        cn = read_counts(sid); if (is.null(cn)) next
        tot = sum((cn$a + cn$b)[sel_idx])
        if (tot > bestn) { bestn = tot; best = list(sid = sid, cn = cn) }
    }
    best
}


# ── B2. Re-arbitrate the linkage holds ────────────────────────────────
# Part B runs on EVERY well that has a candidate, not just the holds.
# The 25 already-reassigned wells (scan_linkage / both_assays /
# rescue_unphasable) were never od-checked, and three of them (C01, G09,
# G10) are the worst on the plate.
lk_rows = wd %>% filter(!is.na(best), nzchar(best))
cat('\n── Part B: linkage re-check (', nrow(lk_rows), ' wells ) ──\n',
    sep = '')
B = list()
for (i in seq_len(nrow(lk_rows))) {
    w   = lk_rows$Wells[i]
    nom = lk_rows$nominal[i]; bst = lk_rows$best[i]
    ds  = deepest_sample(w)
    if (is.null(ds)) next
    fx  = samples$fixed[match(w, samples$Wells)]
    if (!fx %in% zrows) next
    # Candidate or nominal absent from Z (e.g. YBB wells F09/G08): no
    # genotype to phase against, so there is nothing to arbitrate.
    if (!nom %in% zrows || !bst %in% zrows) next
    r   = as.numeric(Z[fx, ]) / vmax
    cc  = cap_depth(ds$cn$a, ds$cn$b, sel_idx)

    fn = linkage_fit(cc$a, cc$b, r, nom)
    fb = linkage_fit(cc$a, cc$b, r, bst)
    dg = geno_dist(nom, bst)[['d']]

    B[[length(B) + 1]] = tibble(
        Wells = w, sample = ds$sid, nominal = nom, best = bst,
        dist_nom_best = dg,
        od_nominal = fn[['od']],   od_best = fb[['od']],
        out_nominal = fn[['frac_out']], out_best = fb[['frac_out']],
        n_sites_nom = fn[['n']],   n_sites_best = fb[['n']],
        depth_nom = fn[['depth']], depth_best = fb[['depth']],
        old_ac_nominal = lk_rows$autocor_nominal[i],
        old_ac_best    = lk_rows$autocor_best[i],
        d_2nd = lk_rows$d_2nd[i], d_nominal = lk_rows$d_nominal[i])
    cat(sprintf('  %-4s %-4s -> %-4s   od %6.2f -> %6.2f   out %5.3f -> %5.3f\n',
        w, nom, bst, fn[['od']], fb[['od']], fn[['frac_out']],
        fb[['frac_out']]))
}
B = bind_rows(B)
if (nrow(B)) {
    # Lower overdispersion = fewer misphased sites = better parent.
    # Absolute test first, because the null is calibrated at 1.00: a
    # candidate under 1.10 is phasing cleanly on its own terms, no
    # comparator needed. The relative test only breaks ties.
    B = B %>% mutate(verdict = case_when(
        is.na(od_nominal) | is.na(od_best) ~ 'no_call_insufficient_sites',
        pmax(depth_nom, depth_best, na.rm = TRUE) < min_od_depth ~
                                             'no_call_low_depth',
        od_nominal < 1.10 & od_best >= 1.20 ~ 'KEEP_nominal',
        od_best    < 1.10 & od_nominal >= 1.20 ~ 'APPLY_best',
        od_nominal < 1.10 & od_best < 1.10  ~ 'both_clean_near_identical',
        od_best < od_nominal * 0.9          ~ 'APPLY_best_marginal',
        od_nominal < od_best * 0.9          ~ 'KEEP_nominal_marginal',
        TRUE                                ~ 'indistinguishable'))
    cat('\n')
    print(as.data.frame(B %>% select(Wells, nominal, best, dist_nom_best,
        od_nominal, od_best, out_nominal, out_best, verdict) %>%
        arrange(verdict, Wells)), row.names = FALSE)
    cat('\nod = residual scatter around a local smooth, over the binomial\n',
        'expectation. The correct parent gives ~1.00 regardless of depth\n',
        '(validated by simulation, see header), so read it absolutely:\n',
        '  < 1.10  phases cleanly\n',
        '  > 1.20  carries misphased sites\n',
        'both_clean_near_identical means neither candidate injects\n',
        'noise — check dist_nom_best; a small distance means the two are\n',
        'near-clones and the choice does not matter.\n', sep = '')
    cat('\nOld lag-1 autocorrelation is carried through for comparison.\n',
        'Where the two disagree, prefer od: the autocorrelation was\n',
        'measured on 4 kb bins and returned 0.08-0.31 for every\n',
        'candidate, which is noise, not signal.\n', sep = '')
}


# ── B3. Write ─────────────────────────────────────────────────────────
out = bind_rows(
    if (nrow(A)) A %>% mutate(part = 'cluster') else NULL,
    if (nrow(B)) B %>% mutate(part = 'linkage') else NULL)
if (nrow(out)) write_csv(out, out.file)

cat('\n══ Summary ══\n')
if (nrow(A)) {
    cat('cluster wells:\n')
    print(as.data.frame(A %>% count(verdict, name = 'n_wells')),
          row.names = FALSE)
}
if (nrow(B)) {
    cat('linkage wells:\n')
    print(as.data.frame(B %>% count(verdict, name = 'n_wells')),
          row.names = FALSE)
}
cat('\nwrote ', out.file, '\n', sep = '')
cat('\nTo act on this: for wells verdicted APPLY_best or\n',
    'clonal_group_pick_any, add their Wells to an override list and\n',
    'set parent_a from `best`. Nothing here edits the keyfile — the\n',
    'point is to tell you which of the 33 are actually resolvable.\n',
    sep = '')
