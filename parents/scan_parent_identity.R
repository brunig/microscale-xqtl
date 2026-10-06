# ══════════════════════════════════════════════════════════════════════
# scan_parent_identity.R — which strain is actually in each well?
# ══════════════════════════════════════════════════════════════════════
# Scores the raw ASEReadCounter counts of every xQTLpaper sample against
# every strain in Z under a three-class binomial model, on a site set
# that is FIXED across candidates so the likelihoods are directly
# comparable. Aggregates by well (the four conditions of a well are the
# same parent, so they are four independent looks at one question) and
# reports, per well, the best-fitting strain, how far it sits above the
# runner-up, and where the keyfile's nominal assignment ranks.
#
# Why not AFD: `afd_failed` collapses low depth, library quality and
# wrong parent into one bit, and the segregating-site set moves with the
# assumed parent, so AFDs computed under different assignments are not
# on a common footing and cannot be ranked. This works on raw counts.
#
# The model. With the partner parent fixed, the expected alt fraction at
# site i is (g_fixed,i + g_cand,i)/2, which takes only three values:
#   both parents ref  -> ~0      (eps)
#   both parents alt  -> ~1      (1 - eps)
#   parents differ    ->  0.5
# so  loglik = SUM_i [ a_i*log(p_i) + (n_i - a_i)*log(1 - p_i) ].
# A wrong candidate puts truly-segregating sites into a fixed class
# (observed ~0.5 where 0 or 1 was predicted, costing ~log(eps) per read)
# and truly-shared sites into the mixed class. Every discordant site
# between the assumed and the true parent contributes.
#
# Because p_i takes only a few values the sum collapses to grouped
# totals of alt and ref reads, and group membership is decided by
# g_cand vs g_fixed — so scoring all ~919 candidates costs TWO sparse
# matrix products per sample over `cbind(a0, b0, a1, b1)`: one against
# the alt-homozygote indicator and one against Z, whose difference
# isolates the heterozygous calls. Exhaustive therefore costs the same
# as a handful of guesses, so there is no reason to restrict the search
# to the strains believed to be on the plate: a re-array error can pull
# any strain from the source plate, and the ~918 wrong answers calibrate
# what a good match looks like.
#
# Memory note: Z is never subset or transposed, but one derived matrix
# (Zu, the alt-homozygote indicator) is held alongside it. Measured on
# Z_1011_filtered.qs — 919 x 1,240,653, nnz 38.3M — that is roughly
# 460 MB each, so ~1 GB for the pair plus a 40 MB dense block per
# sample. -m 8 is enough; -m 16 is comfortable. Nothing is forked, and
# it should stay that way: do not wrap this in pbmclapply, since forking
# a matrix this size is what OOM-killed the §2 workers.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
    library(Matrix)
    library(readr)
    library(dplyr)
})
has_fread = requireNamespace('data.table', quietly = TRUE)


# ── S0. Configuration ─────────────────────────────────────────────────
proj_dir  = '/u/project/kruglyak/bruni/projects/xQTLpaper'
key.file  = file.path(proj_dir, 'files', 'flatkey.csv')
count_dir = file.path(proj_dir, 'count_files')
qc_dir    = file.path(proj_dir, 'qc')
Z.qs.file = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/Z_1011_filtered.qs'

out.sample = file.path(qc_dir, 'parent_identity_scan_sample.csv')
out.well   = file.path(qc_dir, 'parent_identity_scan_well.csv')
# Written in the shape reassign_parents.R expects for its optional join.
out.join   = file.path(qc_dir, 'cross_identity_scan.csv')

# ---- Which parent is being identified -------------------------------
# 'a'     : fix parent_alpha (BBQ), scan candidates for parent_a
# 'alpha' : fix parent_a, scan candidates for parent_alpha — use this to
#           check BBQ itself once parent_a is settled. Run it on the
#           wells that came back clean, where the fixed side is trusted.
scan_role = 'a'

# ---- Model / site selection -----------------------------------------
eps            = 0.005  # per-read error floor; see the eps note in §S9
# Minimum number of alt-homozygous strains for a site to be scored.
# Deliberately low. The instinct is to demand a common variant, but for
# IDENTITY the opposite is true: a variant carried by two strains out of
# 919 is precisely what separates those two from everyone else, and
# common variants are nearly uninformative about which strain you have.
# Raising this throws away exactly the sites that resolve near-clonal
# candidates. Set to 2 rather than 1 so every scored variant is
# corroborated by a second strain, which is the cheap guard against
# single-strain genotyping error. Raise it only if a ranking looks
# driven by a handful of rare calls.
min_alt_strains = 2
# Drop a site if more than this FRACTION of the 919 strains are called
# het there — a quality cut aimed at collapsed paralogs and mapping
# artefacts, not at the het calls themselves, which are modelled in §S4.
# At 0.10 (>91 strains) this removes just over 1% of sites.
max_het_strains = 0.10
min_reads      = 5e4    # below this a sample gets no call of its own
depth_cap_mult = 5      # cap per-site depth at this * median (see below)

# An xQTL pool is rarely exactly 50/50 — uneven mating, bottlenecks and
# genome-wide selection all skew it. Holding the mixed class at p = 0.5
# charges that skew against the TRUE parent. Profiling the parameter
# instead costs nothing (it falls out of the same matmul), is identical
# on a balanced pool, and separates ~50% better on a badly skewed one.
profile_mix    = TRUE

# ---- Verdict thresholds, in nats of log-likelihood ------------------
# Real separations here run to 1e4-1e6, so these are deliberately
# permissive; the raw deltas are in the output, calibrate from them.
ll_decisive = 100       # best must beat 2nd by this to be a clean call
ll_cluster  = 100       # candidates within this of best = same cluster

# ---- Spatial coherence confirmation ---------------------------------
# The correct parent gives allele frequency that varies smoothly along
# the chromosome (linkage blocks); a wrong-but-related parent gives
# discordance scattered at random with respect to position. Computed for
# the best and the nominal candidate on each well's deepest sample.
do_autocor  = TRUE
ac_bin      = 200       # consecutive segregating sites per bin

uchr = paste0('chr', as.character(as.roman(1:16)))


# ── S1. Load Z, detect its coding, derive the site mask ───────────────
# An earlier version of this block kept only sites with zero het calls
# across all strains, so that Z would be {0,2} on the survivors and the
# alt indicator would be exactly Z/2 — one matmul, no companion matrix.
# That filter is untenable: it requires no het among ~919 strains, so
# even a 1% per-strain het rate leaves 0.99^919 ~ 1e-4 of sites, and the
# 1011 collection does contain heterozygous diploids. Hets are now kept
# and modelled (see §S4), at the cost of one extra matrix.
# -Claude generated
cat('loading Z ...\n')
Z.obj = qs2::qs_read(Z.qs.file)
Z     = Matrix::drop0(Z.obj$Z)   # explicit zeros would corrupt the nnz math
stopifnot(inherits(Z, 'dgCMatrix'))
zrows = rownames(Z)
m_all = ncol(Z)
cat('  Z:', nrow(Z), 'strains x', m_all, 'variants, nnz', length(Z@x), '\n')

# ---- Detect the genotype coding rather than assuming it -------------
vmax = max(Z@x); vmin = min(Z@x)
vtab = table(Z@x[seq(1, length(Z@x), length.out = min(2e6, length(Z@x)))])
cat('  stored values (sampled): ',
    paste(sprintf('%s=%s', names(vtab), format(as.integer(vtab),
          big.mark = ',')), collapse = '  '), '\n', sep = '')
if (!(vmax %in% c(1, 2)) || vmin <= 0)
    stop('unexpected genotype coding in Z: range [', vmin, ', ', vmax,
         ']. This scan expects 0/2 (+1 for het) or 0/1.')

# Alt-homozygote indicator. For 0/1/2 coding, x*(x-1)/2 maps 2 -> 1 and
# 1 (het) -> 0; for 0/1 coding it is Z itself and no het class exists.
has_het = (vmax == 2)
Zu = Z
if (has_het) {
    Zu@x = Z@x * (Z@x - 1) / 2
    Zu = Matrix::drop0(Zu)
    cat('  coding: 0=ref, 1=het, 2=alt\n')
} else {
    cat('  coding: 0=ref, 1=alt (no het class)\n')
}

alt_j = Matrix::colSums(Zu)                 # alt-homozygous strains
het_j = if (has_het) Matrix::colSums(Z) - 2 * alt_j else numeric(m_all)
cat('  het calls per site: median ', stats::median(het_j),
    ', mean ', round(mean(het_j), 2), ', sites with none ',
    sum(het_j == 0), ' (', round(100 * mean(het_j == 0), 1), '%)\n', sep = '')

# ---- Site mask ------------------------------------------------------
# Only two requirements now: the site must be polymorphic enough in the
# collection to carry information, and it must not be so het-ridden that
# the genotypes are unreliable. Whether a given CANDIDATE is het at a
# site is handled per candidate in the likelihood, not by dropping the
# site — dropping per candidate would give each candidate a different
# site set and make the likelihoods incomparable, which is the whole
# point of the fixed site set.
sel = (alt_j >= min_alt_strains) &
      (alt_j <= nrow(Z) - min_alt_strains) &
      (het_j <= max_het_strains * nrow(Z))
cat('  sites: ', sum(sel), ' usable of ', m_all,
    '  (', sum(alt_j < min_alt_strains | alt_j > nrow(Z) - min_alt_strains),
    ' non-informative, ', sum(het_j > max_het_strains * nrow(Z)),
    ' too het)\n', sep = '')
if (sum(sel) < 5e4)
    stop('too few usable sites (', sum(sel), '). Relax min_alt_strains ',
         '(now ', min_alt_strains, ') or max_het_strains (now ',
         max_het_strains, '); the value table printed above shows what ',
         'Z actually contains.')

# chr / pos per site, for the spatial coherence check
vid_short = sub(':.*$', '', Z.obj$variant_ids)
site_chr  = sub('_[0-9]+$', '', vid_short)
site_pos  = as.numeric(sub('^.*_', '', vid_short))
# numeric site key, so matching count files never touches strings
chr_i    = match(site_chr, uchr)
site_key = chr_i * 1e9 + site_pos
if (any(is.na(chr_i))) {
    cat('  note:', sum(is.na(chr_i)), 'sites on contigs outside uchr\n')
    # NA would match an NA in a count file's key vector; make it unmatchable
    site_key[is.na(site_key)] = -1
}

sel_idx = which(sel & !is.na(chr_i))
cat('  scoring on', length(sel_idx), 'sites\n')


# ── S2. Keyfile ───────────────────────────────────────────────────────
samples = read_csv(key.file, show_col_types = FALSE, progress = FALSE)
stopifnot(all(c('Sample_ID', 'Wells', 'parent_a', 'parent_alpha',
                'Drug3Let', 'oldLabel_parent_a') %in% names(samples)))
samples = samples %>%
    mutate(parent_a     = trimws(parent_a),
           parent_alpha = trimws(parent_alpha),
           nominal = if (scan_role == 'a') parent_a     else parent_alpha,
           fixed   = if (scan_role == 'a') parent_alpha else parent_a)

# Codes believed to be on this plate — used only to annotate whether the
# winner is one of them, never to restrict the search.
on_plate = sort(unique(samples$parent_a))

bad_fixed = setdiff(unique(samples$fixed), zrows)
if (length(bad_fixed)) {
    cat('\n!! fixed-side parents absent from Z (their samples are skipped):',
        paste(bad_fixed, collapse = ', '), '\n')
    samples = samples %>% filter(!(fixed %in% bad_fixed))
}


# ── S3. Count file reader ─────────────────────────────────────────────
# Returns full-length (m_all) alt / ref vectors, zero where the site had
# no coverage, so every candidate is scored on the same sites.
read_counts = function(sid) {
    f = file.path(count_dir, paste0(sid, '.txt'))
    if (!file.exists(f) || file.info(f)$size < 1000) return(NULL)
    cnt = tryCatch({
        if (has_fread) {
            data.table::fread(f, select = c('contig', 'position',
                                            'refCount', 'altCount'),
                              showProgress = FALSE, data.table = FALSE)
        } else {
            as.data.frame(read_tsv(f, show_col_types = FALSE,
                                   progress = FALSE))[,
                c('contig', 'position', 'refCount', 'altCount')]
        }
    }, error = function(e) NULL)
    if (is.null(cnt) || !nrow(cnt)) return(NULL)

    ckey = match(cnt$contig, uchr) * 1e9 + as.numeric(cnt$position)
    mi   = match(site_key, ckey)          # site -> row of the count file
    a = numeric(m_all); b = numeric(m_all)
    ok = !is.na(mi)
    a[ok] = as.numeric(cnt$altCount[mi[ok]])
    b[ok] = as.numeric(cnt$refCount[mi[ok]])
    a[is.na(a)] = 0; b[is.na(b)] = 0
    list(a = a, b = b)
}

# Collapsed repeats and paralogs show up as a handful of very deep sites
# whose allele ratio is meaningless; left alone they can dominate a sum
# of ~1e6 log terms. Scale those sites back to the cap rather than
# dropping them, which keeps the site set identical across candidates.
cap_depth = function(a, b, idx) {
    n = a + b
    med = median(n[idx][n[idx] > 0])
    if (!is.finite(med) || med <= 0) return(list(a = a, b = b, cap = NA_real_))
    cap = depth_cap_mult * med
    hot = which(n > cap)
    if (length(hot)) {
        f = cap / n[hot]
        a[hot] = a[hot] * f
        b[hot] = b[hot] * f
    }
    list(a = a, b = b, cap = cap)
}


# ── S4. The scan: one sparse matmul scores every candidate ────────────
# r = fixed parent's alt indicator over the scored sites, in {0,1}.
# S0 = sites where the fixed parent is ref, S1 = where it is alt.
#   both-ref  reads: A0 = sum(a over S0) - sum(a over S0 & cand-alt)
#   both-alt  reads: A1 =                  sum(a over S1 & cand-alt)
#   mixed     reads: everything left over
# and sum(over cand-alt) for all four count vectors at once is a single
# product with Z, since the alt indicator is Z/2 on these sites.
# Two products now: Zu gives the alt-homozygous read sums directly, and
# Z gives 2*alt-hom + het, so het = (Z %*% c) - 2*(Zu %*% c). Everything
# else falls out by subtraction from the column totals.
score_sample = function(a, b, r) {
    s0 = numeric(m_all); s1 = numeric(m_all)
    s0[sel_idx] = as.numeric(r[sel_idx] == 0)
    s1[sel_idx] = as.numeric(r[sel_idx] == 1)

    C = cbind(a * s0, b * s0, a * s1, b * s1)   # m_all x 4, dense
    tot = colSums(C)
    if (sum(tot) < 1) return(NULL)

    U = as.matrix(Zu %*% C)                     # cand alt-hom
    H = if (has_het) as.matrix(Z %*% C) - 2 * U  # cand het
        else matrix(0, nrow(U), ncol(U))
    V = sweep(-(U + H), 2, tot, `+`)            # cand ref-hom

    # Columns of U/H/V are (a|S0, b|S0, a|S1, b|S1). With the fixed
    # parent homozygous at every scored site, the candidate's genotype
    # sorts each read into one of five classes:
    #   fixed ref & cand ref  -> pool is fixed ref   (p_alt = eps)
    #   fixed alt & cand alt  -> pool is fixed alt   (p_alt = 1 - eps)
    #   fixed ref & cand het  -> p_alt = 0.25
    #   fixed alt & cand het  -> p_alt = 0.75
    #   parents differ        -> mixed, ratio profiled below
    base = V[, 1] * log(eps)     + V[, 2] * log1p(-eps) +   # both ref
           U[, 3] * log1p(-eps)  + U[, 4] * log(eps)        # both alt

    # Het candidate: half the pool carries a heterozygote, so the pool
    # sits a quarter of the way off fixation rather than at 0.5. Held at
    # a balanced 0.25/0.75 rather than tied to the profiled ratio below,
    # which would cost the closed form for a minority of reads. frac_het
    # reports how much of the data this approximation touches.
    het_ll = H[, 1] * log(0.25) + H[, 2] * log(0.75) +
             H[, 3] * log(0.75) + H[, 4] * log(0.25)
    n_het  = H[, 1] + H[, 2] + H[, 3] + H[, 4]

    # Mixed class, split by which parent carries the alt allele so the
    # pool's parental ratio is identifiable:
    #   Fx = reads carrying the FIXED parent's allele
    #   Kd = reads carrying the CANDIDATE's allele
    Fx = U[, 2] + V[, 3]
    Kd = U[, 1] + V[, 4]
    p_fix = Fx / (Fx + Kd)

    ll = if (profile_mix) {
        pf = pmin(pmax(p_fix, eps), 1 - eps)
        base + het_ll + Fx * log(pf) + Kd * log1p(-pf)
    } else {
        base + het_ll + (Fx + Kd) * log(0.5)
    }
    # A candidate identical to the fixed parent has no mixed sites at all
    # (p_fix is 0/0). That is not a cross, so rule it out rather than
    # letting a NaN propagate into the ranking.
    ll[!is.finite(ll)] = -Inf

    names(ll) = zrows; names(p_fix) = zrows
    list(ll = ll, p_fix = p_fix, n_reads = sum(tot),
         frac_het = n_het / sum(tot))
}


# ── S5. Spatial coherence of a candidate ──────────────────────────────
# Frequency of the FIXED parent's allele at sites segregating between
# the fixed parent and the candidate, binned along each chromosome; the
# statistic is the lag-1 autocorrelation of those binned means. A
# correct parent gives linkage blocks and a high value; a wrong one
# gives position-independent noise and ~0.
autocor_candidate = function(a, b, r, cand) {
    gu = as.numeric(Zu[cand, ])                     # 1 = alt-homozygous
    gh = if (has_het) (as.numeric(Z[cand, ]) - 2 * gu) > 0
         else rep(FALSE, m_all)                     # het at this site
    ok = logical(m_all); ok[sel_idx] = TRUE
    # segregating and phasable: both parents homozygous and different
    seg = (r == 0 & gu == 1) | (r == 1 & gu == 0 & !gh)
    idx = which(ok & seg & (a + b) >= 3)
    if (length(idx) < 10 * ac_bin) return(NA_real_)
    # frequency of the fixed parent's allele
    f = ifelse(r[idx] == 0, b[idx], a[idx]) / (a[idx] + b[idx])
    ch = site_chr[idx]; po = site_pos[idx]
    acs = c()
    for (cc in unique(ch)) {
        j  = which(ch == cc); j = j[order(po[j])]
        if (length(j) < 3 * ac_bin) next
        bin = ceiling(seq_along(j) / ac_bin)
        mu  = tapply(f[j], bin, mean)
        if (length(mu) < 4 || sd(mu) == 0) next
        acs = c(acs, cor(mu[-length(mu)], mu[-1]))
    }
    if (!length(acs)) NA_real_ else mean(acs, na.rm = TRUE)
}


# ── S6. Walk the wells ────────────────────────────────────────────────
# Well-by-well so each count file is read exactly once and a well's
# samples can be aggregated and confirmed while still in memory.
wells = samples %>% distinct(Wells) %>% arrange(Wells) %>% pull(Wells)
cat('\nscanning', length(wells), 'wells /', nrow(samples), 'samples\n')

samp_rows = list(); well_rows = list()
r_cache = list()

for (wi in seq_along(wells)) {
    w   = wells[wi]
    sub = samples %>% filter(Wells == w)
    ll_w = NULL; pf_acc = NULL; fh_acc = NULL; reads_w = 0; nsamp_ok = 0
    keep = list()

    for (k in seq_len(nrow(sub))) {
        sid = sub$Sample_ID[k]
        fx  = sub$fixed[k]
        cn  = read_counts(sid)
        if (is.null(cn)) {
            samp_rows[[length(samp_rows) + 1]] = tibble(
                Sample_ID = sid, Wells = w, Drug3Let = sub$Drug3Let[k],
                nominal = sub$nominal[k], n_reads = 0,
                best = NA_character_, best_ll = NA_real_,
                d_2nd = NA_real_, rank_nominal = NA_integer_,
                d_nominal = NA_real_, note = 'no count file')
            next
        }
        cc = cap_depth(cn$a, cn$b, sel_idx)
        # r in {0, 0.5, 1}: the fixed parent's alt dosage. Only its
        # homozygous sites enter the score (s0/s1 below test == 0 and
        # == 1), and they are excluded identically for every candidate,
        # so comparability holds. Scale by vmax, not a hard-coded 2 —
        # under 0/1 coding dividing by 2 would make r==1 unreachable and
        # silently score every site as if the fixed parent were ref.
        if (is.null(r_cache[[fx]])) r_cache[[fx]] = as.numeric(Z[fx, ]) / vmax
        r  = r_cache[[fx]]

        sc = score_sample(cc$a, cc$b, r)
        if (is.null(sc)) next

        o    = order(sc$ll, decreasing = TRUE)
        nomi = match(sub$nominal[k], names(sc$ll))
        samp_rows[[length(samp_rows) + 1]] = tibble(
            Sample_ID = sid, Wells = w, Drug3Let = sub$Drug3Let[k],
            nominal = sub$nominal[k], n_reads = sc$n_reads,
            best = names(sc$ll)[o[1]], best_ll = unname(sc$ll[o[1]]),
            d_2nd = unname(sc$ll[o[1]] - sc$ll[o[2]]),
            rank_nominal = if (is.na(nomi)) NA_integer_ else match(nomi, o),
            d_nominal = if (is.na(nomi)) NA_real_ else unname(sc$ll[o[1]] - sc$ll[nomi]),
            p_fix = unname(sc$p_fix[o[1]]),
            frac_het = unname(sc$frac_het[o[1]]),
            note = if (sc$n_reads < min_reads) 'low reads' else NA_character_)

        # a well's conditions are the same parent: sum their logliks
        ll_w    = if (is.null(ll_w)) sc$ll else ll_w + sc$ll
        pf_acc  = if (is.null(pf_acc)) sc$p_fix * sc$n_reads
                  else pf_acc + sc$p_fix * sc$n_reads
        fh_acc  = if (is.null(fh_acc)) sc$frac_het * sc$n_reads
                  else fh_acc + sc$frac_het * sc$n_reads
        reads_w = reads_w + sc$n_reads
        nsamp_ok = nsamp_ok + 1
        keep[[sid]] = list(a = cc$a, b = cc$b, r = r, n = sc$n_reads)
    }

    if (is.null(ll_w)) {
        well_rows[[length(well_rows) + 1]] = tibble(
            Wells = w, strain = sub$oldLabel_parent_a[1],
            nominal = sub$nominal[1], n_samples = 0L, n_reads = 0,
            best = NA_character_, verdict = 'no_data')
        next
    }

    o   = order(ll_w, decreasing = TRUE)
    nom = sub$nominal[1]
    ni  = match(nom, names(ll_w))
    d2  = unname(ll_w[o[1]] - ll_w[o[2]])
    dn  = if (is.na(ni)) NA_real_ else unname(ll_w[o[1]] - ll_w[ni])
    clus = names(ll_w)[ll_w >= ll_w[o[1]] - ll_cluster]

    verdict =
        if (reads_w < min_reads)                      'insufficient_data'
        else if (identical(names(ll_w)[o[1]], nom))   'confirms_nominal'
        else if (!is.na(dn) && dn <= ll_cluster)      'nominal_in_cluster'
        else if (d2 < ll_decisive)                    'ambiguous_cluster'
        else                                          'REASSIGN'

    # confirm on the deepest sample of the well
    ac_b = NA_real_; ac_n = NA_real_
    if (do_autocor && length(keep) && verdict != 'insufficient_data') {
        deep = keep[[which.max(vapply(keep, function(x) x$n, 0))]]
        ac_b = autocor_candidate(deep$a, deep$b, deep$r, names(ll_w)[o[1]])
        if (!is.na(ni) && !identical(names(ll_w)[o[1]], nom))
            ac_n = autocor_candidate(deep$a, deep$b, deep$r, nom)
        else ac_n = ac_b
    }

    well_rows[[length(well_rows) + 1]] = tibble(
        Wells = w, strain = sub$oldLabel_parent_a[1], nominal = nom,
        n_samples = nsamp_ok, n_reads = reads_w,
        best = names(ll_w)[o[1]], best_ll = unname(ll_w[o[1]]),
        second = names(ll_w)[o[2]], d_2nd = d2,
        nominal_ll = if (is.na(ni)) NA_real_ else unname(ll_w[ni]),
        d_nominal = dn,
        rank_nominal = if (is.na(ni)) NA_integer_ else match(ni, o),
        d_2nd_per_10k = 1e4 * d2 / reads_w,
        d_nominal_per_10k = if (is.na(dn)) NA_real_ else 1e4 * dn / reads_w,
        n_cluster = length(clus),
        cluster = paste(head(clus, 8), collapse = '/'),
        p_fix_best = unname((pf_acc / reads_w)[o[1]]),
        frac_het_best = unname((fh_acc / reads_w)[o[1]]),
        best_on_plate = names(ll_w)[o[1]] %in% on_plate,
        autocor_best = ac_b, autocor_nominal = ac_n,
        verdict = verdict)

    rm(keep); if (wi %% 8 == 0) invisible(gc())
    cat(sprintf('  [%3d/%3d] %-4s %-5s -> %-9s  d2=%9.0f  d_nom=%10s  %s\n',
        wi, length(wells), w, nom, names(ll_w)[o[1]], d2,
        if (is.na(dn)) 'NA' else sprintf('%.0f', dn), verdict))
}

samp = bind_rows(samp_rows)
well = bind_rows(well_rows)


# ── S7. Write ─────────────────────────────────────────────────────────
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
write_csv(samp, out.sample)
write_csv(well %>% arrange(factor(verdict, levels = c(
    'REASSIGN', 'ambiguous_cluster', 'nominal_in_cluster',
    'confirms_nominal', 'insufficient_data', 'no_data')), desc(d_nominal)),
    out.well)

# Per-sample view in the shape reassign_parents.R joins on.
samples %>%
    select(Sample_ID, Wells) %>%
    left_join(well %>% select(Wells, best, d_2nd, rank_nominal, verdict),
              by = 'Wells') %>%
    transmute(Sample_ID,
              best_p1 = best,
              assigned_rank = rank_nominal,
              depth_flag = ifelse(verdict %in% c('insufficient_data', 'no_data'),
                                  'low', 'ok'),
              well_verdict = verdict,
              well_d_2nd = d_2nd) %>%
    write_csv(out.join)


# ── S8. Report ────────────────────────────────────────────────────────
cat('\n══ Verdicts by well ══\n')
print(as.data.frame(well %>% count(verdict, name = 'n_wells') %>%
    arrange(desc(n_wells))), row.names = FALSE)

if (any(well$verdict == 'REASSIGN')) {
    cat('\n── REASSIGN ──\n')
    print(as.data.frame(well %>% filter(verdict == 'REASSIGN') %>%
        arrange(desc(d_nominal)) %>%
        select(Wells, strain, nominal, best, d_nominal, d_2nd,
               rank_nominal, best_on_plate, autocor_nominal, autocor_best)),
        row.names = FALSE)
}
if (any(well$verdict == 'ambiguous_cluster')) {
    cat('\n── Ambiguous: best not separated from the cluster ──\n')
    print(as.data.frame(well %>% filter(verdict == 'ambiguous_cluster') %>%
        select(Wells, strain, nominal, best, d_2nd, n_cluster, cluster)),
        row.names = FALSE)
}

# Same nominal parent, different verdicts across its wells — the
# signature of a re-array error rather than a bad original assignment.
split_nom = well %>%
    filter(!verdict %in% c('insufficient_data', 'no_data')) %>%
    group_by(nominal) %>%
    filter(n_distinct(best) > 1) %>%
    ungroup()
if (nrow(split_nom)) {
    cat('\n── Wells sharing a nominal parent but resolving differently ──\n')
    cat('   (box sequencing cannot explain these; they are re-array errors)\n')
    print(as.data.frame(split_nom %>%
        select(nominal, strain, Wells, best, d_2nd, verdict) %>%
        arrange(nominal, Wells)), row.names = FALSE)
}

cat('\nwrote ', out.sample, '\n      ', out.well, '\n      ', out.join,
    '\n', sep = '')
cat('\nNext: feed the REASSIGN calls into reassign_parents.R (it joins\n',
    out.join, ' automatically), then re-run §1-§2 for the\n',
    'affected samples.\n', sep = '')


# ── S9. Reading the output ────────────────────────────────────────────
# d_2nd / d_nominal are in nats of log-likelihood, summed over every
# read at every scored site, so they are large. Simulated separation
# between the true parent and a decoy differing from it at the given
# fraction of sites, scaled to 1.2M sites:
#
#   divergence      2x      5x     12x     25x
#   0.03%          918    1619    5315   11953
#   0.1%          2645    7463   17612   36624
#   0.3%          8499   23773   52421  106455
#   1%           27614   71894  175769  359865
#
# So even the tightest near-clonal pairs on this plate separate cleanly
# — M22's top five Hamming hits span 437-465 units, i.e. ~0.03%
# divergence, which is still ~5000 nats apart at 12x. ll_decisive = 100
# is deliberately far below any of this: a call that comes back
# ambiguous at that threshold is genuinely ambiguous, not underpowered.
# Where n_cluster > 1 you still do not need the true name — you need a
# genotype vector that phases correctly, and any cluster member gives
# you that.
#
# frac_het_best is the fraction of scored reads sitting at sites where
# the winning candidate is called heterozygous. Those reads are modelled
# at a balanced 0.25/0.75 rather than at the profiled ratio, so if this
# is large (say >0.2) the approximation is carrying real weight and the
# call deserves a second look — most simply, re-run with
# max_het_strains tightened and check the winner does not move.
#
# p_fix_best is the fitted frequency of the FIXED parent's allele at
# segregating sites — the pool's parental ratio, which the scan
# estimates rather than assuming. ~0.5 is a healthy pool. Values far
# off (say <0.35 or >0.65) mean uneven mating, a bottleneck, or
# genome-wide selection, and are worth knowing about independently of
# identity: that well's AFD is being computed on a lopsided pool.
#
# `eps` sets what a site in the wrong class costs (~log(eps) per read).
# Lower makes the scan sharper but more sensitive to genotyping error
# and collapsed paralogs; depth_cap_mult is the guard on the latter. If
# the ranking looks driven by a few sites, raise eps to 0.02 and re-run
# — a call that survives both is solid.
#
# autocor_* is the independent check: it should be clearly higher for
# `best` than for `nominal` wherever the verdict is REASSIGN. If it is
# not, the likelihood is being driven by something other than genuine
# linkage structure and the call needs a closer look.
