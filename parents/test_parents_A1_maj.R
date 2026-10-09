# ══════════════════════════════════════════════════════════════════════
# test_parents_A1.R — nominal vs reassigned, judged on the actual reads
# ══════════════════════════════════════════════════════════════════════
# Every parent decision so far was made by matching against the 1011
# PANEL. We now have the physically sequenced haploids, so the question
# stops being "which of 919 reference strains is closest" and becomes
# the direct one: for this well's reads, does the A1-sequenced nominal
# or the panel-derived reassignment explain them better?
#
# Two statistics, deliberately different in kind:
#   d_half  mean |AF - 0.5| over segregating sites. A wrong parent
#           invents segregating sites where the pool carries only BBQ,
#           so AF collapses toward 0 and d_half goes large. This is the
#           global-offset signature seen on the fitness plots.
#   od      residual scatter around a local smooth over the binomial
#           expectation (same statistic as recheck_held_wells.R Part B).
#           Catches a related-but-wrong parent that keeps AF near 0.5
#           but misphases individual sites.
# They fail differently, so agreement between them is meaningful and
# disagreement is a flag rather than noise.
#
# n_seg is reported alongside because it is evidence in its own right:
# a wrong parent INVENTS segregating sites (correcting G11 dropped its
# count 50,822 -> 15,845), so a candidate with far more sites than its
# rival is suspect even when the fit statistics are close.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(Matrix) })

base      = '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop'
count_dir = file.path(base, 'notAF/count_files')
z.file    = file.path(base, 'A1ref/files/Z_hybrid_maj.qs')
nom.file  = file.path(base, 'notAF/files/flatkey.csv')            # nominal
cur.file  = file.path(base, 'notAF/files/flatkey_override2.csv')  # after rounds 1-2
out.file  = file.path(base, 'A1ref/qc/parent_test_A1_maj.csv')
MIN_DEPTH = 3      # per-site, to enter the statistics
WIN       = 501    # sites in the local smooth

Z = qs2::qs_read(z.file); zr = rownames(Z$Z)
sid_v = sub(':.*$', '', Z$variant_ids)
chr_v = sub('_[0-9]+$', '', sid_v)
pos_v = as.integer(sub('^.*_', '', sid_v))

tw = function(x) trimws(as.character(x))
nk = read.csv(nom.file, check.names = FALSE); ck = read.csv(cur.file, check.names = FALSE)
key = data.frame(Sample_ID = tw(ck$Sample_ID), Wells = tw(ck$Wells),
                 alpha = tw(ck$parent_alpha), current = tw(ck$parent_a),
                 stringsAsFactors = FALSE)
key$nominal = tw(nk$parent_a)[match(key$Sample_ID, tw(nk$Sample_ID))]

geno = local({ e = new.env(parent = emptyenv())
    function(s) { if (is.null(e[[s]])) e[[s]] = as.numeric(Z$Z[s, ]); e[[s]] } })

read_counts = function(sid) {
    f = file.path(count_dir, paste0(sid, '.txt'))
    if (!file.exists(f)) return(NULL)
    d = try(read.delim(f, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(d, 'try-error') || !nrow(d)) return(NULL)
    j = match(paste0(d$contig, '_', d$position), sid_v)
    ok = !is.na(j)
    list(j = j[ok], ref = d$refCount[ok], alt = d$altCount[ok])
}

# od: variance of (observed - local mean) over the binomial expectation.
od_of = function(af, n, ch, po) {
    num = c(); den = c()
    for (cc in unique(ch)) {
        k = which(ch == cc); k = k[order(po[k])]
        if (length(k) < 2 * WIN) next
        loc = as.numeric(stats::filter(af[k], rep(1/WIN, WIN), sides = 2))
        keep = !is.na(loc) & loc > 0 & loc < 1
        if (!any(keep)) next
        res = af[k][keep] - loc[keep]
        num = c(num, res^2); den = c(den, loc[keep]*(1-loc[keep])/n[k][keep])
    }
    if (!length(num)) return(NA_real_)
    mean(num) / mean(den)
}

score = function(cn, cand, alpha) {
    if (!(cand %in% zr) || !(alpha %in% zr))
        return(c(n_seg = NA, d_half = NA, med_af = NA, od = NA, depth = NA))
    a = geno(cand); b = geno(alpha)
    seg = a != b & a != 1 & b != 1          # no-call is coded 1 and drops out
    s = seg[cn$j]
    if (!any(s)) return(c(n_seg = 0, d_half = NA, med_af = NA, od = NA, depth = NA))
    j = cn$j[s]; r = cn$ref[s]; al = cn$alt[s]
    tot = r + al
    keep = tot >= MIN_DEPTH
    if (sum(keep) < 5 * WIN)
        return(c(n_seg = sum(s), d_half = NA, med_af = NA, od = NA,
                 depth = if (any(keep)) mean(tot[keep]) else 0))
    j = j[keep]; tot = tot[keep]
    # phase: parent_a carries ref where its genotype is 0, else alt
    p1 = ifelse(a[j] == 0, r[keep], al[keep])
    af = p1 / tot
    o = order(chr_v[j], pos_v[j])
    c(n_seg = sum(s), d_half = mean(abs(af - 0.5)), med_af = median(af),
      od = od_of(af[o], tot[o], chr_v[j][o], pos_v[j][o]), depth = mean(tot))
}

wells = unique(key$Wells)
cat('testing ', length(wells), ' wells\n', sep = '')
res = list()
for (i in seq_along(wells)) {
    w = wells[i]
    sub = key[key$Wells == w, , drop = FALSE]
    best = NULL; bn = -1
    for (sid in sub$Sample_ID) {
        cn = read_counts(sid); if (is.null(cn)) next
        tt = sum(cn$ref + cn$alt)
        if (tt > bn) { bn = tt; best = list(sid = sid, cn = cn) }
    }
    if (is.null(best)) next
    nom = sub$nominal[1]; cur = sub$current[1]; alp = sub$alpha[1]
    sn = score(best$cn, nom, alp); sc = score(best$cn, cur, alp)
    res[[length(res)+1]] = data.frame(Wells = w, sample = best$sid,
        nominal = nom, current = cur, same = nom == cur,
        n_seg_nom = sn[['n_seg']], n_seg_cur = sc[['n_seg']],
        dhalf_nom = sn[['d_half']], dhalf_cur = sc[['d_half']],
        medaf_nom = sn[['med_af']], medaf_cur = sc[['med_af']],
        od_nom = sn[['od']], od_cur = sc[['od']], depth = sn[['depth']])
    cat(sprintf('  %-4s %-4s vs %-4s  d_half %5.3f / %5.3f   od %6.2f / %6.2f\n',
        w, nom, cur, sn[['d_half']], sc[['d_half']], sn[['od']], sc[['od']]))
}
R = do.call(rbind, res)
R$winner = ifelse(is.na(R$dhalf_nom) & is.na(R$dhalf_cur), 'no_data',
           ifelse(is.na(R$dhalf_cur) | (!is.na(R$dhalf_nom) & R$dhalf_nom < R$dhalf_cur),
                  'nominal', 'current'))
R$winner[R$same] = 'same'
write.csv(R, out.file, row.names = FALSE)
cat('\n── verdict by d_half (lower = fits the reads better) ──\n')
print(table(R$winner))
cat('\nwrote ', out.file, '\n', sep = '')
