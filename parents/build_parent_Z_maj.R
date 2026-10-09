# ══════════════════════════════════════════════════════════════════════
# build_parent_Z.R — genotype the REAL parental haploids from their
#                    ASEReadCounter count files
# ══════════════════════════════════════════════════════════════════════
# Everything up to now picked each well's parent_a as the nearest match
# among the 919 strains in Z_1011_filtered.qs. That is a reference
# PANEL, not the strains that were actually crossed, so for some wells
# the true parent is not in the candidate set at any distance (C12's two
# candidates differ from each other by 37,966 sites, more than two
# random 1011 strains; YBB is absent from the panel entirely, which is
# why F09/G08/H10 have never been phasable). Genotyping the actual
# haploids removes the whole problem — no nearest-neighbour search, no
# od arbitration, no override rounds.
#
# ── The one decision that matters ─────────────────────────────────────
# makeCountTablesGATK2 finds segregating sites as
#     seg = (p1 != p2) & (p1 != 1) & (p2 != 1)
# so a site with NO COVERAGE must never be coded 0 — that would read as
# a confident reference call, and every site where the other parent is
# alt would become falsely segregating. That is precisely the
# "a wrong parent invents segregating sites" failure. No-call is
# therefore coded 1 (het), which the `!= 1` guards already drop. The
# analysis code needs no change.
#
# ── Output ────────────────────────────────────────────────────────────
# Z_parents.qs, same shape as Z_1011_filtered.qs ($Z, $variant_ids) on
# the IDENTICAL variant_ids, so cross.Z slicing works unchanged. Point
# Z.qs.file at it and nothing else in the pipeline moves.
# -Claude generated
# ══════════════════════════════════════════════════════════════════════

# Same env switch the rest of the pipeline uses. 'mac' runs against the
# local mirror so this can be driven from the laptop -- the panel is
# already there and the parent count files are small enough to copy,
# unlike the FASTQs. -Claude generated
env = 'mac'          # 'mac' | 'hoffman2'

paths_by_env = list(
    hoffman2 = list(
        base       = '/u/project/kruglyak/bruni/projects/xQTLpaper',
        panel.file = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/Z_1011_filtered.qs'),
    mac = list(
        base       = '/Users/kruglyakbws/Desktop/AFxQTL/hoffmanDrop/A1ref',
        panel.file = '/Users/kruglyakbws/Desktop/hoffmanDrop/Z_1011_filtered.qs')
)
stopifnot(env %in% names(paths_by_env))
for (.k in names(paths_by_env[[env]]))
    assign(.k, paths_by_env[[env]][[.k]])

par_dir     = file.path(base, 'count_files')
map.file    = file.path(base, 'files/parent_file_map.csv')  # optional: columns file,strain
out.file    = file.path(base, 'files/Z_parents_maj.qs')
rep.file    = file.path(base, 'qc/parent_genotyping_maj.csv')
dir.create(dirname(rep.file), recursive = TRUE, showWarnings = FALSE)

MIN_DEPTH   = 5      # below this a site is a no-call
HOM_FRAC    = 0.70   # alt frac >= this -> alt hom (2); <= 1-this -> ref hom (0)
                     # in between, with depth, is a real het (1) — for a
                     # haploid that means contamination or a repeat, and
                     # it is dropped by the seg-site guard either way.

suppressPackageStartupMessages({ library(Matrix); library(qs2) })

# ── 1. Which file is which strain ─────────────────────────────────────
# A mapping CSV wins; otherwise the strain is the file basename minus
# extension. The mapping is the one thing that cannot be inferred
# safely, so the script prints what it resolved and stops if anything
# looks wrong rather than silently genotyping the wrong strain.
stopifnot(dir.exists(par_dir))
files = list.files(par_dir, pattern = '\\.(txt|tsv|table|csv)(\\.gz)?$',
                   full.names = TRUE)
stopifnot(length(files) > 0)
if (file.exists(map.file)) {
    mp = read.csv(map.file, stringsAsFactors = FALSE)
    stopifnot(all(c('file','strain') %in% names(mp)))
    strain = mp$strain[match(basename(files), mp$file)]
} else {
    cat('no ', map.file, ' — taking strain from filename\n', sep = '')
    strain = sub('\\.(txt|tsv|table|csv)(\\.gz)?$', '', basename(files))
}
keep = !is.na(strain) & nzchar(strain)
if (any(!keep))
    cat('!! ', sum(!keep), ' file(s) with no strain mapping, skipped: ',
        paste(basename(files)[!keep], collapse = ', '), '\n', sep = '')
files = files[keep]; strain = strain[keep]
# Duplicate strain codes are EXPECTED here -- the A1 plate carries 2-3
# replicate wells for most strains, and they get merged in §3.
if (anyDuplicated(basename(files)))
    stop('duplicate file names: ',
         paste(unique(basename(files)[duplicated(basename(files))]), collapse = ', '))
d = table(strain)
if (any(d > 1))
    cat('  (', sum(d > 1), ' strains have replicate wells, they will be merged)\n', sep = '')
cat('resolved ', length(files), ' parent files:\n', sep = '')
print(data.frame(file = basename(files), strain = strain), row.names = FALSE)

# ── 2. Site universe, taken from the panel so the ids line up exactly ──
cat('\nloading panel for variant_ids ...\n')
P   = qs_read(panel.file)
vid = P$variant_ids
sid = sub(':.*$', '', vid)          # "chrI_27908:C/T" -> "chrI_27908"
m   = length(vid)
pos = new.env(parent = emptyenv())  # id -> index, built once
idx_of = { e = seq_along(sid); names(e) = sid; e }
cat('  ', format(m, big.mark = ','), ' variant sites\n', sep = '')

# ── 3. Genotype each parent ───────────────────────────────────────────
read_counts = function(f) {
    d = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
    need = c('contig','position','refCount','altCount')
    if (!all(need %in% names(d)))
        stop('not ASEReadCounter format (need ', paste(need, collapse = '/'),
             '): ', basename(f))
    d
}

# Replicate wells are MERGED, not picked between: 17 strains have two
# wells on the A1 plate and four have three. Summing ref/alt across them
# before calling raises the depth and so the fraction of sites that
# clear MIN_DEPTH, which is the quantity that actually limits how many
# segregating sites a cross gets.
strains = sort(unique(strain))
cat('\ngenotyping ', length(strains), ' strains from ', length(files),
    ' files\n', sep = '')

i_all = integer(0); j_all = integer(0); x_all = numeric(0)
rep = list()
for (k in seq_along(strains)) {
    st = strains[k]
    fk = files[strain == st]
    ref = numeric(m); alt = numeric(m); seen = logical(m); n_in = 0
    for (f in fk) {
        d  = read_counts(f); n_in = n_in + nrow(d)
        j  = idx_of[paste0(d$contig, '_', d$position)]
        ok = !is.na(j); j = j[ok]
        ref[j] = ref[j] + d$refCount[ok]
        alt[j] = alt[j] + d$altCount[ok]
        seen[j] = TRUE
    }
    tot = ref + alt
    af  = ifelse(tot > 0, alt / tot, NA_real_)
    g = rep_len(NA_integer_, m)
    g[tot >= MIN_DEPTH & !is.na(af) & af >= HOM_FRAC]       = 2L
    g[tot >= MIN_DEPTH & !is.na(af) & af <= (1 - HOM_FRAC)] = 0L
    g[tot >= MIN_DEPTH & is.na(g)]                          = 1L  # real het
    called = !is.na(g)
    # No-call -> 1 (het) so the seg-site guard drops it. See header.
    g[!called] = 1L

    nz = g != 0L                    # 0 is the sparse default
    i_all = c(i_all, rep_len(k, sum(nz)))
    j_all = c(j_all, which(nz)); x_all = c(x_all, g[nz])

    rep[[k]] = data.frame(strain = st, file = paste(basename(fk), collapse = ';'),
        n_files = length(fk), n_in_file = n_in, n_matched = sum(seen),
        n_called = sum(called), pct_called = round(100 * sum(called) / m, 1),
        n_ref = sum(g == 0L & called), n_alt = sum(g == 2L & called),
        n_het = sum(g == 1L & called),
        pct_het_of_called = round(100 * sum(g == 1L & called) /
                                  max(1, sum(called)), 2),
        median_depth = if (any(tot > 0)) median(tot[tot > 0]) else 0)
    cat(sprintf('  %-9s %d file(s)  matched %6.1f%%  called %6.1f%%  het %5.2f%%  med depth %4.0f\n',
        st, length(fk), 100 * sum(seen) / m, rep[[k]]$pct_called,
        rep[[k]]$pct_het_of_called, rep[[k]]$median_depth))
}
strain_rows = strains
R = do.call(rbind, rep)

Zp = sparseMatrix(i = i_all, j = j_all, x = x_all,
                  dims = c(length(strain_rows), m),
                  dimnames = list(strain_rows, NULL))

# ── 4. Sanity check against the panel ─────────────────────────────────
# Any strain code that also exists in the panel should agree with it at
# nearly every confidently-called site. This validates the file->strain
# mapping and the genotyping in one shot: a strain that disagrees with
# its own panel row means the mapping is wrong, not that the panel is.
cat('\n── concordance with the 1011 panel (shared strain codes) ──\n')
shared = intersect(strain_rows, rownames(P$Z))
R$panel_concord = NA_real_; R$panel_n = NA_real_
if (length(shared)) {
    for (s in shared) {
        a = as.numeric(Zp[s, ]); b = as.numeric(P$Z[s, ])
        u = a != 1 & b != 1                       # both confidently called
        R$panel_concord[R$strain == s] = round(100 * mean(a[u] == b[u]), 2)
        R$panel_n[R$strain == s] = sum(u)
        cat(sprintf('  %-8s %6.2f%% over %s sites\n', s,
            round(100 * mean(a[u] == b[u]), 2), format(sum(u), big.mark = ',')))
    }
    bad = R$strain[!is.na(R$panel_concord) & R$panel_concord < 95]
    if (length(bad))
        cat('!! below 95% — check the mapping for: ',
            paste(bad, collapse = ', '), '\n', sep = '')
} else {
    cat('  (no strain codes in common — cannot cross-check the mapping)\n')
}

# ── 5. Carry over any parent we could not sequence, from the panel ────
# parent_alpha (BBQ) and any haploid without a count file still have to
# come from somewhere. Taking them from the panel keeps the object
# complete; the report says which rows are real and which are panel.
missing = setdiff('BBQ', rownames(Zp))
carry = intersect(missing, rownames(P$Z))
if (length(carry)) {
    cat('\ncarrying ', length(carry), ' row(s) from the panel: ',
        paste(carry, collapse = ', '), '\n', sep = '')
    Zp = rbind(Zp, P$Z[carry, , drop = FALSE])
}
R$source = 'sequenced'
if (length(carry))
    R = rbind(R, data.frame(strain = carry, file = NA, n_files = NA,
        n_in_file = NA, n_matched = NA, n_called = NA, pct_called = NA,
        n_ref = NA, n_alt = NA, n_het = NA, pct_het_of_called = NA,
        median_depth = NA, panel_concord = NA, panel_n = NA, source = 'panel'))

out = list(Z = Zp, variant_ids = vid)
qs_save(out, out.file)
write.csv(R, rep.file, row.names = FALSE)
cat('\nwrote ', out.file, ' (', nrow(Zp), ' strains x ',
    format(m, big.mark = ','), ' variants)\n', sep = '')
cat('wrote ', rep.file, '\n', sep = '')
cat('\nNext: point Z.qs.file at Z_parents.qs, clear the cache, re-run.\n')
