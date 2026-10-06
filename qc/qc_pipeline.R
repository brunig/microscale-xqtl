library(ggplot2)
library(dplyr)
library(readr)

# ══════════════════════════════════════════════════════════════════════
# qc_pipeline.R — single pre-flight QC pass for xQTLpaper
# ══════════════════════════════════════════════════════════════════════
# Replaces the scattered qc_levelA.R / qc_bam_parse_plot.R / concordance
# scripts with one file. Run it after alignment and before
# xQTL_analysis4.R; it writes the allowlist that the analysis gates on.
#
# The QC/analysis boundary. Everything here is parent-agnostic: it reads
# BAM summary statistics and raw ASEReadCounter coverage, neither of
# which depends on which strain you think is in a well. That matters,
# because parent assignment is itself under revision — QC has to be
# valid regardless of how the identity question resolves. Anything that
# only exists as a consequence of a pipeline step (a count table that
# would not parse, an AFD that would not converge, a contrast with no
# matched control) is gated inside xQTL_analysis4.R and written back to
# qc/analysis_attrition.csv, which §QC5 below reads for its funnel. So
# the two files hand off in both directions and neither duplicates the
# other's logic.
#
# Order of operations for a fresh run:
#   1. qc_bam_generate.sh          -> qc/samtools_stats/*.stats
#   2. qc_pipeline.R               -> qc/qc_samples.csv, sample_allowlist.csv
#   3. scan_parent_identity.R      -> qc/parent_identity_scan_well.csv
#   4. reassign_parents.R          -> files/flatkey_reassigned.csv
#   5. xQTL_analysis4.R            -> qc/analysis_attrition.csv
#   6. qc_pipeline.R again         -> funnel + consistency now populated
# -Claude generated
# ══════════════════════════════════════════════════════════════════════


# ── QC0. Configuration ────────────────────────────────────────────────
proj_dir   = '/u/project/kruglyak/bruni/projects/xQTLpaper'
key.file   = file.path(proj_dir, 'files', 'flatkey.csv')
count_dir  = file.path(proj_dir, 'count_files')
qc_dir     = file.path(proj_dir, 'qc')
stats_dir  = file.path(qc_dir, 'samtools_stats')   # from qc_bam_generate.sh
fig_dir    = file.path(qc_dir, 'figures')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Which blocks to run. Coverage is the slow one (reads 384 count files);
# it caches to qc_samples_full.csv and is skipped when that is present
# unless named in `force`.
stages = c('bam', 'coverage', 'funnel', 'consistency')
force  = character(0)   # e.g. 'coverage'

condition_col = 'Drug3Let'

# ---- Thresholds ------------------------------------------------------
# The allowlist is deliberately permissive: its job is to remove samples
# that are not data at all, not to make analysis decisions. The real
# depth filter lives in xQTL_analysis4.R §5 and operates on segregating
# sites, which is the depth that actually matters and which cannot be
# computed here without committing to a parent assignment.
min_raw_depth   = 0.5   # mean raw depth over all counted sites
min_frac_covered = 0.10 # fraction of counted sites with >= 1 read
short_insert_bp = 50    # at or below this = adapter/primer-dimer range
low_map_pct     = 40    # flag (not drop) below this mapping rate


# ── QC1. Keyfile ──────────────────────────────────────────────────────
samples = read_csv(key.file, show_col_types = FALSE, progress = FALSE)
stopifnot(all(c('Sample_ID', 'Wells', 'Plate', condition_col) %in%
              names(samples)))
samples$condition = samples[[condition_col]]
samples$cross = paste(samples$parent_a, samples$parent_alpha, sep = 'x')
# Wells look like 'A01' — split for the plate-layout figures.
samples$well_row = substr(samples$Wells, 1, 1)
samples$well_col = as.integer(sub('^[A-Za-z]', '', samples$Wells))

cat('keyfile: ', nrow(samples), ' samples | ',
    length(unique(samples$Wells)), ' wells | conditions: ',
    paste(sort(unique(samples$condition)), collapse = ','), '\n', sep = '')


# ── QC2. BAM-level metrics from `samtools stats` ──────────────────────
# Parses the SN block for the scalar summaries and the IS block for the
# insert-size distribution. Both are stable samtools output formats.
#
# One caveat that matters for interpretation: duplicates are only
# reported if something marked them. If MarkDuplicates never ran,
# `reads duplicated` is 0 for every sample — which reads as "no
# duplicates" when it actually means "not measured". The
# dup_measured flag below keeps that distinction visible instead of
# letting a column of zeros look like a clean result.
# -Claude generated
parse_stats = function(path) {
    ln = readLines(path, warn = FALSE)
    sn = ln[startsWith(ln, 'SN')]
    getn = function(kk) {
        h = grep(paste0('\t', kk), sn, fixed = TRUE, value = TRUE)
        if (!length(h)) return(NA_real_)
        suppressWarnings(as.numeric(strsplit(h[1], '\t')[[1]][3]))
    }
    is_ln = ln[startsWith(ln, 'IS')]
    ins_med = NA_real_; frac_short = NA_real_
    if (length(is_ln)) {
        p = do.call(rbind, strsplit(is_ln, '\t'))
        sz = suppressWarnings(as.numeric(p[, 2]))
        ct = suppressWarnings(as.numeric(p[, 3]))
        okv = !is.na(sz) & !is.na(ct) & ct > 0
        sz = sz[okv]; ct = ct[okv]
        if (length(sz)) {
            cs = cumsum(ct)
            ins_med = sz[which(cs >= cs[length(cs)] / 2)[1]]
            frac_short = sum(ct[sz <= short_insert_bp]) / sum(ct)
        }
    }
    total = getn('raw total sequences:')
    mapped = getn('reads mapped:')
    dup = getn('reads duplicated:')
    paired = getn('reads properly paired:')
    data.frame(
        total_reads = total, n_mapped = mapped, n_dup = dup,
        pct_mapped = 100 * mapped / total,
        pct_proper_pair = 100 * paired / total,
        pct_dup = 100 * dup / total,
        insert_median = ins_med,
        insert_avg = getn('insert size average:'),
        insert_sd  = getn('insert size standard deviation:'),
        frac_short_inserts = frac_short,
        error_rate = getn('error rate:'),
        avg_quality = getn('average quality:'),
        stringsAsFactors = FALSE)
}

bam.file = file.path(qc_dir, 'bam_qc.csv')
if ('bam' %in% stages && dir.exists(stats_dir)) {
    sf = list.files(stats_dir, pattern = '\\.stats$', full.names = TRUE)
    cat('parsing', length(sf), 'samtools stats files\n')
    bam = bind_rows(lapply(sf, function(p) {
        d = tryCatch(parse_stats(p), error = function(e) NULL)
        if (is.null(d)) return(NULL)
        d$Sample_ID = sub('\\.stats$', '', basename(p))
        d
    }))
    if (nrow(bam)) write_csv(bam, bam.file)
} else if (file.exists(bam.file)) {
    bam = read_csv(bam.file, show_col_types = FALSE, progress = FALSE)
    cat('reusing existing', bam.file, '\n')
} else {
    bam = NULL
    cat('!! no BAM metrics — run qc_bam_generate.sh to populate ',
        stats_dir, '\n', sep = '')
}
if (!is.null(bam) && nrow(bam)) {
    bam$dup_measured = !is.na(bam$n_dup) & bam$n_dup > 0
    if (!any(bam$dup_measured))
        cat('!! duplicates unmeasured for every sample (MarkDuplicates',
            'was not run) — treat depth as an upper bound on information\n')
}


# ── QC3. Raw coverage from the count files ────────────────────────────
# Parent-agnostic: counts every site the ASEReadCounter table carries,
# with no reference to which parents were assigned. Definitions are
# stated here rather than inherited, because the old scripts' `frac_
# informative` was ambiguous:
#   n_sites        rows in the count table
#   n_covered      sites with refCount + altCount >= 1
#   frac_covered   n_covered / n_sites
#   mean_depth_raw mean(refCount + altCount) over ALL sites (zeros in)
# -Claude generated
full.file = file.path(qc_dir, 'qc_samples_full.csv')
has_fread = requireNamespace('data.table', quietly = TRUE)

cov_one = function(sid) {
    f = file.path(count_dir, paste0(sid, '.txt'))
    if (!file.exists(f) || file.info(f)$size < 1000)
        return(data.frame(Sample_ID = sid, n_sites = 0L, n_covered = 0L,
                          frac_covered = NA_real_, mean_depth_raw = NA_real_,
                          median_depth_cov = NA_real_,
                          countfile = 'missing_or_empty',
                          stringsAsFactors = FALSE))
    cnt = tryCatch({
        if (has_fread)
            data.table::fread(f, select = c('contig','position','refCount',
                              'altCount'), showProgress = FALSE,
                              data.table = FALSE)
        else as.data.frame(read_tsv(f, show_col_types = FALSE,
                 progress = FALSE))[, c('contig','position','refCount','altCount')]
    }, error = function(e) NULL)
    if (is.null(cnt) || !nrow(cnt))
        return(data.frame(Sample_ID = sid, n_sites = 0L, n_covered = 0L,
                          frac_covered = NA_real_, mean_depth_raw = NA_real_,
                          median_depth_cov = NA_real_,
                          countfile = 'unparseable',
                          stringsAsFactors = FALSE))
    n = cnt$refCount + cnt$altCount
    pc = tapply(n, cnt$contig, mean)
    out = data.frame(Sample_ID = sid, n_sites = nrow(cnt),
        n_covered = sum(n >= 1), frac_covered = mean(n >= 1),
        mean_depth_raw = mean(n),
        median_depth_cov = if (any(n >= 1)) median(n[n >= 1]) else 0,
        countfile = 'ok', stringsAsFactors = FALSE)
    attr(out, 'per_chrom') = pc
    out
}

if ('coverage' %in% stages &&
    (!file.exists(full.file) || 'coverage' %in% force)) {
    cat('reading', nrow(samples), 'count files\n')
    res = lapply(samples$Sample_ID, cov_one)
    covdf = bind_rows(res)
    # per-chromosome depth, long form, for the per-chrom figure
    pcl = bind_rows(lapply(seq_along(res), function(i) {
        pc = attr(res[[i]], 'per_chrom')
        if (is.null(pc)) return(NULL)
        data.frame(Sample_ID = samples$Sample_ID[i], chrom = names(pc),
                   mean_depth = as.numeric(pc), stringsAsFactors = FALSE)
    }))
    write_csv(covdf, full.file)
    write_csv(pcl, file.path(qc_dir, 'qc_per_chrom_depth.csv'))
} else if (file.exists(full.file)) {
    covdf = read_csv(full.file, show_col_types = FALSE, progress = FALSE)
    pcl = if (file.exists(file.path(qc_dir, 'qc_per_chrom_depth.csv')))
        read_csv(file.path(qc_dir, 'qc_per_chrom_depth.csv'),
                 show_col_types = FALSE, progress = FALSE) else NULL
    cat('reusing existing', full.file, '\n')
} else {
    stop('no coverage data and coverage stage not requested')
}


# ── QC4. Assemble, flag, and write the allowlist ──────────────────────
qc = left_join(samples, covdf, by = 'Sample_ID')
if (!is.null(bam) && nrow(bam)) qc = left_join(qc, bam, by = 'Sample_ID')

# Identity-scan verdicts, if the scan has been run. Flagged, never
# dropped — the fix for a mislabelled parent is reassignment, not
# discarding the sample.
scan.file = file.path(qc_dir, 'parent_identity_scan_well.csv')
if (file.exists(scan.file)) {
    sc = read_csv(scan.file, show_col_types = FALSE, progress = FALSE)
    sc_keep = intersect(c('Wells', 'best', 'verdict', 'p_fix_best'), names(sc))
    sc_sub  = sc[, sc_keep, drop = FALSE]
    names(sc_sub)[names(sc_sub) != 'Wells'] =
        paste0('identity_', names(sc_sub)[names(sc_sub) != 'Wells'])
    qc = left_join(qc, sc_sub, by = 'Wells')
    cat('identity verdicts joined:\n')
    print(as.data.frame(sc %>% count(verdict, name = 'n_wells')),
          row.names = FALSE)
}

# Flags built outside mutate(): whether a BAM column exists depends on
# whether qc_bam_generate.sh has run, and testing for that with `.`
# inside mutate() relies on a magrittr substitution rule that is easy to
# get subtly wrong. Explicit is safer here. -Claude generated
.has = function(col) col %in% names(qc)
qc$flag_no_countfile = qc$countfile != 'ok'
qc$flag_low_depth    = !is.na(qc$mean_depth_raw) &
                       qc$mean_depth_raw < min_raw_depth
qc$flag_low_covered  = !is.na(qc$frac_covered) &
                       qc$frac_covered < min_frac_covered
qc$flag_low_mapping  = if (.has('pct_mapped'))
    !is.na(qc$pct_mapped) & qc$pct_mapped < low_map_pct else FALSE
qc$flag_short_insert = if (.has('frac_short_inserts'))
    !is.na(qc$frac_short_inserts) & qc$frac_short_inserts > 0.5 else FALSE

qc$qc_pass = !(qc$flag_no_countfile | qc$flag_low_depth |
               qc$flag_low_covered)
qc$qc_stage = with(qc, ifelse(flag_no_countfile, 'no_countfile',
                       ifelse(flag_low_depth,   'low_depth',
                       ifelse(flag_low_covered, 'low_covered', 'pass'))))

write_csv(qc, file.path(qc_dir, 'qc_samples.csv'))
write_csv(qc %>% select(Sample_ID, qc_pass, qc_stage),
          file.path(qc_dir, 'sample_allowlist.csv'))

cat('\n── pre-flight outcome ──\n')
print(as.data.frame(qc %>% count(qc_stage, name = 'n_samples') %>%
                    arrange(desc(n_samples))), row.names = FALSE)
for (fl in grep('^flag_', names(qc), value = TRUE))
    if (any(qc[[fl]], na.rm = TRUE))
        cat(sprintf('  %-20s %d samples\n', fl, sum(qc[[fl]], na.rm = TRUE)))


# ── QC5. Figures ──────────────────────────────────────────────────────
# One multi-page PDF, so there is a single artefact to look at rather
# than a directory of numbered files. Individual PNGs go to figures/ for
# anything you want to drop into a slide.
# -Claude generated
pdf(file.path(qc_dir, 'qc_report.pdf'), width = 11, height = 8)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

# One bad panel should not abort the report and leave a half-written
# PDF with the device still open. -Claude generated
pp = function(p) tryCatch(print(p),
    error = function(e) cat('  !! panel failed:', conditionMessage(e), '\n'))

# (1) Depth waterfall — samples ranked, split by condition
p1 = qc %>% filter(!is.na(mean_depth_raw)) %>%
    arrange(mean_depth_raw) %>%
    mutate(rank = row_number()) %>%
    ggplot(aes(rank, pmax(mean_depth_raw, 0.01), color = condition)) +
    geom_point(size = 0.7) +
    geom_hline(yintercept = min_raw_depth, linetype = 2, color = 'grey40') +
    scale_y_log10() +
    labs(x = 'samples, ranked', y = 'mean raw depth (log10)',
         title = 'Coverage waterfall',
         subtitle = sprintf('dashed = allowlist floor (%.2fx)', min_raw_depth)) +
    theme_minimal(base_size = 11)
pp(p1)

# (2) Plate layout — where on the plate does coverage fail?
p2 = ggplot(qc, aes(well_col, well_row, fill = log10(pmax(mean_depth_raw, 0.01)))) +
    geom_tile(color = 'white') +
    facet_wrap(~ paste0(Plate, ' — ', condition)) +
    scale_fill_viridis_c(name = 'log10 depth', na.value = 'grey85') +
    scale_x_continuous(breaks = 1:12) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, title = 'Coverage by plate position',
         subtitle = 'edge or column structure = a handling artefact, not biology') +
    theme_minimal(base_size = 10)
pp(p2)

# (3) Covered-site fraction vs depth — separates "thin" from "empty"
p3 = ggplot(qc, aes(pmax(mean_depth_raw, 0.01), frac_covered,
                    color = condition)) +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_x_log10() +
    labs(x = 'mean raw depth (log10)', y = 'fraction of sites covered',
         title = 'Breadth vs depth',
         subtitle = 'low depth + high breadth = thin but even; low both = failed library') +
    theme_minimal(base_size = 11)
pp(p3)

# (4) Per-chromosome depth — aneuploidy and whole-chrom dropout
if (!is.null(pcl) && nrow(pcl)) {
    pc2 = pcl %>% left_join(qc %>% select(Sample_ID, condition, mean_depth_raw),
                            by = 'Sample_ID') %>%
        filter(!is.na(mean_depth_raw), mean_depth_raw > min_raw_depth) %>%
        mutate(rel = mean_depth / mean_depth_raw,
               chrom = factor(chrom, levels = paste0('chr',
                              as.character(as.roman(1:16)))))
    p4 = ggplot(pc2, aes(chrom, rel)) +
        geom_hline(yintercept = 1, color = 'grey60') +
        geom_boxplot(outlier.size = 0.4, fill = 'grey90') +
        coord_cartesian(ylim = c(0, 2.5)) +
        labs(x = NULL, y = 'chromosome depth / sample mean',
             title = 'Per-chromosome relative depth',
             subtitle = 'systematic offsets = aneuploidy or reference bias') +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    pp(p4)
}

# (5)-(7) BAM-level views
if (!is.null(bam) && nrow(bam)) {
    p5 = ggplot(qc, aes(pct_mapped, fill = condition)) +
        geom_histogram(bins = 40) +
        geom_vline(xintercept = low_map_pct, linetype = 2) +
        labs(x = '% reads mapped', y = 'samples',
             title = 'Mapping rate',
             subtitle = sprintf('dashed = flag threshold (%d%%)', low_map_pct)) +
        theme_minimal(base_size = 11)
    pp(p5)

    p6 = ggplot(qc, aes(insert_median, frac_short_inserts,
                        color = condition)) +
        geom_point(size = 0.9, alpha = 0.8) +
        labs(x = 'median insert size (bp)',
             y = sprintf('fraction of pairs <= %d bp', short_insert_bp),
             title = 'Insert size and adapter-dimer load',
             subtitle = 'top-left corner = libraries dominated by dimers') +
        theme_minimal(base_size = 11)
    pp(p6)

    p7 = ggplot(qc, aes(pct_mapped, pmax(mean_depth_raw, 0.01),
                        color = qc_stage)) +
        geom_point(size = 1) + scale_y_log10() +
        labs(x = '% reads mapped', y = 'mean raw depth (log10)',
             title = 'Depth vs mapping rate') +
        theme_minimal(base_size = 11)
    pp(p7)
}

# (8) Attrition funnel — pre-flight plus whatever the analysis reported
att.file = file.path(qc_dir, 'analysis_attrition.csv')
if ('funnel' %in% stages) {
    steps = qc %>% count(stage = qc_stage, name = 'n') %>%
        mutate(source = 'pre-flight')
    if (file.exists(att.file)) {
        att = read_csv(att.file, show_col_types = FALSE, progress = FALSE)
        steps = bind_rows(steps,
            att %>% count(stage = gate, name = 'n') %>%
                mutate(source = 'analysis'))
    } else {
        cat('\n(no ', att.file, ' yet — run xQTL_analysis4.R to populate',
            ' the analysis half of the funnel)\n', sep = '')
    }
    p8 = ggplot(steps, aes(reorder(stage, n), n, fill = source)) +
        geom_col() + coord_flip() +
        geom_text(aes(label = n), hjust = -0.2, size = 3) +
        labs(x = NULL, y = 'samples', title = 'Sample attrition',
             subtitle = 'pre-flight gates from this script; analysis gates from the ledger') +
        theme_minimal(base_size = 11)
    pp(p8)
}

# (9) Well consistency — the four conditions of a well share a parent,
# so a well where some conditions survive and others do not is telling
# you something different from a well that fails outright.
if ('consistency' %in% stages && file.exists(att.file)) {
    att = read_csv(att.file, show_col_types = FALSE, progress = FALSE)
    wellsum = att %>% group_by(Wells, cross) %>%
        summarise(n = n(), n_kept = sum(gate == 'kept'),
                  .groups = 'drop') %>%
        left_join(qc %>% group_by(Wells) %>%
                  summarise(depth_med = median(mean_depth_raw, na.rm = TRUE),
                            .groups = 'drop'), by = 'Wells') %>%
        mutate(outcome = case_when(n_kept == n ~ 'all conditions kept',
                                   n_kept == 0 ~ 'none kept',
                                   TRUE        ~ 'partial'))
    write_csv(wellsum, file.path(qc_dir, 'well_outcome_summary.csv'))

    p9 = ggplot(wellsum, aes(pmax(depth_med, 0.01), n_kept, color = outcome)) +
        geom_jitter(width = 0, height = 0.15, size = 1.6, alpha = 0.85) +
        scale_x_log10() +
        labs(x = 'median raw depth for the well (log10)',
             y = 'conditions surviving to contrasts',
             title = 'Well outcome vs coverage',
             subtitle = 'wells at good depth with nothing kept are identity suspects, not coverage failures') +
        theme_minimal(base_size = 11)
    pp(p9)

    crosssum = att %>% group_by(cross) %>%
        summarise(n = n(), n_kept = sum(gate == 'kept'),
                  pct_kept = round(100 * n_kept / n, 1), .groups = 'drop') %>%
        arrange(pct_kept)
    write_csv(crosssum, file.path(qc_dir, 'cross_outcome_summary.csv'))
    cat('\n── worst crosses by retention ──\n')
    print(as.data.frame(utils::head(crosssum, 12)), row.names = FALSE)

    # The diagnostic worth acting on: good depth, nothing kept.
    suspect = wellsum %>%
        filter(outcome == 'none kept', depth_med >= 5) %>% arrange(-depth_med)
    if (nrow(suspect)) {
        cat('\n!! wells with good coverage and no surviving condition',
            '(identity suspects):\n')
        print(as.data.frame(suspect), row.names = FALSE)
    }
}

dev.off()


# ── QC6. Report ───────────────────────────────────────────────────────
cat('\nwrote:\n')
for (f in c('qc_samples.csv', 'qc_samples_full.csv', 'sample_allowlist.csv',
            'qc_per_chrom_depth.csv', 'bam_qc.csv',
            'well_outcome_summary.csv', 'cross_outcome_summary.csv',
            'qc_report.pdf'))
    if (file.exists(file.path(qc_dir, f)))
        cat('  ', file.path(qc_dir, f), '\n', sep = '')

cat('\nallowlist: ', sum(qc$qc_pass), '/', nrow(qc), ' samples pass\n',
    sep = '')
cat('xQTL_analysis4.R will gate on qc/sample_allowlist.csv.\n')
