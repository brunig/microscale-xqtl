# ══════════════════════════════════════════════════════════════════
# Three-way parent-identity comparison for xQTLpaper
#
# For each sample: compare
#   (a) the parent_a code in the current keyfile,
#   (b) the Hamming-verified code from previously sequenced reference
#       haploids (via strain common name lookup), and
#   (c) the concordance scan's best_p1 (from the read-count scan).
#
# Reads:
#   /u/project/kruglyak/bruni/projects/xQTLpaper/files/flatkey.csv
#   /u/project/kruglyak/bruni/projects/xQTLpaper/files/hamming_reference_haploids.csv
#   /u/project/kruglyak/bruni/projects/xQTLpaper/qc/cross_identity_scan.csv  (optional)
#
# Writes:
#   qc/three_way_identity_check.csv
# ══════════════════════════════════════════════════════════════════

library(dplyr)
library(readr)

# ---- paths --------------------------------------------------------
if (!exists('count_dir'))
    count_dir = '/u/project/kruglyak/bruni/projects/xQTLpaper/count_files/'
if (!exists('key.file'))
    key.file  = '/u/project/kruglyak/bruni/projects/xQTLpaper/files/flatkey.csv'
qc_dir       = file.path(dirname(count_dir), 'qc')
hamming_file = file.path(dirname(key.file),
                         'hamming_reference_haploids.csv')
scan_file    = file.path(qc_dir, 'cross_identity_scan.csv')


# ---- load ---------------------------------------------------------
samples = readr::read_csv(key.file, show_col_types = FALSE)

if (!file.exists(hamming_file))
    stop('hamming_reference_haploids.csv not found at ', hamming_file,
         '\nRun hamming_export.py locally and scp the CSV to Hoffman2.')
ham = readr::read_csv(hamming_file, show_col_types = FALSE)


# ---- build strain -> hamming_code lookup --------------------------
# The reference haploid rows have (strain_common_name, postHam,
# standardized_name, top-3 hamming distances). Some strains have
# multiple replicate rows in the Hamming file; consolidate by
# taking the modal call.

# Prefer standardized_name; fall back to postHam
ham$hamming_code = ifelse(!is.na(ham$standardized_name) &
                           ham$standardized_name != 'NA',
                          ham$standardized_name, ham$postHam)
# Strip trailing '?' from uncertain calls
ham$hamming_code = sub('\\?$', '', ham$hamming_code)

lookup = ham |>
    dplyr::filter(!is.na(strain_common_name)) |>
    dplyr::group_by(strain_common_name) |>
    dplyr::summarise(
        hamming_code    = names(sort(table(hamming_code), decreasing = TRUE))[1],
        n_ref_replicates = dplyr::n(),
        top1_hamming    = suppressWarnings(as.integer(min(top1_hamming, na.rm = TRUE))),
        gap_to_top2     = {
            v1 = suppressWarnings(as.integer(top1_hamming))
            v2 = suppressWarnings(as.integer(top2_hamming))
            if (any(!is.na(v1)) && any(!is.na(v2)))
                as.integer(min(v2 - v1, na.rm = TRUE))
            else NA_integer_
        },
        .groups = 'drop'
    )
cat('Hamming lookup built:', nrow(lookup),
    'unique strains with reference genotypes\n')


# ---- join to xQTLpaper keyfile -----------------------------------
merged = samples |>
    dplyr::left_join(lookup,
                     by = c('oldLabel_parent_a' = 'strain_common_name'))


# ---- optional: bring in concordance scan --------------------------
have_scan = file.exists(scan_file)
if (have_scan) {
    scan = readr::read_csv(scan_file, show_col_types = FALSE)
    merged = merged |>
        dplyr::left_join(scan |>
                         dplyr::select(Sample_ID,
                                       scan_best_p1     = best_p1,
                                       scan_best_conc   = best_conc,
                                       scan_assigned_rank = assigned_rank,
                                       scan_diff        = diff,
                                       scan_depth_flag  = depth_flag),
                         by = 'Sample_ID')
    cat('Concordance scan joined:', sum(!is.na(merged$scan_best_p1)),
        'samples\n')
} else {
    cat('Concordance scan output not found (', scan_file,
        ') - skipping scan join\n')
    merged$scan_best_p1     = NA_character_
    merged$scan_assigned_rank = NA_integer_
    merged$scan_depth_flag  = NA_character_
}


# ---- three-way agreement classification --------------------------
classify = function(keyfile_p1, hamming_code, scan_best_p1, scan_rank,
                    depth_flag) {
    has_ham  = !is.na(hamming_code) && hamming_code != 'NA'
    has_scan = !is.na(scan_best_p1) &&
               !is.na(scan_rank) &&
               depth_flag %in% c('high_confidence', 'medium_confidence')

    if (!has_ham && !has_scan)     return('no_reference')
    if (!has_ham &&  has_scan)     return(
        if (scan_rank == 1) 'scan_agrees_only' else 'scan_disagrees_only')
    if ( has_ham && !has_scan)     return(
        if (keyfile_p1 == hamming_code) 'hamming_agrees_only'
        else                            'hamming_disagrees_only')

    # Both available
    key_eq_ham  = keyfile_p1  == hamming_code
    scan_eq_ham = scan_best_p1 == hamming_code
    key_eq_scan = keyfile_p1  == scan_best_p1

    if (key_eq_ham && key_eq_scan)                  return('all_agree')
    if (!key_eq_ham && scan_eq_ham)                 return('mislabel_confirmed')
    if (!key_eq_ham && !scan_eq_ham && key_eq_scan)
        return('hamming_disagrees_scan_supports_keyfile')
    if (!key_eq_ham && !key_eq_scan && !scan_eq_ham)
        return('three_way_disagree')
    if (key_eq_ham && !scan_eq_ham)
        return('scan_disagrees_but_keyfile_matches_hamming')
    'unclassified'
}

merged$agreement = mapply(classify,
    merged$parent_a, merged$hamming_code,
    merged$scan_best_p1, merged$scan_assigned_rank,
    merged$scan_depth_flag)


# ---- write output -------------------------------------------------
out_cols = c(
    'Sample_ID', 'Plate', 'Wells', 'Drug3Let',
    'oldLabel_parent_a', 'parent_a',
    'hamming_code', 'n_ref_replicates',
    'top1_hamming', 'gap_to_top2',
    'scan_best_p1', 'scan_assigned_rank',
    'scan_depth_flag',
    'agreement'
)
out_cols = intersect(out_cols, colnames(merged))
out = merged[, out_cols]

# Sort: interesting cases first (mislabels), then all_agree, then no_ref
priority = c(
    'mislabel_confirmed'                     = 1,
    'three_way_disagree'                     = 2,
    'hamming_disagrees_only'                 = 3,
    'scan_disagrees_only'                    = 4,
    'hamming_disagrees_scan_supports_keyfile' = 5,
    'scan_disagrees_but_keyfile_matches_hamming' = 6,
    'hamming_agrees_only'                    = 7,
    'scan_agrees_only'                       = 8,
    'all_agree'                              = 9,
    'no_reference'                           = 10
)
out$.priority = priority[out$agreement]
out = out[base::order(out$.priority, out$Sample_ID), ]
out$.priority = NULL

write.csv(out, file.path(qc_dir, 'three_way_identity_check.csv'),
          row.names = FALSE)


# ---- console report -----------------------------------------------
cat('\n==== Three-way identity check ====\n')
cat('\nAgreement classification:\n')
print(table(merged$agreement))

cat('\nMislabels confirmed by BOTH Hamming and scan (highest priority):\n')
mc = out[out$agreement == 'mislabel_confirmed', ]
if (nrow(mc) > 0) {
    print(as.data.frame(mc[, c('Sample_ID', 'Wells', 'Drug3Let',
                               'oldLabel_parent_a', 'parent_a',
                               'hamming_code', 'scan_best_p1',
                               'scan_depth_flag')]))
} else {
    cat('  (none)\n')
}

cat('\nHamming disagrees with keyfile (scan not conclusive):\n')
hd = out[out$agreement == 'hamming_disagrees_only', ]
if (nrow(hd) > 0) {
    print(as.data.frame(hd[, c('Sample_ID', 'Wells', 'Drug3Let',
                               'oldLabel_parent_a', 'parent_a',
                               'hamming_code', 'top1_hamming',
                               'gap_to_top2')]))
} else {
    cat('  (none)\n')
}

cat('\nThree-way disagreements (all three sources call something different):\n')
td = out[out$agreement == 'three_way_disagree', ]
if (nrow(td) > 0) {
    print(as.data.frame(td[, c('Sample_ID', 'Wells', 'Drug3Let',
                               'oldLabel_parent_a', 'parent_a',
                               'hamming_code', 'scan_best_p1',
                               'scan_depth_flag')]))
} else {
    cat('  (none)\n')
}

cat('\nSamples with no reference in either Hamming or scan:\n')
nr = out[out$agreement == 'no_reference', ]
if (nrow(nr) > 0) {
    cat('  n =', nrow(nr), '  unique strains:',
        paste(sort(unique(nr$oldLabel_parent_a)), collapse = ', '), '\n')
}

cat('\nWritten to:\n  ',
    file.path(qc_dir, 'three_way_identity_check.csv'), '\n')
