# ══════════════════════════════════════════════════════════════════════
# reassign_parents.R — parent_a name reassignment for xQTLpaper
# ══════════════════════════════════════════════════════════════════════
# Purpose: reconcile the `parent_a` codes in flatkey.csv against the
# Hamming-distance identity calls made on the individually sequenced
# haploid parents (Box A1 reference run), optionally cross-checked
# against the read-count concordance scan, and emit a corrected keyfile.
#
# Run this BEFORE §1 of xQTL_analysis3.R. Reassigning parent_a changes
# which rows of Z are used to phase a sample, so §1 (cross.Z) and §2
# (makeCountTablesGATK2) must be re-run for every affected sample.
# Alignment / count files are unaffected — this is a relabel, not a
# realignment.
#
# Direction of evidence (important, and the opposite of a first read):
#   `Standardized_name` in the Hamming sheet is the NOMINAL code — the
#   1011 code the strain was *supposed* to be. `col1Name` is what the
#   Hamming scan actually FOUND. So col1Name is the correction, and
#   Standardized_name is what flatkey.csv currently carries. That is
#   confirmed by the earlier run's own relabel columns in the same
#   sheet, which go nominal -> col1Name (ADE->ACS, ADS->ADG, ...).
# -Claude generated
# ══════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})


# ── R0. Configuration ─────────────────────────────────────────────────
proj_dir   = '/u/project/kruglyak/bruni/projects/xQTLpaper'
key.file   = file.path(proj_dir, 'files', 'flatkey.csv')
out.key    = file.path(proj_dir, 'files', 'flatkey_reassigned.csv')
qc_dir     = file.path(proj_dir, 'qc')
out.log    = file.path(qc_dir, 'parent_reassignment_log.csv')

# Optional: read-count concordance / identity scan output from the QC
# work. If the file is absent the script runs on Hamming evidence only.
# Well-level output of scan_parent_identity.R — the richer of its two
# per-sample/per-well files, and the one this script reasons from.
scan.well.file = file.path(qc_dir, 'parent_identity_scan_well.csv')

# Optional: Z matrix, used only to check that proposed codes exist as
# rownames. Set to NA to skip (the check is cheap but loads ~GB).
Z.qs.file  = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/Z_1011_filtered.qs'
check_Z    = TRUE

# ---- What actually gets applied -------------------------------------
# Which evidence tiers to act on is set in §R5 (`apply_tiers`), where
# the tiers are defined. Default is 1-3: apply only where the pool scan
# and a second independent line of evidence agree.
drop_unmapped   = FALSE   # drop samples whose final code is not in Z

# Classification thresholds. Typical clean self-match in this reference
# run is d1 ~ 400-900; a strain whose own best match is in the
# thousands was poorly sequenced and its call means nothing.
d1_clean     = 1500   # <= this: reference sequenced well enough to trust
d1_unusable  = 6000   # >  this: reference is junk, make no call
ratio_clear  = 2.0    # d2/d1 needed to call the top hit unambiguous
cluster_frac = 1.25   # hits within d1 * this are "the same cluster"


# ── R1. Hamming reference: nominal code vs. top-10 observed matches ──
# One row per reference haploid strain, keyed by the common strain name
# used in flatkey.csv$oldLabel_parent_a. Extracted from
# hamming_resultsLongerShuffle_relabelled.csv (reference-haploid rows
# only, i.e. STRAIN != EMPTY). Empty distance fields = the sheet had no
# Hamming result for that strain; its `Standardized_name` was asserted
# with a "?" and carries no supporting evidence. -Claude generated
ham_csv = '
strain,nominal_code,ref_sample,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
CBS2888,ABL,Sample027,ABL,YCB,CCP,AMB,AMC,CDQ,CDB,BEG,AVI,ANK,826,994,46780,48160,48220,48223,48289,48332,48397,48411
CECT10109,ABP,Sample028,ABP,BIR,CPF,AQB,BRH,BRD,BLI,BLG,BQE,CHF,527,11168,11874,12703,12714,12942,13153,13219,13235,13276
CECT10266,AVB,Sample037,,,,,,,,,,,,,,,,,,,,
CLIB192,AAM,Sample001,AAM,AAL,ATR,BCS,BCT,BTF,AMK,APN,ASL,AMG,18642,23362,27830,28308,28556,28590,28662,28721,28786,28819
CLIB219,ACQ,Sample065,ADR,BRD,BSC,BQE,BLI,BLG,AMS,YCP,AIL,CBG,736,6952,7026,7204,7303,7353,7418,7527,7530,7552
CLIB294,ACF,Sample038,ACF,CBG,BRD,BQE,BLI,BLG,ALL,ALM,AIF,AQR,408,6309,6380,6525,6576,6683,6729,6733,7040,7087
CLIB382,AAR,Sample006,AAR,YAG,AEA,BRM,AEQ,YAB,AFB,CHF,BRD,BSL,606,8819,9655,9667,9815,14580,24127,24360,24573,24710
CLIB413,ACG,Sample014,ACK,ADH,ABB,YDA,AKL,ACM,AQS,BHS,CEK,BHP,623,645,648,5447,6573,7409,14721,14916,14960,14990
DBVPG3591,ABS,Sample031,ABS,BRD,BQE,BLI,CBG,BLG,ALL,BKD,BNP,BQH,462,6776,7344,7401,7446,7455,7460,7473,7476,7538
DBVPG4651,ABT,Sample005,ADM,CHG,AKT,BML,CHR,CKS,BMM,BTE,CKK,CBF,16644,26963,27272,28169,28396,29229,29459,29476,29525,29543
IL_01,ADD,Sample063,ACR,BKC,BIH,BLS,BRD,BKD,BNP,BQE,BLI,CBG,698,2914,2915,3151,7599,7715,7774,7855,7983,8008
K1,ADN,Sample041,ACH,ACG,CHV,BRD,CHF,CAI,BSL,BQE,BLI,CAL,5060,5077,10200,19625,19689,19760,19843,19889,19906,19947
M22,ADR,Sample067,ACT,AKV,ACV,AKP,ALE,BRD,BQE,BLI,CBG,BLG,437,462,463,464,465,7420,7658,7670,7725,7743
PW5,ADE,Sample064,ACS,AMA,BFB,CHD,BEP,BBI,BBH,ALS,BBK,ALT,578,22190,22194,22801,22948,23342,23468,23469,23477,23480
T73,ADS,Sample078,ADG,YDI,CDP,CMS,YCO,CDQ,BBR,BEG,CDB,ANF,766,778,836,868,6734,16547,18731,25748,28218,28229
UC1,ADC,Sample062,ADC,YCP,BSB,BRD,BQE,BLI,BLG,CBG,AMS,ALL,2893,2896,7789,8377,8602,8604,8648,8765,8818,8821
UC8,ACT,Sample091,,,,,,,,,,,,,,,,,,,,
WE372,ACV,Sample079,AMB,AMC,CDQ,BFA,BFL,YCW,BEG,CDB,ACD,ANF,816,822,13600,14404,14428,21163,28073,28379,28966,28970
Y10,ACI,Sample042,,,,,,,,,,,,,,,,,,,,
Y12,ACK,Sample043,ACK,ADH,ABB,YDA,AKL,ACM,AQS,BHS,CEK,BHP,659,689,703,5570,6641,7487,14827,15009,15053,15093
Y3,ACL,Sample044,ACL,AHH,AFK,CDA,CCV,CCT,ARG,AVV,BAC,YDE,2160,32024,34692,37806,39011,39364,41754,41782,41802,41812
Y9J,ADB,Sample061,,,,,,,,,,,,,,,,,,,,
YJM145,AAS,Sample007,AAS,YBJ,CHF,BSL,YBW,BRD,CAI,BLI,BQE,CAL,571,36118,36134,36329,36345,36489,36511,36557,36656,36656
YJM269,ACN,Sample045,ACN,ADT,S288C,CBL,AFD,BCR,AQH,CKK,CBM,CBF,22724,23471,23489,25625,25664,26278,26468,26568,26594,26681
YJM280,AAT,Sample008,AAT,YBI,BRV,BSK,CFP,YBW,AMN,CHF,BMD,CAI,742,36935,38885,39581,39778,39848,39910,39966,40036,40127
YJM320,AAV,Sample009,AAV,YAH,BRV,YAP,YAQ,YAR,CHF,ABA,YAI,CAI,737,955,12473,18877,21226,35357,38186,38280,38302,38752
YJM326,ABA,Sample013,,,,,,,,,,,,,,,,,,,,
YJM421,ABC,Sample015,YDL,ABC,YBM,CIL,CRE,CIG,CKN,BBT,AKL,CHN,625,666,17982,38564,38861,39057,39523,40221,40243,40397
YJM428,ABD,Sample017,ABD,YAJ,ATC,CFP,BMD,AMN,AMG,APK,CKK,ATB,625,819,17920,38450,38473,38625,38846,38856,38863,38887
YJM434,ABE,Sample018,,,,,,,,,,,,,,,,,,,,
YJM436,ABF,Sample019,ABF,YAL,ATC,CMG,ATB,CRE,CKN,BPA,BIE,BBT,7772,7778,41776,41783,41839,41860,41918,42028,42040,42143
YJM440,ABG,Sample004,ADL,BDL,CHR,BFE,BFG,BFM,BTD,APL,ADM,AKT,17852,25551,26644,27030,27048,27088,27167,28470,28492,28573
YJM454,ABH,Sample021,ABH,YAN,AMB,AMC,BFA,BFL,ALG,AKS,CDQ,CKN,647,10352,35885,35969,36996,37006,39997,40712,41051,41445
YJM653,ABI,Sample025,,,,,,,,,,,,,,,,,,,,
YJM678,ABK,Sample026,ABK,YAU,ABA,YAI,YAV,BRV,CHF,CAI,CAL,CBF,640,18520,22556,22587,24427,35628,38088,38418,38536,38578
YJM978,SACE_YBB,Sample069,,,,,,,,,,,,,,,,,,,,
YJM981,ADI,Sample029,ABQ,BNV,ABR,BPB,ALM,BQE,ALL,CBG,BRD,BLI,407,424,426,5758,6628,6629,6638,6648,6651,7042
YPS163,ACD,Sample033,ACD,AVI,AND,ANF,ANK,ANI,YCU,ANE,YBS,ANH,750,751,767,767,786,793,793,807,808,811
'
ham = read_csv(I(ham_csv), show_col_types = FALSE, progress = FALSE)

# Independent corroboration: the relabel decisions the earlier
# Fluconazole run actually shipped (parent_a_oldLabel -> parent_a in
# the cross-sample rows of the same sheet, with parent_a_alt). Kept
# separate from the Hamming distances because a human curated these.
# -Claude generated
prior_csv = '
nominal_code,prior_call,prior_alt
ADE,ACS,ACS
ACQ,ADF,ADR/ACK
ADS,ADG,ADG
ADI,ABQ,ACT
ACV,ACV,AMB
'
prior = read_csv(I(prior_csv), show_col_types = FALSE, progress = FALSE)


# ── R2. Classify each reference strain ────────────────────────────────
hcols = paste0('h', 1:10)
dcols = paste0('d', 1:10)

classify_strain = function(row) {
    nm  = row$nominal_code
    hs  = unlist(row[hcols], use.names = FALSE)
    ds  = suppressWarnings(as.numeric(unlist(row[dcols], use.names = FALSE)))

    if (all(is.na(ds))) {
        return(list(class = 'no_hamming_data', call = NA_character_,
                    d1 = NA_real_, ratio = NA_real_, cluster = NA_character_))
    }
    d1    = ds[1]
    d2    = if (!is.na(ds[2])) ds[2] else d1 * 999
    ratio = d2 / d1
    # near-clonal cluster: everything within cluster_frac of the best hit
    in_cl = !is.na(ds) & ds <= d1 * cluster_frac
    clust = hs[in_cl]
    clust = clust[!is.na(clust) & nzchar(clust)]

    cls =
        if (d1 > d1_unusable)                 'unreliable_reference'
        else if (identical(hs[1], nm))        if (ratio >= ratio_clear) 'agree' else 'agree_in_cluster'
        else if (nm %in% clust)               'agree_in_cluster'
        else if (ratio >= ratio_clear && d1 <= d1_clean) 'mislabel_confident'
        else                                  'mislabel_ambiguous'

    list(class = cls, call = hs[1], d1 = d1, ratio = ratio,
         cluster = paste(clust, collapse = '/'))
}

ham_cls = bind_rows(lapply(seq_len(nrow(ham)), function(i) {
    r = classify_strain(ham[i, ])
    tibble(strain       = ham$strain[i],
           nominal_code = ham$nominal_code[i],
           ref_sample   = ham$ref_sample[i],
           ham_call     = r$call,
           ham_class    = r$class,
           ham_d1       = r$d1,
           ham_ratio    = r$ratio,
           ham_cluster  = r$cluster)
})) %>% left_join(prior, by = 'nominal_code')

cat('\n── Hamming reference classification (', nrow(ham_cls), 'strains) ──\n')
print(as.data.frame(ham_cls %>% count(ham_class, name = 'n_strains')), row.names = FALSE)
cat('\n')
print(as.data.frame(ham_cls %>%
    filter(ham_class %in% c('mislabel_confident', 'mislabel_ambiguous',
                            'unreliable_reference', 'no_hamming_data')) %>%
    arrange(factor(ham_class, levels = c('mislabel_confident',
                                         'mislabel_ambiguous',
                                         'unreliable_reference',
                                         'no_hamming_data'))) %>%
    select(strain, nominal_code, ham_call, ham_class, ham_d1, ham_ratio,
           ham_cluster, prior_call, prior_alt)), row.names = FALSE)


# ── R3. Join to the keyfile ───────────────────────────────────────────
samples = read_csv(key.file, show_col_types = FALSE, progress = FALSE)
stopifnot(all(c('Sample_ID', 'oldLabel_parent_a', 'parent_a',
                'parent_alpha') %in% names(samples)))

rea = samples %>%
    mutate(strain = trimws(oldLabel_parent_a),
           parent_a_orig = trimws(parent_a)) %>%
    left_join(ham_cls, by = 'strain') %>%
    mutate(ham_class = ifelse(is.na(ham_class), 'no_reference_row', ham_class))

# Sanity: does the keyfile's parent_a match the sheet's nominal code?
# It should, since the keyfile was built from Standardized_name. Any
# disagreement here is a keyfile bookkeeping error, separate from the
# identity question.
nominal_mismatch = rea %>%
    filter(!is.na(nominal_code), parent_a_orig != nominal_code) %>%
    count(strain, parent_a_orig, nominal_code, name = 'n_samples')
if (nrow(nominal_mismatch)) {
    cat('\n!! Keyfile parent_a disagrees with the sheet nominal code',
        '(bookkeeping, not identity):\n')
    print(as.data.frame(nominal_mismatch), row.names = FALSE)
}


# ── R4. Identity scan, joined per well ────────────────────────────────
# Read the well-level output of scan_parent_identity.R. This is the
# stronger evidence of the two: it is measured on the actual xQTLpaper
# pools, whereas Hamming saw Box A1 as it was at sequencing time and is
# blind to anything the re-array did afterwards. The well is the join
# key because a well carries one parent across all four conditions.
scan_joined = FALSE
if (file.exists(scan.well.file)) {
    scan = read_csv(scan.well.file, show_col_types = FALSE, progress = FALSE)
    want = c('Wells', 'best', 'verdict', 'd_2nd', 'd_nominal',
             'rank_nominal', 'n_cluster', 'cluster', 'autocor_best',
             'autocor_nominal', 'p_fix_best', 'frac_het_best',
             'best_on_plate', 'n_reads')
    have = intersect(want, names(scan))
    if (!all(c('Wells', 'best', 'verdict') %in% have)) {
        cat('\n!! ', scan.well.file, ' lacks Wells/best/verdict; skipping.\n',
            '   columns: ', paste(names(scan), collapse = ', '), '\n', sep = '')
    } else {
        sc = scan[, have, drop = FALSE]
        names(sc)[names(sc) != 'Wells'] =
            paste0('scan_', names(sc)[names(sc) != 'Wells'])
        rea = left_join(rea, sc, by = 'Wells')
        scan_joined = TRUE
        cat('\n── Joined identity scan:', nrow(scan), 'wells ──\n')
        print(as.data.frame(scan %>% count(verdict, name = 'n_wells') %>%
              arrange(desc(n_wells))), row.names = FALSE)
    }
}
if (!scan_joined) {
    rea = rea %>% mutate(scan_best = NA_character_,
                         scan_verdict = NA_character_)
    cat('\n(no scan at ', scan.well.file, ' — Hamming evidence only)\n',
        sep = '')
}
for (cc in c('scan_d_2nd', 'scan_d_nominal', 'scan_autocor_best',
             'scan_autocor_nominal', 'scan_p_fix_best'))
    if (!cc %in% names(rea)) rea[[cc]] = NA_real_


# ── R4b. Does the scan's call agree with Hamming's? ───────────────────
# Hamming's top-10 is embedded in §R1, so "agrees" can mean the scan's
# winner is anywhere in the near-clonal cluster Hamming reported, not
# just its top hit — for a cluster spanning a few hundred Hamming units
# the two assays cannot be expected to pick the same member, and either
# member phases identically.
ham_top = ham[, c('strain', hcols)]
# Also counts as agreement: the relabel the earlier Fluconazole run
# curated by hand (prior_call / prior_alt in §R1). For CLIB219 that call
# was ADF, which is NOT in Hamming's own top-3 (ADR/BRD/BSC) — so
# checking the top-10 alone would miss a corroboration that exists.
in_ham_cluster = function(strain, code, prior, alt) {
    if (is.na(strain) || is.na(code)) return(NA)
    cands = character(0)
    row = ham_top[ham_top$strain == strain, , drop = FALSE]
    if (nrow(row)) {
        hs = unlist(row[1, hcols], use.names = FALSE)
        hs = hs[!is.na(hs) & nzchar(hs)]
        cands = c(cands, head(hs, 3))
    }
    if (!is.na(prior)) cands = c(cands, prior)
    if (!is.na(alt))   cands = c(cands, trimws(strsplit(alt, '/')[[1]]))
    if (!length(cands)) return(NA)
    code %in% cands
}
rea$scan_matches_ham = as.logical(mapply(in_ham_cluster, rea$strain,
    rea$scan_best, rea$prior_call, rea$prior_alt, SIMPLIFY = TRUE))


# ── R5. Decide the final parent_a ─────────────────────────────────────
# The scan now drives the relabel, because it measures the pools rather
# than the box. But a raw scan REASSIGN is not enough on its own: on the
# 2026-08-24 run it reassigned 39 of 96 wells, and a 41% mislabel rate
# is not credible. The likelihood can be won by something other than
# genuine identity — a lopsided pool, a failed cross, a bad reference
# genotype for the nominal strain — so a call is applied only when a
# second, independent line of evidence agrees.
#
# Tiers, most to least trusted:
#   1 both_assays   scan REASSIGN and its winner is in Hamming's top-3
#                   for that strain. Two unrelated measurements agreeing;
#                   on the run above this recovered ACS, ACR, ADF, ADG
#                   and the ACT/ACV, AMB/AMC, ACK/ADH/ABB clusters.
#   2 scan_linkage  scan REASSIGN with no Hamming to check against, but
#                   the winner yields MORE chromosomal structure than
#                   the nominal (autocor_best > autocor_nominal + margin)
#                   and clears d_2nd. Linkage coherence is the second
#                   independent line: a wrong-but-related parent gives
#                   position-independent noise.
#   3 ham_only      no usable scan for the well, Hamming confident.
#   4 hold_linkage  scan REASSIGN but linkage CONTRADICTS it — keep the
#                   nominal and flag. 15 of 39 landed here, concentrated
#                   on strains whose Hamming reference was junk or absent
#                   (CLIB192 d1=18642, Y10, Y9J, YJM145, DBVPG1339).
#   5 hold_cluster  scan could not separate a near-clonal group.
# -Claude generated
apply_tiers = c(1, 2, 3, 7)  # add 4/5/6 only with a reason
ac_margin   = 0.02           # autocor_best must beat nominal by this
ll_decisive = 100            # and clear this many nats over the runner-up
# Higher bar for tier 7, where there is no nominal to compare against
# and the linkage check therefore cannot arbitrate.
ll_rescue   = 10000

rea = rea %>% mutate(
    scan_reassign = !is.na(scan_verdict) & scan_verdict == 'REASSIGN' &
                    !is.na(scan_best) & scan_best != parent_a_orig,
    ac_supports   = !is.na(scan_autocor_best) & !is.na(scan_autocor_nominal) &
                    scan_autocor_best > scan_autocor_nominal + ac_margin,
    ac_contradicts = !is.na(scan_autocor_best) & !is.na(scan_autocor_nominal) &
                    scan_autocor_best <= scan_autocor_nominal + ac_margin,
    # Hamming can only corroborate where Hamming actually measured
    # something. A strain whose own reference self-matched at d1 = 18642
    # (CLIB192) agrees with nothing meaningfully, and letting it vote
    # would launder noise as consensus.
    ham_usable = ham_class %in% c('agree', 'agree_in_cluster',
                                  'mislabel_confident',
                                  'mislabel_ambiguous'),
    # The scan is the pool measurement, so Hamming does not get to
    # override it — tier 3 fires only where the scan had no opinion.
    scan_silent = is.na(scan_verdict) | scan_verdict == 'insufficient_data',
    tier = case_when(
        scan_reassign & ham_usable & !is.na(scan_matches_ham) &
            scan_matches_ham                                    ~ 1L,
        # Nominal absent from Z entirely (rank_nominal is NA because it
        # was never a candidate). There is no nominal to compare
        # against, so linkage cannot arbitrate — but the alternative is
        # losing the well outright, since it cannot be phased as it
        # stands. Rescue it when the call is decisive on its own terms.
        scan_reassign & is.na(scan_rank_nominal) &
            !is.na(scan_d_2nd) & scan_d_2nd >= ll_rescue        ~ 7L,
        # Genuine numerical anomaly: the nominal WAS a candidate and
        # still scored -Inf. Distinct from the NA case above, which an
        # earlier version conflated — that mislabelled the one well we
        # can actually rescue as "numerically suspect".
        scan_reassign & !is.na(scan_rank_nominal) &
            !is.finite(scan_d_nominal)                          ~ 6L,
        scan_reassign & ac_supports &
            !is.na(scan_d_2nd) & scan_d_2nd >= ll_decisive      ~ 2L,
        scan_silent & ham_class == 'mislabel_confident'          ~ 3L,
        scan_reassign & ac_contradicts                          ~ 4L,
        !is.na(scan_verdict) &
            scan_verdict == 'ambiguous_cluster'                 ~ 5L,
        TRUE                                                    ~ NA_integer_),
    tier_label = c('both_assays', 'scan_linkage', 'ham_only',
                   'hold_linkage', 'hold_cluster', 'hold_numeric',
                   'rescue_unphasable')[tier],
    parent_a_new = ifelse(!is.na(tier) & tier %in% apply_tiers,
        ifelse(tier == 3L, ham_call, scan_best), parent_a_orig),
    parent_a_new = ifelse(is.na(parent_a_new), parent_a_orig, parent_a_new),
    changed = parent_a_new != parent_a_orig,
    review_reason = case_when(
        !is.na(tier) & tier == 4L ~ paste0('scan says ', scan_best,
            ' but linkage contradicts it (ac ',
            round(scan_autocor_nominal, 3), ' -> ',
            round(scan_autocor_best, 3), '); kept ', parent_a_orig),
        !is.na(tier) & tier == 5L ~ paste0('unresolved cluster: ',
            scan_cluster),
        !is.na(tier) & tier == 6L ~ paste0('scan says ', scan_best,
            ' but the nominal scored -Inf despite being a candidate, ',
            'which should not happen; numerically suspect, kept ',
            parent_a_orig),
        !is.na(tier) & tier == 7L ~ paste0('nominal ', parent_a_orig,
            ' is not in Z, so the well was unphasable; rescued to ',
            scan_best, ' on a decisive scan call (no linkage check ',
            'possible)'),
        !is.na(scan_verdict) & scan_verdict == 'insufficient_data' ~
            'scan had too few reads to call',
        ham_class == 'unreliable_reference' ~
            'Hamming reference poorly sequenced; unverifiable there',
        ham_class %in% c('no_hamming_data', 'no_reference_row') &
            !scan_reassign ~ 'no Hamming evidence for this strain',
        TRUE ~ NA_character_))


# ── R6. Check final codes against Z rownames ──────────────────────────
rea = rea %>% mutate(in_Z = NA)
if (check_Z && !is.na(Z.qs.file) && file.exists(Z.qs.file)) {
    Z.obj = qs2::qs_read(Z.qs.file)
    zrows = rownames(Z.obj$Z)
    cat('\nZ has', length(zrows), 'strain rows\n')
    rea = rea %>% mutate(in_Z = parent_a_new %in% zrows &
                                trimws(parent_alpha) %in% zrows)
    missing = rea %>% filter(!in_Z) %>%
        count(strain, parent_a_orig, parent_a_new, name = 'n_samples')
    if (nrow(missing)) {
        cat('\n!! Final codes absent from Z (these samples cannot be phased):\n')
        print(as.data.frame(missing), row.names = FALSE)
        # Is the strain in Z under a different case / prefix?
        for (cd in unique(rea$parent_a_new[!rea$in_Z])) {
            stem = sub('^SACE_', '', cd)
            hits = grep(stem, zrows, ignore.case = TRUE, value = TRUE)
            cat('  ', cd, ' -> Z candidates matching "', stem, '": ',
                if (length(hits)) paste(hits, collapse = ', ') else '(none)',
                '\n', sep = '')
        }
    }
    rm(Z.obj); invisible(gc())
}


# ── R7. Collisions: two different strains sharing one code ────────────
# A relabel can collide with another strain's nominal code. Not fatal
# for the pipeline (parent_a is just a Z row index) but it means two
# physically distinct parents are being phased against one genotype —
# at most one of them can be right.
collide = rea %>%
    distinct(strain, parent_a_new) %>%
    count(parent_a_new, name = 'n_strains') %>%
    filter(n_strains > 1)
if (nrow(collide)) {
    cat('\n!! Code collisions after reassignment:\n')
    print(as.data.frame(rea %>%
        filter(parent_a_new %in% collide$parent_a_new) %>%
        count(parent_a_new, strain, parent_a_orig, ham_class,
              name = 'n_samples') %>%
        arrange(parent_a_new)), row.names = FALSE)
}


# ── R8. Write outputs ─────────────────────────────────────────────────
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

log_cols = c('Sample_ID', 'Plate', 'Wells', 'Drug3Let', 'strain',
             'parent_a_orig', 'parent_a_new', 'changed',
             'tier', 'tier_label',
             'ham_class', 'ham_call', 'ham_d1', 'ham_ratio', 'ham_cluster',
             'prior_call', 'prior_alt',
             'scan_best', 'scan_verdict', 'scan_d_2nd', 'scan_d_nominal',
             'scan_rank_nominal', 'scan_cluster', 'scan_autocor_nominal',
             'scan_autocor_best', 'scan_p_fix_best', 'scan_frac_het_best',
             'scan_best_on_plate', 'scan_n_reads', 'scan_matches_ham',
             'in_Z', 'review_reason')
log = rea %>%
    select(any_of(log_cols)) %>%
    arrange(desc(changed), desc(!is.na(review_reason)), strain, Sample_ID)
write_csv(log, out.log)

out = rea %>%
    mutate(parent_a_preReassign = parent_a_orig,
           parent_a             = parent_a_new,
           reassign_tier        = tier_label,
           reassign_class       = ham_class) %>%
    select(all_of(names(samples)), parent_a_preReassign, reassign_tier,
           reassign_class)

if (drop_unmapped && any(!is.na(rea$in_Z))) {
    n0 = nrow(out)
    out = out[which(rea$in_Z), , drop = FALSE]
    cat('\nDropped', n0 - nrow(out), 'samples with codes absent from Z\n')
}
write_csv(out, out.key)


# ── R9. Report ────────────────────────────────────────────────────────
cat('\n══ Summary ══\n')
cat('samples in keyfile:      ', nrow(samples), '\n')
cat('parent_a reassigned:     ', sum(rea$changed), '\n')
cat('flagged for review:      ', sum(!is.na(rea$review_reason)), '\n')
if (any(!is.na(rea$in_Z))) cat('unphasable (not in Z):   ', sum(!rea$in_Z), '\n')
cat('\nwrote keyfile: ', out.key, '\n', sep = '')
cat('wrote log:     ', out.log,  '\n', sep = '')

cat('\n── Evidence tiers (wells) ──\n')
print(as.data.frame(rea %>% distinct(Wells, tier, tier_label) %>%
    count(tier, tier_label, name = 'n_wells') %>% arrange(tier)),
    row.names = FALSE)
cat('applying tiers:', paste(apply_tiers, collapse = ','), '\n')

if (sum(rea$changed)) {
    cat('\n── Applied reassignments ──\n')
    print(as.data.frame(rea %>% filter(changed) %>%
        count(tier_label, Wells, strain, parent_a_orig, parent_a_new,
              scan_d_2nd, scan_autocor_nominal, scan_autocor_best,
              name = 'n_samples') %>%
        arrange(tier_label, strain)), row.names = FALSE)
    cat('\nRe-run §2 and §3 of xQTL_analysis4.R with key.file pointed at\n',
        out.key, ' — the affected samples need new cross.Z slices and\n',
        'new phased count tables. §4 onward then proceeds normally.\n', sep = '')
}

# Applied and flagged are not exclusive — a tier-2 call on a strain with
# a junk Hamming reference is applied AND worth a note. Separating them
# stops the applied rows appearing under a "held back" heading.
held = rea %>% filter(!is.na(review_reason), !changed)
noted = rea %>% filter(!is.na(review_reason), changed)
cat('\n── Held back (not applied) ──\n')
print(as.data.frame(held %>%
    count(tier_label, strain, parent_a_orig, scan_best, review_reason,
          name = 'n_samples') %>%
    arrange(tier_label, strain)), row.names = FALSE)
if (nrow(noted)) {
    cat('\n── Applied, but carrying a caveat ──\n')
    print(as.data.frame(noted %>%
        count(tier_label, strain, parent_a_orig, parent_a_new,
              review_reason, name = 'n_samples') %>%
        arrange(strain)), row.names = FALSE)
}

# A lopsided pool makes identity unidentifiable regardless of depth: if
# nearly all the reads come from one parent there is little information
# about who the other one was. Worth seeing before trusting any call.
if ('scan_p_fix_best' %in% names(rea) && any(!is.na(rea$scan_p_fix_best))) {
    skew = rea %>% distinct(Wells, strain, parent_a_orig, scan_best,
                            scan_verdict, scan_p_fix_best, scan_d_2nd) %>%
        filter(!is.na(scan_p_fix_best),
               scan_p_fix_best < 0.35 | scan_p_fix_best > 0.65) %>%
        arrange(scan_p_fix_best)
    if (nrow(skew)) {
        cat('\n!! Wells with a lopsided pool (fitted parental ratio far',
            'from 0.5).\n   Identity is weakly determined here whatever',
            'the depth:\n')
        print(as.data.frame(skew), row.names = FALSE)
    }
}
