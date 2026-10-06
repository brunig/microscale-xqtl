Microscale X-QTL: processed data for runs v7 and v6
====================================================

Processed data behind "Microscale bulk segregant analysis and genetic mapping"
(Bruni et al.): a 96-well bulk-segregant X-QTL protocol in Saccharomyces
cerevisiae. 36 wild isolates crossed to the common parent BBQ, selected in
glucose, acetic acid, caffeine and sodium chloride; 384 wells.

Code: GitHub repository and its Zenodo code record (see Related identifiers).
Raw reads: NCBI SRA (see Related identifiers).

FILES
-----
runs_v7.tar.gz      Run v7, the run of record. Every number, table and figure
                    in the paper comes from it. Output of xQTL_analysis7.R:
                      qc/          results_key.csv (per-well status; databin
                                   "main" = pools analysed), reproducibility
                      contrasts/   fit_peaks.csv (fitness contrast: allele
                                   frequency against 0.5), stress_peaks.csv
                                   (stress pool against its matched glucose
                                   pool), variant and CNV tables per QTL
                      landscape/   published-locus recovery, enrichment and
                                   cross-level replication tests
                      candidates/  candidate genes per QTL
                      go/          GO enrichment, positional permutation null
                      scripts/     the driver exactly as it last ran
                      *.log        logs of each stage
runs_v6.tar.gz      Run v6, identical except min_med_depth = 5 (v7: 2). The
                    declared sensitivity comparator. Carries its own driver
                    copy in scripts/.
count_files.tar.gz  Per-well allele counts at the 1,011-isolate SNP panel
                    sites: the driver's primary input, one file per well.
driver_inputs.tar.gz
                    flatkey_A1.csv            well -> cross, condition, parents
                                              after reassignment
                    Z_final.qs, variant_lookup.qs
                                              precomputed objects the driver
                                              reads (R, package qs2)
                    saccharomyces_cerevisiae.gff.gz
                                              the GFF the gene models came
                                              from; it differs from SGD's
                                              current release (5,160 genes)
qc.tar.gz           Pre-flight QC output of qc_pipeline.R, including
                    sample_allowlist.csv and unresolved_wells_A1.txt, which
                    the driver gates on.
INPUTS.tsv          Every driver input with MD5 and size. Rows marked "cite
                    source" are published third-party files not redistributed
                    here: the 1,011 Yeast Genomes CNVnator matrix (Peter et
                    al. 2018), the fungAMR supplement, and the SGD GAF
                    R64-5-1 and go-basic.obo.
MD5SUMS.txt         Checksums of the files above.

NOT HERE
--------
Reads (SRA). BAM files (regenerated from the reads by upstream/ in the code
repository). The .qs caches in runs/*/cache (regenerated from count_files by
the driver).

KEY SETTINGS, v7
----------------
min_med_depth 2, skew_cut 0.05, keep_chr7 TRUE, LOD_threshold 5,
LOD_drop_frac 0.1, interval floor 20 kb / cap 200 kb, diverged_frac 0.01,
max_excess 0.12, afd_sample_size 10,000.

LICENSE
-------
CC-BY 4.0.
