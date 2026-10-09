# Microscale X-QTL: code for the v7 analysis

Code behind *Microscale bulk segregant analysis and genetic mapping* (Bruni et
al., G3): a 96-well bulk-segregant X-QTL protocol in *Saccharomyces cerevisiae*,
with 36 wild isolates crossed to the common parent BBQ and selected in glucose,
acetic acid, caffeine and sodium chloride. That is 384 wells in all.

This repository holds the code for **run v7**, the run every number, table and
figure in the paper comes from. It also holds the **Microscale X-QTL Browser**,
the paper's interactive companion site, and the code that builds it.

**Browser:** https://brunig.github.io/microscale-xqtl/
**Data:** reads at SRA <!-- BioProject PRJNA… -->; processed data at Zenodo <!-- DOI -->
**Code DOI:** <!-- Zenodo DOI minted from the GitHub release -->

`RATIONALE.md` explains what is here, what was left out and why.
`MANIFEST.tsv` lists every file with its source and checksum.

## Layout and run order

Steps run top to bottom. Steps 1 to 5 ran on UCLA Hoffman2. Steps 6 to 8 run
on a workstation from a copy of the run directory.

    1  upstream/            reads -> alignments -> per-site allele counts
         run_pipeline.sh      launcher (local or Hoffman2 SGE)
         main.nf, *.config    Nextflow pipeline
         parents_A1/          alignment of the individually sequenced haploid
                              parents, used as the reference panel for step 3
    2  qc/
         qc_bam_generate.sh   samtools stats for every BAM
         qc_pipeline.R        pre-flight QC; writes sample_allowlist.csv,
                              which step 4 reads
    3  parents/             which wild parent is in each well (216 of 384
                            reassigned from the nominal plate map)
         hamming_export.py, hamming_integrate.R
                              genotype distance of each sample to every
                              candidate strain
         scan_parent_identity.R   per-well likelihood over all strains
         reassign_parents.R       applies the calls; holds back the ambiguous
         recheck_held_wells.R     resolves the held wells
         round3_apply.R           round 3: one rule on overdispersion and fit
         build_parent_Z_maj.R     genotypes of the sequenced A1 haploids
         test_parents_A1_maj.R    tests disputed wells' reads against the A1
                                  haploids; rejects candidates, never confirms
         build_final_key.R        applies that test over round 3; writes
                                  flatkey_A1.csv and unresolved_wells_A1.txt,
                                  which step 4 reads
         scan_all_candidates.R, scan_919.R, reconstruct_parents.R,
         reconstruct_cv.R         evidence behind the 15 wells called by hand
         test_row_source.R        per strain, 1011-panel or A1 genotype row
         Z_final_provenance.csv   its choices; Z_final.qs was assembled from
                                  them by hand, no script writes it
         verify_parents.R, verify_final_key.R
                                  independent checks of the final key
         flatkey_A1.csv           THE KEY: sample -> cross, well, condition
    4  pipeline/
         xQTL_analysis7.R     the analysis driver. Allele frequencies,
                              fitness and stress contrasts, peak calling,
                              published-locus tests (§16), candidate genes
                              (§17), GO enrichment (§18). Writes runs/v7/.
    5  (output)             runs/v7/ -> deposited at Zenodo, not in git
                            runs/v6/ (min_med_depth = 5 sensitivity run) likewise
    6  derived/             numbers and tables quoted in the manuscript
         make_numbers.py      NUMBERS.md, from qc/results_key.csv and
                              recorded_from_R.tsv
         estimate_pool_size.R Table S7 (effective segregants per pool)
         make_go_candidates.R Table S6
         make_workbook.py     TABLES.xlsx, every table, one per tab
         make_table_legends.py  lifts the captions out of make_workbook.py
    7  figures/
         make_figures.R       Figures 2 to 4, from runs/v7
         make_fig1.py         Figure 1 blueprint (SVG)
         make_figure_pdf.py   collates the figures into one PDF
    8  browser/scripts/     builds the site in docs/ from runs/v7
       docs/                the site itself, served by GitHub Pages

The order inside `parents/` was reconstructed from the scripts' own headers.

## The run of record: v7

`xQTL_analysis7.R` with `run_tag = "v7"`. Settings, as recorded in the run's
manifest:

| setting | value |
|---|---|
| `min_med_depth` | 2 |
| `skew_cut` | 0.05 |
| `keep_chr7` | TRUE |
| `LOD_threshold` | 5 |
| `LOD_drop_frac` | 0.1 (interval runs until LOD falls by 10% of the summit's) |
| interval floor / cap | 20 kb / 200 kb |
| `diverged_frac` | 0.01 |
| `max_excess` | 0.12 |
| `afd_sample_size` | 10,000 |
| GO annotation | SGD GAF R64-5-1 + go-basic.obo, 1,000 positional permutations |

Of 384 wells, 181 pools pass every filter (`databin == "main"`). They carry
1,829 fitness peaks, of which 1,472 remain after masking MAT and HO, plus 363
stress-contrast peaks.

## Software versions

Steps 1 to 5 ran on Hoffman2 in the conda environment `genomicsV2`. Its
revision history shows one change after it was built (2026-04-15): on
2026-08-21, two days after the alignment run (2026-08-19) and a month before
the analysis runs (2026-09-21/22), samtools 1.24 was added (for QC) and
htslib/bcftools went from 1.23.1 to 1.24. Neither is called by the alignment
pipeline, so the versions below are the ones every step ran under.

| step | software | version |
|---|---|---|
| 1 | Nextflow | 25.10.4 |
| 1 | bwa (mem) | 0.7.19-r1273 |
| 1 | sambamba | 1.0.1 |
| 1 | GATK (ASEReadCounter) | 4.6.2.0, OpenJDK 23.0.2 |
| 2 | samtools | 1.24 |
| 2 to 4 | R | 4.3.3 |
| 4 | xQTLStats | 0.2, github.com/joshsbloom/xQTLStats @ 080a156dcc |
| 4 | data.table, dplyr, tidyr, readr, purrr | 1.17.8, 1.1.4, 1.3.1, 2.1.5, 1.1.0 |
| 4 | ggplot2, ggrepel, patchwork | 3.5.2, 0.9.6, 1.3.2 |
| 4 | qs2, pbmcapply | 0.1.5, 1.5.1 |
| 4 | IRanges, GenomicRanges | 2.36.0, 1.54.1 (Bioconductor 3.18) |
| 4 | ComplexHeatmap, circlize | 2.18.0, 0.4.18 |

Steps 6 to 8 ran on a Mac (conda environment `genomics`): R 4.5.3 with
data.table 1.18.2.1, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, purrr 1.2.2,
ggplot2 4.0.2, patchwork 1.3.2, ragg 1.5.2, qs2 0.1.7, jsonlite 2.0.0, IRanges
2.44.0, GenomicRanges 1.62.1; Python 3.10.20 with openpyxl 3.1.5 and
Pillow 12.3.0.

The reference is sacCer3 with the 1,011-isolate SNP panel (Peter et al. 2018).

## Paths

The steps 6 to 8 scripts were written against a local checkout and contain
absolute paths, e.g. `ROOT <- "/Users/…/AFxQTL"`. Point `ROOT`/`RUN` at your
copy of the Zenodo run directory before running them.

## License

Code: MIT (`LICENSE`). Data at Zenodo: CC-BY 4.0.
