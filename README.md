# Microscale X-QTL: code for the v7 analysis

Code behind *Microscale bulk segregant analysis and genetic mapping* (Bruni et
al., G3): a 96-well bulk-segregant X-QTL protocol in *Saccharomyces cerevisiae*,
with 36 wild isolates crossed to the common parent BBQ and selected in glucose,
acetic acid, caffeine and sodium chloride. That is 384 wells in all.

This repository holds the code for **run v7**, the run every number, table and
figure in the paper comes from. It also holds the **Microscale X-QTL Browser**,
the paper's interactive companion site, and the code that builds it.

**Browser:** <!-- GitHub Pages URL, once enabled -->
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
         qc_pipeline.R        pre-flight QC; writes sample_allowlist.csv and
                              unresolved_wells_A1.txt, which step 4 reads
    3  parents/             which wild parent is in each well (216 of 384
                            reassigned from the nominal plate map)
         hamming_export.py, hamming_integrate.R
                              genotype distance of each sample to every
                              candidate strain
         scan_parent_identity.R   per-well likelihood over all strains
         reassign_parents.R       applies the calls; holds back the ambiguous
         recheck_held_wells.R     resolves the held wells
         round3_apply.R           builds the final key from one rule
         verify_parents.R         independent check of the final key
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

## Requirements

R 4.5 with `data.table`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `patchwork`,
`ragg`, `qs2`, `jsonlite`, `IRanges`, `GenomicRanges`. Python 3 with `Pillow`
and `openpyxl`. Steps 1 to 2 also need Nextflow, samtools and GATK
ASEReadCounter. The reference is sacCer3 with the 1,011-isolate SNP panel
(Peter et al. 2018).

## Paths

The steps 6 to 8 scripts were written against a local checkout and contain
absolute paths, e.g. `ROOT <- "/Users/…/AFxQTL"`. Point `ROOT`/`RUN` at your
copy of the Zenodo run directory before running them.

## License

Code: MIT (`LICENSE`). Data at Zenodo: CC-BY 4.0.
