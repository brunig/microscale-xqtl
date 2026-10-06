# Rationale: what this repository contains, and why

Assembled 2026-10-06. This file records the decisions behind the repository so
a reader, a reviewer or a later lab member can tell what was deliberate.

## 1. One run, not a history

The paper reports run **v7** only, so the repository is a snapshot of the code
that produced v7, not a version history of the analysis. The earlier drivers
(`xQTL_analysis3.R` to `xQTL_analysis6.R`) and the runs made with them (pre-V6,
v6b) are not included, with the one exception of v6 below.

v7 is v6b's analysis with additions, and changes nothing upstream of them. The
acceptance test for v7 was that it reproduce v6b's tables. It does: every
`landscape/` table is byte-identical, and so are `qc/results_key.csv`,
`contrasts/fit_peaks.csv` and `contrasts/stress_peaks.csv`. Two
manuscript-table scripts still read the v6b copies of those files. They were
repointed to v7 when the repository was assembled, and both outputs were
checked unchanged: `NUMBERS.md` differs only in its "generated from" line, and
Table S6 is byte-identical.

**Decided 2026-10-06: v6 is deposited beside v7.** The analysis plan declared
`min_med_depth` a sensitivity analysis, reporting each primary result at both 2
(the run of record) and 5. Run v6 is the `min_med_depth = 5` comparator and
differs from v6b, and so from v7's shared tables, in that one setting only. Its
run directory goes to Zenodo with v7's, carrying its own frozen driver in
`runs/v6/scripts/`. `derived/recorded_from_R.tsv` and the sensitivity column of
`make_workbook.py` read from it. Its code is not repeated here: the driver
that matters for v6 is that frozen copy.

## 2. Hoffman2 is the source of truth for the pipeline

The pipeline ran on UCLA Hoffman2. The working copies on the analysis laptop had
drifted, so every pipeline file was checked by MD5 against an inventory taken on
Hoffman2 (`code_md5.txt`, 2026-10-06). A file was copied from the laptop only
when its checksum matched. Otherwise it is fetched from Hoffman2 by
`tools/hoffman_collect.sh`. `MANIFEST.tsv` records which is which.

The driver is the case that matters. Hoffman2's
`runs/v7/scripts/xQTL_analysis7.R` (211,583 bytes, MD5 `44ffae4…`) matches
neither laptop copy: one is an earlier snapshot (208,942 bytes) and the other
carries later edits (212,278 bytes). The repository takes the Hoffman2 run copy.
That is the file the driver writes into its own run directory as it runs.

**Which driver produced which stage (settled 2026-10-06 from Hoffman2
timestamps, `xqtl_v7_collect/v7_timestamps.txt`).** v7 ran in two sittings:

| stage | ran | driver |
|---|---|---|
| §1-15 allele frequencies, contrasts, peaks | from v6b's cache (warm start); identical by the acceptance test | v6b driver, same §1-15 code |
| variants, CNV tables | 2026-09-21 12:08 | v7 driver, full `stages` |
| §16 landscape, §18 GO | 2026-09-21 14:20-14:48 | v7 driver, `stages = c('landscape','go')` |
| §17 candidates | 2026-09-22 13:54 | the deposited run copy, `stages = c('candidates')` |

The run copy is rewritten each time the driver runs, so it is the driver of the
LAST sitting. It differs from the 09-21 laptop snapshot in four places only:
the `stages` line, the §16 L/R column rename, the §17 NA-score fix and the §18
GO focus word boundaries. The deposited §16 and §18 outputs carry both of
their fixes (`landscape.log` has no `!!` lines and the contrast enrichment
rows are populated; the GO focus table holds no "replication checkpoint" term),
so the code that produced every stage is the deposited driver's code. To
reproduce a stage, set `stages` (line 669) accordingly.

One later edit exists only on the laptop and was never run: a tie rule for the
§17 allele-series self-check (`rank_mid`, average rank across tied genes). It
is not in the repository. If the manuscript adopts it, §17 is re-run on
Hoffman2 and the run copy replaces `pipeline/xQTL_analysis7.R`.

The driver's inputs from outside the run, by its `hoffman2` path block
(lines 560-579): `count_files/`, `files/flatkey_A1.csv`, `files/Z_final.qs`,
`qc/` (allowlist and unresolved wells), the 1,011-isolate `variant_lookup.qs`
and CNVnator matrix, the GFF, the fungAMR supplement, and the SGD GAF and OBO.
The lab-made ones are deposited at Zenodo; the published ones are cited with
their checksums (section 6).

## 3. Included, and why

| directory | why it is here |
|---|---|
| `upstream/` | Reads to allele counts. Without it the counts the driver reads cannot be regenerated. `parents_A1/` aligns the individually sequenced parents that `parents/` uses as its reference panel. |
| `qc/` | `qc_pipeline.R` writes the two files the driver gates on (`sample_allowlist.csv`, `unresolved_wells_A1.txt`). The driver reads them and never writes them, so they come from outside it. |
| `parents/` | 216 of 384 wells hold a different wild parent from the nominal plate map. The reassignment method is spread over these scripts and was never written up as one procedure, so the scripts are the method of record. `flatkey_A1.csv` is their product and the driver's key file. |
| `pipeline/` | The driver. It contains the landscape, candidate and GO analyses inline (V6.1 onward), so no separate child scripts are needed. |
| `derived/` | Every number and table in the manuscript is generated by one of these, not typed. |
| `figures/` | Figures 2 to 4 are generated from the run. Figure 1 is drawn by hand from the SVG blueprint `make_fig1.py` writes. |
| `browser/`, `docs/` | See section 5. |

## 4. Left out, and why

| file(s) | reason |
|---|---|
| `xQTL_analysis3.R` to `6.R`, `v6b/` | Superseded; see section 1. `v6/` goes to Zenodo as the sensitivity comparator, not into git. |
| `qtl_landscape.R`, `qtl_candidates.R`, `qtl_go.R` (standalone) | No longer used. V6.1 moved the first two into the driver, and the standalone copies then drifted from it; `qtl_go.R` became the driver's §18. Including them would invite running code that did not produce the results. |
| `patch_cache.R` | Repointed a V4 cache at a new key without a rebuild. v7's cache was computed from scratch. |
| `detect_skew.R` | The skew screen it prototyped is implemented in the driver (`skew_cut`). |
| `diag_well_A01.R` | One-off diagnostic of a single well's scoring. |
| `qc_bam_parse_plot.R` | Replaced by `qc_pipeline.R`, by that script's own header. |
| `plot_overrides.R` | Styling for the pre-V6 figures. |
| `finemap/` | Exploratory chromosome VII and control analyses; not cited in the manuscript. |
| `misc/interval_test/` | Run against the pre-V6 contrasts, not v7. The statement it supports (the interval rule does not cause the overlap saturation) should be re-run on v7 if it stays in the text. |
| `make_tables.py`, `add_supp_sheets.py` | Superseded by `make_workbook.py`. |
| `.ipynb_checkpoints/` | Jupyter autosaves. |

## 5. The browser lives here

The browser is the paper's companion and is built from the same run, so it ships
with the code that made that run. The site is served by GitHub Pages from
`/docs` at the repository root, which needs no build step and no Actions
workflow. The code that builds the site is in the repository too:

- `browser/scripts/export.R` builds `docs/data.json` from `runs/v7`: traces in
  2 kb bins, peaks, QC status, published loci and genes. Its significance lines
  use the same rule as `figures/make_figures.R`.
- `browser/scripts/export_variants.R` builds `docs/v/`, the per-chromosome
  variant, CNV and candidate-gene tables, from the same run.
- `docs/index.html` is the site's own code, a single page with no dependencies
  beyond Google Fonts.

The gene track uses the SGD R64-5-1 annotation and is for display only. The
pipeline took its gene models from a different GFF on Hoffman2 (5,160 against
5,211 genes). Gene-level numbers therefore come from the pipeline, never from
the browser.

**Nothing goes public before the preprint is posted** (decided 2026-10-06): the
GitHub repository stays private, GitHub Pages stays off, and the Zenodo records
are held as drafts or restricted until then. Before release, check with the
senior authors: the site exposes every pool's trace, including flagged pools,
the parent reassignments and the variant tables.

## 6. Data are not in git

Reads, alignments, caches and run directories are excluded by `.gitignore`.
Raw reads go to NCBI SRA. The v7 and v6 run directories and the inputs the
driver reads go to Zenodo, built and uploaded on Hoffman2 by
`tools/zenodo_package.sh` and `tools/zenodo_upload.sh`. Excluded from Zenodo:
the reads (SRA), the BAMs (regenerated from the reads by step 1) and the `.qs`
caches (regenerated from the count tables by the driver). Third-party
reference files (CNVnator matrix, fungAMR supplement, GAF, OBO) are cited by
source and MD5 in the record's `INPUTS.tsv` rather than redistributed; the GFF
is deposited because the one used differs from SGD's current release. Each
deposit is cited in the README once it has an accession.

## 7. Known rough edges

- Steps 6 to 8 use absolute paths (`ROOT`, `RUN`) from the analysis laptop.
  Left as run rather than rewritten untested; the README says what to change.
- `make_figures.R` resolves the cache by filename pattern and, in `runs/v7/cache`,
  finds both `v6_*` (hardlinked from v6b at warm start) and `v7_*` files. It
  picks the larger file. Per-site values are identical between the two because
  v7 changed nothing upstream of §16, but a clean deposit should carry only
  `v7_*`.
- The `LODfromZ` conversion in the driver uses twice the two-sided P
  (`log(4) + pnorm(-|z|)`). The Methods sentence and the code must agree before
  publication; one of them changes.
