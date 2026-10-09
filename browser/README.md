# Microscale X-QTL Browser

An interactive companion to *Microscale bulk segregant analysis and genetic mapping*
(Bruni et al.). It shows every pool from the 96-well X-QTL experiment: 36 yeast
isolates crossed to the common parent BBQ, selected in glucose, acetate, caffeine
and sodium chloride, 384 wells in all.

**Open it:** https://brunig.github.io/microscale-xqtl/

## What it shows

- **Genome view.** One pool's smoothed allele frequency and its reads, over a
  heatmap of how many crosses carry a QTL in each 50 kb bin per condition.
- **Every well.** The four plates, with each well marked as analysed (main),
  held out with the reason (flagged), or no data. Nothing is hidden: flagged
  pools can be drawn and included in every view.
- **Region view.** All pools' traces over a region, as the fitness contrast
  (allele frequency against 0.5) or the stress contrast (a stress pool minus
  its glucose pool from the same well), with peak intervals, published loci
  (Treusch et al. 2015; Ehrenreich et al. 2010) and genes.
- **Under a peak.** Click a peak to list the variants that differ between the
  two parents inside its interval, ranked by consequence and PROVEAN score, with
  the allele the favoured parent carries and how many detecting crosses share
  it. A second tab gives the pipeline's within-locus candidate-gene ranking.
- **Search.** Genes, coordinates (`chrIV:509600`, `IV:400k-600k`), crosses
  (`AAM`), pools (`Sample302`) and published loci (`Caff.chr10.1`).

The browser shows what the data contain. No gene shown here is claimed as
causal.

## Layout

    docs/               (repo root) the site, served by GitHub Pages from /docs
      index.html        the page; no build step, no dependencies beyond Google Fonts
      data.json         traces (2 kb bins), peaks, QC status, published loci, genes
      v/                variants, CNVs and candidate genes, one file per chromosome,
                        loaded only when a peak on that chromosome is opened
    browser/scripts/
      export.R          builds docs/data.json from a pipeline run
      export_variants.R builds docs/v/ from the same run

## Rebuilding from a pipeline run

Both scripts read a run umbrella written by the analysis driver
(`xQTL_analysis7.R`): `cache/`, `contrasts/`, `qc/`, `landscape/` and
`candidates/`. They need R with `qs2`, `jsonlite` and `data.table`.

    # from the repository root
    Rscript browser/scripts/export.R          <run_dir> docs <annotation.gff.gz>
    Rscript browser/scripts/export_variants.R <run_dir> docs

The gene track uses the SGD S288C annotation, R64-5-1. Variant and
candidate-gene tables come from the pipeline itself.

## Data

Run of record v7: median-depth threshold 2, skew cut 0.05, peaks at LOD 5 with
an interval extending until LOD falls by 10% of the summit's, 20 kb floor,
200 kb cap. Smoothed traces are averaged into 2 kb bins for display.
