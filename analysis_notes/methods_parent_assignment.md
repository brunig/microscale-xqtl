<!-- HANDBACK: proposed Methods text for the user to place. Not the manuscript. Facts verified 2026-10-07 against parents/flatkey_A1.csv, the A1ref scripts and qc/Z_final_provenance.csv. Describes the key the v7 preprint ran on. -->

## Methods — Assigning the wild parent of each well (proposed text)

**Parent assignment.** Each well was seeded from a cross between the common parent BBQ and one wild isolate, named by the plate map. Because isolates can be mislabelled or exchanged during handling, we identified the wild parent of every well from the pools' own sequence. For each pool, we scored every candidate isolate in the 1011 Yeast Genomes panel (919 strains after filtering) on a fixed set of segregating sites. A site was kept if at least two strains were homozygous for the alternative allele and no more than 10% were heterozygous. The score was a binomial likelihood of the observed reads with BBQ fixed as one parent and the mixing proportion profiled at sites where the two parents differ. Log-likelihoods were summed over the four pools of a well, which share their parents. A well was reassigned when another strain beat the plate-map strain by at least 100 log-likelihood units.

Candidate calls were then confirmed on two statistics that do not depend on the ranking. The first was mean |allele frequency − 0.5| (*d*_half), which is near zero across the genome when the parent is right and large when the assumed parent genotype is wrong. The second was the overdispersion of allele frequency around a 501-site running mean, relative to binomial expectation (*od*; sites with depth ≥ 3). A candidate was accepted at *od* ≤ 1.45 when it improved on the plate-map strain. We also sequenced haploid stocks of the parents ourselves (61 libraries, 54 strain labels). These reads tested the two remaining candidates of each disputed well, a candidate being rejected at *d*_half > 0.40 (or > 0.20 with more than 20,000 segregating sites). For nine parents whose own sequence fitted the pools better than their 1011 entry, we used our sequence as that parent's genotype throughout (AAM, ABF, ABG, ACN, ADE, ADI, ADN, ADR, ADS). Each parent's genotype therefore comes from a single source.

Of 384 pools, 216 (54 of 96 wells) were assigned a parent other than the one on the plate map. The 40 plate-map names resolved to 35 strains. Twelve plate-map names were never supported, and eight strains not on the plate map were recovered (AAL, ABB, ACH, ACS, ADA, ADF, ADG, AMB). Six wells (24 pools) could not be resolved to any candidate and are flagged and excluded from analysis.

<!-- For the v8 submission version, add: -->
<!-- "A final check counted, for each well, sites with at least 20 reads at which one parental allele never appeared in any of its pools. It reversed two reassignments (wells A11 and G10, where the plate-map strain fitted better) and left a third well unresolved (A10)." -->

---

### Notes for you, not for the text

- **The "35 strains" and "six wells" are v7's key.** v8b changes the parents of three wells: A11 back to AAM, G10 back to ADD, A10 unresolved in place of G10. It also drops the AAL cross. The commented sentence above covers this for the submission.
- **Undocumented steps:**
  - 15 wells got their final call from targeted rescans (a 47-candidate scan, a 919-strain rescan, and parent reconstruction from pool reads). That's a hand-entered `decision` column with scripts beside it.
  - The `override_source` edits after round 2 have no generating script.
  - Nothing in the repo writes `Z_final.qs`; only its provenance table exists. I haven't seen how the 1011 panel was filtered from 1,011 strains to 919.

  The text above is true of the outcome. A reviewer asking for code would find gaps here.
- **The 100-unit threshold and the round structure** are from `scan_parent_identity.R` and `round3_apply.R`. If you'd rather not quote numbers this fine-grained in a preprint, cut the second paragraph down to its first sentence.
