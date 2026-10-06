#!/usr/bin/env python3
"""Regenerate draft/derived/NUMBERS.md from the authoritative per-well table.

    python3 draft/derived/make_numbers.py

Every number this script prints is computed here, from results_key.csv, at run
time.  Numbers that come out of the R pipeline instead are NOT recomputed --
they are read verbatim from recorded_from_R.tsv, because reimplementing the
mask and peak-calling logic in Python would let the two drift apart silently.
"""

import csv, collections, statistics, subprocess, datetime
from pathlib import Path

ROOT     = Path("/Users/kruglyakbws/Desktop/AFxQTL")
KEY      = ROOT / "hoffmanDrop/notAF/v7/qc/results_key.csv"
RECORDED = Path(__file__).parent / "recorded_from_R.tsv"
OUT      = Path(__file__).parent / "NUMBERS.md"

# Glucose is the control and is listed first everywhere; then alphabetical.
CONDS = ["Glu", "Ace", "Caf", "NaCl"]
LONG  = {"Glu": "glucose", "Ace": "acetic acid", "Caf": "caffeine", "NaCl": "sodium chloride"}

rows = list(csv.DictReader(open(KEY)))
main = [r for r in rows if r["databin"] == "main"]
T    = lambda r, c: r[c].strip().upper() == "TRUE"
med  = statistics.median


def q7(xs, p):
    """R's default quantile (type 7), so this file cannot disagree with the R pipeline."""
    xs = sorted(xs)
    h  = (len(xs) - 1) * p
    lo = int(h)
    return xs[lo] + (h - lo) * (xs[min(lo + 1, len(xs) - 1)] - xs[lo])


def quartiles(xs):
    return q7(xs, 0.25), med(xs), q7(xs, 0.75)


L = []
w = L.append

w("<!-- DERIVED -->")
w("<!-- GENERATED FILE - do not hand-edit.  Regenerate with:")
w("       python3 draft/derived/make_numbers.py")
w("     Hand-maintained inputs live in draft/derived/recorded_from_R.tsv -->")
w("")
w("# NUMBERS")
w("")
w(f"Generated {datetime.date.today().isoformat()} from `{KEY.relative_to(ROOT)}`")
w(f"({len(rows)} wells x {len(rows[0])} columns).")
w("")
w("Every value below is recomputed at generation time.  If a number in the")
w("manuscript disagrees with this file, this file is right.")
w("")

# ---------------------------------------------------------------- design ----
nominal  = sorted(set(r["parent_nominal"] for r in rows))
resolved = sorted(set(r["parent_a"] for r in rows))
changed  = sum(r["parent_a"] != r["parent_a_preReassign"] for r in rows)
crosses  = sorted(set(r["cross"] for r in rows))
percell  = collections.Counter((r["cross"], r["condition"]) for r in rows)
cells    = [percell[(x, c)] for x in crosses for c in CONDS]

w("## Design")
w("")
w("| quantity | value |")
w("|---|---|")
w(f"| wild isolates on the plate map (`parent_nominal`) | **{len(nominal)}** |")
w(f"| wild isolates after sequence-based assignment (`parent_a`) | **{len(resolved)}** |")
w(f"| wells whose parent changed on reassignment | **{changed}** of {len(rows)} "
  f"({100*changed/len(rows):.0f}%) |")
w(f"| distinct crosses | {len(crosses)} (every cross is `<wild> x BBQ`) |")
w(f"| conditions | {len(CONDS)}, one 96-well plate each |")
w(f"| total wells | {len(rows)} |")
w(f"| cross-by-condition cells | {len(cells)}, of which {cells.count(0)} empty |")
w(f"| wells per cross per condition | median {med(cells):.0f}, range {min(cells)} to {max(cells)} |")
w("")
w(f"Dropped on reassignment: {', '.join(sorted(set(nominal) - set(resolved)))}.")
w("")
w(f"Never on the plate map, recovered by assignment: "
  f"{', '.join(sorted(set(resolved) - set(nominal)))}.")
w("")
w("Reassignment class of the wells that changed:")
w("")
w("| `reassign_class` | wells changed |")
w("|---|---|")
ch = [r for r in rows if r["parent_a"] != r["parent_a_preReassign"]]
for k, v in collections.Counter(r["reassign_class"] for r in ch).most_common():
    w(f"| `{k}` | {v} |")
w("")

# ---------------------------------------------------------------- funnel ----
bins  = collections.Counter(r["databin"] for r in rows)
flags = [c for c in rows[0] if c.startswith("flag_")]

w("## Well funnel")
w("")
w("| `databin` | wells |")
w("|---|---|")
for b in ["main", "flagged", "no_data"]:
    w(f"| {b} | **{bins[b]}** |")
w(f"| total | {len(rows)} |")
w("")
w("`main` is the analysed set.  `flagged` wells returned an allele-frequency")
w("estimate but carry a quality flag and are held out.  Flags per bin:")
w("")
w("| flag | no_data | flagged | main |")
w("|---|---|---|---|")
for f in flags:
    counts = {b: sum(T(r, f) for r in rows if r["databin"] == b) for b in ("no_data", "flagged", "main")}
    if sum(counts.values()):
        w(f"| `{f}` | {counts['no_data']} | {counts['flagged']} | {counts['main']} |")
w("")
w("Flags are not mutually exclusive, so columns sum to more than the bin size.")
w("")

# ------------------------------------------------------------- retention ----
w("## Retention by condition")
w("")
w("| condition | wells | main | flagged | no_data | retained |")
w("|---|---|---|---|---|---|")
for c in CONDS:
    sub = [r for r in rows if r["condition"] == c]
    b   = collections.Counter(r["databin"] for r in sub)
    w(f"| {LONG[c]} | {len(sub)} | **{b['main']}** | {b['flagged']} | {b['no_data']} | "
      f"{100*b['main']/len(sub):.0f}% |")
w("")

# ----------------------------------------------------------------- depth ----
dep   = [float(r["mean_depth_seg"]) for r in main if r["mean_depth_seg"] not in ("", "NA")]
sites = [float(r["n_sites"]) for r in main if r["n_sites"] not in ("", "NA")]
q1, q2, q3 = quartiles(dep)

w("## Depth and sites, `main` pools only")
w("")
w("| quantity | value |")
w("|---|---|")
w(f"| pools with a depth value | {len(dep)} of {len(main)} |")
w(f"| mean depth | median **{q2:.1f}x** (IQR {q1:.1f}-{q3:.1f}) |")
w(f"| informative sites | median **{med(sites):,.0f}** |")
w(f"| pools reaching the 30x target | **{sum(d >= 30 for d in dep)}** |")
w(f"| pools below 10x | {sum(d < 10 for d in dep)} ({100*sum(d < 10 for d in dep)/len(dep):.0f}%) |")
w("")
w("| condition | main pools | median depth | below 10x |")
w("|---|---|---|---|")
for c in CONDS:
    d = [float(r["mean_depth_seg"]) for r in main
         if r["condition"] == c and r["mean_depth_seg"] not in ("", "NA")]
    if d:
        w(f"| {LONG[c]} | {len(d)} | {med(d):.1f}x | {100*sum(x < 10 for x in d)/len(d):.0f}% |")
w("")

# ---------------------------------------------------------------- pairing ---
w("## Stress pairing")
w("")
w("| `stress_status` | pools |")
w("|---|---|")
for k, v in collections.Counter(r["stress_status"] for r in main).most_common():
    w(f"| {k or '(empty)'} | {v} |")
w("")
w(f"Crosses represented in `main`: **{len(set(r['cross'] for r in main))}**.")
w("")
w("| condition | main pools | crosses |")
w("|---|---|---|")
for c in CONDS:
    sub = [r for r in main if r["condition"] == c]
    w(f"| {LONG[c]} | {len(sub)} | {len(set(r['cross'] for r in sub))} |")
w("")

# ------------------------------------------------- recorded from R pipeline --
w("## Recorded from the R pipeline")
w("")
w("These are **not** recomputed here.  Masking and peak calling live in")
w("`hoffmanDrop/notAF/qtl_landscape.R` and `qtl_candidates.R`; a second Python")
w("implementation could drift from them without anyone noticing.  Values are")
w("read verbatim from `recorded_from_R.tsv` -- edit that file when the R")
w("pipeline is rerun, then regenerate this one.")
w("")
if RECORDED.exists():
    rec = list(csv.DictReader(open(RECORDED), delimiter="\t"))
    # The depth threshold is a declared sensitivity analysis: the run of record
    # is min_med_depth = 2 (v6b) and v6 is kept as the min_med_depth = 5
    # comparator.  Rows that are reported at both cuts carry a value_mmd5.
    sens = any((r.get("value_mmd5") or "").strip() for r in rec)
    if sens:
        w("The run of record is `min_med_depth = 2`.  The depth threshold is a")
        w("declared sensitivity analysis, so rows reported at both cuts carry the")
        w("`min_med_depth = 5` value from run v6 beside them; a blank means the")
        w("quantity is reported at the run of record only.")
        w("")
        w("| quantity | value (mmd = 2, v6b) | at mmd = 5 (v6) | source |")
        w("|---|---|---|---|")
        for r in rec:
            alt = (r.get("value_mmd5") or "").strip()
            alt = f"**{alt}**" if alt else "--"
            w(f"| {r['quantity']} | **{r['value']}** | {alt} | `{r['source']}` |")
    else:
        w("| quantity | value | source |")
        w("|---|---|---|")
        for r in rec:
            w(f"| {r['quantity']} | **{r['value']}** | `{r['source']}` |")
else:
    w("*(recorded_from_R.tsv not found)*")
w("")

# ------------------------------------------------------------------- open ---
w("## Open placeholders")
w("")
w("| token | meaning | status |")
w("|---|---|---|")
w("| `[N]` | segregants per pool | not yet measured |")
w("| `[$X]` | reagent cost per experiment | not yet costed |")
w("| MATa / MATalpha | which haploid class was selected | `parent_alpha` is BBQ in all "
  f"{len(rows)} wells and `selected_haploid` is `x`; confirm at the bench before it goes in a figure |")
w("")

OUT.write_text("\n".join(L) + "\n")
print(f"wrote {OUT}  ({OUT.stat().st_size} bytes, {len(L)} lines)")
