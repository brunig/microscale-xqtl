#!/usr/bin/env python3
"""Emit draft/table_legends.md from the captions already in TABLES.xlsx, so the
two cannot drift.  The captions are authored in make_workbook.py; this only
lifts them out.

    python3 draft/derived/make_table_legends.py
"""
import re
from pathlib import Path
import openpyxl

HERE = Path(__file__).resolve().parent
WB   = HERE / "TABLES.xlsx"
OUT  = HERE.parent / "table_legends.md"

GROUPS = [
    ("Manuscript tables — design", ["T1", "T2"]),
    ("Manuscript tables — fitness contrast", ["T3", "T4", "T5", "T6", "T7"]),
    ("Manuscript tables — stress contrast", ["T8", "T9", "T10"]),
    ("Supplementary tables", ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]),
]

wb = openpyxl.load_workbook(WB)
caps = {}
for ws in wb:
    tag = ws.title.split()[0]
    rng = list(ws.merged_cells.ranges)[0]
    caps[tag] = (ws.title, ws.cell(row=rng.min_row, column=1).value.strip())

lines = ["<!-- PROSE -->",
         "<!-- Captions are authored in draft/derived/make_workbook.py and lifted",
         "     here by make_table_legends.py.  Edit them THERE, not here, then run",
         "         python3 draft/derived/make_workbook.py",
         "         python3 draft/derived/make_table_legends.py -->",
         "",
         "# Table legends",
         "",
         "One table per tab in `draft/derived/TABLES.xlsx`, in this order.  Each",
         "caption also sits below its own table in the workbook, merged to the",
         "width of that table.",
         ""]

for head, tags in GROUPS:
    lines += [f"## {head}", ""]
    for t in tags:
        tab, cap = caps[t]
        body = re.sub(r"^(Table S?\d+\.\s+[^.]+\.)", r"**\1**", cap, count=1)
        lines += [f"`{tab}`", "", body, ""]

OUT.write_text("\n".join(lines) + "\n")
print(f"wrote {OUT}  ({len(caps)} legends)")
