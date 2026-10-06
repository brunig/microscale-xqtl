#!/usr/bin/env python3
"""Collate the figure PNGs into one US Letter PDF, one figure per page.

    python3 draft/figures/make_figure_pdf.py

Main sequence first, then the alternates that make_figures.R also emits, so a
reader can flip through what would be submitted and then what was considered.
Nothing here regenerates a figure; run make_figures.R for that.
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
GEN  = HERE / "generated"
OUT  = GEN / "AFxQTL_all_figures.pdf"

DPI    = 150
PAGE   = (int(8.5 * DPI), int(11 * DPI))
MARGIN = int(0.6 * DPI)
FOOT   = int(0.32 * DPI)          # room for the running label

PAGES = [
    (HERE / "Figure1_design_v2.png",                      "Figure 1  —  design (blueprint render, not final art)"),
    (GEN / "Figure2.png",                                 "Figure 2"),
    (GEN / "Figure3.png",                                 "Figure 3  —  default encoding (mako, 50 kb) pending choice"),
    (GEN / "Figure4.png",                                 "Figure 4"),
    (GEN / "Figure3_OPTIONS.png",                         "OPTIONS  —  Figure 3, three encodings x three bin sizes"),
]

def font(sz):
    for p in ("/System/Library/Fonts/Supplemental/Arial.ttf",
              "/Library/Fonts/Arial.ttf", "/System/Library/Fonts/Helvetica.ttc"):
        if Path(p).exists():
            try:
                return ImageFont.truetype(p, sz)
            except OSError:
                pass
    return ImageFont.load_default()

F = font(15)
pages, missing = [], []
for path, label in PAGES:
    if not path.exists():
        missing.append(path.name)
        continue
    im = Image.open(path).convert("RGB")
    box = (PAGE[0] - 2 * MARGIN, PAGE[1] - 2 * MARGIN - FOOT)
    s = min(box[0] / im.width, box[1] / im.height)
    im = im.resize((max(1, round(im.width * s)), max(1, round(im.height * s))), Image.LANCZOS)
    page = Image.new("RGB", PAGE, "white")
    page.paste(im, ((PAGE[0] - im.width) // 2, MARGIN + (box[1] - im.height) // 2))
    d = ImageDraw.Draw(page)
    d.text((MARGIN, PAGE[1] - MARGIN - FOOT // 2), label, fill=(90, 90, 90), font=F)
    d.text((PAGE[0] - MARGIN, PAGE[1] - MARGIN - FOOT // 2),
           f"{len(pages) + 1}", fill=(150, 150, 150), font=F, anchor="ra")
    pages.append(page)
    print(f"  {path.name:46s} {im.width}x{im.height} px on page {len(pages)}")

if not pages:
    raise SystemExit("no figures found")
pages[0].save(OUT, "PDF", resolution=DPI, save_all=True, append_images=pages[1:])
print(f"\nwrote {OUT}  ({len(pages)} pages)")
if missing:
    print("MISSING, not in the PDF:", ", ".join(missing))
