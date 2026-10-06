#!/usr/bin/env python3
"""Generate Figure 1 layout blueprint as SVG (editable in Illustrator/BioRender)."""
import math

W, H = 1040, 690
COL = dict(Glu="#5b6b7a", Ace="#0f766e", Caf="#a2560b", NaCl="#6d28d9")
INK, MUTE, LINE, PAPER = "#15181c", "#5c6672", "#b9c0c8", "#ffffff"
BBQ, WILD = "#1f4e79", "#d9a441"

s = []
A = s.append
A(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica Neue, Helvetica, Arial, sans-serif">')
A(f'<rect width="{W}" height="{H}" fill="{PAPER}"/>')
A('<defs>')
A(f'<marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{INK}"/></marker>')
A(f'<marker id="arm" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{MUTE}"/></marker>')
A('</defs>')

def txt(x, y, t, size=12, fill=INK, w="400", anchor="start", style=""):
    A(f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" font-weight="{w}" text-anchor="{anchor}"{style}>{t}</text>')

def panel(letter, x, y, title):
    txt(x, y, letter, 17, INK, "700")
    txt(x + 20, y, title, 13.5, INK, "600")

def arrow(x1, y1, x2, y2, mute=False):
    A(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{MUTE if mute else INK}" stroke-width="1.6" marker-end="url(#{"arm" if mute else "ar"})"/>')

def cell(cx, cy, r, fill, bar=None, barcol=INK):
    A(f'<ellipse cx="{cx}" cy="{cy}" rx="{r}" ry="{r*0.86}" fill="{fill}" stroke="{INK}" stroke-width="1.1"/>')
    if bar:
        for i, bc in enumerate(bar):
            off = (i - (len(bar) - 1) / 2) * 4.2
            A(f'<rect x="{cx+off-1.3}" y="{cy-r*0.5}" width="2.6" height="{r}" fill="{bc}" rx="1"/>')

# ---------------- Panel A : cross construction ----------------
panel("A", 22, 34, "Thirty-six wild isolates, each crossed to the same parent")
ay = 108
stations = [
    (86,  "Parents",            None),
    (218, "Diploid",            "select\ndiploids"),
    (348, "Tetrads",            "sporulate"),
    (478, "Segregant pool",     "select\nhaploids"),
]
# parents
cell(64, ay, 17, "#fdf6e6", [WILD])
cell(112, ay, 17, "#fdf6e6", [BBQ])
txt(88, ay + 6, "&#215;", 16, INK, "600", "middle")
txt(64, ay - 30, "wild", 11, MUTE, "500", "middle")
txt(64, ay - 19, "isolate", 11, MUTE, "500", "middle")
txt(112, ay - 30, "BBQ", 11, MUTE, "600", "middle")
txt(112, ay - 19, "shared", 11, MUTE, "500", "middle")
cell(218, ay, 18, "#fdf6e6", [WILD, BBQ])
for i, (dx, dy) in enumerate([(-13, -11), (13, -11), (-13, 11), (13, 11)]):
    cell(348 + dx, ay + dy, 10, "#fdf6e6", [WILD if i % 2 else BBQ])
A(f'<circle cx="348" cy="{ay}" r="30" fill="none" stroke="{INK}" stroke-width="1.1" stroke-dasharray="3 3"/>')
import random
random.seed(7)
for i in range(6):
    dx = [-26, -6, 14, -26, -6, 14][i]; dy = [-13, -15, -13, 12, 14, 12][i]
    cell(478 + dx, ay + dy, 10, "#fdf6e6", [random.choice([WILD, BBQ])])
for x0, x1, lab in [(142, 192, "select\ndiploids"), (246, 308, "sporulate"), (390, 436, "select\nhaploids")]:
    arrow(x0, ay, x1, ay)
    for j, ln in enumerate(lab.split("\n")):
        txt((x0 + x1) / 2, ay - 32 + j * 11, ln, 10.5, MUTE, "500", "middle")
txt(86, ay + 62, "one wild parent per cross", 11, MUTE, "400", "middle")
txt(478, ay + 62, "segregants, one pool per cross", 11, MUTE, "400", "middle")
# the x36 bracket
A(f'<path d="M20,{ay+78} L20,{ay+90} L536,{ay+90} L536,{ay+78}" fill="none" stroke="{INK}" stroke-width="1.3"/>')
txt(278, ay + 108, "repeated for 36 wild isolates — BBQ is the shared parent of every cross", 12, INK, "500", "middle")

# ---------------- Panel B : four plates ----------------
panel("B", 600, 34, "Four conditions, one 96-well plate each")
plates = [("Glu", "glucose (control)"), ("Ace", "acetic acid"), ("Caf", "caffeine"), ("NaCl", "sodium chloride")]
px0, py0, pw, ph = 608, 64, 196, 88
for i, (key, lab) in enumerate(plates):
    cx = px0 + (i % 2) * 214
    cy = py0 + (i // 2) * 118
    A(f'<rect x="{cx}" y="{cy}" width="{pw}" height="{ph}" rx="7" fill="none" stroke="{INK}" stroke-width="1.4"/>')
    for r in range(8):
        for c in range(12):
            wx = cx + 14 + c * 15.2
            wy = cy + 14 + r * 8.7
            A(f'<circle cx="{wx:.1f}" cy="{wy:.1f}" r="3.3" fill="{COL[key]}" fill-opacity="0.72"/>')
    txt(cx, cy - 6, lab, 11.5, COL[key], "600")
    txt(cx + pw, cy - 6, "96 wells", 10.5, MUTE, "400", "end")
txt(608, py0 + 252, "Every cross appears on every plate: 36 crosses × 4 conditions,", 11.5, MUTE, "400")
txt(608, py0 + 267, "a median of 2 wells per cross per plate (range 1 to 6).", 11.5, MUTE, "400")

# ---------------- Panel C : passage timeline ----------------
panel("C", 22, 380, "Two 48-h passages, one bottleneck")
tx0, tw = 96, 300
rows = [("Glu", "glucose"), ("Caf", "caffeine"), ("NaCl", "sodium chloride"), ("Ace", "acetic acid")]
ty = 418
for i, (key, lab) in enumerate(rows):
    y = ty + i * 40
    two = key != "Ace"
    end = tx0 + (tw if two else tw * 0.5)
    A(f'<line x1="{tx0}" y1="{y}" x2="{end}" y2="{y}" stroke="{COL[key]}" stroke-width="7" stroke-linecap="round" opacity="0.85"/>')
    txt(tx0 - 8, y + 4, lab, 11.5, COL[key], "600", "end")
    A(f'<circle cx="{tx0}" cy="{y}" r="5" fill="{PAPER}" stroke="{INK}" stroke-width="1.4"/>')
    mid = tx0 + tw * 0.5
    A(f'<circle cx="{mid}" cy="{y}" r="5" fill="{PAPER}" stroke="{INK}" stroke-width="1.4"/>')
    if two:
        A(f'<circle cx="{tx0+tw}" cy="{y}" r="5" fill="{INK}"/>')
    else:
        A(f'<circle cx="{mid}" cy="{y}" r="5" fill="{INK}"/>')
        txt(mid + 14, y + 4, "second timepoint not recovered", 10.5, MUTE, "400")
A(f'<line x1="{tx0+tw*0.5}" y1="{ty-16}" x2="{tx0+tw*0.5}" y2="{ty+3*40+16}" stroke="{INK}" stroke-width="1" stroke-dasharray="4 3"/>')
txt(tx0 + tw * 0.5, ty - 24, "bottleneck", 11, INK, "600", "middle")
txt(tx0, ty + 3 * 40 + 34, "0 h", 10.5, MUTE, "400", "middle")
txt(tx0 + tw * 0.5, ty + 3 * 40 + 34, "48 h", 10.5, MUTE, "400", "middle")
txt(tx0 + tw, ty + 3 * 40 + 34, "96 h", 10.5, MUTE, "400", "middle")
txt(22, ty + 3 * 40 + 60, "Growth in HCD medium; pools harvested at the filled timepoint.", 11.5, MUTE, "400")

# ---------------- Panel D : sequencing and the two contrasts ----------------
panel("D", 470, 380, "Allele frequencies give two contrasts per pool")
dx0 = 478
# flow
flow = [("all 384 wells", 540), ("pooled libraries", 672), ("30× target", 800)]
fy = 424
A(f'<rect x="{dx0}" y="{fy-14}" width="86" height="28" rx="4" fill="none" stroke="{INK}" stroke-width="1.2"/>')
txt(dx0 + 43, fy + 4, "384 wells", 11, INK, "500", "middle")
arrow(dx0 + 92, fy, dx0 + 122, fy)
A(f'<rect x="{dx0+126}" y="{fy-14}" width="96" height="28" rx="4" fill="none" stroke="{INK}" stroke-width="1.2"/>')
txt(dx0 + 174, fy + 4, "one library", 11, INK, "500", "middle")
arrow(dx0 + 228, fy, dx0 + 258, fy)
A(f'<rect x="{dx0+262}" y="{fy-14}" width="134" height="28" rx="4" fill="none" stroke="{INK}" stroke-width="1.2"/>')
txt(dx0 + 329, fy + 4, "one run, 30× target", 11, INK, "500", "middle")
arrow(dx0 + 402, fy, dx0 + 432, fy)
txt(dx0 + 470, fy + 4, "allele", 11, INK, "500", "middle")
txt(dx0 + 470, fy + 16, "frequencies", 11, INK, "500", "middle")

# two contrast tracks
def track(x, y, w, h, colour, label, sub, baseline_lab, peaks):
    A(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="none" stroke="{LINE}" stroke-width="1"/>')
    mid = y + h / 2
    A(f'<line x1="{x}" y1="{mid}" x2="{x+w}" y2="{mid}" stroke="{MUTE}" stroke-width="1" stroke-dasharray="4 3"/>')
    txt(x - 8, mid + 4, baseline_lab, 10.5, MUTE, "400", "end")
    pts = []
    import random as rnd
    rnd.seed(11)
    for i in range(0, w + 1, 4):
        v = rnd.uniform(-0.12, 0.12)
        for (pc, pa, pw_) in peaks:
            v += pa * math.exp(-((i - pc) ** 2) / (2 * pw_ ** 2))
        pts.append(f"{x+i:.1f},{mid - v*(h/2*0.88):.1f}")
    A(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{colour}" stroke-width="1.7" stroke-linejoin="round"/>')
    txt(x, y - 7, label, 11.5, INK, "600")
    txt(x + w, y - 7, sub, 10.5, MUTE, "400", "end")

tw2 = 400
track(dx0 + 56, 486, tw2, 62, COL["Glu"],
      "fitness contrast", "each pool against an allele frequency of 0.5", "0.5",
      [(150, 0.92, 22), (300, 0.42, 18)])
track(dx0 + 56, 584, tw2, 62, COL["NaCl"],
      "stress contrast", "each stress pool against its matched glucose pool", "0",
      [(300, 0.72, 18), (150, 0.06, 22)])
A(f'<line x1="{dx0+56+150}" y1="480" x2="{dx0+56+150}" y2="652" stroke="{INK}" stroke-width="0.9" stroke-dasharray="3 3" opacity="0.55"/>')
txt(dx0 + 56 + 150, 474, "chr VII", 10.5, INK, "600", "middle")
txt(dx0 + 56 + tw2, 664, "chromosome position", 10.5, MUTE, "400", "end")
A('</svg>')

out = "/Users/kruglyakbws/Desktop/AFxQTL/draft/figures/Figure1_design_v2.svg"
open(out, "w").write("\n".join(s))
print("wrote", out)
