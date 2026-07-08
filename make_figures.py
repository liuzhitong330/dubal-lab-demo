"""
Generate two SVG figures for the Dubal lab demo:
1. aging_heatmap.svg  — aging suppressor genes across brain regions (row-normalized)
2. xescape_heatmap.svg — X-chromosome escapee genes across brain regions
"""
import csv, math
from pathlib import Path

OUT = Path(__file__).parent
STRUCTURES = ["Midbrain","Isocortex","Hippocampus","Cerebellum",
              "Hypothalamus","Striatum","Olfact_bulb","Thalamus"]
STRUCT_LABELS = ["Midbrain","Isocortex","Hippocampus","Cerebellum",
                 "Hypothalamus","Striatum","Olf. bulb","Thalamus"]

# ── Load data ──────────────────────────────────────────────────────────────────
aging_rows, xescape_rows = [], []
with open(OUT / "expression.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        vals = [float(r[s]) for s in STRUCTURES]
        entry = {"gene": r["gene"], "vals": vals}
        if r["group"] == "aging":
            aging_rows.append(entry)
        elif r["group"] == "x_escape" and max(vals) > 0.01:
            xescape_rows.append(entry)


def row_norm(vals):
    mn, mx = min(vals), max(vals)
    if mx == mn:
        return [0.0] * len(vals)
    return [(v - mn) / (mx - mn) for v in vals]


def blue_red(v):
    """0→deep blue, 0.5→white, 1→deep red"""
    if v < 0.5:
        t = v * 2
        r = int(255 * t + 230 * (1 - t))
        g = int(255 * t + 230 * (1 - t))
        b = int(255 * (1 - t) + 255 * t) if False else int(255 - 100 * t)
        # low end: white(255,255,255) → blue(40,80,200)
        r = int(230 + (40 - 230) * (1 - t * 2 + 1 - 1) )  # redo cleanly
        # v=0 → (220,235,255), v=0.5 → (255,255,255)
        frac = v / 0.5
        r = int(220 + (255 - 220) * frac)
        g = int(235 + (255 - 235) * frac)
        b = 255
    else:
        frac = (v - 0.5) / 0.5
        r = 255
        g = int(255 + (60 - 255) * frac)
        b = int(255 + (60 - 255) * frac)
    return f"rgb({r},{g},{b})"


def warm_cold(v):
    """0→cool blue, 1→warm red via white midpoint"""
    if v <= 0.5:
        t = v / 0.5
        r = int(180 + (255 - 180) * t)
        g = int(210 + (255 - 210) * t)
        b = 255
    else:
        t = (v - 0.5) / 0.5
        r = 255
        g = int(255 + (80 - 255) * t)
        b = int(255 + (80 - 255) * t)
    return f"rgb({r},{g},{b})"


# ── Figure 1: aging suppressor heatmap ────────────────────────────────────────
CELL_W = 62
CELL_H = 28
LEFT   = 72
TOP    = 95
n_genes = len(aging_rows)
n_structs = len(STRUCTURES)
W = LEFT + n_structs * CELL_W + 80
H = TOP + n_genes * CELL_H + 80

cells = ""
for gi, row in enumerate(aging_rows):
    nv = row_norm(row["vals"])
    for si, v in enumerate(nv):
        x = LEFT + si * CELL_W
        y = TOP + gi * CELL_H
        color = warm_cold(v)
        raw = row["vals"][si]
        cells += (
            f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" '
            f'fill="{color}" stroke="white" stroke-width="1.5"/>'
        )
        # annotate value if notable
        if raw > 0.5 or (raw > 0 and v > 0.7):
            txt_col = "#333" if v < 0.8 else "white"
            cells += (
                f'<text x="{x+CELL_W/2:.0f}" y="{y+CELL_H/2+4:.0f}" '
                f'text-anchor="middle" font-size="8.5" fill="{txt_col}">{raw:.2f}</text>'
            )

# column headers
col_hdrs = ""
for si, lbl in enumerate(STRUCT_LABELS):
    x = LEFT + si * CELL_W + CELL_W // 2
    lines = lbl.split("\n")
    y0 = TOP - 8 - (len(lines) - 1) * 11
    # highlight cerebellum column (Kl peak)
    CB_SI = STRUCTURES.index("Cerebellum")
    if si == CB_SI:
        col_hdrs += (
            f'<rect x="{LEFT}" y="{TOP-40}" width="{CELL_W}" height="40" '
            f'fill="#fff8e1" rx="3"/>'
            f'<text x="{x}" y="{TOP-2}" text-anchor="middle" font-size="9.5" '
            f'fill="#b07800" font-weight="bold">'
        )
    else:
        col_hdrs += f'<text x="{x}" y="{TOP-2}" text-anchor="middle" font-size="9.5" fill="#444">'
    col_hdrs += lbl.replace("\n", "</text>")
    if "\n" in lbl:
        parts = lbl.split("\n")
        col_hdrs = col_hdrs.rsplit("<text", 1)[0]  # redo
        for li, part in enumerate(parts):
            y = TOP - 4 - (len(parts) - 1 - li) * 12
            fw = "bold" if si == CB_SI else "normal"
            fc = "#b07800" if si == CB_SI else "#444"
            col_hdrs += (f'<text x="{x}" y="{y}" text-anchor="middle" font-size="9.5" '
                         f'fill="{fc}" font-weight="{fw}">{part}</text>')
    else:
        if si == CB_SI:
            col_hdrs += f'<text x="{x}" y="{TOP-8}" text-anchor="middle" font-size="9.5" fill="#b07800" font-weight="bold">{lbl}</text>'
        else:
            col_hdrs += f'<text x="{x}" y="{TOP-8}" text-anchor="middle" font-size="9.5" fill="#444">{lbl}</text>'

# row labels
row_lbls = ""
gene_colors = {
    "Kl": "#c0392b", "Fgf21": "#e67e22", "Igf1": "#27ae60",
    "Foxo3": "#8e44ad", "Sirt1": "#2980b9", "Tert": "#16a085",
    "Gdf11": "#d35400", "Gdf15": "#7f8c8d"
}
for gi, row in enumerate(aging_rows):
    y = TOP + gi * CELL_H + CELL_H // 2 + 4
    color = gene_colors.get(row["gene"], "#333")
    row_lbls += (f'<text x="{LEFT-6}" y="{y}" text-anchor="end" '
                 f'font-size="11" fill="{color}" font-weight="600" '
                 f'font-family="ui-monospace,Menlo,monospace">'
                 f'{row["gene"]}</text>')

# colorbar
cb_x = LEFT + n_structs * CELL_W + 12
cb_y = TOP
cb_h = n_genes * CELL_H
cb_w = 14
colorbar = f'<defs><linearGradient id="cb1" x1="0" y1="1" x2="0" y2="0">'
colorbar += f'<stop offset="0%" stop-color="{warm_cold(0)}"/>'
colorbar += f'<stop offset="50%" stop-color="{warm_cold(0.5)}"/>'
colorbar += f'<stop offset="100%" stop-color="{warm_cold(1)}"/>'
colorbar += f'</linearGradient></defs>'
colorbar += f'<rect x="{cb_x}" y="{cb_y}" width="{cb_w}" height="{cb_h}" fill="url(#cb1)" stroke="#ccc" stroke-width="0.5"/>'
colorbar += f'<text x="{cb_x+cb_w/2:.0f}" y="{cb_y-4}" text-anchor="middle" font-size="8" fill="#555">high</text>'
colorbar += f'<text x="{cb_x+cb_w/2:.0f}" y="{cb_y+cb_h+10}" text-anchor="middle" font-size="8" fill="#555">low</text>'

svg1 = f"""<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif; background:white;">
  <text x="{W//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Aging Suppressor Genes — ISH Expression Across Mouse Brain Regions
  </text>
  <text x="{W//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas · expression energy (row-normalized per gene) · values annotated where &gt;0.5
  </text>
  <text x="{W//2}" y="52" text-anchor="middle" font-size="10" fill="#b07800">
    ▶ Cerebellum highlighted: highest Kl (klotho) expression site in brain
  </text>
  {colorbar}
  {col_hdrs}
  {cells}
  {row_lbls}
</svg>"""

with open(OUT / "aging_heatmap.svg", "w") as f:
    f.write(svg1)
print("Wrote aging_heatmap.svg")


# ── Figure 2: X-escapee gene heatmap + bar inset ──────────────────────────────
# Show raw (not normalized) values to convey actual expression levels
# Also add an inset showing the sex-dosage concept

CELL_W2 = 68
CELL_H2 = 32
LEFT2   = 80
TOP2    = 90
n2 = len(xescape_rows)
W2 = LEFT2 + n_structs * CELL_W2 + 100
H2 = TOP2 + n2 * CELL_H2 + 110

# global max for shared color scale
all_vals2 = [v for row in xescape_rows for v in row["vals"]]
gmax = max(all_vals2)

def blue_scale(v, gmax):
    """white → deep teal, shared scale"""
    t = min(v / gmax, 1.0)
    r = int(255 - 180 * t)
    g = int(255 - 100 * t)
    b = int(255 - 40 * t)
    return f"rgb({r},{g},{b})"

cells2 = ""
for gi, row in enumerate(xescape_rows):
    for si, v in enumerate(row["vals"]):
        x = LEFT2 + si * CELL_W2
        y = TOP2 + gi * CELL_H2
        color = blue_scale(v, gmax)
        cells2 += (
            f'<rect x="{x}" y="{y}" width="{CELL_W2}" height="{CELL_H2}" '
            f'fill="{color}" stroke="white" stroke-width="1.5"/>'
        )
        if v > 1.0:
            cells2 += (
                f'<text x="{x+CELL_W2/2:.0f}" y="{y+CELL_H2/2+4:.0f}" '
                f'text-anchor="middle" font-size="8.5" fill="#333">{v:.1f}</text>'
            )

# col headers
col_hdrs2 = ""
for si, lbl in enumerate(STRUCT_LABELS):
    x = LEFT2 + si * CELL_W2 + CELL_W2 // 2
    parts = lbl.split("\n")
    for li, part in enumerate(parts):
        y = TOP2 - 4 - (len(parts) - 1 - li) * 12
        col_hdrs2 += (f'<text x="{x}" y="{y}" text-anchor="middle" '
                      f'font-size="9.5" fill="#444">{part}</text>')

# row labels
xescape_colors = {"Kdm6a": "#8e44ad", "Kdm5c": "#2980b9",
                  "Eif2s3x": "#27ae60", "Usp9x": "#c0392b"}
row_lbls2 = ""
for gi, row in enumerate(xescape_rows):
    y = TOP2 + gi * CELL_H2 + CELL_H2 // 2 + 4
    color = xescape_colors.get(row["gene"], "#333")
    row_lbls2 += (f'<text x="{LEFT2-6}" y="{y}" text-anchor="end" '
                  f'font-size="11" fill="{color}" font-weight="600" '
                  f'font-family="ui-monospace,Menlo,monospace">'
                  f'{row["gene"]}</text>')

# Colorbar
cb2_x = LEFT2 + n_structs * CELL_W2 + 12
cb2_h = n2 * CELL_H2
colorbar2 = f'<defs><linearGradient id="cb2" x1="0" y1="1" x2="0" y2="0">'
colorbar2 += f'<stop offset="0%" stop-color="{blue_scale(0, gmax)}"/>'
colorbar2 += f'<stop offset="100%" stop-color="{blue_scale(gmax, gmax)}"/>'
colorbar2 += f'</linearGradient></defs>'
colorbar2 += f'<rect x="{cb2_x}" y="{TOP2}" width="14" height="{cb2_h}" fill="url(#cb2)" stroke="#ccc" stroke-width="0.5"/>'
colorbar2 += f'<text x="{cb2_x+7}" y="{TOP2-4}" text-anchor="middle" font-size="8" fill="#555">{gmax:.1f}</text>'
colorbar2 += f'<text x="{cb2_x+7}" y="{TOP2+cb2_h+10}" text-anchor="middle" font-size="8" fill="#555">0</text>'

# Sex-dosage concept box at bottom
box_y = TOP2 + n2 * CELL_H2 + 20
box_x = LEFT2
box_w = n_structs * CELL_W2
SVG_NOTE = (
    f'<rect x="{box_x}" y="{box_y}" width="{box_w}" height="68" '
    f'fill="#f0f4ff" rx="5" stroke="#c3d0ee" stroke-width="1"/>'
    f'<text x="{box_x+12}" y="{box_y+17}" font-size="10.5" fill="#2c3e7a" font-weight="600">'
    f'Sex-chromosome dosage hypothesis</text>'
    f'<text x="{box_x+12}" y="{box_y+32}" font-size="9.5" fill="#444">'
    f'XY individuals: 1 copy of each escapee gene  →  baseline expression</text>'
    f'<text x="{box_x+12}" y="{box_y+46}" font-size="9.5" fill="#444">'
    f'XX individuals: 2 copies (both X active for escapees)  →  ~2× dosage</text>'
    f'<text x="{box_x+12}" y="{box_y+60}" font-size="9.5" fill="#555">'
    f'KDM6A (H3K27 demethylase) &amp; KDM5C (H3K4 demethylase) are chromatin regulators — '
    f'dosage differences may underlie sex-biased resilience to brain aging.</text>'
)

svg2 = f"""<svg viewBox="0 0 {W2} {H2+10}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif; background:white;">
  <text x="{W2//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    X-Chromosome Escapee Genes — Expression Across Mouse Brain Regions
  </text>
  <text x="{W2//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas · raw expression energy (shared color scale) · genes that escape X-inactivation
  </text>
  <text x="{W2//2}" y="53" text-anchor="middle" font-size="10" fill="#555">
    In XX individuals these genes are expressed from both X chromosomes — a potential molecular basis for sex-linked longevity
  </text>
  {colorbar2}
  {col_hdrs2}
  {cells2}
  {row_lbls2}
  {SVG_NOTE}
</svg>"""

with open(OUT / "xescape_heatmap.svg", "w") as f:
    f.write(svg2)
print("Wrote xescape_heatmap.svg")
