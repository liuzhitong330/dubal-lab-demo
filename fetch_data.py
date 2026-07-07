"""
Pull Allen Mouse Brain Atlas ISH expression data for:
1. Aging suppressor / longevity genes: Kl, Fgf21, Igf1, Foxo3, Sirt1, Tert, Gdf11, Gdf15
2. X-chromosome escapee genes (expressed from both X copies):
   Kdm6a, Kdm5c, Ddx3x, Eif2s3x, Usp9x, Xist

Structures chosen to span the brain and include choroid plexus
(the primary brain source of klotho):
  - Lateral ventricle choroid plexus (CP-lat, id 220)
  - Isocortex (315)
  - Hippocampal formation (HPF, 1089)
  - Cerebellum (CB, 512)
  - Hypothalamus (HY, 1097)
  - Striatum / Caudoputamen (CP-str, 672)
  - Olfactory bulb (MOB, 507)
  - Thalamus (TH, 549)
"""
import requests, csv, time
from pathlib import Path

OUT = Path(__file__).parent
BASE = "http://api.brain-map.org/api/v2/data/query.json"

AGING_GENES = ["Kl", "Fgf21", "Igf1", "Foxo3", "Sirt1", "Tert", "Gdf11", "Gdf15"]
XESCAPE_GENES = ["Kdm6a", "Kdm5c", "Ddx3x", "Eif2s3x", "Usp9x"]  # Xist is nuclear RNA, skip

STRUCTURES = {
    "Choroid_plexus": 220,
    "Isocortex":      315,
    "Hippocampus":    1089,
    "Cerebellum":     512,
    "Hypothalamus":   1097,
    "Striatum":       672,
    "Olfact_bulb":    507,
    "Thalamus":       549,
}

def get_dataset(gene):
    """Find best ISH dataset for a gene (sagittal preferred)."""
    url = (f"{BASE}?criteria=model::SectionDataSet,"
           f"rma::criteria,genes[acronym$eq'{gene}'],"
           f"[failed$eqfalse],"
           f"products[abbreviation$eqMouse]"
           f"&include=genes&num_rows=5&start_row=0")
    r = requests.get(url, timeout=20).json()
    if not r.get("msg"):
        return None
    # prefer sagittal plane (plane_of_section_id=2)
    for row in r["msg"]:
        if row.get("plane_of_section_id") == 2:
            return row["id"]
    return r["msg"][0]["id"] if r["msg"] else None

def get_expression(dataset_id, struct_ids):
    """Get expression_energy per structure for a dataset."""
    id_str = ",".join(str(i) for i in struct_ids)
    url = (f"{BASE}?criteria=model::StructureUnionize,"
           f"rma::criteria,[section_data_set_id$eq{dataset_id}],"
           f"[structure_id$in{id_str}]"
           f"&num_rows=50")
    r = requests.get(url, timeout=20).json()
    out = {}
    for row in r.get("msg", []):
        out[row["structure_id"]] = row.get("expression_energy", 0.0) or 0.0
    return out

all_genes = AGING_GENES + XESCAPE_GENES
struct_ids = list(STRUCTURES.values())
id_to_name = {v: k for k, v in STRUCTURES.items()}

rows = []
for gene in all_genes:
    print(f"  {gene}...", end=" ", flush=True)
    ds = get_dataset(gene)
    if ds is None:
        print("no dataset")
        continue
    expr = get_expression(ds, struct_ids)
    row = {"gene": gene, "dataset_id": ds, "group": "aging" if gene in AGING_GENES else "x_escape"}
    for sid, name in id_to_name.items():
        row[name] = round(expr.get(sid, 0.0), 4)
    rows.append(row)
    print(f"dataset {ds} ✓")
    time.sleep(0.3)

# Write TSV
fieldnames = ["gene", "group", "dataset_id"] + list(STRUCTURES.keys())
with open(OUT / "expression.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"\nSaved expression.tsv ({len(rows)} genes)")
for r in rows:
    vals = {k: r[k] for k in STRUCTURES}
    top = max(vals, key=vals.get)
    print(f"  {r['gene']:12s}  peak={top} ({vals[top]:.3f})")
