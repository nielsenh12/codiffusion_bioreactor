#!/usr/bin/env python
"""
backfill_taxonomy.py — fill each ASV's midas_taxonomy + rel_ab from the FULL new_data taxonomy
(new_data/taxonomy_codif_all.csv, 100% of the 3,276 ASVs) into the alignment-hits outputs.

The alignment run annotated MiDAS taxonomy from model_inputs/taxonomy.csv (~91% coverage); this patches
those two annotation fields to 100% coverage. The BV-BRC alignment hits themselves are unchanged.
"""
import csv, json
from pathlib import Path

DIR = Path("/home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits")
TAX = Path("/home/freiburger/Documents/codiffusion_bioreactor/new_data/taxonomy_codif_all.csv")
RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

def main():
    full = {}
    for r in csv.DictReader(open(TAX)):
        lineage = " ".join(r[k] for k in RANKS if r.get(k))
        try:
            rel = float(r.get("rel_ab", "") or 0.0)
        except ValueError:
            rel = 0.0
        full[r["seq"]] = (lineage, rel)

    jpath = DIR / "asv_top20_alignment_hits.json"
    mapping = json.load(open(jpath))
    before = sum(1 for a in mapping if mapping[a].get("midas_taxonomy"))
    filled = 0
    for asv, rec in mapping.items():
        if asv in full:
            lineage, rel = full[asv]
            rec["midas_taxonomy"] = lineage
            rec["rel_ab"] = rel
            filled += 1
    json.dump(mapping, open(jpath, "w"))
    after = sum(1 for a in mapping if mapping[a].get("midas_taxonomy"))
    print(f"mapping: {len(mapping)} ASVs | midas_taxonomy non-empty {before} -> {after} "
          f"| in new_data taxonomy {filled} ({100*filled/len(mapping):.1f}%)")

    # summary CSV
    cpath = DIR / "asv_alignment_summary.csv"
    rows = list(csv.DictReader(open(cpath)))
    for row in rows:
        if row["asv"] in full:
            row["midas_taxonomy"], rel = full[row["asv"]]
            row["rel_ab"] = rel
    with open(cpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"summary.csv: {len(rows)} rows back-filled")

if __name__ == "__main__":
    main()
