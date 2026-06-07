#!/usr/bin/env python
"""
build_new_data_taxonomy.py — parse new_data/taxonomy_codif_all.qza (QIIME2) + table_rel_codif_all.csv
into new_data/taxonomy_codif_all.csv (seq, Kingdom..Species, rel_ab) for all 3,276 codif_all ASVs.

rel_ab = each ASV's MAX relative abundance across the 29 samples, in PERCENT (table_rel columns already
sum to 100 per sample — do NOT rescale). Used to annotate the alignment-hits outputs at 100% coverage.
"""
import zipfile, csv
import pandas as pd
from pathlib import Path

NEW = Path("/home/freiburger/Documents/codiffusion_bioreactor/new_data")
RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
LET = {"k": "Kingdom", "d": "Kingdom", "p": "Phylum", "c": "Class",
       "o": "Order", "f": "Family", "g": "Genus", "s": "Species"}


def main():
    z = zipfile.ZipFile(NEW / "taxonomy_codif_all.qza")
    tsv = z.read([x for x in z.namelist() if x.endswith("taxonomy.tsv")][0]).decode().splitlines()
    rows = {}
    for line in tsv[1:]:
        p = line.split("\t")
        if len(p) < 2:
            continue
        d = {r: "" for r in RANKS}
        for tok in p[1].split(";"):
            tok = tok.strip()
            if len(tok) >= 3 and tok[1:3] == "__":
                rk = LET.get(tok[0].lower())
                if rk:
                    d[rk] = tok[3:].strip()
        rows[p[0]] = d
    rel = pd.read_csv(NEW / "table_rel_codif_all.csv")
    seqcol, samp = rel.columns[0], rel.columns[1:]
    relmax = {r[seqcol]: float(r[list(samp)].max()) for _, r in rel.iterrows()}   # already % per sample
    with open(NEW / "taxonomy_codif_all.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["seq"] + RANKS + ["rel_ab"]); w.writeheader()
        for asv, d in rows.items():
            w.writerow({"seq": asv, **d, "rel_ab": round(relmax.get(asv, 0.0), 6)})
    print(f"wrote new_data/taxonomy_codif_all.csv : {len(rows)} ASVs")


if __name__ == "__main__":
    main()
