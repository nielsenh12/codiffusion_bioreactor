"""Rebuild model_inputs/* for ALL ASVs from the canonical new_data/ dataset.

Replaces the old >0.1%-abundance-filtered inputs (261 ASVs, 16 samples) with the
full new_data codif_all set (3,276 ASVs, 29 samples), so every ASV receives an
iterativeID.

Inputs (canonical new_data):
  new_data/table_rel_codif_all.csv     wide rel-abundance (seq x 29 samples, PERCENT)
  new_data/taxonomy_codif_all.csv      taxonomy per seq (Kingdom..Species, + rel_ab)
  new_data/16s_metadata_codif_all.xlsx sample -> Day, Description, Media

Outputs (overwrite model_inputs/):
  taxonomy.json            { seq: { Kingdom..Species } }              (3,276)
  iterativeIDs.json        { seq: "<lowest_named_taxon>.<N>" }        (3,276)
  iterativeID_levels.json  { "<taxon>.<N>": "<rank>" }                (3,276)
  abundances.json          { sample: { seq: rel_ab_fraction } }       (29 samples)
  total.csv                long: seq, sample, rel_ab(%), taxonomy, date, media, timepoint, notes

Design decisions:
  * iterativeID naming replicates data_processing.py's committed logic exactly:
    walk Kingdom->...->Genus, keep the lowest rank in an UNBROKEN string lineage
    (never Species); collisions get an incrementing numeric suffix.
  * STABILITY: the 261 seq->iterativeID mappings already published (and cited in the
    manuscript, e.g. Methanobacterium.1) are PRESERVED verbatim; new ASVs continue
    the per-taxon numbering from where the existing IDs leave off. All 261 existing
    seqs are present in new_data, so none are orphaned.
"""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("model_inputs")
NEW = Path("new_data")
WIDE = NEW / "table_rel_codif_all.csv"
TAXO_CSV = NEW / "taxonomy_codif_all.csv"
META = NEW / "16s_metadata_codif_all.xlsx"

TAXO_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# The analysis sample scope is the original mature-community set (days 0, 140-315):
# inoculum (5-10) plus the post-retention-screen samples. The 13 new_data early
# ramp-up samples (A-Ga, days 21-133) are EXCLUDED so correlation/network results
# reflect the mature community, matching the manuscript. iterativeIDs still cover
# ALL 3,276 ASVs; only the per-sample abundance tables are scoped to these columns.
MATURE_SAMPLES = ["5-10", "H", "Ha", "I", "Ia", "J", "Ja", "Jb",
                  "K", "L", "M", "N", "Na", "P", "Q", "Qa"]


def clean(v):
    """Normalise a taxonomy cell to a python str or None."""
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    s = str(v).strip()
    if s in ("", "nan", "NaN", "NA", "None", "unclassified"):
        return None
    return s


def best_level_and_name(taxa):
    """data_processing.py naming logic: lowest rank in an unbroken string
    lineage from Kingdom, never descending to Species."""
    best_level = "Kingdom"
    for level in TAXO_LEVELS:
        taxon = taxa[level]
        if not isinstance(taxon, str) or level == "Species":
            break
        best_level = level
    return best_level, taxa[best_level]


def main():
    # ---- 1. Wide rel-abundance (percent -> fraction) ----
    wide = pd.read_csv(WIDE).set_index("seq")
    wide.columns = [str(c) for c in wide.columns]
    # Scope to the mature-community samples (drop the new early ramp-up samples).
    missing_samp = [s for s in MATURE_SAMPLES if s not in wide.columns]
    if missing_samp:
        print(f"WARNING: mature samples absent from wide table: {missing_samp}")
    sample_cols = [s for s in MATURE_SAMPLES if s in wide.columns]
    dropped = [c for c in wide.columns if c not in sample_cols]
    wide = wide[sample_cols]
    print(f"sample scope: keeping {len(sample_cols)} mature samples; "
          f"dropped {len(dropped)} early samples ({dropped})")
    max_val = float(wide.values.max())
    is_percent = max_val > 1.5
    if is_percent:
        wide_frac = wide / 100.0
    else:
        wide_frac = wide
    print(f"wide: {wide.shape[0]} seqs x {len(sample_cols)} samples; "
          f"max={max_val:.3f} -> treated as {'PERCENT (/100)' if is_percent else 'FRACTION'}")

    # ---- 2. Taxonomy ----
    tax_df = pd.read_csv(TAXO_CSV).drop_duplicates(subset=["seq"]).set_index("seq")
    # Ensure every wide seq has a taxonomy row
    missing = set(wide.index) - set(tax_df.index)
    if missing:
        print(f"WARNING: {len(missing)} wide seqs lack taxonomy; padding with None")
        pad = pd.DataFrame({c: [None] * len(missing) for c in TAXO_LEVELS},
                           index=list(missing))
        tax_df = pd.concat([tax_df, pad])

    taxonomy_dict = {}
    for seq in wide.index:
        row = tax_df.loc[seq]
        taxonomy_dict[seq] = {l: clean(row[l]) for l in TAXO_LEVELS}

    # Phylum backfill (matches the committed pipeline): any ASV with a Kingdom but
    # no Phylum gets Phylum = "<Kingdom> (unclassified)" so it still receives a
    # phylum-level iterativeID and a phylum colour (e.g. "Bacteria (unclassified)").
    backfilled = 0
    for seq, taxa in taxonomy_dict.items():
        if taxa["Phylum"] is None and taxa["Kingdom"] is not None:
            taxa["Phylum"] = f"{taxa['Kingdom']} (unclassified)"
            backfilled += 1
    print(f"phylum-backfilled {backfilled} ASVs to '<Kingdom> (unclassified)'")

    json.dump(taxonomy_dict, open(ROOT / "taxonomy.json", "w"), allow_nan=False, indent=4)
    print(f"wrote taxonomy.json ({len(taxonomy_dict)} seqs)")

    # ---- 3. iterativeIDs: preserve the published 261, extend to all ----
    # Seed from the committed 261-ASV stable seed (NOT the live iterativeIDs.json,
    # which this script overwrites -- reading that back would re-preserve prior runs
    # and defeat idempotency).
    seed_ids = ROOT / "iterativeIDs_stable_seed.json"
    seed_levels = ROOT / "iterativeID_levels_stable_seed.json"
    existing_ids = json.load(open(seed_ids if seed_ids.exists() else ROOT / "iterativeIDs.json"))
    existing_levels = json.load(open(seed_levels if seed_levels.exists() else ROOT / "iterativeID_levels.json"))

    iterativeIDs = {}          # seq -> ID
    iterativeID_levels = {}    # ID -> rank
    used = set()
    counters = defaultdict(int)  # taxon name -> max suffix used

    # 3a. lock in preserved assignments (only for seqs still present)
    preserved = 0
    for seq, ID in existing_ids.items():
        if seq not in taxonomy_dict:
            continue
        iterativeIDs[seq] = ID
        iterativeID_levels[ID] = existing_levels.get(ID, best_level_and_name(taxonomy_dict[seq])[0])
        used.add(ID)
        name, suf = ID.rsplit(".", 1)
        counters[name] = max(counters[name], int(suf))
        preserved += 1
    print(f"preserved {preserved} existing iterativeIDs")

    # 3b. diagnostic: how many preserved IDs would be renamed under fresh logic?
    renamed = 0
    for seq, ID in list(iterativeIDs.items()):
        _, fresh_name = best_level_and_name(taxonomy_dict[seq])
        if fresh_name is not None and ID.rsplit(".", 1)[0] != fresh_name:
            renamed += 1
    print(f"  ({renamed}/{preserved} preserved IDs have a taxon name that differs from "
          f"fresh new_data taxonomy; kept old for manuscript stability)")

    # 3c. assign new IDs to remaining seqs, in wide-table order (deterministic)
    new_assigned = 0
    for seq in wide.index:
        if seq in iterativeIDs:
            continue
        level, name = best_level_and_name(taxonomy_dict[seq])
        if name is None:
            name, level = "Unknown", "Unknown"
        counters[name] += 1
        ID = f"{name}.{counters[name]}"
        while ID in used:
            counters[name] += 1
            ID = f"{name}.{counters[name]}"
        iterativeIDs[seq] = ID
        iterativeID_levels[ID] = level
        used.add(ID)
        new_assigned += 1

    assert len(iterativeIDs) == wide.shape[0] == len(taxonomy_dict)
    assert len(set(iterativeIDs.values())) == len(iterativeIDs), "duplicate iterativeIDs!"
    json.dump(iterativeIDs, open(ROOT / "iterativeIDs.json", "w"))
    json.dump(iterativeID_levels, open(ROOT / "iterativeID_levels.json", "w"))
    print(f"wrote iterativeIDs.json ({len(iterativeIDs)}: {preserved} preserved + "
          f"{new_assigned} new) and iterativeID_levels.json ({len(iterativeID_levels)})")

    # ---- 4. abundances.json (sample -> seq -> fraction) ----
    # Scope the abundance table to ASVs actually present (>0) in >=1 mature sample;
    # the ~700 ASVs seen only in the excluded early samples keep their iterativeID
    # (full labelling) but are dropped here so they don't inflate the mature-
    # community figures.
    present = wide.index[(wide_frac[sample_cols].fillna(0) > 0).any(axis=1)]
    dropped_absent = len(wide.index) - len(present)
    abundances = {}
    for sample in sample_cols:
        abundances[sample] = {seq: float(wide_frac.at[seq, sample]) for seq in present
                              if pd.notna(wide_frac.at[seq, sample])}
    json.dump(abundances, open(ROOT / "abundances.json", "w"), indent=2)
    print(f"wrote abundances.json ({len(abundances)} samples x {len(present)} ASVs; "
          f"{dropped_absent} labelled-but-absent ASVs excluded from abundance table)")

    # ---- 5. metadata for total.csv ----
    meta_df = pd.read_excel(META, sheet_name="16s_metadata_2026")
    meta_df["sample"] = meta_df["sample"].astype(str)
    sample_to_day = dict(zip(meta_df["sample"], meta_df["Day"]))
    sample_to_desc = dict(zip(meta_df["sample"], meta_df["Description"]))
    sample_to_media = dict(zip(meta_df["sample"], meta_df["Media"]))

    # ---- 6. total.csv (long) --- present ASVs only, mature samples ----
    rows = []
    for seq in present:
        tax = taxonomy_dict[seq]
        for sample in sample_cols:
            v = wide.at[seq, sample]  # percent units, like the old total.csv
            rows.append({
                "seq": seq,
                "sample": sample,
                "rel_ab": float(v) if pd.notna(v) else 0.0,
                **{l: tax[l] for l in TAXO_LEVELS},
                "date": "",
                "media": sample_to_media.get(sample, ""),
                "timepoint": sample_to_day.get(sample, np.nan),
                "notes": sample_to_desc.get(sample, ""),
            })
    cols = ["seq", "sample", "rel_ab"] + TAXO_LEVELS + ["date", "media", "timepoint", "notes"]
    pd.DataFrame(rows)[cols].to_csv(ROOT / "total.csv", index=False)
    print(f"wrote total.csv ({len(rows)} rows = {len(present)} present seqs x {len(sample_cols)} samples)")

    # ---- 7. summary ----
    print("\n=== iterativeID coverage summary ===")
    by_level = defaultdict(int)
    for ID, lvl in iterativeID_levels.items():
        by_level[lvl] += 1
    for lvl in TAXO_LEVELS + ["Unknown"]:
        if by_level.get(lvl):
            print(f"  {lvl}: {by_level[lvl]}")
    # spot-check a few manuscript-cited labels remain on their expected seqs
    check = ["Methanobacterium.1", "Methanobacterium.2", "Methanobacteriaceae.1",
             "Enterobacteriaceae.1", "Burkholderiales.1", "Lentimicrobium.1",
             "Petrimonas.1", "Proteiniphilum.1", "Thermovirga.1", "Desulfovibrio.1"]
    id_to_seq = {v: k for k, v in iterativeIDs.items()}
    old_id_to_seq = {v: k for k, v in existing_ids.items()}
    print("  manuscript label stability:")
    for lbl in check:
        ok = id_to_seq.get(lbl) == old_id_to_seq.get(lbl) and lbl in id_to_seq
        print(f"    {lbl}: {'STABLE' if ok else 'CHANGED/absent'}")


if __name__ == "__main__":
    main()
