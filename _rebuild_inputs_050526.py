"""Rebuild canonical model_inputs/* files from the 050526 raw data.

Inputs:
  model_inputs/table_rel_codif_output_050526.csv   wide rel-abundance (seq x sample)
  model_inputs/table1_codif_output_050526.csv      long w/ taxonomy   (Q-only is fine; we use it solely for taxonomy per seq)
  model_inputs/16s_metadata_codif_050526.xlsx      sample -> Day, Description, Media

Outputs (overwrite):
  model_inputs/abundances.json               { sample: { seq: rel_ab_fraction } }
  model_inputs/total.csv                     long: seq, sample, rel_ab, taxonomy..., Day, media, notes
  model_inputs/taxonomy.json                 { seq: { Kingdom, Phylum, ..., Species } }
  model_inputs/iterativeIDs.json             { seq: "<lowest_named_taxon>.<N>" }
  model_inputs/iterativeID_levels.json       { "<lowest_named_taxon>.<N>": "<rank>" }
  model_inputs/measurements/summary_samples.json   re-keyed to new sample names where Day matches
"""
import pandas as pd
import numpy as np
from json import load, dump
from collections import defaultdict
from pathlib import Path

ROOT = Path("model_inputs")
WIDE = ROOT / "table_rel_codif_output_050526.csv"
LONG = ROOT / "table1_codif_output_050526.csv"
META = ROOT / "16s_metadata_codif_050526.xlsx"

TAXO_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def main():
    # --- 1. Load metadata ---
    meta_df = pd.read_excel(META, sheet_name="16s_metadata_2026")
    sample_to_day = dict(zip(meta_df["sample"].astype(str), meta_df["Day"]))
    sample_to_desc = dict(zip(meta_df["sample"].astype(str), meta_df["Description"]))
    sample_to_media = dict(zip(meta_df["sample"].astype(str), meta_df["Media"]))
    print(f"metadata: {len(meta_df)} samples, days = {sorted(sample_to_day.values())}")

    # --- 2. Wide rel-abundance ---
    wide_df = pd.read_csv(WIDE)
    wide_df = wide_df.set_index("seq")
    wide_df.columns = [str(c) for c in wide_df.columns]
    sample_cols = list(wide_df.columns)
    print(f"wide: {wide_df.shape[0]} seqs x {wide_df.shape[1]} samples ({sample_cols})")

    # New input has rel_ab as a fraction or percentage? Check magnitude.
    max_val = float(wide_df.values.max())
    print(f"wide max value = {max_val}")
    is_percentage = max_val > 1.5
    print("treating wide table as", "PERCENT (will divide by 100)" if is_percentage else "FRACTION")
    if is_percentage:
        wide_df = wide_df / 100.0

    # --- 3. Taxonomy from long file (1 row per seq, Q-only is fine) ---
    long_df = pd.read_csv(LONG)
    tax_per_seq = (
        long_df.drop_duplicates(subset=["seq"])
               .set_index("seq")[TAXO_LEVELS]
    )
    # Ensure all wide-table seqs have taxonomy (missing → all None)
    missing = set(wide_df.index) - set(tax_per_seq.index)
    if missing:
        print(f"WARNING: {len(missing)} seqs in wide but absent from long taxonomy; padding with None")
        pad = pd.DataFrame({c: [None] * len(missing) for c in TAXO_LEVELS}, index=list(missing))
        tax_per_seq = pd.concat([tax_per_seq, pad])

    # Replace pandas NA / "NA" string with python None for JSON-friendliness
    tax_per_seq = tax_per_seq.where(tax_per_seq.notna() & (tax_per_seq != "NA"), None)

    # --- 4. Build taxonomy.json ---
    taxonomy_dict = {}
    for seq, row in tax_per_seq.iterrows():
        taxonomy_dict[seq] = {l: row[l] for l in TAXO_LEVELS}
    dump(taxonomy_dict, open(ROOT / "taxonomy.json", "w"))
    print(f"wrote taxonomy.json ({len(taxonomy_dict)} seqs)")

    # --- 5. Build iterativeIDs.json + iterativeID_levels.json ---
    # For each seq, pick the lowest taxon (Species → ... → Kingdom) that is non-null/non-"NA",
    # then assign a unique numeric suffix per name.
    counters = defaultdict(int)
    iterativeIDs = {}
    iterativeID_levels = {}
    for seq, row in tax_per_seq.iterrows():
        chosen_name = None
        chosen_level = "Unknown"
        for lvl in reversed(TAXO_LEVELS):  # Species, Genus, ..., Kingdom
            val = row[lvl]
            if val is None or (isinstance(val, float) and np.isnan(val)) or val == "NA" or val == "":
                continue
            chosen_name = str(val)
            chosen_level = lvl
            break
        if chosen_name is None:
            chosen_name = "Unknown"
            chosen_level = "Unknown"
        counters[chosen_name] += 1
        unique_id = f"{chosen_name}.{counters[chosen_name]}"
        iterativeIDs[seq] = unique_id
        iterativeID_levels[unique_id] = chosen_level
    dump(iterativeIDs, open(ROOT / "iterativeIDs.json", "w"))
    dump(iterativeID_levels, open(ROOT / "iterativeID_levels.json", "w"))
    print(f"wrote iterativeIDs.json + iterativeID_levels.json ({len(iterativeIDs)} entries)")

    # --- 6. Build abundances.json (sample -> seq -> rel_ab fraction) ---
    abundances = {}
    for sample in sample_cols:
        abundances[sample] = {seq: float(v) for seq, v in wide_df[sample].items() if pd.notna(v)}
    dump(abundances, open(ROOT / "abundances.json", "w"))
    print(f"wrote abundances.json ({len(abundances)} samples)")

    # --- 7. Build total.csv (long: seq, sample, rel_ab, taxonomy, Day, media, notes) ---
    total_rows = []
    for seq, tax_row in tax_per_seq.iterrows():
        for sample in sample_cols:
            rel_ab = float(wide_df.at[seq, sample]) if seq in wide_df.index else 0.0
            row = {
                "seq": seq,
                "sample": sample,
                "rel_ab": rel_ab * 100.0 if not is_percentage else rel_ab * 100.0,  # keep % units like the old total.csv
                **{l: tax_row[l] for l in TAXO_LEVELS},
                "date": "",
                "media": sample_to_media.get(sample, ""),
                "timepoint": sample_to_day.get(sample, np.nan),
                "notes": sample_to_desc.get(sample, ""),
            }
            total_rows.append(row)
    total_df = pd.DataFrame(total_rows)
    # Match the OLD total.csv column order
    cols = ["seq", "sample", "rel_ab"] + TAXO_LEVELS + ["date", "media", "timepoint", "notes"]
    total_df = total_df[cols]
    total_df.to_csv(ROOT / "total.csv", index=False)
    print(f"wrote total.csv ({len(total_df)} rows = {wide_df.shape[0]} seqs x {len(sample_cols)} samples)")

    # --- 8. Re-key summary_samples.json by Day match ---
    # Load existing operational measurements (keyed by old sample names)
    old_summary = load(open(ROOT / "_backup_pre_050526" / "summary_samples.json"))
    old_meta_lookup = {}  # old sample name -> day (need a separate source)
    # Old day mapping was hardcoded; reconstruct it from inspection
    old_sample_days = {
        "10AB": 0, "A34": 22, "B34": 35, "C3": 69, "D34": 83, "E34": 91,
        "F34": 105, "G12": 126, "G3": 126, "H34": 140, "I34": 149,
        "J12": 161, "J34": 126, "K12": 182, "L12": 189, "M12": 210,
        "N12": 238, "P12": 283, "Q": 300,
    }
    day_to_old = {}
    for s, d in old_sample_days.items():
        day_to_old.setdefault(d, s)  # first old sample at each day wins
    new_summary = {}
    for new_sample, day in sample_to_day.items():
        old_sample = day_to_old.get(int(day))
        if old_sample and old_sample in old_summary:
            new_summary[new_sample] = old_summary[old_sample]
            print(f"  {new_sample} (day {day}) <- {old_sample}")
        else:
            print(f"  {new_sample} (day {day}) -- no matching operational data, skipping")
    dump(new_summary, open(ROOT / "measurements" / "summary_samples.json", "w"))
    print(f"wrote summary_samples.json ({len(new_summary)} samples with operational data)")

    # --- 9. Print suggested sample_days dict for the notebook ---
    print("\nSuggested sample_days dict for notebook:")
    print("sample_days = {")
    for sample, day in sample_to_day.items():
        commented = "" if sample in new_summary else "    # "
        print(f'    {commented}"{sample}": {int(day)},')
    print("}")


if __name__ == "__main__":
    main()
