"""
Correlate every parameter in summary_samples.json (across H34..Qa, day >= 140) with each ASV's
relative abundance at the iterativeID level. Writes one ASV_correlations_all_relevant_samples_<param>.json
per numeric parameter.

Day-matching is the linker: abundance samples (new naming: H, Ha, I, ..., Qa) and summary samples (mostly
old naming after consolidation: H34, Ha, I34, ..., Qa) are joined by Day.
"""
import os, re, glob
from json import load, dump
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ABUND_PATH = "model_inputs/abundances.json"
ITID_PATH = "model_inputs/iterativeIDs.json"
SUMMARY_JSON_PATH = "model_inputs/measurements/summary_samples.json"
SUMMARY_CSV_PATH = "model_inputs/measurements/Summary_interpolated.csv"
META_PATH = "model_inputs/16s_metadata_codif_050526.xlsx"
OUT_DIR = "modeling_files/correlations"
MIN_DAY = 140        # H34..Qa = day 140 onward
MIN_PAIRS = 5        # require at least 5 paired (day, value, ASV) samples

os.makedirs(OUT_DIR, exist_ok=True)


def slugify(s: str) -> str:
    """Make a filename-safe slug. Keep '%', alphanumerics, _, -. Collapse runs of separators."""
    s = re.sub(r"[\s/]+", "_", s.strip())
    s = re.sub(r"[^A-Za-z0-9%_\-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "param"


# --- Load everything ---
abundances = load(open(ABUND_PATH))
iterativeIDs = load(open(ITID_PATH))
summary_samples = load(open(SUMMARY_JSON_PATH))

new_meta = pd.read_excel(META_PATH, sheet_name="16s_metadata_2026")
abund_sample_to_day = dict(zip(new_meta["sample"].astype(str), new_meta["Day"].astype(int)))

# Summary key -> Day from Summary_interpolated.csv
summary_df = pd.read_csv(SUMMARY_CSV_PATH)
summary_day_map = {}
for _, row in summary_df.dropna(subset=["Biomass Sample ID"]).iterrows():
    biom = str(row["Biomass Sample ID"]).strip()
    day = int(round(float(row["Days of Operation"])))
    summary_day_map[biom] = day

# Restrict to relevant samples (Day >= 140)
abund_relevant = {s: d for s, d in abund_sample_to_day.items() if d >= MIN_DAY}
summary_relevant = {s: d for s, d in summary_day_map.items() if d >= MIN_DAY}

print(f"Relevant abundance samples (day >= {MIN_DAY}): {len(abund_relevant)}")
print(f"  {sorted(abund_relevant.items(), key=lambda x: x[1])}")
print(f"Relevant summary samples (day >= {MIN_DAY}): {len(summary_relevant)}")
print(f"  {sorted(summary_relevant.items(), key=lambda x: x[1])}")

# --- Build seq_hash -> day -> rel_ab nested dict ---
seq_to_day_abund = {}
for sample, inner in abundances.items():
    if sample not in abund_relevant:
        continue
    d = abund_relevant[sample]
    for seq, val in inner.items():
        seq_to_day_abund.setdefault(seq, {})[d] = val

# Filter to ASVs present (>0) in at least 5 samples
ASV_min_samples = 5
seq_to_day_abund = {
    seq: dd for seq, dd in seq_to_day_abund.items()
    if sum(1 for v in dd.values() if v and v > 0) >= ASV_min_samples
}
print(f"ASVs surviving '>0 in >={ASV_min_samples} samples' filter: {len(seq_to_day_abund)}")

# --- Determine which parameters are numeric (correlation is meaningful) ---
all_params = list(next(iter(summary_samples.values())).keys())
numeric_params = []
for p in all_params:
    vals = [content.get(p) for s, content in summary_samples.items() if s in summary_relevant]
    finite = [v for v in vals if isinstance(v, (int, float)) and not (isinstance(v, float) and np.isnan(v))]
    if len(finite) >= MIN_PAIRS:
        numeric_params.append(p)
print(f"Numeric parameters with >= {MIN_PAIRS} non-null values across relevant samples: {len(numeric_params)} / {len(all_params)}")

# --- Optional: clear stale ASV_correlations files from previous run before regenerating ---
stale = glob.glob(os.path.join(OUT_DIR, "ASV_correlations_all_relevant_samples_*.json"))
for f in stale:
    os.remove(f)
print(f"Removed {len(stale)} stale ASV_correlations files")

# --- Correlate each parameter against each ASV ---
metrics_summary = []
for param in numeric_params:
    metric_by_day = {}
    for s, content in summary_samples.items():
        if s not in summary_relevant:
            continue
        v = content.get(param)
        if v is None:
            continue
        if isinstance(v, float) and np.isnan(v):
            continue
        try:
            metric_by_day[summary_relevant[s]] = float(v)
        except (TypeError, ValueError):
            pass
    if len(metric_by_day) < MIN_PAIRS:
        continue
    metric_series = pd.Series(metric_by_day)

    correlations = {}
    constant_skipped = 0
    for seq, day_abund in seq_to_day_abund.items():
        abund_series = pd.Series(day_abund)
        a, b = abund_series.align(metric_series, join="inner")
        if len(a) < MIN_PAIRS:
            continue
        if len(set(a)) <= 1 or len(set(b)) <= 1:
            constant_skipped += 1
            continue
        rho, p = spearmanr(a, b)
        if np.isnan(rho):
            continue
        ID = iterativeIDs.get(seq, seq)
        correlations[ID] = {"correlation": float(rho), "p_value": float(p)}

    if not correlations:
        continue
    correlations = dict(sorted(correlations.items(), key=lambda x: -abs(x[1]["correlation"])))
    slug = slugify(param)
    out_path = os.path.join(OUT_DIR, f"ASV_correlations_all_relevant_samples_{slug}.json")
    dump(correlations, open(out_path, "w"))
    metrics_summary.append((slug, len(correlations), constant_skipped, len(metric_by_day)))

# --- Print summary table ---
print(f"\nWrote {len(metrics_summary)} ASV_correlations_all_relevant_samples_*.json files to {OUT_DIR}/")
print(f"\n{'param_slug':45s}  {'n_ASVs':>8s}  {'n_paired_samples':>16s}  {'const_skipped':>14s}")
print("-" * 90)
for slug, n_corr, n_const, n_samples in sorted(metrics_summary, key=lambda x: -x[1]):
    print(f"{slug[:45]:45s}  {n_corr:8d}  {n_samples:16d}  {n_const:14d}")
