# Workflow: Populating the Node/Edge Tables for the Co-occurrence Plot

This documents the step-by-step pipeline that turned the Louvain module
assignment and z-P analysis into the `nodes1.csv` / `edges1.csv` (originally
`nodes.csv`/`edges.csv`) tables consumed by the plotting scripts
(`plot_circos_signed.py` and earlier force-directed versions).

## Prerequisite outputs (from earlier steps)

- `asv_community_assignments.csv` — Louvain module per ASV (see the Louvain
  module README)
- `zP_analysis.csv` — z-score, participation coefficient, and role per ASV
  (see the z-P analysis README)
- `iterativeIDs.json` — a dict mapping ASV hash → a human-readable
  taxonomic/iterative ID (261 entries; not every ASV in the network has one)
- `Rel_ab_codif.csv` — taxonomy (Kingdom, Phylum, Class, Order, Family,
  Genus, Species) and per-sample relative abundance, one row per ASV
- `Sample_ID_to_day_of_operation.xlsx` — maps each sample column name to its
  day of bioreactor operation (used to identify and exclude the Day 0
  sample)

## Step 1 — Compute average relative abundance, excluding Day 0

`Rel_ab_codif.csv` has one relative-abundance column per sample:
`46152, H, Ha, I, Ia, J, Ja, Jb, K, L, M, N, Na, P, Q, Qa`. Cross-referencing
`Sample_ID_to_day_of_operation.xlsx` confirms column `46152` corresponds to
Day 0 (the `5-10` sample ID, mangled into a date-like value on export); the
other 15 columns correspond to Days 140–315.

```python
day_samples = ['H','Ha','I','Ia','J','Ja','Jb','K','L','M','N','Na','P','Q','Qa']
rel_ab['avg_rel_ab_excl_day0'] = rel_ab[day_samples].mean(axis=1)
```

## Step 2 — Merge z-P results with taxonomy and abundance

`zP_analysis.csv` was merged with `Rel_ab_codif.csv` (join key: ASV hash /
`seq`) to attach `Kingdom`, `Phylum`, `Genus`, and `avg_rel_ab_excl_day0` to
every ASV.

## Step 3 — Build a display name for each node (`Node` column)

Not every ASV has an `iterative_ID`. The naming fallback, in priority order:

1. **`iterative_ID`** from `iterativeIDs.json`, if present (261 of 266 ASVs
   in earlier stages; ultimately 121 of 266 in the final network after
   cross-referencing — see note below).
2. **Genus**, if `iterative_ID` is missing but a genus-level taxonomic call
   is available.
3. **Next available higher taxonomic rank** (Family → Order → Class →
   Phylum → Kingdom), for the small number of ASVs with neither an
   `iterative_ID` nor a genus-level classification (23 ASVs fell into this
   case — e.g. `Christensenellaceae.1`, `Bacteria.1`).

Duplicate names (e.g. multiple ASVs falling back to the same genus or
family) were disambiguated with a `.N` numeric suffix, checked against
already-used names (including existing `iterative_ID` values like
`Methanobacterium.4`) to guarantee every node label in the final network is
unique.

## Step 4 — Assemble `nodes.csv`

Final columns: `Node, Phylum, Domain, Relative_abundance, Role, Module`

| Column | Source |
|---|---|
| `Node` | Step 3 |
| `Phylum` | `Rel_ab_codif.csv` |
| `Domain` | `Kingdom` column from `Rel_ab_codif.csv` |
| `Relative_abundance` | `avg_rel_ab_excl_day0` from Step 1 |
| `Role` | `role` from `zP_analysis.csv` |
| `Module` | `community` from `asv_community_assignments.csv` |

266 rows — one per ASV in the positive co-occurrence network.

## Step 5 — Build `edges.csv`

The positive-only significant SparCC pairs (`sparcc_fdr_significant_pairs_r03.csv`,
filtered to `correlation > 0`) were translated from ASV hash IDs to the
`Node` display names built in Step 3, using the same ASV→Node mapping.

Final columns: `Source, Target, rho` — 1,269 rows, `rho` = the (always
positive) SparCC correlation coefficient.

For the signed circos plot specifically, a second edge table
(`edges_negative.csv`) was built the same way from the negative-correlation
pairs (949 pairs, 920 successfully mapped to nodes in the network; 29
involved ASVs outside the positive-only node set and were dropped), and used
as an additional, separately-styled edge layer (see the circos plot script).

## Step 6 — Manual correction: `nodes.csv` → `nodes1.csv`

The auto-generated `nodes.csv` did not carry a "Dominant methanogen"
designation — this label existed in an earlier/legacy version of the node
table (used by the original plotting code) but had no equivalent field in
the z-P `role` output, since it isn't a z-P category — it's a
manually-curated abundance/taxonomy-based call-out for the most abundant
methanogens. The user re-added this designation for 3 ASVs directly in the
node table and re-uploaded it as **`nodes1.csv`**, which is the file the
plotting scripts read from. `edges.csv` required no equivalent manual step
and is used as-is (referenced as `edges_new.csv` in the plotting scripts).

## Final files consumed by the plotting scripts

| File | Rows | Used by |
|---|---|---|
| `nodes1.csv` | 266 | All plotting scripts (`plot_circos_signed.py`, etc.) |
| `edges_new.csv` | 1,269 | Positive-edge layer |
| `edges_negative.csv` | 920 | Negative-edge layer (signed circos plot only) |

## Reproducing this pipeline end-to-end

```python
import pandas as pd
import json

# Step 1: relative abundance excluding Day 0
rel_ab = pd.read_csv("Rel_ab_codif.csv")
day_samples = ['H','Ha','I','Ia','J','Ja','Jb','K','L','M','N','Na','P','Q','Qa']
rel_ab["avg_rel_ab_excl_day0"] = rel_ab[day_samples].mean(axis=1)

# Step 2: merge
zp = pd.read_csv("zP_analysis.csv")
merged = zp.merge(
    rel_ab[["seq", "Kingdom", "Phylum", "Genus", "avg_rel_ab_excl_day0"]],
    left_on="ASV", right_on="seq", how="left"
)

# Step 3: build Node names (fallback chain + disambiguation) -- see full
# implementation in the project notebook/chat history for the disambiguation
# logic; omitted here for brevity.

# Step 4: assemble nodes.csv
nodes_out = pd.DataFrame({
    "Node": merged["Node"],
    "Phylum": merged["Phylum"],
    "Domain": merged["Kingdom"],
    "Relative_abundance": merged["avg_rel_ab_excl_day0"],
    "Role": merged["role"],
    "Module": merged["community"],
})

# Step 5: assemble edges.csv (positive-only, translated to Node names)
sig = pd.read_csv("sparcc_fdr_significant_pairs_r03.csv")
pos = sig[sig["correlation"] > 0].copy()
asv_to_node = dict(zip(merged["ASV"], merged["Node"]))
pos["Source"] = pos["ASV_1"].map(asv_to_node)
pos["Target"] = pos["ASV_2"].map(asv_to_node)
edges_out = pos[["Source", "Target", "correlation"]].rename(columns={"correlation": "rho"})
```

## Notes / caveats

- 23 of 266 nodes are labeled by a taxonomic rank above genus (Family or
  higher), since they had no `iterative_ID` and no genus-level
  classification. These labels are less specific than the rest and should
  be read accordingly in any figure.
- The "Dominant methanogen" designation in `Role` (3 ASVs) is a manual
  addition, not derived from the z-P calculation — it coexists in the same
  column as the 7 z-P role categories but comes from a different source.
  Downstream code that filters or colors by `Role` should be aware these are
  two different classification schemes sharing one column.
