# Manuscript / correlation reconciliation

**No manuscript edits are applied.** An earlier draft of this file proposed number
updates; they were reverted after verifying they were based on a flawed comparison
(details below). The manuscript `04.13.26.docx` is back to its original state.

## Does relabeling / the new_data switch change defined organisms' results? NO.

A pairwise Spearman ρ between two defined organisms depends only on their own two
abundance vectors — adding or labeling other ASVs cannot change it. Verified directly:
at an identical sample scope, the old (manuscript-era) and new_data abundance tables
give **bit-identical** correlations.

| pair | old-data, 15 samp | new-data, 15 samp | old, 12 samp (d149–300) | new, 12 samp |
|---|---|---|---|---|
| Methanobacterium.2 vs .1 | −0.75 | −0.75 | −0.83 | −0.83 |
| Methanobacterium.2 vs Methanobacteriaceae.1 | −0.31 | −0.31 | −0.29 | −0.29 |

The abundances of defined organisms are also essentially identical across datasets
(e.g. Methanobacterium.1 at day 300: 49.75% [050526] vs 49.66% [codif_all]; Thermovirga.1
is bit-identical). The only differences are <0.5% compositional shifts from a few extra
ASVs in the denominator. **The correlation *value* depends on the sample scope, not the
ASV set or the labels.**

## Why my earlier report looked like things changed (my error)

1. **Different sample window.** The manuscript's correlations came from specific
   day-windowed analyses (there were `ASV_correlations_Days 150-300_*` and
   `Days 189-238_*` files). I compared them against the regenerated pipeline's scope
   (15 post-inoculum samples), so scope differences looked like data changes. E.g.
   Methanobacterium.2 vs .1 is −0.75 at 15 samples but −0.83 at the 12-sample
   days-149–300 window; the manuscript's −0.78 sits in between.
2. **"Only 2 significant" → "135 significant" for H2** is purely a consequence of
   testing all 2,567 ASVs instead of 261 — it does **not** change any individual
   organism's ρ. Thermovirga.1's H2 correlation is ≈ +0.65 on both the old committed
   data (+0.67) and new_data (+0.63); the manuscript's stated −0.93 does not reproduce
   on either and came from a different windowed analysis.

## Genuine old_data → new_data differences (pre-date this work)

These are NOT caused by labeling; they reflect the different DADA2 runs, and they
already existed between the manuscript and the committed 050526 pipeline before any
of this work:

| quantity | old_data (3,215) ≈ manuscript | committed 050526 (2,571) | new_data (3,276) |
|---|---|---|---|
| Methanobacteriaceae, day 300 | 81.6% (ms: 80.2%) | 77.8% | 77.8% |
| Shannon, day 300 | 1.89 (ms: "low of 1.9") | 2.11 | 2.11 |
| Methanobacteriaceae, inoculum | 13.1% (ms: 11.3%) | 13.1% | 13.1% |
| total ASVs | 3,215 (ms: 3,216) | 2,571 | 3,276 |

So the committed pipeline the figures already used (050526) disagreed with the
manuscript's headline numbers *before* this task; new_data matches 050526. If you want
the manuscript to track the current/canonical data, the day-300 enrichment would read
77.8% and Shannon low 2.1 — **but that is a separate decision about which dataset the
paper reports, independent of the iterativeID work.** Tell me and I'll apply it.

## What this task actually delivered (unchanged, correct)

- iterativeIDs for **all 3,276 ASVs** (261 published labels preserved; manuscript-cited
  labels verified stable).
- All figures regenerated so every ASV shows a readable label (no raw md5 hashes),
  on the preserved mature-community sample scope.
- The organism-level science in the figures is the same as before — only the labeling
  coverage and (for the co-occurrence network) the ASV membership changed.
