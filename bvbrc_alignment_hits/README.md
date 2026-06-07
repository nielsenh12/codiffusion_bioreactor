# BV-BRC alignment hits for the codiffusion 16S ASVs

Two-stage **alignment** mapping of all **3,276** codiffusion V4–V5 ASVs
(`new_data/dna-sequences_codif_all.fasta`) against the deduplicated BV-BRC 16S reference DB
(**459,301** unique sequences ≥200 bp from `model_inputs/16S_md5_seq.json`):

1. **Stage 1 — edlib** infix (`HW`) edit-distance scan with an adaptive k-band over every reference;
   keep the top **K1 = 500** lowest-distance candidates per ASV.
2. **Stage 2 — Biopython** local Smith–Waterman (match +2 / mismatch −3 / gap −5,−2) rescoring of those
   500; keep the top **K2 = 20** by alignment score.

Fully parallel (`edlib_biopython_hits.py`, 60 workers): the run did **~12.8 CPU-hours** (edlib 37,952 s +
Biopython 8,086 s) in **788 s wall** — a **59.5× speedup** on the 64-thread Threadripper.

## Outputs

| file | contents |
|---|---|
| **`asv_top20_alignment_hits.json`** | per-ASV mapping: `asv_len, midas_taxonomy, rel_ab, best_align_score, best_identity, n_edlib_candidates`, and `top20[]` each with **`align_score`** (raw Smith–Waterman), **`identity`** (normalized %, matches/aligned_len), `n_matches`, `aligned_len`, `edlib_distance`, `edlib_identity`, `organism`, `genome_id`, `taxon_id`, `feature_id`, `md5`, `ref_seq_len`, `ref_aln_start/end`, and `lineage{Kingdom…Species}` (NCBI via taxopy). |
| `asv_alignment_summary.csv` | one row per ASV: best hit + key fields |
| `run_stats.json` | parameters, counts, per-stage CPU time, wall time, parallel speedup |

## Result summary

Best-hit % identity (Biopython, normalized): mean **0.974**, median **0.987**, p5 0.914.
**69.8%** of ASVs have a best hit ≥97% identity, **43.1%** ≥99%, only **3.0%** below 90%. The dominant
*Methanobacterium* ASV matches *Methanobacterium* sp. at 100% identity (edit distance 0).

Reproduce: `python edlib_biopython_hits.py --workers 60 --chunksize 1 --k1 500 --k2 20`
(reference DB + ASV/taxonomy paths are set at the top of the script).
