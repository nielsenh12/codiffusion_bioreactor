# SparCC Co-occurrence Network Analysis

This repository documents the full analysis pipeline used to identify statistically significant, strongly-correlated ASV pairs from 16S amplicon count data using SparCC (Sparse Correlations for Compositional data; Friedman & Alm, 2012).

## Method attribution

The correlation algorithm (`sparcc_core.py`) is a faithful Python 3 translation of the original SparCC method. It is **not** a reimplementation from the paper's equations — it is a direct, line-for-line syntax port (Python 2 → 3: print statements, `xrange` → `range`, tuple-based fancy indexing) of the original authors' own code, sourced from:

- [`bio-developer/sparcc`](https://github.com/bio-developer/sparcc) — a mirror of the original `yonatanf/sparcc` (Friedman & Alm) repository
- [`shafferm/fast_sparCC`](https://github.com/shafferm/fast_sparCC) — whose core functions are explicitly attributed by that repo's own comments as copied from `https://bitbucket.org/yonatanf/pysurvey/` (the original authors' source)

No algorithmic changes were made to the correlation estimation, exclusion-iteration logic, or bootstrap resampling procedure. The only functional addition (see `sparcc_core.py`) is clipping minor floating-point overshoot past ±1 on the correlation estimate, which occurs occasionally with small sample sizes and is a known, expected numerical artifact rather than an algorithmic change.

Reference: Friedman J, Alm EJ (2012) Inferring Correlation Networks from Genomic Survey Data. *PLoS Comput Biol* 8(9): e1002687.

## Repository contents

| File | Description |
|---|---|
| `sparcc_core.py` | Core algorithm module (correlation estimation + bootstrap resampling) |
| `sparcc_counts.tsv` | Original input: raw ASV counts, 2,571 ASVs × 16 samples (rows=ASVs, columns=day-of-operation) |
| `sparcc_counts_final.tsv` | Preprocessed input actually used for analysis: 309 ASVs × 13 samples |
| `cor_sparcc_final.tsv` | Real-data SparCC correlation matrix (309 × 309) |
| `permutations_final.tar.gz` | 100 bootstrap-resampled datasets (null model for significance testing) |
| `perm_cor_final.tar.gz` | SparCC correlation matrix computed on each of the 100 permutations |
| `sparcc_final_results_v2.csv` | All ASV pairs with real correlation + permutation-based p-value |
| `sparcc_final_results_with_fdr.csv` | Same, with Benjamini-Hochberg FDR-adjusted q-value added |
| `sparcc_fdr_significant_pairs.csv` | Final significant pairs: FDR q ≤ 0.01 (2,240 pairs) |
| `sparcc_fdr_significant_pairs_r03.csv` | Final reported network: FDR q ≤ 0.01 **and** \|r\| ≥ 0.3 (2,218 pairs) |

## Requirements

```
python >= 3.9
numpy
pandas
scipy   # only needed to independently verify the FDR correction; not required to run the pipeline
```

## Pipeline

### Step 0: Preprocessing

Three adjustments were made to the raw count data before correlation analysis, based on diagnostic checks described in the Notes section below:

1. **Dropped Day 0** — temporally distant from all other timepoints (baseline/pre-shift community, fundamentally different composition from Days 140-315).
2. **Dropped Days 147 and 174** — the two lowest-depth samples (28,376 and 22,677 reads respectively), roughly 10-15x lower depth than most other samples.
3. **Rarefied** all 13 remaining samples to 31,614 reads (the next-lowest depth after removing the two outliers above), to equalize detection probability across samples.
4. **Prevalence-filtered** to ASVs present (nonzero) in at least 3 of the 13 remaining samples, removing singleton/near-absent ASVs that cannot support a meaningful correlation estimate.

```python
import pandas as pd
import numpy as np

df = pd.read_csv('sparcc_counts.tsv', sep='\t', index_col=0)
df.columns = [int(c) for c in df.columns]

# Drop Day 0 and the two lowest-depth samples
df_dropped = df.drop(columns=[0, 174, 147])

# Rarefy remaining samples to the next-lowest depth
RAREFY_DEPTH = 31614
rng = np.random.default_rng(42)

rarefied = pd.DataFrame(index=df_dropped.index, columns=df_dropped.columns, dtype=int)
for col in df_dropped.columns:
    counts = df_dropped[col].values.astype(np.int64)
    sampled = rng.multivariate_hypergeometric(counts, RAREFY_DEPTH)
    rarefied[col] = sampled

# Prevalence filter: keep ASVs present in >= 3 of 13 samples
n_samples_present = (rarefied > 0).sum(axis=1)
filtered = rarefied[n_samples_present >= 3].copy()

filtered.index.name = 'ASV_ID'
filtered.to_csv('sparcc_counts_final.tsv', sep='\t')
```

**Input:** `sparcc_counts.tsv` (2,571 ASVs × 16 samples)
**Output:** `sparcc_counts_final.tsv` (309 ASVs × 13 samples)

### Step 1: Compute real SparCC correlations

```python
import pandas as pd
import sys
sys.path.insert(0, '.')
from sparcc_core import sparcc

df = pd.read_csv('sparcc_counts_final.tsv', sep='\t', index_col=0)
counts = df.T  # sparcc() expects rows=samples, columns=ASVs

cor_df, cov_df = sparcc(counts, iters=20, th=0.1, xiter=10, verbose=False)
cor_df.to_csv('cor_sparcc_final.tsv', sep='\t')
```

**Input:** `sparcc_counts_final.tsv` + `sparcc_core.py`
**Output:** `cor_sparcc_final.tsv`

### Step 2: Generate bootstrap permutations

Following the original method's default (100 permutations is the tool's built-in default and documented recommendation for empirical p-value estimation):

```python
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '.')
from sparcc_core import permute_w_replacement

df = pd.read_csv('sparcc_counts_final.tsv', sep='\t', index_col=0)
counts = df.T

N_PERM = 100
for i in range(N_PERM):
    np.random.seed(i)
    perm = permute_w_replacement(counts)
    perm_out = perm.T
    perm_out.index.name = 'ASV_ID'
    perm_out.to_csv(f'permutations_final/permutation_{i}.tsv', sep='\t')
```

**Input:** `sparcc_counts_final.tsv`
**Output:** `permutations_final/permutation_0.tsv` ... `permutation_99.tsv`

### Step 3: Run SparCC on each permutation

```python
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '.')
from sparcc_core import sparcc

iter_idx = int(sys.argv[1])  # run once per index, 0-99

df = pd.read_csv(f'permutations_final/permutation_{iter_idx}.tsv', sep='\t', index_col=0)
counts = df.T

cor_df, cov_df = sparcc(counts, iters=20, th=0.1, xiter=10, verbose=False)
np.save(f'iters_final/perm_cor_{iter_idx}.npy', cor_df.values.astype(np.float32))
```

**Input:** `permutations_final/permutation_<i>.tsv` (each) + `sparcc_core.py`
**Output:** `perm_cor_final.tar.gz` (containing `perm_cor_0.npy` ... `perm_cor_99.npy`)

### Step 4: Pseudo p-values and Benjamini-Hochberg FDR correction

```python
import pandas as pd
import numpy as np
import glob

# --- Pseudo p-values ---
cor_df = pd.read_csv('cor_sparcc_final.tsv', sep='\t', index_col=0)
asvs = cor_df.index.tolist()
n = len(asvs)

cor_vals = cor_df.values
iu = np.triu_indices(n, k=1)
real_cor = cor_vals[iu]

perm_files = sorted(glob.glob('iters_final/perm_cor_*.npy'),
                     key=lambda x: int(x.split('_')[-1].split('.')[0]))
n_perm = len(perm_files)
perm_array = np.empty((n_perm, len(real_cor)), dtype=np.float32)
for i, f in enumerate(perm_files):
    perm_array[i] = np.load(f)[iu]

# Two-sided pseudo p-value: fraction of permutations at least as extreme as the real correlation
abs_real = np.abs(real_cor)
abs_perm = np.abs(perm_array)
ge_count = (abs_perm >= abs_real[np.newaxis, :]).sum(axis=0)
pval = ge_count / n_perm

asv_index = pd.DataFrame({'index': range(n), 'ASV_ID': asvs})
idx_to_id = dict(zip(asv_index['index'], asv_index['ASV_ID']))

edge_df = pd.DataFrame({
    'ASV_1_idx': iu[0], 'ASV_2_idx': iu[1],
    'correlation': real_cor.round(4), 'pvalue': pval
})
edge_df['ASV_1'] = edge_df['ASV_1_idx'].map(idx_to_id)
edge_df['ASV_2'] = edge_df['ASV_2_idx'].map(idx_to_id)
final = edge_df[['ASV_1', 'ASV_2', 'correlation', 'pvalue']].sort_values('correlation', key=abs, ascending=False)
final.to_csv('sparcc_final_results_v2.csv', index=False)

# --- Benjamini-Hochberg FDR correction ---
df = pd.read_csv('sparcc_final_results_v2.csv')
n_tests = len(df)

df_sorted = df.sort_values('pvalue').reset_index(drop=True)
ranks = np.arange(1, n_tests + 1)
bh_q = df_sorted['pvalue'] * n_tests / ranks

# Enforce monotonicity (required step of the BH procedure)
bh_q_adj = bh_q.values[::-1]
bh_q_adj = np.minimum.accumulate(bh_q_adj)[::-1]
bh_q_adj = np.clip(bh_q_adj, 0, 1)

df_sorted['fdr_qvalue'] = bh_q_adj
df_sorted.to_csv('sparcc_final_results_with_fdr.csv', index=False)

# --- Final filtering ---
sig = df_sorted[df_sorted['fdr_qvalue'] <= 0.01].copy()
sig.to_csv('sparcc_fdr_significant_pairs.csv', index=False)

sig_r03 = sig[sig['correlation'].abs() >= 0.3].copy()
sig_r03.to_csv('sparcc_fdr_significant_pairs_r03.csv', index=False)
```

**Input:** `cor_sparcc_final.tsv` + all `iters_final/perm_cor_<i>.npy`
**Output:** `sparcc_final_results_v2.csv` → `sparcc_final_results_with_fdr.csv` → `sparcc_fdr_significant_pairs.csv` → `sparcc_fdr_significant_pairs_r03.csv`

## Validation

The manual Benjamini-Hochberg implementation above was independently cross-checked against `scipy.stats.false_discovery_control` (method='bh') on the full 47,586-pair result set. The two implementations produced identical q-values (max absolute difference: 0.0) and identical significant-pair counts at every threshold tested.

## Final result

**2,218 ASV pairs** met both criteria: FDR q ≤ 0.01 and \|r\| ≥ 0.3 (the correlation magnitude threshold used in the original SparCC publication's example networks). Of these, 1,269 are positive correlations and 949 are negative.

## Notes on preprocessing decisions

- **Day 0 exclusion:** diagnostic checks showed 1,876 of the original 2,571 ASVs (73%) were detected only at Day 0 and never again afterward, indicating a substantial community shift between Day 0 and Day 140. Retaining Day 0 caused several ASVs with a shared "absent-then-blooming" trajectory to appear strongly (and likely spuriously) anti-correlated with whatever was dominant at baseline.
- **Rarefaction:** raw sample depths ranged from 22,677 to 631,616 reads (~28x difference). Several apparent correlation "hub" ASVs were found to be undetectable specifically on the lowest-depth samples regardless of day, consistent with a sequencing-depth detection artifact rather than true biological absence.
- **Correlation magnitude threshold:** given the small sample size (n=13), the standard error of a correlation estimate is approximately 1/√(n-3) ≈ 0.32, comparable in magnitude to the field-standard 0.3 cutoff itself. The combined FDR + magnitude filter was used to avoid over-relying on either criterion alone.
