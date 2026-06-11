# Large files (git-ignored)

The files below exceed **GitHub's 100 MB per-file limit**, so they are **git-ignored** (never committed)
to keep the repo pushable. They are large data artifacts, all **regenerable** from the notebooks/scripts in
this repo. Listed here for reference so their absence from a fresh clone is expected, not a mistake.

| file | size | what it is | regenerate |
|---|---|---|---|
| `model_inputs/BV_BRC_16S.json` | 2.0 GB | all BV-BRC 16S rRNA sequences `{fasta_header: sequence}` | `BV_BRC-Copy1.ipynb` cells 11–12 (BV-BRC API pull → merge) |
| `model_inputs/16S_metadata.json` | 957 MB | BV-BRC 16S feature metadata `{feature_id: record}` | `BV_BRC-Copy1.ipynb` cells 4–5 |
| `model_inputs/16S_md5_seq.json` | 496 MB | deduplicated unique 16S sequences `{md5: sequence}` | derived from `BV_BRC_16S.json` |

Ignored via `.gitignore` patterns: `model_inputs/BV_BRC_16S*.json`, `16S_metadata*.json`, `*md5*json`.
(`model_inputs/16S_md5_ID.json`, 81 MB, is under 100 MB but also `*md5*json`-ignored as a related artifact.)

> Threshold note: this excludes files **> 100 MB** (GitHub's hard file-size limit). No file in the repo is
> anywhere near 100 GB.
