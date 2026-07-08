"""Produce _regen_newdata.ipynb from _exec_run.ipynb, migrated to the new_data
ASV set (all 3,276 ASVs labelled) while PRESERVING the original 16-sample
mature-community scope (5-10, H..Qa; days 0,140-315).

Edits:
  * cells 1,3,5  -> load-only (consume the already-built model_inputs foundation
                   instead of rebuilding from the 050526 table + >0.1% filter).
                   cell 3 scopes taxonomy_js to the mature-present ASVs so the
                   taxonomic bar charts reflect the mature community.
  * cells 19,23,25 -> pandas-2.x fix: groupby(axis=1) -> .T.groupby().T
  * cell 14      -> harden the asv_IDs table against stale/absent labels.
  * cell 23      -> define taxonomy_series locally (cross-cell dep in marimo).
The 16-sample sample_days/_days dicts in cells 10,13,15,18 are LEFT UNCHANGED
(the mature scope). Everything else reads model_inputs/* from disk.
"""
import json
import ast
from pathlib import Path

SRC = Path("_exec_run.ipynb")
OUT = Path("_regen_newdata.ipynb")

CELL1 = '''# --- MIGRATED: operational bridge already rebuilt for the 16 mature samples ---
from json import load
summary_samples = load(open("model_inputs/measurements/summary_samples.json", "r"))
print(f"summary_samples.json: {len(summary_samples)} mature samples")
'''

CELL3 = '''# --- MIGRATED: load the new_data foundation (all 3,276 ASVs labelled; abundance
# table scoped to the 16 mature samples and the ASVs present in them) ---
from json import load, dump
from pandas import DataFrame
import numpy as np
ab_js = load(open("model_inputs/abundances.json", "r"))          # {sample: {seq: fraction}}  16 mature samples
full_taxonomy = load(open("model_inputs/taxonomy.json", "r"))    # {seq: {Kingdom..Species}}  all 3,276
present = set(next(iter(ab_js.values())).keys())                 # ASVs present in the mature community
taxonomy_js = {seq: full_taxonomy[seq] for seq in present}       # scope taxonomic charts to the mature community
meaningful_ASVs = DataFrame(ab_js)                               # index=seq, cols=mature samples
print(f"abundances: {meaningful_ASVs.shape[0]} present ASVs x {meaningful_ASVs.shape[1]} mature samples; "
      f"taxonomy(all): {len(full_taxonomy)} ASVs")
'''

CELL5 = '''# --- MIGRATED: load iterativeIDs (261 published labels preserved + extended to all 3,276) ---
from json import load
iterativeIDs = load(open("model_inputs/iterativeIDs.json", "r"))            # {seq: ID}  all 3,276
iterativeID_levels = load(open("model_inputs/iterativeID_levels.json", "r"))  # {ID: rank}
print(f"iterativeIDs: {len(iterativeIDs)} ASVs labelled")
'''

GROUPBY_OLD = 'df = df.groupby(df.columns, axis=1).sum()'
GROUPBY_NEW = 'df = df.T.groupby(df.columns).sum().T'


def code_cells(nb):
    return [i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code"]


def harden_cell14(src):
    return src.replace(
        'for iterativeID, taxonomy in all_top10.items():',
        'for iterativeID, taxonomy in all_top10.items():\n'
        '    if iterativeID not in IDasvs:  # stale label absent from current map\n'
        '        continue')


def comment_orphaned_dicts(src):
    """Some marimo cells have a commented dict-open (`# sample_days = {`) followed
    by an UN-commented body -- a syntax error once we stop rewriting those dicts.
    Comment the orphaned body lines through the closing brace."""
    out, commenting, depth = [], False, 0
    for ln in src.split('\n'):
        s = ln.lstrip()
        if not commenting:
            if s.startswith('#') and '=' in ln and ln.count('{') > ln.count('}'):
                commenting = True
                depth = ln.count('{') - ln.count('}')
            out.append(ln)
        else:
            out.append(ln if s.startswith('#') else '# ' + ln)
            depth += ln.count('{') - ln.count('}')
            if depth <= 0:
                commenting = False
    return '\n'.join(out)


def fix_cell23_taxseries(src):
    lines = src.split('\n')
    lines = [l for l in lines if l.strip() != 'taxonomy_series = Series(taxonomies)']
    end = next(j for j, l in enumerate(lines) if 'for col in df.columns}' in l)
    lines = lines[:end + 1] + ['taxonomy_series = Series(taxonomies)'] + lines[end + 1:]
    return '\n'.join(lines)


def main():
    nb = json.load(open(SRC))
    idx = code_cells(nb)

    def setcell(ordinal, text):
        nb["cells"][idx[ordinal]]["source"] = text.splitlines(keepends=True)

    setcell(1, CELL1)
    setcell(3, CELL3)
    setcell(5, CELL5)

    # groupby fix (cells 19, 23, 25 contain the pattern; replace wherever present)
    ng = 0
    for i in idx:
        s = "".join(nb["cells"][i]["source"])
        if GROUPBY_OLD in s:
            nb["cells"][i]["source"] = s.replace(GROUPBY_OLD, GROUPBY_NEW).splitlines(keepends=True)
            ng += 1

    # cell 14 harden
    c14 = "".join(nb["cells"][idx[14]]["source"])
    nb["cells"][idx[14]]["source"] = harden_cell14(c14).splitlines(keepends=True)

    # cell 23 taxonomy_series (after groupby replace)
    c23 = "".join(nb["cells"][idx[23]]["source"])
    nb["cells"][idx[23]]["source"] = fix_cell23_taxseries(c23).splitlines(keepends=True)

    # clear outputs for a clean execute; fix orphaned commented dicts; ast-check
    bad = 0
    for c in nb["cells"]:
        if c["cell_type"] == "code":
            c["outputs"] = []
            c["execution_count"] = None
            fixed = comment_orphaned_dicts("".join(c["source"]))
            c["source"] = fixed.splitlines(keepends=True)
            try:
                ast.parse("".join(c["source"]))
            except SyntaxError as e:
                bad += 1
                print(f"SYNTAX ERROR: {e}")

    json.dump(nb, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}: neutered cells 1/3/5, groupby-fixed {ng} cells, "
          f"hardened cell14, fixed cell23; sample_days LEFT at 16 mature; syntax-bad={bad}")


if __name__ == "__main__":
    main()
