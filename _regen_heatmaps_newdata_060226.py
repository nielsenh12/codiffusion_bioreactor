"""
Rebuild BOTH abundance heatmaps entirely from the NEW dataset (new_data/codif_all):
  - heatmap cells  : new-data relative abundances (table_rel_codif_all.csv, all 29 samples)
  - row labels      : iterativeIDs generated from the new ASV taxonomy (notebook logic, lines 203-220)
  - row colors      : phylum palette regenerated from the new taxonomy (notebook logic, lines 1810-1834)
  - Shannon track   : new-data H' (normalized), now 29 points (days 0-315)

Nothing in the canonical model_inputs/ is overwritten -- all new-data inputs are built in memory.
`taxonomy_linkage` and `_create_heatmap` are faithful copies from
data_processing_codiffusion_bioreactor.py, with the only changes being that the color map and
iterativeID levels are taken from in-memory new-data dicts instead of the old-data files.
"""
import matplotlib
matplotlib.use('Agg')  # never pop figures on screen

from json import load
from math import log
from numpy import log10, nanmean, where, isnan, inf, nan
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
from scipy.cluster import hierarchy
import sigfig

TAX = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# ================================================================= 1. load NEW data
_tr = read_csv('new_data/table_rel_codif_all.csv', index_col=0)          # seq x sample, PERCENT
_meta = pd.read_excel('new_data/16s_metadata_codif_all.xlsx')
_sample_days = {str(s): int(d) for s, d in zip(_meta['sample'], _meta['Day'])}
_sample_days = dict(sorted(_sample_days.items(), key=lambda kv: kv[1]))   # day-ordered, 29 samples
# clip timepoints: keep only Day 0 before Day 147 (drop intermediate early days 21-140)
_sample_days = {s: d for s, d in _sample_days.items() if d == 0 or d >= 147}
print(f'{_tr.shape[0]} ASVs x {len(_sample_days)} samples (clipped); days {list(_sample_days.values())}')

_t1 = read_csv('new_data/table1_codif_all.csv')
_taxa_df = _t1.drop_duplicates('seq').set_index('seq')[TAX]


def _clean(v):
    if isinstance(v, str) and v.strip() not in ('', 'NA', 'nan', 'None', 'NaN'):
        return v
    return None


taxonomy = {}
for _seq in _tr.index:                      # table_rel order drives iterativeID counters
    if _seq in _taxa_df.index:
        _row = _taxa_df.loc[_seq]
        taxonomy[_seq] = {l: _clean(_row[l]) for l in TAX}
    else:
        taxonomy[_seq] = {l: None for l in TAX}

# ================================================================= 2. iterativeIDs (notebook lines 203-220)
_iterativeIDs_inv, iterativeID_levels = ({}, {})
for _ASV, _taxa in taxonomy.items():
    best_level = 'Kingdom'
    for level, taxon in _taxa.items():
        if not isinstance(taxon, str) or level == 'Species':
            break
        best_level = level
    name = _taxa[best_level] if isinstance(_taxa[best_level], str) else 'Unknown'
    cand = name + '.1'
    while cand in _iterativeIDs_inv:
        base, cnt = cand.rsplit('.', 1)
        cand = f'{base}.{int(cnt) + 1}'
    _iterativeIDs_inv[cand] = _ASV
    iterativeID_levels[cand] = best_level
_iterativeIDs = {_v: _k for _k, _v in _iterativeIDs_inv.items()}          # seq -> iterativeID
_ID_LEVELS_NEW = iterativeID_levels
print(f'generated {len(_iterativeIDs)} iterativeIDs')

# ================================================================= 3. phylum color map (notebook lines 1810-1834)
_iid_phylum = {_iid: taxonomy[_seq].get('Phylum') for _seq, _iid in _iterativeIDs.items()}
_arch_markers = ['euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'thermoplasmatota',
                 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota']
archaea_phyla = sorted({p for p in _iid_phylum.values() if p and ('archaeo' in p.lower() or any(a in p.lower() for a in _arch_markers))})
bacteria_phyla = sorted({p for p in _iid_phylum.values() if p and p not in archaea_phyla})
_na, _nb = len(archaea_phyla), len(bacteria_phyla)
_taxa_color_map = {}
for _i, _p in enumerate(archaea_phyla):
    _taxa_color_map[_p] = plt.cm.turbo(_i / max(_na, 1) * 0.15)
for _i, _p in enumerate(bacteria_phyla):
    _taxa_color_map[_p] = plt.cm.turbo(0.2 + _i / max(_nb, 1) * 0.8)
iterativeID_color_map_1 = {_iid: _taxa_color_map[_p] for _iid, _p in _iid_phylum.items() if _p}
genera_color_map_1 = {_iid.split('.')[0]: _v for _iid, _v in iterativeID_color_map_1.items()}
print(f'{_na} archaea phyla, {_nb} bacteria phyla')

# ================================================================= 4. Shannon track (NEW data, normalized, 29 pts)
def shannon_index(abundances):
    _vals = [a for a in abundances if a > 0]
    _tot = sum(_vals)
    return -sum((a / _tot) * log(a / _tot) for a in _vals)


shannon_indices = {_sample_days[_s]: shannon_index(_tr[_s].values) for _s in _tr.columns if _s in _sample_days}


# ================================================================= taxonomy_linkage (verbatim)
def taxonomy_linkage(taxonomy_series):
    from collections import defaultdict
    import numpy as np
    from numpy import array
    parsed = taxonomy_series.fillna('unknown').astype(str).str.split('[;|,]', regex=True).apply(
        lambda lst: [x.strip() for x in lst or [] if x.strip()])
    ranks_n = max((len(_p) for _p in parsed))
    if ranks_n == 0:
        ranks_n = 1
    parsed = parsed.apply(lambda lst: (lst + [''] * ranks_n)[:ranks_n])
    _n = len(parsed)
    if _n <= 1:
        return np.empty((0, 4), dtype=float)
    idx_to_pos = {idx: pos for pos, idx in enumerate(parsed.index)}
    linkage_rows = []
    next_cluster_id = _n

    def build_subtree(indices, depth):
        nonlocal next_cluster_id
        if len(indices) == 1:
            return (idx_to_pos[indices[0]], 1)
        if depth >= ranks_n:
            merge_distance = 0.5
            cluster_id, total_count = (idx_to_pos[indices[0]], 1)
            for idx in indices[1:]:
                new_id = next_cluster_id
                next_cluster_id = next_cluster_id + 1
                linkage_rows.append([cluster_id, idx_to_pos[idx], merge_distance, total_count + 1])
                cluster_id = new_id
                total_count = total_count + 1
            return (cluster_id, total_count)
        groups = defaultdict(list)
        for idx in indices:
            rank_value = parsed.loc[idx][depth]
            if rank_value == '':
                rank_value = f'__unclassified_{idx}'
            groups[rank_value].append(idx)
        subclusters = []
        for rank_value, group_indices in groups.items():
            sub_id, sub_count = build_subtree(group_indices, depth + 1)
            subclusters.append((sub_id, sub_count))
        if len(subclusters) == 1:
            return subclusters[0]
        merge_distance = float(ranks_n - depth)
        cluster_id, total_count = subclusters[0]
        for sub_id, sub_count in subclusters[1:]:
            new_id = next_cluster_id
            next_cluster_id = next_cluster_id + 1
            linkage_rows.append([cluster_id, sub_id, merge_distance, total_count + sub_count])
            cluster_id = new_id
            total_count = total_count + sub_count
        return (cluster_id, total_count)
    build_subtree(list(parsed.index), depth=0)
    Z = array(linkage_rows, dtype=float)
    if len(Z) > 0:
        for _i in range(1, len(Z)):
            if Z[_i, 2] < Z[_i - 1, 2]:
                Z[_i, 2] = Z[_i - 1, 2]
    return Z


# ================================================================= _create_heatmap (verbatim except color/levels source)
def _create_heatmap(df, taxonomies, title, inlayed_data=None, linkage=None):
    new_cmap = LinearSegmentedColormap.from_list('NewMap', [(0.0, 'aliceblue'), (0.25, 'lightblue'), (1.0, 'navy')])
    new_cmap.set_bad('aliceblue')
    vmin = df.min().min()
    vmax = df.max().max()
    vcenter = log10(0.1)
    vcenter = min(max(vcenter, vmin + 1e-06), vmax - 1e-06)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    DEFAULT_COLOR = 'lightgray'

    def lookup(idx):
        if idx in iterativeID_color_map_1:
            return iterativeID_color_map_1[idx]
        elif idx in genera_color_map_1:
            return genera_color_map_1[idx]
        return DEFAULT_COLOR
    _row_colors = Series({idx: lookup(idx) for idx in df.index}, name='Phylum')
    _clusterMap = sns.clustermap(df, row_colors=_row_colors, cmap=new_cmap, norm=norm, clip_on=True, col_cluster=False, row_cluster=True, figsize=(20, 20), row_linkage=taxonomy_linkage(taxonomies), dendrogram_ratio=(0.2, 0.15))
    for tick in _clusterMap.ax_row_colors.get_xticklabels():
        tick.set_fontsize(22)
    _cbar = _clusterMap.ax_cbar
    log_ticks = np.array([vmin, log10(0.01), vcenter, log10(0.2), vmax])
    _cbar.set_yticks(log_ticks)
    original_ticks = 10 ** log_ticks
    original_labels = [f'{x * 100:.1e}'.replace('e-0', 'E-') if x * 100 < 0.1 else round(x * 100) for x in original_ticks]
    original_labels[0] = '0'
    _cbar.set_yticklabels(original_labels, fontsize=28)
    _cbar.set_yticks(log_ticks)
    _cbar.set_xlabel('')
    _cbar.set_ylabel('Rel. Abundance %', fontsize=32, labelpad=20, rotation=90)
    _cbar.tick_params(labelsize=28, length=10, width=2)
    if inlayed_data:
        top_ax = _clusterMap.ax_col_dendrogram
        top_ax.clear()
        _day_to_idx = {day: i + 0.5 for i, day in enumerate(df.columns)}
        _xs = [_day_to_idx[d] for d in inlayed_data if d in _day_to_idx]
        _ys = [inlayed_data[d] for d in inlayed_data if d in _day_to_idx]
        top_ax.plot(_xs, _ys, color='black', marker='o', linewidth=3)
        top_ax.set_xlim(0, len(df.columns))
        top_ax.set_ylabel('Shannon Diversity', fontsize=20)
        top_ax.yaxis.tick_right()
        top_ax.yaxis.set_label_position('right')
        top_ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        top_ax.set_xticks([])
        top_ax.tick_params(axis='y', labelsize=16)
    for _i in range(_clusterMap.data2d.shape[0]):
        for _j in range(_clusterMap.data2d.shape[1]):
            cell = _clusterMap.data2d.iloc[_i, _j]
            if isnan(cell):
                continue
            _value = 10 ** cell * 100
            if _value < 0.1:
                continue
            _color = 'black' if _value < 10 else 'white'
            if _value >= 1:
                str_val = str(sigfig.round(_value, 2))
            else:
                str_val = str(round(_value, 1))
            if str_val.endswith('.0'):
                str_val = str_val[:-2]
            _clusterMap.ax_heatmap.text(_j + 0.5, _i + 0.5, str_val, ha='center', va='center', color=_color, fontsize=20)
    _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=24, rotation=0)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=24, rotation=70, ha='right', rotation_mode='anchor')
    _clusterMap.ax_heatmap.set_xlabel('Days after inoculation', fontsize=32, labelpad=20)
    _clusterMap.ax_row_dendrogram.xaxis.set_visible(False)
    _clusterMap.ax_row_dendrogram.text(0.65, -0.02, 'Taxonomical tree', fontsize=24, ha='center', transform=_clusterMap.ax_row_dendrogram.transAxes)
    _clusterMap.ax_heatmap.yaxis.tick_left()
    _clusterMap.ax_heatmap.yaxis.set_label_position('left')
    label_right = 0.22
    heatmap_right = 0.68
    dendro_left = 0.72
    dendro_right = 0.85
    hm_pos = _clusterMap.ax_heatmap.get_position()
    dend_pos = _clusterMap.ax_row_dendrogram.get_position()
    _clusterMap.ax_heatmap.set_position([label_right, hm_pos.y0, heatmap_right - label_right, hm_pos.height])
    _clusterMap.ax_row_dendrogram.set_position([dendro_left, dend_pos.y0, dendro_right - dendro_left, dend_pos.height])
    _cbar_w = 0.025
    _cbar_h = hm_pos.height * 0.6
    _clusterMap.ax_cbar.set_position([
        dendro_right + 0.05,
        hm_pos.y0 + (hm_pos.height - _cbar_h) / 2,
        _cbar_w,
        _cbar_h,
    ])
    _clusterMap.ax_row_dendrogram.invert_xaxis()
    gap = 0.03
    col_dend_pos = _clusterMap.ax_col_dendrogram.get_position()
    _clusterMap.ax_col_dendrogram.set_position([label_right, col_dend_pos.y0 + gap, heatmap_right - label_right, col_dend_pos.height])
    _organisms_to_highlight = ['Methanobacterium', 'Methanosarcina', 'Methanobacteriaceae']
    _iterativeID_levels = _ID_LEVELS_NEW
    _ID_levels = {}
    for _k, _v in _iterativeID_levels.items():
        _ID_levels[_k] = _v
        _ID_levels.setdefault(_k.split('.')[0], _v)
    for _label in _clusterMap.ax_heatmap.get_yticklabels():
        _text = _label.get_text()
        if any((x in _text for x in _organisms_to_highlight)):
            _label.set_fontweight('bold')
        if _ID_levels.get(_text) == 'Genus':
            _label.set_fontstyle('italic')
    hm_pos = _clusterMap.ax_heatmap.get_position()
    fig_w = _clusterMap.figure.get_figwidth()
    fig_h = _clusterMap.figure.get_figheight()
    strip_w = 0.015
    strip_h = 0.015
    _clusterMap.ax_row_colors.set_position([hm_pos.x1 + 0.005, hm_pos.y0, strip_w, hm_pos.height])
    y_pad_pts = strip_w * fig_w * 72 + 12
    _clusterMap.ax_heatmap.tick_params(axis='y', pad=y_pad_pts)
    if isinstance(_row_colors, Series):
        rc = _row_colors
    else:
        rc = Series(list(_row_colors), index=df.index)
    phylum_color = {}
    for idx in df.index:
        parts = str(taxonomies.get(idx, '')).split('|')
        if len(parts) < 2:
            continue
        phylum = parts[1]
        if phylum in ('None', '', 'Unknown'):
            continue
        _color = rc.get(idx)
        if _color is None:
            continue
        phylum_color.setdefault(phylum, _color)
    archaea_markers = ('archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'thermoplasmatota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota')

    def is_archaea(p):
        return any((m in p.lower() for m in archaea_markers))
    archaea = sorted((_p for _p in phylum_color if is_archaea(_p)))
    bacteria = sorted((_p for _p in phylum_color if not is_archaea(_p)))
    handles = []
    if archaea:
        handles.append(Patch(color='none', label='$\\bf{Archaea}$'))
        handles = handles + [Patch(facecolor=phylum_color[_p], label=_p) for _p in archaea]
    if bacteria:
        handles.append(Patch(color='none', label='$\\bf{Bacteria}$'))
        handles = handles + [Patch(facecolor=phylum_color[_p], label=_p) for _p in bacteria]
    _clusterMap.figure.legend(handles=handles, title='Phylum', title_fontsize=20, fontsize=16, loc='upper right', bbox_to_anchor=(0.07, 1.0), frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6)
    for spine in _clusterMap.ax_heatmap.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(1)
    _clusterMap.figure.savefig(f"abundance_heatmaps/{title.lower().replace(' ', '_')}.png", bbox_inches='tight', dpi=800)
    _clusterMap.figure.savefig(f"abundance_heatmaps/{title.lower().replace(' ', '_')}.svg", bbox_inches='tight')
    print('  saved', f"abundance_heatmaps/{title.lower().replace(' ', '_')}.(png|svg)  [{df.shape[0]} rows x {df.shape[1]} cols]")
    plt.close(_clusterMap.figure)


def _change_columns(df):
    return df[[c for c in _sample_days.values() if c in df.columns]]


_frac = _tr / 100.0   # seq x sample fraction

# ================================================================= 5. Methanogens figure (Archaea, genus-aggregated)
print('\n[1/2] Methanogens Abundances (% abundance)')
_dic, taxonomies_m = ({}, {})
for _seq, _iid in _iterativeIDs.items():
    _t = taxonomy[_seq]
    if _t['Kingdom'] == 'Bacteria':
        continue
    _gid = _iid.split('.')[0]
    taxonomies_m.setdefault(_gid, '|'.join([str(_t[l]) for l in TAX[:TAX.index('Genus') + 1]]))
    _col = _frac.loc[_seq]
    for _s, _day in _sample_days.items():
        _v = float(_col[_s])
        if _v > 0:
            _dic.setdefault(_day, {}).setdefault(_gid, 0.0)
            _dic[_day][_gid] += _v
_new_dic = {_day: {g: log10(a) for g, a in c.items()} for _day, c in _dic.items()}
df_m = _change_columns(DataFrame(_new_dic)).astype(float).replace([inf, -inf], nan)
tax_m = Series({idx: taxonomies_m.get(idx, f'Unknown|{idx}') for idx in df_m.index})
_create_heatmap(df_m, tax_m, 'Methanogens Abundances (% abundance)', shannon_indices)

# ================================================================= 6. Top-10 ASV figure (per-ASV, top10/day union)
print('\n[2/2] Top 10 ASVs (% abundance)')
_nonzero_per_day, taxonomies_t = ({}, {})
for _seq, _iid in _iterativeIDs.items():
    _t = taxonomy[_seq]
    taxonomies_t.setdefault(_iid, '|'.join([str(_t[l]) for l in TAX[:TAX.index('Genus') + 1]]))
for _s, _day in _sample_days.items():
    _col = _frac[_s]
    _nz = {_iterativeIDs[_seq]: float(_v) for _seq, _v in _col.items() if _v > 0}
    _nonzero_per_day[_day] = dict(sorted(_nz.items(), key=lambda kv: kv[1], reverse=True))
_topNum = 10
_top10_union = set()
for _day, _d in _nonzero_per_day.items():
    _top10_union.update(list(_d.keys())[:_topNum])
_top_per_day = {_day: {_iid: log10(_v) for _iid, _v in _d.items() if _iid in _top10_union} for _day, _d in _nonzero_per_day.items()}
taxonomies_t = {_iid: _tx for _iid, _tx in taxonomies_t.items() if _iid in _top10_union}
df_t = _change_columns(DataFrame(_top_per_day)).astype(float).replace([inf, -inf], nan)
tax_t = Series({idx: taxonomies_t.get(idx, f'Unknown|{idx}') for idx in df_t.index})
_create_heatmap(df_t, tax_t, f'Top {_topNum} ASVs (% abundance)', shannon_indices)

print('\nDone.')
