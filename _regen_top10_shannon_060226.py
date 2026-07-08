"""
Regenerate the Top-10 ASV (and methanogens) % abundance heatmaps with the
Shannon diversity track computed from the NEW dataset
(new_data/table_rel_codif_all.csv), proportions normalized to sum to 1.

Faithful copy of `taxonomy_linkage` and `_create_heatmap` from
data_processing_codiffusion_bioreactor.py (cells at lines 336 and 588). The
heatmap cells remain the old-data top-10 ASV abundances; only the Shannon
Diversity track now follows the new-data values -- which notably raises the
inoculum (day 0, sample '5-10') to 5.583 (vs old-data 4.488, stored 3.633).
"""
import matplotlib
matplotlib.use('Agg')  # never pop figures on screen

from json import load, dump
from math import log
from numpy import log10, nanmean, where, isnan, inf, nan, logspace
import numpy as np
import seaborn as sns
from pandas import DataFrame, Series, read_csv
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
from scipy.cluster import hierarchy
import sigfig

_taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
_sample_days = {'5-10': 0, 'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168,
                'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283,
                'Q': 300, 'Qa': 315}


# ---------------------------------------------------------------- taxonomy_linkage (verbatim)
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


# ---------------------------------------------------------------- inputs (verbatim)
iterativeID_color_map_1 = load(open('iterativeID_color_map.json', 'r'))
genera_color_map_1 = {_ID.split('.')[0]: _v for _ID, _v in iterativeID_color_map_1.items()}
_total_df = read_csv('model_inputs/total.csv').set_index('seq')
_total_df['rel_ab'] = _total_df['rel_ab'] / 100


# ---------------------------------------------------------------- _create_heatmap (verbatim)
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
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
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
    archaea_markers = ('archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota')

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
    print('  saved', f"abundance_heatmaps/{title.lower().replace(' ', '_')}.png / .svg")


# ---------------------------------------------------------------- CORRECTED Shannon (normalized)
_iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))


def shannon_index(abundances):
    """H' = -sum(p ln p) with p normalized to sum to 1 (proper Shannon, natural log)."""
    _vals = [a for a in abundances if a > 0]
    _tot = sum(_vals)
    return -sum((a / _tot) * log(a / _tot) for a in _vals)


# Shannon track follows the NEW dataset (new_data/table_rel_codif_all.csv).
# shannon_index() normalizes internally, so passing percentage columns is fine.
_tr_new = read_csv('new_data/table_rel_codif_all.csv', index_col=0)
shannon_indices = {_sample_days[_s]: shannon_index(_tr_new[_s].values)
                   for _s in _tr_new.columns if _s in _sample_days}
print('New-data Shannon track (day -> H\'):')
for _d in sorted(shannon_indices):
    print(f'  day {_d:>3}: {shannon_indices[_d]:.4f}')
# keep the canonical record in sync with what the figure now plots
dump({str(_d): shannon_indices[_d] for _d in sorted(shannon_indices)},
     open('shannon_indices.json', 'w'), indent=2)


def _change_columns(df):
    return df[[c for c in _sample_days.values() if c in df.columns]]


# ---------------------------------------------------------------- methanogens figure (verbatim prep)
print('\n[1/2] Methanogens Abundances (% abundance)')
_dic, taxonomies_1 = ({}, {})
level_2 = 'Genus'
for _seq, _row in _total_df.iterrows():
    _day = _sample_days.get(_row['sample'])
    if _row['Kingdom'] == 'Bacteria' or _day is None:
        continue
    _abundance = _row['rel_ab']
    if _abundance == 0:
        continue
    _iterativeID = _iterativeIDs.get(_seq) or 'Archaea.99'
    if isinstance(_iterativeID, str):
        _iterativeID = _iterativeID.split('.')[0]
    taxonomies_1.setdefault(_iterativeID, '|'.join([str(_row[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index(level_2)]))
    _dic.setdefault(_day, {})
    _dic[_day].setdefault(_iterativeID, 0)
    _dic[_day][_iterativeID] = _dic[_day][_iterativeID] + _abundance
new_dic = {_day: {genus: log10(abund) for genus, abund in _content.items()} for _day, _content in _dic.items()}
df_2 = _change_columns(DataFrame(new_dic)).astype(float).replace([inf, -inf], nan)
taxonomy_series_1 = Series({idx: taxonomies_1.get(idx, f'Unknown|{idx}') for idx in df_2.index})
_create_heatmap(df_2, taxonomy_series_1, 'Methanogens Abundances (% abundance)', shannon_indices)


# ---------------------------------------------------------------- top-10 ASV figure (verbatim prep)
print('\n[2/2] Top 10 ASVs (% abundance)')
_dic, taxonomies_1 = ({}, {})
for _seq, _row in _total_df.iterrows():
    _day = _sample_days.get(_row['sample'])
    if _day is None:
        continue
    _iterativeID = _iterativeIDs.get(_seq)
    taxonomies_1.setdefault(_iterativeID, '|'.join([str(_row[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index('Genus')]))
    _dic.setdefault(_day, {})
    _dic[_day].setdefault(_iterativeID, 0)
    _dic[_day][_iterativeID] = _dic[_day][_iterativeID] + _row['rel_ab']
_nonzero_per_day = {_day: dict(sorted({_k: _v for _k, _v in _org_dict.items() if _v > 0}.items(), key=lambda item: item[1], reverse=True)) for _day, _org_dict in _dic.items()}
_topNum = 10
top10_all_days = []
for _day, _org_dict in _nonzero_per_day.items():
    top10_all_days.extend(list(_org_dict.keys())[:_topNum])
top10_all_days = list(set(top10_all_days))
_top_per_day = {}
for _day, _org_dict in _nonzero_per_day.items():
    _top_per_day[_day] = {_org: log10(_v) for _org, _v in _org_dict.items() if _org in top10_all_days}
taxonomies_1 = {_org: _taxa for _org, _taxa in taxonomies_1.items() if _org in top10_all_days}
df_2 = _change_columns(DataFrame(_top_per_day)).astype(float).replace([inf, -inf], nan)
taxonomy_series_1 = Series({idx: taxonomies_1.get(idx, f'Unknown|{idx}') for idx in df_2.index})
_create_heatmap(df_2, taxonomy_series_1, f'Top {_topNum} ASVs (% abundance)', shannon_indices)

print('\nDone.')
