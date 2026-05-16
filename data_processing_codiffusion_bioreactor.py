import marimo

__generated_with = "0.23.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # load files
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Measurement data
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### summarizing all of the non-sequencing data in `measurements/*.csv`
    - VFA.csv: Contains Acetate, Propionate, Butyrate, and COD for the reactor, waste, and media
    - Ammonia.csv: Ammonia concentrations in the reactor, waste, and media
    - TSS_VSS.csv: Abundances of TSS & VSS in the system
    - Pre-Summary.csv: Defining sampling codes as days of operation
    - GC_TCD.csv: gaseous components of various reactors
    - Summary.csv: `main content`.  Complete gaseous ins/outs, carbon accounting, occasional media VFA and TSS measurements
    """)
    return


@app.cell
def _(sheetname):
    from pandas import read_excel

    measurements = read_excel("model_inputs/measurements.xlsx", engine='openpyxl', sheet_name=None)
    print(len(measurements))
    for sheetName, csv in measurements.items():
        csv.to_csv(f"model_inputs/{sheetname}.csv")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Process Summary.csv
    """)
    return


@app.cell
def _(display):
    from pandas import api, read_csv, set_option
    from json import load, dump
    summary = read_csv('model_inputs/measurements/Summary.csv', header=1).set_index('Days of Operation').drop(['0', '2024-05-14 12:42:00'], axis=1)
    for _col in summary.columns:
        if not api.types.is_numeric_dtype(summary[_col]):
            continue
        series = summary[_col]  # Skip non-numeric columns
        if not series.isna().any():
            continue
        if series.isna().all():
            continue
        summary[_col] = series.interpolate(method='linear', limit_direction='forward', limit_area='inside')  # Skip if no missing values
    summary.iloc[0:4, 5:9] = summary.iloc[5, 5:9].values
    summary.loc[300, 'Biomass Sample ID'] = 'Q'
    _summary_samples = summary[summary['Biomass Sample ID'].notna()]  # Skip if all values are missing (can't interpolate)
    _summary_samples.set_index('Biomass Sample ID', inplace=True)
    summary.to_csv('model_inputs/measurements/Summary_interpolated.csv')
    _dic = _summary_samples.T.to_dict()
    display(list(_dic.items())[0])
    # display(summary_samples) #summary.loc["300"])
    dump(_dic, open('model_inputs/measurements/summary_samples.json', 'w'), indent=3)  # Only fill NaN values inside valid values (not at edges)
    return dump, load, read_csv, set_option


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### ModelSEED compounds for these column names
    """)
    return


@app.cell
def _():
    summary_key_to_modelseed = {

        # ---- MEDIA COMPOUNDS ----
        "Media Acetate (mol/L)": "cpd00029",        # acetate
        "Media Propionate (mol/L)": "cpd00141",     # propionate
        "Media Butyrate (mol/L)": "cpd00211",       # butyrate
        "Media Ammonia (mol/L)": "cpd00013",        # NH3/NH4+
        "Media Dissolved Inorganic Carbon (mol C/L)": "cpd00011",  # CO2 (aq)

        # ---- GAS FEED ----
        "H2 delivery rate (mol/min)": "cpd00067",   # hydrogen
        "CO2 delivery rate (mol/min)": "cpd00011",  # carbon dioxide

        # ---- WASTE / EFFLUENT COMPOUNDS ----
        "Waste Effluent Acetate (mol/L)": "cpd00029",
        "Waste Effluent Propionate (mol/L)": "cpd00141",
        "Waste Effluent Butyrate (mol/L)": "cpd00211",
        "Waste Effluent Ammonia (mol/L)": "cpd00013",
        "Waste Effluent Dissolved Inorganic Carbon (mol C/L)": "cpd00011",
        "Waste Effluent Dissolved Organic Carbon (mol C/L)": "cpd00027",  # generic glucose-equivalent DOC proxy

        # ---- GAS PHASE ----
        "Gas Composition (% CH4)": "cpd00104",      # methane
        "Gas Composition (% H2)": "cpd00067",
        "Gas Composition (% CO2)": "cpd00011",
        "CH4 Production (mol/min)": "cpd00104",
        "H2 breakthrough (mol/min)": "cpd00067",
        "CO2 breakthrough (mol/min)": "cpd00011",

        # ---- AQUEOUS CARBONATE SPECIATION ----
        "[CO2aq]": "cpd00011",
        "[HCO3-]": "cpd00099",   # bicarbonate
        "[CO3-2]": "cpd00060",  # carbonate

    }
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Taxonomy and Abundance data
    """)
    return


@app.cell
def _(display, dump, read_csv):
    ab = read_csv('model_inputs/abundances.csv').set_index('seq')
    meaningful_ASVs = ab[(ab > 0.001).any(axis=1)]
    print(len(ab) - len(meaningful_ASVs), 'negligible abundance ASVs')
    display(meaningful_ASVs)
    ab_js = meaningful_ASVs.to_dict()
    dump(ab_js, open('model_inputs/abundances.json', 'w'), indent=2)
    df = read_csv('model_inputs/taxonomy.csv').set_index('seq').drop('sample', axis=1)
    df = df.where(df.notna(), None)
    taxonomy_js = {}
    for _i, _row in df.iterrows():
        if _i not in meaningful_ASVs.index:
            continue
        taxonomy_js.setdefault(_i, {})
        taxonomy_js[_i]['Kingdom'] = _row['Kingdom']
        taxonomy_js[_i]['Phylum'] = _row['Phylum']
        taxonomy_js[_i]['Class'] = _row['Class']
        taxonomy_js[_i]['Order'] = _row['Order']
        taxonomy_js[_i]['Family'] = _row['Family']
        taxonomy_js[_i]['Genus'] = _row['Genus']
        taxonomy_js[_i]['Species'] = _row['Species']
    display(list(taxonomy_js.items())[:3])
    dump(taxonomy_js, open('model_inputs/taxonomy.json', 'w'), allow_nan=False, indent=4)
    return ab_js, taxonomy_js


@app.cell
def _():
    # from json import load, dump
    # IDmd5 = load(open("model_inputs/ID_md5.json", 'r'))
    # dump({v: k.split()[1].split(".rna")[0] for k,v in IDmd5.items()}, open("model_inputs/md5_ID.json", 'w'), indent=2)
    # list(IDmd5.items())[:2]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Taxonomy
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Create iterative IDs
    """)
    return


@app.cell
def _(display, dump, load):
    from collections import Counter
    midas_taxonomy = load(open('model_inputs/taxonomy.json', 'r'))
    genera_counts = Counter([x['Genus'] for x in midas_taxonomy.values()])
    display({_k: _v for _k, _v in genera_counts.most_common() if _v > 2})
    _iterativeIDs, _iterativeID_levels = ({}, {})
    for _ASV, _taxa in midas_taxonomy.items():
        best_level = 'Kingdom'
        for level, taxon in _taxa.items():
            if not isinstance(taxon, str) or level == 'Species':
                break
            best_level = level
        taxon = _taxa[best_level]
        taxon = taxon + '.1'
        while taxon in _iterativeIDs:
            taxon, _count = taxon.split('.')
            taxon = f'{taxon}.{int(_count) + 1}'
        _iterativeIDs[taxon] = _ASV
        _iterativeID_levels[taxon] = best_level
    _iterativeIDs = {_v: _k for _k, _v in _iterativeIDs.items()}
    display(list(_iterativeIDs.items())[:5])
    dump(_iterativeIDs, open('model_inputs/iterativeIDs.json', 'w'))
    dump(_iterativeID_levels, open('model_inputs/iterativeID_levels.json', 'w'))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## taxonomic representation
    """)
    return


@app.cell
def _(ab_js, taxonomy_js):
    import matplotlib.pyplot as plt
    from numpy import nan

    def taxa_bar_chart(taxa_asvs, x_label, title, export_name):
        limit_size = 40  # Create the bar chart
        sorted_data = dict(sorted({_k: len(x) for _k, x in taxa_asvs.items()}.items(), key=lambda x: x[1], reverse=True))
        sorted_data = dict(list(sorted_data.items())[:limit_size])
        plt.bar(list(sorted_data.keys()), list(sorted_data.values()), color='lightblue', edgecolor='black')
        for _i, _value in enumerate(list(sorted_data.values())):
            plt.text(_i, _value + max(list(sorted_data.values())) / 30, str(_value), ha='center', fontsize=12, color='black', rotation=90)
        plt.xlabel(x_label)  # Add labels on top of each bar
        plt.ylabel(f'Frequency ({sum(list(sorted_data.values()))} total)')
        plt.title(title)
        plt.xticks(rotation=90)
        plt.tight_layout()  # Add labels and title
        plt.savefig(export_name)

    def define_taxa_asvs(taxa, full_taxonomy):
        taxa_asvs = {}
        for asv, taxonomy in full_taxonomy.items():  # Rotate x-axis labels if needed
            if str(taxonomy.get(_taxa, 'none')) in [nan, 'nan', 'NaN', 'none', 'None'] or 'midas' in taxonomy.get(_taxa, 'none'):
                continue
            if taxonomy.get(_taxa, 'none') not in taxa_asvs:  # Show the plot
                taxa_asvs[taxonomy.get(_taxa, 'none')] = []
            taxa_asvs[taxonomy.get(_taxa, 'none')].append(asv)
        return taxa_asvs  # plt.show()
    for _taxa in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
        taxa_asvs = define_taxa_asvs(_taxa, taxonomy_js)
        taxa_bar_chart(taxa_asvs, _taxa, f'{_taxa} ({len(taxa_asvs)} total) ASVs', f'taxonomic_representations/{_taxa}_asvs.png')
    for _taxa in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
        for _sample, _content in ab_js.items():
            full_taxonomy = {_k: taxonomy_js[_k] for _k, _v in _content.items() if _v > 0}
            taxa_asvs = define_taxa_asvs(_taxa, full_taxonomy)
    # from deepdiff import DeepDiff
            taxa_bar_chart(taxa_asvs, _taxa, f'{_sample} {_taxa} ({len(taxa_asvs)} total) ASVs', f'taxonomic_representations/{_sample}_{_taxa}_asvs.png')  # display(DeepDiff(taxonomy_js, full_taxonomy))  # break
    return nan, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Abundance
    """)
    return


@app.cell
def _(ab_js, plt):
    from numpy import array
    for _sample, _content in ab_js.items():
        print(_sample)
        _abundances = array([x for x in list(_content.values()) if x > 0])
        plt.hist(_abundances, bins=200)
        plt.xlabel('% Relative Abundance')
        plt.ylabel(f'ASVs ({len(_abundances)} total)')
        plt.yscale('log')
        plt.title(f'{_sample} ASV Abundance')
        plt.savefig(f'abundance_distributions/{_sample}_abundance.png')
        plt.show()
    return (array,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Shannon Diversity Index
    """)
    return


@app.cell
def _(array, load, plt):
    from numpy import log
    _asvSet_abundances = load(open(f'modeling_files/ASVset_abundances.json', 'r'))
    SDI = {}
    for _sample, asvAbund in _asvSet_abundances.items():
        _abundances = array(list(asvAbund.values()))
        _abundances = _abundances / sum(_abundances)
        SDI[_sample] = -sum(log(_abundances) * _abundances)
    plt.figure(figsize=(10, 6))
    plt.plot(SDI.keys(), SDI.values(), marker='o', linestyle='-', color='blue', linewidth=2, markersize=8)
    plt.xlabel('Samples')
    plt.ylabel('Shannon Index', fontsize=20)
    plt.title('Shannon Diversity Index for all samples', fontsize=20)
    plt.grid(True)
    plt.yticks(fontsize=20)
    plt.xticks(rotation=60, fontsize=20)
    plt.tight_layout()
    plt.show()
    return (log,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Abundance heatmaps
    """)
    return


@app.cell
def _():
    def taxonomy_linkage(taxonomy_series):
        """
        Build a linkage matrix that EXACTLY follows taxonomy hierarchy.
        Auto-detects depth from the taxonomy strings.
        """
        from collections import defaultdict
        import numpy as np
        from numpy import array
        parsed = taxonomy_series.fillna('unknown').astype(str).str.split('[;|,]', regex=True).apply(lambda lst: [x.strip() for x in lst or [] if x.strip()])
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

    return (taxonomy_linkage,)


@app.cell
def _(display, nan, read_csv, taxonomy_linkage):
    from pandas import DataFrame, Series, concat
    import seaborn as sns
    from matplotlib import pyplot
    from numpy import inf, empty, logspace
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage
    from collections import defaultdict

    def concatenate_taxonomy(df, level):
        _dic = {}
        for index, _row in df.iterrows():
            _dic['|'.join([str(_row[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index(level)])] = index
        return Series(_dic)

    def _change_columns(df):
        return df[[c for c in _sample_days.values() if c in df.columns]]
    from matplotlib.colors import LinearSegmentedColormap
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    _sample_days = {'5-10': 0, 'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}
    _total_df = read_csv('model_inputs/total.csv').set_index('seq')
    _max_abundance = _total_df['rel_ab'].max()
    print('max', _max_abundance)
    _total_df['rel_ab'] = _total_df['rel_ab'] / 100
    _max_abundance = _total_df['rel_ab'].max()
    print('max', _max_abundance)
    _intervals = [0] + list(logspace(-2, -0.3, 5))
    print(_intervals)

    def _create_heatmap(df, level, abundance_range, taxonomies, cmap=None, log=False):
        cmap = cmap or LinearSegmentedColormap.from_list('NewMap', ['skyblue', 'red'])
        _clusterMap = sns.clustermap(df, cmap=cmap, vmin=0, vmax=max(list(df.to_numpy().flatten())), col_cluster=False, figsize=(20, 20), row_linkage=taxonomy_linkage(Series(taxonomies)), dendrogram_ratio=(0.1, 0.2))
        _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
        _clusterMap.figure.suptitle(f'{level.capitalize()} {abundance_range} Abundances', fontsize=30)
        _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=20, rotation=0)
        _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=24, rotation=70, ha='right', rotation_mode='anchor')
        _clusterMap.ax_heatmap.set_xlabel('Days after inoculation', fontsize=30)
        _clusterMap.figure.tight_layout()
        _clusterMap.figure.savefig(f'abundance_heatmaps/{level}_{abundance_range}_abundances.png')
    for _i, abundance_limit in enumerate(_intervals):
        if abundance_limit == 0:
            continue
        previous_limit = _intervals[_i - 1]
        print(previous_limit, abundance_limit)
        abundance_range = f'{round(previous_limit, 2)}-{round(abundance_limit, 2)}'
        for level_1 in _taxonomic_levels:
            _dic = {}
            taxonomies = {}
            for _seq, _row in _total_df.iterrows():
                _day = _sample_days.get(_row['sample'])
                if _day is None:
                    continue
                _abundance = _row['rel_ab']
                if _abundance <= previous_limit or _abundance > abundance_limit:
                    continue
                taxonomies.setdefault(_row[level_1], '|'.join([str(_row[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index(level_1)]))
                _dic.setdefault(_day, {})
                _dic[_day].setdefault(_row[level_1], 0)
                _dic[_day][_row[level_1]] = _dic[_day][_row[level_1]] + _abundance
            if _dic == {}:
                print(f'The {level_1} does not have organisms in {abundance_range}')
                continue
            elif len(taxonomies) == 1:
                print(f'The {level_1} only has one organism in {abundance_range}:  {_dic} and {taxonomies}')
                continue
            if level_1 == 'Phylum':
                display(taxonomies)
            df_1 = _change_columns(DataFrame(_dic))
            df_1 = df_1.astype(float).replace([inf, -inf], nan).fillna(0)
            taxonomy_series = Series({idx: taxonomies.get(idx, f'Unknown|{idx}') for idx in df_1.index})
            print(f"Detected depth: {max((len(t.split('|')) for t in taxonomy_series))}")
            print(f'Sample entries:\n{taxonomy_series.head()}')
            _create_heatmap(df_1, level_1, abundance_range, taxonomy_series)
    return (
        DataFrame,
        LinearSegmentedColormap,
        Series,
        defaultdict,
        df_1,
        inf,
        logspace,
        sns,
        taxonomies,
    )


@app.cell
def _(taxonomies):
    taxonomies
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### Making Heather's requested figures
    """)
    return


@app.cell
def _(
    G,
    archaea_phyla,
    bacteria_phyla,
    df_1,
    display,
    dump,
    iterativeID_level,
    load,
    mpatches,
    plt,
    taxonomy,
):
    from pathlib import Path

    def rgba_to_hex(rgba):
        return '#{:02x}{:02x}{:02x}'.format(*(int(c * 255) for c in rgba[:3]))

    def build_phylum_colors(out_path='model_inputs/phylum_colors.json'):
        phyla = sorted({t.get('Phylum') for t in taxonomy.values() if t.get('Phylum')})
        archaea_markers = ['archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota']
        archaea = [_p for _p in phyla if any((m in _p.lower() for m in archaea_markers))]
        bacteria = [_p for _p in phyla if _p not in archaea]
        colors = {}
        for _i, _p in enumerate(archaea):
            colors[_p] = rgba_to_hex(plt.cm.turbo(_i / max(len(archaea), 1) * 0.15))
        for _i, _p in enumerate(bacteria):
            colors[_p] = rgba_to_hex(plt.cm.turbo(0.2 + _i / max(len(bacteria), 1) * 0.8))
        dump(colors, open(out_path, 'w'), indent=2)
        return colors
    build_phylum_colors()
    phylum_colors = load(open('model_inputs/phylum_colors.json', 'r'))
    iterativeID_color_map = load(open(f'iterativeID_color_map.json', 'r'))
    display(iterativeID_color_map)
    DEFAULT = 'lightgray'
    _node_colors = [phylum_colors.get(iterativeID_level.get(_n, 'Unknown'), DEFAULT) for _n in G.nodes()]
    print(df_1.index)
    if '.' in df_1.index[0]:
        _row_colors = [iterativeID_color_map[_org] for _org in df_1.index]
    else:
        genera_color_map = {_ID.split('.')[0]: _v for _ID, _v in iterativeID_color_map.items()}
        _row_colors = [genera_color_map.get(_org, DEFAULT) for _org in df_1.index]
    print(len(_row_colors), _row_colors)
    _archaea_patches = [mpatches.Patch(color=phylum_colors[_p], label=_p) for _p in archaea_phyla]
    _bacteria_patches = [mpatches.Patch(color=phylum_colors[_p], label=_p) for _p in bacteria_phyla]
    return


@app.cell
def _(
    DataFrame,
    LinearSegmentedColormap,
    Series,
    display,
    dump,
    inf,
    load,
    log,
    logspace,
    nan,
    np,
    read_csv,
    sns,
    taxonomy_linkage,
):
    from matplotlib.colors import LogNorm, Normalize, BoundaryNorm, TwoSlopeNorm
    from matplotlib.patches import Patch
    from matplotlib.transforms import Bbox
    from scipy.cluster import hierarchy
    import sigfig
    from numpy import log10, delete, nanmean, where, isnan
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    _sample_days = {'5-10': 0, 'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}
    iterativeID_color_map_1 = load(open(f'iterativeID_color_map.json', 'r'))
    genera_color_map_1 = {_ID.split('.')[0]: _v for _ID, _v in iterativeID_color_map_1.items()}
    _total_df = read_csv('model_inputs/total.csv').set_index('seq')
    _max_abundance = _total_df['rel_ab'].max()
    print('max', _max_abundance)
    _total_df['rel_ab'] = _total_df['rel_ab'] / 100
    _max_abundance = _total_df['rel_ab'].max()
    print('max', _max_abundance)
    _intervals = [0] + list(logspace(-2, -0.3, 5))
    print(_intervals)

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
        print('phylum_color', phylum_color)
        archaea_markers = ('archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota')

        def is_archaea(p):
            return any((m in _p.lower() for m in archaea_markers))
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
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))

    def shannon_index(abundances):
        return -sum([_abundance * log(_abundance) for _abundance in _abundances if _abundance > 0])
    abundances_js = load(open('model_inputs/abundances.json', 'r'))
    shannon_indices = {_sample_days[_s]: shannon_index(list(_abundances.values())) for _s, _abundances in abundances_js.items() if _s in _sample_days}
    display(shannon_indices)
    _dic = {}
    taxonomies_1 = {}
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
    new_dic = {}
    for _day, _content in _dic.items():
        new_dic[_day] = {genus: log10(abund) for genus, abund in _content.items()}

    def _change_columns(df):
        return df[[c for c in _sample_days.values() if c in df.columns]]
    df_2 = _change_columns(DataFrame(new_dic))
    display(df_2)
    df_2 = df_2.astype(float).replace([inf, -inf], nan)
    taxonomy_series_1 = Series({idx: taxonomies_1.get(idx, f'Unknown|{idx}') for idx in df_2.index})
    _arr = df_2.values.astype(float)
    row_means = nanmean(_arr, axis=1, keepdims=True)
    arr_filled = where(isnan(_arr), row_means, _arr)
    row_linkage = hierarchy.linkage(arr_filled, method='average', metric='euclidean')
    _create_heatmap(df_2, taxonomy_series_1, 'Methanogens Abundances (% abundance)', shannon_indices, row_linkage)
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
    dump(_nonzero_per_day, open('model_inputs/nonzero_per_day.json', 'w'), indent=2)
    _top_per_day = {}
    _all_orgs = {}
    _topNum = 10
    top10_all_days = []
    for _day, _org_dict in _nonzero_per_day.items():
        top10_all_days.extend(list(_org_dict.keys())[:_topNum])
    top10_all_days = list(set(top10_all_days))
    for _day, _org_dict in _nonzero_per_day.items():
        _top_per_day.setdefault(_day, {})
        _top_per_day[_day] = {_org: log10(_v) for _org, _v in _org_dict.items() if _org in top10_all_days}
    taxonomies_1 = {_org: _taxa for _org, _taxa in taxonomies_1.items() if _org in top10_all_days}
    display(_top_per_day)
    display(taxonomies_1)
    df_2 = _change_columns(DataFrame(_top_per_day))
    df_2 = df_2.astype(float).replace([inf, -inf], nan)
    taxonomy_series_1 = Series({idx: taxonomies_1.get(idx, f'Unknown|{idx}') for idx in df_2.index})
    display(df_2)
    _arr = df_2.values.astype(float)
    row_means = nanmean(_arr, axis=1, keepdims=True)
    arr_filled = where(isnan(_arr), row_means, _arr)
    row_linkage = hierarchy.linkage(arr_filled, method='average', metric='euclidean')
    _create_heatmap(df_2, taxonomy_series_1, f'Top {_topNum} ASVs (% abundance)', shannon_indices, row_linkage)
    return genera_color_map_1, iterativeID_color_map_1


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### table of iterative IDs
    """)
    return


@app.cell
def _(DataFrame, load):
    all_top10 = {'Methanobacterium.1': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobacterium', 'Methanobacterium.2': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobacterium', 'Methanobacteriaceae.1': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|nan', 'Mesotoga.1': 'Bacteria|Thermotogota|Thermotogae|Kosmotogales|Kosmotogaceae|Mesotoga', 'Lentimicrobium.1': 'Bacteria|Bacteroidota|Bacteroidia|Sphingobacteriales|Lentimicrobiaceae|Lentimicrobium', 'Aminivibrio.1': 'Bacteria|Synergistota|Synergistia|Synergistales|Synergistaceae|Aminivibrio', 'midas_g_94288.1': 'Bacteria|Firmicutes|Clostridia|Gracilibacteraceae|Lutispora|midas_g_94288', 'Desulfovibrio.1': 'Bacteria|Desulfobacterota|Desulfovibrionia|Desulfovibrionales|Desulfovibrionaceae|Desulfovibrio', 'Methanobacterium.3': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobacterium', 'Petrimonas.1': 'Bacteria|Bacteroidota|Bacteroidia|Bacteroidales|Dysgonomonadaceae|Petrimonas', 'Burkholderiales.1': 'Bacteria|Proteobacteria|Gammaproteobacteria|Burkholderiales|nan|nan', 'midas_g_9269.1': 'Bacteria|Bacteroidota|Bacteroidia|Bacteroidales|Paludibacteraceae|midas_g_9269', 'Methanobacterium.5': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobacterium', 'Proteiniphilum.1': 'Bacteria|Bacteroidota|Bacteroidia|Bacteroidales|Dysgonomonadaceae|Proteiniphilum', 'midas_g_36215.1': 'Bacteria|Bacteroidota|Bacteroidia|Bacteroidales|Rikenellaceae|midas_g_36215', 'Methanobacterium.6': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobacterium', 'Enterobacteriaceae.1': 'Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|nan', 'EBM-39.1': 'Bacteria|Synergistota|Synergistia|Synergistales|Synergistaceae|EBM-39', 'Dysgonomonas.1': 'Bacteria|Bacteroidota|Bacteroidia|Bacteroidales|Dysgonomonadaceae|Dysgonomonas', 'Methanobacterium.78': 'Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobacterium', 'midas_g_25828.1': 'Bacteria|Bacteroidota|Bacteroidia|Bacteroidales|Prolixibacteraceae|midas_g_25828', 'Azoarcus.4': 'Bacteria|Proteobacteria|Gammaproteobacteria|Burkholderiales|Rhodocyclaceae|Azoarcus', 'Desulfovibrio.3': 'Bacteria|Desulfobacterota|Desulfovibrionia|Desulfovibrionales|Desulfovibrionaceae|Desulfovibrio', 'Desulfomicrobium.1': 'Bacteria|Desulfobacterota|Desulfovibrionia|Desulfovibrionales|Desulfomicrobiaceae|Desulfomicrobium', 'Stappia.1': 'Bacteria|Proteobacteria|Alphaproteobacteria|Rhizobiales|Stappiaceae|Stappia', 'midas_g_1799.7': 'Bacteria|Firmicutes|Clostridia|Proteinivoracales|midas_f_1799|midas_g_1799', 'midas_g_7.1': 'Bacteria|Firmicutes|Dethiobacteria|Dethiobacterales|Dethiobacteraceae|midas_g_7', 'Dethiobacteraceae.4': 'Bacteria|Firmicutes|Dethiobacteria|Dethiobacterales|Dethiobacteraceae|nan', 'Soehngenia.2': 'Bacteria|Firmicutes|Clostridia|Peptostreptococcales-Tissierellales|Family_XI|Soehngenia'}
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))
    IDasvs = {}
    for _k, _v in _iterativeIDs.items():
        IDasvs.setdefault(_v, []).append(_k)
    total_dic = {'ASV': [], 'Iterative ID': [], 'MiDAS Taxonomy': []}
    for _iterativeID, taxonomy in all_top10.items():
        total_dic['Iterative ID'].append(_iterativeID)
        total_dic['MiDAS Taxonomy'].append(taxonomy.replace('|', ' '))
        total_dic['ASV'].append(IDasvs.get(_iterativeID)[0])
    df_3 = DataFrame(total_dic).set_index('ASV')
    df_3.to_csv('model_inputs/asv_IDs.csv')
    df_3
    return (taxonomy,)


@app.cell
def _():
    # from pandas import read_csv, DataFrame, Series, concat
    # import seaborn as sns
    # from matplotlib import pyplot
    # from json import load
    # from numpy import inf, nan, empty, logspace
    # from scipy.spatial.distance import squareform
    # from scipy.cluster.hierarchy import linkage
    # from collections import Counter

    # # taxonomy = load(open("/Users/andrewfreiburger/Documents/Research/MicrobiomeNotebooks/digestor/new_abundances.json", 'r'))

    # def taxonomy_linkage(taxonomy_series, ranks=('kingdom','phylum','class','order','family','genus','species'),
    #                      method: str = 'average'):
    #     """
    #     Build a SciPy linkage matrix from a taxonomy Series by defining
    #     a simple distance: (#ranks) - (#equal ranks at the same positions).

    #     Parameters
    #     ----------
    #     taxonomy_series : pd.Series
    #         Index order must match df.index (rows). Values are taxonomy strings
    #         like "k__Bacteria;p__Firmicutes;c__Bacilli;...".
    #     ranks : tuple or int
    #         Number of ranks considered. If tuple, its length is used.
    #     method : str
    #         Linkage method for scipy.cluster.hierarchy.linkage.

    #     Returns
    #     -------
    #     Z : np.ndarray
    #         Linkage matrix suitable for seaborn.clustermap(row_linkage=Z).
    #     """
    #     # Normalize rank count
    #     ranks_n = len(ranks) if not isinstance(ranks, int) else int(ranks)

    #     # Split and standardize to fixed-length lists (pad/truncate)
    #     parts = (
    #         taxonomy_series.fillna("")
    #         .astype(str)
    #         .str.split(r"[;|,]", regex=True)
    #         .apply(lambda lst: [x.strip() for x in (lst or [])])
    #         .apply(lambda lst: (lst + [""] * ranks_n)[:ranks_n])
    #     )

    #     n = len(parts)
    #     if n <= 1:
    #         # Nothing to cluster; return empty linkage
    #         return empty((0, 4), dtype=float)

    #     # Build condensed distance vector: distance = ranks_n - (#equal positions)
    #     D_condensed = []
    #     for i in range(n - 1):
    #         ai = parts.iloc[i]  # POSitional access fixes KeyError
    #         for j in range(i + 1, n):
    #             aj = parts.iloc[j]
    #             same = sum(a == b and a != "" and b != "" for a, b in zip(ai, aj))
    #             D_condensed.append(ranks_n - same)

    #     # Convert to linkage (hierarchical clustering tree)
    #     Z = linkage(D_condensed, method=method)
    #     return Z


    # def concatenate_taxonomy(df, level):
    #     dic = {}
    #     for index, row in df.iterrows():
    #         dic["|".join([str(row[l]) for l in taxonomic_levels
    #                     if taxonomic_levels.index(l) <= taxonomic_levels.index(level)])] = index
    #     return Series(dic)

    # def change_columns(df):
    #     return df[[c for c in sample_days.values() if c in df.columns]]
    #     # return df[[c for c in sample_days if c in df.columns]]


    # from matplotlib.colors import LinearSegmentedColormap

    # # Define a single-color gradient colormap (white → blue)
    # new_cmap = LinearSegmentedColormap.from_list("NewMap", [(0, "skyblue"), (1, "red")])


    # from pandas import read_excel
    # from json import dump
    # # ab = read_csv("GAME Sequencing/ASV_abundance_wideform.csv").set_index("seq")
    # taxonomic_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    # # sample_order = ["10AB", "A34", "B34", "C3", "D34", "E34", "F34", "G12", "G3", "H34", "I34",
    # #                 "J12", #"J34",
    # #                 "K12", "L12", "M12", "N12", "P12"]
    # sample_days = {
    #     "5-10": 0,
    #     # "H": 140,
    #     "I": 149,
    #     "J": 161,
    #     "K": 182,
    #     "L": 189,
    #     "M": 210,
    #     "N": 238,
    #     "P": 283,
    #     "Q": 300,
    #     # additional new time points (no operational data)
    #     # "Ha": 147,
    #     # "Ia": 153,
    #     # "Ja": 168,
    #     # "Jb": 174,
    #     # "Na": 262,
    #     # "Qa": 315,
    # }

    # og = read_csv("model_inputs/").set_index("seq")
    # sample_dfs = []
    # for sample in sample_days:
    #     # print(sample)
    #     sample_og = og[og["sample"] == sample]
    #     # display(sample_og)
    #     sample_og["rel_ab"] = sample_og["rel_ab"] / sample_og["rel_ab"].sum()
    #     sample_dfs.append(sample_og)
    # og = concat(sample_dfs)
    # # og["rel_ab"] /= og["rel_ab"].sum()
    # # max_abundance = og["rel_ab"].max()
    # # print(max_abundance)
    # # display(og)
    # intervals = [0] + list(logspace(-2, -.3, 5))
    # print(intervals)
    # for i, abundance_limit in enumerate(intervals):
    #     if abundance_limit == 0:  continue
    #     previous_limit = intervals[i-1]
    #     print(previous_limit, abundance_limit)
    #     abundance_range = f"{round(previous_limit, 2)}-{round(abundance_limit, 2)}"
    #     for level in taxonomic_levels:
    #         ab = og.copy()
    #         # ab = ab.set_index(level).drop([l for l in taxonomic_levels if l != level], axis=1)#.groupby(level).sum()
    #         # display(ab)

    #         dic = {}
    #         taxonomies = {}
    #         for seq, row in ab.iterrows():
    #             s = row["sample"]
    #             if s == "J34":  continue

    #             day = sample_days[s]
    #             abundance = row["rel_ab"]
    #             # print(abundance)
    #             if abundance <= previous_limit or abundance > abundance_limit:  continue
    #             # asvs_used.append(row["seq"])
    #             taxonomies[row[level]] = "|".join([str(row[l]) for l in taxonomic_levels
    #                                 if taxonomic_levels.index(l) <= taxonomic_levels.index(level)])
    #             # taxonomies.append("|".join([str(row[l]) for l in taxonomic_levels
    #             #                     if taxonomic_levels.index(l) <= taxonomic_levels.index(level)]))
    #             dic.setdefault(day, {})
    #             dic[day].setdefault(row[level], 0)
    #             dic[day][row[level]] += abundance

    #         # display(dic)
    #         if dic == {}:
    #             print(f"The {level} does not have organisms in {abundance_range}")
    #             continue
    #         elif len(taxonomies) == 1:
    #             print(f"The {level} only has one organism in {abundance_range}:  {dic} and {taxonomies}")
    #             continue

    #         # display(dic)

    #         df = change_columns(DataFrame(dic))
    #         df = df.astype(float).replace([inf, -inf], nan).fillna(0)
    #         # display(df)
    #         # display(df.to_numpy())
    #         max_abundance = max(list(df.to_numpy().flatten()))
    #         # display(df)
    #         # clustermap = ClusterHeatmap(df, level, ) # Series(dict(Counter(taxonomies))))
    #         clusterMap = sns.clustermap(df, 
    #                                     # row_colors=row_colors,
    #                                     cmap=new_cmap,
    #                                     vmin=0,
    #                                     vmax=max_abundance,
    #                                     col_cluster=False, figsize=(20, 20),
    #                                     row_linkage=taxonomy_linkage(Series(taxonomies)),
    #                                     dendrogram_ratio=(.1, .2))
    #         clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)  # Adjust these values as needed to fit labels
    #         clusterMap.figure.suptitle(f"{level.capitalize()} {abundance_range} Abundances", fontsize=30)
    #         clusterMap.ax_heatmap.set_yticklabels(clusterMap.ax_heatmap.get_yticklabels(), fontsize=20, rotation=0)
    #         clusterMap.ax_heatmap.set_xticklabels(clusterMap.ax_heatmap.get_xticklabels(), fontsize=24)
    #         clusterMap.ax_heatmap.set_xlabel("Days after inoculation", fontsize=30)
    #         clusterMap.figure.tight_layout()
    #         clusterMap.figure.savefig(f"abundance_heatmaps/{level}_{abundance_range}_abundances.png")
    #     # break
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Correlate member abundances with methane
    """)
    return


@app.cell
def _(merged_normalized_abundances, normalized_methane):
    merged_correlations = {}
    for index, _row in merged_normalized_abundances.iterrows():
        corr = _row.corr(normalized_methane)
        merged_correlations[index] = corr
    merged_correlations = dict(sorted(merged_correlations.items(), key=lambda item: item[1], reverse=True))  # print(f"{index} correlates {corr} with methane")
    merged_correlations
    return


@app.cell
def _(DataFrame, Series, dump, load):
    from pandas import isna
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))
    _abundances = load(open('model_inputs/abundances.json', 'r'))
    _asvSet_abundances = load(open(f'modeling_files/ASVset_abundances.json', 'r'))
    omitted_columns = {'all_relevant_samples': ['5-10']}
    from scipy.stats import spearmanr

    def _correlations(ser1, ser2):
        aligned1, aligned2 = ser1.align(ser2, join='inner')
        if len(set(aligned1)) == 1 or len(set(aligned2)) == 1:
            global constant_vals
            constant_vals = constant_vals + 1
            return (float('nan'), float('nan'))
        return spearmanr(aligned1, aligned2)
    abundances_df = DataFrame(_abundances).drop(ommitted, axis=1)
    abundances_df = abundances_df.loc[:, (abundances_df.fillna(0) > 0).sum() >= 5]
    ASVset_abundances_df = DataFrame(_asvSet_abundances).drop(ommitted, axis=1).fillna(0)
    ASVset_abundances_df = ASVset_abundances_df.loc[:, ASVset_abundances_df.notna().sum() >= 5]
    print(abundances_df.shape)
    _summary_samples = load(open('model_inputs/measurements/summary_samples.json', 'r'))
    print([_summary_samples[s]['CH4 Production (mol/min)'] for s in _summary_samples.keys()])
    print([_summary_samples[s]['H2 delivery rate (mol/min)'] for s in _summary_samples.keys()])
    for dataName, term in {'CH4%': 'Gas Composition (% CH4)', 'CH4_mol': 'CH4 Production (mol/min)', 'CO2_BT': 'CO2 breakthrough (mol/min)', 'CO2_feed': 'CO2 delivery rate (mol/min)', 'H2_BT': 'H2 breakthrough (mol/min)', 'H2_feed': 'H2 delivery rate (mol/min)'}.items():
        for _name, ommitted in omitted_columns.items():
            constant_vals = 0
            sample_values = {_sample: _content[term] for _sample, _content in _summary_samples.items() if _sample not in ommitted}
            print(dataName)
            ASV_correlations = {}
            constant_vals = 0
            for _ASV, sampleAbun in abundances_df.iterrows():
                _ID = _iterativeIDs.get(_ASV, _ASV)
                correlation, _p = _correlations(sampleAbun, Series(sample_values))
                if isna(correlation):
                    continue
                ASV_correlations[_ID] = {'correlation': correlation, 'p_value': _p}
            ASV_correlations = dict(sorted(ASV_correlations.items(), key=lambda item: item[1]['correlation'], reverse=True))
            print(f'Constant values: {constant_vals}')
            dump(ASV_correlations, open(f'modeling_files/correlations/ASV_correlations_{_name}_{dataName}.json', 'w'))
            IterativeID_correlations = {}
            for _ASV, sampleAbun in ASVset_abundances_df.iterrows():
                _ID = _iterativeIDs.get(_ASV, _ASV)
                correlation, _p = _correlations(sampleAbun, Series(sample_values))
                if isna(correlation):
                    continue
                IterativeID_correlations[_ID] = {'correlation': correlation, 'p_value': _p}
            IterativeID_correlations = dict(sorted(IterativeID_correlations.items(), key=lambda item: item[1]['correlation'], reverse=True))
            dump(IterativeID_correlations, open(f'modeling_files/correlations/IterativeID_correlations_{_name}_{dataName}.json', 'w'))
            print(IterativeID_correlations)
            ASVset_abundances_df.index = [_iterativeIDs.get(_ASV, _ASV).split('.')[0] for _ASV in ASVset_abundances_df.index]
            ASVset_abundances_df.groupby(ASVset_abundances_df.index).sum()
            ASVset_correlations = {}
            for genus, sampleAbun in ASVset_abundances_df.iterrows():
                correlation, _p = _correlations(sampleAbun, Series(sample_values))
                if isna(correlation):
                    continue
                ASVset_correlations[genus] = {'correlation': correlation, 'p_value': _p}
            ASVset_correlations = dict(sorted(ASVset_correlations.items(), key=lambda item: item[1]['correlation'], reverse=True))
            dump(ASVset_correlations, open(f'modeling_files/correlations/ASVset_correlations_{_name}_{dataName}.json', 'w'))
    return ommitted, spearmanr


@app.cell
def _(DataFrame, Series, display, load, read_csv, set_option, sns):
    from matplotlib import colors, patches
    from glob import glob
    from statsmodels.stats.multitest import multipletests
    set_option('display.max_rows', None)
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    _sample_days = {'5-10': 0, 'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}
    order = ['all samples $H_2$ feed', 'all samples $H_2$ BT', 'all samples $CO_2$ feed', 'all samples $CO_2$ BT', 'all samples $CH_4$%', 'all samples $CH_4$ mol']
    _correlations, _pvals = ({}, {})
    for _k in order:
        _correlations.setdefault(_k, {})
        _pvals.setdefault(_k, {})
    for cor in glob('modeling_files/correlations/ASV_correlations_all_relevant_samples_*.json'):
        print(cor)
        _name = cor.split('/')[-1].split('.')[0].split('_correlations_')[1].replace('_', ' ')
        _name = _name.replace('relevant ', '').replace('phase ', '')
        _name = _name.replace('H2', '$H_2$').replace('CO2', '$CO_2$').replace('CH4', '$CH_4$')
        print(_name)
        _content = load(open(cor, 'r'))
        _correlations[_name].update({_k: _v['correlation'] for _k, _v in _content.items()})
        _pvals[_name].update({_k: _v['p_value'] for _k, _v in _content.items()})
    q_vals = {}
    for _name, _content in _pvals.items():
        pvals_list = list(_content.values())
        q_vals[_name] = multipletests(pvals_list, alpha=0.05, method='fdr_bh')
        print(sum(q_vals[_name][0]))
    new_correlations = {_name: dict(inner) for _name, inner in _correlations.items()}
    new_pvals = {_name: dict(inner) for _name, inner in _pvals.items()}
    for _name, _content in _pvals.items():
        for _i, (_k, _v) in enumerate(_content.items()):
            if not q_vals[_name][0][_i]:
                new_correlations[_name].pop(_k)
                new_pvals[_name].pop(_k)
    display(new_correlations)
    display(new_pvals)
    df_4 = DataFrame(new_correlations).fillna(0)
    pval_matrix = DataFrame(new_pvals)
    _rename_full = {
        'all samples $H_2$ feed':   '$H_2$ loading $\\left(\\frac{mol}{min}\\right)$',
        'all samples $H_2$ BT':     '$H_2$ breakthrough $\\left(\\frac{mol}{min}\\right)$',
        'all samples $CO_2$ feed':  '$CO_2$ loading $\\left(\\frac{mol}{min}\\right)$',
        'all samples $CO_2$ BT':    '$CO_2$ breakthrough $\\left(\\frac{mol}{min}\\right)$',
        'all samples $CH_4$%':      '$CH_4$ composition $\\left(\\%\\right)$',
        'all samples $CH_4$ mol':   '$CH_4$ production $\\left(\\frac{mol}{min}\\right)$',
    }
    _rename_short = {k: v.split(' $\\left')[0] for k, v in _rename_full.items()}
    df_4.rename(columns=_rename_full, inplace=True)
    pval_matrix.rename(columns=_rename_short, inplace=True)
    reduced = True
    if reduced:
        df_4 = df_4[(df_4 != float(0)).any(axis=1)]
    _total_df = read_csv('model_inputs/total.csv').set_index('seq')
    _orgs = {}
    for _i in df_4.index:
        _orgs.setdefault(_i.split('.')[0], []).append(_i)
    taxonomies_2 = {}
    level_3 = 'Genus'
    for _seq, _row in _total_df.iterrows():
        _day = _sample_days.get(_row['sample'])
        taxonomy_1 = []
        for l in reversed(_taxonomic_levels):
            _taxa = str(_row[l])
            IDs = _orgs.get(_taxa)
            if IDs is not None:
                break
        if _day is None or _row['rel_ab'] == 0 or IDs is None:
            continue
        taxonomy_string = '|'.join([str(_row[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index(level_3)])
        for _ID in IDs:
            taxonomies_2.setdefault(_ID, taxonomy_string)
    taxonomy_series_2 = Series({_ID: taxonomies_2.get(_ID, f"Unknown|{_ID.split('.')[0]}") for _ID in df_4.index})
    min_max = (round(df_4.min().min(), 1), round(df_4.max().max(), 1))
    _clusterMap = sns.clustermap(df_4, cmap='coolwarm_r', norm=colors.TwoSlopeNorm(vmin=min_max[0], vcenter=0, vmax=min_max[1]), col_cluster=False, clip_on=True, figsize=(20, 80), cbar_kws={'label': 'Correlation'})
    _labelsize = 70
    _clusterMap.ax_heatmap.set_xlabel('Operational metric', fontsize=_labelsize, labelpad=40)
    _clusterMap.ax_heatmap.set_ylabel('ASV', fontsize=_labelsize, labelpad=40)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_heatmap.yaxis.tick_left()
    _clusterMap.ax_heatmap.yaxis.set_label_position('left')
    hm = _clusterMap.ax_heatmap.get_position()
    rd = _clusterMap.ax_row_dendrogram.get_position()
    _clusterMap.ax_heatmap.set_position([rd.x0, hm.y0, hm.x1 - rd.x0, hm.height])
    _cbar = _clusterMap.ax_cbar
    _ticks = _cbar.get_yticks()[::2]
    _ticks[-1] = min_max[1]
    _ticks[0] = min_max[0]
    _ticks = [round(t, 1) for t in _ticks]
    _cbar.set_yticklabels(_ticks, fontsize=40)
    print(min_max, _ticks)
    _cbar.set_yticks(_ticks)
    _cbar.set_xlabel('')
    _cbar.set_ylabel('Spearman $\\rho$', fontsize=60, labelpad=30, rotation=90)
    _cbar.tick_params(labelsize=40, length=12, width=2)
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _cbar_w = 0.03
    _cbar_h = _heatmap_pos.height * 0.6
    _clusterMap.ax_cbar.set_position([
        _heatmap_pos.x1 + 0.05,
        _heatmap_pos.y0 + (_heatmap_pos.height - _cbar_h) / 2,
        _cbar_w,
        _cbar_h,
    ])
    _organisms_to_highlight = ['Methanobacterium', 'Methanosarcina', 'Methanobacteriaceae']
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {_k.split('.')[0]: _v for _k, _v in _iterativeID_levels.items()}
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    for _label in _ylabels:
        if any([x in _label.get_text() for x in _organisms_to_highlight]):
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
        _label.set_rotation(90)
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = list(range(pval_matrix.shape[1]))
    _pvals_reordered = pval_matrix.iloc[_dendrogram_row, _dendrogram_col]
    _ax = _clusterMap.ax_heatmap
    for _i in range(_pvals_reordered.shape[0]):
        for _j in range(_pvals_reordered.shape[1]):
            _p = _pvals_reordered.iloc[_i, _j]
            _value = _clusterMap.data2d.iloc[_i, _j]
            if _value == 0:
                continue
            marker = ''
            _color = 'black' if abs(_value) < 0.7 else 'white'
            _clusterMap.ax_heatmap.text(_j + 0.5, _i + 0.5, f'{_value:.2f}{marker}', ha='center', va='center', color=_color, fontsize=60)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=50, rotation=80)
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=50, rotation=0)
    _clusterMap.figure.savefig(f"correlations_heatmap{('' if not reduced else '_reduced')}.png", bbox_inches='tight', dpi=300)
    return level_3, multipletests, patches


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### All v All correlation matrix
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### Abundance
    """)
    return


@app.cell
def _(DataFrame, display, load, plt, sns):
    from locale import D_FMT
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))
    df_5 = DataFrame(load(open('model_inputs/abundances.json', 'r'))).T
    columns = []
    remove_columns = []
    for _col in df_5.columns:
        if df_5[_col].max() < 0.001:
            remove_columns.append(_col)
            continue
        columns.append(_iterativeIDs[_col])
    df_5.drop(remove_columns, axis=1, inplace=True)
    df_5.columns = [_ID.split('.')[0] for _ID in columns]
    df_5 = df_5.groupby(df_5.columns, axis=1).sum()
    display(df_5)
    _corr_matrix = df_5.corr('spearman')
    display(_corr_matrix)
    plt.figure(figsize=(150, 120))
    sns.heatmap(_corr_matrix, annot=False, cmap='coolwarm', center=0)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.title('Correlation Matrix')
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(
    DataFrame,
    Series,
    colorsys,
    display,
    genera_color_map_1,
    iterativeID_color_map_1,
    level_3,
    load,
    mcolors,
    nan,
    patches,
    sns,
    spearmanr,
):
    from numpy import ones, triu, ones_like
    _days = {'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}

    def _corr_pvalues(df):
        _n = df.shape[1]
        _pvals = DataFrame(ones((_n, _n)), index=df.columns, columns=df.columns)
        for _i in range(_n):
            for _j in range(_i + 1, _n):
                _, _p = spearmanr(df.iloc[:, _i], df.iloc[:, _j])
                _pvals.iloc[_i, _j] = _p
                _pvals.iloc[_j, _i] = _p
        return _pvals
    _significantly_connected_organisms = list(load(open('significantly_connected_organisms.json', 'r')))
    _significantly_connected_organisms.append('Methanobacteriaceae.1')
    _taxa_color_map = load(open(f'iterativeID_color_map.json', 'r'))
    print(_significantly_connected_organisms)
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))
    df_6 = DataFrame(load(open('model_inputs/abundances.json', 'r'))).T
    df_6.drop(df_6.index.difference(_days.keys()), inplace=True)
    df_6.columns = [_iterativeIDs.get(_ID, _ID) for _ID in df_6.columns]
    df_6.drop([_col for _col in df_6.columns if _col not in _significantly_connected_organisms], axis=1, inplace=True)
    df_6.drop([_col for _col in df_6.columns if df_6[_col].max() < 0.0005], axis=1, inplace=True)
    total_captured = df_6.sum(axis=1)
    print(sum(total_captured), total_captured)
    display(df_6)
    _corr_matrix = df_6.corr('spearman')
    display(_corr_matrix)
    pval_matrix_1 = _corr_pvalues(df_6)
    from statsmodels.stats.multitest import multipletests as _multipletests
    from numpy import triu_indices as _triu_indices, full_like as _full_like, isfinite as _isfinite, ones_like as _ones_like
    _p_arr = pval_matrix_1.values
    _n = _p_arr.shape[0]
    _iu, _ju = _triu_indices(_n, k=1)
    _flat_p = _p_arr[_iu, _ju]
    _finite = _isfinite(_flat_p)
    _q_arr = _full_like(_flat_p, float("nan"))
    if _finite.any():
        _, _q_finite, _, _ = _multipletests(_flat_p[_finite], alpha=0.05, method='fdr_bh')
        _q_arr[_finite] = _q_finite
    _q_matrix = _ones_like(_p_arr)
    for _k, (_i_q, _j_q) in enumerate(zip(_iu, _ju)):
        _q_matrix[_i_q, _j_q] = _q_arr[_k]
        _q_matrix[_j_q, _i_q] = _q_arr[_k]
    pval_matrix_1 = DataFrame(_q_matrix, index=pval_matrix_1.index, columns=pval_matrix_1.columns)
    print(f"BH-FDR survivors at q<0.05: {(_q_arr < 0.05).sum()} / {_finite.sum()} pairs")

    taxonomy_2 = {_iterativeIDs.get(_k, _k): _v for _k, _v in load(open('model_inputs/taxonomy.json', 'r')).items()}
    taxonomies_3 = {}
    for _col in df_6.columns:
        taxonomies_3[_col] = '|'.join([_v for _k, _v in taxonomy_2.get(_col, 'Unknown').items() if _k != 'Species' and _v is not None])
    taxonomy_series_3 = Series({idx: taxonomies_3.get(idx, f'Unknown|{idx}') for idx in df_6.columns})
    print(f"Detected depth: {max((len(t.split('|')) for t in taxonomy_series_3))}")
    print(f'df rows: {len(df_6.columns)}')
    print(f'taxonomy_series length: {len(taxonomy_series_3)}')
    print(f'Sample entries:\n{taxonomy_series_3.head()}')
    DEFAULT_COLOR = 'lightgray'

    def _lighten(c, factor):
        h, l, s = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
        return colorsys.hls_to_rgb(h, max(0.0, min(1.0, l * factor)), s)
    proteo_base = iterativeID_color_map_1.get('Proteobacteria', 'tab:purple')
    proteo_classes = sorted({str(taxonomies_3.get(_i, '')).split('|')[2] for _i in df_6.index if len(str(taxonomies_3.get(_i, '')).split('|')) >= 3 and str(taxonomies_3[_i]).split('|')[1] == 'Proteobacteria'})
    _n = max(len(proteo_classes), 1)
    proteo_class_color = {cls: _lighten(proteo_base, 0.6 + 0.8 * _i / max(_n - 1, 1)) for _i, cls in enumerate(proteo_classes)}

    def lookup(idx):
        parts = str(taxonomies_3.get(idx, '')).split('|')
        if len(parts) >= 3 and parts[1] == 'Proteobacteria':
            return proteo_class_color.get(parts[2], proteo_base)
        if idx in iterativeID_color_map_1:
            return iterativeID_color_map_1[idx]
        elif idx in genera_color_map_1:
            return genera_color_map_1[idx]
        return DEFAULT_COLOR
    _row_colors = Series({idx: lookup(idx) for idx in _corr_matrix.index}, name='Phylum')
    col_colors = Series({idx: lookup(idx) for idx in _corr_matrix.columns}, name='Phylum')
    _clusterMap = sns.clustermap(_corr_matrix, row_colors=_row_colors, col_colors=col_colors, cbar_pos=None, cmap='coolwarm_r', center=0, figsize=(60, 70), dendrogram_ratio=(0.1, 0.2))
    _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    secondary_labels = [tick.get_text() for tick in _clusterMap.ax_heatmap.get_yticklabels()]
    _clusterMap.ax_heatmap.yaxis.set_ticks_position('left')
    _clusterMap.ax_heatmap.yaxis.set_label_position('left')
    _clusterMap.ax_heatmap.set_yticklabels(secondary_labels, rotation=0)
    hm_pos = _clusterMap.ax_heatmap.get_position()
    fig_w = _clusterMap.figure.get_figwidth()
    fig_h = _clusterMap.figure.get_figheight()
    strip_w = 0.015
    strip_h = 0.015
    _clusterMap.ax_row_colors.set_position([hm_pos.x0 - strip_w, hm_pos.y0, strip_w, hm_pos.height])
    _clusterMap.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0 - strip_h, hm_pos.width, strip_h])
    y_pad = strip_w * fig_w * 72 + 15
    _clusterMap.ax_heatmap.tick_params(axis='y', pad=y_pad)
    x_pad = strip_h * fig_h * 72 + 15
    _clusterMap.ax_heatmap.tick_params(axis='x', pad=x_pad)
    _labelsize = 40
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=_labelsize, rotation=0)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=_labelsize, rotation=60, ha='right', rotation_mode='anchor')
    _organisms_to_highlight = ['Methanobacterium', 'Methanosarcina', 'Methanobacteriaceae']
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    _orgs = set()
    for _label in _ylabels:
        if any([x in _label.get_text() for x in _organisms_to_highlight]):
            _orgs.add(_label.get_text())
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text()
        if _iterativeID_levels.get(_text) == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels)
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {_k.split('.')[0]: _v for _k, _v in _iterativeID_levels.items()}
    _xlabels = _clusterMap.ax_heatmap.get_xticklabels()
    for _label in _xlabels:
        if any([x in _label.get_text() for x in _organisms_to_highlight]):
            _orgs.add(_label.get_text())
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text()
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_xticklabels(_xlabels)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = _clusterMap.dendrogram_col.reordered_ind
    _one_triangle = True
    if _one_triangle:
        _df_reordered = _corr_matrix.iloc[_dendrogram_row, _dendrogram_col]
        _mask = triu(ones_like(_df_reordered, dtype=bool), k=1)
        _mesh = _clusterMap.ax_heatmap.collections[0]
        _arr = _mesh.get_array().reshape(_df_reordered.shape)
        _arr[_mask] = nan
        _mesh.set_array(_arr.ravel())
    for _org in _orgs:
        _orgIx = _corr_matrix.index.get_loc(_org)
        if _orgIx not in _dendrogram_row:
            continue
        _row_pos = _dendrogram_row.index(_orgIx)
        _rect = patches.Rectangle((0, _row_pos), len(_corr_matrix.columns) if not _one_triangle else _row_pos + 1, 1, linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
        _colIx = _corr_matrix.columns.get_loc(_org)
        _col_pos = _dendrogram_col.index(_colIx)
        _rect = patches.Rectangle((_col_pos, 0) if not _one_triangle else (_col_pos, len(_corr_matrix.index)), 1, len(_corr_matrix.index) if not _one_triangle else -(len(_corr_matrix.index) - _col_pos), linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
    _pvals_reordered = pval_matrix_1.iloc[_dendrogram_row, _dendrogram_col]
    _ax = _clusterMap.ax_heatmap
    for _i in range(_pvals_reordered.shape[0]):
        for _j in range(_pvals_reordered.shape[1]):
            if _pvals_reordered.iloc[_i, _j] < 0.05 and (not _one_triangle or _mask[_i, _j] == 0):
                _ax.text(_j + 0.5, _i + 0.5, '*', color='lightgreen', fontsize=60, fontweight='bold', ha='center', va='center')
    _clusterMap.ax_heatmap.set_xlabel('Member ASVs', fontsize=50)
    _clusterMap.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)
    fig_width = _clusterMap.figure.get_figwidth()
    phylum_color = {}
    for idx in _corr_matrix.index:
        parts = str(taxonomies_3.get(idx, '')).split('|')
        if len(parts) >= 2 and parts[1] not in ('None', '', 'Unknown', 'nan'):
            phylum_color.setdefault(parts[1], _row_colors.get(idx))
    shown_phyla = set()
    for idx in _corr_matrix.index:
        parts = str(taxonomies_3.get(idx, '')).split('|')
        if len(parts) >= 2 and parts[1] not in ('None', '', 'Unknown', 'nan'):
            shown_phyla.add(parts[1])
    archaea_markers = ('archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota')

    def is_archaea(p):
        return any((m in _p.lower() for m in archaea_markers))
    archaea_phyla = sorted((_p for _p in shown_phyla if is_archaea(_p)))
    bacteria_phyla = sorted((_p for _p in shown_phyla if not is_archaea(_p)))
    handle = []
    if len(archaea_phyla) > 0:
        _archaea_patches = [patches.Patch(color=phylum_color.get(_p, DEFAULT_COLOR), label=_p) for _p in archaea_phyla]
        handle = [header_patch('Archaea')] + _archaea_patches
    if len(bacteria_phyla) > 0:
        _bacteria_patches = [patches.Patch(color=phylum_color.get(_p, DEFAULT_COLOR), label=_p) for _p in bacteria_phyla]
        handle = handle + ([header_patch('Bacteria')] + _bacteria_patches)

    def header_patch(title):
        """Section header: invisible swatch with a bold label."""
        return patches.Patch(color='none', label=f'$\\bf{{{title}}}$')
    _legend_handles = handle
    _clusterMap.ax_heatmap.legend(handles=_legend_handles, title=level_3, title_fontsize=8 * (fig_width / 10), loc='lower left', bbox_to_anchor=(0.5, 0.6), fontsize=7 * (fig_width / 10), frameon=True)
    _clusterMap.figure.savefig(f"abundance_heatmaps/abundance_correlatons{('_one_triangle' if _one_triangle else '')}.png", dpi=300, bbox_inches='tight')
    return (
        archaea_phyla,
        bacteria_phyla,
        header_patch,
        ones,
        ones_like,
        pval_matrix_1,
        taxonomy_series_3,
        triu,
    )


@app.cell
def _(pval_matrix_1):
    pval_matrix_1[pval_matrix_1 < 0.05]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
 
    """)
    return


@app.cell
def _(KMeans, best_k, corr_df, pd, sig_mask):
    # sig_mask is a boolean DataFrame (True = p < threshold)
    corr_masked = corr_df.where(sig_mask, 0)
    km = KMeans(n_clusters=best_k, n_init=20, random_state=42)
    cluster_labels = pd.Series(km.fit_predict(corr_masked), index=corr_df.index)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### methanogen correlations
    """)
    return


@app.cell
def _(
    DataFrame,
    display,
    load,
    ones,
    patches,
    sns,
    spearmanr,
    taxonomy_series_3,
):
    _days = {'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}

    def _corr_pvalues(df):
        _n = df.shape[1]
        _pvals = DataFrame(ones((_n, _n)), index=df.columns, columns=df.columns)
        for _i in range(_n):
            for _j in range(_i + 1, _n):
                _, _p = spearmanr(df.iloc[:, _i], df.iloc[:, _j])
                _pvals.iloc[_i, _j] = _p
                _pvals.iloc[_j, _i] = _p
        return _pvals
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json', 'r'))
    df_7 = DataFrame(load(open('model_inputs/abundances.json', 'r'))).T
    df_7.drop(df_7.index.difference(_days.keys()), inplace=True)
    df_7.columns = [_iterativeIDs.get(_ID, _ID) for _ID in df_7.columns]
    df_7 = df_7.groupby(df_7.columns, axis=1).sum()
    taxonomy_3 = {_iterativeIDs.get(_k, _k): _v for _k, _v in load(open('model_inputs/taxonomy.json', 'r')).items()}
    taxonomies_4 = {_col: '|'.join([_v for _k, _v in taxonomy_3.get(_col, {}).items() if _k != 'Species' and _v is not None]) for _col in df_7.columns}
    df_7.drop([_col for _col in df_7.columns if taxonomies_4.get(_col, '|').split('|')[0] != 'Archaea' or 'midas' in _col], axis=1, inplace=True)
    _corr_matrix = df_7.corr('spearman').dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    display(_corr_matrix)
    pval_matrix_2 = _corr_pvalues(df_7).dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    from statsmodels.stats.multitest import multipletests as _multipletests
    from numpy import triu_indices as _triu_indices, full_like as _full_like, isfinite as _isfinite, ones_like as _ones_like
    _p_arr = pval_matrix_2.values
    _n = _p_arr.shape[0]
    _iu, _ju = _triu_indices(_n, k=1)
    _flat_p = _p_arr[_iu, _ju]
    _finite = _isfinite(_flat_p)
    _q_arr = _full_like(_flat_p, float("nan"))
    if _finite.any():
        _, _q_finite, _, _ = _multipletests(_flat_p[_finite], alpha=0.05, method='fdr_bh')
        _q_arr[_finite] = _q_finite
    _q_matrix = _ones_like(_p_arr)
    for _k, (_i_q, _j_q) in enumerate(zip(_iu, _ju)):
        _q_matrix[_i_q, _j_q] = _q_arr[_k]
        _q_matrix[_j_q, _i_q] = _q_arr[_k]
    pval_matrix_2 = DataFrame(_q_matrix, index=pval_matrix_2.index, columns=pval_matrix_2.columns)
    print(f"BH-FDR survivors at q<0.05: {(_q_arr < 0.05).sum()} / {_finite.sum()} pairs")

    pval_matrix_2 = pval_matrix_2[_corr_matrix.columns]
    pval_matrix_2 = pval_matrix_2.loc[_corr_matrix.index]
    display(pval_matrix_2)
    print(f'df rows: {len(df_7.columns)}')
    print(f'taxonomy_series length: {len(taxonomy_series_3)}')
    print(f'Sample entries:\n{taxonomy_series_3.head()}')
    _clusterMap = sns.clustermap(_corr_matrix, cmap='coolwarm_r', center=0, figsize=(60, 60), dendrogram_ratio=(0.1, 0.2))
    _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = _clusterMap.dendrogram_col.reordered_ind
    _organisms_to_highlight = ['Methanobacterium.2', 'Methanobacteriaceae.1', 'Methanobacterium.1']
    for _org in _organisms_to_highlight:
        _orgIx = _corr_matrix.index.get_loc(_org)
        if _orgIx not in _dendrogram_row:
            continue
        _row_pos = _dendrogram_row.index(_orgIx)
        _rect = patches.Rectangle((0, _row_pos), len(_corr_matrix.columns), 1, linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
        _colIx = _corr_matrix.columns.get_loc(_org)
        _col_pos = _dendrogram_col.index(_colIx)
        _rect = patches.Rectangle((_col_pos, 0), 1, len(_corr_matrix.index), linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
    _pvals_reordered = pval_matrix_2.iloc[_dendrogram_row, _dendrogram_col]
    _ax = _clusterMap.ax_heatmap
    for _i in range(_pvals_reordered.shape[0]):
        for _j in range(_pvals_reordered.shape[1]):
            if _pvals_reordered.iloc[_i, _j] < 0.05:
                _ax.add_patch(patches.Rectangle((_j, _i), 1, 1, fill=False, edgecolor='lightgreen', linewidth=8))
    _cbar = _clusterMap.ax_cbar
    _ticks = _cbar.get_yticks()
    _ticks[-1] = round(_corr_matrix.max().max())
    _ticks[0] = round(_corr_matrix.min().min(), 2)
    _cbar.set_yticklabels(_ticks, fontsize=35)
    _cbar.set_yticks(_ticks)
    _cbar.set_xlabel('Spearman', fontsize=40, labelpad=20)
    _labelsize = 40
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=_labelsize, rotation=0)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=_labelsize, rotation=60, ha='right', rotation_mode='anchor')
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    for _label in _ylabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels)
    _xlabels = _clusterMap.ax_heatmap.get_xticklabels()
    for _label in _xlabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
    _clusterMap.ax_heatmap.set_xticklabels(_xlabels)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.96, _heatmap_pos.y0 + 0.5, 0.06, _heatmap_pos.height / 5])
    _clusterMap.ax_heatmap.set_xlabel('Member ASVs', fontsize=50)
    _clusterMap.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)
    _clusterMap.figure.savefig('abundance_heatmaps/methanogens_abundance_correlatons.png', bbox_inches='tight', dpi=300)
    return


@app.cell
def _(
    DataFrame,
    display,
    inf,
    load,
    nan,
    ones,
    ones_like,
    patches,
    sns,
    spearmanr,
    triu,
):
    _days = {'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}

    def _corr_pvalues(df):
        _n = df.shape[1]
        _pvals = DataFrame(ones((_n, _n)), index=df.columns, columns=df.columns)
        for _i in range(_n):
            for _j in range(_i + 1, _n):
                _, _p = spearmanr(df.iloc[:, _i], df.iloc[:, _j])
                _pvals.iloc[_i, _j] = _p
                _pvals.iloc[_j, _i] = _p
        return _pvals
    _nonzero_per_day = load(open('model_inputs/nonzero_per_day.json', 'r'))
    _top_per_day = {}
    _all_orgs = {}
    _topNum = 10
    for _day, _org_dict in _nonzero_per_day.items():
        _orgs = dict(list(_org_dict.items())[:_topNum])
        _top_per_day[_day] = _orgs
        _all_orgs.update(_orgs)
    df_8 = DataFrame(_top_per_day).T
    df_8 = df_8.astype(float).replace([inf, -inf], nan)
    df_8 = df_8.loc[:, df_8.notna().sum() >= 5]
    display(df_8)
    _corr_matrix = df_8.corr('spearman').dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    display(_corr_matrix)
    pval_matrix_3 = _corr_pvalues(df_8).dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    from statsmodels.stats.multitest import multipletests as _multipletests
    from numpy import triu_indices as _triu_indices, full_like as _full_like, isfinite as _isfinite, ones_like as _ones_like
    _p_arr = pval_matrix_3.values
    _n = _p_arr.shape[0]
    _iu, _ju = _triu_indices(_n, k=1)
    _flat_p = _p_arr[_iu, _ju]
    _finite = _isfinite(_flat_p)
    _q_arr = _full_like(_flat_p, float("nan"))
    if _finite.any():
        _, _q_finite, _, _ = _multipletests(_flat_p[_finite], alpha=0.05, method='fdr_bh')
        _q_arr[_finite] = _q_finite
    _q_matrix = _ones_like(_p_arr)
    for _k, (_i_q, _j_q) in enumerate(zip(_iu, _ju)):
        _q_matrix[_i_q, _j_q] = _q_arr[_k]
        _q_matrix[_j_q, _i_q] = _q_arr[_k]
    pval_matrix_3 = DataFrame(_q_matrix, index=pval_matrix_3.index, columns=pval_matrix_3.columns)
    print(f"BH-FDR survivors at q<0.05: {(_q_arr < 0.05).sum()} / {_finite.sum()} pairs")

    pval_matrix_3 = pval_matrix_3[_corr_matrix.columns]
    pval_matrix_3 = pval_matrix_3.loc[_corr_matrix.index]
    display(pval_matrix_3)
    _clusterMap = sns.clustermap(_corr_matrix, cmap='coolwarm_r', center=0, figsize=(60, 60), col_cluster=True, row_cluster=True, dendrogram_ratio=(0.1, 0.2))
    _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = _clusterMap.dendrogram_col.reordered_ind
    _one_triangle = True
    if _one_triangle:
        _df_reordered = _corr_matrix.iloc[_dendrogram_row, _dendrogram_col]
        _mask = triu(ones_like(_df_reordered, dtype=bool), k=1)
        _mesh = _clusterMap.ax_heatmap.collections[0]
        _arr = _mesh.get_array().reshape(_df_reordered.shape)
        _arr[_mask] = nan
        _mesh.set_array(_arr.ravel())
    _organisms_to_highlight = ['Methanobacterium.2', 'Methanobacteriaceae.1', 'Methanobacterium.1']
    for _org in _organisms_to_highlight:
        _orgIx = _corr_matrix.index.get_loc(_org)
        if _orgIx not in _dendrogram_row:
            continue
        _row_pos = _dendrogram_row.index(_orgIx)
        _rect = patches.Rectangle((0, _row_pos), len(_corr_matrix.columns) if not _one_triangle else _row_pos + 1, 1, linewidth=16, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
        _colIx = _corr_matrix.columns.get_loc(_org)
        _col_pos = _dendrogram_col.index(_colIx)
        _rect = patches.Rectangle((_col_pos, 0) if not _one_triangle else (_col_pos, len(_corr_matrix.index)), 1, len(_corr_matrix.index) if not _one_triangle else -(len(_corr_matrix.index) - _col_pos), linewidth=16, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
    _cbar = _clusterMap.ax_cbar
    _ticks = _cbar.get_yticks()
    _ticks[-1] = round(_corr_matrix.max().max())
    _ticks[0] = round(_corr_matrix.min().min(), 2)
    _cbar.set_yticklabels(_ticks, fontsize=50)
    _cbar.set_yticks(_ticks)
    _cbar.set_ylabel('Spearman $\\rho$', fontsize=60, labelpad=15, rotation=90)
    _labelsize = 60
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=_labelsize, rotation=0)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=_labelsize, rotation=60, ha='right', rotation_mode='anchor')
    for _i in range(_clusterMap.data2d.shape[0]):
        for _j in range(_clusterMap.data2d.shape[1]):
            if _mask[_i, _j] == 1:
                continue
            _value = _clusterMap.data2d.iloc[_i, _j]
            if _value == 0:
                continue
            _color = 'black' if abs(_value) < 0.7 else 'white'
            _clusterMap.ax_heatmap.text(_j + 0.5, _i + 0.5, f'{_value:.2f}', ha='center', va='center', color=_color, fontsize=60)
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {_k.split('.')[0]: _v for _k, _v in _iterativeID_levels.items()}
    _clusterMap.ax_heatmap.yaxis.set_ticks_position('left')
    _clusterMap.ax_heatmap.yaxis.set_label_position('left')
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    for _label in _ylabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels, rotation=0)
    _xlabels = _clusterMap.ax_heatmap.get_xticklabels()
    for _label in _xlabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_xticklabels(_xlabels)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 1.01, _heatmap_pos.y0 + 0.5, 0.06, _heatmap_pos.height / 5])
    if _one_triangle:
        _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.3, _heatmap_pos.y0 + 0.4, 0.06, _heatmap_pos.height / 5])
    _clusterMap.ax_heatmap.set_xlabel('ASVs', fontsize=100, labelpad=20)
    _clusterMap.ax_heatmap.set_ylabel('ASVs', fontsize=100, labelpad=20)
    _clusterMap.figure.savefig(f"abundance_heatmaps/Top_{_topNum}_ASVs_abundance_correlation{('_one_triangle' if _one_triangle else '')}.png", bbox_inches='tight', dpi=300)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # co-occurrence figure
    """)
    return


@app.cell
def _(
    DataFrame,
    defaultdict,
    display,
    dump,
    load,
    multipletests,
    order_color_map,
    plt,
    spearmanr,
):
    import networkx as nx
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches
    import math
    import numpy as np
    from itertools import combinations
    _days = {'H': 140, 'Ha': 147, 'I': 149, 'Ia': 153, 'J': 161, 'Ja': 168, 'Jb': 174, 'K': 182, 'L': 189, 'M': 210, 'N': 238, 'Na': 262, 'P': 283, 'Q': 300, 'Qa': 315}
    _iterativeIDs = load(open('model_inputs/iterativeIDs.json'))
    taxonomy_4 = load(open('model_inputs/taxonomy.json'))
    iterativeID_taxonomy = {_iterativeIDs.get(_ID, _ID): taxonomy_4.get(_ID, 'Unknown') for _ID in _iterativeIDs}
    level_4 = 'Phylum'
    iterativeID_level = {_ID: _content.get(level_4, 'Unknown') for _ID, _content in iterativeID_taxonomy.items()}
    genera_level = False
    if genera_level:
        iterativeID_level = iterativeID_level = {_ID.split('.')[0]: _k for _ID, _k in iterativeID_level.items()}
    none_keys = [_k for _k, _v in iterativeID_level.items() if _v is None]
    print('None taxonomy entries:', none_keys)
    archaea_phyla_1 = sorted({_v for _k, _v in iterativeID_level.items() if _v is not None and ('archaeo' in _v.lower() or any((a in _v.lower() for a in ['euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'candidatus thermoplasmatota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota'])))})
    bacteria_phyla_1 = sorted({_v for _k, _v in iterativeID_level.items() if _v is not None and _v not in archaea_phyla_1})
    n_archaea = len(archaea_phyla_1)
    archaea_colors = [plt.cm.turbo(_i / max(n_archaea, 1) * 0.15) for _i in range(n_archaea)]
    n_bacteria = len(bacteria_phyla_1)
    bacteria_colors = [plt.cm.turbo(0.2 + _i / max(n_bacteria, 1) * 0.8) for _i in range(n_bacteria)]
    _taxa_color_map = {}
    for phylum, _color in zip(archaea_phyla_1, archaea_colors):
        _taxa_color_map[phylum] = _color
    for phylum, _color in zip(bacteria_phyla_1, bacteria_colors):
        _taxa_color_map[phylum] = _color
    dump(_taxa_color_map, open(f'{level_4}_color_map.json', 'w'))
    iterativeID_color_map_2 = {_ID: _taxa_color_map[phylum] for _ID, phylum in iterativeID_level.items() if phylum}
    dump(iterativeID_color_map_2, open(f'iterativeID_color_map.json', 'w'))
    df_9 = DataFrame(load(open('model_inputs/abundances.json', 'r'))).T
    df_9.drop(df_9.index.difference(_days.keys()), inplace=True)
    df_9 = df_9.loc[:, (df_9.fillna(0) > 0).sum() >= 3]
    df_9.columns = [_iterativeIDs.get(_ID, _ID) for _ID in df_9.columns]
    if genera_level:
        df_9.columns = [_col.split('.')[0] for _col in df_9.columns]
        df_9 = df_9.groupby(df_9.columns, axis=1).sum()
    display(df_9)
    relative_abundance = df_9.div(df_9.sum(axis=1), axis=0)
    mean_rel_abund = relative_abundance.mean(axis=0)
    presence = (df_9 > 0).astype(int)
    cooccurrence = defaultdict(int)
    for _sample in presence.itertuples(index=False):
        present = [_col for _col, val in zip(presence.columns, _sample) if val]
        for pair in combinations(sorted(present), 2):
            cooccurrence[pair] = cooccurrence[pair] + 1
    pair_data = []
    for (a, b), _count in cooccurrence.items():
        rho, p_value = spearmanr(df_9[a], df_9[b])
        if np.isnan(rho):
            continue
        pair_data.append((a, b, rho, p_value, _count))
    _pvals = np.array([t[3] for t in pair_data])
    reject, qvals, _, _ = multipletests(_pvals, alpha=0.05, method='fdr_bh')
    G = nx.Graph()
    for (a, b, rho, _p, _count), q, keep in zip(pair_data, qvals, reject):
        if keep:
            G.add_edge(a, b, weight=abs(rho), rho=rho, pvalue=_p, qvalue=q, cooccurrence=_count)
    print(f'Tests run: {len(pair_data)}')
    print(f'FDR-significant pairs (q < 0.05): {int(reject.sum())}')
    print(f'Edges after FDR correction: {G.number_of_edges()}')
    print(f'Nodes in graph: {G.number_of_nodes()}')
    print('nodes:', G.number_of_nodes(), 'edges:', G.number_of_edges())
    print('max cooccurrence:', max(cooccurrence.values()) if cooccurrence else 0)
    pos = nx.spring_layout(G, seed=42, iterations=300, k=4.0, scale=1.5)   # expanded layout
    edges = G.edges(data=True)
    rho_values = [d['rho'] for _, _, d in edges]
    edge_widths = [3 * d['weight'] for _, _, d in edges]
    print(sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10])
    norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = plt.get_cmap('coolwarm_r')
    edge_colors = [cmap(norm(r)) for r in rho_values]
    width = 40
    height = 30
    fig, _ax = plt.subplots(figsize=(width, height))
    scale = 5000 * (width / 10) * 2
    compressor = np.sqrt
    node_sizes = [scale * compressor(mean_rel_abund.get(_n, 0)) for _n in G.nodes()]
    _node_colors = [_taxa_color_map.get(iterativeID_level.get(_n, 'Unknown'), 'lightgray') for _n in G.nodes()]
    print('unique node colors:', len(set(_node_colors)))
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=_node_colors, alpha=0.9, ax=_ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors, alpha=0.85, ax=_ax)
    _significantly_connected_organisms = list(pos.keys())
    print(_significantly_connected_organisms)
    dump(_significantly_connected_organisms, open('significantly_connected_organisms.json', 'w'))
    node_connectivity = dict(sorted({_n: int(G.degree(_n)) for _n in G.nodes()}.items(), key=lambda kv: kv[1], reverse=True))
    dump(node_connectivity, open('node_connectivity.json', 'w'), indent=2)
    import matplotlib.patheffects as path_effects
    LABEL_MIN_ABUND = 0.002
    for _n, (x, y) in pos.items():
        abund = mean_rel_abund.get(_n, 0)
        if abund < LABEL_MIN_ABUND:
            continue
        font_size = 6 + 14 * compressor(abund) / compressor(mean_rel_abund.max()) * (width / 10)
        font_size = max(5, min(font_size, 48))
        _txt = _ax.text(x, y, str(_n), fontsize=font_size, color='black', fontweight='bold', ha='center', va='center')
        _txt.set_path_effects([path_effects.Stroke(linewidth=max(1.5, font_size / 6), foreground='white'), path_effects.Normal()])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    _cbar = plt.colorbar(sm, ax=_ax, shrink=0.6, pad=0.02)
    _cbar.set_label('Spearman $\\rho$', fontsize=10 * (width / 10))
    _cbar.ax.tick_params(labelsize=8 * (width / 10), length=8, width=2)
    legend_entries = [(0.002, '0.2%'), (0.02, '2%'), (0.2, '20%')]
    sizes_pt2 = [scale * compressor(a) for a, _ in legend_entries]
    diameters_pt = [2 * np.sqrt(s / np.pi) for s in sizes_pt2]
    breather_pt = 8
    offsets_pt = [0.0]
    for _i in range(1, len(legend_entries)):
        offsets_pt.append(offsets_pt[-1] + (diameters_pt[_i - 1] + diameters_pt[_i]) / 2 + breather_pt)
    fig = _ax.figure
    fig_h_pts = fig.get_figheight() * 72
    top_y = 0.9
    circle_x = 0.85
    label_x = 0.88
    _ax.text(circle_x, top_y + 0.025, 'Mean rel. abundance', transform=fig.transFigure, va='bottom', fontweight='bold', fontsize=6 * (width / 10), clip_on=False)
    for (abund, _label), s, off in zip(legend_entries, sizes_pt2, offsets_pt):
        y = top_y - off / fig_h_pts
        _ax.scatter([circle_x], [y], s=s, color='slategray', alpha=0.9, transform=fig.transFigure, clip_on=False)
        _ax.text(label_x, y, _label, transform=fig.transFigure, va='center', fontsize=5 * (width / 10), clip_on=False)
    _archaea_patches = [mpatches.Patch(color=order_color_map[_p], label=_p) for _p in archaea_phyla_1 if _p in order_color_map]
    _bacteria_patches = [mpatches.Patch(color=order_color_map[_p], label=_p) for _p in bacteria_phyla_1 if _p in order_color_map]

    def header_patch_1(title):
        return mpatches.Patch(color='none', label=f'$\\bf{{{title}}}$')
    _legend_handles = [header_patch_1('Archaea')] + _archaea_patches + [header_patch_1('Bacteria')] + _bacteria_patches
    _ax.legend(handles=_legend_handles, title='Taxonomic ' + level_4, title_fontsize=8 * (width / 10), loc='lower left', bbox_to_anchor=(-0.24, 0.1), fontsize=7 * (width / 10), frameon=True)
    _ax.axis('off')
    plt.tight_layout()
    plt.savefig(f'cooccurrence_network_p_value_FDR.png', dpi=300, bbox_inches='tight')
    plt.show()
    return G, iterativeID_level, mcolors, mpatches, np


if __name__ == "__main__":
    app.run()
