"""Regenerate the phylum / iterativeID colour maps for ALL 3,276 ASVs.

Replicates the colour logic embedded in data_processing.py so the maps match the
published figures, but runs it up front (standalone) over the complete, phylum-
backfilled taxonomy so every iterativeID has a colour BEFORE any heatmap cell
reads the map from disk (fixes the end-of-notebook ordering gotcha).

Outputs (overwrite):
  model_inputs/phylum_colors.json   { phylum: "#hex" }     (build_phylum_colors)
  Phylum_color_map.json             { phylum: [r,g,b,a] }  (co-occurrence cell)
  iterativeID_color_map.json        { iterativeID: [r,g,b,a] }  (per-phylum colour, all IDs)
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path("model_inputs")

# markers used by data_processing.py:build_phylum_colors (line 521)
BPC_ARCHAEA = ['archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota',
               'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota']
# markers used by the co-occurrence cell (line 1738)
COOC_ARCHAEA = ['euryarchaeota', 'crenarchaeota', 'thaumarchaeota',
                'candidatus thermoplasmatota', 'halobacterota', 'methanobacteriota',
                'micrarchaeota', 'nanoarchaeota']


def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(*(int(c * 255) for c in rgba[:3]))


def split_archaea(phyla, markers, allow_archaeo_substring):
    def is_arch(p):
        pl = p.lower()
        if allow_archaeo_substring and 'archaeo' in pl:
            return True
        return any(m in pl for m in markers)
    archaea = sorted(p for p in phyla if is_arch(p))
    bacteria = sorted(p for p in phyla if p not in archaea)
    return archaea, bacteria


def turbo_palette(archaea, bacteria):
    """archaea in turbo[0,0.15], bacteria in turbo[0.2,1.0] (data_processing.py)."""
    colors = {}
    for i, p in enumerate(archaea):
        colors[p] = plt.cm.turbo(i / max(len(archaea), 1) * 0.15)
    for i, p in enumerate(bacteria):
        colors[p] = plt.cm.turbo(0.2 + i / max(len(bacteria), 1) * 0.8)
    return colors


def main():
    taxonomy = json.load(open(ROOT / "taxonomy.json"))
    iterativeIDs = json.load(open(ROOT / "iterativeIDs.json"))   # seq -> ID

    # --- phylum_colors.json (hex), build_phylum_colors logic ---
    phyla = sorted({t.get('Phylum') for t in taxonomy.values() if t.get('Phylum')})
    a1, b1 = split_archaea(phyla, BPC_ARCHAEA, allow_archaeo_substring=False)
    pal_hex = {p: rgba_to_hex(c) for p, c in turbo_palette(a1, b1).items()}
    json.dump(pal_hex, open(ROOT / "phylum_colors.json", "w"), indent=2)
    print(f"wrote model_inputs/phylum_colors.json ({len(pal_hex)} phyla)")

    # --- Phylum_color_map.json (rgba) + iterativeID_color_map.json, co-occurrence logic ---
    # per-iterativeID phylum
    iterativeID_phylum = {}
    for seq, ID in iterativeIDs.items():
        iterativeID_phylum[ID] = taxonomy[seq].get('Phylum')
    present_phyla = {p for p in iterativeID_phylum.values() if p}
    a2, b2 = split_archaea(present_phyla, COOC_ARCHAEA, allow_archaeo_substring=True)
    pal_rgba = {p: list(c) for p, c in turbo_palette(a2, b2).items()}
    json.dump(pal_rgba, open("Phylum_color_map.json", "w"))
    print(f"wrote Phylum_color_map.json ({len(pal_rgba)} phyla; {len(a2)} archaea / {len(b2)} bacteria)")

    iterativeID_color_map = {ID: pal_rgba[ph] for ID, ph in iterativeID_phylum.items() if ph}
    json.dump(iterativeID_color_map, open("iterativeID_color_map.json", "w"))
    missing = [ID for ID, ph in iterativeID_phylum.items() if not ph]
    print(f"wrote iterativeID_color_map.json ({len(iterativeID_color_map)}/{len(iterativeIDs)} IDs)")
    if missing:
        print(f"  WARNING: {len(missing)} IDs have no phylum colour (e.g. {missing[:3]})")
    else:
        print("  every iterativeID has a colour (no bare-index KeyError risk)")


if __name__ == "__main__":
    main()
