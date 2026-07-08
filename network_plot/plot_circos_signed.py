import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as pe
from matplotlib.patches import Wedge, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

phylum_colors = {
    'Euryarchaeota': '#4B006E', 'Halobacterota': '#A77AC5', 'Thermoplasmatota': '#C9A0DC',
    'Actinobacteriota': '#A6E3E9', 'Armatimonadota': '#7FD8BE', 'Firmicutes': '#F4D03F',
    'Bacteroidota': '#00C853', 'Caldatribacteriota': '#B6F36C', 'Campylobacterota': '#9DE548',
    'Chloroflexi': '#C8F04D', 'Cloacimonadota': '#B5D334', 'Desulfobacterota': '#B4A800',
    'Planctomycetota': '#FF7F7F', 'Proteobacteria': '#FF5C5C', 'SAR324_cladeMarine_group_B': '#E08283',
    'Spirochaetota': '#FFB3B3', 'Synergistota': '#B22222', 'Thermotogota': '#8B0000',
    'WPS-2': '#D9AF6B', 'WS1': '#D62728', 'Unknown': '#BDBDBD'
}

SMALL_THRESHOLD = 0.5
MEDIUM_THRESHOLD = 5
LARGE_THRESHOLD = 15
SMALL_SIZE = 30
MEDIUM_SIZE = 90
LARGE_SIZE = 260

#####################################
# LOAD DATA
#####################################

nodes = pd.read_csv("/home/claude/network_analysis/nodes1.csv")
pos_edges_df = pd.read_csv("/home/claude/network_analysis/edges_new.csv")
neg_edges_df = pd.read_csv("/home/claude/network_analysis/edges_negative.csv")

# keep only negative edges between nodes that are actually in our node set/module structure
valid_nodes = set(nodes["Node"])
neg_edges_df = neg_edges_df[neg_edges_df["Source"].isin(valid_nodes) & neg_edges_df["Target"].isin(valid_nodes)]

module_ids = sorted(nodes["Module"].unique())
module_cmap = mpl.colormaps["tab20"].resampled(len(module_ids))
module_colors = {m: module_cmap(i) for i, m in enumerate(module_ids)}

#####################################
# ORDER NODES AROUND CIRCLE (clockwise, grouped by module)
#####################################

nodes_sorted = nodes.sort_values(["Module", "Relative_abundance"], ascending=[True, False]).reset_index(drop=True)

GAP_DEG = 3.0
n_modules = len(module_ids)
total_gap = GAP_DEG * n_modules
usable_deg = 360 - total_gap

module_counts = nodes_sorted.groupby("Module").size().to_dict()
total_n = len(nodes_sorted)
deg_per_node = {m: usable_deg * (module_counts[m] / total_n) for m in module_ids}

angles = {}
module_arc_range = {}
cursor = 90.0
for m in module_ids:
    span = deg_per_node[m]
    sub = nodes_sorted[nodes_sorted["Module"] == m]
    n_m = len(sub)
    step = span / n_m if n_m > 0 else 0
    start = cursor
    for i, (_, row) in enumerate(sub.iterrows()):
        angles[row["Node"]] = cursor - step * (i + 0.5)
    end = cursor - span
    module_arc_range[m] = (end, start)
    cursor = end - GAP_DEG

R = 1.0

def polar(angle_deg, r=R):
    a = np.radians(angle_deg)
    return r * np.cos(a), r * np.sin(a)

node_pos = {row["Node"]: polar(angles[row["Node"]]) for _, row in nodes_sorted.iterrows()}
node_info = nodes_sorted.set_index("Node").to_dict("index")

def node_size(abundance):
    if abundance < SMALL_THRESHOLD:
        return SMALL_SIZE
    elif abundance < LARGE_THRESHOLD:
        return MEDIUM_SIZE
    else:
        return LARGE_SIZE

#####################################
# DRAW
#####################################

fig, ax = plt.subplots(figsize=(15, 15), subplot_kw={"aspect": "equal"})

BAND_INNER = 1.04
BAND_OUTER = 1.14

for m, (start, end) in module_arc_range.items():
    wedge = Wedge(center=(0, 0), r=BAND_OUTER, theta1=start, theta2=end,
                   width=BAND_OUTER - BAND_INNER, facecolor=module_colors[m],
                   edgecolor="white", linewidth=1.0, zorder=2)
    ax.add_patch(wedge)
    mid = np.radians((start + end) / 2)
    lx, ly = (BAND_OUTER + 0.09) * np.cos(mid), (BAND_OUTER + 0.09) * np.sin(mid)
    rot = np.degrees(mid) % 360
    if 90 < rot < 270:
        rot += 180
    ax.text(lx, ly, f"Module {m}", ha="center", va="center", rotation=rot,
             rotation_mode="anchor", fontsize=13, fontweight="bold", color=module_colors[m],
             path_effects=[pe.withStroke(linewidth=3, foreground="white")])

#####################################
# POSITIVE EDGES (co-occurrence) -- green, as before
#####################################

pos_rhos = pos_edges_df["rho"].values
pos_min, pos_max = pos_rhos.min(), pos_rhos.max()
coolwarm = mpl.colormaps["coolwarm"]
pos_cmap = lambda t: coolwarm(0.35 - 0.35 * t)   # blue half of coolwarm; t=1 -> strongest blue

pos_edges_df = pos_edges_df.copy()
pos_edges_df["same_module"] = pos_edges_df.apply(
    lambda r: node_info.get(r["Source"], {}).get("Module") == node_info.get(r["Target"], {}).get("Module"), axis=1)
pos_edges_sorted = pos_edges_df.sort_values("same_module")

for _, row in pos_edges_sorted.iterrows():
    u, v, rho = row["Source"], row["Target"], row["rho"]
    if u not in node_pos or v not in node_pos:
        continue
    x1, y1 = node_pos[u]
    x2, y2 = node_pos[v]
    same_mod = row["same_module"]
    pull = 0.15 if same_mod else 0.85
    cx, cy = (x1 + x2) / 2 * (1 - pull), (y1 + y2) / 2 * (1 - pull)
    norm = (rho - pos_min) / (pos_max - pos_min) if pos_max > pos_min else 0.5
    color = pos_cmap(norm)
    width = 0.15 + 1.1 * rho if not same_mod else 0.1 + 0.6 * rho
    alpha = 0.45 if not same_mod else 0.22
    path = Path([(x1, y1), (cx, cy), (x2, y2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=width, alpha=alpha, zorder=2))

#####################################
# NEGATIVE EDGES (mutual exclusion) -- red/orange, drawn on top, pulled deep
# through the center so cross-module antagonism reads as bold sweeping chords.
#####################################

neg_rhos = neg_edges_df["rho"].abs().values
neg_min, neg_max = neg_rhos.min(), neg_rhos.max()
neg_cmap = lambda t: coolwarm(0.65 + 0.35 * t)   # red half of coolwarm; t=1 -> strongest red

for _, row in neg_edges_df.iterrows():
    u, v, rho = row["Source"], row["Target"], abs(row["rho"])
    if u not in node_pos or v not in node_pos:
        continue
    x1, y1 = node_pos[u]
    x2, y2 = node_pos[v]
    pull = 0.9
    cx, cy = (x1 + x2) / 2 * (1 - pull), (y1 + y2) / 2 * (1 - pull)
    norm = (rho - neg_min) / (neg_max - neg_min) if neg_max > neg_min else 0.5
    color = neg_cmap(norm)
    width = 0.2 + 1.3 * rho
    path = Path([(x1, y1), (cx, cy), (x2, y2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=width, alpha=0.35, zorder=1.5))

#####################################
# NODES
#####################################

for _, row in nodes_sorted.iterrows():
    node = row["Node"]
    x, y = node_pos[node]
    phylum = row["Phylum"] if pd.notna(row["Phylum"]) else "Unknown"
    fill = phylum_colors.get(phylum, "#BDBDBD")
    size = node_size(row["Relative_abundance"])
    role = row["Role"].strip() if isinstance(row["Role"], str) else ""

    is_hub = role in ("Provincial hub", "Connector hub", "Kinless hub")
    is_connector = role == "Non-hub connector"

    z = 5 if (is_hub or is_connector) else 3
    ax.scatter([x], [y], s=size, c=[fill], edgecolors="black", linewidths=0.6, zorder=z)

    if is_hub:
        ax.scatter([x], [y], s=size + 220, facecolors="none", edgecolors="black", linewidths=2.2, zorder=z - 1)
    elif is_connector:
        ax.scatter([x], [y], s=size + 160, facecolors="none", edgecolors="#1f6f8b",
                   linewidths=1.8, linestyle="dashed", zorder=z - 1)

    if is_hub or is_connector or role == "Dominant methanogen":
        ang = angles[node]
        lx, ly = polar(ang, r=BAND_OUTER + 0.16)
        rot = ang
        ha = "left" if -90 < ((ang + 180) % 360 - 180) < 90 else "right"
        if ha == "right":
            rot += 180
        ax.plot([x, (BAND_OUTER + 0.02) * np.cos(np.radians(ang))],
                 [y, (BAND_OUTER + 0.02) * np.sin(np.radians(ang))], color="black", lw=0.8, zorder=4)
        ax.text(lx, ly, node, ha=ha, va="center", rotation=rot, rotation_mode="anchor",
                 fontsize=11, fontstyle="italic",
                 path_effects=[pe.withStroke(linewidth=3, foreground="white")], zorder=6)

#####################################
# LEGEND
#####################################

legend_elements = [Patch(facecolor='none', edgecolor='none', label="NETWORK FEATURES")]
legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
                               markeredgecolor='black', markeredgewidth=2.0, markersize=14,
                               linestyle='None', label='Hub (z \u2265 2.5)'))
legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
                               markeredgecolor='#1f6f8b', markeredgewidth=1.6, markersize=13,
                               linestyle='None', label='Connector (non-hub, high P)'))
legend_elements.append(Line2D([0], [0], color=pos_cmap(1.0), lw=2.5, label='Positive edge'))
legend_elements.append(Line2D([0], [0], color=neg_cmap(1.0), lw=2.5, label='Negative edge'))

legend_elements.append(Patch(facecolor='none', edgecolor='none', label="RELATIVE ABUNDANCE (%)"))
legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                               markersize=np.sqrt(SMALL_SIZE), label=f'< {SMALL_THRESHOLD}'))
legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                               markersize=np.sqrt(MEDIUM_SIZE), label=f'{SMALL_THRESHOLD}-{LARGE_THRESHOLD}'))
legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                               markersize=np.sqrt(LARGE_SIZE), label=f'> {LARGE_THRESHOLD}'))

legend_elements.append(Patch(facecolor='none', edgecolor='none', label=" "))
legend_elements.append(Patch(facecolor='none', edgecolor='none', label="PHYLUM (node color)"))
used_phyla = set(nodes['Phylum'].dropna().unique()) | {"Unknown"}
for phylum, color in phylum_colors.items():
    if phylum not in used_phyla:
        continue
    legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                                   markeredgecolor='white', markersize=13, label=phylum))

leg = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1.0),
                 frameon=True, framealpha=0.95, fontsize=12)
for text in leg.get_texts():
    if text.get_text() in ["NETWORK FEATURES", "RELATIVE ABUNDANCE (%)", "PHYLUM (node color)"]:
        text.set_fontweight('bold')

#####################################
# FINAL
#####################################

lim = BAND_OUTER + 0.9
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.axis("off")
plt.tight_layout()
plt.savefig("/home/claude/network_analysis/cooccurrence_circos_signed.png", dpi=300, bbox_inches="tight")
print("Saved plot")
