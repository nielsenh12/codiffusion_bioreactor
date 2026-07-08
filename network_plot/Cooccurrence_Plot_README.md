# Co-occurrence Network Plot — README

This documents `plot_circos_signed.py`, the circos-style visualization of the
positive co-occurrence network and its Louvain modules, with negative
(mutual exclusion) edges overlaid.

## Inputs

| File | Role |
|---|---|
| `nodes1.csv` | One row per ASV: `Node` (display name), `Phylum`, `Domain`, `Relative_abundance`, `Role` (z-P category or "Dominant methanogen"), `Module` (Louvain module 0–10). See `Node_Edge_Table_Workflow.md` for how this was built. |
| `edges_new.csv` | Positive-correlation edges: `Source, Target, rho`. 1,269 rows. |
| `edges_negative.csv` | Negative-correlation edges: `Source, Target, rho`. 920 rows (of 949 total negative significant pairs; 29 involved ASVs outside the positive-only node set and are excluded). |

## Why a circos layout

Earlier iterations of this figure used a force-directed (spring) layout,
which struggles at this network's size (266 nodes) — modules either had to
be forced apart on a grid (clean but arbitrary) or left to a physics
simulation (organic but prone to messy overlap between modules with heavy
cross-connectivity). A circos layout sidesteps this entirely: module
membership is encoded by position (a contiguous arc), not by anything the
layout algorithm has to "discover" or negotiate, so it stays legible
regardless of how densely connected modules are to each other. It also
scales better to overlaying a second edge type (negative correlations)
without the two edge sets fighting over layout space.

## Encodings

| Visual channel | Encodes |
|---|---|
| Position around the circle | Louvain module (contiguous arc per module, modules ordered clockwise starting at 12 o'clock) |
| Outer colored band (ideogram) | Module identity (color from `tab20`, one color per module) |
| Node fill color | Phylum |
| Node size | Relative abundance (excluding Day 0), in 3 tiers: <0.5%, 0.5–15%, >15% |
| Node black outline | All nodes (uniform) |
| Solid black ring around a node | Hub (z ≥ 2.5: Provincial/Connector/Kinless hub) |
| Dashed teal ring around a node | Non-hub connector (z < 2.5, high participation coefficient) |
| Blue chords | Positive correlation (co-occurrence); intensity scales with correlation strength |
| Red chords | Negative correlation (mutual exclusion); intensity scales with correlation strength |
| Italic name label + leader line | Hub, connector, or "Dominant methanogen" nodes only (not all 266 — labeling every node would be unreadable) |

Node ordering within each module's arc is by descending relative abundance,
so the largest/most abundant members of a module sit closest to the "front"
(the module's starting edge going clockwise).

## Chord (edge) styling

- **Positive edges**: drawn as quadratic Bézier curves. Within-module edges
  are pulled only slightly toward the center (`pull = 0.15`), so they read
  as short, faint arcs hugging the rim — visually confirming local density
  without cluttering the middle. Cross-module edges are pulled deep toward
  the center (`pull = 0.85`), so any module-to-module bridging shows up as a
  bold sweeping chord.
- **Negative edges**: always pulled deep toward center (`pull = 0.9`),
  regardless of module, since the negative-edge layer exists specifically to
  surface cross-module exclusion patterns (see the exclusion-module
  analysis, where Module 0 ↔ Module 3 stood out sharply).
- **Draw order**: positive edges render on top of negative edges (`zorder=2`
  vs `zorder=1.5`), so the blue co-occurrence layer isn't visually buried
  under the (typically denser) red exclusion layer.
- **Color**: both edge types sample from the `coolwarm` diverging colormap
  — positive edges from its blue half, negative edges from its red half —
  rather than independent colormaps, so the two are visually paired as two
  ends of one spectrum. Labeled simply "Positive edge" / "Negative edge" in
  the legend, without implying a specific causal mechanism (co-occurrence
  vs. exclusion) behind either.

## Layout parameters

| Parameter | Value | Meaning |
|---|---|---|
| `GAP_DEG` | 3.0° | Angular gap between adjacent module arcs |
| `BAND_INNER` / `BAND_OUTER` | 1.04 / 1.14 | Inner/outer radius of the module ideogram band |
| Node label radius | `BAND_OUTER + 0.16` | Distance of hub/connector/methanogen name labels from the circle |
| `SMALL_SIZE` / `MEDIUM_SIZE` / `LARGE_SIZE` | 30 / 90 / 260 | Marker sizes (matplotlib `s` units) for the 3 abundance tiers |
| Figure size / resolution | 15×15 in, 300 dpi | — |

## Known limitations / things to watch

- With 920 negative edges all drawn through the center, dense regions can
  visually compete with each other — the intended Module 0↔3 signal should
  still be the densest red bundle, but weaker pairwise negative
  relationships elsewhere may be harder to individually distinguish.
- Module label rotation is normalized (`% 360`) before deciding whether to
  flip the text right-side-up; if module count or ordering changes
  significantly, double-check labels aren't upside-down on the left half of
  the circle.
- Only hub/connector/"Dominant methanogen" nodes are labeled by name — this
  is a deliberate legibility choice, not a limitation of the data (all role
  and taxonomy information is still present in `nodes1.csv` for any node,
  labeled or not).

## Reproducing

```bash
python3 plot_circos_signed.py
```

Reads `nodes1.csv`, `edges_new.csv`, and `edges_negative.csv` from the
working directory; writes `cooccurrence_circos_signed.png`.
