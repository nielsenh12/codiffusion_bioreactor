# Summary Report: Louvain Module Detection on the Positive Co-occurrence Network

## Objective

Identify community structure ("modules") within a SparCC-derived microbial
co-occurrence network, to support downstream topological-role analysis
(z-P / Guimerà-Amaral classification) and visualization of the network.

## Data source

`sparcc_fdr_significant_pairs_r03.csv` — pairwise SparCC correlations between
ASVs, restricted to pairs passing FDR-corrected significance filtering
(2,218 total pairs: 1,269 positive, 949 negative).

## Methodological decision: why positive-only

The significant-pairs table mixes co-occurrence (positive correlation) and
mutual-exclusion (negative correlation) relationships. Before building a
network for community detection, we had to decide how to handle the negative
correlations, since they carry a different biological meaning (organisms
that tend *not* to co-occur) than positive ones.

We considered three options: using positive correlations only, using the
absolute value of all correlations, or an unweighted network ignoring sign
and magnitude. **We chose positive-only**, for two reasons:

1. **Algorithmic**: standard modularity — and therefore the Louvain
   algorithm — is built around comparing observed edge density within a
   group to the density expected under random rewiring. That comparison
   assumes non-negative weights; mixing in negative weights (or using their
   absolute value) changes what "community structure" even means and can
   produce misleading partitions if used naively with an unsigned modularity
   formula.
2. **Convention**: this is also standard practice in the 16S/microbiome
   co-occurrence network literature. Purpose-built tools such as SCNIC
   default to Louvain modularity maximization on a positive-only (or
   magnitude-thresholded) network; where negative correlations are analyzed
   at all, they are typically treated as a separate, complementary question
   (e.g., "co-exclusion" or "structural balance" analysis) rather than
   combined into a single signed community-detection score.

We did explore an alternative — implementing signed modularity (Gómez,
Jensen & Arenas, 2009) from scratch, since no available Python library
supported it — and validated that implementation against this positive-only
result (it reproduces the same partition when negative edges are excluded).
However, we ultimately kept the positive-only Louvain modules documented
here as the primary module assignment, with negative-edge structure examined
separately (see the companion exclusion-module analysis) rather than folded
into this clustering.

## Network construction

Only positive-correlation pairs (1,269 of 2,218) were used. This produced an
undirected, weighted graph of **266 ASVs** connected by **1,269 edges**,
where edge weight = the SparCC correlation coefficient (all positive by
construction). The resulting network has two connected components: one
giant component of 264 nodes, and one isolated pair.

## Method

Community detection was performed with the **Louvain algorithm** (Blondel et
al., 2008), using networkx's built-in `louvain_communities` implementation
(networkx 3.6.1). Louvain optimizes modularity — a measure of how much more
densely connected a proposed group of nodes is compared to what would be
expected by chance — via an iterative two-phase process: nodes are first
greedily reassigned to whichever neighboring group most increases
modularity, then groups are collapsed into single "super-nodes" and the
process repeats on the coarser graph. This continues until no further
reassignment improves modularity.

**Parameters used:**
- `weight = 'weight'` — the SparCC correlation coefficient was used directly
  as edge weight, so more strongly correlated pairs exert proportionally
  more influence on the grouping.
- `resolution = 1.0` — the default resolution; no attempt was made to tune
  this toward coarser or finer modules.
- `seed = 42` — fixed for reproducibility. Louvain's local-moving phase has
  some randomness in the order nodes are visited, which can change the exact
  partition found (though not usually the overall quality) between runs
  without a fixed seed.

## Results

The algorithm converged on **11 modules**, with an overall **modularity of
Q = 0.6065**. As a rule of thumb, modularity above ~0.4 is considered
evidence of real, non-random community structure — 0.6065 indicates the
positive co-occurrence network has a clean, well-separated modular
structure.

**Module sizes:** 63, 57, 52, 39, 25, 8, 8, 7, 3, 2, 2 (266 ASVs total).

The four largest modules (63, 57, 52, and 39 ASVs) account for 79% of all
ASVs in the network and represent the primary community structure. The three
smallest modules (3, 2, and 2 ASVs) are small, loosely connected fragments
rather than substantial communities in their own right, and should be
interpreted cautiously in any downstream analysis that assumes modules
reflect meaningful ecological groupings.

## Output

`asv_community_assignments.csv` — the per-ASV module assignment (module IDs
0–10, ordered from largest to smallest), along with each ASV's iterative ID
and unweighted degree in the positive network. This assignment was used
as-is in all subsequent analyses in this project: the z-P topological role
classification, the node/edge tables built for network visualization, and
the module-level mutual exclusion analysis.

## Limitations and notes for future work

- Modularity optimization via Louvain is a greedy heuristic, not an exact
  solver — it finds a good partition, not provably the best one. Different
  seeds can yield slightly different module boundaries, though the overall
  structure (module count, approximate sizes, modularity value) was
  consistent across the checks performed.
- No resolution sweep was performed; a different `resolution` parameter
  could merge or split these 11 modules differently. This may be worth
  exploring if a coarser or finer module structure is desired.
- This module assignment reflects positive co-occurrence only. It does not
  by itself indicate anything about mutual exclusion between ASVs or
  modules — that relationship was examined separately and is documented in
  the companion exclusion-module analysis.
