"""Regenerate the derived co-occurrence data files from the freshly rebuilt
network (network_edges.csv + network_metrics.csv/json), keep them consistent
with the new 3,276-ASV / 29-sample analysis, and sync into cooccurence/."""
import csv
import json
import shutil
from pathlib import Path

COOC = Path("cooccurence")
COOC.mkdir(exist_ok=True)

# ---- inputs (fresh, at repo root) ----
metrics = json.load(open("network_metrics.json"))
tax = json.load(open("model_inputs/taxonomy.json"))
itids = json.load(open("model_inputs/iterativeIDs.json"))          # seq -> ID
id2seq = {v: k for k, v in itids.items()}
ab = json.load(open("model_inputs/abundances.json"))              # sample -> seq -> frac
# per-node rows
node_rows = list(csv.DictReader(open("network_metrics.csv")))
edge_rows = list(csv.DictReader(open("network_edges.csv")))

# ---- node_connectivity.json  {node: degree}  (sorted desc) ----
node_deg = {r["node"]: int(r["degree"]) for r in node_rows}
node_deg = dict(sorted(node_deg.items(), key=lambda kv: kv[1], reverse=True))
json.dump(node_deg, open("node_connectivity.json", "w"), indent=2)

# ---- network_module_membership.json  from positive communities ----
comms = sorted(metrics["communities_positive"], key=len, reverse=True)
membership = {}
for i, members in enumerate(comms, 1):
    membership[f"module_{i}"] = {"size": len(members), "type": "louvain",
                                 "members": sorted(members)}
json.dump(membership, open("network_module_membership.json", "w"), indent=2)

# ---- Gephi nodes.csv  (Node,Phylum,Domain,Relative_abundance,Role,Module) ----
node_to_comm = {r["node"]: r["community"] for r in node_rows}
node_to_phylum = {r["node"]: r["phylum"] for r in node_rows}
# max rel abundance (%) across samples for each node's ASV
def max_ab(ID):
    seq = id2seq.get(ID)
    if seq is None:
        return ""
    return round(100 * max((ab[s].get(seq, 0) for s in ab), default=0), 4)
def domain(ID):
    seq = id2seq.get(ID)
    return tax.get(seq, {}).get("Kingdom", "") if seq else ""
with open("nodes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Node", "Phylum", "Domain", "Relative_abundance", "Role", "Module"])
    for r in node_rows:
        n = r["node"]
        w.writerow([n, node_to_phylum.get(n, ""), domain(n), max_ab(n), "", node_to_comm.get(n, "")])

# ---- Gephi edges.csv  (Source,Target,rho) ----
with open("edges.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Source", "Target", "rho"])
    for e in edge_rows:
        w.writerow([e["node_a"], e["node_b"], e["rho"]])

# ---- sync fresh files into cooccurence/ (the committed location) ----
for name in ["network_edges.csv", "network_metrics.csv", "network_metrics.json",
             "network_module_membership.json", "nodes.csv", "edges.csv"]:
    shutil.copy(name, COOC / name)

print("node_connectivity.json:", len(node_deg), "nodes")
print("network_module_membership.json:", len(membership), "modules; sizes",
      [membership[k]["size"] for k in membership])
print("nodes.csv rows:", len(node_rows), "| edges.csv rows:", len(edge_rows))
print("synced to cooccurence/:", ", ".join(sorted(p.name for p in COOC.glob("network_*")) +
      ["nodes.csv", "edges.csv"]))
