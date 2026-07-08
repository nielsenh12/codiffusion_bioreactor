"""Re-execute specific self-contained cells of _regen_newdata.ipynb (each loads
its inputs from disk), used to regenerate the correlation figures that failed in
the full run (pandas-2.x groupby fix + cleared stale correlation files)."""
import sys
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TARGETS = [int(x) for x in sys.argv[1:]] or [17, 18, 19, 23]

nb = json.load(open("_regen_newdata.ipynb"))
cells = [c for c in nb["cells"] if c["cell_type"] == "code"]


def _display(*a, **k):
    pass


for i in TARGETS:
    src = "".join(cells[i]["source"])
    ns = {"__name__": "__main__", "display": _display, "get_ipython": lambda: None}
    print(f"\n===== running cell {i} =====")
    try:
        exec(compile(src, f"<cell{i}>", "exec"), ns)
        print(f"cell {i}: OK")
    except Exception as e:
        import traceback
        print(f"cell {i}: FAILED -> {type(e).__name__}: {e}")
        traceback.print_exc()
    finally:
        plt.close("all")
