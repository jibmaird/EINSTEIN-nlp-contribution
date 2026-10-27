"""
Helpers for running the pipeline without a database.

In production, evaluate.py fetches the test-id -> order-type map over ODBC.
Every other code path in both algorithm modules is pure CSV, so with that one
map supplied from a file the whole evaluation runs offline and unmodified.
"""

import csv
import importlib.util
import math
import os
import sys
import types
from pathlib import Path

IFIS_DIR = Path(__file__).resolve().parents[1] / "ifi_pipeline"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_evaluate():
    """Import evaluate.py with pyodbc stubbed if it is not installed.

    evaluate.py imports pyodbc at module level but only read_test_type_map()
    uses it, and we replace that with load_test_type_map() below.
    """
    try:
        import pyodbc  # noqa: F401
    except ImportError:
        sys.modules["pyodbc"] = types.ModuleType("pyodbc")
    return _load("ifis_evaluate", IFIS_DIR / "evaluate.py")


def load_combine_tests():
    """Return a callable run(outdir) that executes episode detection.

    combine_tests.py reads <outdir>/tmp_merged2.csv and <outdir>/all_db_tests.csv
    and writes <outdir>/epi_table.csv. It touches no database: its --password,
    --database and --server arguments are accepted only for CLI compatibility.
    """
    mod = _load("ifis_combine_tests", IFIS_DIR / "combine_tests.py")

    def run(outdir):
        cwd = os.getcwd()
        try:
            mod.group_episode_org(outdir=outdir)
            mod.group_episode(outdir=outdir)
        finally:
            os.chdir(cwd)

    return run


def load_test_type_map(data_dir):
    """Rebuild read_test_type_map()'s (map, orders) from the two CSVs."""
    data_dir = Path(data_dir)
    ref, tst = data_dir / "ref_order_type.csv", data_dir / "test_order_type.csv"
    for p in (ref, tst):
        if not p.exists():
            sys.exit(f"Missing {p}")
    with open(ref, newline="", encoding="utf-8") as f:
        names = {r["order_type_id"]: r["order_type"] for r in csv.DictReader(f)}
    tmap, orders = {}, []
    with open(tst, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            name = names.get(r["order_type_id"])
            if name is None:
                continue
            tmap[r["test_id"]] = name
            if name not in orders:
                orders.append(name)
    return tmap, sorted(orders)


def wilson_ci(k, n, confidence=0.95):
    """Wilson score interval for a binomial proportion (Section 2.6)."""
    if n == 0:
        return 0.0, 0.0
    z = {0.90: 1.6449, 0.95: 1.9600, 0.99: 2.5758}.get(
        confidence, 1.9600)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, centre - half), min(1.0, centre + half)
