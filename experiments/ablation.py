#!/usr/bin/env python3
"""
Ablation study: remove each NLP test type, then all of them, and measure the
effect on episode detection against the gold standard.

Reproduces Tables 3, 4, 5 and 6 of the paper, and the super-additivity check
of Section 3.3.

For each condition the positive tests of the ablated modality are dropped from
the input, the unmodified episode-detection algorithm is re-run, and the
resulting episodes are matched against the gold standard. Nothing about the
algorithm or the matching changes between conditions -- only the input.

    order_type_id 9  = Histopathology
    order_type_id 14 = CT scan
    order_type_id 15 = PET scan

Usage:
    python -m src.experiments.ablation --data synthetic/out --outdir results/
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.offline import (  # noqa: E402
    load_combine_tests, load_evaluate, load_test_type_map, wilson_ci,
)

CONDITIONS = [
    ("baseline",  "Baseline",           []),
    ("no_histo",  "No Histopathology",  ["9"]),
    ("no_ct",     "No CT",              ["14"]),
    ("no_pet",    "No PET",             ["15"]),
    ("no_nlp",    "No NLP",             ["9", "14", "15"]),
]

MONTHS_IN_STUDY = 60


def filter_tests(src, dst, exclude_type_ids):
    """Write src -> dst, dropping positive tests of the excluded modalities."""
    kept = dropped = 0
    with open(src, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows, cols = [], reader.fieldnames
        for r in reader:
            if str(r["t1.order_type_id"]) in set(exclude_type_ids):
                dropped += 1
            else:
                rows.append(r)
                kept += 1
    with open(dst, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return kept, dropped


def run_condition(key, exclude, data, workroot, combine, ev, tmap, orders, load_gold):
    work = workroot / key
    (work / "eval").mkdir(parents=True, exist_ok=True)

    kept, dropped = filter_tests(data / "tmp_merged2.csv",
                                 work / "tmp_merged2.csv", exclude)
    shutil.copy(data / "all_db_tests.csv", work / "all_db_tests.csv")

    combine(str(work))

    pred = ev.date_filter(
        ev.read_dict("t1.patient_id", str(work / "epi_table.csv")), "%d-%m-%Y")
    # eval() normalises dates on the rows it is given, so each condition gets
    # a freshly parsed copy of the gold standard rather than a shared one.
    ev.eval(load_gold(), pred, ["species", "genus", "class_specific", "class_broad"],
            tmap, orders, evaldir=str(work / "eval"))

    def count(name):
        p = work / "eval" / name
        return sum(1 for _ in csv.DictReader(open(p, encoding="utf-8"))) \
            if p.exists() else 0

    tp, fn, fp = count("tp.csv"), count("fn.csv"), count("fp.csv")
    episodes = sum(len(v) for v in pred.values())
    return dict(key=key, excluded=dropped, input=kept, episodes=episodes,
                tp=tp, fn=fn, fp=fp)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True,
                    help="directory holding tmp_merged2.csv, all_db_tests.csv, "
                         "gold_standard.csv, ref_order_type.csv, test_order_type.csv")
    ap.add_argument("--outdir", default="results",
                    help="where to write ablation_results.csv and per-condition output")
    args = ap.parse_args()

    data = Path(args.data).resolve()
    out = Path(args.outdir).resolve()
    work = out / "conditions"
    work.mkdir(parents=True, exist_ok=True)

    combine = load_combine_tests()
    ev = load_evaluate()
    tmap, orders = load_test_type_map(data)
    def load_gold():
        return ev.date_filter(
            ev.read_dict("patient_id", str(data / "gold_standard.csv")),
            "%d/%m/%Y")

    n_gold = sum(len(v) for v in load_gold().values())

    print(f"Gold standard: {n_gold} episodes\n")
    rows = []
    for key, label, exclude in CONDITIONS:
        print(f"  {label} ...", flush=True)
        r = run_condition(key, exclude, data, work, combine, ev, tmap, orders,
                          load_gold)
        r["label"] = label
        rows.append(r)

    base = rows[0]
    for r in rows:
        r["sensitivity"] = 100.0 * r["tp"] / n_gold if n_gold else 0.0
        det = r["tp"] + r["fp"]
        r["precision"] = 100.0 * r["tp"] / det if det else 0.0
        s, p = r["sensitivity"] / 100, r["precision"] / 100
        r["f1"] = 2 * s * p / (s + p) if (s + p) else 0.0
        r["d_sens"] = r["sensitivity"] - base["sensitivity"]
        r["d_episodes"] = r["episodes"] - base["episodes"]
        r["per_month"] = r["episodes"] / MONTHS_IN_STUDY
        lo, hi = wilson_ci(r["tp"], n_gold)
        r["sens_lo"], r["sens_hi"] = 100 * lo, 100 * hi

    w = [20, 9, 7, 9, 5, 4, 5, 22, 10, 7]
    hdr = ["Condition", "Excluded", "Input", "Episodes", "TP", "FN", "FP",
           "Sensitivity % (95% CI)", "Precision", "F1"]
    print("\n" + "  ".join(h.ljust(x) for h, x in zip(hdr, w)))
    print("  ".join("-" * x for x in w))
    for r in rows:
        print("  ".join([
            r["label"].ljust(w[0]), str(r["excluded"]).rjust(w[1]),
            str(r["input"]).rjust(w[2]), str(r["episodes"]).rjust(w[3]),
            str(r["tp"]).rjust(w[4]), str(r["fn"]).rjust(w[5]),
            str(r["fp"]).rjust(w[6]),
            f'{r["sensitivity"]:.2f} ({r["sens_lo"]:.1f}-{r["sens_hi"]:.1f})'.rjust(w[7]),
            f'{r["precision"]:.2f}'.rjust(w[8]), f'{r["f1"]:.4f}'.rjust(w[9]),
        ]))

    by = {r["key"]: r for r in rows}
    ind = sum(-by[k]["d_sens"] for k in ("no_ct", "no_pet", "no_histo"))
    comb = -by["no_nlp"]["d_sens"]
    print(f"\nSuper-additivity check (Section 3.3)")
    print(f"  drop(no NLP)                              = {comb:.2f} pp")
    print(f"  drop(no CT) + drop(no PET) + drop(no Histo) = {ind:.2f} pp")
    print(f"  -> combined drop {'exceeds' if comb > ind else 'does NOT exceed'} "
          f"the sum of individual drops"
          f"{f', by {comb - ind:.2f} pp' if comb > ind else ''}")
    redundant = by["no_nlp"]["tp"] and (
        base["tp"] - by["no_nlp"]["tp"]) - sum(
        base["tp"] - by[k]["tp"] for k in ("no_ct", "no_pet", "no_histo"))
    print(f"  episodes supported by more than one modality: {redundant}")

    print(f"\nDiagnostic yield (TP lost per 1,000 positive tests excluded)")
    for k, name in (("no_ct", "CT"), ("no_pet", "PET"),
                    ("no_histo", "Histopathology"), ("no_nlp", "Combined NLP")):
        r, lost = by[k], base["tp"] - by[k]["tp"]
        y = 1000.0 * lost / r["excluded"] if r["excluded"] else 0.0
        print(f"  {name:<16} {r['excluded']:>6} tests  {lost:>3} TP lost  "
              f"yield {y:6.2f}")

    cols = ["key", "label", "excluded", "input", "episodes", "per_month",
            "tp", "fn", "fp", "sensitivity", "sens_lo", "sens_hi",
            "precision", "f1", "d_sens", "d_episodes"]
    out_csv = out / "ablation_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    print(f"\nWrote {out_csv}")


if __name__ == "__main__":
    main()
