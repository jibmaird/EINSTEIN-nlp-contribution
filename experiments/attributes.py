#!/usr/bin/env python3
"""
Contribution of NLP to episode classification and diagnosis date.

Reproduces Table 9, Figure 2 and Section 3.4 of the paper.

The ablation and upper-bound analyses measure whether an episode is detected.
The system also assigns each detected episode an EORTC/MSG classification and
a diagnosis date, and NLP evidence can change those independently of whether
the episode is found.

To isolate that, every comparison here is restricted to gold episodes detected
under *both* the baseline and the ablated condition. A difference therefore
reflects the attribute assigned, not which episodes were found. Denominators
differ between conditions because the common set does.

Two mechanisms in ifi_pipeline/combine_tests.py drive the classification
result, and both are worth knowing when reading the output:

  - Imaging tests are skipped when a classification is assigned, so an
    imaging-only episode falls back to Possible.
  - Possible/Pending is upgraded to Probable only when the episode holds a CT
    (order type 14) AND a non-imaging supporting test (type 2, 7, 8, 9 or 11).

Probable is therefore unreachable without CT, which is why removing CT alone
collapses Probable agreement rather than merely reducing it.

Organism and anatomical site are not reported. NLP-derived tests identify no
organism, and site is carried from the structured tests, so neither attribute
is affected by these ablations.

Usage:
    python -m experiments.attributes --data synthetic/out --outdir results/
"""

import argparse
import csv
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.offline import (  # noqa: E402
    load_combine_tests, load_evaluate, load_test_type_map,
)

CONDITIONS = [
    ("baseline", "Baseline",          []),
    ("no_histo", "No Histopathology", ["9"]),
    ("no_ct",    "No CT",             ["14"]),
    ("no_pet",   "No PET",            ["15"]),
    ("no_nlp",   "No NLP",            ["9", "14", "15"]),
]

CLASS_ORDER = ["Proven", "Probable", "Possible", "Pending"]
DATE_FORMATS = ("%d-%m-%Y", "%d/%m/%Y", "%Y-%m-%d")


def parse_date(s):
    for fmt in DATE_FORMATS:
        try:
            return datetime.strptime(s, fmt)
        except (ValueError, TypeError):
            continue
    return None


def filter_tests(src, dst, exclude):
    with open(src, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        cols = r.fieldnames
        rows = [x for x in r if str(x["t1.order_type_id"]) not in set(exclude)]
    with open(dst, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def run_condition(key, exclude, data, work_root, combine, ev, tmap, orders,
                  load_gold):
    """Run one condition and return its true positives keyed by gold episode."""
    work = work_root / key
    (work / "eval").mkdir(parents=True, exist_ok=True)
    filter_tests(data / "tmp_merged2.csv", work / "tmp_merged2.csv", exclude)
    shutil.copy(data / "all_db_tests.csv", work / "all_db_tests.csv")
    combine(str(work))
    pred = ev.date_filter(
        ev.read_dict("t1.patient_id", str(work / "epi_table.csv")), "%d-%m-%Y")
    ev.eval(load_gold(), pred,
            ["species", "genus", "class_specific", "class_broad"],
            tmap, orders, evaldir=str(work / "eval"))

    tp_path = work / "eval" / "tp.csv"
    tps = {}
    if tp_path.exists():
        with open(tp_path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                gid = (str(row.get("patient_id", "")),
                       str(row.get("ifi_episode_id", "")))
                tps[gid] = row
    return tps


def agreement(rows):
    """Overall and per-gold-category classification agreement."""
    same = 0
    by_cat = defaultdict(lambda: [0, 0])     # gold category -> [agree, total]
    for r in rows:
        g = (r.get("ifi_classification") or "").strip()
        p = (r.get("pred_ifi_classification") or "").strip()
        ok = (g == p)
        same += ok
        by_cat[g][0] += ok
        by_cat[g][1] += 1
    return same, len(rows), by_cat


def confusion(rows):
    c = Counter()
    for r in rows:
        g = (r.get("ifi_classification") or "").strip() or "(none)"
        p = (r.get("pred_ifi_classification") or "").strip() or "(none)"
        c[(g, p)] += 1
    return c


def date_stats(rows):
    """Signed offset (system date - gold date), in days.

    A negative offset means the system dated the episode before the gold
    standard. Because a matched episode is within the +/-30-day window by
    construction, an earlier date is lead time rather than error: the NLP
    evidence anchoring the episode precedes microbiological confirmation.
    """
    offsets = []
    for r in rows:
        g = parse_date(r.get("system_diagnosis_date", ""))
        p = parse_date(r.get("pred_system_diagnosis_date", ""))
        if g and p:
            offsets.append((p - g).days)
    earlier = [o for o in offsets if o < 0]
    return dict(
        n=len(offsets),
        mean=sum(offsets) / len(offsets) if offsets else 0.0,
        n_earlier=len(earlier),
        pct_earlier=100.0 * len(earlier) / len(offsets) if offsets else 0.0,
        mean_lead=-sum(earlier) / len(earlier) if earlier else 0.0,
    )


def print_confusion(title, c, labels):
    print(f"\n  {title}  (n = {sum(c.values())})")
    head = "    gold \\ system   " + "".join(f"{l:>10}" for l in labels)
    print(head)
    print("    " + "-" * (len(head) - 4))
    for g in labels:
        row = "".join(f"{c.get((g, p), 0):>10}" for p in labels)
        print(f"    {g:<16}{row}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="directory of input CSVs")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    data, out = Path(args.data).resolve(), Path(args.outdir).resolve()
    work_root = out / "attributes"
    work_root.mkdir(parents=True, exist_ok=True)

    combine, ev = load_combine_tests(), load_evaluate()
    tmap, orders = load_test_type_map(data)

    def load_gold():
        return ev.date_filter(
            ev.read_dict("patient_id", str(data / "gold_standard.csv")),
            "%d/%m/%Y")

    tps = {}
    for key, label, exclude in CONDITIONS:
        print(f"  {label} ...", flush=True)
        tps[key] = run_condition(key, exclude, data, work_root, combine, ev,
                                 tmap, orders, load_gold)

    base_ids = set(tps["baseline"])

    # ---- Table 9: agreement by gold category -------------------------
    print("\n\nEORTC classification agreement, restricted to episodes "
          "detected under both conditions")
    labels = [k for k, _, _ in CONDITIONS]
    names = {k: l for k, l, _ in CONDITIONS}
    common, summary = {}, {}
    for key in labels:
        ids = base_ids & set(tps[key])
        common[key] = ids
        rows = [tps[key][i] for i in sorted(ids)]
        same, total, by_cat = agreement(rows)
        summary[key] = dict(same=same, total=total, by_cat=by_cat,
                            pct=100.0 * same / total if total else 0.0)

    w = 20
    print("\n  " + "Gold category".ljust(16)
          + "".join(names[k].rjust(w) for k in labels))
    print("  " + "-" * (16 + w * len(labels)))
    cats = [c for c in CLASS_ORDER
            if any(c in summary[k]["by_cat"] for k in labels)]
    for cat in cats:
        cells = []
        for k in labels:
            a, t = summary[k]["by_cat"].get(cat, [0, 0])
            cells.append((f"{a}/{t} ({100.0*a/t:.0f}%)" if t else "-").rjust(w))
        print(f"  {cat:<16}" + "".join(cells))
    print("  " + "-" * (16 + w * len(labels)))
    print("  " + "Overall".ljust(16) + "".join(
        f'{summary[k]["same"]}/{summary[k]["total"]} '
        f'({summary[k]["pct"]:.0f}%)'.rjust(w) for k in labels))
    print("\n  Note: each column is restricted to gold episodes detected in "
          "both the\n  baseline and that condition, so denominators differ "
          "between columns.")

    # ---- Figure 2: confusion matrices on one common set ---------------
    shared = base_ids & set(tps["no_nlp"])
    print(f"\n\nClassification assigned, baseline vs no NLP")
    print(f"  Both panels use the same {len(shared)} episodes, detected "
          f"under both conditions.")
    cb = confusion([tps["baseline"][i] for i in sorted(shared)])
    cn = confusion([tps["no_nlp"][i] for i in sorted(shared)])
    seen = {g for g, _ in list(cb) + list(cn)} | \
           {p for _, p in list(cb) + list(cn)}
    lab = [c for c in CLASS_ORDER if c in seen] + \
          sorted(seen - set(CLASS_ORDER))
    print_confusion("Baseline (all NLP)", cb, lab)
    print_confusion("No NLP", cn, lab)

    changed = [(i, tps["baseline"][i].get("pred_ifi_classification"),
                tps["no_nlp"][i].get("pred_ifi_classification"))
               for i in sorted(shared)
               if tps["baseline"][i].get("pred_ifi_classification")
               != tps["no_nlp"][i].get("pred_ifi_classification")]
    print(f"\n  {len(changed)} of {len(shared)} episodes change "
          f"classification when NLP is removed:")
    for (before, after), n in Counter(
            (c[1], c[2]) for c in changed).most_common():
        print(f"    {n:>3}  {before} -> {after}")
    assignable = {k: sorted({r.get("pred_ifi_classification")
                             for r in tps[k].values()} - {None, ""})
                  for k in ("baseline", "no_ct", "no_nlp")}
    for k in ("baseline", "no_ct", "no_nlp"):
        print(f"\n  categories the system can assign, {names[k]:<18}: "
              f"{', '.join(assignable[k])}")

    # ---- Diagnosis date ----------------------------------------------
    print("\n\nDiagnosis date, relative to the gold standard")
    print("\n  " + "Condition".ljust(20) + "Matched".rjust(9)
          + "Dated earlier".rjust(16) + "Mean lead (d)".rjust(15)
          + "Mean offset (d)".rjust(17))
    print("  " + "-" * 77)
    dates = {}
    for key in labels:
        ids = common[key]
        st = date_stats([tps[key][i] for i in sorted(ids)])
        dates[key] = st
        print("  " + names[key].ljust(20) + str(st["n"]).rjust(9)
              + f'{st["n_earlier"]} ({st["pct_earlier"]:.0f}%)'.rjust(16)
              + f'{st["mean_lead"]:.1f}'.rjust(15)
              + f'{st["mean"]:+.2f}'.rjust(17))
    print("\n  A negative offset means the system dated the episode before "
          "the gold\n  standard. Such an episode is inside the +/-30-day "
          "matching window by\n  construction, so an earlier date is lead "
          "time, not error.")

    # ---- Output files -------------------------------------------------
    f1 = out / "classification_agreement.csv"
    with open(f1, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["condition", "gold_category", "agree", "total", "pct"])
        for k in labels:
            for cat in cats:
                a, t = summary[k]["by_cat"].get(cat, [0, 0])
                wr.writerow([k, cat, a, t,
                             f"{100.0*a/t:.1f}" if t else ""])
            wr.writerow([k, "OVERALL", summary[k]["same"],
                         summary[k]["total"], f'{summary[k]["pct"]:.1f}'])

    f2 = out / "classification_confusion.csv"
    with open(f2, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["condition", "gold_classification",
                     "system_classification", "n"])
        for name, c in (("baseline", cb), ("no_nlp", cn)):
            for (g, p), n in sorted(c.items()):
                wr.writerow([name, g, p, n])

    f3 = out / "diagnosis_date.csv"
    with open(f3, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["condition", "matched", "n_earlier", "pct_earlier",
                     "mean_lead_days", "mean_signed_offset_days"])
        for k in labels:
            st = dates[k]
            wr.writerow([k, st["n"], st["n_earlier"],
                         f'{st["pct_earlier"]:.1f}',
                         f'{st["mean_lead"]:.2f}', f'{st["mean"]:.2f}'])

    print(f"\nWrote {f1}\n      {f2}\n      {f3}")


if __name__ == "__main__":
    main()
