#!/usr/bin/env python3
"""
Upper-bound analysis: the sensitivity ceiling under perfect NLP recall.

Reproduces Table 8 and the persistent-false-negative breakdown of the paper
(Section 3.3).

The ablation measures what the NLP classifiers currently contribute. This
measures what they could contribute at best. For each modality, every report
of that modality is treated as classifier-positive -- a classifier with 100%
recall -- and the pipeline is re-run.

Reports are added only for patients who have a gold episode the baseline
missed. Adding them for every patient would inflate the false-positive count
with episodes a perfect classifier would not be credited for, and could
re-date existing true positives and lose them. Restricting to currently
missed patients keeps the baseline true positives intact, so the result is a
true ceiling: no choice of classifier threshold can do better.

A gold episode still missed under this condition cannot be recovered by
improving NLP at all, and is attributed to one of three causes:

    no NLP evidence     the patient has no histopathology, CT or PET report
    episode date shift  the episode is anchored to an earlier test, placing
                        it outside the +/-30-day matching window
    imaging-only        the episode is imaging-only and lacks corroborating
                        structured evidence within +/-28 days

Usage:
    python -m experiments.upper_bound --data synthetic/out --outdir results/
"""

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.offline import (  # noqa: E402
    load_combine_tests, load_evaluate, load_test_type_map, wilson_ci,
)

HISTO, CT, PET = "9", "14", "15"

CONDITIONS = [
    ("baseline", "Baseline", []),
    ("ub_histo", "All Histo Positive", [HISTO]),
    ("ub_ct",    "All CT Positive",    [CT]),
    ("ub_pet",   "All PET Positive",   [PET]),
    ("ub_all",   "All NLP Positive",   [HISTO, CT, PET]),
]

# Order types that corroborate an imaging-only episode
# (ifi_pipeline/combine_tests.py: filter_episodes_nlp).
CORROBORATING = {"1", "2", "5", "6", "8", "11"}

MATCH_WINDOW_DAYS = 30
CORROBORATION_WINDOW_DAYS = 28


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        return list(r), r.fieldnames


def template_row(positives, order_type_id):
    """A representative positive row of this modality, to copy fields from.

    Added tests need the column values the pipeline expects. Everything
    patient-specific is overwritten by the caller; the template only supplies
    the modality-constant fields (an NLP-derived test identifies no organism,
    carries no specimen, and is scored from a report).
    """
    for r in positives:
        if str(r["t1.order_type_id"]) == order_type_id:
            return dict(r)
    return None


def make_upper_bound_input(data, work, target_ots, fn_patients, positives,
                           cols, all_tests):
    """Write tmp_merged2.csv with every target-modality report of the missed
    patients marked positive, alongside the unchanged baseline positives."""
    existing = {(str(r["t1.patient_id"]), str(r["t1.test_id"]))
                for r in positives}
    rows = [dict(r) for r in positives]
    added = 0

    templates = {ot: template_row(positives, ot) for ot in target_ots}
    for t in all_tests:
        ot, pid = str(t["order_type"]), str(t["patient_id"])
        if ot not in target_ots or pid not in fn_patients:
            continue
        if (pid, str(t["test_id"])) in existing:
            continue        # already positive; leave it exactly as it is
        tmpl = templates.get(ot)
        if tmpl is None:
            continue
        r = dict(tmpl)
        r["t1.test_id"] = t["test_id"]
        r["t1.patient_id"] = pid
        r["t1.collection_date"] = t["test_date"]
        r["t1.order_type_id"] = ot
        rows.append(r)
        added += 1

    with open(work / "tmp_merged2.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    shutil.copy(data / "all_db_tests.csv", work / "all_db_tests.csv")
    return added


def eval_condition(work, ev, tmap, orders, load_gold):
    pred = ev.date_filter(
        ev.read_dict("t1.patient_id", str(work / "epi_table.csv")), "%d-%m-%Y")
    (work / "eval").mkdir(parents=True, exist_ok=True)
    ev.eval(load_gold(), pred,
            ["species", "genus", "class_specific", "class_broad"],
            tmap, orders, evaldir=str(work / "eval"))

    def rows_of(name):
        p = work / "eval" / name
        return list(csv.DictReader(open(p, encoding="utf-8"))) \
            if p.exists() else []

    tp, fn, fp = rows_of("tp.csv"), rows_of("fn.csv"), rows_of("fp.csv")
    episodes = sum(len(v) for v in pred.values())
    return dict(tp=len(tp), fn=len(fn), fp=len(fp), episodes=episodes,
                fn_rows=fn, pred=pred)


def classify_persistent_fn(fn_rows, all_tests, positives):
    """Attribute each still-missed gold episode to one of the three causes."""
    by_patient = defaultdict(list)
    for t in all_tests:
        by_patient[str(t["patient_id"])].append(t)
    pos_by_patient = defaultdict(list)
    for r in positives:
        pos_by_patient[str(r["t1.patient_id"])].append(r)

    buckets = defaultdict(list)
    for row in fn_rows:
        pid = str(row.get("patient_id") or row.get("t1.patient_id") or "")
        gold_date = None
        for fmt in ("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d"):
            try:
                gold_date = datetime.strptime(row["system_diagnosis_date"], fmt)
                break
            except (ValueError, KeyError):
                continue

        nlp = [t for t in by_patient.get(pid, [])
               if str(t["order_type"]) in {HISTO, CT, PET}]
        if not nlp:
            buckets["no NLP evidence"].append(pid)
            continue

        # Is there corroborating structured evidence near the gold date?
        corroborated = False
        if gold_date:
            for t in by_patient.get(pid, []):
                if str(t["order_type"]) not in CORROBORATING or not t["test_date"]:
                    continue
                try:
                    d = datetime.strptime(t["test_date"], "%Y-%m-%d")
                except ValueError:
                    continue
                if abs((gold_date - d).days) <= CORROBORATION_WINDOW_DAYS:
                    corroborated = True
                    break

        # Does the patient have a positive test far from the gold date that
        # could have anchored the episode outside the matching window?
        anchored_early = False
        if gold_date:
            for r in pos_by_patient.get(pid, []):
                try:
                    d = datetime.strptime(r["t1.collection_date"], "%Y-%m-%d")
                except (ValueError, KeyError):
                    continue
                if abs((gold_date - d).days) > MATCH_WINDOW_DAYS:
                    anchored_early = True
                    break

        if anchored_early:
            buckets["episode date shift"].append(pid)
        elif not corroborated:
            buckets["imaging-only episode filtered"].append(pid)
        else:
            buckets["other"].append(pid)
    return buckets


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="directory of input CSVs")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    data, out = Path(args.data).resolve(), Path(args.outdir).resolve()
    work_root = out / "upper_bound"
    work_root.mkdir(parents=True, exist_ok=True)

    combine, ev = load_combine_tests(), load_evaluate()
    tmap, orders = load_test_type_map(data)

    def load_gold():
        return ev.date_filter(
            ev.read_dict("patient_id", str(data / "gold_standard.csv")),
            "%d/%m/%Y")

    n_gold = sum(len(v) for v in load_gold().values())
    positives, cols = read_rows(data / "tmp_merged2.csv")
    all_tests, _ = read_rows(data / "all_db_tests.csv")

    print(f"Gold standard: {n_gold} episodes\n")

    # Baseline first: its false negatives decide which patients get extra
    # reports in every other condition.
    base_work = work_root / "baseline"
    base_work.mkdir(parents=True, exist_ok=True)
    shutil.copy(data / "tmp_merged2.csv", base_work / "tmp_merged2.csv")
    shutil.copy(data / "all_db_tests.csv", base_work / "all_db_tests.csv")
    print("  Baseline ...", flush=True)
    combine(str(base_work))
    base = eval_condition(base_work, ev, tmap, orders, load_gold)

    fn_patients = set()
    for row in base["fn_rows"]:
        pid = row.get("patient_id") or row.get("t1.patient_id")
        if pid:
            fn_patients.add(str(pid))
    print(f"    {base['tp']} TP, {base['fn']} FN "
          f"({len(fn_patients)} patients with a missed episode)\n")

    results = [dict(key="baseline", label="Baseline", added=0, **{
        k: base[k] for k in ("tp", "fn", "fp", "episodes")})]

    for key, label, ots in CONDITIONS[1:]:
        work = work_root / key
        work.mkdir(parents=True, exist_ok=True)
        added = make_upper_bound_input(data, work, set(ots), fn_patients,
                                       positives, cols, all_tests)
        print(f"  {label} ... (+{added} reports treated as positive)", flush=True)
        combine(str(work))
        r = eval_condition(work, ev, tmap, orders, load_gold)
        results.append(dict(key=key, label=label, added=added, **{
            k: r[k] for k in ("tp", "fn", "fp", "episodes")}))
        if key == "ub_all":
            persistent = classify_persistent_fn(r["fn_rows"], all_tests,
                                                positives)

    for r in results:
        r["sensitivity"] = 100.0 * r["tp"] / n_gold if n_gold else 0.0
        det = r["tp"] + r["fp"]
        r["precision"] = 100.0 * r["tp"] / det if det else 0.0
        r["d_tp"] = r["tp"] - results[0]["tp"]
        lo, hi = wilson_ci(r["tp"], n_gold)
        r["sens_lo"], r["sens_hi"] = 100 * lo, 100 * hi
        # Paired Wilson interval on the change, as in Section 2.6.
        if r["d_tp"]:
            plo, phi = wilson_ci(abs(r["d_tp"]), n_gold)
            r["d_tp_lo"], r["d_tp_hi"] = 100 * plo, 100 * phi
        else:
            r["d_tp_lo"] = r["d_tp_hi"] = 0.0

    w = [20, 9, 5, 16, 4, 5, 22, 10]
    hdr = ["Condition", "Episodes", "TP", "dTP (95% CI)", "FN", "FP",
           "Sensitivity % (95% CI)", "Precision"]
    print("\n" + "  ".join(h.ljust(x) for h, x in zip(hdr, w)))
    print("  ".join("-" * x for x in w))
    for r in results:
        dtp = "---" if r["key"] == "baseline" else \
            f'+{r["d_tp"]} ({r["d_tp_lo"]:.1f}-{r["d_tp_hi"]:.1f}pp)'
        print("  ".join([
            r["label"].ljust(w[0]), str(r["episodes"]).rjust(w[1]),
            str(r["tp"]).rjust(w[2]), dtp.rjust(w[3]),
            str(r["fn"]).rjust(w[4]), str(r["fp"]).rjust(w[5]),
            f'{r["sensitivity"]:.2f} ({r["sens_lo"]:.1f}-{r["sens_hi"]:.1f})'.rjust(w[6]),
            f'{r["precision"]:.2f}'.rjust(w[7]),
        ]))

    ceiling = results[-1]
    print(f"\nSensitivity ceiling under perfect NLP recall: "
          f"{ceiling['tp']}/{n_gold} = {ceiling['sensitivity']:.2f}%")
    print(f"Recovered {ceiling['d_tp']} of the {results[0]['fn']} "
          f"baseline false negatives; {ceiling['fn']} remain.\n")

    print("Persistent false negatives by cause (Table 8)")
    total = 0
    for cause in ("no NLP evidence", "episode date shift",
                  "imaging-only episode filtered", "other"):
        pids = persistent.get(cause, [])
        if pids:
            print(f"  {cause:<32} {len(pids):>3}   patients {sorted(pids)}")
            total += len(pids)
    algorithmic = len(persistent.get("episode date shift", [])) + \
        len(persistent.get("imaging-only episode filtered", []))
    print(f"  {'':<32} {total:>3} total")
    print(f"\n  {algorithmic} of {total} are addressable by changing the "
          f"episode-detection rules rather than NLP.")

    cols_out = ["key", "label", "added", "episodes", "tp", "fn", "fp",
                "sensitivity", "sens_lo", "sens_hi", "precision",
                "d_tp", "d_tp_lo", "d_tp_hi"]
    out_csv = out / "upper_bound_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols_out, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(results)

    fn_csv = out / "persistent_false_negatives.csv"
    with open(fn_csv, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["cause", "patient_id"])
        for cause, pids in persistent.items():
            for pid in sorted(pids):
                wr.writerow([cause, pid])
    print(f"\nWrote {out_csv}\n      {fn_csv}")


if __name__ == "__main__":
    main()
