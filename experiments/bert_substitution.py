"""
full_substitution.py
---------------------
Offline BERT-vs-production substitution that can ADD tests, not only remove
them.

build_histo_ct_pet_csv.py can only keep rows already present in the reference
tmp_merged2.csv (production positives), so a BERT classifier could reject
production positives but never contribute its own extra positives. This
script does both. For each swapped modality it
  1. drops reference rows of that order type whose test is not BERT-positive,
  2. adds a row for every BERT-positive test that was not production-positive.

Added rows need patient_id and collection_date. Histopathology tests get them
from all_db_tests.csv; CT/PET tests get them from the export written by
export_ct_pet_test_dates.py (all_db_tests.csv has no CT/PET rows). Organism
fields are constant for histo/CT/PET tests ("Unspecified organism") and are
copied from a template row of the same order type. Fields that do not affect
episode grouping or gold-standard matching (site, order name, encounter,
specimen classification) are copied from that template too.

Known approximation: production histo rows classified "Colonisation" (urine
specimens, ~3.5% of histo positives) are skipped by combine_tests.py; an added
histo row inherits the template's classification, so a few added urine
specimens would be included that the real pipeline would skip.

Usage (from the IFI_algorithm directory):
    python Experiments/Bert_runs/full_substitution.py \\
        --ref-merged Data/runs/Histo_Bert/tmp_merged.csv \\
        --all-db-tests Data/runs/Histo_Bert/all_db_tests.csv \\
        --dates Experiments/NLP_contribution/tmp.csv \\
        --swap pet=Data/runs/Histo_Bert/pet_scans_processed.csv,Data/runs/Histo_Bert/pet_scans_processed_bert.csv \\
        --outdir /tmp/pet_full

    --swap is repeatable (histo=..., ct=..., pet=...) as NAME=PROD_CSV,BERT_CSV.
    --veto-only reproduces the earlier reject-only design for comparison.
    With no --swap the reference is rebuilt unchanged (baseline).
"""

import argparse
import importlib.util
import os
import shutil
import subprocess
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
IFI_DIR = os.path.normpath(os.path.join(_HERE, "..", ".."))
ORDER_TYPE = {"histo": "9", "ct": "14", "pet": "15"}


def positive_ids(csv_path):
    d = pd.read_csv(csv_path, usecols=["test_id", "y_pred"], dtype={"test_id": str})
    return set(d[d.y_pred == 1].test_id)


def load_dates(all_db_tests, dates_export):
    a = pd.read_csv(all_db_tests, dtype=str).dropna(subset=["order_type"])
    a = a.drop_duplicates("test_id").set_index("test_id")[["patient_id", "test_date"]]
    a.columns = ["patient_id", "date"]
    parts = [a]
    if dates_export:
        e = pd.read_csv(dates_export, sep=None, engine="python", dtype=str)
        e["date"] = pd.to_datetime(e.collection_date, format="%m/%d/%Y", errors="coerce").dt.strftime("%Y-%m-%d")
        e = e.drop_duplicates("test_id").set_index("test_id")[["patient_id", "date"]]
        parts.append(e)
    return pd.concat(parts)[lambda d: ~d.index.duplicated()]


def evaluate(outdir):
    sys.argv = ["x"]
    spec = importlib.util.spec_from_file_location(
        "abl", os.path.join(IFI_DIR, "Experiments", "NLP_contribution", "ablation_from_runs.py"))
    abl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(abl)
    gold = abl.load_gold(abl.GOLD_PATH)
    TP, FN, FP, s, p, f = abl.evaluate(abl.load_predictions(os.path.join(outdir, "epi_table.csv")), gold)
    return TP + FP, TP, FN, FP, s, p, f


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-merged", required=True, help="Reference tmp_merged(2).csv (production positives + structured tests)")
    ap.add_argument("--all-db-tests", required=True)
    ap.add_argument("--dates", default=None, help="CT/PET patient/date export (export_ct_pet_test_dates.py)")
    ap.add_argument("--swap", action="append", default=[], help="NAME=PROD_CSV,BERT_CSV (histo|ct|pet)")
    ap.add_argument("--veto-only", action="store_true", help="Only reject production positives; do not add BERT-only positives")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    ref = pd.read_csv(args.ref_merged, dtype=str, keep_default_na=False)
    dates = load_dates(args.all_db_tests, args.dates)
    os.makedirs(args.outdir, exist_ok=True)
    shutil.copy(args.all_db_tests, os.path.join(args.outdir, "all_db_tests.csv"))

    out = ref
    for spec_ in args.swap:
        name, paths = spec_.split("=", 1)
        prod_csv, bert_csv = paths.split(",")
        ot = ORDER_TYPE[name]
        prod, bert = positive_ids(prod_csv), positive_ids(bert_csv)
        is_type = out["t1.order_type_id"] == ot
        n_before = int(is_type.sum())
        out = out[~(is_type & ~out["t1.test_id"].isin(bert))]
        rejected = n_before - int((out["t1.order_type_id"] == ot).sum())
        added_rows = []
        missing = 0
        if not args.veto_only:
            tmpl = ref[ref["t1.order_type_id"] == ot]
            tmpl = tmpl[tmpl["system_ifi_classification"] == tmpl["system_ifi_classification"].mode()[0]].iloc[0].to_dict()
            for tid in sorted(bert - prod):
                if tid not in dates.index or pd.isna(dates.loc[tid, "date"]):
                    missing += 1
                    continue
                r = dict(tmpl)
                r.update({"t1.test_id": tid, "t1.patient_id": dates.loc[tid, "patient_id"],
                          "t1.collection_date": dates.loc[tid, "date"],
                          "t1.encounter_id": "", "t1.order_id": "", "t7.order_name": ""})
                added_rows.append(r)
            out = pd.concat([out, pd.DataFrame(added_rows, columns=ref.columns)], ignore_index=True)
        print(f"{name}: production+ {len(prod):,}, BERT+ {len(bert):,}; rejected {rejected:,} reference rows, "
              f"added {len(added_rows):,} BERT-only tests ({missing} without a date)")

    out.to_csv(os.path.join(args.outdir, "tmp_merged2.csv"), index=False)
    r = subprocess.run([sys.executable, os.path.join(IFI_DIR, "combine_tests.py"), "--outdir", args.outdir],
                       capture_output=True, text=True, cwd=IFI_DIR)
    if r.returncode != 0:
        sys.exit("combine_tests.py failed:\n" + r.stderr[-800:])
    e, TP, FN, FP, s, p, f = evaluate(args.outdir)
    print(f"RESULT episodes={e} TP={TP} FN={FN} FP={FP} sensitivity={100*s:.2f}% precision={100*p:.2f}% F1={f:.4f}")


if __name__ == "__main__":
    main()
