#!/usr/bin/env python3
"""
Inference cost: production classifiers vs their BERT counterparts.

The substitution experiment (experiments/bert_substitution.py) shows that
replacing the production classifiers with BERT does not improve episode-level
sensitivity. This measures what the replacement would cost to run.

Three configurations are timed on the same reports and the same hardware:

  production        TF-IDF + gradient boosting over the whole report, as
                    deployed for CT
  BERT whole-report one forward pass per report, as deployed for
                    histopathology and PET
  BERT sentence     one forward pass per sentence, max-pooled to a report
                    label, as deployed for CT

Timing depends on the architecture and the number of tokens, not on the
values of the trained weights, so a fine-tuned checkpoint and the public
base checkpoint of the same architecture cost the same to run. The public
checkpoint is used here so the benchmark is reproducible.

With --texts, reports are read from the `text` column of a CSV. Without it,
synthetic reports of a comparable length are generated, which is sufficient
because cost tracks token count rather than content.

Usage:
    python -m experiments.benchmark_cost
    python -m experiments.benchmark_cost --texts reports.csv --n 300
"""

import argparse
import csv
import os
import random
import statistics
import sys
import time
from pathlib import Path

# Thread count must be fixed before torch and the BLAS backends are imported,
# or the two sides of the comparison get different amounts of CPU.
DEFAULT_THREADS = 8


def set_threads(n):
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ[var] = str(n)


# Median CT report in the study corpus is ~164 words; this vocabulary only
# has to produce text of a realistic length and token distribution.
VOCAB = """chest computed tomography contrast ground glass opacity nodule
consolidation bronchial wall thickening pleural effusion lymphadenopathy
mediastinal hilar cavitation halo sign reverse nodular infiltrate segmental
lobe right left upper lower lingula bilateral diffuse patchy peripheral
interstitial septal thickening traction bronchiectasis fibrosis scarring
atelectasis no evidence of previous comparison stable increased decreased
resolved new interval change clinical correlation recommended follow up
imaging findings consistent with infection inflammation malignancy
appearances suggest further assessment radiology report impression""".split()


def synthetic_reports(n, median_words=164, seed=0):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        k = max(40, int(rng.gauss(median_words, median_words * 0.35)))
        out.append(" ".join(rng.choice(VOCAB) for _ in range(k)))
    return out


def load_reports(path, n):
    csv.field_size_limit(10 ** 9)
    texts = []
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f):
            t = (r.get("text") or "").strip()
            if t:
                texts.append(t)
            if len(texts) >= n:
                break
    if not texts:
        sys.exit(f"No `text` column content found in {path}")
    return texts


def split_sentences(text):
    parts = [s.strip() for s in text.replace("\n", " ").split(".")]
    return [s for s in parts if s] or [text]


def timed(fn, repeats):
    """Best-of-N wall clock, to suppress scheduling noise."""
    runs = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        runs.append(time.perf_counter() - t0)
    return min(runs), statistics.median(runs)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--texts", default=None,
                    help="CSV with a `text` column of real reports")
    ap.add_argument("--n", type=int, default=300, help="reports to time")
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--model", default="emilyalsentzer/Bio_ClinicalBERT")
    ap.add_argument("--production-model", default=None,
                    help="joblib Pipeline(TfidfVectorizer, classifier); "
                         "without it, an equivalent pipeline is fitted here")
    ap.add_argument("--corpus-size", type=int, default=47971,
                    help="reports scored per full pipeline run "
                         "(14,741 histopathology + 18,629 PET + 14,601 CT)")
    args = ap.parse_args()

    set_threads(args.threads)
    import joblib                      # noqa: E402
    import torch                       # noqa: E402
    torch.set_num_threads(args.threads)

    texts = (load_reports(args.texts, args.n) if args.texts
             else synthetic_reports(args.n))
    wc = [len(t.split()) for t in texts]
    n_sent = sum(len(split_sentences(t)) for t in texts)
    print(f"{len(texts)} reports "
          f"({'real' if args.texts else 'synthetic'}), "
          f"median {int(statistics.median(wc))} words, "
          f"{n_sent / len(texts):.1f} sentences each")
    print(f"{args.threads} CPU threads, best of {args.repeats}\n")

    results = {}

    # ---------------- production ----------------
    if args.production_model:
        pipe = joblib.load(args.production_model)
    else:
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn.pipeline import make_pipeline
        fit_texts = texts + synthetic_reports(200, seed=1)
        y = [i % 2 for i in range(len(fit_texts))]
        pipe = make_pipeline(TfidfVectorizer(min_df=2),
                             GradientBoostingClassifier())
        pipe.fit(fit_texts, y)

    pipe.predict_proba(texts[:20])
    best, med = timed(lambda: pipe.predict_proba(texts), args.repeats)
    results["production"] = best
    print(f"production  TF-IDF + gradient boosting")
    print(f"  {best * 1000 / len(texts):9.3f} ms/report   "
          f"{len(texts) / best:10.1f} reports/s")

    # ---------------- BERT ----------------
    from transformers import AutoModelForSequenceClassification, AutoTokenizer
    import logging
    logging.getLogger("transformers").setLevel(logging.ERROR)
    tok = AutoTokenizer.from_pretrained(args.model)
    mdl = AutoModelForSequenceClassification.from_pretrained(
        args.model, num_labels=2)
    mdl.eval()
    nparam = sum(p.numel() for p in mdl.parameters())

    @torch.no_grad()
    def run(items):
        for i in range(0, len(items), args.batch_size):
            batch = tok(items[i:i + args.batch_size], truncation=True,
                        max_length=512, padding=True, return_tensors="pt")
            mdl(**batch)

    run(texts[:args.batch_size])
    best, med = timed(lambda: run(texts), args.repeats)
    results["bert_report"] = best
    print(f"\nBERT whole-report  {args.model} ({nparam / 1e6:.0f}M parameters)")
    print(f"  {best * 1000 / len(texts):9.3f} ms/report   "
          f"{len(texts) / best:10.1f} reports/s")

    sentences = [s for t in texts for s in split_sentences(t)]
    run(sentences[:args.batch_size])
    best, med = timed(lambda: run(sentences), args.repeats)
    results["bert_sentence"] = best
    print(f"\nBERT sentence-level  ({len(sentences)} sentences)")
    print(f"  {best * 1000 / len(texts):9.3f} ms/report   "
          f"{len(texts) / best:10.1f} reports/s")

    # ---------------- summary ----------------
    prod = results["production"]
    print(f"\n{'':22}{'x production':>14}{'per full run':>16}")
    print("-" * 52)
    for key, label in (("production", "production"),
                       ("bert_report", "BERT whole-report"),
                       ("bert_sentence", "BERT sentence-level")):
        t = results[key]
        total = t / len(texts) * args.corpus_size
        ratio = t / prod
        total_s = (f"{total:.1f} s" if total < 120
                   else f"{total / 60:.0f} min" if total < 7200
                   else f"{total / 3600:.1f} h")
        print(f"{label:<22}{ratio:>13.0f}x{total_s:>16}")
    print(f"\nPer full pipeline run over {args.corpus_size:,} reports, "
          f"on {args.threads} CPU threads.")


if __name__ == "__main__":
    main()
