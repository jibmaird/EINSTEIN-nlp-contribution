# EINSTEIN: NLP Contribution Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

Analysis code for *Impact of NLP Modules in a Real-World Hospital Surveillance
System for Invasive Fungal Infection*.

## Overview

EINSTEIN is a deployed invasive fungal infection (IFI) surveillance system. It
combines structured laboratory results with three NLP classifiers that read
histopathology, CT and PET reports, then groups the positive tests into
episodes of infection per patient. This repository measures what those three
classifiers contribute to the finished surveillance record.

Four experiments, each against the same clinician-adjudicated gold standard:

| Experiment | Question | Paper |
|---|---|---|
| `experiments/ablation.py` | What is lost when a modality is removed? | Tables 3–6, §3.3 |
| `experiments/upper_bound.py` | What is the ceiling under perfect NLP recall, and what stays missed? | Table 8, §3.3 |
| `experiments/attributes.py` | Beyond detection, how do classification and diagnosis date change? | Table 9, Figure 2, §3.4 |
| `experiments/bert_substitution.py` | Would transformer classifiers do better? | §3.5 |
| `experiments/benchmark_cost.py` | What would a transformer cost to run? | §3.5 |

The ablation is the core design: the positive tests of one modality are
dropped from the input, the unmodified pipeline is re-run, and the resulting
episodes are rematched against the gold standard. Nothing about the algorithm
or the matching changes between conditions — only the input.

## The data are synthetic

The study data are clinical records from two hospitals and **cannot be
shared**. This repository ships a synthetic dataset instead, with the same
schema, so every experiment runs end to end.

The synthetic data are **not** a statistical imitation of the real cohort and
**do not reproduce the published numbers**. Each synthetic patient is instead
designed to exercise one mechanism the paper reports, so the experiments
reproduce the *qualitative* findings. `synthetic/out/SCENARIOS.csv` maps every
patient to the claim it exercises:

| Label | Mechanism |
|---|---|
| `NEC-CT`, `NEC-PET`, `NEC-HISTO` | one modality is necessary for an episode |
| `RED-CTPET` | an episode supported by both CT and PET, lost only when both go |
| `STRUCT`, `STRUCT-LEAD` | structured-only episode; NLP dates one earlier |
| `CLASS-PROBABLE` | Probable depends structurally on CT |
| `UB-RECOVERABLE-*` | a report the classifier missed, recovered at perfect recall |
| `DATESHIFT`, `IMGONLY-FILTERED`, `NO-NLP-EVIDENCE` | the three persistent false negatives |
| `FP-*` | detected episodes with no gold counterpart |

## Quick start

```bash
pip install -r requirements.txt

python synthetic/generate.py --outdir synthetic/out
python -m experiments.ablation    --data synthetic/out --outdir results/
python -m experiments.upper_bound --data synthetic/out --outdir results/
python -m experiments.attributes  --data synthetic/out --outdir results/

python figures/generate_figures.py --results results/
```

No database is needed. Expected ablation output on the shipped data:

```
Condition             Excluded   Input   Episodes   TP   FN   FP   Sensitivity %   Precision
Baseline                     0      42         34   14    6   20           70.00       41.18
No Histopathology            2      40         32   13    7   19           65.00       40.62
No CT                       16      26         26   13    7   13           65.00       50.00
No PET                      10      32         26   13    7   13           65.00       50.00
No NLP                      28      14         14    9   11    5           45.00       64.29

Super-additivity check
  drop(no NLP)                               = 25.00 pp
  drop(no CT) + drop(no PET) + drop(no Histo) = 15.00 pp
  -> combined drop exceeds the sum of individual drops, by 10.00 pp
  episodes supported by more than one modality: 2
```

Every qualitative finding of the paper is visible: NLP is needed to hold
sensitivity up, precision moves the other way, removing a modality costs more
in combination than individually, histopathology has the highest per-test
yield despite the lowest volume, the ceiling under perfect recall is still
short of 100%, and `Probable` becomes unassignable without CT.

## Repository structure

```
ifi_pipeline/            the deployed pipeline (see Provenance)
  combine_tests.py       episode detection: groups positive tests into episodes
  evaluate.py            gold-standard matching and metrics
experiments/
  offline.py             runs the pipeline with no database
  ablation.py            remove each modality, then all of them
  upper_bound.py         perfect-recall ceiling, persistent false negatives
  attributes.py          classification and diagnosis date
  bert_substitution.py   transformer comparison (see Limitations)
  benchmark_cost.py      inference cost, production vs BERT
figures/
  generate_figures.py    Figures 1 and 2, from a run or from the paper
synthetic/
  generate.py            builds the synthetic dataset
  out/                   the generated CSVs, committed for convenience
data/README.md           schema of the five input files
```

## Using your own data

The experiments read five CSVs from the directory given as `--data`, described
in [data/README.md](data/README.md):

| File | Rows are |
|---|---|
| `tmp_merged2.csv` | one IFI-positive test |
| `all_db_tests.csv` | every test ordered, positive or not |
| `gold_standard.csv` | one adjudicated IFI episode |
| `ref_order_type.csv` | `order_type_id` → name |
| `test_order_type.csv` | `test_id` → `order_type_id` |

Point `--data` at a directory holding files in that schema. Nothing in
`ifi_pipeline/` or `experiments/` requires a database. `evaluate.py` normally
fetches the test-type map over ODBC; `experiments/offline.py` supplies it from
the two order-type CSVs instead. If you do connect to a database, its identity
is read from the environment:

```bash
export IFI_DB_SERVER=your-server
export IFI_DB_DATABASE=your-database
export IFI_DB_USER=your-username
```

## Provenance of `ifi_pipeline/`

`combine_tests.py` and `evaluate.py` are the production episode-detection
algorithm and gold-standard matcher. They are **not** a contribution of this
paper — the algorithm, the development of the gold standard and the system's
clinical performance are described in a companion clinical paper. They are
included because the experiments cannot be reproduced without them.

They are copied from the deployed system with no functional change. Internal
database hostnames, database names, usernames and gold-standard filenames are
replaced by placeholders or read from the environment.

The deployed system also carries an older copy of the episode-detection code
(`combine_tests_noNLP.py`) which the original ablation driver invoked. It
differs in two code paths: the most-recent-species selection in
`prev_test_recent`, and part of the species-upgrade logic. On the study
snapshot the two produce byte-identical episode tables, so only the current
version is shipped here.

## Figures

`figures/generate_figures.py` plots either a run or the published values, and
labels each figure with which it used:

```bash
python figures/generate_figures.py --results results/   # a run (default)
python figures/generate_figures.py --published          # the paper's values
```

With `--results` it reads `ablation_results.csv` and
`classification_confusion.csv` from that directory, so the figures always
match the numbers the experiments just produced. With `--published` it plots
the values in the paper, which come from clinical data not included here. The
two will not agree, by design.

## Citation

Under review. Citation details will be added on publication.

```bibtex
@article{einstein_nlp_contribution,
  title   = {Impact of NLP Modules in a Real-World Hospital Surveillance
             System for Invasive Fungal Infection},
  author  = {Martinez, D. and Khanina, A. and Rozova, V. and Neoh, C. F. and
             Han, W. and Teng, J. and Cavedon, L. and Slavin, M. A. and
             Thursky, K. A. and Worth, L. J. and Verspoor, K.},
  note    = {Manuscript under review},
  year    = {2026}
}
```

## Limitations

- **The synthetic data do not reproduce the published numbers.** They
  reproduce the qualitative findings only. Figures in the paper come from
  clinical data that cannot be released.
- **`experiments/bert_substitution.py` cannot be run from this repository.**
  It needs per-report predictions from the fine-tuned ClinicalBERT and
  Bio_ClinicalBERT models, which are not included. It is provided as the
  record of what produced §3.5.
- **The 2018–2022 evaluation window is hardcoded** in `evaluate.py`'s
  `date_filter`. Data outside it is silently dropped; synthetic dates are
  anchored in 2020.

## License

MIT — see [LICENSE](LICENSE).
