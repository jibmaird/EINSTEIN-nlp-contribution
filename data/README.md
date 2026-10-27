# Input data format

The experiments read five CSVs from the directory given as `--data`. To run
against your own data, produce files in these schemas.

| File | Rows are | Used by |
|---|---|---|
| `tmp_merged2.csv` | one IFI-positive test | episode detection |
| `all_db_tests.csv` | every test ordered, positive or not | the imaging-only filter |
| `gold_standard.csv` | one adjudicated IFI episode | evaluation |
| `ref_order_type.csv` | `order_type_id` -> name | evaluation |
| `test_order_type.csv` | `test_id` -> `order_type_id` | evaluation |

## `tmp_merged2.csv`

Tests the pipeline considers IFI-positive, whether from a structured result or
from an NLP classifier. These are the rows an ablation removes. Column names
keep the `tN.` prefixes of the production query.

| Column | Notes |
|---|---|
| `t1.test_id` | unique |
| `t1.encounter_id` | |
| `t1.patient_id` | |
| `t1.order_name_id`, `t7.order_name` | the specific assay |
| `t1.order_id` | |
| `t3.organism_id` | |
| `t4.species`, `t4.genus`, `t4.class_specific`, `t4.class_broad` | organism hierarchy; NLP-derived tests carry `Unspecified*` values, since the classifiers do not identify an organism |
| `t6.specimen_source`, `t5.specimen_type` | may be empty for imaging |
| `t1.order_type_id` | **9** histopathology, **14** CT, **15** PET; see `ref_order_type.csv` for the rest |
| `t1.is_sterile` | 1 sterile, 3 non-sterile |
| `t1.collection_date` | `YYYY-MM-DD` |
| `final_infection_site_id`, `infection_site_name` | |
| `system_ifi_classification_id`, `system_ifi_classification` | `Proven`, `Probable`, `Possible`, `Pending`, `Colonisation` |

Two behaviours are easy to trip over:

- **`Colonisation` tests are skipped** during episode grouping, so they never
  create an episode.
- **Imaging tests are skipped when a classification is assigned.** An
  imaging-only episode falls back to `Possible`, and is upgraded to `Probable`
  only if the episode also holds a non-imaging test of type 2, 7, 8, 9 or 11.
  This is why removing CT makes `Probable` unreachable (Section 3.4).

## `all_db_tests.csv`

Every test ordered for a patient, regardless of result. Only used to decide
whether an imaging-only episode has corroborating evidence.

| Column | Notes |
|---|---|
| `test_id` | |
| `patient_id` | |
| `test_date` | `YYYY-MM-DD` |
| `order_type` | the `order_type_id` value |

An episode whose tests are *all* CT or PET is kept only if the patient has
some test of type **1, 2, 5, 6, 8 or 11** (antifungal sensitivity, Aspergillus
PCR, cryptococcal antigen, galactomannan, panfungal PCR) within **±28 days**
of the episode date. Otherwise it is discarded.

Note this corroboration comes from *ordered* tests, not positive ones. That is
what lets a CT or PET be the sole positive test behind a true-positive
episode — and therefore what makes a modality *necessary*.

## `gold_standard.csv`

One row per clinician-adjudicated episode.

| Column | Notes |
|---|---|
| `gold_doc_id` | |
| `patient_id` | must match `tmp_merged2.csv` |
| `ifi_episode_id`, `if_episode_organism_id` | |
| `system_diagnosis_date` | **`DD/MM/YYYY`** — note this differs from the test files |
| `organism_id`, `species`, `genus`, `class_specific`, `class_broad` | |
| `ifi_infection_site_id`, `infection_site` | |
| `system_ifi_classification_id`, `ifi_classification` | |
| `hospital` | |

## Matching

A detected episode matches a gold episode when it is the same patient and the
diagnosis dates are within **±30 days**. Organisms are compared
hierarchically — species, then genus, then specific class, then broad class —
accepting the first level at which a date-compatible episode is found. An
episode of unspecified organism matches on patient and date alone. Each
detected episode matches at most one gold episode.

A gold episode with a match is a true positive, one without a false negative,
and a detected episode matching nothing a false positive. True negatives are
undefined at episode level, so specificity is not computed.

`evaluate.py` restricts both sides to diagnosis dates in **2018–2022**; the
range is hardcoded in `date_filter`. Synthetic dates are anchored in 2020 so
they fall inside it.
