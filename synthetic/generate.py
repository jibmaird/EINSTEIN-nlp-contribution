"""
Generate a synthetic dataset for the NLP-contribution experiments.

The real study data are clinical records and cannot be shared. This script
builds a small, fully synthetic substitute with the same schema, so that every
experiment in the paper can be run end to end.

The data are not a statistical imitation of the real cohort. Each synthetic
patient is hand-designed to exercise one mechanism that the paper reports, so
that running the experiments on this dataset reproduces the *qualitative*
findings (which modality is necessary for which episode, super-additivity,
the date-shift and imaging-only false negatives, the classification collapse)
rather than the published numbers.

Scenario labels map to the paper as follows:

    NEC-CT, NEC-PET, NEC-HISTO   a modality is necessary for one episode
                                 -> individual ablation drops (Table 5)
    RED-CTPET                    one episode supported by both CT and PET,
                                 lost only when both go
                                 -> super-additive effect (Section 3.3)
    STRUCT                       structured-only episode, survives every
                                 ablation -> the non-NLP floor
    DATESHIFT                    episode anchored to an earlier test, falls
                                 outside the +/-30d window
                                 -> "episode date shift" FN (Table 8)
    IMGONLY-FILTERED             imaging-only episode with no corroborating
                                 order within +/-28d, excluded by the filter
                                 -> "imaging-only episode filtered" FN
    NO-NLP-EVIDENCE              gold episode with no NLP-eligible test at all
                                 -> irreducible FN
    FP-*                         detected episodes with no gold counterpart
                                 -> the precision cost

Usage:
    python synthetic/generate.py --outdir synthetic/out
"""

import argparse
import csv
import os
from datetime import date, timedelta

# --------------------------------------------------------------------------
# Reference vocabulary (order_type_id values are those the pipeline uses)
# --------------------------------------------------------------------------

ORDER_TYPES = [
    (1, "Antifungal Sensitivity"),
    (2, "Aspergillus PCR"),
    (3, "Biofire"),
    (4, "CMV"),
    (5, "Cryptococcal antigen, Blood"),
    (6, "Cryptococcal antigen, CSF"),
    (7, "Culture"),
    (8, "Galactomannan"),
    (9, "Histopathology"),
    (10, "Meningitis PCR"),
    (11, "Panfungal PCR"),
    (12, "PJP"),
    (13, "TDM"),
    (14, "CT scan"),
    (15, "PET scan"),
    (16, "Full blood count"),
    (17, "Other"),
]

HISTO, CT, PET = 9, 14, 15
NLP_TYPES = {HISTO, CT, PET}

# Order types that count as corroborating evidence for an imaging-only
# episode (combine_tests.filter_episodes_nlp).
CORROBORATING = {1, 2, 5, 6, 8, 11}

ORGANISMS = {
    "aspergillus_fumigatus": dict(
        organism_id=101, species="Aspergillus fumigatus", genus="Aspergillus",
        class_specific="Aspergillus", class_broad="Mould"),
    "aspergillus_sp": dict(
        organism_id=102, species="Aspergillus species", genus="Aspergillus",
        class_specific="Aspergillus", class_broad="Mould"),
    "candida_albicans": dict(
        organism_id=201, species="Candida albicans", genus="Candida",
        class_specific="Candida and related species", class_broad="Yeast"),
    "rhizopus": dict(
        organism_id=301, species="Rhizopus arrhizus", genus="Rhizopus",
        class_specific="Mucormycota", class_broad="Mould"),
    # NLP-derived tests never identify an organism.
    "unspecified_mould": dict(
        organism_id=901, species="Unspecified mould", genus="Unspecified",
        class_specific="Unspecified", class_broad="Mould"),
    "unspecified": dict(
        organism_id=903, species="Unspecified organism", genus="Unspecified",
        class_specific="Unspecified", class_broad="Unspecified"),
}

CLASSIFICATIONS = {
    "Proven": 1, "Probable": 2, "Possible": 3, "Pending": 4, "Colonisation": 5,
}

SITES = {"Lung": (1, "Lung"), "Blood": (2, "Blood"),
         "Sinus": (3, "Sinus"), "Unknown": (99, "Unknown")}

ANCHOR = date(2020, 6, 1)   # mid-study reference date


def d(offset):
    """Test collection date, ISO -- matches t1.collection_date in the real export."""
    return (ANCHOR + timedelta(days=offset)).isoformat()


def d_gold(offset):
    """Gold-standard date. evaluate.py parses these with %d/%m/%Y."""
    return (ANCHOR + timedelta(days=offset)).strftime("%d/%m/%Y")


# --------------------------------------------------------------------------
# Scenario definitions
#
# Each scenario is one patient. "positive" lists tests that the pipeline
# treats as IFI-positive (they enter episode detection). "ordered" lists
# additional tests that were merely requested -- they never create an episode
# but they can corroborate an imaging-only one.
# --------------------------------------------------------------------------

def scenarios():
    S = []

    def add(label, gold, positive, ordered=(), note=""):
        S.append(dict(label=label, gold=gold, positive=list(positive),
                      ordered=list(ordered), note=note))

    # ---- A modality is necessary for an episode -----------------------
    # Imaging-only episode, corroborated by a galactomannan ORDER (not a
    # positive result) 5 days earlier. Removing the modality leaves no
    # positive test at all, so the episode disappears.
    add("NEC-CT",
        gold=dict(day=0, org="aspergillus_fumigatus", site="Lung",
                  classification="Probable"),
        positive=[dict(day=0, type=CT, org="unspecified_mould",
                       site="Lung", classification="Probable")],
        ordered=[dict(day=-5, type=8)],
        note="CT is the only positive test; removing CT loses the episode.")

    add("NEC-PET",
        gold=dict(day=0, org="aspergillus_fumigatus", site="Lung",
                  classification="Probable"),
        positive=[dict(day=0, type=PET, org="unspecified_mould",
                       site="Lung", classification="Probable")],
        ordered=[dict(day=3, type=2)],
        note="PET is the only positive test; removing PET loses the episode.")

    # Histopathology is not subject to the imaging-only filter, so it needs
    # no corroborating order.
    add("NEC-HISTO",
        gold=dict(day=0, org="rhizopus", site="Sinus",
                  classification="Proven"),
        positive=[dict(day=0, type=HISTO, org="unspecified_mould",
                       site="Sinus", classification="Proven")],
        note="Histopathology is the only positive test.")

    # ---- Redundant support: the super-additive effect ------------------
    # CT and PET 15 days apart merge into one episode (28-day imaging merge
    # window). Either alone still dates the episode inside the gold
    # standard's +/-30-day window, so neither single ablation loses it --
    # only removing both does.
    add("RED-CTPET",
        gold=dict(day=10, org="aspergillus_fumigatus", site="Lung",
                  classification="Probable"),
        positive=[dict(day=0, type=CT, org="unspecified_mould",
                       site="Lung", classification="Probable"),
                  dict(day=15, type=PET, org="unspecified_mould",
                       site="Lung", classification="Probable")],
        ordered=[dict(day=-3, type=8)],
        note="Counted in no single-modality ablation, only in the combined one.")

    add("RED-CTPET-2",
        gold=dict(day=12, org="aspergillus_sp", site="Lung",
                  classification="Probable"),
        positive=[dict(day=5, type=PET, org="unspecified_mould",
                       site="Lung", classification="Probable"),
                  dict(day=20, type=CT, org="unspecified_mould",
                       site="Lung", classification="Probable")],
        ordered=[dict(day=2, type=11)],
        note="Second redundant episode, PET first.")

    # ---- Structured-only episodes: the non-NLP floor ------------------
    for i, (org, site, cls) in enumerate([
            ("aspergillus_fumigatus", "Lung", "Proven"),
            ("candida_albicans", "Blood", "Proven"),
            ("aspergillus_sp", "Lung", "Probable"),
            ("candida_albicans", "Blood", "Proven"),
            ("rhizopus", "Sinus", "Proven"),
    ]):
        add(f"STRUCT-{i+1}",
            gold=dict(day=0, org=org, site=site, classification=cls),
            positive=[dict(day=0, type=7, org=org, site=site,
                           classification=cls)],
            note="Culture-positive; unaffected by any NLP ablation.")

    # Structured episode that NLP additionally dates earlier (lead time).
    add("STRUCT-LEAD",
        gold=dict(day=20, org="aspergillus_fumigatus", site="Lung",
                  classification="Probable"),
        positive=[dict(day=8, type=CT, org="unspecified_mould",
                       site="Lung", classification="Probable"),
                  dict(day=20, type=7, org="aspergillus_fumigatus",
                       site="Lung", classification="Proven")],
        ordered=[dict(day=6, type=8)],
        note="CT anchors the episode 12 days before confirmation: lead time.")

    # ---- Classification: Probable depends structurally on CT ----------
    # combine_tests.py skips imaging tests when assigning a classification,
    # then upgrades Possible/Pending to Probable only when the episode holds
    # a CT (order_type 14) AND a non-imaging supporting test from
    # {2,7,8,9,11}. Removing CT therefore makes Probable unreachable, which
    # is the mechanism behind the classification collapse in Section 3.4.
    for i, org in enumerate(["aspergillus_fumigatus", "aspergillus_sp",
                             "aspergillus_fumigatus"]):
        add(f"CLASS-PROBABLE-{i+1}",
            gold=dict(day=0, org=org, site="Lung",
                      classification="Probable"),
            positive=[dict(day=0, type=CT, org="unspecified_mould",
                           site="Lung", classification="Probable"),
                      dict(day=2, type=11, org="unspecified_mould",
                           site="Lung", classification="Pending")],
            note="Upgraded to Probable at baseline; falls back to Pending without CT.")

    # ---- Persistent false negatives -----------------------------------
    # Episode anchored to a culture 60 days before the gold date; the CT near
    # the gold date merges into that open episode rather than starting a new
    # one, so the episode date stays outside the +/-30-day window.
    add("DATESHIFT",
        gold=dict(day=60, org="aspergillus_fumigatus", site="Lung",
                  classification="Probable"),
        positive=[dict(day=0, type=7, org="aspergillus_fumigatus",
                       site="Lung", classification="Proven"),
                  dict(day=55, type=CT, org="unspecified_mould",
                       site="Lung", classification="Probable")],
        ordered=[dict(day=54, type=8)],
        note="Index date anchored early: FN despite NLP evidence near the gold date.")

    # Imaging-only episode with no corroborating order inside +/-28 days.
    add("IMGONLY-FILTERED",
        gold=dict(day=0, org="aspergillus_fumigatus", site="Lung",
                  classification="Probable"),
        positive=[dict(day=0, type=CT, org="unspecified_mould",
                       site="Lung", classification="Probable")],
        ordered=[dict(day=-90, type=8)],
        note="Excluded by the imaging-only evidentiary filter.")

    # Gold episode with no NLP-eligible test at all: irreducible even under
    # perfect NLP recall. Culture is 100 days before the gold date, so it
    # cannot match either.
    add("NO-NLP-EVIDENCE",
        gold=dict(day=100, org="candida_albicans", site="Blood",
                  classification="Proven"),
        positive=[dict(day=0, type=7, org="candida_albicans",
                       site="Blood", classification="Proven")],
        note="No histopathology, CT or PET exists for this patient.")

    # ---- Recoverable under perfect NLP recall -------------------------
    # The classifier missed an NLP report that does exist for the patient:
    # the test appears in all_db_tests.csv but not among the positives. The
    # upper-bound analysis treats every report of the modality as positive,
    # so these episodes are recovered -- unlike the three false negatives
    # above, which no amount of classifier recall can reach.
    for i, (mod, name) in enumerate([(CT, "CT"), (PET, "PET"),
                                     (HISTO, "histopathology")]):
        add(f"UB-RECOVERABLE-{name.upper()}",
            gold=dict(day=0, org="aspergillus_fumigatus", site="Lung",
                      classification="Probable"),
            positive=[],
            ordered=[dict(day=0, type=mod), dict(day=-4, type=8)],
            note=f"{name} report exists but was scored negative: "
                 f"recovered only under perfect recall.")

    # ---- False positives: the precision cost --------------------------
    # Imaging-positive episodes, corroborated so they survive the filter,
    # with no gold episode for the patient.
    for i in range(14):
        mod = CT if i % 2 == 0 else PET
        add(f"FP-IMG-{i+1}",
            gold=None,
            positive=[dict(day=i * 7, type=mod, org="unspecified_mould",
                           site="Lung", classification="Probable")],
            ordered=[dict(day=i * 7 - 2, type=8)],
            note="Imaging-positive, no gold episode: false positive.")

    # A couple of structured false positives, so precision does not become
    # purely an NLP artefact. Note "Colonisation" would NOT work here:
    # combine_tests.py skips colonisation tests during episode grouping.
    for i in range(3):
        add(f"FP-STRUCT-{i+1}",
            gold=None,
            positive=[dict(day=i * 11, type=7, org="candida_albicans",
                           site="Blood", classification="Pending")],
            note="Structured-only false positive.")

    # Histopathology false positive.
    add("FP-HISTO-1", gold=None,
        positive=[dict(day=0, type=HISTO, org="unspecified",
                       site="Unknown", classification="Pending")],
        note="Histopathology false positive.")

    return S


# --------------------------------------------------------------------------
# Emit the five CSVs the experiments consume
# --------------------------------------------------------------------------

TMP_MERGED_COLS = [
    "t1.test_id", "t1.encounter_id", "t1.patient_id", "t1.order_name_id",
    "t7.order_name", "t1.order_id", "t3.organism_id", "t4.species",
    "t4.genus", "t4.class_broad", "t4.class_specific", "t6.specimen_source",
    "t5.specimen_type", "t1.order_type_id", "t1.is_sterile",
    "t1.collection_date", "final_infection_site_id", "infection_site_name",
    "system_ifi_classification_id", "system_ifi_classification",
]

GOLD_COLS = [
    "gold_doc_id", "patient_id", "ifi_episode_id", "if_episode_organism_id",
    "system_diagnosis_date", "organism_id", "species", "genus",
    "class_specific", "class_broad", "ifi_infection_site_id",
    "infection_site", "system_ifi_classification_id", "ifi_classification",
    "hospital",
]

ORDER_NAMES = {
    7: "Fungal Culture", 8: "Aspergillus Antigen (Galactomannan)",
    2: "Aspergillus Nucleic Acid Detection", 11: "Panfungal PCR",
    9: "Histopathology Report", 14: "CT Chest", 15: "FDG PET-CT",
    1: "Antifungal Sensitivity", 5: "Cryptococcal Antigen, Blood",
    6: "Cryptococcal Antigen, CSF",
}


def build(outdir):
    os.makedirs(outdir, exist_ok=True)
    S = scenarios()

    merged, all_tests, gold, test_type = [], [], [], []
    test_id, patient_id, gold_doc = 1000, 100, 1
    manifest = []

    for sc in S:
        patient_id += 1
        pid_tests = []

        for t in sc["positive"]:
            org = ORGANISMS[t["org"]]
            site_id, site_name = SITES[t["site"]]
            merged.append({
                "t1.test_id": test_id,
                "t1.encounter_id": 50000 + patient_id,
                "t1.patient_id": patient_id,
                "t1.order_name_id": 400 + t["type"],
                "t7.order_name": ORDER_NAMES.get(t["type"], "Other"),
                "t1.order_id": 70000 + test_id,
                "t3.organism_id": org["organism_id"],
                "t4.species": org["species"],
                "t4.genus": org["genus"],
                "t4.class_broad": org["class_broad"],
                "t4.class_specific": org["class_specific"],
                "t6.specimen_source": "" if t["type"] in NLP_TYPES else "Blood",
                "t5.specimen_type": "" if t["type"] in NLP_TYPES else "Blood",
                "t1.order_type_id": t["type"],
                "t1.is_sterile": 1 if t["type"] in (7, 9) else 3,
                "t1.collection_date": d(t["day"]),
                "final_infection_site_id": site_id,
                "infection_site_name": site_name,
                "system_ifi_classification_id":
                    CLASSIFICATIONS[t["classification"]],
                "system_ifi_classification": t["classification"],
            })
            test_type.append({"test_id": test_id, "order_type_id": t["type"]})
            pid_tests.append((test_id, d(t["day"]), t["type"]))
            test_id += 1

        for t in sc["ordered"]:
            test_type.append({"test_id": test_id, "order_type_id": t["type"]})
            pid_tests.append((test_id, d(t["day"]), t["type"]))
            test_id += 1

        for tid, dt, ty in pid_tests:
            all_tests.append({"test_id": tid, "patient_id": patient_id,
                              "test_date": dt, "order_type": ty})

        if sc["gold"]:
            g = sc["gold"]
            org = ORGANISMS[g["org"]]
            site_id, site_name = SITES[g["site"]]
            gold.append({
                "gold_doc_id": gold_doc,
                "patient_id": patient_id,
                "ifi_episode_id": 8000 + gold_doc,
                "if_episode_organism_id": 9000 + gold_doc,
                "system_diagnosis_date": d_gold(g["day"]),
                "organism_id": org["organism_id"],
                "species": org["species"],
                "genus": org["genus"],
                "class_specific": org["class_specific"],
                "class_broad": org["class_broad"],
                "ifi_infection_site_id": site_id,
                "infection_site": site_name,
                "system_ifi_classification_id":
                    CLASSIFICATIONS[g["classification"]],
                "ifi_classification": g["classification"],
                "hospital": "Hospital A" if patient_id % 2 else "Hospital B",
            })
            gold_doc += 1

        manifest.append({
            "label": sc["label"], "patient_id": patient_id,
            "gold_episode": "yes" if sc["gold"] else "no",
            "positive_tests": ";".join(
                f'{ORDER_NAMES.get(t["type"], t["type"])}@{t["day"]:+d}'
                for t in sc["positive"]),
            "note": sc["note"],
        })

    def write(name, rows, cols):
        path = os.path.join(outdir, name)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"  {name:28s} {len(rows):5d} rows")
        return path

    print(f"Writing synthetic dataset to {outdir}/")
    write("tmp_merged2.csv", merged, TMP_MERGED_COLS)
    write("all_db_tests.csv", all_tests,
          ["test_id", "patient_id", "test_date", "order_type"])
    write("gold_standard.csv", gold, GOLD_COLS)
    write("test_order_type.csv", test_type, ["test_id", "order_type_id"])
    write("ref_order_type.csv",
          [{"order_type_id": i, "order_type": n} for i, n in ORDER_TYPES],
          ["order_type_id", "order_type"])
    write("SCENARIOS.csv", manifest,
          ["label", "patient_id", "gold_episode", "positive_tests", "note"])

    n_nlp = sum(1 for r in merged if int(r["t1.order_type_id"]) in NLP_TYPES)
    print(f"\n  patients        {len(S)}")
    print(f"  positive tests  {len(merged)} "
          f"({n_nlp} NLP-derived, {len(merged) - n_nlp} structured)")
    print(f"  gold episodes   {len(gold)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="synthetic/out")
    build(ap.parse_args().outdir)
