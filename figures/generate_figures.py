"""
Figures 1 and 2 of the paper.

Two sources of data, and the figure is labelled with whichever was used:

    --results DIR   plot the output of an actual run, reading
                    ablation_results.csv and classification_confusion.csv
                    written by experiments/ablation.py and
                    experiments/attributes.py

    --published     plot the values published in the paper, which come from
                    clinical data not included in this repository

Default is --results results/. Plotting a synthetic run and plotting the
published figures give different numbers, by design: the synthetic data
reproduce the paper's qualitative findings, not its values.

Usage:
    python figures/generate_figures.py --results results/
    python figures/generate_figures.py --published
"""

import argparse
import csv
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Published values (clinical data, not reproducible from this repository).
# Used only with --published.
# ---------------------------------------------------------------------------

PUBLISHED_CONDITIONS  = ["Baseline", "No Histo", "No PET", "No CT", "No NLP"]
PUBLISHED_SENSITIVITY = [92.57, 92.00, 92.00, 90.29, 86.86]   # percent
PUBLISHED_PRECISION   = [29.40, 30.67, 31.38, 33.33, 41.76]   # percent

# ---------------------------------------------------------------------------
# Figure: Sensitivity-Precision Trade-off (grouped bar chart, dual y-axis)
# ---------------------------------------------------------------------------

def make_tradeoff_figure(out_dir, conditions, sensitivity, precision,
                         source="") -> None:
    os.makedirs(out_dir, exist_ok=True)

    x = np.arange(len(conditions))
    bar_width = 0.35

    # Colour palette: dark blue for sensitivity, muted orange for precision
    sens_colour = "#2c6fad"
    prec_colour = "#e08030"
    threshold_colour = "#cc3333"

    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax2 = ax1.twinx()

    bars1 = ax1.bar(x - bar_width / 2, sensitivity, bar_width,
                    color=sens_colour, alpha=0.88, label="Sensitivity (%)", zorder=3)
    bars2 = ax2.bar(x + bar_width / 2, precision, bar_width,
                    color=prec_colour, alpha=0.88, label="Precision (%)", zorder=3)

    # 90 % sensitivity reference line, drawn only when it falls inside the
    # plotted range -- a synthetic run sits well below it, and a legend entry
    # for a line that is not on the axes would be misleading.
    show_reference = min(sensitivity) <= 90 <= max(sensitivity) + 5
    if show_reference:
        ax1.axhline(y=90, color=threshold_colour, linestyle="--",
                    linewidth=1.4, zorder=4,
                    label="90% sensitivity (reference)")

    # Axis formatting
    ax1.set_ylabel("Sensitivity (%)", color=sens_colour, fontsize=11)
    ax2.set_ylabel("Precision (%)", color=prec_colour, fontsize=11)
    ax1.tick_params(axis="y", labelcolor=sens_colour)
    ax2.tick_params(axis="y", labelcolor=prec_colour)

    # Follow the data: a synthetic run sits well outside the published range.
    # More headroom above than below, so the value labels on the tallest bars
    # clear the legend.
    def limits(values, headroom=0.55, footroom=0.25):
        lo, hi = min(values), max(values)
        span = hi - lo
        return (max(0.0, lo - max(2.0, span * footroom)),
                min(100.0, hi + max(4.0, span * headroom)))

    ax1.set_ylim(*limits(sensitivity))
    ax2.set_ylim(*limits(precision))

    ax1.set_xticks(x)
    ax1.set_xticklabels(conditions, fontsize=10)
    ax1.set_xlabel("Ablation Condition", fontsize=11)

    ax1.yaxis.grid(True, linestyle=":", alpha=0.5, zorder=0)
    ax1.set_axisbelow(True)

    # Value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2, height + 0.15,
                 f"{height:.2f}", ha="center", va="bottom",
                 fontsize=8, color=sens_colour)

    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2, height + 0.15,
                 f"{height:.2f}", ha="center", va="bottom",
                 fontsize=8, color=prec_colour)

    # Legend
    legend_handles = [
        mpatches.Patch(color=sens_colour, alpha=0.88, label="Sensitivity (%)"),
        mpatches.Patch(color=prec_colour, alpha=0.88, label="Precision (%)"),
    ]
    if show_reference:
        legend_handles.append(
            plt.Line2D([0], [0], color=threshold_colour, linestyle="--",
                       linewidth=1.4, label="90% sensitivity (reference)"))
    ax1.legend(handles=legend_handles, loc="best", fontsize=9,
               framealpha=0.9)

    plt.title("Sensitivity–Precision Trade-off Across Ablation Conditions",
              fontsize=12, pad=10)
    if source:
        fig.text(0.99, 0.01, source, ha="right", va="bottom",
                 fontsize=7, color="#666666")
    fig.tight_layout()

    for ext in ("pdf", "png"):
        path = os.path.join(out_dir, f"fig_tradeoff.{ext}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved: {path}")

    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure: EORTC classification confusion matrices, baseline vs no-NLP
#
# Published counts, from the study snapshot. Both panels use the 151 episodes
# detected under BOTH conditions, so they are directly comparable and differ
# only in classification, not detection. Used only with --published;
# experiments/attributes.py computes the same matrices for an actual run.
# ---------------------------------------------------------------------------

CLASS_ORDER = ["Proven", "Probable", "Possible", "Pending"]

# rows = gold standard, columns = system-assigned
PUBLISHED_CONFUSION_BASELINE = [
    [75,  0, 0,  1],   # gold Proven
    [ 4, 26, 0, 33],   # gold Probable
    [ 0,  4, 0,  1],   # gold Possible
    [ 0,  2, 0,  5],   # gold Pending
]
PUBLISHED_CONFUSION_NO_NLP = [
    [74, 0, 0,  2],    # gold Proven
    [ 4, 0, 0, 59],    # gold Probable
    [ 0, 0, 0,  5],    # gold Possible
    [ 0, 0, 0,  7],    # gold Pending
]


def make_classification_figure(out_dir, baseline, no_nlp,
                               labels=None, source="") -> None:
    os.makedirs(out_dir, exist_ok=True)
    labels = labels or CLASS_ORDER
    panels = [("Baseline (all NLP)", np.array(baseline)),
              ("No NLP", np.array(no_nlp))]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    vmax = max(m.max() for _, m in panels)

    for ax, (title, m) in zip(axes, panels):
        ax.imshow(m, cmap="Blues", vmin=0, vmax=vmax)
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("System-assigned classification", fontsize=10)
        ax.set_ylabel("Gold standard classification", fontsize=10)
        ax.set_title(f"{title}  (n = {int(m.sum())})", fontsize=11, pad=8)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = int(m[i, j])
                ax.text(j, i, str(v), ha="center", va="center", fontsize=10,
                        color="white" if v > vmax * 0.55 else "#222222",
                        fontweight="bold" if i == j else "normal")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.5)
        ax.tick_params(which="minor", length=0)

    if source:
        fig.text(0.99, 0.01, source, ha="right", va="bottom",
                 fontsize=7, color="#666666")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(out_dir, f"fig_classification_confusion.{ext}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved: {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Reading an actual run
# ---------------------------------------------------------------------------

CONDITION_LABELS = {"baseline": "Baseline", "no_histo": "No Histo",
                    "no_pet": "No PET", "no_ct": "No CT", "no_nlp": "No NLP"}
PLOT_ORDER = ["baseline", "no_histo", "no_pet", "no_ct", "no_nlp"]


def load_ablation(results_dir):
    """Read ablation_results.csv -> (conditions, sensitivity, precision)."""
    path = os.path.join(results_dir, "ablation_results.csv")
    if not os.path.exists(path):
        sys.exit(f"Missing {path}\n"
                 f"Run: python -m experiments.ablation --data synthetic/out "
                 f"--outdir {results_dir}")
    with open(path, newline="", encoding="utf-8") as f:
        rows = {r["key"]: r for r in csv.DictReader(f)}
    keys = [k for k in PLOT_ORDER if k in rows] or list(rows)
    return ([CONDITION_LABELS.get(k, k) for k in keys],
            [float(rows[k]["sensitivity"]) for k in keys],
            [float(rows[k]["precision"]) for k in keys])


def load_confusion(results_dir):
    """Read classification_confusion.csv -> (baseline, no_nlp, labels)."""
    path = os.path.join(results_dir, "classification_confusion.csv")
    if not os.path.exists(path):
        sys.exit(f"Missing {path}\n"
                 f"Run: python -m experiments.attributes --data synthetic/out "
                 f"--outdir {results_dir}")
    cells = {"baseline": {}, "no_nlp": {}}
    seen = set()
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g, s = r["gold_classification"], r["system_classification"]
            cells.setdefault(r["condition"], {})[(g, s)] = int(r["n"])
            seen.add(g)
            seen.add(s)
    labels = [c for c in CLASS_ORDER if c in seen] + \
             sorted(seen - set(CLASS_ORDER))

    def matrix(cond):
        return [[cells.get(cond, {}).get((g, s), 0) for s in labels]
                for g in labels]

    return matrix("baseline"), matrix("no_nlp"), labels


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--results", metavar="DIR", default=None,
                     help="plot the output of a run in DIR (default: results/)")
    src.add_argument("--published", action="store_true",
                     help="plot the values published in the paper")
    ap.add_argument("--outdir", default=None,
                    help="where to write the figures "
                         "(default: <results>/figures, or figures/ "
                         "with --published)")
    args = ap.parse_args()

    if args.published:
        out = args.outdir or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "published")
        note = "published values (clinical data)"
        print(f"Plotting {note}")
        make_tradeoff_figure(out, PUBLISHED_CONDITIONS, PUBLISHED_SENSITIVITY,
                             PUBLISHED_PRECISION, source=note)
        make_classification_figure(out, PUBLISHED_CONFUSION_BASELINE,
                                   PUBLISHED_CONFUSION_NO_NLP,
                                   CLASS_ORDER, source=note)
    else:
        results = args.results or "results"
        out = args.outdir or os.path.join(results, "figures")
        note = f"computed from {os.path.normpath(results)}"
        print(f"Plotting {note}")
        conds, sens, prec = load_ablation(results)
        base, no_nlp, labels = load_confusion(results)
        make_tradeoff_figure(out, conds, sens, prec, source=note)
        make_classification_figure(out, base, no_nlp, labels, source=note)


if __name__ == "__main__":
    main()
