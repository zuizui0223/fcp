#!/usr/bin/env python3
"""Render the verified five-cohort robustness results without fitting new models."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    result = json.loads((args.results / "result.json").read_text())
    table_path = args.results / "group_statistic_summary.csv"
    digest = hashlib.sha256(table_path.read_bytes()).hexdigest()
    if digest != result["output_sha256"][table_path.name]:
        raise ValueError("The plotted table does not match the verified result receipt")
    data = pd.read_csv(table_path).set_index(["group", "statistic"])
    teal, grey, ink = "#167D8D", "#959DA5", "#252D34"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelcolor": ink, "text.color": ink,
                         "xtick.color": ink, "ytick.color": ink,
                         "svg.fonttype": "none", "svg.hashsalt": "fcp447500"})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.0),
                             gridspec_kw={"width_ratios": [1.05, 1]})
    fig.subplots_adjust(left=.18, right=.965, top=.74, bottom=.29, wspace=.75)
    fig.text(.035, .945, "FCP | 447,500 measured photographs", fontsize=19, weight="bold")
    fig.text(.035, .90, "4,475 species measured; 1,761 species eligible for the original 50-km test", fontsize=11)
    handles = [Line2D([0], [0], marker="o", linestyle="", color=teal, markersize=7, label="Observed"),
               Line2D([0], [0], marker="D", color=grey, markersize=5, label="Matched null mean and central 95% range")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.03, .868),
               ncol=2, frameon=False, fontsize=10)

    def point_pair(ax, row, y):
        lo, hi, middle, observed = [100 * row[k] for k in
                                    ("null_q025", "null_q975", "null_mean", "observed")]
        ax.hlines(y, lo, hi, color=grey, linewidth=2, zorder=2)
        ax.plot(middle, y, "D", color=grey, markersize=5, zorder=3)
        ax.plot(observed, y, "o", color=teal, markersize=7, zorder=4)
        return observed

    ax = axes[0]
    groups = ["discovery", "validation", "prospective", "wave4", "new_U100", "pooled"]
    names = ["Discovery", "Validation", "Prospective", "Fourth cohort", "Additional U100", "All five pooled"]
    labels = []
    for i, (group, name) in enumerate(zip(groups, names)):
        row = data.loc[(group, "fraction_positive")]
        observed = point_pair(ax, row, 5-i)
        ax.text(observed + 1.3, 5-i, f"{observed:.1f}%", va="center", fontsize=10, color=teal)
        labels.append(f"{name}  (n={int(row.n_species):,})")
    ax.set(yticks=np.arange(6), yticklabels=labels[::-1], ylim=(-.6, 5.6), xlim=(30, 67),
           xticks=[30, 40, 50, 60], xlabel="Species with positive spatial effect (%)")
    ax.set_title("A  Prevalence across species", loc="left", fontsize=12, weight="bold", pad=14)
    ax.axhline(.5, color="#D7DBDE", linewidth=.8)

    ax = axes[1]
    methods = ["mean", "median", "drop_top_10pct", "drop_top_20pct"]
    labels = ["Mean (original)", "Median", "Remove highest 10%\n177 species", "Remove highest 20%\n353 species"]
    for i, method in enumerate(methods):
        observed = point_pair(ax, data.loc[("pooled", method)], 3-i)
        ax.text(observed + .16, 3-i, f"{observed:+.2f}", va="center", fontsize=10, color=teal)
    ax.set(yticks=np.arange(4), yticklabels=labels[::-1], ylim=(-.6, 3.6), xlim=(-2.6, 2.55),
           xticks=[-2, -1, 0, 1, 2], xlabel="Global minus local mismatch (percentage points)")
    ax.set_title("B  Robustness to high-effect species", loc="left", fontsize=12, weight="bold", pad=14)
    ax.axvline(0, color=ink, linewidth=.8, linestyle="--", alpha=.6, zorder=1)
    for ax in axes:
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_color("#B7BEC4")
        ax.tick_params(axis="y", length=0, pad=9)
        ax.grid(axis="x", color="#E7EAEC", linewidth=.7, zorder=0)
    fig.text(.035, .175, "All four pooled robustness tests: raw P = 0.005; Holm P = 0.020 within this four-test family.", fontsize=10, weight="bold")
    fig.text(.035, .125, "Grey ranges come from 199 original within-species vertex permutations; they are not confidence intervals on observed effects.", fontsize=9.3)
    fig.text(.035, .088, "Highest-effect species are selected again in every null world. A negative trimmed value can still exceed its matched null.", fontsize=9.3)
    fig.text(.035, .051, "The 50-km radius and added tests are post hoc. Photographic sorting does not identify inherited morphs or a selection mechanism.", fontsize=9.3)
    args.outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.outdir / "fcp_447500_robust_matched_null.png", dpi=180, facecolor="white")
    svg_path = args.outdir / "fcp_447500_robust_matched_null.svg"
    fig.savefig(svg_path, facecolor="white", metadata={"Date": None})
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
    plt.close(fig)
    print("FCP_447500_ROBUST_FIGURE_COMPLETE")


if __name__ == "__main__":
    main()
