#!/usr/bin/env python3
"""Plot authenticated descriptive outputs; intervals are species IQRs, not CIs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    axis_dir = root / "results/fcp_real_447500_white_axis_decomposition_20261018"
    env_dir = root / "results/fcp_real_447500_environment_description_20261018"
    out = root / "docs/figures/fcp_447500_axis_environment_20261018"
    axis = json.loads((axis_dir / "result.json").read_text())
    receipt = json.loads((env_dir / "execution_receipt.json").read_text())
    for folder, hashes in [(axis_dir, axis["output_sha256"]), (env_dir, receipt["output_sha256"])]:
        for name, expected in hashes.items():
            if hashlib.sha256((folder / name).read_bytes()).hexdigest() != expected:
                raise ValueError(f"Output hash mismatch: {folder / name}")
    env = pd.read_csv(env_dir / "group_environment_summary.csv")
    env = env[(env["group"] == "pooled") & (env["metric"] == "standardized_difference")]
    primary = env[env["panel"] == "primary_ge5"].set_index("variable")
    sensitivity = env[env["panel"] == "sensitivity_ge10"].set_index("variable")
    pooled = axis["pooled_spatial"]
    global_share = pooled["white_nonwhite"]["observed_mean_global_mismatch"] / pooled["total"]["observed_mean_global_mismatch"]
    effect_share = axis["pooled_effect_decomposition"]["white_nonwhite_over_positive_pooled_total"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 7.8), gridspec_kw={"width_ratios": [1.0, 1.35]})
    fig.subplots_adjust(left=.085, right=.98, top=.81, bottom=.30, wspace=.43)
    fig.suptitle("A recurring visible colour axis; heterogeneous climate separation", x=.06, ha="left", y=.97, fontsize=16, fontweight="bold")
    fig.text(.06, .91, "447,500 source photographs across five cohorts; eligible subsets shown below", fontsize=10.5, color="#45536a")
    white_color, chromatic_color = "#49669d", "#c0cbdc"
    shares = np.array([global_share, effect_share]) * 100
    ys = np.array([1., 0.])
    left.barh(ys, shares, color=white_color, height=.40, label="White–chromatic")
    left.barh(ys, 100-shares, left=shares, color=chromatic_color, height=.40, label="Chromatic–chromatic")
    for y, share in zip(ys, shares):
        left.text(share/2, y, f"{share:.1f}%", ha="center", va="center", color="white", fontweight="bold")
        left.text(share+(100-share)/2, y, f"{100-share:.1f}%", ha="center", va="center", fontsize=9)
    left.set(yticks=ys, yticklabels=["Global colour\nmismatch", "50-km\nsorting effect"], xlim=(0, 100), ylim=(-.55, 1.7), xlabel="Composition of the pooled statistic (%)")
    left.set_title("A  Exact pairwise partition · 1,761 species", loc="left", fontsize=11, pad=14, fontweight="bold")
    left.legend(loc="upper left", bbox_to_anchor=(.085, .185), bbox_transform=fig.transFigure, borderaxespad=0, frameon=False, fontsize=9)
    left.text(0, 1.42, "Pooled Δ = 0.014763 + 0.002852 = 0.017615", fontsize=9, color="#45536a")
    variables = ["BIO1", "BIO4", "BIO5", "BIO12", "BIO14", "BIO15", "SRADmean"]
    names = ["Annual temperature", "Temperature seasonality", "Warmest-month maximum", "Annual precipitation", "Driest-month precipitation", "Precipitation seasonality", "Solar radiation"]
    for i, (v, name) in enumerate(zip(variables, names)):
        p, s = primary.loc[v], sensitivity.loc[v]
        right.plot([p.q25, p.q75], [i, i], color="#b2bac7", lw=3, solid_capstyle="round")
        right.scatter([p["median"]], [i], color="#49669d", s=42, zorder=3, label="≥5 photos per colour: 1,574 species" if i == 0 else None)
        right.scatter([s["median"]], [i+.15], marker="^", color="#b26137", s=28, zorder=3, label="≥10 per colour: 876 species" if i == 0 else None)
    right.axvline(0, color="#7d8697", ls="--", lw=.9, zorder=0)
    right.set(yticks=np.arange(len(names)), yticklabels=names, ylim=(6.65, -.65), xlim=(-.39, .42), xlabel="White − chromatic environment (within-species SD)")
    right.set_title("B  Unadjusted climate separation", loc="left", fontsize=11, pad=14, fontweight="bold")
    right.legend(loc="upper left", bbox_to_anchor=(.557, .185), bbox_transform=fig.transFigure, borderaxespad=0, frameon=False, fontsize=9)
    fig.text(.06, .045, "A: Ratios of species-equal means, not causal variance explained. Baseline composition is already 79.6% white–chromatic.\nB: Dots are species medians; grey lines are species interquartile ranges, not confidence intervals. No environmental p-tests;\n    spatial, observer, background and exposure effects are not adjusted in this description.", fontsize=9, color="#45536a", linespacing=1.45)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "fcp_447500_axis_environment.png", dpi=180)
    svg = out / "fcp_447500_axis_environment.svg"
    fig.savefig(svg)
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)
    print(out / "fcp_447500_axis_environment.png")


if __name__ == "__main__":
    main()
