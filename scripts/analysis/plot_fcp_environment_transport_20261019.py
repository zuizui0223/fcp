#!/usr/bin/env python3
"""Draw the executed common-unit and held-out prediction results from CSVs.

Each plotted estimate/interval is also saved in plotted_values.csv. No inferred
or illustrative data enter these scientific figures.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
NATIVE = ROOT / "results/fcp_environment_native_units_20261019"
TRANSFER = ROOT / "results/fcp_environment_transport_20261019"
VARIABLES = ["BIO5", "BIO14", "SRADmean"]
AXIS_TITLES = ["BIO5\nWarmest-month maximum", "BIO14\nDriest-month precipitation", "Mean solar radiation"]
UNIT_LABELS = ["Log odds per 1°C", "Log odds per 10 mm", "Log odds per 1,000 kJ m⁻² day⁻¹"]
TEAL, BLUE, GRAY = "#147D79", "#305F9A", "#4D5661"
LIGHT, GRID = "#B1C9C8", "#E9EDF0"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def style_axis(ax):
    ax.axvline(0., color="#89949D", linewidth=1., zorder=0)
    ax.grid(axis="x", color=GRID, linewidth=.6, zorder=-1)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#B9C1C8")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.tick_params(axis="x", color="#B9C1C8")


def record_and_draw(records, ax, *, panel, source, row_key, estimate, low, high,
                    y, color=TEAL, family_low=None, family_high=None, multiplier=1.,
                    n_species=None, n_photos=None):
    vals = np.array([estimate, low, high], float)
    if not np.isfinite(vals).all() or not low <= estimate <= high:
        raise ValueError(f"Invalid interval in {panel}/{row_key}")
    if family_low is not None:
        if not family_low <= low <= high <= family_high:
            raise ValueError(f"Family interval does not enclose nominal interval: {row_key}")
        ax.plot(np.array([family_low, family_high]) * multiplier, [y, y],
                color=LIGHT, lw=6., solid_capstyle="round", zorder=2)
    ax.plot(np.array([low, high]) * multiplier, [y, y], color=color,
            lw=2., solid_capstyle="round", zorder=3)
    ax.plot(estimate * multiplier, y, "o", color=color, markeredgecolor="white",
            markeredgewidth=.8, markersize=7., zorder=4)
    records.append({"panel":panel, "source_csv":str(source.relative_to(ROOT)),
                    "source_row_key":row_key, "estimate":estimate,
                    "nominal_low":low, "nominal_high":high,
                    "family_low":family_low, "family_high":family_high,
                    "display_multiplier":multiplier, "n_species":n_species,
                    "n_photos":n_photos})


def write_figure(fig, outdir, stem):
    png, svg = outdir / (stem + ".png"), outdir / (stem + ".svg")
    fig.savefig(png, dpi=180, facecolor="white", metadata={"Software":"FCP source-backed plotting"})
    fig.savefig(svg, facecolor="white", metadata={"Date":None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)
    return [png, svg]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path,
                        default=ROOT / "results/fcp_environment_transport_synthesis_20261019")
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    sources = [NATIVE / "primary_summary.csv", TRANSFER / "genus_summary.csv",
               TRANSFER / "spatial_summary.csv", NATIVE / "slope_modifier_contrasts.csv"]
    native, genus, spatial, traits = [pd.read_csv(p) for p in sources]
    if len(native) != 3 or len(genus) != 2 or len(spatial) != 4 or len(traits) != 9:
        raise ValueError("Unexpected executed-result comparison count")
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10,
                         "axes.titlesize":11, "axes.labelsize":10,
                         "svg.hashsalt":"FCP-20261019-environment-transport"})
    records = []
    fig = plt.figure(figsize=(13., 11.6))
    grid = fig.add_gridspec(3, 3, left=.18, right=.94, bottom=.13, top=.875,
                           height_ratios=[1., 1.05, 1.8], hspace=.73, wspace=.42)
    fig.text(.07, .965, "FCP · environmental association and transfer", fontsize=18,
             fontweight="bold", ha="left", color="#1F303B")
    fig.text(.07, .934, "1,574 species · 86,114 measured photographs · 727 genera",
             fontsize=11, color=GRAY, ha="left")
    for j, variable in enumerate(VARIABLES):
        ax = fig.add_subplot(grid[0, j])
        r = native.loc[native.variable.eq(variable)].iloc[0]
        record_and_draw(records, ax, panel="A", source=sources[0], row_key=variable + "|equal_species",
                        estimate=r.equal_species_unshrunk_mean, low=r.genus_cluster_percentile_low,
                        high=r.genus_cluster_percentile_high, y=1., color=TEAL, n_species=1574, n_photos=86114)
        record_and_draw(records, ax, panel="A", source=sources[0], row_key=variable + "|REML",
                        estimate=r.reml_mean, low=r.reml_mean_ci_low,
                        high=r.reml_mean_ci_high, y=0., color=BLUE, n_species=1574, n_photos=86114)
        ax.set_yticks([1., 0.], ["Equal species", "REML"] if j == 0 else ["", ""])
        ax.set_ylim(-.6, 1.6)
        ax.set_xlabel(UNIT_LABELS[j], labelpad=7)
        ax.set_title(("A  " if j == 0 else "") + AXIS_TITLES[j], loc="left", pad=14)
        style_axis(ax)
        ax.margins(x=.15)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)

    ax = fig.add_subplot(grid[1, :])
    for y, model, label in [(1., "BIO5_component", "BIO5 component of joint fit"),
                            (0., "joint_three", "Joint three-variable score")]:
        r = genus.loc[genus.model.eq(model)].iloc[0]
        record_and_draw(records, ax, panel="B", source=sources[1], row_key=model,
                        estimate=r.pair_logloss_gain, low=r.pair_logloss_gain_ci95_low,
                        high=r.pair_logloss_gain_ci95_high, y=y,
                        family_low=r.pair_logloss_gain_family2_low,
                        family_high=r.pair_logloss_gain_family2_high, multiplier=1000.,
                        n_species=r.n_species, n_photos=r.n_photos)
        ax.text(.995, y, f"AUC {r.pair_auc:.3f}", transform=ax.get_yaxis_transform(),
                ha="right", va="center", color=GRAY, fontsize=10,
                bbox={"facecolor":"white", "edgecolor":"none", "pad":2.})
    ax.set_yticks([1., 0.], ["BIO5 joint-fit\ncomponent", "Three climate\ncomponents"])
    ax.set_ylim(-.6, 1.6)
    ax.set_title("B  Discrimination within species from held-out genera", loc="left", pad=14)
    ax.set_xlabel("Conditional pair log-loss gain (10⁻³ nats per white–chromatic pair)", labelpad=8)
    style_axis(ax)
    ax.margins(x=.14)

    ax = fig.add_subplot(grid[2, :])
    comparisons = [(0., "climate_vs_intercept", "Climate vs intercept\nNo distance buffer"),
                   (0., "geoclimate_vs_geography", "Climate + geography\nvs geography · no buffer"),
                   (50., "climate_vs_intercept", "Climate vs intercept\n50 km buffer"),
                   (50., "geoclimate_vs_geography", "Climate + geography\nvs geography · 50 km buffer")]
    positions = [3.4, 2.4, 1., 0.]
    for y, (buffer, comparison, label) in zip(positions, comparisons):
        r = spatial.loc[spatial.buffer_km.eq(buffer) & spatial.comparison.eq(comparison)].iloc[0]
        record_and_draw(records, ax, panel="C", source=sources[2], row_key=f"{buffer:g}|{comparison}",
                        estimate=r.logloss_gain, low=r.logloss_gain_ci95_low,
                        high=r.logloss_gain_ci95_high, y=y,
                        family_low=r.logloss_gain_family4_low,
                        family_high=r.logloss_gain_family4_high,
                        n_species=r.n_species, n_photos=r.n_photos)
        ax.text(.99, y, f"{int(r.n_species):,} species · {int(r.n_photos):,} photos",
                transform=ax.get_yaxis_transform(), ha="right", va="center", color=GRAY,
                fontsize=9, bbox={"facecolor":"white", "edgecolor":"none", "pad":2.})
    ax.set_yticks(positions, [row[2] for row in comparisons])
    ax.set_ylim(-.6, 4.)
    ax.set_title("C  Prediction in held-out 1° cells within known species", loc="left", pad=14)
    ax.set_xlabel("Binary log-loss gain (nats per held-out photograph)", labelpad=8)
    style_axis(ax)
    ax.margins(x=.15)
    # Keep coverage labels in a distinct right-hand area, including when every
    # climate model performs worse than its paired baseline.
    left, right = ax.get_xlim()
    ax.set_xlim(left, max(right, .32 * (right - left)))
    fig.text(.07, .063, "Positive gains favour the climate score/model. Geography is a limited linear latitude/longitude model.",
             color=GRAY, fontsize=9)
    fig.text(.07, .043, "A: 95% genus-bootstrap intervals (equal species) and modified-HKSJ intervals (REML).",
             color=GRAY, fontsize=9)
    fig.text(.07, .023, "B/C: dark 95% and pale approximate family-adjusted genus-bootstrap intervals. Post hoc internal validation.",
             color=GRAY, fontsize=9)
    artifacts = write_figure(fig, outdir, "fcp_environment_transfer")

    fig, axes = plt.subplots(1, 3, figsize=(13., 5.2))
    fig.subplots_adjust(left=.19, right=.965, bottom=.24, top=.75, wspace=.25)
    fig.text(.06, .925, "GIFT traits · differences in common-unit environmental slopes", fontsize=16,
             fontweight="bold", color="#1F303B")
    fig.text(.06, .855, "Same categories and species as the SD-scale comparison · all nine intervals include zero",
             fontsize=10, color=GRAY)
    trait_order = ["self_fertilization", "lifecycle", "woodiness"]
    if set(traits.trait) != set(trait_order):
        raise ValueError("Unexpected GIFT trait keys")
    if ((traits.percentile_ci_low > 0) | (traits.percentile_ci_high < 0)).any():
        raise ValueError("GIFT caption no longer matches the executed nominal intervals")
    trait_labels = ["Self-fertilization record\npresent − absent", "Annual − perennial", "Woody − non-woody"]
    for j, (ax, variable) in enumerate(zip(axes, VARIABLES)):
        for y, trait in zip([2., 1., 0.], trait_order):
            r = traits.loc[traits.variable.eq(variable) & traits.trait.eq(trait)].iloc[0]
            record_and_draw(records, ax, panel="S1", source=sources[3], row_key=trait + "|" + variable,
                            estimate=r.difference_first_minus_reference, low=r.percentile_ci_low,
                            high=r.percentile_ci_high, y=y,
                            family_low=r.family9_percentile_ci_low, family_high=r.family9_percentile_ci_high,
                            n_species=r.n_first + r.n_reference)
        ax.set_yticks([2., 1., 0.], trait_labels if j == 0 else ["", "", ""])
        ax.set_ylim(-.5, 2.5)
        ax.set_title(variable, pad=14)
        ax.set_xlabel(UNIT_LABELS[j], labelpad=9)
        style_axis(ax)
        ax.margins(x=.12)
    fig.text(.06, .105, "Dark: nominal 95% intervals. Pale: approximate family9 intervals. 5,000 shared-genus draws per trait.", fontsize=9, color=GRAY)
    fig.text(.06, .060, "Equal-species contrasts of fitted slopes; missing traits remain unknown. Physical increments are slope units, not extrapolation validation.", fontsize=9, color=GRAY)
    artifacts += write_figure(fig, outdir, "fcp_native_gift_traits")
    values_path = outdir / "plotted_values.csv"
    pd.DataFrame(records).to_csv(values_path, index=False)
    if len(records) != 21:
        raise ValueError("Expected exactly 21 plotted estimates")
    artifacts.append(values_path)
    manifest = {"schema":"fcp_environment_transfer_figure_v1", "source_sha256":{
        str(p.relative_to(ROOT)):sha(p) for p in sources}, "script_sha256":sha(Path(__file__)),
        "n_plotted_estimates":len(records), "no_synthetic_or_imputed_data":True,
        "output_sha256":{p.name:sha(p) for p in artifacts}}
    (outdir / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status":"complete", "n_plotted_estimates":len(records),
                      "figures":[str(p) for p in artifacts if p.suffix == ".png"]}, indent=2))


if __name__ == "__main__":
    main()
