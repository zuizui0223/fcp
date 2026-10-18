#!/usr/bin/env python3
"""Generate the canonical flower-colour polymorphism manuscript figures.

This is a reporting-only pipeline. It reads frozen repository results and the
artifact-derived reporting source freeze. It does not refit biological models,
change thresholds, rerun null models, or alter any scientific decision.
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


PRIMARY = "#D55E00"
SECONDARY = "#0072B2"
SUPPORT = "#009E73"
NEUTRAL = "#666666"
LIGHT = "#D9D9D9"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def configure_matplotlib() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.5,
            "figure.titlesize": 13,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


def save_pair(fig: plt.Figure, output_dir: Path, stem: str) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(
        pdf,
        bbox_inches="tight",
        facecolor="white",
        metadata={
            "Title": stem,
            "Author": "fcp frozen polymorphism figure pipeline",
            "Creator": "matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)
    return {
        "png": str(png.name),
        "png_sha256": sha256(png),
        "pdf": str(pdf.name),
        "pdf_sha256": sha256(pdf),
    }


def validate_reporting_source(root: Path) -> dict:
    source_dir = root / "results" / "polymorphism_publication_figure_source_20260918"
    manifest = load_json(source_dir / "manifest.json")
    if manifest.get("schema") != "polymorphism_publication_figure_source_v1":
        raise ValueError("unexpected reporting-source manifest schema")
    for name, record in manifest["files"].items():
        path = source_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        if sha256(path) != record["sha256"]:
            raise ValueError(f"reporting-source hash mismatch: {name}")
    return manifest


def figure1(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    d_path = root / "results" / "polymorphism_publication_figure_source_20260918" / "h3a_species_D.csv"
    d = pd.read_csv(d_path)
    expected = {"discovery": 369, "reserve": 363}
    if d.groupby("cohort").size().to_dict() != expected:
        raise ValueError("unexpected H1/H3 D cohort sizes in reporting source")

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.3), gridspec_kw={"width_ratios": [1.0, 1.22]})
    fig.suptitle("From outcome-blind opportunity frame to prospective confirmation", y=1.02, fontweight="bold")

    ax = axes[0]
    bins = np.linspace(0, 0.72, 25)
    for cohort, colour, label in [
        ("discovery", SECONDARY, "Discovery (n=369)"),
        ("reserve", PRIMARY, "Validation (n=363)"),
    ]:
        vals = d.loc[d["cohort"].eq(cohort), "D"].to_numpy(dtype=float)
        ax.hist(vals, bins=bins, density=True, histtype="step", linewidth=2.0, color=colour, label=label)
        ax.axvline(np.median(vals), color=colour, linestyle="--", linewidth=1.2)
    ax.set_xlabel(r"Species-level polymorphism, $D = 1-\sum p_k^2$")
    ax.set_ylabel("Density")
    ax.set_xlim(0, 0.72)
    ax.legend(frameon=False)
    ax.text(
        0.98,
        0.96,
        "High-depth cohorts\nnot a prevalence sample",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color=NEUTRAL,
    )
    panel_label(ax, "A")

    ax = axes[1]
    ax.axis("off")
    box_style = {"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": "#999999"}
    secondary_box_style = {
        "boxstyle": "round,pad=0.40",
        "facecolor": "white",
        "edgecolor": "#999999",
        "linestyle": "--",
    }
    boxes = [
        (0.50, 0.87, "Stage 0 | Outcome-blind frame\niNaturalist metadata: 42,111 species; 4,730 high-depth-capable\nWHY: define candidates before colour"),
        (0.50, 0.70, "Stage 1 | Discovery + validation\niNaturalist images: 500 × 100 photos each; D eligible 369 + 363\nWHY: estimate distributions and hold out species"),
        (0.50, 0.49, "Freeze recurrent target after discovery\nwhite ↔ non-white; W + structured null\nWHY: discovery is not confirmation"),
        (0.50, 0.27, "Stage 2 | Prospective confirmation\n3,230 unused candidates → 500 frozen → 499 × 100 new photos\n377 evaluable; primary H2 n=158\nWHY: untouched species/photo test"),
    ]
    for x, y, text in boxes:
        ax.text(
            x,
            y,
            text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7.9,
            bbox=box_style,
        )

    for source, target in [
        ((0.50, 0.79), (0.50, 0.76)),
        ((0.50, 0.61), (0.50, 0.57)),
        ((0.50, 0.40), (0.50, 0.36)),
    ]:
        ax.annotate(
            "",
            xy=target,
            xytext=source,
            xycoords=ax.transAxes,
            arrowprops={"arrowstyle": "->", "color": "#777777", "lw": 1.2},
        )

    ax.text(
        0.50,
        0.065,
        "Stage 3 | Post-H2 annotations\nimage highlight + WorldClim 2.1 + V.PhyloMaker2\nWHY: bound measurement artifacts and alternative explanations",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=7.8,
        color=NEUTRAL,
        bbox=secondary_box_style,
    )
    ax.annotate(
        "",
        xy=(0.50, 0.125),
        xytext=(0.50, 0.17),
        xycoords=ax.transAxes,
        arrowprops={"arrowstyle": "->", "color": "#777777", "lw": 1.1, "linestyle": "--"},
    )
    ax.text(
        0.98,
        0.99,
        "Inferential sequence — cohorts are not nested samples",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.6,
        color=NEUTRAL,
    )
    panel_label(ax, "C")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure1_measurement_frame")
    meta = {
        "D_source": str(d_path.relative_to(root)),
        "discovery_n": 369,
        "reserve_n": 363,
        "global_frame_species": 42111,
        "high_depth_capable_species": 4730,
        "prospective_candidate_species": 3230,
        "third_cohort_species": 499,
        "third_cohort_measurement_evaluable_species": 377,
        "layout_contract": {
            "cohort_topology": "inferential_sequence_not_nested_samples",
            "arrow_direction": "top_to_bottom",
            "stage_necessity": "shown_in_each_stage_box",
            "secondary_followup": "dashed_post_h2_annotations",
        },
    }
    return files, meta

def figure2(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    repeated = load_json(root / "results" / "polymorphism_h1_observer_disjoint_reliability_20260913" / "result.json")
    strict = load_json(root / "results" / "polymorphism_h1_observer_disjoint_D_reliability_20260913" / "result.json")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.35))
    fig.suptitle("H1: observer-disjoint reproducibility of the polymorphism score", y=1.03, fontweight="bold")

    ax = axes[0]
    cohorts = ["discovery", "reserve"]
    # Keep the two cohorts visually close: this is an interval comparison,
    # not a two-level vertical scale.
    ys = np.array([0.42, 0.58])
    for y, cohort, colour in zip(ys, cohorts, [SECONDARY, PRIMARY], strict=True):
        s = repeated[cohort]["primary20"]
        ax.hlines(y, s["rho_q05"], s["rho_q95"], color=colour, linewidth=3)
        ax.scatter(s["rho_median"], y, s=80, color=colour, edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(s["rho_q95"] + 0.005, y, f"median {s['rho_median']:.3f}", va="center", fontsize=8.5)
    primary_floor = repeated["primary_gate"]["minimum_median_split_rho"]
    ax.axvline(primary_floor, color=NEUTRAL, linestyle="--", linewidth=1.2)
    ax.text(
        primary_floor + 0.002,
        0.50,
        "primary floor = 2/3",
        rotation=90,
        ha="left",
        va="center",
        fontsize=7.9,
        color=NEUTRAL,
    )
    reserve_q05 = float(repeated["reserve"]["primary20"]["rho_q05"])
    ax.text(
        reserve_q05,
        ys[1] + 0.045,
        f"q05={reserve_q05:.3f}",
        ha="center",
        va="top",
        fontsize=7.8,
        color=NEUTRAL,
    )
    ax.set_yticks(ys, ["Discovery", "Validation"])
    ax.set_ylim(0.68, 0.32)
    ax.set_xlim(0.63, 0.86)
    ax.set_xlabel("Observer-disjoint split Spearman rho")
    ax.set_title("200 frozen partitions: q05–median–q95")
    panel_label(ax, "A")

    ax = axes[1]
    rows = []
    for cohort in cohorts:
        s = strict[cohort]
        rows.append((cohort, s["spearman_D_A_D_B"], *s["spearman_bootstrap_95_percentile_ci"], s["lin_ccc"]))
    for y, (cohort, rho, low, high, ccc), colour in zip(ys, rows, [SECONDARY, PRIMARY], strict=True):
        ax.hlines(y, low, high, color=colour, linewidth=3)
        ax.scatter(rho, y, s=80, color=colour, edgecolor="white", linewidth=0.8, zorder=3)
        label_y = y - 0.035 if cohort == "discovery" else y + 0.035
        ax.text(
            0.878,
            label_y,
            f"rho={rho:.3f}; CCC={ccc:.3f}",
            ha="right",
            va="center",
            fontsize=8.2,
        )
    ax.axvline(0.80, color=NEUTRAL, linestyle="--", linewidth=1.3)
    ax.text(
        0.802,
        0.50,
        "strict floor = 0.80",
        rotation=90,
        ha="left",
        va="center",
        fontsize=8.0,
        color=NEUTRAL,
    )
    ax.set_yticks(ys, ["Discovery", "Validation"])
    ax.set_ylim(0.68, 0.32)
    ax.set_xlim(0.72, 0.88)
    ax.set_xlabel("Deterministic split Spearman rho")
    ax.set_title("Later prespecified stress test: bootstrap 95% CI")
    panel_label(ax, "B")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure2_h1_reproducibility")
    meta = {
        "primary_verdict": repeated["decision"]["verdict"],
        "reserve_primary_rho_median": repeated["reserve"]["primary20"]["rho_median"],
        "reserve_primary_rho_q05": repeated["reserve"]["primary20"]["rho_q05"],
        "stress_verdict": strict["decision"]["verdict"],
        "stress_reserve_rho": strict["reserve"]["spearman_D_A_D_B"],
        "stress_floor": 0.80,
        "layout_contract": {
            "stress_annotations": "offset_no_legend_overlap",
        },
    }
    return files, meta


def figure3(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    result = load_json(root / "results" / "polymorphism_white_axis_targeted_test_20260912" / "result.json")
    palette_order = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
    loadings = [float(result["fixed_axis_palette_loadings"][name]) for name in palette_order]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.7), gridspec_kw={"width_ratios": [1.0, 1.35]})
    fig.suptitle("H2 target localization in the original cohorts", y=1.02, fontweight="bold")

    ax = axes[0]
    x = np.arange(len(palette_order))
    bars = ax.bar(x, loadings, color=[PRIMARY if name == "white" else LIGHT for name in palette_order], edgecolor="#777777", linewidth=0.5)
    ax.axhline(0, color="#777777", linewidth=0.8)
    ax.set_xticks(x, palette_order, rotation=45, ha="right")
    ax.set_ylabel("Fixed q_white loading")
    ax.set_title("White versus equal mean of eight non-white coordinates")
    ax.text(
        0.98,
        0.95,
        "Named axis isolated after\noriginal broad H2 was opened",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.3,
        color=NEUTRAL,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#DDDDDD", "alpha": 0.95},
    )
    panel_label(ax, "A")

    ax = axes[1]
    entries = []
    for tier_key, tier_label in [("primary_0_10", "Primary 0.10"), ("strict_0_20", "Strict 0.20")]:
        tier = result["thresholds"][tier_key]
        for cohort in ["discovery", "reserve"]:
            s = tier[cohort]
            entries.append(
                {
                    "label": f"{tier_label}\n{("Validation" if cohort == "reserve" else "Discovery")}",
                    "observed": float(s["observed_mean_squared_white_axis_alignment"]),
                    "median": float(s["structured_null_summary"]["q50"]),
                    "low": float(s["structured_null_summary"]["q025"]),
                    "high": float(s["structured_null_summary"]["q975"]),
                    "p": float(s["structured_null_upper_p"]),
                    "n": int(s["species"]),
                }
            )
    y = np.arange(len(entries))
    for i, e in enumerate(entries):
        ax.hlines(i, e["low"], e["high"], color="#AAAAAA", linewidth=5, zorder=1)
        ax.scatter(e["median"], i, s=55, color=SECONDARY, label="Null median" if i == 0 else None, zorder=2)
        ax.scatter(e["observed"], i, s=85, marker="D", color=PRIMARY, label="Observed W" if i == 0 else None, zorder=3)
        ax.text(e["observed"] + 0.008, i, f"n={e['n']}, p={e['p']:.3f}", va="center", fontsize=8.2)
    ax.set_yticks(y, [e["label"] for e in entries])
    ax.invert_yaxis()
    ax.set_xlabel("Mean squared alignment with fixed white axis, W")
    ax.set_title("Observed alignment versus structured-null 95% interval")
    ax.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=2,
        borderaxespad=0,
    )
    panel_label(ax, "B")

    fig.tight_layout(rect=[0, 0.10, 1, 1])
    files = save_pair(fig, output_dir, "polymorphism_figure3_h2_target_localization")
    meta = {
        "status": result["status"],
        "verdict": result["decision"]["verdict"],
        "structured_null_replicates": result["structured_null_replicates"],
        "layout_contract": {
            "legend": "outside_below_axis",
        },
    }
    return files, meta


def _h2_hist_panel(ax: plt.Axes, null: np.ndarray, observed: float, species: int, p: float, title: str) -> None:
    ax.hist(null, bins=36, density=True, color=LIGHT, edgecolor="white", linewidth=0.5)
    median = float(np.median(null))
    ax.axvline(median, color=SECONDARY, linestyle="--", linewidth=2.0, label=f"Null median = {median:.3f}")
    ax.axvline(observed, color=PRIMARY, linewidth=2.6, label=f"Observed W = {observed:.3f}")
    ax.set_xlabel("Structured-null W")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.text(
        0.03,
        0.95,
        f"n = {species} species\n999 frozen null worlds\nupper-tail p = {p:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.95},
    )


def figure4(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    result_path = root / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
    result = load_json(result_path)
    base = result_path.parent
    primary_null = pd.read_csv(base / "primary_0_10_structured_null.csv")["W"].to_numpy(dtype=float)
    strict_null = pd.read_csv(base / "strict_0_20_structured_null.csv")["W"].to_numpy(dtype=float)
    if len(primary_null) != 999 or len(strict_null) != 999:
        raise ValueError("prospective H2 null arrays must each contain 999 rows")

    primary = result["thresholds"]["primary_0_10"]
    strict = result["thresholds"]["strict_0_20"]

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))
    fig.suptitle("Prospective confirmation of the pre-frozen achromatic–chromatic axis", y=1.02, fontweight="bold")
    _h2_hist_panel(
        axes[0],
        primary_null,
        float(primary["observed_W"]),
        int(primary["species"]),
        float(primary["structured_null_upper_p"]),
        "Primary admissibility tier (0.10)",
    )
    _h2_hist_panel(
        axes[1],
        strict_null,
        float(strict["observed_W"]),
        int(strict["species"]),
        float(strict["structured_null_upper_p"]),
        "Strict admissibility tier (0.20)",
    )
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    # The same-universe / non-independent-source boundary belongs in the
    # manuscript caption; keeping it out of the plotting field improves
    # legibility after journal-width reduction.
    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure4_prospective_h2")
    meta = {
        "verdict": result["decision"]["verdict"],
        "support_evaluable_species": int(result["measurement_support"]["evaluable_species"]),
        "primary": {
            "species": int(primary["species"]),
            "W": float(primary["observed_W"]),
            "null_median": float(primary["structured_null_summary"]["q50"]),
            "p": float(primary["structured_null_upper_p"]),
        },
        "strict": {
            "species": int(strict["species"]),
            "W": float(strict["observed_W"]),
            "null_median": float(strict["structured_null_summary"]["q50"]),
            "p": float(strict["structured_null_upper_p"]),
        },
        "layout_contract": {
            "same_universe_nonreplication_note": "caption_not_plot_field",
        },
    }
    return files, meta


def verify_five_cohort_original_reporting(root: Path) -> dict:
    """Figure reporting only: authenticate the final original-photo extension."""
    result = load_json(root / "results" / "fcp_original_447500_five_cohort_verified_reporting_20261018.json")
    source = result.get("source", {})
    frame = result.get("study", {})
    spatial = result.get("spatial", {})
    if (result.get("schema") != "fcp_five_cohort_real_447500_original_method_reporting_receipt_v1"
        or source.get("original_pixel_roi_run_id") != 38997169400
        or source.get("source_aggregation_run_id") != 39049598596
        or source.get("measured_archive_sha256") !=
            "99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a8de07edd38341b27889310bc25"
        or frame.get("genuine_independent_100photo_species") != 4475
        or frame.get("actual_original_photo_terminals") != 447500
        or frame.get("four_state_classifiable_photo_rows") != 227107
        or frame.get("min_40_classifiable_species") != 3297
        or spatial.get("primary_radius_post_hoc") is not True
        or spatial.get("original_within_species_vertex_null_replicates") != 199
        or spatial.get("pooled_50km", {}).get("n_species") != 1761
        or abs(spatial.get("pooled_50km", {}).get("mean_depletion", -1) -
               0.017614902168541096) > 1e-12
        or spatial.get("pooled_50km", {}).get("mean_depletion_p_upper") != 0.005):
        raise ValueError("final five-cohort FCP actual original-photo evidence identity drift")
    expected = {
        "discovery": (166, 0.020529254583812922, 0.005),
        "validation": (181, 0.018672971642749295, 0.005),
        "prospective": (204, 0.01468491968437185, 0.005),
        "wave4": (190, 0.017293596294403593, 0.005),
        "new_U100": (1020, 0.017598697724511534, 0.005),
    }
    for cohort, (n, depletion, p) in expected.items():
        row = spatial["all_sources"][cohort]["50"]
        if (row.get("n_species") != n
            or abs(row.get("mean_depletion", -1) - depletion) > 1e-12
            or row.get("p_upper") != p):
            raise ValueError("frozen fifth-cohort source effect changed: " + cohort)
    return result


def figure5(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    distributed = load_json(root / "results" / "polymorphism_distributed_polymorphism_posthoc_20261007" / "result.json")
    robustness = load_json(root / "results" / "polymorphism_distributed_polymorphism_robustness_20261007" / "result.json")
    ibd_ibe = load_json(root / "results" / "polymorphism_phenotypic_IBD_IBE_posthoc_20261007" / "result.json")

    cohorts = ["discovery", "validation", "third"]
    labels = ["Discovery", "Validation", "Third cohort"]
    five = verify_five_cohort_original_reporting(root)

    fig, axes = plt.subplots(1, 3, figsize=(17.1, 5.1),
                             gridspec_kw={"width_ratios": [1.50, 1.25, 1.0]})
    fig.suptitle("Species-wide flower-colour variation is geographically partitioned",
                 y=1.01, fontweight="bold")

    # Panel A: five independent actual-pixel species cohorts and one pooled mean.
    # The pooled bar must NEVER be labelled or interpreted as a sixth cohort.
    ax = axes[0]
    final_groups = ("discovery", "validation", "prospective", "wave4", "new_U100")
    final_labels = ["Discovery", "Validation", "Prospective", "Wave 4", "New U100", "Pooled (not cohort)"]
    final_rows = [five["spatial"]["all_sources"][co]["50"] for co in final_groups]
    final_rows.append(five["spatial"]["pooled_50km"])
    effect = np.array([float(r["mean_depletion"]) for r in final_rows])
    axis_y = np.arange(len(final_labels))
    colours = [SECONDARY, PRIMARY, SUPPORT, "#9673B9", "#9A6324", "#444444"]
    ax.barh(axis_y, effect, color=colours, height=0.62)
    ax.axvline(0, color="#777777", linewidth=1.0)
    for yy, rr, val in zip(axis_y, final_rows, effect, strict=True):
        ax.text(val+0.00045, yy, f"n={int(rr['n_species']):,}; p={float(rr['mean_depletion_p_upper'] if 'mean_depletion_p_upper' in rr else rr['p_upper']):.3f}",
                ha="left", va="center", fontsize=7.9)
    ax.set_yticks(axis_y, final_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, max(effect)*1.98)
    ax.set_xlabel("Global minus local pair dissimilarity")
    ax.set_title("Five same-provider cohorts at 50 km (post hoc)")
    ax.text(0.015, 0.015, "Pooled = species-equal mean; not a new cohort",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.4, color=NEUTRAL)
    panel_label(ax, "A")

    # Panel B: falsification tests.
    ax = axes[1]
    tests = [
        ("Different\nobserver", "R1_different_observer"),
        ("Nonwhite\nonly", "R2_nonwhite_only"),
        ("Continuous\n9-colour", "R3_continuous_nine_colour"),
    ]
    offsets = [-0.22, 0.0, 0.22]
    for offset, co, label, colour in zip(offsets, cohorts, labels, [SECONDARY, PRIMARY, SUPPORT], strict=True):
        vals = [float(robustness[co][key]["mean_depletion"]) for _, key in tests]
        xs = np.arange(len(tests)) + offset
        ax.scatter(xs, vals, s=70, color=colour, edgecolor="white", linewidth=0.8, label=label, zorder=3)
        for xx, val, (_, key) in zip(xs, vals, tests, strict=True):
            pval = float(robustness[co][key]["p_upper"])
            ax.text(xx, val + 0.0012, f"{pval:.3f}", ha="center", va="bottom", fontsize=6.9)
    ax.axhline(0, color="#777777", linewidth=1.0)
    ax.set_xticks(np.arange(len(tests)), [t[0] for t in tests])
    ax.set_ylabel("Mean local depletion")
    ax.set_ylim(0, 0.0255)
    ax.set_title("Partitioning survives measurement falsifications")
    ax.legend(frameon=False, loc="upper right", fontsize=7.2)
    ax.text(0.03, 0.04, "Labels above points are matched-null p", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7.2, color=NEUTRAL)
    panel_label(ax, "B")

    # Panel C: phenotypic IBD versus IBE-like structure.
    ax = axes[2]
    width = 0.34
    x = np.arange(3)
    ibd = np.array([float(ibd_ibe["cohorts"][co]["IBD_like"]["mean_partial_rho"]) for co in cohorts])
    ibe = np.array([float(ibd_ibe["cohorts"][co]["IBE_like"]["mean_partial_rho"]) for co in cohorts])
    ax.bar(x - width/2, ibd, width=width, color=PRIMARY, label="IBD-like")
    ax.bar(x + width/2, ibe, width=width, color=SECONDARY, label="BIO5 IBE-like")
    for xi, co, a, b in zip(x, cohorts, ibd, ibe, strict=True):
        p_ibd = float(ibd_ibe["cohorts"][co]["IBD_like"]["p"])
        p_ibe = float(ibd_ibe["cohorts"][co]["IBE_like"]["p"])
        ax.text(xi-width/2, a+0.0010, f"p={p_ibd:.3f}", ha="center", va="bottom", fontsize=6.7, rotation=90)
        ax.text(xi+width/2, b+0.0010, f"p={p_ibe:.3f}", ha="center", va="bottom", fontsize=6.7, rotation=90)
    ax.axhline(0, color="#777777", linewidth=1.0)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Equal-species mean partial rho")
    ax.set_ylim(0, max(ibd) * 1.38)
    ax.set_title("Geographic distance dominates a smaller BIO5 residual")
    ax.legend(frameon=False, loc="upper right", fontsize=7.2)
    panel_label(ax, "C")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure5_explanatory_boundaries")
    meta = {
        "distributed_polymorphism": {
            "receipt": "results/polymorphism_distributed_polymorphism_posthoc_20261007/result.json",
            "primary_radius_km": 50,
            "discovery_depletion": float(distributed["primary"]["discovery"]["mean_local_depletion"]),
            "discovery_p": float(distributed["primary"]["discovery"]["p_upper"]),
            "validation_depletion": float(distributed["primary"]["validation"]["mean_local_depletion"]),
            "validation_p": float(distributed["primary"]["validation"]["p_upper"]),
            "third_depletion": float(distributed["primary"]["third"]["mean_local_depletion"]),
            "third_p": float(distributed["primary"]["third"]["p_upper"]),
            "multiscale_all_three": bool(all(distributed["multiscale_local_depletion_supported_all_three"].values())),
        },
        "five_cohort_original_method": {
            "receipt": "results/fcp_original_447500_five_cohort_verified_reporting_20261018.json",
            "measured_species": 4475,
            "measured_photo_rows": 447500,
            "frozen_posthoc_radius_km": 50,
            "spatial_eligible_species": 1761,
            "pooled_mean_depletion": float(five["spatial"]["pooled_50km"]["mean_depletion"]),
            "pooled_p_upper": float(five["spatial"]["pooled_50km"]["mean_depletion_p_upper"]),
            "cohorts": {
                name: {
                    "n_species": int(five["spatial"]["all_sources"][name]["50"]["n_species"]),
                    "mean_depletion": float(five["spatial"]["all_sources"][name]["50"]["mean_depletion"]),
                    "p_upper": float(five["spatial"]["all_sources"][name]["50"]["p_upper"]),
                } for name in ("discovery", "validation", "prospective", "wave4", "new_U100")
            },
            "pooled_not_sixth_cohort": True,
            "falsification_only_original_three": True,
            "ibd_ibe_only_original_three": True,
        },
        "robustness": {
            "receipt": "results/polymorphism_distributed_polymorphism_robustness_20261007/result.json",
            "different_observer_all_three": bool(robustness["cross_cohort"]["different_observer_supported_all_three"]),
            "nonwhite_only_all_three": bool(robustness["cross_cohort"]["nonwhite_only_supported_all_three"]),
            "continuous_nine_colour_all_three": bool(robustness["cross_cohort"]["continuous_nine_colour_supported_all_three"]),
        },
        "ibd_ibe": {
            "receipt": "results/polymorphism_phenotypic_IBD_IBE_posthoc_20261007/result.json",
            "IBE_like_all_three": bool(ibd_ibe["replication"]["IBE_like_all_three"]),
            "IBD_like_all_three": bool(ibd_ibe["replication"]["IBD_like_all_three"]),
            "IBE_stronger_than_IBD_all_three": bool(ibd_ibe["replication"]["IBE_stronger_than_IBD_all_three"]),
            "cohorts": {
                co: {
                    "IBD_like": float(ibd_ibe["cohorts"][co]["IBD_like"]["mean_partial_rho"]),
                    "IBE_like": float(ibd_ibe["cohorts"][co]["IBE_like"]["mean_partial_rho"]),
                }
                for co in cohorts
            },
        },
        "layout_contract": {
            "panel_A": "five_real_source_cohorts_plus_species_equal_pooled_50km_depletion",
            "panel_B": "three_falsification_tests_across_three_cohorts",
            "panel_C": "paired_IBD_IBE_effects",
        },
    }
    return files, meta



def figure6(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    """Show FCP commonality, not global proportions of particular flower hues.

    Source is a SHA-linked summary of 447500 genuinely measured original photo
    terminals, conditionally sampled rather than a botanical population census.
    The figure is entirely post hoc and does not modify the original null/W.
    """
    source = load_json(
        root / "results" / "fcp_original_447500_polymorphism_commonality_20261018.json"
    )
    if source.get("schema") != "fcp_verified_original_447500_intraspecific_polymorphism_commonality_summary_v1":
        raise ValueError("missing source-verified five-cohort FCP commonality receipt")
    sampling = source["sampling"]
    comp = source["composition"]
    structure = source["space"]
    checks = source["falsification"]
    if (sampling["genuinely_measured_taxa"] != 4475 or
        sampling["species_with_ge40_photographs"] != 3297 or
        structure["original_four_state"]["n"] != 1761 or
        abs(structure["original_four_state"]["effect"] - 0.017614902168541096) > 1e-10 or
        abs(structure["species_majority_versus_other"]["effect"] - 0.016485145475518544) > 1e-10):
        raise ValueError("FCP commonality figure input does not match authentic original study")

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.6))
    fig.suptitle(
        "Within-species floral-colour polymorphism: skew, allocation and limits",
        y=1.005, fontweight="bold",
    )

    # A: State number measured per species, NOT incidence of biological morphs.
    ax = axes[0, 0]
    states = comp["observed_coarse_state_counts"]
    xs = np.arange(1, 5)
    ys = np.array([int(states[str(x)]) for x in xs])
    ax.bar(xs, ys, color=[NEUTRAL, SECONDARY, SUPPORT, PRIMARY], width=.67)
    for x, v in zip(xs, ys, strict=True):
        ax.text(x, v+35, f"{v:,}", ha="center", fontsize=9)
    ax.set_xticks(xs)
    ax.set_xlabel("Number of photographed coarse states / species")
    ax.set_ylabel("Species with >=40 classifiable photos")
    ax.set_ylim(0, max(ys)*1.20)
    ax.text(.03, .94,
            f"Median species D = {comp['median_Gini_Simpson_D']:.3f}\n"
            f"Median dominant-state photo share = "
            f"{100*comp['median_majority_photo_state_fraction']:.1f}%",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox={"facecolor":"white","alpha":.94,"edgecolor":"none","pad":5})
    ax.set_title("Uneven within-species photographed morph composition")
    panel_label(ax, "A")

    # B: Group by whichever photographed morph is modal within that species.
    ax = axes[0, 1]
    groups = structure["dominant_state_categories"]
    order = ("white", "yellow_orange", "red_pink", "blue_purple")
    label = {"white":"White","yellow_orange":"Yellow/orange",
             "red_pink":"Red/pink","blue_purple":"Blue/purple"}
    colors = ["#999999", "#E2AC43", "#C95D7B", "#7466A5"]
    ordered = [next(t for t in groups if t["state"]==k) for k in order]
    p = np.array([t["effect"] for t in ordered])
    ax.bar(np.arange(4), p, color=colors, width=.68)
    for x, row in enumerate(ordered):
        ax.text(x, row["effect"]+.00125, f"n={row['n']}\np={row['p_upper']:.3f}",
                ha="center", fontsize=8.0)
    ax.axhline(0, color=NEUTRAL, lw=.9)
    ax.set_xticks(np.arange(4), [label[k] for k in order])
    ax.set_ylabel("Species-wide minus within-50-km dissimilarity")
    ax.set_ylim(0, .037)
    ax.set_title("Geographic allocation spans all photographed majorities")
    panel_label(ax, "B")

    # C: Effect changes with the amount of within-species D, and is not
    # a uniform biological species-level coefficient.
    ax = axes[1, 0]
    quartiles = structure["diversity_quartiles"]
    effects = [q["effect"] for q in quartiles]
    ax.bar(np.arange(4), effects,
           color=[LIGHT, "#AFC7DA", SECONDARY, PRIMARY], width=.7)
    for x, q in enumerate(quartiles):
        ax.text(x, q["effect"]+.0014, f"n={q['n']}\np={q['p_upper']:.3f}",
                ha="center", fontsize=8)
    ax.set_xticks(np.arange(4), ["Lowest\nD quartile", "Q2", "Q3", "Highest\nD quartile"])
    ax.set_ylabel("Mean original four-state 50-km depletion")
    ax.set_ylim(0, max(effects)*1.32)
    ax.set_title("Stronger raw depletion where species D is larger")
    panel_label(ax, "C")

    # D: Every point has its own eligibility and sampling denominator.
    ax = axes[1, 1]
    records = [
        ("Original 4-state", structure["original_four_state"], SECONDARY),
        ("Majority vs alternatives", structure["species_majority_versus_other"], SUPPORT),
        ("Equal-weight genera", {"effect":structure["genus_balanced_four_state"]["effect"],
                                "n":structure["genus_balanced_four_state"]["n_genera"],
                                "p_upper":structure["genus_balanced_four_state"]["p_upper"]}, PRIMARY),
        ("Different observers", checks["different_observer"], SECONDARY),
        ("Different years", checks["different_year"], SECONDARY),
        ("No white records", checks["nonwhite_only_pooled"], NEUTRAL),
    ]
    for i, (_, result, color) in enumerate(records):
        ax.barh(i, result["effect"], color=color, height=.67)
        ax.text(result["effect"] + .0007, i,
                f"{result['effect']:.4f}  n={result['n']}",
                va="center", fontsize=8)
    ax.set_yticks(range(len(records)), [z[0] for z in records])
    ax.invert_yaxis()
    ax.axvline(0, lw=.9, color=NEUTRAL)
    ax.set_xlim(0, .027)
    ax.set_xlabel("Mean photographed local depletion (50 km)")
    ax.set_title("Same-source robustness, with nonwhite exceptions")
    ax.text(.98, .96, "Nonwhite p=0.14 in prospective;\np=0.72 in fourth cohort",
            transform=ax.transAxes, fontsize=7.4, ha="right", va="top")
    panel_label(ax, "D")

    fig.tight_layout(h_pad=2.9, w_pad=2.5)
    files = save_pair(fig, output_dir, "polymorphism_figure6_fcp_commonality")
    meta = {
        "source":"results/fcp_original_447500_polymorphism_commonality_20261018.json",
        "actual_measured_photos":sampling["photographs"],
        "descriptively_eligible_species":sampling["species_with_ge40_photographs"],
        "spatially_eligible_species":structure["original_four_state"]["n"],
        "majority_fraction_median":comp["median_majority_photo_state_fraction"],
        "four_state_50km_depletion":structure["original_four_state"]["effect"],
        "majority_binary_50km_depletion":structure["species_majority_versus_other"]["effect"],
        "original_null_worlds":structure["original_complete_vertex_null_replicates"],
        "post_hoc":True,
        "nonwhite_support_not_all_cohorts":True,
        "not_true_population_morph_frequency":True,
        "files":files,
    }
    return files, meta


def supplementary_figure9(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    env = load_json(root / "results" / "polymorphism_white_environment_mechanism_20260925" / "result.json")
    transport = load_json(root / "results" / "polymorphism_legacy_white_bio5_replication_20260925" / "result.json")
    observer = load_json(root / "results" / "polymorphism_white_environment_observer_sensitivity_20260925" / "result.json")
    molecular = load_json(root / "results" / "polymorphism_silene_molecular_anchor_20260928" / "result.json")
    silene = load_json(root / "results" / "polymorphism_silene_decoupling_persistence_20260925" / "result.json")
    cross = load_json(root / "results" / "polymorphism_crossspecies_pal_wal_frequency_20260925" / "result.json")

    bio5 = next(x for x in env["results"] if x["variable"] == "bio5")
    paired_bio5 = next(x for x in observer["paired_results"] if x["variable"] == "bio5")
    balanced_bio5 = next(x for x in observer["observer_balanced_results"] if x["variable"] == "bio5")

    fig, axes = plt.subplots(1, 3, figsize=(17.2, 4.8), gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})
    fig.suptitle("Secondary molecular, persistence and environmental evidence", y=1.02, fontweight="bold")

    ax = axes[0]
    labels = ["Prospective\nprimary", "Observer-\nbalanced", "Observer-\npaired", "Discovery\ntransport", "Validation\ntransport"]
    deltas = [
        float(bio5["median_delta_white_minus_nonwhite_SD"]),
        float(balanced_bio5["median_delta"]),
        float(paired_bio5["median_delta"]),
        float(transport["discovery"]["species_level"]["median_delta_white_minus_nonwhite_SD"]),
        float(transport["reserve"]["species_level"]["median_delta_white_minus_nonwhite_SD"]),
    ]
    pvals = [
        float(bio5["wilcoxon_holm_p"]),
        float(balanced_bio5["wilcoxon_two_sided_p"]),
        float(paired_bio5["wilcoxon_two_sided_p"]),
        float(transport["discovery"]["species_level"]["wilcoxon_two_sided_p"]),
        float(transport["reserve"]["species_level"]["wilcoxon_two_sided_p"]),
    ]
    ns = [
        int(bio5["n_species_delta"]),
        int(balanced_bio5["n_species"]),
        int(paired_bio5["n_species"]),
        int(transport["discovery"]["eligible_species"]),
        int(transport["reserve"]["eligible_species"]),
    ]
    colours = [SUPPORT, SECONDARY, NEUTRAL, LIGHT, LIGHT]
    x = np.arange(5)
    ax.axhline(0, color="#777777", linewidth=1.0)
    ax.bar(x, deltas, color=colours, width=0.62)
    for xi, delta, pval, n in zip(x, deltas, pvals, ns, strict=True):
        offset = 0.007 if delta >= 0 else -0.007
        va = "bottom" if delta >= 0 else "top"
        ax.text(xi, delta + offset, f"n={n}\np={pval:.3f}", ha="center", va=va, fontsize=7.4)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Median within-species white − nonwhite BIO5 (SD)")
    ax.set_ylim(-0.045, 0.105)
    ax.set_title("BIO5 is positive in the frozen primary analysis\nbut observer-sensitive and non-transporting")
    ax.text(
        0.02,
        0.04,
        f"Prospective conditional OR = {float(bio5['OR_per_within_species_SD']):.3f}\n"
        f"Holm-adjusted p = {float(bio5['wilcoxon_holm_p']):.4f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=NEUTRAL,
    )
    panel_label(ax, "A")

    ax = axes[1]
    molecular_rows = [
        ("F3h1\nD/W", molecular["bud_expression"]["dark_vs_white_significant"][0]),
        ("C4h2\nD/W", molecular["bud_expression"]["dark_vs_white_significant"][1]),
        ("Myb1a\nD/W", molecular["bud_expression"]["dark_vs_white_significant"][2]),
        ("F3h1\nL/W", molecular["bud_expression"]["light_vs_white_significant"][0]),
        ("F3′h\nL/W", molecular["bud_expression"]["light_vs_white_significant"][1]),
    ]
    x = np.arange(len(molecular_rows))
    fold = np.array([float(row["fold_change_pigmented_over_white"]) for _, row in molecular_rows])
    p_mol = np.array([float(row["p"]) for _, row in molecular_rows])
    bar_colours = [PRIMARY, SECONDARY, SUPPORT, PRIMARY, SECONDARY]
    ax.bar(x, fold, color=bar_colours, width=0.65)
    for xi, value, pval in zip(x, fold, p_mol, strict=True):
        ax.text(xi, value * 1.08, f"{value:g}×\np={pval:.3f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x, [label for label, _ in molecular_rows])
    ax.set_ylabel("Published bud-stage fold change\n(pigmented / white)")
    ax.set_yscale("log")
    ax.set_ylim(3, 90)
    ax.set_title("Published Silene RNA-seq localizes a\npetal regulatory anchor near F3h1/Myb1a")
    ax.text(
        0.02,
        0.04,
        "All source-reported p < 0.05\npigmented-vs-white bud contrasts",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color=NEUTRAL,
    )
    panel_label(ax, "B")

    ax = axes[2]
    x = np.array([0, 1, 2.4, 3.4])
    y = np.array([
        float(silene["phenotypes"]["PAL"]["positive_frequency_median_percent"]),
        float(silene["phenotypes"]["WAL"]["positive_frequency_median_percent"]),
        float(cross["PAL"]["lower_bound_median_percent"]),
        float(cross["WAL"]["numeric_upper_bound_median_percent"]),
    ])
    labels = ["Silene\nPAL", "Silene\nWAL", "Cross-system\nPAL lower", "Cross-system\nWAL upper"]
    colours = [PRIMARY, SECONDARY, PRIMARY, SECONDARY]
    ax.scatter(x, y, s=90, color=colours, edgecolor="white", linewidth=0.8, zorder=3)

    pal_min = float(silene["phenotypes"]["PAL"]["positive_frequency_min_percent"])
    pal_max = float(silene["phenotypes"]["PAL"]["positive_frequency_max_percent"])
    wal_min = float(silene["phenotypes"]["WAL"]["positive_frequency_min_percent"])
    wal_max = float(silene["phenotypes"]["WAL"]["positive_frequency_max_percent"])
    ax.vlines(0, pal_min, pal_max, color=PRIMARY, linewidth=3, alpha=0.65)
    ax.vlines(1, wal_min, wal_max, color=SECONDARY, linewidth=3, alpha=0.65)

    wal_cross_max = float(cross["WAL"]["numeric_upper_bound_max_percent"])
    ax.axhline(wal_cross_max, color=NEUTRAL, linestyle="--", linewidth=1.0)
    ax.text(3.55, wal_cross_max, f"max quantified WAL upper bound = {wal_cross_max:.1f}%", ha="right", va="bottom", fontsize=7.8, color=NEUTRAL)

    for xi, yi in zip(x, y, strict=True):
        ax.text(xi, yi * 1.25, f"{yi:g}%", ha="center", va="bottom", fontsize=8.3)

    ax.set_yscale("log")
    ax.set_ylim(0.03, 60)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Reported natural frequency (%)")
    ax.set_title("Biochemically anchored flower-restricted loss reaches\nhigher reported natural frequencies")
    ax.text(
        0.02,
        0.04,
        "Silene PAL/WAL anchor: HPLC-DAD-MS$^n$ tissue profiles;\ncross-system entries retain source-table classifications",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=NEUTRAL,
    )
    panel_label(ax, "B")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figureS9_secondary_mechanism_evidence")
    meta = {
        "bio5": {
            "prospective_median_delta_SD": float(bio5["median_delta_white_minus_nonwhite_SD"]),
            "prospective_holm_p": float(bio5["wilcoxon_holm_p"]),
            "prospective_conditional_OR": float(bio5["OR_per_within_species_SD"]),
            "observer_balanced_median_delta_SD": float(balanced_bio5["median_delta"]),
            "observer_balanced_wilcoxon_p": float(balanced_bio5["wilcoxon_two_sided_p"]),
            "observer_paired_median_delta_SD": float(paired_bio5["median_delta"]),
            "observer_paired_wilcoxon_p": float(paired_bio5["wilcoxon_two_sided_p"]),
            "observer_paired_conditional_OR": float(paired_bio5["OR_per_within_species_SD"]),
            "discovery_median_delta_SD": float(transport["discovery"]["species_level"]["median_delta_white_minus_nonwhite_SD"]),
            "discovery_p": float(transport["discovery"]["species_level"]["wilcoxon_two_sided_p"]),
            "validation_median_delta_SD": float(transport["reserve"]["species_level"]["median_delta_white_minus_nonwhite_SD"]),
            "validation_p": float(transport["reserve"]["species_level"]["wilcoxon_two_sided_p"]),
        },
        "molecular_anchor": {
            "source_doi": molecular["source"]["doi"],
            "F3h1_dark_vs_white_fold": float(molecular["bud_expression"]["shared_F3h1_fold_changes"]["dark_vs_white"]),
            "F3h1_light_vs_white_fold": float(molecular["bud_expression"]["shared_F3h1_fold_changes"]["light_vs_white"]),
            "F3h1_shared_significant": bool(molecular["derived_synthesis"]["F3h1_is_only_locus_significant_in_both_pigmented_vs_white_bud_contrasts"]),
            "Myb1a_dark_vs_white_fold": float(molecular["bud_expression"]["dark_vs_white_significant"][2]["fold_change_pigmented_over_white"]),
            "expanded_sequence_survey_individuals": int(molecular["sequence_evidence"]["expanded_sequence_survey_individuals"]),
            "causal_variant_identified": bool(molecular["derived_synthesis"]["causal_variant_identified"]),
            "role": molecular["role"],
        },
        "pal_wal": {
            "silene_PAL_median_percent": float(silene["phenotypes"]["PAL"]["positive_frequency_median_percent"]),
            "silene_WAL_median_percent": float(silene["phenotypes"]["WAL"]["positive_frequency_median_percent"]),
            "cross_PAL_lower_bound_median_percent": float(cross["PAL"]["lower_bound_median_percent"]),
            "cross_WAL_upper_bound_median_percent": float(cross["WAL"]["numeric_upper_bound_median_percent"]),
            "cross_WAL_max_upper_bound_percent": float(cross["WAL"]["numeric_upper_bound_max_percent"]),
            "phenotype_anchor": "silene_HPLC_DAD_MSn_tissue_profiles_crosssystem_source_classified",
        },
        "claim_boundary": "bounded_secondary_empirical_evidence_not_universal_causation",
    }
    return files, meta

def generate_all(root: Path, output_dir: Path) -> Path:
    root = Path(root).resolve()
    output_dir = Path(output_dir).resolve()
    configure_matplotlib()
    source_manifest = validate_reporting_source(root)

    f1_files, f1_meta = figure1(root, output_dir)
    f2_files, f2_meta = figure2(root, output_dir)
    f3_files, f3_meta = figure3(root, output_dir)
    f4_files, f4_meta = figure4(root, output_dir)
    f5_files, f5_meta = figure5(root, output_dir)
    f6_files, f6_meta = figure6(root, output_dir)
    s9_files, s9_meta = supplementary_figure9(root, output_dir)

    manifest = {
        "schema": "polymorphism_manuscript_figure_manifest_v1",
        "date_jst": "2026-09-18",
        "status": "generated_from_frozen_results",
        "scientific_claims_changed": False,
        "reporting_source_manifest_sha256": sha256(
            root / "results" / "polymorphism_publication_figure_source_20260918" / "manifest.json"
        ),
        "source_artifacts": source_manifest["source_artifacts"],
        "figures": {
            "figure1": {**f1_meta, "files": f1_files},
            "figure2": {**f2_meta, "files": f2_files},
            "figure3": {**f3_meta, "files": f3_files},
            "figure4": {**f4_meta, "files": f4_files},
            "figure5": {**f5_meta, "files": f5_files},
            "figure6": {**f6_meta, "files": f6_files},
            "supplementary_figure9": {**s9_meta, "files": s9_files},
        },
        "hard_nonclaims": [
            "high-depth cohorts do not estimate global flower-colour polymorphism prevalence",
            "prospective H2 is species-disjoint within the same iNaturalist opportunity universe, not independent-source replication",
            "the white/nonwhite axis does not identify pigment chemistry, transition direction, or adaptive mechanism",
            "the D-spatial association is structural and does not identify a causal maintenance mechanism",
            "H3a non-support does not imply absence of all phylogenetic structure",
            "sampled photographic span is not true biological range size",
        ],
    }
    manifest_path = output_dir / "polymorphism_figure_manifest_20260918.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def regenerate_figure5_only(root: Path, output_dir: Path) -> Path:
    """Refresh Figure 5 framing while preserving the frozen numerical bundle."""
    root = Path(root).resolve()
    output_dir = Path(output_dir).resolve()
    configure_matplotlib()
    validate_reporting_source(root)

    manifest_path = output_dir / "polymorphism_figure_manifest_20260918.json"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    manifest = load_json(manifest_path)
    if manifest.get("schema") != "polymorphism_manuscript_figure_manifest_v1":
        raise ValueError("unexpected canonical figure-manifest schema")

    files, meta = figure5(root, output_dir)
    manifest["figures"]["figure5"] = {**meta, "files": files}
    manifest["scientific_claims_changed"] = False
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="repository root",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/figures/polymorphism_20260918"),
    )
    parser.add_argument(
        "--figure5-only",
        action="store_true",
        help="regenerate only Figure 5 and update its canonical manifest entry",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.figure5_only:
        manifest = regenerate_figure5_only(args.root, args.output_dir)
    else:
        manifest = generate_all(args.root, args.output_dir)
    print(f"Wrote frozen polymorphism publication figures and manifest: {manifest}")


if __name__ == "__main__":
    main()
