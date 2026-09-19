#!/usr/bin/env python3
"""Generate the canonical disttrait methods-paper displays from frozen receipts.

Reporting only:
- reads frozen result JSON/CSV files;
- never reruns simulation worlds or empirical inference;
- writes figures, Table 1 and a source manifest.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "figures" / "disttrait_methods_20260919"

SOURCES = {
    "conditioning": ROOT / "results" / "disttrait_species_conditioning_benchmark_v0_2_20260918" / "result.json",
    "performance": ROOT / "results" / "disttrait_performance_surface_v0_3_20260919" / "cells.csv",
    "mnar": ROOT / "results" / "disttrait_mnar_observation_v0_10_20260919" / "cells.csv",
    "model": ROOT / "results" / "disttrait_model_comparator_surface_v0_5_20260919" / "cells.csv",
    "direction": ROOT / "results" / "disttrait_direction_heterogeneity_v0_7_20260919" / "cells.csv",
    "meta": ROOT / "results" / "disttrait_random_slope_meta_v0_8_20260919" / "cells.csv",
    "nonlinear": ROOT / "results" / "disttrait_nonlinear_curvature_v0_11_20260919" / "cells.csv",
    "multivariate": ROOT / "results" / "disttrait_multivariate_orientation_v0_12_20260919" / "cells.csv",
    "sf": ROOT / "results" / "disttrait_sf_street_tree_empirical_v0_6_20260919" / "result.json",
    "gammarus": ROOT / "results" / "disttrait_gammarus_empirical_v0_9_20260919" / "result.json",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(key: str) -> dict:
    return json.loads(SOURCES[key].read_text(encoding="utf-8"))


def load_csv(key: str) -> pd.DataFrame:
    return pd.read_csv(SOURCES[key])


def save_figure(fig: plt.Figure, stem: str) -> list[str]:
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    fig.savefig(png, dpi=240, bbox_inches="tight")
    fig.savefig(
        pdf,
        bbox_inches="tight",
        metadata={
            "Title": stem,
            "Creator": "disttrait reporting-only figure generator",
        },
    )
    plt.close(fig)
    return [png.name, pdf.name]


def figure1() -> list[str]:
    result = load_json("conditioning")
    labels = ["Null", "Signal"]
    naive = [
        result["null"]["naive_pooled_false_positive_fraction"],
        result["signal"]["naive_pooled_detection_fraction"],
    ]
    conditioned = [
        result["null"]["species_conditioned_false_positive_fraction"],
        result["signal"]["species_conditioned_detection_fraction"],
    ]

    x = np.arange(len(labels))
    width = 0.34
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.bar(x - width / 2, naive, width, label="Naive pooled")
    ax.bar(x + width / 2, conditioned, width, label="Species-conditioned")
    ax.axhline(0.05, linewidth=1, linestyle="--", label="0.05")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Rejection / detection fraction")
    ax.set_ylim(0, 1.08)
    ax.set_title("Pooling confounds between-species turnover with within-species organization")
    ax.legend(frameon=False)
    return save_figure(fig, "disttrait_figure1_pooling")


def figure2() -> list[str]:
    perf = load_csv("performance")
    mnar = load_csv("mnar")

    effect_cols = ["Effect 0", "Effect 0.8", "Effect 1.6", "Effect 2.4"]
    selection_cols = ["Selection 0.8", "Selection 1.6"]
    columns = effect_cols + selection_cols

    perf_rows = []
    perf_labels = []
    for imbalance in (1.0, 4.0):
        for missing in (0.0, 0.25, 0.5):
            values = []
            for effect in (0.0, 0.8, 1.6, 2.4):
                row = perf[
                    (perf["effect_size"] == effect)
                    & (perf["imbalance_ratio"] == imbalance)
                    & (perf["missing_fraction"] == missing)
                ].iloc[0]
                values.append(float(row["conditioned_detection"]))
            perf_rows.append(values + [np.nan, np.nan])
            perf_labels.append(f"{int(imbalance)}x, MCAR {int(missing * 100)}%")

    mnar_rows = []
    mnar_labels = []
    mechanism_order = ("mcar", "trait_only", "position_only", "joint_trait_position")
    display_name = {
        "mcar": "MCAR",
        "trait_only": "Trait-only",
        "position_only": "Position-only",
        "joint_trait_position": "Trait x position",
    }
    for mechanism in mechanism_order:
        values = []
        for strength in (0.8, 1.6):
            row = mnar[
                (mnar["mechanism"] == mechanism)
                & (mnar["strength"] == strength)
            ].iloc[0]
            values.append(float(row["equal_rejection_fraction"]))
        mnar_rows.append([np.nan] * 4 + values)
        mnar_labels.append(display_name[mechanism])

    matrix = np.asarray(perf_rows + mnar_rows, dtype=float)
    labels = perf_labels + mnar_labels

    fig, ax = plt.subplots(figsize=(9.0, 6.2))
    image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(columns)), columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.set_title("Calibration, power and the observation-process boundary")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center")
    fig.colorbar(image, ax=ax, label="Rejection / detection fraction")
    return save_figure(fig, "disttrait_figure2_calibration_boundary")


def _cell(
    frame: pd.DataFrame,
    **conditions: float | str,
) -> pd.Series:
    keep = np.ones(len(frame), dtype=bool)
    for column, value in conditions.items():
        keep &= frame[column] == value
    rows = frame.loc[keep]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one row for {conditions}, found {len(rows)}")
    return rows.iloc[0]


def figure3() -> list[str]:
    model = load_csv("model")
    direction = load_csv("direction")
    nonlinear = load_csv("nonlinear")
    multivariate = load_csv("multivariate")

    scenarios = [
        "Binary\nshared",
        "Linear\nshared",
        "Linear\n50% reverse",
        "Quadratic\nshared",
        "Quadratic\n50% reverse",
        "Multivariate\nshared",
        "Multivariate\nfull spread",
    ]

    m0 = _cell(model, effect_size=0.8, imbalance_ratio=1.0, missing_fraction=0.0)
    d0 = _cell(direction, effect_size=0.4, reversal_fraction=0.0, missing_fraction=0.0)
    d5 = _cell(direction, effect_size=0.4, reversal_fraction=0.5, missing_fraction=0.0)
    q0 = _cell(nonlinear, effect_size=0.4, reversal_fraction=0.0, missing_fraction=0.0)
    q5 = _cell(nonlinear, effect_size=0.4, reversal_fraction=0.5, missing_fraction=0.0)
    v0 = _cell(multivariate, effect_size=0.4, orientation_spread=0.0, missing_fraction=0.0)
    v1 = _cell(multivariate, effect_size=0.4, orientation_spread=1.0, missing_fraction=0.0)

    matched = [
        m0["equal_matched_detection"],
        d0["equal_matched_detection"],
        d5["equal_matched_detection"],
        q0["equal_matched_detection"],
        q5["equal_matched_detection"],
        v0["equal_matched_detection"],
        v1["equal_matched_detection"],
    ]
    response_model = [
        m0["fixed_effect_logit_detection"],
        d0["common_slope_detection"],
        d5["common_slope_detection"],
        q0["quadratic_detection"],
        q5["quadratic_detection"],
        v0["common_vector_detection"],
        v1["common_vector_detection"],
    ]

    x = np.arange(len(scenarios))
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    ax.plot(x, matched, marker="o", label="Matched-null organization")
    ax.plot(x, response_model, marker="o", label="Shared-response model")
    ax.axhline(0.05, linewidth=1, linestyle="--")
    ax.set_xticks(x, scenarios)
    ax.set_ylabel("Detection fraction")
    ax.set_ylim(0, 1.08)
    ax.set_title("Performance follows the target estimand")
    ax.legend(frameon=False)
    return save_figure(fig, "disttrait_figure3_estimand_alignment")


def figure4() -> list[str]:
    meta = load_csv("meta")
    null = meta[meta["effect_size"] == 0.0]
    weak_reversed = _cell(
        meta,
        effect_size=0.4,
        reversal_fraction=0.5,
        missing_fraction=0.0,
    )

    labels = [
        "Null\nanalytic mean",
        "Null\nanalytic Q",
        "Null\nanalytic omnibus",
        "Null\ncalibrated mean",
        "Null\ncalibrated Q",
        "Null\ncalibrated omnibus",
        "Weak reversed\ncommon slope",
        "Weak reversed\nmatched-null",
        "Weak reversed\ncalibrated Q",
        "Weak reversed\ncalibrated omnibus",
    ]
    values = [
        float(null["meta_mean_detection"].max()),
        float(null["meta_heterogeneity_detection"].max()),
        float(null["meta_omnibus_detection"].max()),
        float(null["calibrated_meta_mean_detection"].max()),
        float(null["calibrated_meta_heterogeneity_detection"].max()),
        float(null["calibrated_meta_omnibus_detection"].max()),
        float(weak_reversed["common_slope_detection"]),
        float(weak_reversed["equal_matched_detection"]),
        float(weak_reversed["calibrated_meta_heterogeneity_detection"]),
        float(weak_reversed["calibrated_meta_omnibus_detection"]),
    ]

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(11.0, 4.8))
    ax.bar(x, values)
    ax.axhline(0.05, linewidth=1, linestyle="--", label="0.05")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_ylabel("Rejection / detection fraction")
    ax.set_ylim(0, 1.08)
    ax.set_title("Flexible species slopes still require calibration")
    ax.legend(frameon=False)
    return save_figure(fig, "disttrait_figure4_meta_calibration")


def figure5() -> list[str]:
    sf = load_json("sf")
    g = load_json("gammarus")

    labels = [
        "Street trees\nequal-taxon mean",
        "Gammarus\nlog rate",
        "Gammarus\nraw rate",
    ]
    rho = [
        float(sf["analysis"]["species_conditioned_mean_rho"]),
        float(g["primary_log_rate"]["observed_rho"]),
        float(g["raw_rate_sensitivity"]["observed_rho"]),
    ]
    p = [
        float(sf["analysis"]["species_conditioned_p_upper"]),
        float(g["primary_log_rate"]["p_upper"]),
        float(g["raw_rate_sensitivity"]["p_upper"]),
    ]

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    ax.axhline(0, linewidth=1)
    ax.scatter(x, rho, s=70)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Observed spatial rho")
    ax.set_title("External transport can be positive or non-supporting")
    for xi, yi, pi in zip(x, rho, p, strict=True):
        ax.annotate(f"p={pi:.3g}", (xi, yi), xytext=(0, 10), textcoords="offset points", ha="center")
    return save_figure(fig, "disttrait_figure5_empirical_transport")


def table1() -> str:
    rows = [
        {
            "method": "Naive pooled pairwise",
            "species_conditioned": "No",
            "trait_representation": "Pairwise dissimilarity",
            "direction_sensitive": "No",
            "assumption": "Pairs pooled across species",
            "calibration": "Nominal correlation test",
            "primary_estimand": "Pooled turnover + within-species structure",
            "benchmark_role": "Confounded baseline",
            "main_limitation": "Cannot isolate within-species organization",
        },
        {
            "method": "Equal-species matched-null",
            "species_conditioned": "Yes",
            "trait_representation": "Categorical, scalar or multivariate dissimilarity",
            "direction_sensitive": "No",
            "assumption": "Within-species exchangeability under null",
            "calibration": "Matched vertex permutation",
            "primary_estimand": "Mean strength of within-species spatial organization",
            "benchmark_role": "Assumption-light default",
            "main_limitation": "Conditional on observed rows; metric choice matters",
        },
        {
            "method": "Pair-weighted matched-null",
            "species_conditioned": "Yes",
            "trait_representation": "Pairwise dissimilarity",
            "direction_sensitive": "No",
            "assumption": "Within-species exchangeability under null",
            "calibration": "Matched vertex permutation",
            "primary_estimand": "Pair-weighted spatial organization",
            "benchmark_role": "Aggregation comparator",
            "main_limitation": "High-depth species receive greater inferential weight",
        },
        {
            "method": "Species-rho one-sample test",
            "species_conditioned": "Yes",
            "trait_representation": "Species rho",
            "direction_sensitive": "No",
            "assumption": "Parametric sampling of species summaries",
            "calibration": "One-sample t reference",
            "primary_estimand": "Mean species spatial rho",
            "benchmark_role": "Simple species-summary comparator",
            "main_limitation": "Reference distribution may be fragile",
        },
        {
            "method": "Fixed-intercept logistic common slope",
            "species_conditioned": "Yes",
            "trait_representation": "Binary response",
            "direction_sensitive": "Yes",
            "assumption": "Common logistic response form and slope",
            "calibration": "Profile likelihood ratio",
            "primary_estimand": "Shared signed binary response",
            "benchmark_role": "Correct-model efficiency comparator",
            "main_limitation": "Can lose target alignment under heterogeneous directions",
        },
        {
            "method": "Fixed-intercept continuous common slope",
            "species_conditioned": "Yes",
            "trait_representation": "Scalar continuous",
            "direction_sensitive": "Yes",
            "assumption": "Common linear slope",
            "calibration": "OLS/t reference",
            "primary_estimand": "Shared signed linear response",
            "benchmark_role": "Direction comparator",
            "main_limitation": "Opposing slopes cancel",
        },
        {
            "method": "Common quadratic curvature",
            "species_conditioned": "Yes",
            "trait_representation": "Scalar continuous",
            "direction_sensitive": "Yes",
            "assumption": "Shared quadratic response",
            "calibration": "OLS/t reference",
            "primary_estimand": "Shared curvature sign and magnitude",
            "benchmark_role": "Nonlinear comparator",
            "main_limitation": "Opposing curvature cancels",
        },
        {
            "method": "Common multivariate response vector",
            "species_conditioned": "Yes",
            "trait_representation": "Continuous multivariate",
            "direction_sensitive": "Yes",
            "assumption": "Shared trait-space response orientation",
            "calibration": "Multivariate Wald test",
            "primary_estimand": "Shared signed response vector",
            "benchmark_role": "Multivariate comparator",
            "main_limitation": "Different orientations cancel",
        },
        {
            "method": "Species-specific slope random effects",
            "species_conditioned": "Yes",
            "trait_representation": "Scalar continuous signed slopes",
            "direction_sensitive": "Yes",
            "assumption": "Meta-analytic slope distribution",
            "calibration": "Analytic large-sample reference",
            "primary_estimand": "Mean slope and slope heterogeneity",
            "benchmark_role": "Flexible-model diagnostic",
            "main_limitation": "Analytic calibration was anti-conservative in small-sample benchmark",
        },
        {
            "method": "Permutation-calibrated slope/meta omnibus",
            "species_conditioned": "Yes",
            "trait_representation": "Scalar continuous signed slopes",
            "direction_sensitive": "Yes",
            "assumption": "Within-species exchangeability under null",
            "calibration": "Matched within-species trait permutation",
            "primary_estimand": "Average signed response and/or directional heterogeneity",
            "benchmark_role": "Calibrated flexible comparator",
            "main_limitation": "Different estimand from distance-dissimilarity; computationally heavier",
        },
    ]
    frame = pd.DataFrame(rows)
    path = OUT / "disttrait_table1_method_estimand_map.csv"
    frame.to_csv(path, index=False)
    return path.name


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    outputs: list[str] = []
    outputs += figure1()
    outputs += figure2()
    outputs += figure3()
    outputs += figure4()
    outputs += figure5()
    outputs.append(table1())

    manifest = {
        "schema": "disttrait_methods_display_manifest_v1",
        "date_jst": "2026-09-19",
        "package_version": "0.12.0",
        "role": "reporting_only_from_frozen_receipts",
        "source_sha256": {
            key: sha256(path)
            for key, path in SOURCES.items()
        },
        "outputs": outputs,
        "claim_boundary": [
            "no simulation or biological inference is rerun by this generator",
            "finite benchmark rejection fractions are not universal type-I error guarantees",
            "model-based and distance-dissimilarity methods can target different estimands",
            "external empirical examples are transport validations, not causal analyses",
        ],
    }
    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PASS", "outputs": len(outputs)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
