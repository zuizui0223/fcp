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

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.7), gridspec_kw={"width_ratios": [1.25, 1.0]})
    fig.suptitle("A continuous species-level flower-colour polymorphism phenotype", y=1.02, fontweight="bold")

    ax = axes[0]
    bins = np.linspace(0, 0.72, 25)
    for cohort, colour, label in [
        ("discovery", SECONDARY, "Discovery (n=369)"),
        ("reserve", PRIMARY, "Reserve (n=363)"),
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
        "High-depth validation cohorts\nnot a prevalence sample",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color=NEUTRAL,
    )
    panel_label(ax, "A")

    ax = axes[1]
    ax.axis("off")
    boxes = [
        (0.5, 0.86, "42,111-species global frame\nSampling / opportunity universe"),
        (0.5, 0.62, "Original high-depth source\n500 discovery + 500 reserve\n100 photos per species"),
        (0.5, 0.38, "D inference after ≥40 classifiable\n369 discovery + 363 reserve"),
        (0.5, 0.14, "Prospective H2 third cohort\n499 species × 100 rows\n377 measurement-evaluable"),
    ]
    for idx, (x, y, text) in enumerate(boxes):
        ax.text(
            x,
            y,
            text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "white", "edgecolor": "#999999"},
        )
        if idx < len(boxes) - 1:
            ax.annotate(
                "",
                xy=(0.5, y - 0.105),
                xytext=(0.5, y - 0.175),
                xycoords=ax.transAxes,
                arrowprops={"arrowstyle": "->", "color": "#777777", "lw": 1.2},
            )
    panel_label(ax, "B")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure1_measurement_frame")
    meta = {
        "D_source": str(d_path.relative_to(root)),
        "discovery_n": 369,
        "reserve_n": 363,
        "global_frame_species": 42111,
        "third_cohort_species": 499,
        "third_cohort_measurement_evaluable_species": 377,
    }
    return files, meta


def figure2(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    repeated = load_json(root / "results" / "polymorphism_h1_observer_disjoint_reliability_20260913" / "result.json")
    strict = load_json(root / "results" / "polymorphism_h1_observer_disjoint_D_reliability_20260913" / "result.json")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    fig.suptitle("H1: observer-disjoint reproducibility of the polymorphism score", y=1.02, fontweight="bold")

    ax = axes[0]
    cohorts = ["discovery", "reserve"]
    ys = np.arange(2)
    for y, cohort, colour in zip(ys, cohorts, [SECONDARY, PRIMARY], strict=True):
        s = repeated[cohort]["primary20"]
        ax.hlines(y, s["rho_q05"], s["rho_q95"], color=colour, linewidth=3)
        ax.scatter(s["rho_median"], y, s=80, color=colour, edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(s["rho_q95"] + 0.005, y, f"median {s['rho_median']:.3f}", va="center", fontsize=8.5)
    ax.axvline(repeated["primary_gate"]["minimum_median_split_rho"], color=NEUTRAL, linestyle="--", linewidth=1.2)
    ax.set_yticks(ys, ["Discovery", "Reserve"])
    ax.invert_yaxis()
    ax.set_xlim(0.63, 0.86)
    ax.set_xlabel("Observer-disjoint split Spearman rho")
    ax.set_title("200 frozen partitions: q05–median–q95")
    ax.text(
        0.02,
        0.05,
        "Primary median floor = 2/3\nReserve q05 = 0.765",
        transform=ax.transAxes,
        fontsize=8.5,
        color=NEUTRAL,
    )
    panel_label(ax, "A")

    ax = axes[1]
    rows = []
    for cohort in cohorts:
        s = strict[cohort]
        rows.append((cohort, s["spearman_D_A_D_B"], *s["spearman_bootstrap_95_percentile_ci"], s["lin_ccc"]))
    for y, (cohort, rho, low, high, ccc), colour in zip(ys, rows, [SECONDARY, PRIMARY], strict=True):
        ax.hlines(y, low, high, color=colour, linewidth=3)
        ax.scatter(rho, y, s=80, color=colour, edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(high + 0.004, y, f"rho={rho:.3f}; CCC={ccc:.3f}", va="center", fontsize=8.2)
    ax.axvline(0.80, color=NEUTRAL, linestyle="--", linewidth=1.3, label="Strict floor = 0.80")
    ax.set_yticks(ys, ["Discovery", "Reserve"])
    ax.invert_yaxis()
    ax.set_xlim(0.72, 0.88)
    ax.set_xlabel("Deterministic split Spearman rho")
    ax.set_title("Later prespecified stress test: bootstrap 95% CI")
    ax.legend(frameon=False, loc="lower right")
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
        0.02,
        0.96,
        "Named axis isolated after\noriginal broad H2 was opened",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color=NEUTRAL,
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
                    "label": f"{tier_label}\n{cohort.capitalize()}",
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
    ax.legend(frameon=False, loc="lower right")
    panel_label(ax, "B")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure3_h2_target_localization")
    meta = {
        "status": result["status"],
        "verdict": result["decision"]["verdict"],
        "structured_null_replicates": result["structured_null_replicates"],
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
    fig.text(
        0.5,
        -0.01,
        "Species-disjoint prospective confirmation within the same iNaturalist opportunity universe; not an independent-source replication.",
        ha="center",
        fontsize=8.5,
        color=NEUTRAL,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
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
    }
    return files, meta


def figure5(root: Path, output_dir: Path) -> tuple[dict[str, str], dict]:
    h3a = load_json(root / "results" / "polymorphism_h3a_phylogenetic_signal_20260912" / "frozen_result_manifest.json")
    source_dir = root / "results" / "polymorphism_publication_figure_source_20260918"
    h3b = pd.read_csv(source_dir / "h3b_span_summary.csv")
    h3b_freeze = (root / "docs" / "POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md").read_text(encoding="utf-8")
    h3b_verdict = "H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED"
    if h3b_verdict not in h3b_freeze:
        raise ValueError("H3b frozen verdict not found in canonical result freeze")

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6))
    fig.suptitle("H3: broad explanatory tests do not replicate in reserve", y=1.02, fontweight="bold")

    ax = axes[0]
    scenarios = ["S1", "S2", "S3"]
    k = [float(h3a["reserve_primary"][s]["K"]) for s in scenarios]
    p = [float(h3a["reserve_primary"][s]["p_K"]) for s in scenarios]
    x = np.arange(3)
    ax.scatter(x, k, s=100, color=SECONDARY, edgecolor="white", linewidth=0.8)
    ax.plot(x, k, color=SECONDARY, linewidth=1.2, alpha=0.7)
    for xi, kval, pval in zip(x, k, p, strict=True):
        ax.text(xi, kval + 0.004, f"p={pval:.4f}", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x, scenarios)
    ax.set_ylabel("Reserve Blomberg K")
    ax.set_xlabel("Frozen tree-placement scenario")
    ax.set_ylim(0, max(k) + 0.035)
    ax.set_title("Broad tree-wide phylogenetic signal")
    ax.text(
        0.03,
        0.93,
        "0/3 raw-D scenarios p < 0.05",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color=NEUTRAL,
    )
    panel_label(ax, "A")

    ax = axes[1]
    order = ["discovery", "reserve"]
    x = np.arange(2)
    vals = [float(h3b.loc[h3b["cohort"].eq(c), "rho_D_span"].iloc[0]) for c in order]
    ps = [float(h3b.loc[h3b["cohort"].eq(c), "p_D_span"].iloc[0]) for c in order]
    ax.axhline(0, color="#777777", linewidth=1.0)
    ax.bar(x, vals, color=[SECONDARY, PRIMARY], width=0.55)
    for xi, val, pval in zip(x, vals, ps, strict=True):
        offset = 0.012 if val >= 0 else -0.012
        ax.text(xi, val + offset, f"rho={val:.3f}\np={pval:.4f}", ha="center", va="bottom" if val >= 0 else "top", fontsize=8.5)
    ax.set_xticks(x, ["Discovery\ncalibration", "Reserve\nreplication"])
    ax.set_ylabel("Spearman rho(D, sampled span)")
    ax.set_ylim(-0.07, 0.23)
    ax.set_title("Sampled photographic span")
    ax.text(
        0.03,
        0.06,
        "Predictor is sampled span,\nnot true biological range size",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=NEUTRAL,
    )
    panel_label(ax, "B")

    fig.tight_layout()
    files = save_pair(fig, output_dir, "polymorphism_figure5_explanatory_boundaries")
    meta = {
        "h3a": {
            "verdict": h3a["decision"]["verdict"],
            "reserve_scenarios": {
                s: {"K": float(h3a["reserve_primary"][s]["K"]), "p": float(h3a["reserve_primary"][s]["p_K"])}
                for s in scenarios
            },
        },
        "h3b": {
            "verdict": h3b_verdict,
            "discovery_rho": vals[0],
            "discovery_p": ps[0],
            "reserve_rho": vals[1],
            "reserve_p": ps[1],
        },
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
        },
        "hard_nonclaims": [
            "high-depth cohorts do not estimate global flower-colour polymorphism prevalence",
            "third-cohort H2 is species-disjoint within the same iNaturalist opportunity universe, not independent-source replication",
            "the white/nonwhite axis does not identify pigment chemistry, transition direction, or adaptive mechanism",
            "H3a non-support does not imply absence of all phylogenetic structure",
            "sampled photographic span is not true biological range size",
        ],
    }
    manifest_path = output_dir / "polymorphism_figure_manifest_20260918.json"
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = generate_all(args.root, args.output_dir)
    print(f"Wrote frozen polymorphism publication figures and manifest: {manifest}")


if __name__ == "__main__":
    main()
