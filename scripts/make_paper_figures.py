#!/usr/bin/env python3
"""Generate canonical JBI figures from frozen 34-species inputs and analysis outputs.

Figures 2--4 are generated directly from the checksum-locked 34-species
statistical freeze. Figure 5 and Supporting Figure S2 use outputs produced by
the canonical 34-species workflow. Figure 1 and Supporting Figure S1 use the
broader exact GBIF occurrence subset only as geographic context; those
occurrence points are not represented as the exact occurrence sample that
created the frozen climatic metrics.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from fcp_pipeline.constants import METRICS
from fcp_pipeline.models import fit_model, zscore
from fcp_pipeline.validation import validate_frozen_file

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except Exception as exc:  # pragma: no cover - dependency checked in figure CI
    raise RuntimeError("cartopy is required for geographic figures") from exc


CATEGORY_ORDER = ["within_population", "among_population"]
CATEGORY_LABEL = {
    "within_population": "Within-population coexistence",
    "among_population": "Geographic differentiation",
}
CATEGORY_SHORT = {
    "within_population": "Within",
    "among_population": "Among",
}
CATEGORY_COLOR = {
    "within_population": "#0072B2",
    "among_population": "#D55E00",
}
CATEGORY_MARKER = {
    "within_population": "o",
    "among_population": "^",
}
METRIC_LABELS = {
    "temperature_breadth": "Temperature breadth",
    "moisture_breadth": "Moisture breadth",
    "climatic_heterogeneity": "Climatic heterogeneity",
    "pca_dispersion": "PCA dispersion",
    "pca_hull_area": "PCA hull area",
}
SENSITIVITY_METHODS = [
    "Primary family-clustered",
    "CR2 / Satterthwaite",
    "Open Tree / Grafen",
    "Dated phylogeny S1-S3",
]
SENSITIVITY_MARKERS = ["o", "s", "D", "^"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="data/frozen/frozen_34species_five_metric_dataset.csv")
    p.add_argument(
        "--occurrences",
        default="docs/supporting/jbi_gbif_doi_bundle/jbi_gbif_exact_occurrence_subset.csv.gz",
    )
    p.add_argument("--analysis-dir", default="analysis_outputs/34species_paper")
    p.add_argument("--outdir", default="docs/figures")
    return p.parse_args()


def save_figure(fig: plt.Figure, outbase: Path) -> None:
    outbase.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outbase.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def read_required_csv(path: Path, required: set[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    d = pd.read_csv(path)
    missing = required - set(d.columns)
    if missing:
        raise RuntimeError(f"{path} missing columns: {sorted(missing)}")
    return d


def load_occurrences(path: Path, frozen: pd.DataFrame) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    occ = pd.read_csv(path, compression="gzip", low_memory=False)
    needed = {"canonical_name", "decimalLatitude", "decimalLongitude"}
    missing = needed - set(occ.columns)
    if missing:
        raise RuntimeError(f"Occurrence archive missing columns: {sorted(missing)}")
    occ = occ.loc[occ["canonical_name"].isin(frozen["canonical_name"])].copy()
    occ["decimalLatitude"] = pd.to_numeric(occ["decimalLatitude"], errors="coerce")
    occ["decimalLongitude"] = pd.to_numeric(occ["decimalLongitude"], errors="coerce")
    occ = occ.dropna(subset=["decimalLatitude", "decimalLongitude"])
    occ = occ.loc[
        occ["decimalLatitude"].between(-90, 90)
        & occ["decimalLongitude"].between(-180, 180)
    ].copy()
    scale_map = dict(zip(frozen["canonical_name"], frozen["spatial_scale"]))
    occ["spatial_scale"] = occ["canonical_name"].map(scale_map)
    if occ["canonical_name"].nunique() != 34:
        raise RuntimeError(
            f"Expected geographic context for 34 species; found {occ['canonical_name'].nunique()}"
        )
    return occ


def add_world_base(ax) -> None:
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor="0.94", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
    ax.coastlines(linewidth=0.35, color="0.35")
    ax.add_feature(cfeature.BORDERS, linewidth=0.18, edgecolor="0.7")


def figure1_global_context(frozen: pd.DataFrame, occ: pd.DataFrame, outdir: Path) -> None:
    fig = plt.figure(figsize=(11.2, 7.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.3, 1.0], hspace=0.16)
    ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    add_world_base(ax)
    pc = ccrs.PlateCarree()
    for scale in CATEGORY_ORDER:
        d = occ.loc[occ["spatial_scale"] == scale]
        ax.scatter(
            d["decimalLongitude"],
            d["decimalLatitude"],
            s=4,
            alpha=0.16,
            color=CATEGORY_COLOR[scale],
            marker=CATEGORY_MARKER[scale],
            linewidths=0,
            transform=pc,
            rasterized=True,
            label=CATEGORY_LABEL[scale],
        )
    ax.set_title("Geographic context of the 34 focal species", fontsize=12, pad=8)
    handles = [
        Line2D([0], [0], marker=CATEGORY_MARKER[s], color="none",
               markerfacecolor=CATEGORY_COLOR[s], markeredgecolor="none",
               markersize=7, label=CATEGORY_LABEL[s])
        for s in CATEGORY_ORDER
    ]
    ax.legend(handles=handles, loc="lower left", frameon=True, fontsize=8)

    ax2 = fig.add_subplot(gs[1])
    d = frozen.sort_values(["spatial_scale", "n_climate_cells", "canonical_name"]).reset_index(drop=True)
    y = np.arange(len(d))
    for scale in CATEGORY_ORDER:
        m = d["spatial_scale"].eq(scale)
        ax2.scatter(
            d.loc[m, "n_climate_cells"], y[m], s=22,
            color=CATEGORY_COLOR[scale], marker=CATEGORY_MARKER[scale], zorder=3,
        )
    ax2.set_yticks(y)
    ax2.set_yticklabels(d["canonical_name"].str.replace(" ", "\u00a0"), fontsize=5.7)
    ax2.set_xlabel("Occupied climate cells in the frozen analysis")
    ax2.grid(axis="x", linewidth=0.35, alpha=0.35)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.5, 0.01,
        "Map points show the broader exact GBIF occurrence subset for geographic context; "
        "the climatic analysis itself uses the checksum-locked 34-species frozen summaries.",
        ha="center", va="bottom", fontsize=7.5,
    )
    save_figure(fig, outdir / "figure1_geographic_context")


def figure_s1_species_maps(frozen: pd.DataFrame, occ: pd.DataFrame, outdir: Path) -> None:
    species = frozen.sort_values(["spatial_scale", "canonical_name"])["canonical_name"].tolist()
    scale_map = dict(zip(frozen["canonical_name"], frozen["spatial_scale"]))
    fig, axes = plt.subplots(
        6, 6, figsize=(17.5, 15.8),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    pc = ccrs.PlateCarree()
    for ax, sp in zip(axes.flat, species):
        d = occ.loc[occ["canonical_name"] == sp]
        scale = scale_map[sp]
        lon = d["decimalLongitude"].to_numpy()
        lat = d["decimalLatitude"].to_numpy()
        lon_min, lon_max = float(np.min(lon)), float(np.max(lon))
        lat_min, lat_max = float(np.min(lat)), float(np.max(lat))
        lon_span = lon_max - lon_min
        lat_span = lat_max - lat_min
        if lon_span > 300:
            extent = [-180, 180, max(-90, lat_min - 5), min(90, lat_max + 5)]
        else:
            pad_lon = max(2.5, lon_span * 0.10)
            pad_lat = max(2.5, lat_span * 0.10)
            extent = [
                max(-180, lon_min - pad_lon), min(180, lon_max + pad_lon),
                max(-90, lat_min - pad_lat), min(90, lat_max + pad_lat),
            ]
        ax.set_extent(extent, crs=pc)
        ax.add_feature(cfeature.LAND, facecolor="0.95")
        ax.coastlines(linewidth=0.35, color="0.35")
        ax.add_feature(cfeature.BORDERS, linewidth=0.18, edgecolor="0.7")
        ax.scatter(
            lon, lat, s=3.2, alpha=0.32, color=CATEGORY_COLOR[scale],
            marker=CATEGORY_MARKER[scale], linewidths=0, transform=pc, rasterized=True,
        )
        ax.set_title(f"{sp}\n{CATEGORY_SHORT[scale]} · n={len(d):,}", fontsize=6.8)
    for ax in axes.flat[len(species):]:
        ax.set_visible(False)
    handles = [
        Line2D([0], [0], marker=CATEGORY_MARKER[s], color="none",
               markerfacecolor=CATEGORY_COLOR[s], markeredgecolor="none",
               markersize=7, label=CATEGORY_LABEL[s])
        for s in CATEGORY_ORDER
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(
        "Supporting Figure S1. Geographic occurrence context for each focal species",
        fontsize=13, y=0.995,
    )
    fig.text(
        0.5, 0.008,
        "Panels use the broader exact GBIF occurrence subset as distribution context/QC and are not "
        "the morph-labelled data or the exact occurrence sample used to construct the frozen climate metrics.",
        ha="center", fontsize=7.5,
    )
    save_figure(fig, outdir / "figureS1_34_species_distribution_context")


def figure2_forest(frozen: pd.DataFrame, outdir: Path) -> None:
    rows = []
    for metric in METRICS:
        fit, _ = fit_model(frozen, metric, clustered=True)
        if fit is None:
            raise RuntimeError(f"Could not fit {metric}")
        beta = float(fit.params["metric_z"])
        se = float(fit.bse["metric_z"])
        rows.append({
            "metric": metric,
            "or": float(np.exp(beta)),
            "low": float(np.exp(beta - 1.96 * se)),
            "high": float(np.exp(beta + 1.96 * se)),
        })
    d = pd.DataFrame(rows)
    d["label"] = d["metric"].map(METRIC_LABELS)
    d = d.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    y = np.arange(len(d))
    xerr = np.vstack([d["or"] - d["low"], d["high"] - d["or"]])
    ax.errorbar(d["or"], y, xerr=xerr, fmt="o", capsize=3, linewidth=1.4, markersize=5)
    ax.axvline(1.0, linestyle="--", linewidth=1, color="0.45")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(d["label"])
    ax.set_xlabel("Odds ratio for geographic differentiation per 1 SD increase")
    ax.set_title("Five climatic-niche metrics show the same effect direction")
    ax.grid(axis="x", which="both", linewidth=0.35, alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    for yi, row in d.iterrows():
        ax.text(row["high"] * 1.03, yi, f"{row['or']:.2f}", va="center", fontsize=8)
    fig.text(
        0.5, 0.01,
        "OR < 1 indicates that broader occupied climatic breadth is associated with lower odds of "
        "geographically structured rather than within-population variation.",
        ha="center", fontsize=7.5,
    )
    save_figure(fig, outdir / "figure2_five_metric_forest")


def figure3_raw_species(frozen: pd.DataFrame, outdir: Path) -> None:
    rng = np.random.default_rng(20260821)
    fig, axes = plt.subplots(1, 5, figsize=(13.5, 4.8), sharey=True)
    for ax, metric in zip(axes, METRICS):
        z = zscore(frozen[metric])
        for xi, scale in enumerate(CATEGORY_ORDER):
            m = frozen["spatial_scale"].eq(scale)
            vals = z[m].to_numpy()
            jitter = rng.normal(0, 0.055, size=len(vals))
            ax.scatter(
                np.full(len(vals), xi) + jitter,
                vals,
                s=28,
                alpha=0.78,
                color=CATEGORY_COLOR[scale],
                marker=CATEGORY_MARKER[scale],
                edgecolors="white",
                linewidths=0.35,
            )
            med = float(np.median(vals))
            ax.plot([xi - 0.20, xi + 0.20], [med, med], color="black", linewidth=2)
        ax.axhline(0, linewidth=0.6, color="0.7")
        ax.set_xticks([0, 1], ["Within", "Among"], rotation=35, ha="right")
        ax.set_title(METRIC_LABELS[metric], fontsize=9.2)
        ax.grid(axis="y", linewidth=0.3, alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Species-level climatic metric (z score)")
    fig.suptitle("Raw 34-species climatic breadth by spatial organization", fontsize=12.5)
    fig.text(
        0.5, 0.01,
        "Each point is one species (20 within-population, 14 geographically structured); horizontal bars are medians.",
        ha="center", fontsize=7.5,
    )
    save_figure(fig, outdir / "figure3_raw_species_metrics")


def figure4_family_deletion(frozen: pd.DataFrame, outdir: Path) -> None:
    rows: list[dict] = []
    for metric in METRICS:
        full_fit, model_data = fit_model(frozen, metric, clustered=False)
        if full_fit is None:
            raise RuntimeError(f"Could not fit full model for {metric}")
        full_or = float(np.exp(full_fit.params["metric_z"]))
        families = sorted(model_data["family"].dropna().astype(str).unique())
        for family in families:
            subset = frozen.loc[frozen["family"].astype(str) != family].copy()
            fit, _ = fit_model(subset, metric, clustered=False)
            if fit is None:
                raise RuntimeError(f"Leave-one-family-out fit failed for {metric}: {family}")
            rows.append({
                "metric": metric,
                "omitted_family": family,
                "odds_ratio": float(np.exp(fit.params["metric_z"])),
                "full_odds_ratio": full_or,
            })

    d = pd.DataFrame(rows)
    metric_order = list(reversed(METRICS))
    fig, ax = plt.subplots(figsize=(8.0, 5.1))
    y_base = np.arange(len(metric_order))
    for yi, metric in enumerate(metric_order):
        x = d.loc[d["metric"] == metric, "odds_ratio"].to_numpy()
        full_or = float(d.loc[d["metric"] == metric, "full_odds_ratio"].iloc[0])
        jitter = np.linspace(-0.13, 0.13, len(x))
        ax.hlines(yi, float(np.min(x)), float(np.max(x)), linewidth=2.2, color="0.72", zorder=1)
        ax.scatter(x, np.full(len(x), yi) + jitter, s=17, alpha=0.58, color="0.45", zorder=2)
        ax.scatter([full_or], [yi], s=54, marker="D", color="black", zorder=3)

    ax.axvline(1.0, linestyle="--", linewidth=1, color="0.45")
    ax.set_xscale("log")
    ax.set_yticks(y_base)
    ax.set_yticklabels([METRIC_LABELS[m] for m in metric_order])
    ax.set_xlabel("Odds ratio after omitting one plant family")
    ax.set_title("The shared effect direction is not concentrated in one family")
    ax.grid(axis="x", which="both", linewidth=0.35, alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="0.45",
               markeredgecolor="none", alpha=0.7, label="One family omitted"),
        Line2D([0], [0], marker="D", linestyle="none", markerfacecolor="black",
               markeredgecolor="black", label="Full 34-species estimate"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8, frameon=True)
    fig.text(
        0.5, 0.01,
        "Each grey point is one unclustered leave-one-family-out refit. All 125 refits remain below OR = 1; "
        "family deletion tests concentration, not phylogenetic independence.",
        ha="center", fontsize=7.5,
    )
    save_figure(fig, outdir / "figure4_leave_one_family_out")


def build_method_sensitivity_table(analysis_dir: Path) -> pd.DataFrame:
    primary = read_required_csv(
        analysis_dir / "models" / "environmental_niche_five_metric_models.csv",
        {"metric", "odds_ratio", "odds_ratio_ci_low", "odds_ratio_ci_high"},
    )
    opentree = read_required_csv(
        analysis_dir / "phylogeny" / "environmental_niche_opentree_summary.csv",
        {"metric", "median_odds_ratio", "median_ci_low", "median_ci_high"},
    )
    dated = read_required_csv(
        analysis_dir / "phylogeny" / "environmental_niche_dated_phyloglm.csv",
        {"metric", "scenario", "odds_ratio", "ci_low", "ci_high", "status"},
    )
    cr2 = read_required_csv(
        analysis_dir / "finite_sample" / "cr2" / "jbi_34species_cr2_satterthwaite_summary.csv",
        {"metric", "odds_ratio", "cr2_ci_low", "cr2_ci_high"},
    )

    rows: list[dict] = []
    for metric in METRICS:
        p = primary.loc[primary["metric"].eq(metric)].iloc[0]
        rows.append({
            "metric": metric,
            "method": "Primary family-clustered",
            "or": float(p["odds_ratio"]),
            "low": float(p["odds_ratio_ci_low"]),
            "high": float(p["odds_ratio_ci_high"]),
        })

        c = cr2.loc[cr2["metric"].eq(metric)].iloc[0]
        rows.append({
            "metric": metric,
            "method": "CR2 / Satterthwaite",
            "or": float(c["odds_ratio"]),
            "low": float(c["cr2_ci_low"]),
            "high": float(c["cr2_ci_high"]),
        })

        o = opentree.loc[opentree["metric"].eq(metric)].iloc[0]
        rows.append({
            "metric": metric,
            "method": "Open Tree / Grafen",
            "or": float(o["median_odds_ratio"]),
            "low": float(o["median_ci_low"]),
            "high": float(o["median_ci_high"]),
        })

        q = dated.loc[dated["metric"].eq(metric) & dated["status"].eq("complete")].copy()
        if q.empty:
            raise RuntimeError(f"No complete dated-phylogeny results for {metric}")
        rows.append({
            "metric": metric,
            "method": "Dated phylogeny S1-S3",
            "or": float(q["odds_ratio"].median()),
            "low": float(q["ci_low"].min()),
            "high": float(q["ci_high"].max()),
        })
    return pd.DataFrame(rows)


def figure5_method_sensitivity(analysis_dir: Path, outdir: Path) -> None:
    d = build_method_sensitivity_table(analysis_dir)
    metric_order = list(reversed(METRICS))
    fig, ax = plt.subplots(figsize=(9.2, 6.1))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    offsets = np.linspace(-0.27, 0.27, len(SENSITIVITY_METHODS))
    y_base = np.arange(len(metric_order))

    for mi, method in enumerate(SENSITIVITY_METHODS):
        xs, ys, lows, highs = [], [], [], []
        for yi, metric in enumerate(metric_order):
            row = d.loc[d["metric"].eq(metric) & d["method"].eq(method)].iloc[0]
            xs.append(float(row["or"]))
            lows.append(float(row["low"]))
            highs.append(float(row["high"]))
            ys.append(float(yi + offsets[mi]))
        xs = np.asarray(xs)
        lows = np.asarray(lows)
        highs = np.asarray(highs)
        xerr = np.vstack([xs - lows, highs - xs])
        ax.errorbar(
            xs, ys, xerr=xerr,
            fmt=SENSITIVITY_MARKERS[mi],
            color=colors[mi % len(colors)],
            capsize=2.5,
            linewidth=1.1,
            markersize=4.8,
            label=method,
        )

    ax.axvline(1.0, linestyle="--", linewidth=1, color="0.45")
    ax.set_xscale("log")
    ax.set_yticks(y_base)
    ax.set_yticklabels([METRIC_LABELS[m] for m in metric_order])
    ax.set_xlabel("Odds ratio for geographic differentiation per 1 SD increase")
    ax.set_title("Direction persists while uncertainty expands under stricter treatments")
    ax.grid(axis="x", which="both", linewidth=0.35, alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", fontsize=7.6, frameon=True)
    fig.text(
        0.5, 0.01,
        "Open Tree points are medians across 100 polytomy resolutions; dated-phylogeny points are medians across "
        "S1-S3 with the interval envelope shown. These are sensitivity analyses, not four independent tests.",
        ha="center", fontsize=7.2,
    )
    save_figure(fig, outdir / "figure5_inference_method_sensitivity")


def figure_s2_power_precision(analysis_dir: Path, outdir: Path) -> None:
    sim = read_required_csv(
        analysis_dir / "finite_sample" / "power_precision" / "jbi_34species_power_precision_simulation.csv",
        {"metric", "scenario", "odds_ratio_true", "prob_estimate_negative", "prob_p_lt_0_05"},
    )
    grid = sim.loc[sim["scenario"].astype(str).str.startswith("OR_")].copy()
    grid["odds_ratio_true"] = pd.to_numeric(grid["odds_ratio_true"], errors="coerce")
    grid = grid.dropna(subset=["odds_ratio_true"]).sort_values("odds_ratio_true")

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), sharex=True, sharey=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for mi, metric in enumerate(METRICS):
        d = grid.loc[grid["metric"].eq(metric)].sort_values("odds_ratio_true")
        axes[0].plot(
            d["odds_ratio_true"], d["prob_estimate_negative"],
            marker="o", markersize=3.5, linewidth=1.1,
            color=colors[mi % len(colors)], label=METRIC_LABELS[metric],
        )
        axes[1].plot(
            d["odds_ratio_true"], d["prob_p_lt_0_05"],
            marker="o", markersize=3.5, linewidth=1.1,
            color=colors[mi % len(colors)],
        )
    for ax in axes:
        ax.axhline(0.05, linestyle=":", linewidth=0.8, color="0.55")
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Simulated true odds ratio")
        ax.grid(linewidth=0.3, alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Probability across 3,000 simulations")
    axes[0].set_title("Estimated coefficient is negative")
    axes[1].set_title("Family-clustered Wald p < 0.05")
    axes[0].legend(fontsize=7.2, frameon=True, loc="lower left")
    fig.suptitle("Supporting Figure S2. Finite-sample design diagnostic", fontsize=12)
    fig.text(
        0.5, 0.01,
        "Simulations retain the observed 34-species predictor, effort and family structure. They describe design "
        "precision under specified effects and are not evidence for the ecological hypothesis.",
        ha="center", fontsize=7.2,
    )
    save_figure(fig, outdir / "figureS2_power_precision_design")


def main() -> None:
    args = parse_args()
    dataset = Path(args.dataset)
    occ_path = Path(args.occurrences)
    analysis_dir = Path(args.analysis_dir)
    outdir = Path(args.outdir)
    frozen = validate_frozen_file(dataset)
    if len(frozen) != 34:
        raise RuntimeError(f"Expected 34 frozen species; found {len(frozen)}")
    occ = load_occurrences(occ_path, frozen)
    figure1_global_context(frozen, occ, outdir)
    figure_s1_species_maps(frozen, occ, outdir)
    figure2_forest(frozen, outdir)
    figure3_raw_species(frozen, outdir)
    figure4_family_deletion(frozen, outdir)
    figure5_method_sensitivity(analysis_dir, outdir)
    figure_s2_power_precision(analysis_dir, outdir)
    print(f"Wrote canonical paper figures to {outdir}")


if __name__ == "__main__":
    main()
