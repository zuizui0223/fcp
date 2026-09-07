#!/usr/bin/env python3
"""Render discovery-only publication figures without reopening any inference."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE_COMMIT = "29584f3ad7ae0cd99a1d8f43459252af38f615da"
PREFIX = "global_rgfca_within_species_spatial_omnibus"
SOURCES = {
    "result": f"docs/supporting/{PREFIX}_result_v1.json",
    "species": f"data/derived/{PREFIX}_species_v1.csv",
    "null": f"data/derived/{PREFIX}_null_v1.csv",
    "measured": "data/derived/global_monte_carlo_measured_photos_v1.csv",
}
COLOURS = ["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]
ANCHORS = ["#F4F1E8", "#DBAD34", "#B84A78", "#655AB0"]
INK, BLUE = "#30343A", "#477A9F"


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_sources(root: Path = ROOT) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    # Exact committed bytes, not implicit CRLF repair of a working-tree file.
    raw = {key: subprocess.check_output(
        ["git", "show", f"{SOURCE_COMMIT}:{path}"], cwd=root
    ) for key, path in SOURCES.items()}
    result = json.loads(raw["result"])
    verification = result["verification"]
    for key, field in [("species", "species_results_sha256"), ("null", "global_null_sha256"),
                       ("measured", "measured_table_sha256")]:
        if digest(raw[key]) != verification[field]:
            raise ValueError(f"Frozen {key} input digest differs")
    frames = {key: pd.read_csv(io.BytesIO(raw[key])) for key in ("species", "null", "measured")}
    lineage = {SOURCES[key]: digest(value) for key, value in raw.items()}
    return result, frames["species"], frames["null"], frames["measured"], lineage


def validate_inputs(result: dict, species: pd.DataFrame, null: pd.DataFrame,
                    measured: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    if result["status"] != "complete_verified_rgfca_within_species_spatial_omnibus":
        raise ValueError("Incomplete discovery result")
    if result["six_species_used"] or result["thirty_four_species_used"]:
        raise ValueError("Legacy inference cannot enter these figures")
    if len(species) != 369 or species.species.duplicated().any():
        raise ValueError("Species census differs")
    if len(null) != 999 or sorted(null.permutation_index.tolist()) != list(range(999)):
        raise ValueError("Null-index census differs")
    if len(measured) != 50000 or measured.photo_id.duplicated().any():
        raise ValueError("Measured photo census differs")
    if measured.species.nunique() != 500 or not measured.groupby("species").size().eq(100).all():
        raise ValueError("Measured species census differs")
    flags = measured.global_classifiable.astype(str).str.casefold()
    if not flags.isin(["true", "false"]).all():
        raise ValueError("Unknown classifiable flag")
    classified = measured.loc[flags.eq("true")].copy()
    counts = classified.groupby("species").size()
    eligible = counts[counts >= 40]
    if len(classified) != 25377 or set(eligible.index) != set(species.species):
        raise ValueError("Classifiable/eligible census differs")
    pool = classified.loc[classified.species.isin(eligible.index)].copy()
    if len(pool) != 21424:
        raise ValueError("Eligible photo count differs")
    reported = species.set_index("species").sort_index()
    if not reported.photos.eq(eligible.sort_index()).all():
        raise ValueError("Species photo counts disagree")
    if not reported.pairs.eq(reported.photos * (reported.photos - 1) // 2).all():
        raise ValueError("Species pair counts disagree")
    vectors = pool[COLOURS].to_numpy(float)
    if not np.isfinite(vectors).all() or (vectors < 0).any() or not np.allclose(vectors.sum(axis=1), 1, atol=1e-10, rtol=0):
        raise ValueError("Invalid colour vectors")
    coords = pool[["latitude", "longitude"]].to_numpy(float)
    if not np.isfinite(coords).all() or (np.abs(coords[:, 0]) > 90).any() or (np.abs(coords[:, 1]) > 180).any():
        raise ValueError("Invalid coordinates")
    values = species.observed_rho.to_numpy(float)
    null_values = null.mean_rho.to_numpy(float)
    if not np.isfinite(values).all() or not np.isfinite(null_values).all() or (np.abs(values) > 1).any():
        raise ValueError("Invalid effect/null values")
    observed = float(np.mean(values))
    # Literal frozen upper-tail rule; no new randomization or tie adjustment.
    p = float((1 + np.count_nonzero(null_values >= observed)) / (len(null_values) + 1))
    checks = {
        "observed_mean_rho": observed, "p_upper": p,
        "null_mean": float(np.mean(null_values)), "null_sd": float(np.std(null_values, ddof=1)),
        "null_q025": float(np.quantile(null_values, .025)),
        "null_q975": float(np.quantile(null_values, .975)),
    }
    for key, value in checks.items():
        if not np.isclose(value, result["primary"][key], atol=1e-12, rtol=0):
            raise ValueError(f"Recomputed {key} differs")
    checks.update(species=len(species), eligible_photos=len(pool), permutations=len(null))
    return pool, checks


def display_rgb(pool: pd.DataFrame) -> np.ndarray:
    return pool[COLOURS].to_numpy(float) @ np.array([to_rgb(c) for c in ANCHORS])


def configure() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "savefig.facecolor": "white"})


def save_pair(fig: plt.Figure, output: Path, stem: str) -> list[dict]:
    files = []
    for suffix in ("png", "pdf"):
        path = output / f"{stem}.{suffix}"
        metadata = {"Creator": "FCP RGFCA", "CreationDate": None, "ModDate": None} if suffix == "pdf" else {"Software": "FCP RGFCA"}
        fig.savefig(path, dpi=240, metadata=metadata)
        files.append({"name": path.name, "sha256": digest(path.read_bytes()), "bytes": path.stat().st_size})
    plt.close(fig)
    return files


def draw_atlas(pool: pd.DataFrame, output: Path) -> tuple[list[dict], dict]:
    import geopandas as gpd
    if gpd.__version__ != "0.14.4":
        raise ValueError("Display-only basemap requires GeoPandas 0.14.4")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        path = Path(gpd.datasets.get_path("naturalearth_lowres"))
    basemap = {"source": "Natural Earth low resolution, bundled with GeoPandas 0.14.4",
               "use": "display only; no inference", "sha256": {
                   p.name: digest(p.read_bytes()) for p in sorted(path.parent.glob("naturalearth_lowres.*"))}}
    world = gpd.read_file(path)
    fig, ax = plt.subplots(figsize=(11, 5.6))
    fig.subplots_adjust(left=.07, right=.98, top=.82, bottom=.21)
    world.plot(ax=ax, facecolor="#EEEEEE", edgecolor="#B8B8B8", linewidth=.3, zorder=0)
    order = pool.photo_id.astype(str).map(lambda p: digest(p.encode())).sort_values().index
    points = pool.loc[order]
    ax.scatter(points.longitude, points.latitude, c=display_rgb(points), s=5,
               linewidths=.12, edgecolors="#545454", alpha=1, rasterized=True, zorder=2)
    ax.set(xlim=(-180, 180), ylim=(-90, 90), xlabel="Longitude (degrees)", ylabel="Latitude (degrees)")
    ax.set_aspect("auto")
    ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-90, 91, 30))
    ax.grid(color="#DDDDDD", linewidth=.4, zorder=-1)
    fig.text(.07, .94, "Discovery photo-derived colour map", fontsize=15)
    fig.text(.07, .885, "21,424 photographs | 369 species | one point per eligible photo | discovery data only", fontsize=10)
    handles = [Line2D([], [], marker="o", linestyle="none", markerfacecolor=c,
        markeredgecolor="#545454", markersize=6, label=label)
        for c, label in zip(ANCHORS, ["White", "Yellow / orange", "Red / pink", "Blue / purple"])]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.52, .085), ncol=4, frameon=False)
    fig.text(.07, .05, "Point colours mix display anchors, not calibrated reflectance. Uneven coverage and overlap are not biological abundance.", fontsize=9)
    fig.text(.07, .018, "Species-free display; all inference remains species-conditioned. Natural Earth land outlines: reference only.", fontsize=9)
    return save_pair(fig, output, "rgfca_figure1_discovery_atlas"), basemap


def draw_omnibus(species: pd.DataFrame, null: pd.DataFrame, checks: dict, output: Path) -> list[dict]:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    fig.subplots_adjust(left=.075, right=.975, bottom=.26, top=.70, wspace=.29)
    mean = checks["observed_mean_rho"]
    axes[0].hist(species.observed_rho, bins=24, color=BLUE, edgecolor="white", linewidth=.7)
    axes[0].axvline(0, color=INK, ls=":", lw=1)
    axes[0].axvline(mean, color=INK, lw=1.5, label=f"Species-equal mean = {mean:.4f}")
    axes[0].set(title="a   Observed species effects", xlabel="Within-species Spearman rho", ylabel="Number of species")
    axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    axes[1].axvspan(checks["null_q025"], checks["null_q975"], color="#EEEEEE", zorder=0)
    axes[1].hist(null.mean_rho, bins=30, color=BLUE, edgecolor="white", linewidth=.7)
    axes[1].axvline(0, color=INK, ls=":", lw=1)
    axes[1].axvline(mean, color=INK, lw=1.5)
    axes[1].set(title="b   Species-conditioned null", xlabel="Species-equal mean rho", ylabel="Number of permutations")
    axes[1].text(.88, .95, f"Observed = {mean:.4f}\nUpper-tail p = {checks['p_upper']:.3f}",
                 ha="right", va="top", transform=axes[1].transAxes, fontsize=9)
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.tick_params(axis="both", labelsize=9)
    fig.text(.075, .94, "Within-species distance-colour association", fontsize=15)
    fig.text(.075, .875, "369 equally weighted species | 21,424 photos | 999 within-species permutations", fontsize=10)
    fig.text(.075, .82, "Small positive exploratory association; independent replication and flower-specific controls remain pending.", fontsize=10)
    fig.text(.075, .115, "All species are shown, not only significant ones. Pairwise distances are not independent replicates.", fontsize=9)
    fig.text(.075, .065, "Gray band: central 95% of the null distribution, NOT a confidence interval on the observed effect.", fontsize=9)
    fig.text(.075, .015, "A distance association does not establish a shared geographic boundary, environmental driver or pollinator mechanism.", fontsize=9)
    return save_pair(fig, output, "rgfca_figure2_discovery_omnibus")


def build(output: Path, manifest_path: Path) -> dict:
    result, species, null, measured, lineage = load_sources()
    pool, checks = validate_inputs(result, species, null, measured)
    output.mkdir(parents=True, exist_ok=True)
    configure()
    atlas_files, basemap = draw_atlas(pool, output)
    files = atlas_files + draw_omnibus(species, null, checks, output)
    manifest = {"protocol": "rgfca-publication-figures-v1",
        "status": "discovery_figures_reconstructed_not_submission_ready",
        "source_commit": SOURCE_COMMIT, "input_sha256": lineage,
        "source_bytes": "exact git blobs; no newline normalization or input edits",
        "reconstructed": checks, "outputs": files, "basemap": basemap,
        "reserve_outcomes_read": False, "legacy_results_read": False,
        "new_inferential_tests_run": False, "photo_bar_present": False,
        "independent_replication_claim_allowed": False, "shared_boundary_claim_allowed": False,
        "display_palette": dict(zip(COLOURS, ANCHORS)),
        "renderer": {"matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__}}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs/figures")
    parser.add_argument("--manifest", type=Path, default=ROOT / "docs/supporting/rgfca_publication_figure_manifest_v1.json")
    args = parser.parse_args()
    print(json.dumps(build(args.output_dir, args.manifest), indent=2))
