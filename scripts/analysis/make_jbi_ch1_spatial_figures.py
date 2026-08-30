#!/usr/bin/env python3
"""Generate the canonical Chapter 1 spatial result figures.

The script reads only committed, frozen evaluation products.  It does not refit models,
change spatial support, inspect environmental layers or create new inferential tests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

SPECIES = [
    "Antirrhinum majus",
    "Dactylorhiza sambucina",
    "Gentiana lutea",
    "Ipomoea purpurea",
    "Lysimachia arvensis",
    "Raphanus sativus",
]
SPECIES_ABBREVIATED = {
    "Antirrhinum majus": r"$\it{A.\ majus}$",
    "Dactylorhiza sambucina": r"$\it{D.\ sambucina}$",
    "Gentiana lutea": r"$\it{G.\ lutea}$",
    "Ipomoea purpurea": r"$\it{I.\ purpurea}$",
    "Lysimachia arvensis": r"$\it{L.\ arvensis}$",
    "Raphanus sativus": r"$\it{R.\ sativus}$",
}
SPECIES_SLUG = {species: species.lower().replace(" ", "_") for species in SPECIES}
CAP_COLOURS = {500: "#0072B2", 1000: "#D55E00", 2000: "#009E73"}
PRIMARY_COLOUR = "#D55E00"
NULL_COLOUR = "#0072B2"
SUPPORT_COLOUR = "#009E73"
NEUTRAL = "#666666"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def require_columns(frame: pd.DataFrame, columns: list[str], *, label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def configure_matplotlib() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 14,
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


def save_pair(fig: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
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
            "Author": "fcp Chapter 1 frozen figure pipeline",
            "Creator": "matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)
    return [png, pdf]


def validate_inputs(stage_a: dict, stage_b: dict, stage_a_null: pd.DataFrame, surface: pd.DataFrame, split: pd.DataFrame) -> None:
    if stage_a.get("status") != "stage_a_evaluation_complete":
        raise ValueError("Stage A result is not complete")
    if stage_a.get("n_evaluation_records") != 720:
        raise ValueError("Stage A must contain 720 evaluation records")
    if stage_a.get("primary_k") != 5 or stage_a.get("sensitivity_k") != [3, 8]:
        raise ValueError("unexpected Stage A graph contract")
    if stage_a.get("primary_rejects_random_labelling_at_0_05") is not True:
        raise ValueError("Stage A gate was not passed")
    if stage_b.get("status") != "stage_b_evaluation_complete":
        raise ValueError("Stage B result is not complete")
    if stage_b.get("n_evaluation_records") != 720:
        raise ValueError("Stage B must contain 720 evaluation records")
    if stage_b.get("geometry_selection_used_colour_values") is not False:
        raise ValueError("Stage B geometry was not label-blind")
    if stage_b.get("environment_used") is not False or stage_b.get("geographic_reference_library_used") is not False:
        raise ValueError("environment/geographic-reference layers entered the frozen analysis")
    if stage_b.get("primary_rejects_shared_concentration_null_at_0_05") is not False:
        raise ValueError("unexpected Stage B primary decision")

    require_columns(
        stage_a_null,
        ["permutation", "global_equal_species_mean_q", *SPECIES],
        label="Stage A null CSV",
    )
    if len(stage_a_null) != 9999:
        raise ValueError("Stage A null CSV must contain 9,999 permutations")
    require_columns(
        surface,
        [
            "cell_id",
            "latitude",
            "longitude",
            "opportunity_A",
            "evaluable_A_ge_minimum",
            "shared_transition_intensity",
        ],
        label="Stage B surface CSV",
    )
    selected = stage_b["selected_primary_configuration"]
    if len(surface) != int(selected["grid"]["n_cells"]):
        raise ValueError("Stage B surface row count does not match the selected grid")
    require_columns(split, ["species", "photo_id", "latitude", "longitude", "split"], label="frozen split")
    evaluation = split.loc[split["split"].astype(str).eq("evaluation")]
    if len(evaluation) != 720 or evaluation["photo_id"].astype(str).nunique() != 720:
        raise ValueError("frozen split must contain 720 unique evaluation rows")


def make_stage_a_global(stage_a: dict, stage_a_null: pd.DataFrame, output_dir: Path) -> list[Path]:
    primary = stage_a["primary_global_result"]
    null_values = stage_a_null["global_equal_species_mean_q"].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1.0]})
    fig.suptitle("Stage A: continuous flower colour is locally organized within species", y=1.02, fontweight="bold")

    ax = axes[0]
    ax.hist(null_values, bins=42, density=True, color="#D9D9D9", edgecolor="white", linewidth=0.5)
    ax.axvline(primary["null_mean"], color=NULL_COLOUR, linestyle="--", linewidth=2.0, label="Null mean")
    ax.axvline(primary["observed"], color=PRIMARY_COLOUR, linewidth=2.5, label="Observed")
    ax.set_xlabel("Equal-species mean edge discontinuity, Q")
    ax.set_ylabel("Permutation density")
    ax.legend(frameon=False)
    ax.text(
        0.03,
        0.95,
        f"9,999 within-species permutations\nLower-tail p = {primary['p_lower_tail']:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.95},
    )
    panel_label(ax, "A")

    ax = axes[1]
    ks = [3, 5, 8]
    z = [stage_a["analyses_by_k"][str(k)]["global_equal_species_mean_q"]["standardized_clustering_deficit"] for k in ks]
    p = [stage_a["analyses_by_k"][str(k)]["global_equal_species_mean_q"]["p_lower_tail"] for k in ks]
    ax.plot(ks, z, color=NULL_COLOUR, linewidth=1.8, zorder=1)
    for k_value, z_value, p_value in zip(ks, z, p, strict=True):
        if k_value == 5:
            ax.scatter(k_value, z_value, marker="*", s=240, color=PRIMARY_COLOUR, edgecolor="black", linewidth=0.6, zorder=3)
        else:
            ax.scatter(k_value, z_value, marker="o", s=90, color=NULL_COLOUR, edgecolor="white", linewidth=0.8, zorder=2)
        ax.text(k_value, z_value + 0.12, f"p={p_value:.4f}", ha="center", va="bottom", fontsize=9)
    ax.axhline(0, color="#888888", linewidth=1.0)
    ax.set_xticks(ks, ["k=3", "k=5\nprimary", "k=8"])
    ax.set_ylabel("Standardized clustering deficit")
    ax.set_xlabel("Spherical nearest-neighbour degree")
    ax.set_ylim(min(-0.2, min(z) - 0.4), max(z) + 0.55)
    ax.text(
        0.03,
        0.05,
        "Positive values = neighbouring\nobservations are more colour-similar\nthan the species-conditioned null",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=NEUTRAL,
    )
    panel_label(ax, "B")

    fig.tight_layout()
    return save_pair(fig, output_dir, "jbi_ch1_figure_c1_stage_a_global")


def make_stage_a_species(stage_a: dict, output_dir: Path) -> list[Path]:
    primary_species = stage_a["analyses_by_k"]["5"]["species"]
    records = []
    for species in SPECIES:
        species_result = primary_species[species]
        q = species_result["q"]
        records.append(
            {
                "species": species,
                "z": float(q["standardized_clustering_deficit"]),
                "p": float(q["p_lower_tail"]),
                "edges": int(species_result["n_graph_edges"]),
            }
        )
    records.sort(key=lambda record: record["z"])

    z_values = np.asarray([record["z"] for record in records], dtype=float)
    edge_values = np.asarray([record["edges"] for record in records], dtype=float)
    if np.ptp(edge_values) > 0:
        sizes = 90 + 170 * (edge_values - edge_values.min()) / np.ptp(edge_values)
    else:
        sizes = np.full(len(records), 160.0)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    y = np.arange(len(records))
    ax.axvline(0, color="#888888", linewidth=1.0)
    global_z = float(stage_a["primary_global_result"]["standardized_clustering_deficit"])
    ax.axvline(global_z, color=NULL_COLOUR, linestyle="--", linewidth=1.4, alpha=0.8, label=f"Global equal-species deficit = {global_z:.2f}")

    for index, (record, size) in enumerate(zip(records, sizes, strict=True)):
        significant = record["p"] <= 0.05
        colour = SUPPORT_COLOUR if significant else ("#999999" if record["z"] < 0 else NULL_COLOUR)
        ax.scatter(record["z"], index, s=size, color=colour, edgecolor="white", linewidth=0.9, zorder=3)
        ax.text(
            max(z_values.max(), global_z) + 0.28,
            index,
            f"p={record['p']:.4f}",
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold" if significant else "normal",
        )

    ax.set_yticks(y, [SPECIES_ABBREVIATED[record["species"]] for record in records])
    ax.set_xlabel("Standardized clustering deficit")
    ax.set_title("Stage A species heterogeneity in the primary k=5 graph", fontweight="bold")
    ax.set_xlim(min(-1.25, z_values.min() - 0.45), max(z_values.max(), global_z) + 1.12)
    ax.legend(frameon=False, loc="lower right")
    ax.text(
        0.01,
        0.99,
        "Point size = retained graph edges\nGreen = species-level p ≤ 0.05",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color=NEUTRAL,
    )
    fig.tight_layout()
    return save_pair(fig, output_dir, "jbi_ch1_figure_c2_stage_a_species")


def draw_world_background(ax: plt.Axes) -> str:
    """Draw bundled Natural Earth coastlines when GeoPandas is available."""

    try:
        import geopandas as gpd

        path = gpd.datasets.get_path("naturalearth_lowres")
        world = gpd.read_file(path)
        world.plot(ax=ax, facecolor="#F5F5F5", edgecolor="#777777", linewidth=0.35, zorder=0)
        return "geopandas-naturalearth-lowres"
    except Exception:
        ax.set_facecolor("#FAFAFA")
        return "graticule-only-fallback"


def make_stage_b_surface(stage_b: dict, surface: pd.DataFrame, split: pd.DataFrame, output_dir: Path) -> tuple[list[Path], str]:
    selected = stage_b["selected_primary_configuration"]
    n_lon = int(selected["grid"]["n_lon"])
    n_sinlat = int(selected["grid"]["n_sinlat"])
    evaluation = split.loc[split["split"].astype(str).eq("evaluation")]

    evaluable = surface["evaluable_A_ge_minimum"].astype(str).str.lower().eq("true")
    cells = surface.loc[evaluable].copy()
    if len(cells) != int(stage_b["primary_result"]["surface"]["n_cells_evaluable_A_ge_minimum"]):
        raise ValueError("evaluable Stage B cell count differs between JSON and CSV")

    fig = plt.figure(figsize=(13.5, 6.4))
    grid_spec = fig.add_gridspec(2, 4, width_ratios=[1.25, 1.25, 1.25, 0.92], height_ratios=[1, 1])
    ax_map = fig.add_subplot(grid_spec[:, :3])
    ax_opportunity = fig.add_subplot(grid_spec[0, 3])
    ax_species = fig.add_subplot(grid_spec[1, 3])
    fig.suptitle("Stage B: observed shared-transition surface under the frozen primary support", y=0.995, fontweight="bold")

    basemap_source = draw_world_background(ax_map)
    ax_map.scatter(
        evaluation["longitude"],
        evaluation["latitude"],
        s=6,
        color="#666666",
        alpha=0.22,
        linewidths=0,
        zorder=1,
        label="720 evaluation locations",
    )

    cmap = plt.get_cmap("viridis")
    norm = Normalize(vmin=0.0, vmax=1.0)
    lon_width = 360.0 / n_lon
    for row in cells.itertuples(index=False):
        cell_id = int(row.cell_id)
        row_index = cell_id // n_lon
        col_index = cell_id % n_lon
        lon0 = -180.0 + col_index * lon_width
        sin0 = -1.0 + row_index * (2.0 / n_sinlat)
        sin1 = -1.0 + (row_index + 1) * (2.0 / n_sinlat)
        lat0 = float(np.rad2deg(np.arcsin(np.clip(sin0, -1.0, 1.0))))
        lat1 = float(np.rad2deg(np.arcsin(np.clip(sin1, -1.0, 1.0))))
        intensity = float(row.shared_transition_intensity)
        ax_map.add_patch(
            Rectangle(
                (lon0, lat0),
                lon_width,
                lat1 - lat0,
                facecolor=cmap(norm(intensity)),
                edgecolor="#222222",
                linewidth=0.45,
                alpha=0.88,
                zorder=2,
            )
        )

    ax_map.scatter(
        cells["longitude"],
        cells["latitude"],
        s=18 + 16 * cells["opportunity_A"].to_numpy(dtype=float),
        facecolors="none",
        edgecolors="black",
        linewidths=0.7,
        zorder=3,
    )
    ax_map.set_xlim(-180, 180)
    ax_map.set_ylim(-72, 82)
    ax_map.set_xticks(np.arange(-180, 181, 60))
    ax_map.set_yticks(np.arange(-60, 81, 30))
    ax_map.grid(color="#CCCCCC", linewidth=0.45, linestyle=":", zorder=-1)
    ax_map.set_xlabel("Longitude")
    ax_map.set_ylabel("Latitude")
    ax_map.legend(frameon=False, loc="lower left")
    ax_map.text(
        0.99,
        0.02,
        "Cell fill = shared transition intensity S(x)\nCircle size = opportunity A(x); absent cells are not evaluable",
        transform=ax_map.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.93},
    )
    colourbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax_map, fraction=0.03, pad=0.015)
    colourbar.set_label("Shared transition intensity, S(x)")
    panel_label(ax_map, "A")

    surface_summary = stage_b["primary_result"]["surface"]
    opportunity_labels = ["A≥2", "A≥3", "A≥4"]
    opportunity_counts = [
        int(surface_summary["n_cells_A_ge_2"]),
        int(surface_summary["n_cells_A_ge_3"]),
        int(surface_summary["n_cells_A_ge_4"]),
    ]
    bars = ax_opportunity.bar(opportunity_labels, opportunity_counts, color=["#56B4E9", "#0072B2", "#003B5C"])
    for bar, count in zip(bars, opportunity_counts, strict=True):
        ax_opportunity.text(bar.get_x() + bar.get_width() / 2, count + 0.5, str(count), ha="center", va="bottom")
    ax_opportunity.set_ylabel("Number of cells")
    ax_opportunity.set_title("Shared opportunity")
    ax_opportunity.set_ylim(0, max(opportunity_counts) + 5)
    ax_opportunity.text(
        0.98,
        0.95,
        f"max A = {int(surface_summary['maximum_A'])}",
        transform=ax_opportunity.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )
    panel_label(ax_opportunity, "B")

    species_counts = selected["detectable_cells_by_species"]
    y = np.arange(len(SPECIES))
    values = [int(species_counts[species]) for species in SPECIES]
    ax_species.barh(y, values, color="#999999")
    ax_species.set_yticks(y, [SPECIES_ABBREVIATED[species] for species in SPECIES])
    ax_species.invert_yaxis()
    ax_species.set_xlabel("Detectable cells")
    ax_species.set_title("Geometry support by species")
    for index, value in enumerate(values):
        ax_species.text(value + 0.25, index, str(value), va="center", ha="left", fontsize=8.5)
    ax_species.set_xlim(0, max(values) + 4)
    panel_label(ax_species, "C")

    fig.tight_layout(rect=[0, 0, 1, 0.975])
    return save_pair(fig, output_dir, "jbi_ch1_figure_c3_stage_b_surface"), basemap_source


def parse_configuration_key(key: str) -> tuple[int, str]:
    match = re.fullmatch(r"cap_(\d+(?:p\d+)?)km_grid_(\d+)x(\d+)", key)
    if match is None:
        raise ValueError(f"unrecognized Stage B configuration key: {key}")
    cap = int(float(match.group(1).replace("p", ".")))
    return cap, f"{match.group(2)}×{match.group(3)}"


def stage_b_configuration_records(stage_b: dict) -> list[dict[str, object]]:
    selected = stage_b["selected_primary_configuration"]
    selected_key = str(selected["configuration"])
    selected_grid = selected["grid"]
    records: list[dict[str, object]] = [
        {
            "key": selected_key,
            "cap": int(float(selected["max_edge_km"])),
            "grid": f"{int(selected_grid['n_lon'])}×{int(selected_grid['n_sinlat'])}",
            "p": float(stage_b["primary_result"]["global_concentration"]["p_upper_tail"]),
            "z": float(stage_b["primary_result"]["global_concentration"]["standardized_concentration_excess"]),
            "primary": True,
        }
    ]
    for key, result in stage_b["sensitivity_results"].items():
        cap, grid = parse_configuration_key(key)
        global_result = result["global_concentration"]
        records.append(
            {
                "key": key,
                "cap": cap,
                "grid": grid,
                "p": float(global_result["p_upper_tail"]),
                "z": float(global_result["standardized_concentration_excess"]),
                "primary": False,
            }
        )
    if len(records) != 9 or sum(bool(record["primary"]) for record in records) != 1:
        raise ValueError("expected exactly nine Stage B configurations with one primary")
    return records


def make_stage_b_sensitivity(stage_b: dict, output_dir: Path) -> list[Path]:
    records = stage_b_configuration_records(stage_b)
    grid_order = ["36×18", "24×12", "18×9"]
    grid_x = {label: index for index, label in enumerate(grid_order)}

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8))
    fig.suptitle("Stage B: shared-transition concentration depends on spatial support", y=1.02, fontweight="bold")

    for cap in [500, 1000, 2000]:
        cap_records = sorted((record for record in records if record["cap"] == cap), key=lambda record: grid_x[str(record["grid"])])
        x = [grid_x[str(record["grid"])] for record in cap_records]
        p = [float(record["p"]) for record in cap_records]
        z = [float(record["z"]) for record in cap_records]
        axes[0].plot(x, p, marker="o", linewidth=1.8, color=CAP_COLOURS[cap], label=f"{cap:,} km cap")
        axes[1].plot(x, z, marker="o", linewidth=1.8, color=CAP_COLOURS[cap], label=f"{cap:,} km cap")
        for x_value, p_value in zip(x, p, strict=True):
            axes[0].text(x_value, p_value + 0.018, f"{p_value:.4f}", ha="center", va="bottom", fontsize=7.8, color=CAP_COLOURS[cap])

    primary = next(record for record in records if record["primary"])
    primary_x = grid_x[str(primary["grid"])]
    axes[0].scatter(primary_x, float(primary["p"]), marker="*", s=260, color=PRIMARY_COLOUR, edgecolor="black", linewidth=0.7, zorder=5)
    axes[1].scatter(primary_x, float(primary["z"]), marker="*", s=260, color=PRIMARY_COLOUR, edgecolor="black", linewidth=0.7, zorder=5)

    axes[0].axhline(0.05, color="#444444", linestyle="--", linewidth=1.1)
    axes[0].set_xticks(range(len(grid_order)), grid_order)
    axes[0].set_xlabel("Equal-area grid")
    axes[0].set_ylabel("Monte Carlo upper-tail p")
    axes[0].set_ylim(0, max(float(record["p"]) for record in records) + 0.12)
    axes[0].text(0.02, 0.96, "Star = selected primary", transform=axes[0].transAxes, ha="left", va="top", fontsize=8.5)
    panel_label(axes[0], "A")

    axes[1].axhline(0, color="#888888", linewidth=1.0)
    axes[1].set_xticks(range(len(grid_order)), grid_order)
    axes[1].set_xlabel("Equal-area grid")
    axes[1].set_ylabel("Standardized concentration excess")
    axes[1].legend(frameon=False, loc="best")
    axes[1].text(
        0.02,
        0.04,
        "Positive values indicate more concentrated shared intensity\nthan the species-conditioned null.",
        transform=axes[1].transAxes,
        ha="left",
        va="bottom",
        fontsize=8.3,
        color=NEUTRAL,
    )
    panel_label(axes[1], "B")

    fig.tight_layout()
    return save_pair(fig, output_dir, "jbi_ch1_figure_c4_stage_b_sensitivity")


def draw_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    *,
    facecolor: str,
    edgecolor: str = "#555555",
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=1.2,
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.70, title, transform=ax.transAxes, ha="center", va="center", fontweight="bold", fontsize=10)
    ax.text(x + width / 2, y + height * 0.34, body, transform=ax.transAxes, ha="center", va="center", fontsize=8.7, linespacing=1.25)


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        transform=ax.transAxes,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.2,
        color="#555555",
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(arrow)


def make_workflow_figure(stage_a: dict, stage_b: dict, output_dir: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(13.5, 6.8))
    ax.set_axis_off()
    ax.set_title("Frozen Chapter 1 workflow and inferential gates", fontsize=15, fontweight="bold", pad=12)

    top_y = 0.62
    top_w = 0.19
    top_h = 0.22
    top_x = [0.02, 0.265, 0.51, 0.755]
    draw_box(ax, top_x[0], top_y, top_w, top_h, "Acquisition", "1,200 photographs\n6 species × 200", facecolor="#F2F2F2")
    draw_box(ax, top_x[1], top_y, top_w, top_h, "Outcome-blind split", "480 calibration\n720 held-out evaluation", facecolor="#E8F1F8")
    draw_box(ax, top_x[2], top_y, top_w, top_h, "Representation frozen", "Species-specific continuous\ncolour vectors; calibration scaling", facecolor="#E8F1F8")
    draw_box(ax, top_x[3], top_y, top_w, top_h, "Evaluation opened", "720/720 measured\n0 localization failures", facecolor="#E7F4EC")
    for index in range(3):
        draw_arrow(ax, (top_x[index] + top_w, top_y + top_h / 2), (top_x[index + 1], top_y + top_h / 2))

    bottom_y = 0.18
    bottom_w = 0.22
    bottom_h = 0.25
    bottom_x = [0.12, 0.39, 0.66]
    stage_a_p = float(stage_a["primary_global_result"]["p_lower_tail"])
    stage_b_p = float(stage_b["primary_result"]["global_concentration"]["p_upper_tail"])
    draw_box(
        ax,
        bottom_x[0],
        bottom_y,
        bottom_w,
        bottom_h,
        "Stage A gate: passed",
        f"Within-species local organization\nk=5, 9,999 permutations\np={stage_a_p:.4f}",
        facecolor="#DDF0E4",
        edgecolor="#2E7D4B",
    )
    draw_box(
        ax,
        bottom_x[1],
        bottom_y,
        bottom_w,
        bottom_h,
        "Stage B gate: not passed",
        f"Shared-transition concentration\n500 km / 36×18 primary\np={stage_b_p:.4f}",
        facecolor="#FCE8D8",
        edgecolor="#B85C00",
    )
    draw_box(
        ax,
        bottom_x[2],
        bottom_y,
        bottom_w,
        bottom_h,
        "Claim stops here",
        "Retain repeated local organization\nDo not promote one shared boundary\nor geographic-cause correspondence",
        facecolor="#F2F2F2",
    )
    draw_arrow(ax, (top_x[3] + top_w / 2, top_y), (bottom_x[0] + bottom_w / 2, bottom_y + bottom_h))
    draw_arrow(ax, (bottom_x[0] + bottom_w, bottom_y + bottom_h / 2), (bottom_x[1], bottom_y + bottom_h / 2))
    draw_arrow(ax, (bottom_x[1] + bottom_w, bottom_y + bottom_h / 2), (bottom_x[2], bottom_y + bottom_h / 2))

    ax.text(
        0.5,
        0.06,
        "Complete colour vectors were permuted strictly within species; locations and all label-blind geometry remained fixed.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=9.5,
        color=NEUTRAL,
    )
    fig.tight_layout()
    return save_pair(fig, output_dir, "jbi_ch1_figure_cs1_workflow")


def make_detectability_figure(stage_b: dict, surface: pd.DataFrame, output_dir: Path) -> list[Path]:
    evaluable = surface["evaluable_A_ge_minimum"].astype(str).str.lower().eq("true")
    cells = surface.loc[evaluable].copy().sort_values(["latitude", "longitude", "cell_id"]).reset_index(drop=True)
    if len(cells) == 0:
        raise ValueError("primary Stage B surface has no evaluable cells")

    matrix = np.empty((len(SPECIES), len(cells)), dtype=float)
    for species_index, species in enumerate(SPECIES):
        slug = SPECIES_SLUG[species]
        require_columns(
            cells,
            [f"{slug}_detectable", f"{slug}_transition_intensity"],
            label="Stage B surface",
        )
        detectable = cells[f"{slug}_detectable"].astype(str).str.lower().eq("true").to_numpy()
        values = pd.to_numeric(cells[f"{slug}_transition_intensity"], errors="coerce").to_numpy(dtype=float)
        values[~detectable] = np.nan
        matrix[species_index] = values

    labels = [
        f"{float(row.longitude):g}°\n{float(row.latitude):g}°"
        for row in cells.itertuples(index=False)
    ]
    opportunity = cells["opportunity_A"].to_numpy(dtype=int)

    fig = plt.figure(figsize=(14.2, 6.5))
    grid_spec = fig.add_gridspec(2, 1, height_ratios=[1.0, 4.2], hspace=0.08)
    ax_top = fig.add_subplot(grid_spec[0, 0])
    ax_heat = fig.add_subplot(grid_spec[1, 0], sharex=ax_top)
    fig.suptitle("Supporting Figure: detectability behind the primary shared-transition surface", y=0.995, fontweight="bold")

    x = np.arange(len(cells))
    ax_top.bar(x, opportunity, color="#777777", width=0.85)
    ax_top.set_ylabel("A")
    ax_top.set_ylim(0, max(6.7, opportunity.max() + 0.7))
    ax_top.set_title("Opportunity denominator by geographically ordered evaluable cell", fontsize=11)
    ax_top.tick_params(axis="x", labelbottom=False)
    for index, value in enumerate(opportunity):
        ax_top.text(index, value + 0.1, str(value), ha="center", va="bottom", fontsize=7.5)
    panel_label(ax_top, "A")

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#E5E5E5")
    masked = np.ma.masked_invalid(matrix)
    image = ax_heat.imshow(masked, aspect="auto", interpolation="none", vmin=0, vmax=1, cmap=cmap)
    ax_heat.set_yticks(np.arange(len(SPECIES)), [SPECIES_ABBREVIATED[species] for species in SPECIES])
    ax_heat.set_xticks(x, labels, rotation=90, ha="center", fontsize=6.8)
    ax_heat.set_xlabel("Cell centre (longitude / latitude); ordered by latitude then longitude")
    ax_heat.set_title("Species-cell transition intensity; grey = species not detectable in that cell", fontsize=11)
    colourbar = fig.colorbar(image, ax=ax_heat, fraction=0.025, pad=0.012)
    colourbar.set_label("Within-species ranked transition intensity")
    panel_label(ax_heat, "B")

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return save_pair(fig, output_dir, "jbi_ch1_figure_cs2_detectability")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-a", type=Path, default=Path("docs/supporting/jbi_ch1_stage_a_continuous_graph_v1.json"))
    parser.add_argument("--stage-a-null", type=Path, default=Path("data/evaluation/jbi_ch1_stage_a_primary_null_v1.csv"))
    parser.add_argument("--stage-b", type=Path, default=Path("docs/supporting/jbi_ch1_stage_b_shared_transition_concentration_v1.json"))
    parser.add_argument("--stage-b-surface", type=Path, default=Path("data/evaluation/jbi_ch1_stage_b_shared_transition_surface_v1.csv"))
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("docs/figures"))
    parser.add_argument("--manifest", type=Path, default=Path("docs/supporting/jbi_ch1_figure_manifest_v1.json"))
    args = parser.parse_args()

    configure_matplotlib()
    stage_a = load_json(args.stage_a)
    stage_b = load_json(args.stage_b)
    stage_a_null = pd.read_csv(args.stage_a_null)
    surface = pd.read_csv(args.stage_b_surface)
    split = pd.read_csv(args.split)
    validate_inputs(stage_a, stage_b, stage_a_null, surface, split)

    outputs: list[Path] = []
    outputs.extend(make_stage_a_global(stage_a, stage_a_null, args.output_dir))
    outputs.extend(make_stage_a_species(stage_a, args.output_dir))
    surface_outputs, basemap_source = make_stage_b_surface(stage_b, surface, split, args.output_dir)
    outputs.extend(surface_outputs)
    outputs.extend(make_stage_b_sensitivity(stage_b, args.output_dir))
    outputs.extend(make_workflow_figure(stage_a, stage_b, args.output_dir))
    outputs.extend(make_detectability_figure(stage_b, surface, args.output_dir))

    for path in outputs:
        if not path.is_file() or path.stat().st_size < 1_000:
            raise ValueError(f"figure output missing or implausibly small: {path}")

    manifest = {
        "protocol": "jbi-ch1-spatial-figure-manifest-v1",
        "status": "canonical_figures_generated_from_frozen_stage_a_and_stage_b",
        "inputs": {
            str(args.stage_a): sha256(args.stage_a),
            str(args.stage_a_null): sha256(args.stage_a_null),
            str(args.stage_b): sha256(args.stage_b),
            str(args.stage_b_surface): sha256(args.stage_b_surface),
            str(args.split): sha256(args.split),
        },
        "results": {
            "stage_a_primary_lower_tail_p": float(stage_a["primary_global_result"]["p_lower_tail"]),
            "stage_a_primary_rejects_random_labelling_at_0_05": bool(stage_a["primary_rejects_random_labelling_at_0_05"]),
            "stage_b_primary_upper_tail_p": float(stage_b["primary_result"]["global_concentration"]["p_upper_tail"]),
            "stage_b_primary_rejects_shared_concentration_at_0_05": bool(stage_b["primary_rejects_shared_concentration_null_at_0_05"]),
        },
        "basemap_source": basemap_source,
        "outputs": [
            {
                "path": str(path),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in sorted(outputs, key=lambda item: str(item))
        ],
        "environment_used_for_inference": False,
        "geographic_reference_library_used_for_inference": False,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
