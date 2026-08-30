#!/usr/bin/env python3
"""Render the canonical Chapter 1 figures with visual-QA layout corrections.

The underlying figure pipeline and all frozen numerical inputs remain unchanged.  This
wrapper replaces only three plotting functions whose first render had overlapping text:
Stage-A species labels, Stage-B sensitivity annotations, and the detectability-panel
heading.  The base module still performs all input validation, output hashing, and
manifest generation.
"""
from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

BASE_PATH = Path(__file__).with_name("make_jbi_ch1_spatial_figures.py")
SPEC = spec_from_file_location("jbi_ch1_spatial_figures_base", BASE_PATH)
assert SPEC is not None and SPEC.loader is not None
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def make_stage_a_species(stage_a: dict, output_dir: Path) -> list[Path]:
    """Render species heterogeneity without a legend/p-value collision."""

    primary_species = stage_a["analyses_by_k"]["5"]["species"]
    records = []
    for species in BASE.SPECIES:
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

    fig, ax = plt.subplots(figsize=(8.4, 5.3))
    y = np.arange(len(records))
    ax.axvline(0, color="#888888", linewidth=1.0)
    global_z = float(stage_a["primary_global_result"]["standardized_clustering_deficit"])
    ax.axvline(global_z, color=BASE.NULL_COLOUR, linestyle="--", linewidth=1.4, alpha=0.8)

    p_text_x = max(z_values.max(), global_z) + 0.28
    for index, (record, size) in enumerate(zip(records, sizes, strict=True)):
        significant = record["p"] <= 0.05
        colour = BASE.SUPPORT_COLOUR if significant else ("#999999" if record["z"] < 0 else BASE.NULL_COLOUR)
        ax.scatter(record["z"], index, s=size, color=colour, edgecolor="white", linewidth=0.9, zorder=3)
        ax.text(
            p_text_x,
            index,
            f"p={record['p']:.4f}",
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold" if significant else "normal",
        )

    ax.set_yticks(y, [BASE.SPECIES_ABBREVIATED[record["species"]] for record in records])
    ax.set_xlabel("Standardized clustering deficit")
    ax.set_title("Stage A species heterogeneity in the primary k=5 graph", fontweight="bold")
    ax.set_xlim(min(-1.25, z_values.min() - 0.45), p_text_x + 0.72)
    ax.text(
        0.01,
        0.99,
        f"Point size = retained graph edges\nGreen = species-level p ≤ 0.05\nBlue dashed = global equal-species deficit ({global_z:.2f})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.4,
        color=BASE.NEUTRAL,
        linespacing=1.25,
    )
    fig.tight_layout()
    return BASE.save_pair(fig, output_dir, "jbi_ch1_figure_c2_stage_a_species")


def make_stage_b_sensitivity(stage_b: dict, output_dir: Path) -> list[Path]:
    """Render scale sensitivity with separated annotations and an external note."""

    records = BASE.stage_b_configuration_records(stage_b)
    grid_order = ["36×18", "24×12", "18×9"]
    grid_x = {label: index for index, label in enumerate(grid_order)}

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))
    fig.suptitle("Stage B: shared-transition concentration depends on spatial support", y=1.01, fontweight="bold")

    for cap in [500, 1000, 2000]:
        cap_records = sorted(
            (record for record in records if record["cap"] == cap),
            key=lambda record: grid_x[str(record["grid"])],
        )
        x = [grid_x[str(record["grid"])] for record in cap_records]
        p = [float(record["p"]) for record in cap_records]
        z = [float(record["z"]) for record in cap_records]
        axes[0].plot(x, p, marker="o", linewidth=1.8, color=BASE.CAP_COLOURS[cap], label=f"{cap:,} km cap")
        axes[1].plot(x, z, marker="o", linewidth=1.8, color=BASE.CAP_COLOURS[cap], label=f"{cap:,} km cap")

        for record, x_value, p_value in zip(cap_records, x, p, strict=True):
            offset = 0.018
            va = "bottom"
            # The 500- and 2,000-km points at 18×9 are nearly coincident; put
            # their labels on opposite sides of the points.
            if record["grid"] == "18×9" and cap == 500:
                offset = -0.025
                va = "top"
            elif record["grid"] == "18×9" and cap == 2000:
                offset = 0.030
            axes[0].text(
                x_value,
                p_value + offset,
                f"{p_value:.4f}",
                ha="center",
                va=va,
                fontsize=7.8,
                color=BASE.CAP_COLOURS[cap],
            )

    primary = next(record for record in records if record["primary"])
    primary_x = grid_x[str(primary["grid"])]
    axes[0].scatter(
        primary_x,
        float(primary["p"]),
        marker="*",
        s=260,
        color=BASE.PRIMARY_COLOUR,
        edgecolor="black",
        linewidth=0.7,
        zorder=5,
    )
    axes[1].scatter(
        primary_x,
        float(primary["z"]),
        marker="*",
        s=260,
        color=BASE.PRIMARY_COLOUR,
        edgecolor="black",
        linewidth=0.7,
        zorder=5,
    )

    axes[0].axhline(0.05, color="#444444", linestyle="--", linewidth=1.1)
    axes[0].set_xticks(range(len(grid_order)), grid_order)
    axes[0].set_xlabel("Equal-area grid")
    axes[0].set_ylabel("Monte Carlo upper-tail p")
    axes[0].set_ylim(0, max(float(record["p"]) for record in records) + 0.12)
    axes[0].text(0.02, 0.96, "Star = selected primary", transform=axes[0].transAxes, ha="left", va="top", fontsize=8.5)
    BASE.panel_label(axes[0], "A")

    axes[1].axhline(0, color="#888888", linewidth=1.0)
    axes[1].set_xticks(range(len(grid_order)), grid_order)
    axes[1].set_xlabel("Equal-area grid")
    axes[1].set_ylabel("Standardized concentration excess")
    axes[1].legend(frameon=False, loc="best")
    axes[1].text(
        0.02,
        -0.16,
        "Positive values indicate more concentrated shared intensity\nthan the species-conditioned null.",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=8.3,
        color=BASE.NEUTRAL,
        clip_on=False,
    )
    BASE.panel_label(axes[1], "B")

    fig.tight_layout(rect=[0, 0.08, 1, 1])
    return BASE.save_pair(fig, output_dir, "jbi_ch1_figure_c4_stage_b_sensitivity")


def make_detectability_figure(stage_b: dict, surface, output_dir: Path) -> list[Path]:
    """Render the opportunity/heatmap panels with a larger inter-panel gap."""

    evaluable = surface["evaluable_A_ge_minimum"].astype(str).str.lower().eq("true")
    cells = surface.loc[evaluable].copy().sort_values(["latitude", "longitude", "cell_id"]).reset_index(drop=True)
    if len(cells) == 0:
        raise ValueError("primary Stage B surface has no evaluable cells")

    matrix = np.empty((len(BASE.SPECIES), len(cells)), dtype=float)
    for species_index, species in enumerate(BASE.SPECIES):
        slug = BASE.SPECIES_SLUG[species]
        BASE.require_columns(
            cells,
            [f"{slug}_detectable", f"{slug}_transition_intensity"],
            label="Stage B surface",
        )
        detectable = cells[f"{slug}_detectable"].astype(str).str.lower().eq("true").to_numpy()
        values = BASE.pd.to_numeric(cells[f"{slug}_transition_intensity"], errors="coerce").to_numpy(dtype=float)
        values[~detectable] = np.nan
        matrix[species_index] = values

    labels = [
        f"{float(row.longitude):g}°\n{float(row.latitude):g}°"
        for row in cells.itertuples(index=False)
    ]
    opportunity = cells["opportunity_A"].to_numpy(dtype=int)

    fig = plt.figure(figsize=(14.2, 7.1))
    grid_spec = fig.add_gridspec(2, 1, height_ratios=[1.0, 4.2], hspace=0.30)
    ax_top = fig.add_subplot(grid_spec[0, 0])
    ax_heat = fig.add_subplot(grid_spec[1, 0], sharex=ax_top)
    fig.suptitle(
        "Supporting Figure: detectability behind the primary shared-transition surface",
        y=0.995,
        fontweight="bold",
    )

    x = np.arange(len(cells))
    ax_top.bar(x, opportunity, color="#777777", width=0.85)
    ax_top.set_ylabel("A")
    ax_top.set_ylim(0, max(6.7, opportunity.max() + 0.7))
    ax_top.set_title("Opportunity denominator by geographically ordered evaluable cell", fontsize=11, pad=8)
    ax_top.tick_params(axis="x", labelbottom=False)
    for index, value in enumerate(opportunity):
        ax_top.text(index, value + 0.1, str(value), ha="center", va="bottom", fontsize=7.5)
    BASE.panel_label(ax_top, "A")

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#E5E5E5")
    masked = np.ma.masked_invalid(matrix)
    image = ax_heat.imshow(masked, aspect="auto", interpolation="none", vmin=0, vmax=1, cmap=cmap)
    ax_heat.set_yticks(np.arange(len(BASE.SPECIES)), [BASE.SPECIES_ABBREVIATED[species] for species in BASE.SPECIES])
    ax_heat.set_xticks(x, labels, rotation=90, ha="center", fontsize=6.8)
    ax_heat.set_xlabel("Cell centre (longitude / latitude); ordered by latitude then longitude")
    ax_heat.set_title(
        "Species-cell transition intensity; grey = species not detectable in that cell",
        fontsize=11,
        pad=10,
    )
    colourbar = fig.colorbar(image, ax=ax_heat, fraction=0.025, pad=0.012)
    colourbar.set_label("Within-species ranked transition intensity")
    BASE.panel_label(ax_heat, "B")

    fig.tight_layout(rect=[0, 0, 1, 0.965], h_pad=2.2)
    return BASE.save_pair(fig, output_dir, "jbi_ch1_figure_cs2_detectability")


BASE.make_stage_a_species = make_stage_a_species
BASE.make_stage_b_sensitivity = make_stage_b_sensitivity
BASE.make_detectability_figure = make_detectability_figure


if __name__ == "__main__":
    raise SystemExit(BASE.main())
