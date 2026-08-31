"""Colour-blind environmental boundary surfaces for the FCP atlas."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np


ENVIRONMENT_PROTOCOL = "jbi-atlas-environmental-overlay-freeze-v1"
CLIMATE_FIELDS = ("bio1", "bio4", "bio12", "bio15")
TERRAIN_FIELDS = ("elevation", "slope", "terrain_ruggedness")
LAND_COVER_CODES = (10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100)
PRIMARY_FAMILIES = ("macroclimate", "terrain", "land_cover", "ecoregion")
BOUNDARY_FIELD_BY_FAMILY = {
    "macroclimate": "macroclimate_boundary",
    "terrain": "terrain_boundary",
    "land_cover": "land_cover_boundary",
    "ecoregion": "ecoregion_boundary",
}


def validate_environment_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != ENVIRONMENT_PROTOCOL:
        raise ValueError("unexpected environmental overlay protocol")
    if contract.get("scaleout_colour_opened") is not False:
        raise ValueError("environmental overlay was not frozen before colour")
    if contract.get("grid", {}).get("scales_km") != [100, 250, 500]:
        raise ValueError("environmental overlay scales changed")
    coverage = contract.get("coverage_gate", {})
    if coverage.get("minimum_atlas_opportunity_cell_fraction_per_family") != 0.7:
        raise ValueError("environmental family coverage changed")
    if coverage.get("minimum_evaluable_primary_families") != 2:
        raise ValueError("environmental family minimum changed")
    if coverage.get("macroclimate_must_be_evaluable") is not True:
        raise ValueError("macroclimate must remain required")
    if contract.get("inference", {}).get("randomizations") != 9999:
        raise ValueError("final environmental randomization count changed")
    if contract.get("worldcover_sampling", {}).get("minimum_boundary_spearman") != 0.9:
        raise ValueError("WorldCover sampling-stability threshold changed")


def evaluate_environmental_coverage_gate(
    geometry_audit: Sequence[Mapping[str, Any]],
    selected_taxon_ids: Sequence[str | int],
    boundary_rows_by_scale: Mapping[int, Sequence[Mapping[str, Any]]],
    contract: Mapping[str, Any],
    *,
    families: Sequence[str] = ("macroclimate", "land_cover", "ecoregion"),
) -> dict[str, Any]:
    """Evaluate pre-colour environmental coverage on transition-opportunity cells.

    Opportunity is the union of detectable cells for the frozen species at each
    scale.  Missing environmental cells remain missing and never become zero.
    """

    validate_environment_contract(contract)
    selected_values = tuple(str(value) for value in selected_taxon_ids)
    selected = set(selected_values)
    if not selected or len(selected) != len(selected_values):
        raise ValueError("selected taxon IDs must be non-empty and unique")
    families = tuple(str(value) for value in families)
    if (
        not families
        or len(families) != len(set(families))
        or not set(families).issubset(BOUNDARY_FIELD_BY_FAMILY)
    ):
        raise ValueError("environmental coverage families are invalid")
    geometry_by_taxon = {str(row.get("taxon_id", "")): row for row in geometry_audit}
    if len(geometry_by_taxon) != len(geometry_audit) or "" in geometry_by_taxon:
        raise ValueError("geometry audit contains duplicate taxon IDs")
    missing_taxa = selected - set(geometry_by_taxon)
    if missing_taxa:
        raise ValueError(f"selected taxa lack geometry audit rows: {sorted(missing_taxa)}")

    minimum_fraction = float(
        contract["coverage_gate"]["minimum_atlas_opportunity_cell_fraction_per_family"]
    )
    required_scales = tuple(int(value) for value in contract["grid"]["scales_km"])
    results_by_scale: list[dict[str, Any]] = []
    family_pass_by_scale: dict[str, list[bool]] = {family: [] for family in families}
    for scale in required_scales:
        opportunity: set[int] = set()
        for taxon_id in selected:
            matches = [
                row
                for row in geometry_by_taxon[taxon_id].get("scale_results", ())
                if int(row.get("scale_km", -1)) == scale
            ]
            if len(matches) != 1:
                raise ValueError(f"taxon {taxon_id} lacks one geometry row at {scale} km")
            cell_ids = [int(value) for value in matches[0].get("detectable_cell_ids", ())]
            if len(cell_ids) != len(set(cell_ids)):
                raise ValueError("detectable cell IDs must be unique within a species-scale")
            opportunity.update(cell_ids)
        if not opportunity:
            raise ValueError(f"no transition-opportunity cells at {scale} km")

        rows = list(boundary_rows_by_scale.get(scale, ()))
        by_cell: dict[int, Mapping[str, Any]] = {}
        for row in rows:
            cell = int(row["cell_id"])
            if cell in by_cell:
                raise ValueError(f"duplicate environmental cell {cell} at {scale} km")
            by_cell[cell] = row
        family_results: dict[str, Any] = {}
        for family in families:
            field = BOUNDARY_FIELD_BY_FAMILY[family]
            finite_cells = 0
            for cell in opportunity:
                value = by_cell.get(cell, {}).get(field, "")
                try:
                    finite = bool(value != "" and np.isfinite(float(value)))
                except (TypeError, ValueError):
                    finite = False
                finite_cells += int(finite)
            fraction = finite_cells / len(opportunity)
            passed = fraction >= minimum_fraction
            family_pass_by_scale[family].append(passed)
            family_results[family] = {
                "finite_opportunity_cells": finite_cells,
                "coverage_fraction": fraction,
                "passes": passed,
            }
        results_by_scale.append(
            {
                "scale_km": scale,
                "opportunity_cells": len(opportunity),
                "families": family_results,
            }
        )

    evaluable_families = [
        family for family in families if all(family_pass_by_scale[family])
    ]
    minimum_families = int(contract["coverage_gate"]["minimum_evaluable_primary_families"])
    macroclimate_required = bool(
        contract["coverage_gate"]["macroclimate_must_be_evaluable"]
    )
    passed = (
        len(evaluable_families) >= minimum_families
        and (not macroclimate_required or "macroclimate" in evaluable_families)
    )
    return {
        "protocol": ENVIRONMENT_PROTOCOL,
        "status": (
            "pass_precolour_environmental_coverage"
            if passed
            else "not_evaluable_precolour_environmental_coverage"
        ),
        "scaleout_colour_opened": False,
        "selected_species": len(selected),
        "minimum_coverage_fraction": minimum_fraction,
        "evaluable_families": evaluable_families,
        "not_evaluable_families": [
            family for family in families if family not in evaluable_families
        ],
        "scales": results_by_scale,
    }


def rook_adjacency_without_repair(
    cell_ids: Sequence[int] | np.ndarray, *, n_lon: int, n_sinlat: int
) -> np.ndarray:
    """Return global-grid rook adjacency without joining distant islands."""

    ids = np.asarray(cell_ids, dtype=int)
    if ids.ndim != 1 or not ids.size or len(np.unique(ids)) != ids.size:
        raise ValueError("cell_ids must be non-empty and unique")
    if n_lon < 2 or n_sinlat < 2 or np.any((ids < 0) | (ids >= n_lon * n_sinlat)):
        raise ValueError("cell IDs or grid dimensions are invalid")
    lookup = {int(cell): index for index, cell in enumerate(ids)}
    adjacency = np.zeros((ids.size, ids.size), dtype=bool)
    for index, cell in enumerate(ids):
        row, column = divmod(int(cell), n_lon)
        neighbours = [
            row * n_lon + (column - 1) % n_lon,
            row * n_lon + (column + 1) % n_lon,
        ]
        if row > 0:
            neighbours.append((row - 1) * n_lon + column)
        if row + 1 < n_sinlat:
            neighbours.append((row + 1) * n_lon + column)
        for neighbour in neighbours:
            other = lookup.get(neighbour)
            if other is not None:
                adjacency[index, other] = True
                adjacency[other, index] = True
    return adjacency


def _edges(adjacency: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    adjacency = np.asarray(adjacency, dtype=bool)
    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError("adjacency must be square")
    if not np.array_equal(adjacency, adjacency.T) or np.any(np.diag(adjacency)):
        raise ValueError("adjacency must be symmetric without self edges")
    return np.nonzero(np.triu(adjacency, k=1))


def cell_mean_edge_values(adjacency: np.ndarray, edge_values: np.ndarray) -> np.ndarray:
    """Average undirected edge changes at each incident cell."""

    first, second = _edges(adjacency)
    edge_values = np.asarray(edge_values, dtype=float)
    if edge_values.shape != first.shape or not np.isfinite(edge_values).all():
        raise ValueError("edge values must be finite and match adjacency edges")
    total = np.zeros(adjacency.shape[0], dtype=float)
    count = np.zeros(adjacency.shape[0], dtype=int)
    np.add.at(total, first, edge_values)
    np.add.at(total, second, edge_values)
    np.add.at(count, first, 1)
    np.add.at(count, second, 1)
    output = np.full(adjacency.shape[0], np.nan, dtype=float)
    valid = count > 0
    output[valid] = total[valid] / count[valid]
    return output


def continuous_boundary(
    values: np.ndarray, adjacency: np.ndarray
) -> tuple[np.ndarray, dict[str, list[float]]]:
    """Standardized RMS neighbour change and its frozen land-cell scaling."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] != adjacency.shape[0]:
        raise ValueError("continuous values must have one row per cell")
    if not np.isfinite(values).all():
        raise ValueError("continuous environmental values must be complete and finite")
    mean = values.mean(axis=0)
    sd = values.std(axis=0)
    if np.any(sd <= 0):
        raise ValueError("continuous environmental variables need positive variation")
    standardized = (values - mean) / sd
    first, second = _edges(adjacency)
    distance = np.sqrt(np.mean((standardized[first] - standardized[second]) ** 2, axis=1))
    return cell_mean_edge_values(adjacency, distance), {
        "mean": mean.tolist(),
        "sd": sd.tolist(),
    }


def composition_boundary(composition: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Hellinger neighbour change for a land-cover composition estimate."""

    composition = np.asarray(composition, dtype=float)
    if composition.ndim != 2 or composition.shape[0] != adjacency.shape[0]:
        raise ValueError("composition must have one row per cell")
    if not np.isfinite(composition).all() or np.any(composition < 0):
        raise ValueError("land-cover composition must be finite and non-negative")
    totals = composition.sum(axis=1)
    if np.any(totals <= 0):
        raise ValueError("land-cover composition contains an empty cell")
    normalized = composition / totals[:, None]
    first, second = _edges(adjacency)
    distance = np.sqrt(
        0.5 * np.sum((np.sqrt(normalized[first]) - np.sqrt(normalized[second])) ** 2, axis=1)
    )
    return cell_mean_edge_values(adjacency, distance)


def average_ranks(values: np.ndarray) -> np.ndarray:
    """Return deterministic one-based average ranks, including exact ties."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not values.size or not np.isfinite(values).all():
        raise ValueError("rank input must be a non-empty finite vector")
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(values.size, dtype=float)
    start = 0
    while start < values.size:
        stop = start + 1
        while stop < values.size and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def spearman_rank_correlation(first: np.ndarray, second: np.ndarray) -> float:
    """Spearman correlation with deterministic average ranks for ties."""

    first_rank = average_ranks(first)
    second_rank = average_ranks(second)
    if first_rank.shape != second_rank.shape:
        raise ValueError("Spearman inputs must have the same length")
    first_sd = float(first_rank.std())
    second_sd = float(second_rank.std())
    if first_sd <= 0 or second_sd <= 0:
        raise ValueError("Spearman inputs need positive rank variation")
    return float(
        np.mean(
            ((first_rank - first_rank.mean()) / first_sd)
            * ((second_rank - second_rank.mean()) / second_sd)
        )
    )


def categorical_boundary(labels: Sequence[object], adjacency: np.ndarray) -> np.ndarray:
    """Mean fraction of incident neighbours assigned to another category."""

    labels = np.asarray([str(value) for value in labels], dtype=object)
    if labels.shape != (adjacency.shape[0],) or np.any(labels == ""):
        raise ValueError("categorical labels must be non-empty and match cells")
    first, second = _edges(adjacency)
    return cell_mean_edge_values(adjacency, (labels[first] != labels[second]).astype(float))


def environmental_boundary_surfaces(
    rows: Sequence[Mapping[str, Any]],
    adjacency: np.ndarray,
    *,
    families: Sequence[str] = PRIMARY_FAMILIES,
) -> dict[str, Any]:
    """Build four primary boundary surfaces without any flower-colour input."""

    forbidden = ("colour", "color", "flower", "transition")
    if not rows or any(
        any(token in str(key).casefold() for token in forbidden)
        for row in rows
        for key in row
    ):
        raise ValueError("environment table is empty or contains flower-colour outcomes")
    families = tuple(families)
    if len(set(families)) != len(families) or not set(families).issubset(PRIMARY_FAMILIES):
        raise ValueError("unknown or duplicate environmental family")
    output: dict[str, Any] = {"scaling": {}}
    if "macroclimate" in families:
        climate = np.asarray([[float(row[key]) for key in CLIMATE_FIELDS] for row in rows])
        surface, scaling = continuous_boundary(climate, adjacency)
        output["macroclimate"] = surface
        output["scaling"]["macroclimate"] = scaling
    if "terrain" in families:
        terrain = np.asarray([[float(row[key]) for key in TERRAIN_FIELDS] for row in rows])
        surface, scaling = continuous_boundary(terrain, adjacency)
        output["terrain"] = surface
        output["scaling"]["terrain"] = scaling
    if "land_cover" in families:
        land_cover = np.asarray(
            [[float(row[f"worldcover_{code}"]) for code in LAND_COVER_CODES] for row in rows]
        )
        output["land_cover"] = composition_boundary(land_cover, adjacency)
    if "ecoregion" in families:
        output["ecoregion"] = categorical_boundary(
            [row["ecoregion"] for row in rows], adjacency
        )
        output["realm_sensitivity"] = categorical_boundary(
            [row["realm"] for row in rows], adjacency
        )
        output["biome_sensitivity"] = categorical_boundary(
            [row["biome"] for row in rows], adjacency
        )
    return output


def weighted_cell_means(
    cell_ids: np.ndarray,
    weights: np.ndarray,
    values: np.ndarray,
    *,
    n_cells: int,
) -> np.ndarray:
    """Area-weight finite source pixels into equal-area grid cells."""

    cell_ids = np.asarray(cell_ids, dtype=int)
    weights = np.asarray(weights, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if (
        cell_ids.ndim != 1
        or weights.shape != cell_ids.shape
        or values.shape[0] != cell_ids.size
        or n_cells < 1
        or np.any((cell_ids < 0) | (cell_ids >= n_cells))
        or np.any(weights <= 0)
        or not np.isfinite(weights).all()
        or not np.isfinite(values).all()
    ):
        raise ValueError("invalid weighted source pixels")
    denominator = np.bincount(cell_ids, weights=weights, minlength=n_cells)
    result = np.full((n_cells, values.shape[1]), np.nan, dtype=float)
    present = denominator > 0
    for column in range(values.shape[1]):
        numerator = np.bincount(
            cell_ids, weights=weights * values[:, column], minlength=n_cells
        )
        result[present, column] = numerator[present] / denominator[present]
    return result


def weighted_dominant_labels(
    cell_ids: np.ndarray,
    weights: np.ndarray,
    labels: np.ndarray,
) -> dict[int, int]:
    """Choose maximum weighted label per cell, breaking exact ties numerically."""

    cell_ids = np.asarray(cell_ids, dtype=int)
    weights = np.asarray(weights, dtype=float)
    labels = np.asarray(labels, dtype=int)
    if (
        cell_ids.ndim != 1
        or weights.shape != cell_ids.shape
        or labels.shape != cell_ids.shape
        or np.any(cell_ids < 0)
        or np.any(labels <= 0)
        or np.any(weights <= 0)
        or not np.isfinite(weights).all()
    ):
        raise ValueError("invalid weighted categorical source pixels")
    order = np.lexsort((labels, cell_ids))
    ordered_cells = cell_ids[order]
    ordered_labels = labels[order]
    ordered_weights = weights[order]
    starts = np.r_[
        0,
        np.flatnonzero(
            (np.diff(ordered_cells) != 0) | (np.diff(ordered_labels) != 0)
        )
        + 1,
    ]
    pair_weight = np.add.reduceat(ordered_weights, starts)
    pair_cell = ordered_cells[starts]
    pair_label = ordered_labels[starts]
    result: dict[int, int] = {}
    best_weight: dict[int, float] = {}
    for cell, label, weight in zip(pair_cell, pair_label, pair_weight, strict=True):
        cell = int(cell)
        label = int(label)
        if cell not in result or weight > best_weight[cell]:
            result[cell] = label
            best_weight[cell] = float(weight)
    return result


def weighted_class_composition(
    cell_ids: np.ndarray,
    weights: np.ndarray,
    classes: np.ndarray,
    *,
    class_codes: Sequence[int],
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate a systematic categorical sample into cell compositions."""

    cell_ids = np.asarray(cell_ids, dtype=int)
    weights = np.asarray(weights, dtype=float)
    classes = np.asarray(classes, dtype=int)
    codes = tuple(int(code) for code in class_codes)
    if (
        not codes
        or len(set(codes)) != len(codes)
        or cell_ids.ndim != 1
        or weights.shape != cell_ids.shape
        or classes.shape != cell_ids.shape
        or np.any((cell_ids < 0) | (cell_ids >= n_cells))
        or np.any(weights <= 0)
        or not np.isfinite(weights).all()
        or not set(np.unique(classes)).issubset(codes)
    ):
        raise ValueError("invalid categorical samples")
    weighted = np.zeros((n_cells, len(codes)), dtype=float)
    counts = np.zeros(n_cells, dtype=int)
    lookup = {code: index for index, code in enumerate(codes)}
    np.add.at(counts, cell_ids, 1)
    for code, column in lookup.items():
        keep = classes == code
        np.add.at(weighted[:, column], cell_ids[keep], weights[keep])
    denominator = weighted.sum(axis=1)
    composition = np.full_like(weighted, np.nan)
    present = denominator > 0
    composition[present] = weighted[present] / denominator[present, None]
    return composition, counts
