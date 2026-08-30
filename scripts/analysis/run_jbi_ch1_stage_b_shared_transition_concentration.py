#!/usr/bin/env python3
"""Run frozen Stage-B shared continuous transition concentration analysis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.continuous_colour_boundaries import (
    average_rank_intensity,
    edge_colour_discontinuity,
    opportunity_weighted_concentration,
    shared_boundary_intensity,
)
from fcp_pipeline.shared_transition_surface import (
    EdgeCellGeometry,
    EqualAreaGrid,
    build_edge_cell_geometry,
    cell_mean_intensity,
    equal_area_cell_centers,
    geometry_opportunity_summary,
)
from fcp_pipeline.spatial_graph import spherical_knn_edges

PROTOCOL = "jbi-ch1-stage-b-shared-transition-concentration-v1"
EXPECTED_SPECIES = [
    "Antirrhinum majus",
    "Dactylorhiza sambucina",
    "Gentiana lutea",
    "Ipomoea purpurea",
    "Lysimachia arvensis",
    "Raphanus sativus",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def species_slug(species: str) -> str:
    return species.lower().replace(" ", "_")


def configuration_key(max_edge_km: float, grid: EqualAreaGrid) -> str:
    cap = format(float(max_edge_km), ".12g").replace(".", "p")
    return f"cap_{cap}km_grid_{grid.n_lon}x{grid.n_sinlat}"


def mc_upper_summary(observed: float, null: np.ndarray) -> dict[str, object]:
    null = np.asarray(null, dtype=float)
    if null.ndim != 1 or len(null) < 1 or not np.isfinite(null).all():
        raise ValueError("null must be a non-empty finite vector")
    mean = float(null.mean())
    sd = float(null.std(ddof=1)) if len(null) > 1 else 0.0
    upper_p = float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))
    observed_dev = abs(observed - mean)
    two_p = float((1 + np.count_nonzero(np.abs(null - mean) >= observed_dev)) / (len(null) + 1))
    q = np.quantile(null, [0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995])
    return {
        "observed": float(observed),
        "null_mean": mean,
        "null_sd": sd,
        "concentration_excess": float(observed - mean),
        "standardized_concentration_excess": float((observed - mean) / sd) if sd > 0 else None,
        "p_upper_tail": upper_p,
        "p_two_sided_descriptive": two_p,
        "null_quantiles": {
            "p005": float(q[0]),
            "p025": float(q[1]),
            "p05": float(q[2]),
            "p50": float(q[3]),
            "p95": float(q[4]),
            "p975": float(q[5]),
            "p995": float(q[6]),
        },
    }


def build_species_inputs(
    rows: list[dict],
    coord_by_photo: dict[str, tuple[str, float, float]],
    representation: dict,
    *,
    k: int,
    minimum_rows: int,
) -> dict[str, dict[str, np.ndarray]]:
    """Validate evaluation rows and build the colour-blind base graph per species."""

    out: dict[str, dict[str, np.ndarray]] = {}
    for species in EXPECTED_SPECIES:
        group = [
            row
            for row in rows
            if row.get("species") == species
            and row.get("feature_status") == "ok"
            and row.get("continuous_colour_vector_z") is not None
        ]
        group = sorted(group, key=lambda row: str(row["photo_id"]))
        if len(group) < minimum_rows:
            raise ValueError(f"{species}: only {len(group)} evaluable rows; minimum is {minimum_rows}")

        expected_dimension = len(representation["per_species"][species]["feature_names"])
        values = np.asarray([row["continuous_colour_vector_z"] for row in group], dtype=float)
        if values.shape != (len(group), expected_dimension) or not np.isfinite(values).all():
            raise ValueError(f"{species}: invalid continuous colour vectors")

        photo_id = np.asarray([str(row["photo_id"]) for row in group], dtype=object)
        coordinates = []
        for pid in photo_id:
            if pid not in coord_by_photo:
                raise ValueError(f"missing frozen coordinate for photo {pid}")
            coord_species, latitude, longitude = coord_by_photo[pid]
            if coord_species != species:
                raise ValueError(f"species mismatch for photo {pid}")
            coordinates.append((latitude, longitude))
        coordinates_array = np.asarray(coordinates, dtype=float)
        edges, edge_distance = spherical_knn_edges(
            coordinates_array[:, 0],
            coordinates_array[:, 1],
            k=k,
            max_edge_km=None,
        )
        out[species] = {
            "photo_id": photo_id,
            "values": values,
            "latitude": coordinates_array[:, 0],
            "longitude": coordinates_array[:, 1],
            "base_edges": edges,
            "base_edge_distance_km": edge_distance,
        }
    return out


def freeze_geometry_candidates(
    species_inputs: dict[str, dict[str, np.ndarray]],
    contract: dict,
) -> tuple[dict[str, object], dict[str, tuple[EqualAreaGrid, list[EdgeCellGeometry]]], str | None]:
    """Evaluate and select geometry configurations without reading any colour values."""

    selection = contract["geometry_only_primary_selection"]
    min_edges_per_cell = int(selection["species_cell_detectable_rule"].split("at least ")[1].split()[0])
    min_detectable_species = int(contract["shared_surface"]["minimum_detectable_species"])
    criteria = selection["passing_criteria"]

    candidate_audit: list[dict[str, object]] = []
    geometry_cache: dict[str, tuple[EqualAreaGrid, list[EdgeCellGeometry]]] = {}
    selected_key: str | None = None
    rank = 0

    for max_edge_km_raw in selection["candidate_max_edge_km_in_priority_order"]:
        max_edge_km = float(max_edge_km_raw)
        for grid_spec in selection["candidate_equal_area_grids_in_priority_order_within_edge_cap"]:
            rank += 1
            grid = EqualAreaGrid(n_lon=int(grid_spec["n_lon"]), n_sinlat=int(grid_spec["n_sinlat"]))
            key = configuration_key(max_edge_km, grid)
            geometries: list[EdgeCellGeometry] = []
            failure: str | None = None
            for species in EXPECTED_SPECIES:
                data = species_inputs[species]
                try:
                    geometry = build_edge_cell_geometry(
                        data["latitude"],
                        data["longitude"],
                        data["base_edges"],
                        data["base_edge_distance_km"],
                        grid=grid,
                        max_edge_km=max_edge_km,
                        min_edges_per_cell=min_edges_per_cell,
                    )
                except ValueError as exc:
                    failure = f"{species}: {exc}"
                    break
                geometries.append(geometry)

            if failure is not None:
                candidate_audit.append(
                    {
                        "selection_rank": rank,
                        "configuration": key,
                        "max_edge_km": max_edge_km,
                        "grid": {
                            "n_lon": grid.n_lon,
                            "n_sinlat": grid.n_sinlat,
                            "n_cells": grid.n_cells,
                            "cell_area_km2": grid.cell_area_km2,
                        },
                        "passes_geometry_only_criteria": False,
                        "failure": failure,
                    }
                )
                continue

            summary = geometry_opportunity_summary(
                geometries,
                min_detectable_species=min_detectable_species,
            )
            retained = {
                species: int(value)
                for species, value in zip(EXPECTED_SPECIES, summary.pop("retained_edges_per_species"), strict=True)
            }
            detectable_cells = {
                species: int(value)
                for species, value in zip(EXPECTED_SPECIES, summary.pop("detectable_cells_per_species"), strict=True)
            }
            checks = {
                "minimum_retained_edges_per_species": min(retained.values())
                >= int(criteria["minimum_retained_edges_per_species"]),
                "minimum_cells_A_ge_2": int(summary["n_cells_A_ge_2"])
                >= int(criteria["minimum_cells_A_ge_2"]),
                "minimum_cells_A_ge_3": int(summary["n_cells_A_ge_3"])
                >= int(criteria["minimum_cells_A_ge_3"]),
                "minimum_species_with_any_A_ge_2_opportunity": int(summary["species_with_any_shared_opportunity"])
                >= int(criteria["minimum_species_with_any_A_ge_2_opportunity"]),
            }
            passed = bool(all(checks.values()))
            candidate = {
                "selection_rank": rank,
                "configuration": key,
                "max_edge_km": max_edge_km,
                "grid": {
                    "n_lon": grid.n_lon,
                    "n_sinlat": grid.n_sinlat,
                    "n_cells": grid.n_cells,
                    "cell_area_km2": grid.cell_area_km2,
                },
                "min_edges_per_species_cell": min_edges_per_cell,
                "minimum_detectable_species_for_shared_cell": min_detectable_species,
                "geometry_summary": summary,
                "retained_edges_by_species": retained,
                "detectable_cells_by_species": detectable_cells,
                "criteria_checks": checks,
                "passes_geometry_only_criteria": passed,
                "failure": None,
            }
            candidate_audit.append(candidate)
            geometry_cache[key] = (grid, geometries)
            if passed and selected_key is None:
                selected_key = key

    audit: dict[str, object] = {
        "protocol": PROTOCOL,
        "status": "geometry_selection_complete_before_observed_colour_scoring",
        "selection_used_colour_values": False,
        "selection_rule": selection["priority_rule"],
        "passing_criteria": criteria,
        "selected_primary_configuration": selected_key,
        "candidate_configurations": candidate_audit,
    }
    return audit, geometry_cache, selected_key


def compute_surface(
    species_inputs: dict[str, dict[str, np.ndarray]],
    geometries: list[EdgeCellGeometry],
    *,
    permuted_values: list[np.ndarray] | None,
    min_detectable_species: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Recompute the complete transition-intensity surface for one label assignment."""

    species_cell_intensity = []
    detectable = []
    for index, species in enumerate(EXPECTED_SPECIES):
        data = species_inputs[species]
        geometry = geometries[index]
        values = data["values"] if permuted_values is None else permuted_values[index]
        scores = edge_colour_discontinuity(values, geometry.retained_edges)
        edge_intensity = average_rank_intensity(scores)
        species_cell_intensity.append(cell_mean_intensity(edge_intensity, geometry))
        detectable.append(geometry.detectable)

    intensity_matrix = np.vstack(species_cell_intensity)
    detectable_matrix = np.vstack(detectable)
    shared, opportunity = shared_boundary_intensity(
        intensity_matrix,
        detectable_matrix,
        min_detectable_species=min_detectable_species,
    )
    return intensity_matrix, detectable_matrix, shared, opportunity


def run_configuration(
    species_inputs: dict[str, dict[str, np.ndarray]],
    grid: EqualAreaGrid,
    geometries: list[EdgeCellGeometry],
    *,
    min_detectable_species: int,
    n_permutations: int,
    seed: int,
) -> tuple[dict[str, object], np.ndarray, dict[str, np.ndarray]]:
    """Run the observed and complete within-species permutation pipeline."""

    intensity, detectable, shared, opportunity = compute_surface(
        species_inputs,
        geometries,
        permuted_values=None,
        min_detectable_species=min_detectable_species,
    )
    observed = opportunity_weighted_concentration(
        shared,
        opportunity,
        min_opportunity=min_detectable_species,
    )

    rng_by_species = [np.random.default_rng(seed + 1009 * (index + 1)) for index in range(len(EXPECTED_SPECIES))]
    null = np.empty(n_permutations, dtype=float)
    for permutation in range(n_permutations):
        permuted = [
            species_inputs[species]["values"][rng_by_species[index].permutation(len(species_inputs[species]["values"]))]
            for index, species in enumerate(EXPECTED_SPECIES)
        ]
        _, _, shared_b, opportunity_b = compute_surface(
            species_inputs,
            geometries,
            permuted_values=permuted,
            min_detectable_species=min_detectable_species,
        )
        if not np.array_equal(opportunity_b, opportunity):
            raise RuntimeError("label-independent opportunity changed after colour permutation")
        null[permutation] = opportunity_weighted_concentration(
            shared_b,
            opportunity_b,
            min_opportunity=min_detectable_species,
        )

    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    evaluable = (opportunity >= min_detectable_species) & np.isfinite(shared)
    top_order = np.flatnonzero(evaluable)
    top_order = top_order[np.lexsort((-opportunity[top_order], -shared[top_order]))]
    top_cells = [
        {
            "cell_id": int(cell_id[cell]),
            "latitude": float(latitude[cell]),
            "longitude": float(longitude[cell]),
            "A": int(opportunity[cell]),
            "shared_transition_intensity": float(shared[cell]),
        }
        for cell in top_order[:20]
    ]
    species_summary = {
        species: {
            "retained_edges": int(len(geometries[index].retained_edges)),
            "detectable_cells": int(np.count_nonzero(detectable[index])),
            "mean_intensity_over_detectable_cells": float(np.nanmean(intensity[index, detectable[index]])),
            "maximum_cell_intensity": float(np.nanmax(intensity[index, detectable[index]])),
        }
        for index, species in enumerate(EXPECTED_SPECIES)
    }
    result: dict[str, object] = {
        "n_permutations": int(n_permutations),
        "random_seed": int(seed),
        "global_concentration": mc_upper_summary(observed, null),
        "surface": {
            "n_cells_total": int(grid.n_cells),
            "n_cells_evaluable_A_ge_minimum": int(np.count_nonzero(evaluable)),
            "n_cells_A_ge_2": int(np.count_nonzero(opportunity >= 2)),
            "n_cells_A_ge_3": int(np.count_nonzero(opportunity >= 3)),
            "n_cells_A_ge_4": int(np.count_nonzero(opportunity >= 4)),
            "maximum_A": int(opportunity.max(initial=0)),
            "shared_intensity_mean": float(np.nanmean(shared[evaluable])),
            "shared_intensity_sd": float(np.nanstd(shared[evaluable], ddof=1)) if np.count_nonzero(evaluable) > 1 else 0.0,
            "top_observed_cells": top_cells,
        },
        "species": species_summary,
    }
    arrays = {
        "cell_id": cell_id,
        "latitude": latitude,
        "longitude": longitude,
        "species_intensity": intensity,
        "detectable": detectable,
        "shared": shared,
        "opportunity": opportunity,
    }
    return result, null, arrays


def write_surface_csv(
    path: Path,
    arrays: dict[str, np.ndarray],
    geometries: list[EdgeCellGeometry],
    *,
    min_detectable_species: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "cell_id",
        "latitude",
        "longitude",
        "opportunity_A",
        "evaluable_A_ge_minimum",
        "shared_transition_intensity",
    ]
    for species in EXPECTED_SPECIES:
        slug = species_slug(species)
        header.extend([f"{slug}_detectable", f"{slug}_edge_count", f"{slug}_transition_intensity"])

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for cell in range(len(arrays["cell_id"])):
            shared_value = arrays["shared"][cell]
            row: list[object] = [
                int(arrays["cell_id"][cell]),
                format(float(arrays["latitude"][cell]), ".12g"),
                format(float(arrays["longitude"][cell]), ".12g"),
                int(arrays["opportunity"][cell]),
                bool(arrays["opportunity"][cell] >= min_detectable_species),
                "" if not np.isfinite(shared_value) else format(float(shared_value), ".12g"),
            ]
            for index, _species in enumerate(EXPECTED_SPECIES):
                value = arrays["species_intensity"][index, cell]
                row.extend(
                    [
                        bool(arrays["detectable"][index, cell]),
                        int(geometries[index].cell_edge_count[cell]),
                        "" if not np.isfinite(value) else format(float(value), ".12g"),
                    ]
                )
            writer.writerow(row)


def write_null_csv(path: Path, null: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["permutation", "opportunity_weighted_shared_transition_concentration"])
        for index, value in enumerate(null, start=1):
            writer.writerow([index, format(float(value), ".12g")])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/evaluation/jbi_ch1_florence_evaluation_features_v1.jsonl"))
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument("--representation", type=Path, default=Path("docs/supporting/jbi_ch1_continuous_colour_representation_v1.json"))
    parser.add_argument("--stage-a", type=Path, default=Path("docs/supporting/jbi_ch1_stage_a_continuous_graph_v1.json"))
    parser.add_argument("--contract", type=Path, default=Path("docs/supporting/jbi_ch1_stage_b_shared_transition_contract_v1.json"))
    parser.add_argument("--geometry-audit", type=Path, default=Path("docs/supporting/jbi_ch1_stage_b_geometry_audit_v1.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/supporting/jbi_ch1_stage_b_shared_transition_concentration_v1.json"))
    parser.add_argument("--surface-csv", type=Path, default=Path("data/evaluation/jbi_ch1_stage_b_shared_transition_surface_v1.csv"))
    parser.add_argument("--primary-null-csv", type=Path, default=Path("data/evaluation/jbi_ch1_stage_b_primary_null_v1.csv"))
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    representation = json.loads(args.representation.read_text(encoding="utf-8"))
    stage_a = json.loads(args.stage_a.read_text(encoding="utf-8"))
    if contract.get("protocol") != PROTOCOL or contract.get("status") != "frozen_after_stage_a_before_any_shared_boundary_result":
        raise ValueError("Stage B contract is not frozen")
    gate = contract["gate"]
    if stage_a.get("protocol") != gate["required_stage_a_protocol"] or stage_a.get("status") != gate["required_stage_a_status"]:
        raise ValueError("Stage A gate is not complete")
    if bool(stage_a.get("primary_rejects_random_labelling_at_0_05")) is not bool(gate["required_stage_a_rejects_random_labelling_at_0_05"]):
        raise ValueError("Stage A gate decision does not match the frozen Stage B requirement")
    if representation.get("status") != "frozen_before_evaluation_values_inspected":
        raise ValueError("continuous representation is not frozen")

    rows = load_jsonl(args.features)
    if len(rows) != int(gate["required_evaluation_records"]) or len({str(row["photo_id"]) for row in rows}) != int(gate["required_evaluation_records"]):
        raise ValueError("expected exactly 720 unique evaluation feature records")
    if any(row.get("evaluation_row") is not True or row.get("final_label") is not False for row in rows):
        raise ValueError("evaluation feature contract violation")
    if sorted({str(row["species"]) for row in rows}) != EXPECTED_SPECIES:
        raise ValueError("unexpected evaluation species set")
    if stage_a.get("evaluation_feature_sha256") != sha256(args.features):
        raise ValueError("Stage B feature input differs from the completed Stage A input")

    split = pd.read_csv(args.split)
    evaluation_split = split.loc[split["split"].astype(str).eq("evaluation")].copy()
    if len(evaluation_split) != int(gate["required_evaluation_records"]):
        raise ValueError("frozen split does not contain 720 evaluation rows")
    coord_by_photo = {
        str(row.photo_id): (str(row.species), float(row.latitude), float(row.longitude))
        for row in evaluation_split.itertuples(index=False)
    }
    if set(coord_by_photo) != {str(row["photo_id"]) for row in rows}:
        raise ValueError("evaluation feature IDs do not exactly match frozen coordinates")

    minimum_rows = int(contract["input"]["minimum_evaluable_rows_per_species"])
    base_k = int(contract["species_transition_intensity"]["base_graph_k"])
    species_inputs = build_species_inputs(
        rows,
        coord_by_photo,
        representation,
        k=base_k,
        minimum_rows=minimum_rows,
    )

    # This entire selection block uses only species, coordinates, graph edges and edge
    # distances.  No observed colour score is computed before selected_key is fixed.
    geometry_audit, geometry_cache, selected_key = freeze_geometry_candidates(species_inputs, contract)
    geometry_audit.update(
        {
            "evaluation_feature_sha256": sha256(args.features),
            "frozen_split_sha256": sha256(args.split),
            "frozen_contract_sha256": sha256(args.contract),
            "base_graph_k": base_k,
        }
    )
    args.geometry_audit.parent.mkdir(parents=True, exist_ok=True)
    args.geometry_audit.write_text(
        json.dumps(geometry_audit, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    common = {
        "protocol": PROTOCOL,
        "primary_representation": "species-specific continuous colour vector standardized from frozen calibration parameters",
        "n_evaluation_records": len(rows),
        "species": EXPECTED_SPECIES,
        "base_graph_k": base_k,
        "geometry_selection_used_colour_values": False,
        "geometry_audit_sha256": sha256(args.geometry_audit),
        "evaluation_feature_sha256": sha256(args.features),
        "frozen_split_sha256": sha256(args.split),
        "frozen_representation_sha256": sha256(args.representation),
        "completed_stage_a_sha256": sha256(args.stage_a),
        "frozen_contract_sha256": sha256(args.contract),
        "environment_used": False,
        "geographic_reference_library_used": False,
        "interpretation_limit": "Stage B tests cross-species geographic concentration of relative continuous transition intensity; it does not identify environmental, historical, or mechanistic causes.",
    }
    if selected_key is None:
        result = {
            **common,
            "status": "stage_b_not_estimable_under_frozen_geometry_support",
            "selected_primary_configuration": None,
            "primary_result": None,
            "primary_rejects_shared_concentration_null_at_0_05": None,
            "sensitivity_results": {},
            "next_gate": "retain the supported Stage A local-organization result and report that shared-boundary concentration was not estimable under the frozen geometry-support criteria",
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0

    min_detectable_species = int(contract["shared_surface"]["minimum_detectable_species"])
    primary_permutations = int(contract["null"]["primary_permutations"])
    sensitivity_permutations = int(contract["null"]["sensitivity_permutations"])
    seed = int(contract["null"]["random_seed"])
    primary_grid, primary_geometries = geometry_cache[selected_key]
    primary_result, primary_null, primary_arrays = run_configuration(
        species_inputs,
        primary_grid,
        primary_geometries,
        min_detectable_species=min_detectable_species,
        n_permutations=primary_permutations,
        seed=seed,
    )
    write_surface_csv(
        args.surface_csv,
        primary_arrays,
        primary_geometries,
        min_detectable_species=min_detectable_species,
    )
    write_null_csv(args.primary_null_csv, primary_null)

    passing_keys = [
        candidate["configuration"]
        for candidate in geometry_audit["candidate_configurations"]
        if candidate.get("passes_geometry_only_criteria") is True
    ]
    sensitivity_results: dict[str, object] = {}
    for key in passing_keys:
        if key == selected_key:
            continue
        grid, geometries = geometry_cache[key]
        result_s, _null_s, _arrays_s = run_configuration(
            species_inputs,
            grid,
            geometries,
            min_detectable_species=min_detectable_species,
            n_permutations=sensitivity_permutations,
            seed=seed,
        )
        sensitivity_results[key] = result_s

    selected_audit = next(
        candidate
        for candidate in geometry_audit["candidate_configurations"]
        if candidate["configuration"] == selected_key
    )
    primary_summary = primary_result["global_concentration"]
    result = {
        **common,
        "status": "stage_b_evaluation_complete",
        "selected_primary_configuration": selected_audit,
        "primary_result": primary_result,
        "primary_permutations": primary_permutations,
        "primary_rejects_shared_concentration_null_at_0_05": bool(primary_summary["p_upper_tail"] <= 0.05),
        "sensitivity_permutations_per_configuration": sensitivity_permutations,
        "sensitivity_results": sensitivity_results,
        "surface_csv_sha256": sha256(args.surface_csv),
        "primary_null_csv_sha256": sha256(args.primary_null_csv),
        "next_gate": (
            "freeze a geographic reference library and test post-discovery correspondence without changing the discovered transition surface"
            if primary_summary["p_upper_tail"] <= 0.05
            else "retain Stage A local organization and report no confirmatory evidence for shared-boundary concentration at the frozen Stage B supports"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
