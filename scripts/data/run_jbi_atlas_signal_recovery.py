#!/usr/bin/env python3
"""Qualify the exact atlas shared-boundary statistic on synthetic colour fields."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from fcp_pipeline.atlas_expansion import validate_expansion_contract
from fcp_pipeline.atlas_signal_recovery import (
    AtlasSpeciesGeometry,
    permutation_p_value,
    synthetic_colour_vectors,
)
from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    build_edge_cell_geometry,
    equal_area_cell_centers,
)
from fcp_pipeline.spatial_graph import spherical_knn_edges


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def numeric_seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--expansion-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json"),
    )
    parser.add_argument(
        "--atlas-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=Path("data/atlas/jbi_image_first_atlas_observation_manifest_v1.csv"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repetitions", type=int)
    parser.add_argument("--permutations", type=int)
    return parser.parse_args()


def build_geometries(
    rows: list[dict[str, str]], atlas_contract: dict[str, Any], scale_km: int
) -> tuple[list[AtlasSpeciesGeometry], EqualAreaGrid, np.ndarray]:
    geometry_contract = atlas_contract["geometry_only_scale_selection"]
    scale = next(
        row for row in geometry_contract["candidates"] if int(row["scale_km"]) == scale_km
    )
    grid = EqualAreaGrid(n_lon=int(scale["n_lon"]), n_sinlat=int(scale["n_sinlat"]))
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["species"], []).append(row)

    items: list[AtlasSpeciesGeometry] = []
    for species in sorted(grouped):
        species_rows = grouped[species]
        latitude = np.asarray([float(row["latitude"]) for row in species_rows])
        longitude = np.asarray([float(row["longitude"]) for row in species_rows])
        edges, distance = spherical_knn_edges(
            latitude,
            longitude,
            k=int(geometry_contract["knn_k"]),
        )
        geometry = build_edge_cell_geometry(
            latitude,
            longitude,
            edges,
            distance,
            grid=grid,
            max_edge_km=scale_km,
            min_edges_per_cell=int(geometry_contract["minimum_edges_per_species_cell"]),
        )
        items.append(
            AtlasSpeciesGeometry(
                species=species,
                latitude=latitude,
                longitude=longitude,
                geometry=geometry,
            )
        )
    detectable = np.vstack([item.geometry.detectable for item in items])
    return items, grid, detectable.sum(axis=0)


def main() -> None:
    args = parse_args()
    expansion = json.loads(args.expansion_contract.read_text(encoding="utf-8"))
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    validate_expansion_contract(expansion)
    recovery = expansion["estimator_qualification"]["signal_recovery"]
    repetitions = args.repetitions or int(recovery["simulation_repetitions"])
    permutations = args.permutations or int(recovery["permutations_per_repetition"])
    formal = (
        repetitions == int(recovery["simulation_repetitions"])
        and permutations == int(recovery["permutations_per_repetition"])
    )
    scale_km = int(expansion["spatial_design"]["primary_scale_km"])
    geometries, grid, opportunity = build_geometries(
        read_csv(args.observations), atlas, scale_km
    )
    cell_id, center_lat, center_lon = equal_area_cell_centers(grid)
    anchor_id = int(cell_id[np.flatnonzero(opportunity == opportunity.max())[0]])
    anchor_lat = float(center_lat[anchor_id])
    anchor_lon = float(center_lon[anchor_id])
    min_species = int(
        atlas["geometry_only_scale_selection"]["passing_criteria"][
            "minimum_detectable_species_per_shared_cell"
        ]
    )

    jobs: list[tuple[str, float]] = [("null_stationary", 0.0)]
    jobs.append(("within_species_heterogeneous_boundaries", 2.0))
    jobs.extend(
        ("shared_geographic_boundary", float(effect))
        for effect in recovery["effect_sizes"]
        if float(effect) > 0
    )
    rows: list[dict[str, Any]] = []
    for scenario, effect in jobs:
        for repetition in range(repetitions):
            generator = np.random.default_rng(
                numeric_seed("signal", scenario, effect, repetition)
            )
            permuter = np.random.default_rng(
                numeric_seed("permutation", scenario, effect, repetition)
            )
            vectors = synthetic_colour_vectors(
                geometries,
                effect_size=effect,
                scenario=scenario,
                shared_anchor_latitude=anchor_lat,
                shared_anchor_longitude=anchor_lon,
                rng=generator,
            )
            statistic, p_value = permutation_p_value(
                vectors,
                geometries,
                min_detectable_species=min_species,
                permutations=permutations,
                rng=permuter,
            )
            rows.append(
                {
                    "scenario": scenario,
                    "effect_size": effect,
                    "repetition": repetition + 1,
                    "statistic": statistic,
                    "p_value": p_value,
                    "detected": p_value <= float(recovery["alpha"]),
                }
            )
            print(
                f"scenario={scenario} effect={effect:g} repetition={repetition + 1}/{repetitions}",
                flush=True,
            )

    rates: dict[str, float] = {}
    for scenario, effect in jobs:
        subset = [
            row
            for row in rows
            if row["scenario"] == scenario and float(row["effect_size"]) == effect
        ]
        rates[f"{scenario}|{effect:g}"] = sum(bool(row["detected"]) for row in subset) / len(
            subset
        )
    shared_rates = [
        rates[f"shared_geographic_boundary|{effect:g}"]
        for effect in (0.5, 1.0, 2.0)
    ]
    checks = {
        "maximum_null_false_positive_rate": rates["null_stationary|0"]
        <= float(recovery["maximum_null_false_positive_rate"]),
        "maximum_heterogeneous_false_shared_rate": rates[
            "within_species_heterogeneous_boundaries|2"
        ]
        <= float(recovery["maximum_heterogeneous_false_shared_rate"]),
        "minimum_shared_boundary_power_at_effect_2": rates[
            "shared_geographic_boundary|2"
        ]
        >= float(recovery["minimum_shared_boundary_power_at_effect_2"]),
        "monotone_power_required": shared_rates == sorted(shared_rates),
    }
    passed = formal and all(checks.values())
    result = {
        "status": (
            "pass_exact_geometry_signal_recovery"
            if passed
            else "smoke_only_not_qualification"
            if not formal
            else "stop_signal_recovery_failed"
        ),
        "formal_contract_counts": formal,
        "species": len(geometries),
        "scale_km": scale_km,
        "anchor_cell_id": anchor_id,
        "anchor_latitude": anchor_lat,
        "anchor_longitude": anchor_lon,
        "repetitions": repetitions,
        "permutations_per_repetition": permutations,
        "detection_rates": rates,
        "checks": checks,
        "atlas_pixels_permitted_by_signal_gate": passed,
    }
    write_csv(args.output_dir / "signal_recovery_repetitions.csv", rows)
    write_json(args.output_dir / "signal_recovery_result.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if formal and not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
