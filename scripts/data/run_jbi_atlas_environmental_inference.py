#!/usr/bin/env python3
"""Join the sealed atlas only after completeness, then run frozen inference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_colour_inference import (
    SEASON_CONFIGURATIONS,
    build_species_transition_surface,
    equal_species_cohort_surface,
    joint_equal_cohort_spectral_test,
    prepare_spectral_cohort_test,
    validate_colour_inference_contract,
)
from fcp_pipeline.atlas_measurement import (
    validate_inference_contract,
    validate_measurement_result_rows,
)
from fcp_pipeline.flower_roi_v4_runtime import (
    file_sha256,
    validate_scaleout_authorization,
)
from scripts.data.validate_jbi_atlas_roi_v4_gate_evidence import (
    load_committed_locked_scaleout_result,
)
from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    equal_area_cell_centers,
)


ROLE_FIELDS = {
    "flower": ("flower_L_mean", "flower_a_mean", "flower_b_mean"),
    "background": (
        "background_L_mean",
        "background_a_mean",
        "background_b_mean",
    ),
}
OVERLAY_FIELDS = (
    "macroclimate_boundary",
    "land_cover_boundary",
    "ecoregion_boundary",
    "realm_sensitivity_boundary",
    "biome_sensitivity_boundary",
)
PRIMARY_OVERLAYS = (
    "macroclimate_boundary",
    "land_cover_boundary",
    "ecoregion_boundary",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    union = {key for row in rows for key in row}
    fields.extend(sorted(union - set(fields)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    folded = str(value).strip().casefold()
    if folded == "true":
        return True
    if folded == "false":
        return False
    raise ValueError(f"expected a boolean, got {value!r}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-results-dir", type=Path, required=True)
    parser.add_argument("--measurement-gate", type=Path, required=True)
    parser.add_argument("--sealed-coordinate-key", type=Path, required=True)
    parser.add_argument("--environmental-coverage-result", type=Path, required=True)
    parser.add_argument("--environment-dir", type=Path, required=True)
    parser.add_argument("--roi-evidence-dir", type=Path, required=True)
    parser.add_argument("--trained-weight", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_colour_surface_contract_v1.json",
    )
    parser.add_argument(
        "--inference-contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v3.json",
    )
    return parser.parse_args()


def load_location_free_results(directory: Path) -> list[dict[str, Any]]:
    paths = sorted(directory.glob("measurement_shard_*.csv"))
    if not paths:
        raise RuntimeError("no location-free measurement shards were found")
    rows: list[dict[str, Any]] = []
    for path in paths:
        for row in read_csv(path):
            row["background_features_available"] = parse_bool(
                row["background_features_available"]
            )
            rows.append(row)
    return validate_measurement_result_rows(rows)


def finite_lab(row: Mapping[str, Any], role: str) -> np.ndarray | None:
    try:
        values = np.asarray([float(row[field]) for field in ROLE_FIELDS[role]], dtype=float)
    except (KeyError, TypeError, ValueError):
        return None
    return values if np.isfinite(values).all() else None


def season_labels(rows: Sequence[Mapping[str, Any]], configuration: str) -> np.ndarray | None:
    if configuration == "all_dates":
        return None
    if configuration == "same_calendar_month_edges":
        return np.asarray([int(row["observed_month"]) for row in rows], dtype=int)
    if configuration == "same_local_solar_quarter_edges":
        return np.asarray([int(row["local_solar_quarter"]) for row in rows], dtype=int)
    raise ValueError(f"unknown season configuration: {configuration}")


def environment_vectors(
    path: Path, *, grid: EqualAreaGrid
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    rows = read_csv(path)
    cells, latitude, longitude = equal_area_cell_centers(grid)
    overlays = {
        field: np.full(grid.n_cells, np.nan, dtype=float) for field in OVERLAY_FIELDS
    }
    seen: set[int] = set()
    for row in rows:
        cell = int(row["cell_id"])
        if cell in seen or cell < 0 or cell >= grid.n_cells:
            raise ValueError(f"invalid duplicate environmental cell: {cell}")
        seen.add(cell)
        for field in OVERLAY_FIELDS:
            value = row.get(field, "")
            if value not in (None, ""):
                parsed = float(value)
                if math.isfinite(parsed):
                    overlays[field][cell] = parsed
    return cells, latitude, longitude, overlays


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_colour_inference_contract(contract)
    inference_contract = json.loads(
        args.inference_contract.read_text(encoding="utf-8")
    )
    validate_inference_contract(inference_contract)

    # All authorization and completeness evidence is checked before the coordinate
    # key is read.  Failure above this line leaves the sealed join unopened.
    trained_weight_sha = file_sha256(args.trained_weight)
    roi_result = load_committed_locked_scaleout_result(args.roi_evidence_dir)
    validate_scaleout_authorization(
        roi_result, trained_weight_sha256=trained_weight_sha
    )
    measurement_gate = json.loads(
        args.measurement_gate.read_text(encoding="utf-8")
    )
    if (
        measurement_gate.get("status") != "pass_scaleout_measurement_completeness"
        or measurement_gate.get("coordinate_join_permitted") is not True
        or measurement_gate.get("coordinates_opened") is not False
        or measurement_gate.get("frozen_measurements") != 60_000
    ):
        raise RuntimeError("measurement completeness does not authorize coordinates")
    coverage = json.loads(
        args.environmental_coverage_result.read_text(encoding="utf-8")
    )
    if (
        coverage.get("status") != "pass_precolour_environmental_coverage"
        or coverage.get("scaleout_colour_opened") is not False
        or "macroclimate" not in coverage.get("evaluable_families", [])
        or len(coverage.get("evaluable_families", [])) < 2
    ):
        raise RuntimeError("pre-colour environmental coverage did not pass")
    measurements = load_location_free_results(args.measurement_results_dir)
    measurement_by_id = {str(row["measurement_id"]): row for row in measurements}
    if len(measurements) != 60_000 or len(measurement_by_id) != 60_000:
        raise RuntimeError("location-free measurement denominator changed")

    # First and only protected join in this runner.
    coordinates = read_csv(args.sealed_coordinate_key)
    coordinate_by_id = {str(row["measurement_id"]): row for row in coordinates}
    if (
        len(coordinates) != 60_000
        or len(coordinate_by_id) != 60_000
        or set(coordinate_by_id) != set(measurement_by_id)
    ):
        raise RuntimeError("sealed coordinate key does not match measurements")
    joined: list[dict[str, Any]] = []
    for measurement_id in sorted(measurement_by_id):
        measured = measurement_by_id[measurement_id]
        coordinate = coordinate_by_id[measurement_id]
        if measured["species_blind_id"] != coordinate["species_blind_id"]:
            raise RuntimeError("species-blind identity changed at protected join")
        joined.append({**coordinate, **measured})

    expected_cohorts = [f"C{index:02d}" for index in range(1, 9)]
    cohorts = sorted({str(row["cohort_id"]) for row in joined})
    if cohorts != expected_cohorts:
        raise RuntimeError("atlas cohort denominator changed")
    species_to_cohort: dict[str, str] = {}
    for row in joined:
        species_id = str(row["species_blind_id"])
        prior = species_to_cohort.setdefault(species_id, str(row["cohort_id"]))
        if prior != row["cohort_id"]:
            raise RuntimeError("one species crossed frozen cohorts")
    if len(species_to_cohort) != 200:
        raise RuntimeError("atlas species denominator changed")

    transition = contract["transition_surface"]
    groups: dict[str, list[Any]] = {}
    species_status: list[dict[str, Any]] = []
    cohort_status: list[dict[str, Any]] = []
    environment_cache: dict[int, tuple[Any, ...]] = {}
    rows_by_species = {
        species_id: [row for row in joined if row["species_blind_id"] == species_id]
        for species_id in sorted(species_to_cohort)
    }
    for scale in transition["scales_km"]:
        n_lon, n_sinlat = transition["grids"][str(scale)]
        grid = EqualAreaGrid(int(n_lon), int(n_sinlat))
        environment_cache[scale] = environment_vectors(
            args.environment_dir / f"environmental_boundary_cells_{scale}km.csv",
            grid=grid,
        )
        cells, cell_latitude, cell_longitude, overlays = environment_cache[scale]
        for configuration in SEASON_CONFIGURATIONS:
            for role in ("flower", "background"):
                cohort_surfaces: dict[str, tuple[np.ndarray, np.ndarray]] = {}
                for cohort in cohorts:
                    surfaces = []
                    cohort_species = sorted(
                        species
                        for species, assigned in species_to_cohort.items()
                        if assigned == cohort
                    )
                    for species_id in cohort_species:
                        raw_rows = rows_by_species[species_id]
                        analysis_rows = []
                        lab_rows = []
                        for row in raw_rows:
                            if row["automated_colour_state_status"] != (
                                "automated_colour_state_admitted"
                            ):
                                continue
                            if role == "background" and row[
                                "background_features_available"
                            ] is not True:
                                continue
                            lab = finite_lab(row, role)
                            if lab is None:
                                continue
                            analysis_rows.append(row)
                            lab_rows.append(lab)
                        if analysis_rows:
                            surface = build_species_transition_surface(
                                [float(row["latitude"]) for row in analysis_rows],
                                [float(row["longitude"]) for row in analysis_rows],
                                np.vstack(lab_rows),
                                grid=grid,
                                scale_km=int(scale),
                                season_labels=season_labels(
                                    analysis_rows, configuration
                                ),
                                knn_k=int(transition["knn_k"]),
                                minimum_edges_per_cell=int(
                                    transition["minimum_edges_per_species_cell"]
                                ),
                                minimum_retained_edges=int(
                                    transition[
                                        "minimum_retained_edges_per_species_configuration"
                                    ]
                                ),
                                minimum_detectable_cells=int(
                                    transition[
                                        "minimum_detectable_cells_per_species_configuration"
                                    ]
                                ),
                            )
                        else:
                            surface = build_species_transition_surface(
                                [],
                                [],
                                np.empty((0, 3)),
                                grid=grid,
                                scale_km=int(scale),
                            )
                        surfaces.append(surface)
                        species_status.append(
                            {
                                "species_blind_id": species_id,
                                "species": raw_rows[0]["species"],
                                "cohort_id": cohort,
                                "surface_role": role,
                                "scale_km": scale,
                                "season_configuration": configuration,
                                "analysis_observations": len(analysis_rows),
                                "retained_edges": surface.retained_edges,
                                "detectable_cells": surface.detectable_cells,
                                "nonconstant_components": surface.nonconstant_components,
                                "status": surface.status,
                            }
                        )
                    evaluable_species = sum(
                        item.status == "evaluable" for item in surfaces
                    )
                    cohort_record = {
                        "cohort_id": cohort,
                        "surface_role": role,
                        "scale_km": scale,
                        "season_configuration": configuration,
                        "evaluable_species": evaluable_species,
                        "status": "not_evaluable",
                    }
                    try:
                        cohort_surface = equal_species_cohort_surface(
                            surfaces,
                            minimum_species=int(
                                transition["minimum_evaluable_species_per_primary_cohort"]
                            ),
                        )
                    except ValueError:
                        cohort_status.append(cohort_record)
                        continue
                    cohort_record["status"] = "evaluable"
                    cohort_status.append(cohort_record)
                    cohort_surfaces[cohort] = cohort_surface
                if len(cohort_surfaces) != 8:
                    continue
                for overlay_name in OVERLAY_FIELDS:
                    tests = []
                    for cohort in cohorts:
                        flower_surface, opportunity = cohort_surfaces[cohort]
                        try:
                            test = prepare_spectral_cohort_test(
                                cohort,
                                flower_surface,
                                opportunity,
                                cells,
                                cell_latitude,
                                cell_longitude,
                                {overlay_name: overlays[overlay_name]},
                                n_lon=grid.n_lon,
                                n_sinlat=grid.n_sinlat,
                                minimum_cells=int(
                                    contract["environmental_families"][
                                        "minimum_test_cells"
                                    ]
                                ),
                            )
                        except ValueError:
                            tests = []
                            break
                        tests.append(test)
                    if len(tests) == 8:
                        group = f"{role}|{scale}|{configuration}|{overlay_name}"
                        groups[group] = tests

    primary_groups = [
        f"flower|100|all_dates|{overlay}" for overlay in PRIMARY_OVERLAYS
    ]
    if not groups:
        joint: dict[str, Any] = {"groups": {}, "randomizations": 9999}
        outcome = "not_evaluable"
    else:
        joint = joint_equal_cohort_spectral_test(
            groups,
            randomizations=int(contract["joint_inference"]["randomizations"]),
            rng=np.random.default_rng(int(contract["joint_inference"]["seed"])),
        )
        if any(name not in joint["groups"] for name in primary_groups):
            outcome = "not_evaluable"
        elif any(
            joint["groups"][name]["familywise_adjusted_p"]
            <= float(contract["joint_inference"]["alpha"])
            for name in primary_groups
        ):
            outcome = "supported"
        else:
            outcome = "not_supported"

    background_primary = [
        f"background|100|all_dates|{overlay}" for overlay in PRIMARY_OVERLAYS
    ]
    if any(name not in joint["groups"] for name in background_primary):
        background_wording = "background_not_evaluable"
    elif any(
        joint["groups"][name]["familywise_adjusted_p"]
        <= float(contract["joint_inference"]["alpha"])
        for name in background_primary
    ):
        background_wording = "background_primary_supported"
    else:
        background_wording = "background_primary_not_supported"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    species_status_path = args.output_dir / "species_configuration_status.csv"
    cohort_status_path = args.output_dir / "cohort_configuration_status.csv"
    map_path = args.output_dir / "species_free_atlas_points.csv"
    write_csv(species_status_path, species_status)
    write_csv(cohort_status_path, cohort_status)
    map_rows = []
    for row in joined:
        if row["automated_colour_state_status"] != "automated_colour_state_admitted":
            continue
        lab = finite_lab(row, "flower")
        if lab is None:
            continue
        map_rows.append(
            {
                "measurement_id": row["measurement_id"],
                "cohort_id": row["cohort_id"],
                "latitude": row["latitude"],
                "longitude": row["longitude"],
                "flower_L_mean": float(lab[0]),
                "flower_a_mean": float(lab[1]),
                "flower_b_mean": float(lab[2]),
            }
        )
    write_csv(map_path, map_rows)
    result = {
        "protocol": contract["protocol"],
        "status": f"environmental_concordance_{outcome}",
        "environmental_branch_outcome": outcome,
        "geographic_branch_outcome": "not_evaluable",
        "pollinator_branch_outcome": "not_evaluable",
        "pollinator_reason": "frozen precolour Bombus source-access gate",
        "background_wording_state": background_wording,
        "background_wording": contract["background_wording"][background_wording],
        "frozen_measurements": len(joined),
        "admitted_map_points": len(map_rows),
        "species": len(species_to_cohort),
        "cohorts": len(cohorts),
        "primary_groups": primary_groups,
        "joint_inference": joint,
        "source_sha256": {
            "measurement_gate": sha256(args.measurement_gate),
            "environmental_coverage_result": sha256(
                args.environmental_coverage_result
            ),
            "roi_evidence_manifest": sha256(
                args.roi_evidence_dir / "gate_evidence_manifest.json"
            ),
            "trained_weight": trained_weight_sha,
            "sealed_coordinate_key": sha256(args.sealed_coordinate_key),
            "contract": sha256(args.contract),
        },
        "output_sha256": {
            species_status_path.name: sha256(species_status_path),
            cohort_status_path.name: sha256(cohort_status_path),
            map_path.name: sha256(map_path),
        },
        "coordinate_colour_key_joined": True,
        "claim_ceiling": contract["claim_ceiling"],
    }
    result_path = args.output_dir / "environmental_inference_result.json"
    result_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
