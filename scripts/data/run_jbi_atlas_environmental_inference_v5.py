#!/usr/bin/env python3
"""Run the reached v5 environmental branch on the already-frozen flower Lab field.

The primary flower field is never remeasured or restandardized here.  This runner
uses the exact standardized flower table emitted by ``run_jbi_atlas_real_inference_v5.py``.
The background camera/scene diagnostic is standardized once from the same frozen
measurement bundle, as required by the pre-pixel colour-surface contract, and is
included in the joint familywise null without being able to create or erase the
primary flower outcome.
"""

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
    SpeciesTransitionSurface,
    equal_species_cohort_surface,
    joint_equal_cohort_spectral_test,
    prepare_spectral_cohort_test,
    robust_standardize_lab,
    validate_colour_inference_contract,
)
from fcp_pipeline.atlas_inference_cascade import next_step, validate_contract as validate_inference_v5
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from fcp_pipeline.atlas_real_inference_v5 import validate_real_inference_amendment
from fcp_pipeline.continuous_colour_boundaries import average_rank_intensity, edge_colour_discontinuity
from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    build_edge_cell_geometry,
    cell_mean_intensity,
    equal_area_cell_centers,
)
from fcp_pipeline.spatial_graph import spherical_knn_edges
from scripts.data.evaluate_jbi_atlas_measurement_gate_v5 import load_complete_v5_measurement_bundle


PROTOCOL = "jbi-atlas-environmental-concordance-v5-v1"
FLOWER_FIELDS = (
    "flower_L_standardized",
    "flower_a_standardized",
    "flower_b_standardized",
)
BACKGROUND_RAW_FIELDS = (
    "background_L_mean",
    "background_a_mean",
    "background_b_mean",
)
BACKGROUND_STD_FIELDS = (
    "background_L_standardized",
    "background_a_standardized",
    "background_b_standardized",
)
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
EXPECTED_COHORTS = tuple(f"C{index:02d}" for index in range(1, 9))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"cannot infer columns for empty table: {path}")
        ordered = list(rows[0])
        union = {key for row in rows for key in row}
        ordered.extend(sorted(union - set(ordered)))
    else:
        ordered = list(fields)
        union = {key for row in rows for key in row}
        ordered.extend(sorted(union - set(ordered)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def finite_vector(row: Mapping[str, Any], fields: Sequence[str]) -> np.ndarray | None:
    try:
        values = np.asarray([float(row[field]) for field in fields], dtype=float)
    except (KeyError, TypeError, ValueError):
        return None
    return values if values.shape == (3,) and np.isfinite(values).all() else None


def season_labels(rows: Sequence[Mapping[str, Any]], configuration: str) -> np.ndarray | None:
    if configuration == "all_dates":
        return None
    if configuration == "same_calendar_month_edges":
        return np.asarray([int(row["observed_month"]) for row in rows], dtype=int)
    if configuration == "same_local_solar_quarter_edges":
        return np.asarray([int(row["local_solar_quarter"]) for row in rows], dtype=int)
    raise ValueError(f"unknown season configuration: {configuration}")


def _season_edge_subset(
    edges: np.ndarray,
    distances: np.ndarray,
    labels: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    if labels is None:
        return edges, distances
    if labels.ndim != 1 or len(labels) <= int(edges.max(initial=-1)):
        raise ValueError("season labels must match observations")
    keep = labels[edges[:, 0]] == labels[edges[:, 1]]
    return edges[keep], distances[keep]


def build_surface_from_standardized_lab(
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    *,
    grid: EqualAreaGrid,
    scale_km: int,
    configuration: str,
    knn_k: int = 5,
    minimum_edges_per_cell: int = 2,
    minimum_retained_edges: int = 100,
    minimum_detectable_cells: int = 10,
) -> SpeciesTransitionSurface:
    """Build one transition surface without any downstream Lab restandardization."""

    empty = np.full(grid.n_cells, np.nan, dtype=float)
    detectable_empty = np.zeros(grid.n_cells, dtype=bool)
    try:
        if len(rows) < 2:
            raise ValueError("fewer_than_two_standardized_rows")
        lab = np.vstack([finite_vector(row, fields) for row in rows])
        if lab.shape != (len(rows), 3) or not np.isfinite(lab).all():
            raise ValueError("standardized_lab_not_finite")
        q25, q75 = np.quantile(lab, [0.25, 0.75], axis=0)
        variable = (q75 - q25) > 0
        if not np.any(variable):
            raise ValueError("all_standardized_components_constant")
        latitude = np.asarray([float(row["latitude"]) for row in rows], dtype=float)
        longitude = np.asarray([float(row["longitude"]) for row in rows], dtype=float)
        edges, distances = spherical_knn_edges(latitude, longitude, k=knn_k)
        edges, distances = _season_edge_subset(
            edges,
            distances,
            season_labels(rows, configuration),
        )
        geometry = build_edge_cell_geometry(
            latitude,
            longitude,
            edges,
            distances,
            grid=grid,
            max_edge_km=float(scale_km),
            min_edges_per_cell=minimum_edges_per_cell,
        )
    except (TypeError, ValueError):
        return SpeciesTransitionSurface(
            status="not_evaluable",
            surface=empty,
            detectable=detectable_empty,
            retained_edges=0,
            detectable_cells=0,
            nonconstant_components=0,
            geometry=None,
        )

    retained_edges = len(geometry.retained_edges)
    detectable_cells = int(np.count_nonzero(geometry.detectable))
    if retained_edges < minimum_retained_edges or detectable_cells < minimum_detectable_cells:
        return SpeciesTransitionSurface(
            status="not_evaluable",
            surface=empty,
            detectable=geometry.detectable,
            retained_edges=retained_edges,
            detectable_cells=detectable_cells,
            nonconstant_components=int(np.count_nonzero(variable)),
            geometry=geometry,
        )
    scores = edge_colour_discontinuity(lab, geometry.retained_edges)
    surface = cell_mean_intensity(average_rank_intensity(scores), geometry)
    return SpeciesTransitionSurface(
        status="evaluable",
        surface=surface,
        detectable=geometry.detectable,
        retained_edges=retained_edges,
        detectable_cells=detectable_cells,
        nonconstant_components=int(np.count_nonzero(variable)),
        geometry=geometry,
    )


def environment_vectors(
    path: Path,
    *,
    grid: EqualAreaGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    rows = read_csv(path)
    cells, latitude, longitude = equal_area_cell_centers(grid)
    overlays = {field: np.full(grid.n_cells, np.nan, dtype=float) for field in OVERLAY_FIELDS}
    seen: set[int] = set()
    for row in rows:
        cell = int(row["cell_id"])
        if cell in seen or cell < 0 or cell >= grid.n_cells:
            raise ValueError(f"invalid duplicate environmental cell: {cell}")
        seen.add(cell)
        for field in OVERLAY_FIELDS:
            raw = row.get(field, "")
            if raw not in (None, ""):
                value = float(raw)
                if math.isfinite(value):
                    overlays[field][cell] = value
    return cells, latitude, longitude, overlays


def freeze_background_standardized_field(
    flower_rows: Sequence[Mapping[str, Any]],
    measurement_by_id: Mapping[str, Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Freeze the predeclared background diagnostic once on the flower-field support."""

    rows_by_species: dict[str, list[tuple[Mapping[str, Any], np.ndarray]]] = {}
    for flower in flower_rows:
        measurement_id = str(flower["measurement_id"])
        measured = measurement_by_id.get(measurement_id)
        if measured is None:
            raise RuntimeError("standardized flower row is absent from measurement bundle")
        if measured.get("automated_colour_state_status") != "automated_colour_state_admitted":
            raise RuntimeError("standardized flower row is not admitted in frozen measurements")
        if measured.get("background_features_available") is not True:
            continue
        vector = finite_vector(measured, BACKGROUND_RAW_FIELDS)
        if vector is None:
            continue
        rows_by_species.setdefault(str(flower["species_blind_id"]), []).append((flower, vector))

    output: list[dict[str, Any]] = []
    status: list[dict[str, Any]] = []
    all_species = sorted({str(row["species_blind_id"]) for row in flower_rows})
    for species_id in all_species:
        pairs = rows_by_species.get(species_id, [])
        if len(pairs) < 2:
            status.append(
                {
                    "species_blind_id": species_id,
                    "status": "not_evaluable",
                    "reason": "fewer_than_two_background_rows",
                    "background_rows": len(pairs),
                }
            )
            continue
        raw = np.vstack([vector for _, vector in pairs])
        try:
            standardized, variable = robust_standardize_lab(raw)
        except ValueError as exc:
            status.append(
                {
                    "species_blind_id": species_id,
                    "status": "not_evaluable",
                    "reason": str(exc),
                    "background_rows": len(pairs),
                }
            )
            continue
        for (flower, _), z in zip(pairs, standardized, strict=True):
            output.append(
                {
                    "measurement_id": str(flower["measurement_id"]),
                    "species_blind_id": species_id,
                    "cohort_id": str(flower["cohort_id"]),
                    "latitude": float(flower["latitude"]),
                    "longitude": float(flower["longitude"]),
                    "observed_month": int(flower["observed_month"]),
                    "local_solar_quarter": int(flower["local_solar_quarter"]),
                    "background_L_standardized": float(z[0]),
                    "background_a_standardized": float(z[1]),
                    "background_b_standardized": float(z[2]),
                }
            )
        status.append(
            {
                "species_blind_id": species_id,
                "status": "background_field_evaluable",
                "background_rows": len(pairs),
                "nonzero_iqr_components": int(np.count_nonzero(variable)),
            }
        )
    return output, status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-results-dir", type=Path, required=True)
    parser.add_argument("--measurement-gate", type=Path, required=True)
    parser.add_argument("--primary-inference-dir", type=Path, required=True)
    parser.add_argument("--environmental-coverage-result", type=Path, required=True)
    parser.add_argument("--environment-dir", type=Path, default=ROOT / "data/atlas/environment")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--measurement-contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json",
    )
    parser.add_argument(
        "--inference-v5",
        type=Path,
        default=ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json",
    )
    parser.add_argument(
        "--real-inference",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_real_colour_inference_amendment_v5.json",
    )
    parser.add_argument(
        "--colour-surface-contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_colour_surface_contract_v1.json",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    measurement_contract = load_json(args.measurement_contract)
    inference = load_json(args.inference_v5)
    real_contract = load_json(args.real_inference)
    surface_contract = load_json(args.colour_surface_contract)
    validate_inference_v5(inference)
    validate_real_inference_amendment(real_contract)
    validate_colour_inference_contract(surface_contract)
    validate_measurement_execution_contract(measurement_contract, inference)

    measurement_gate = load_json(args.measurement_gate)
    if (
        measurement_gate.get("status") != "pass_scaleout_measurement_completeness"
        or measurement_gate.get("coordinate_join_permitted") is not True
        or measurement_gate.get("coordinates_opened") is not False
    ):
        raise RuntimeError("environmental branch cannot run without a passed v5 measurement gate")
    measurements, bundle = load_complete_v5_measurement_bundle(
        args.measurement_results_dir,
        contract=measurement_contract,
        inference=inference,
    )
    if bundle != measurement_gate.get("measurement_bundle"):
        raise RuntimeError("measurement bundle changed after completeness decision")
    measurement_by_id = {str(row["measurement_id"]): row for row in measurements}

    primary_result_path = args.primary_inference_dir / "real_colour_primary_inference_v5.json"
    primary_manifest_path = args.primary_inference_dir / "real_colour_primary_inference_v5_manifest.json"
    standardized_path = args.primary_inference_dir / "frozen_standardized_flower_lab_v5.csv"
    primary = load_json(primary_result_path)
    primary_manifest = load_json(primary_manifest_path)
    if (
        primary.get("protocol") != real_contract["protocol"]
        or primary.get("status") != "complete_reached_primary_spatial_branches_v5"
        or primary.get("next_confirmatory_branch") != "environmental_concordance"
        or primary.get("same_standardized_colour_field_required_for_all_downstream_branches") is not True
        or primary.get("standardized_colour_field_sha256") != sha256(standardized_path)
        or primary_manifest.get("files", {}).get(standardized_path.name) != sha256(standardized_path)
        or primary_manifest.get("files", {}).get(primary_result_path.name) != sha256(primary_result_path)
    ):
        raise RuntimeError("primary v5 result does not authorize environmental concordance")

    coverage = load_json(args.environmental_coverage_result)
    if (
        coverage.get("status") != "pass_precolour_environmental_coverage"
        or coverage.get("source_stage") != "final-source-v5"
        or coverage.get("scaleout_colour_opened") is not False
        or not {"macroclimate", "land_cover", "ecoregion"}.issubset(
            set(coverage.get("evaluable_families", []))
        )
    ):
        raise RuntimeError("frozen pre-colour environmental coverage did not pass")

    flower_rows = read_csv(standardized_path)
    if not flower_rows:
        raise RuntimeError("frozen standardized flower field is empty")
    if "species" in flower_rows[0] or "inat_taxon_id" in flower_rows[0]:
        raise RuntimeError("standardized inferential field leaked taxon identity")
    species_to_cohort: dict[str, str] = {}
    flower_by_species: dict[str, list[dict[str, Any]]] = {}
    for raw in flower_rows:
        species_id = str(raw["species_blind_id"])
        cohort = str(raw["cohort_id"])
        prior = species_to_cohort.setdefault(species_id, cohort)
        if prior != cohort:
            raise RuntimeError("one blinded species crossed frozen cohorts")
        flower_by_species.setdefault(species_id, []).append(raw)
    if tuple(sorted(set(species_to_cohort.values()))) != EXPECTED_COHORTS:
        raise RuntimeError("standardized flower field does not contain all eight frozen cohorts")

    background_rows, background_field_status = freeze_background_standardized_field(
        flower_rows,
        measurement_by_id,
    )
    background_by_species: dict[str, list[dict[str, Any]]] = {}
    for row in background_rows:
        background_by_species.setdefault(str(row["species_blind_id"]), []).append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    background_path = args.output_dir / "frozen_standardized_background_lab_v5.csv"
    write_csv(
        background_path,
        background_rows,
        fields=(
            "measurement_id",
            "species_blind_id",
            "cohort_id",
            "latitude",
            "longitude",
            "observed_month",
            "local_solar_quarter",
            *BACKGROUND_STD_FIELDS,
        ),
    )
    background_status_path = args.output_dir / "background_field_status_v5.csv"
    write_csv(background_status_path, background_field_status)

    transition = surface_contract["transition_surface"]
    groups: dict[str, list[Any]] = {}
    species_status: list[dict[str, Any]] = []
    cohort_status: list[dict[str, Any]] = []
    roles = {
        "flower": (flower_by_species, FLOWER_FIELDS),
        "background": (background_by_species, BACKGROUND_STD_FIELDS),
    }
    for scale in transition["scales_km"]:
        n_lon, n_sinlat = transition["grids"][str(scale)]
        grid = EqualAreaGrid(int(n_lon), int(n_sinlat))
        cells, cell_latitude, cell_longitude, overlays = environment_vectors(
            args.environment_dir / f"environmental_boundary_cells_{scale}km.csv",
            grid=grid,
        )
        for configuration in SEASON_CONFIGURATIONS:
            for role, (rows_by_species, fields) in roles.items():
                cohort_surfaces: dict[str, tuple[np.ndarray, np.ndarray]] = {}
                for cohort in EXPECTED_COHORTS:
                    surfaces: list[SpeciesTransitionSurface] = []
                    cohort_species = sorted(
                        species_id
                        for species_id, assigned in species_to_cohort.items()
                        if assigned == cohort
                    )
                    for species_id in cohort_species:
                        rows = rows_by_species.get(species_id, [])
                        surface = build_surface_from_standardized_lab(
                            rows,
                            fields,
                            grid=grid,
                            scale_km=int(scale),
                            configuration=configuration,
                            knn_k=int(transition["knn_k"]),
                            minimum_edges_per_cell=int(transition["minimum_edges_per_species_cell"]),
                            minimum_retained_edges=int(
                                transition["minimum_retained_edges_per_species_configuration"]
                            ),
                            minimum_detectable_cells=int(
                                transition["minimum_detectable_cells_per_species_configuration"]
                            ),
                        )
                        surfaces.append(surface)
                        species_status.append(
                            {
                                "species_blind_id": species_id,
                                "cohort_id": cohort,
                                "surface_role": role,
                                "scale_km": int(scale),
                                "season_configuration": configuration,
                                "analysis_observations": len(rows),
                                "retained_edges": surface.retained_edges,
                                "detectable_cells": surface.detectable_cells,
                                "nonconstant_components": surface.nonconstant_components,
                                "status": surface.status,
                            }
                        )
                    evaluable_species = sum(item.status == "evaluable" for item in surfaces)
                    record = {
                        "cohort_id": cohort,
                        "surface_role": role,
                        "scale_km": int(scale),
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
                        cohort_status.append(record)
                        continue
                    record["status"] = "evaluable"
                    cohort_status.append(record)
                    cohort_surfaces[cohort] = cohort_surface
                if len(cohort_surfaces) != 8:
                    continue
                for overlay_name in OVERLAY_FIELDS:
                    tests = []
                    for cohort in EXPECTED_COHORTS:
                        cohort_surface, opportunity = cohort_surfaces[cohort]
                        try:
                            test = prepare_spectral_cohort_test(
                                cohort,
                                cohort_surface,
                                opportunity,
                                cells,
                                cell_latitude,
                                cell_longitude,
                                {overlay_name: overlays[overlay_name]},
                                n_lon=grid.n_lon,
                                n_sinlat=grid.n_sinlat,
                                minimum_cells=int(surface_contract["environmental_families"]["minimum_test_cells"]),
                            )
                        except ValueError:
                            tests = []
                            break
                        tests.append(test)
                    if len(tests) == 8:
                        groups[f"{role}|{scale}|{configuration}|{overlay_name}"] = tests

    primary_groups = [f"flower|100|all_dates|{name}" for name in PRIMARY_OVERLAYS]
    if not groups:
        joint: dict[str, Any] = {"groups": {}, "randomizations": int(real_contract["environmental_concordance"]["randomizations"])}
        outcome = "not_evaluable"
    else:
        joint = joint_equal_cohort_spectral_test(
            groups,
            randomizations=int(real_contract["environmental_concordance"]["randomizations"]),
            rng=np.random.default_rng(int(real_contract["environmental_concordance"]["seed"])),
        )
        if any(name not in joint["groups"] for name in primary_groups):
            outcome = "not_evaluable"
        elif any(
            float(joint["groups"][name]["familywise_adjusted_p"])
            <= float(surface_contract["joint_inference"]["alpha"])
            for name in primary_groups
        ):
            outcome = "supported"
        else:
            outcome = "unsupported"

    background_primary = [f"background|100|all_dates|{name}" for name in PRIMARY_OVERLAYS]
    if any(name not in joint.get("groups", {}) for name in background_primary):
        background_wording_state = "background_not_evaluable"
    elif any(
        float(joint["groups"][name]["familywise_adjusted_p"])
        <= float(surface_contract["joint_inference"]["alpha"])
        for name in background_primary
    ):
        background_wording_state = "background_primary_supported"
    else:
        background_wording_state = "background_primary_not_supported"

    species_status_path = args.output_dir / "environmental_species_configuration_status_v5.csv"
    cohort_status_path = args.output_dir / "environmental_cohort_configuration_status_v5.csv"
    write_csv(species_status_path, species_status)
    write_csv(cohort_status_path, cohort_status)
    result = {
        "protocol": PROTOCOL,
        "inference_version": inference["version"],
        "status": "complete_environmental_concordance_v5",
        "branch": "environmental_concordance",
        "outcome": outcome,
        "next_confirmatory_branch": next_step("environmental_concordance", outcome, inference),
        "standardized_flower_colour_field_sha256": sha256(standardized_path),
        "standardized_flower_field_remeasured": False,
        "standardized_flower_field_restandardized": False,
        "standardized_background_diagnostic_sha256": sha256(background_path),
        "background_standardized_once_from_frozen_measurements": True,
        "primary_groups": primary_groups,
        "joint_inference": joint,
        "background_wording_state": background_wording_state,
        "background_wording": surface_contract["background_wording"][background_wording_state],
        "parents_sha256": {
            "measurement_gate": sha256(args.measurement_gate),
            "primary_inference_result": sha256(primary_result_path),
            "primary_standardized_flower_field": sha256(standardized_path),
            "environmental_coverage_result": sha256(args.environmental_coverage_result),
            "measurement_contract": sha256(args.measurement_contract),
            "inference_v5": sha256(args.inference_v5),
            "real_inference": sha256(args.real_inference),
            "colour_surface_contract": sha256(args.colour_surface_contract),
        },
        "claim_ceiling": real_contract["publication_ceiling"],
    }
    result_path = args.output_dir / "environmental_concordance_v5.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "status": result["status"],
        "branch": result["branch"],
        "outcome": outcome,
        "next_confirmatory_branch": result["next_confirmatory_branch"],
        "files": {
            result_path.name: sha256(result_path),
            background_path.name: sha256(background_path),
            background_status_path.name: sha256(background_status_path),
            species_status_path.name: sha256(species_status_path),
            cohort_status_path.name: sha256(cohort_status_path),
        },
    }
    (args.output_dir / "environmental_concordance_v5_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
