#!/usr/bin/env python3
"""Freeze the terminal Lab field once, then execute reached v5 primary spatial branches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_colour_inference import robust_standardize_lab
from fcp_pipeline.atlas_inference_cascade import next_step, validate_contract as validate_inference_v5
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from fcp_pipeline.atlas_real_inference_v5 import (
    COHORTS,
    FrozenSpeciesColourState,
    run_shared_transition_test,
    run_spatial_organization_test,
    validate_real_inference_amendment,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, build_edge_cell_geometry
from fcp_pipeline.spatial_graph import spherical_knn_edges
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs
from scripts.data.evaluate_jbi_atlas_measurement_gate_v5 import load_complete_v5_measurement_bundle


DEFAULT_MEASUREMENT_CONTRACT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
DEFAULT_INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"
DEFAULT_REAL = ROOT / "docs/supporting/jbi_atlas_real_colour_inference_amendment_v5.json"
DEFAULT_SHARED = ROOT / "docs/supporting/jbi_atlas_shared_transition_v5_signal_recovery_result.json"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    fields = list(rows[0])
    union = {key for row in rows for key in row}
    fields.extend(sorted(union - set(fields)))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_flower_lab(row: Mapping[str, Any]) -> np.ndarray | None:
    try:
        values = np.asarray(
            [float(row["flower_L_mean"]), float(row["flower_a_mean"]), float(row["flower_b_mean"])],
            dtype=float,
        )
    except (KeyError, TypeError, ValueError):
        return None
    return values if np.isfinite(values).all() else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-results-dir", type=Path, required=True)
    parser.add_argument("--measurement-gate", type=Path, required=True)
    parser.add_argument("--measurement-firewall", type=Path, required=True)
    parser.add_argument("--sealed-coordinate-key", type=Path, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT_CONTRACT)
    parser.add_argument("--inference-v5", type=Path, default=DEFAULT_INFERENCE)
    parser.add_argument("--real-inference", type=Path, default=DEFAULT_REAL)
    parser.add_argument("--shared-qualification", type=Path, default=DEFAULT_SHARED)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def _authorize_coordinate_join(
    *,
    measurement_gate: Mapping[str, Any],
    firewall: Mapping[str, Any],
    coordinate_key: Path,
    contract: Mapping[str, Any],
    inference: Mapping[str, Any],
) -> None:
    if (
        measurement_gate.get("status") != "pass_scaleout_measurement_completeness"
        or measurement_gate.get("protocol") != contract["protocol"]
        or measurement_gate.get("inference_version") != inference["version"]
        or measurement_gate.get("coordinate_join_permitted") is not True
        or measurement_gate.get("coordinates_opened") is not False
        or measurement_gate.get("frozen_measurements") != 60000
        or measurement_gate.get("superseded_v3_ordered_inference_used") is not False
    ):
        raise RuntimeError("v5 measurement completeness does not authorize coordinate join")
    if (
        firewall.get("status") != "pass_scaleout_measurement_firewall_v5"
        or firewall.get("protocol") != contract["protocol"]
        or firewall.get("inference_version") != inference["version"]
        or firewall.get("candidate_image_pixels_opened") is not False
        or firewall.get("superseded_v3_ordered_inference_used") is not False
        or firewall.get("sealed_keys", {}).get(coordinate_key.name) != sha256(coordinate_key)
    ):
        raise RuntimeError("v5 firewall does not authorize the sealed coordinate key")


def main() -> int:
    args = parse_args()
    measurement_contract = load_json(args.measurement_contract)
    inference = load_json(args.inference_v5)
    real = load_json(args.real_inference)
    shared_qualification = load_json(args.shared_qualification)
    validate_inference_v5(inference)
    validate_real_inference_amendment(real)
    validate_measurement_execution_contract(measurement_contract, inference)
    verify_repo_parent_blobs(measurement_contract)
    if (
        shared_qualification.get("status") != "pass"
        or shared_qualification.get("scope", {}).get("method_gate_only") is not True
        or shared_qualification.get("scope", {}).get("biological_support_claimed") is not False
    ):
        raise RuntimeError("shared-transition preimage qualification is not a method-only pass")

    measurement_gate = load_json(args.measurement_gate)
    firewall = load_json(args.measurement_firewall)
    _authorize_coordinate_join(
        measurement_gate=measurement_gate,
        firewall=firewall,
        coordinate_key=args.sealed_coordinate_key,
        contract=measurement_contract,
        inference=inference,
    )
    measurements, measured_bundle = load_complete_v5_measurement_bundle(
        args.measurement_results_dir,
        contract=measurement_contract,
        inference=inference,
    )
    if measured_bundle != measurement_gate.get("measurement_bundle"):
        raise RuntimeError("measurement bundle changed after completeness decision")
    measurement_by_id = {str(row["measurement_id"]): row for row in measurements}

    # First protected coordinate read occurs only after every validation above.
    coordinates = read_csv(args.sealed_coordinate_key)
    coordinate_by_id = {str(row["measurement_id"]): row for row in coordinates}
    if (
        len(coordinates) != 60000
        or len(coordinate_by_id) != 60000
        or set(coordinate_by_id) != set(measurement_by_id)
    ):
        raise RuntimeError("sealed coordinate denominator differs from v5 measurements")

    evaluable_species = {
        str(row["species_blind_id"])
        for row in measurement_gate.get("species_results", [])
        if row.get("status") == "measurement_evaluable"
    }
    if len(evaluable_species) < 160:
        raise RuntimeError("measurement gate says pass but fewer than 160 species are evaluable")

    raw_rows: list[dict[str, Any]] = []
    by_species: dict[str, list[dict[str, Any]]] = {species: [] for species in evaluable_species}
    for measurement_id in sorted(measurement_by_id):
        measured = measurement_by_id[measurement_id]
        coordinate = coordinate_by_id[measurement_id]
        if str(measured["species_blind_id"]) != str(coordinate["species_blind_id"]):
            raise RuntimeError("species-blind ID changed at protected coordinate join")
        species_id = str(measured["species_blind_id"])
        if species_id not in evaluable_species:
            continue
        if measured.get("automated_colour_state_status") != "automated_colour_state_admitted":
            continue
        lab = finite_flower_lab(measured)
        if lab is None:
            continue
        row = {
            "measurement_id": measurement_id,
            "species_blind_id": species_id,
            "cohort_id": str(coordinate["cohort_id"]),
            "latitude": float(coordinate["latitude"]),
            "longitude": float(coordinate["longitude"]),
            "observed_month": int(coordinate["observed_month"]),
            "local_solar_quarter": int(coordinate["local_solar_quarter"]),
            "flower_L_mean": float(lab[0]),
            "flower_a_mean": float(lab[1]),
            "flower_b_mean": float(lab[2]),
        }
        raw_rows.append(row)
        by_species[species_id].append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / "frozen_admitted_flower_lab_v5.csv"
    write_csv(raw_path, raw_rows)

    grid = EqualAreaGrid(320, 160)
    standardized_rows: list[dict[str, Any]] = []
    states: list[FrozenSpeciesColourState] = []
    state_ledger: list[dict[str, Any]] = []
    for species_id in sorted(evaluable_species):
        rows = by_species.get(species_id, [])
        if not rows:
            state_ledger.append(
                {"species_blind_id": species_id, "cohort_id": "", "status": "not_evaluable", "reason": "no_finite_admitted_flower_lab"}
            )
            continue
        cohorts = {str(row["cohort_id"]) for row in rows}
        if len(cohorts) != 1:
            raise RuntimeError("one blinded species crossed frozen cohorts after join")
        cohort = next(iter(cohorts))
        labs = np.asarray(
            [[row["flower_L_mean"], row["flower_a_mean"], row["flower_b_mean"]] for row in rows],
            dtype=float,
        )
        try:
            standardized, variable = robust_standardize_lab(labs)
            latitude = np.asarray([row["latitude"] for row in rows], dtype=float)
            longitude = np.asarray([row["longitude"] for row in rows], dtype=float)
            edges, distances = spherical_knn_edges(latitude, longitude, k=5)
            geometry = build_edge_cell_geometry(
                latitude,
                longitude,
                edges,
                distances,
                grid=grid,
                max_edge_km=100.0,
                min_edges_per_cell=2,
            )
            if len(geometry.retained_edges) < 100:
                raise ValueError("fewer_than_100_retained_edges")
            detectable_cells = int(np.count_nonzero(geometry.detectable))
            if detectable_cells < 10:
                raise ValueError("fewer_than_10_detectable_cells")
        except ValueError as exc:
            state_ledger.append(
                {
                    "species_blind_id": species_id,
                    "cohort_id": cohort,
                    "status": "not_evaluable",
                    "reason": str(exc),
                }
            )
            continue
        for row, z in zip(rows, standardized, strict=True):
            standardized_rows.append(
                {
                    **row,
                    "flower_L_standardized": float(z[0]),
                    "flower_a_standardized": float(z[1]),
                    "flower_b_standardized": float(z[2]),
                }
            )
        states.append(
            FrozenSpeciesColourState(
                species_id=species_id,
                cohort_id=cohort,
                standardized_lab=standardized,
                geometry=geometry,
            )
        )
        state_ledger.append(
            {
                "species_blind_id": species_id,
                "cohort_id": cohort,
                "status": "spatial_primary_evaluable",
                "admitted_rows": len(rows),
                "nonzero_iqr_components": int(np.count_nonzero(variable)),
                "retained_edges_100km": int(len(geometry.retained_edges)),
                "detectable_cells_100km": detectable_cells,
            }
        )

    standardized_path = args.output_dir / "frozen_standardized_flower_lab_v5.csv"
    write_csv(standardized_path, standardized_rows)
    state_path = args.output_dir / "species_primary_state_v5.csv"
    write_csv(state_path, state_ledger)
    colour_field_hash = sha256(standardized_path)

    cohort_counts = {cohort: sum(item.cohort_id == cohort for item in states) for cohort in COHORTS}
    terminal_evaluable = len(states) >= 160 and all(cohort_counts[cohort] >= 20 for cohort in COHORTS)
    branches: list[dict[str, Any]] = []
    if not terminal_evaluable:
        spatial = {
            "branch": "species_conditioned_spatial_organization",
            "outcome": "not_evaluable",
            "reason": "primary_spatial_evaluability_below_frozen_160_and_20_per_cohort_gate",
            "evaluable_species": len(states),
            "evaluable_species_by_cohort": cohort_counts,
        }
        branches.append(spatial)
        next_branch = "STOP_CONFIRMATORY"
    else:
        spatial_rule = real["species_conditioned_spatial_organization"]
        spatial = run_spatial_organization_test(
            states,
            inference_v5=inference,
            randomizations=int(spatial_rule["randomizations"]),
            seed=int(spatial_rule["seed"]),
            require_terminal=True,
        )
        branches.append(spatial)
        next_branch = next_step(spatial["branch"], spatial["outcome"], inference)
        if next_branch == "shared_transition":
            shared_rule = real["shared_transition"]
            shared = run_shared_transition_test(
                states,
                inference_v5=inference,
                qualification_passed=True,
                randomizations=int(shared_rule["randomizations"]),
                seed=int(shared_rule["seed"]),
                high_transition_quantile=float(shared_rule["high_transition_quantile"]),
                min_detectable_species=int(shared_rule["minimum_detectable_species_per_tested_cell"]),
                require_terminal=True,
            )
            branches.append(shared)
            next_branch = next_step(shared["branch"], shared["outcome"], inference)

    result = {
        "protocol": real["protocol"],
        "status": "complete_reached_primary_spatial_branches_v5",
        "candidate_image_pixels_opened": True,
        "coordinate_colour_join_performed": True,
        "superseded_v3_ordered_inference_used": False,
        "raw_admitted_colour_field_sha256": sha256(raw_path),
        "standardized_colour_field_sha256": colour_field_hash,
        "same_standardized_colour_field_required_for_all_downstream_branches": True,
        "spatial_primary_evaluable_species": len(states),
        "spatial_primary_evaluable_species_by_cohort": cohort_counts,
        "branches": branches,
        "next_confirmatory_branch": next_branch,
        "parents_sha256": {
            "measurement_gate": sha256(args.measurement_gate),
            "measurement_firewall": sha256(args.measurement_firewall),
            "sealed_coordinate_key": sha256(args.sealed_coordinate_key),
            "measurement_contract": sha256(args.measurement_contract),
            "inference_v5": sha256(args.inference_v5),
            "real_inference": sha256(args.real_inference),
            "shared_qualification": sha256(args.shared_qualification),
        },
        "claim_ceiling": real["publication_ceiling"],
    }
    result_path = args.output_dir / "real_colour_primary_inference_v5.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "status": result["status"],
        "next_confirmatory_branch": next_branch,
        "files": {
            raw_path.name: sha256(raw_path),
            standardized_path.name: sha256(standardized_path),
            state_path.name: sha256(state_path),
            result_path.name: sha256(result_path),
        },
    }
    (args.output_dir / "real_colour_primary_inference_v5_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
