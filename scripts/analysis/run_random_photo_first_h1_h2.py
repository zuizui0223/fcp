#!/usr/bin/env python3
"""Run the prospectively frozen random photo-first H1 and hierarchical H2.

This entrypoint is deliberately post-measurement. It refuses partial or mismatched
measurement denominators, executes H1 with the predeclared cached null, retains
the matched 999 edge-persistence maps for H2, and only then evaluates the frozen
macroclimate concordance. H2 can never rescue a non-significant or not-evaluable
H1.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_atlas import prepare_photo_grid, species_capped_sampling_capacity
from fcp_pipeline.photo_first_h1_fast import persistence_null_test_cached
from fcp_pipeline.photo_first_h1_null_maps import persistence_null_maps_cached
from fcp_pipeline.photo_first_h2_climate import (
    aggregate_climate_to_h1_grid,
    build_edge_climate_contrasts,
    climate_concordance_test,
    climate_driver_decomposition,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid


ROOT = Path(__file__).resolve().parents[2]
MEASURED = ROOT / "data/derived/random_photo_first_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_measurement_result_v1.json"
CANDIDATE_MANIFEST = ROOT / "docs/supporting/random_photo_first_candidate_pool_manifest_v1.json"
H1_CONTRACT = ROOT / "docs/supporting/random_photo_first_boundary_persistence_contract_v1.json"
H2_CONTRACT = ROOT / "docs/supporting/random_photo_first_h2_climate_contract_v1.json"
EXECUTION_CONTRACT = ROOT / "docs/supporting/random_photo_first_inference_execution_v1.json"

H1_EDGES = ROOT / "data/derived/random_photo_first_h1_edges_v1.csv"
H1_NULL_CONCENTRATIONS = ROOT / "data/derived/random_photo_first_h1_null_concentrations_v1.csv"
H1_NULL_MAPS = ROOT / "data/derived/random_photo_first_h1_null_maps_v1.npz"
H1_QC_CELLS = ROOT / "data/derived/random_photo_first_h1_measurement_qc_cells_v1.csv"
H1_SENSITIVITIES = ROOT / "data/derived/random_photo_first_h1_sensitivities_v1.csv"
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"

H2_CLIMATE_CELLS = ROOT / "data/derived/random_photo_first_h2_climate_cells_250km_v1.csv"
H2_CLIMATE_EDGES = ROOT / "data/derived/random_photo_first_h2_climate_edges_250km_v1.csv"
H2_DRIVERS = ROOT / "data/derived/random_photo_first_h2_driver_decomposition_v1.csv"
H2_SENSITIVITIES = ROOT / "data/derived/random_photo_first_h2_sensitivities_v1.csv"
H2_RESULT = ROOT / "docs/supporting/random_photo_first_h2_result_v1.json"

BIOLOGICAL_MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
ALLOWED_MORPHS = frozenset((*BIOLOGICAL_MORPHS, "mixed_uncertain"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def validate_inputs() -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    execution = load_json(EXECUTION_CONTRACT)
    h1_contract = load_json(H1_CONTRACT)
    h2_contract = load_json(H2_CONTRACT)
    candidate = load_json(CANDIDATE_MANIFEST)
    measurement = load_json(MEASUREMENT_RESULT)

    if execution.get("protocol") != "random-photo-first-inference-execution-v1":
        raise RuntimeError("unexpected inference execution contract")
    if execution.get("status") != "prospectively_frozen_before_complete_fresh_measurement_outcome":
        raise RuntimeError("inference execution was not frozen before measurement outcome")
    if any(execution.get("outcome_firewall_at_freeze", {}).values()):
        raise RuntimeError("inference execution contract contains opened outcomes")

    frame = execution["immutable_input_frame"]
    expected_rows = int(frame["frozen_candidate_rows"])
    expected_candidate_sha = str(frame["candidate_table_sha256"])
    if candidate.get("candidate_table_sha256") != expected_candidate_sha:
        raise RuntimeError("candidate table identity changed before inference")
    if int(candidate.get("counts", {}).get("observations", -1)) != expected_rows:
        raise RuntimeError("candidate denominator changed before inference")
    if candidate.get("outcome_firewall", {}).get("legacy_pr21_terminal_records_used") is not False:
        raise RuntimeError("legacy PR21 records entered the new candidate frame")

    if measurement.get("status") != frame["measurement_manifest_required_status"]:
        raise RuntimeError("complete fresh measurement is required before H1")
    for key in ("terminal_result_rows", "joined_rows", "frozen_candidate_rows"):
        if int(measurement.get(key, -1)) != expected_rows:
            raise RuntimeError(f"measurement denominator changed: {key}")
    if measurement.get("coordinate_colour_join_opened_after_complete_measurement") is not True:
        raise RuntimeError("measurement metadata join was not opened after complete measurement")
    if measurement.get("legacy_pr21_terminal_records_used") is not False:
        raise RuntimeError("legacy PR21 terminal records entered fresh measurement")
    if measurement.get("h1_run") is not False or measurement.get("h2_run") is not False:
        raise RuntimeError("measurement artifact already contains downstream inference")
    if measurement.get("firewall_candidate_table_sha256") != expected_candidate_sha:
        raise RuntimeError("measurement firewall candidate identity changed")
    if measurement.get("measurement_table_sha256") != sha256(MEASURED):
        raise RuntimeError("measured photo table hash does not match its manifest")

    photos = pd.read_csv(MEASURED)
    required = {"measurement_id", "species", "latitude", "longitude", "morph"}
    missing = sorted(required.difference(photos.columns))
    if missing:
        raise RuntimeError(f"measured table lacks H1 fields: {missing}")
    if len(photos) != expected_rows or photos["measurement_id"].astype(str).nunique() != expected_rows:
        raise RuntimeError("measured table does not exactly cover the frozen denominator")
    observed_morphs = set(photos["morph"].astype(str))
    if not observed_morphs.issubset(ALLOWED_MORPHS):
        raise RuntimeError(f"unexpected measured morph labels: {sorted(observed_morphs - ALLOWED_MORPHS)}")

    primary = execution["h1_primary"]
    if (
        int(primary["photos_per_replicate"]) != int(h1_contract["primary_sampling"]["photos_per_replicate"])
        or int(primary["replicates"]) != int(h1_contract["primary_sampling"]["replicates"])
        or int(primary["species_cap_per_cell_per_replicate"]) != int(h1_contract["primary_sampling"]["species_cap_per_cell_per_replicate"])
        or int(primary["minimum_classifiable_photos_per_cell"]) != int(h1_contract["primary_sampling"]["minimum_classifiable_photos_per_cell"])
        or float(primary["transition_quantile"]) != float(h1_contract["transition"]["transition_quantile"])
        or int(primary["permutations"]) != int(h1_contract["primary_null"]["permutations"])
    ):
        raise RuntimeError("H1 execution parameters drifted from the frozen scientific contract")
    if h2_contract.get("hierarchical_role", {}).get("h2_must_not_rescue_a_non_significant_or_not_evaluable_h1") is not True:
        raise RuntimeError("H2 hierarchy firewall changed")
    return photos, execution, h1_contract, h2_contract, measurement


def build_qc(photos: pd.DataFrame, grid: EqualAreaGrid) -> dict[str, Any]:
    work = prepare_photo_grid(photos, grid=grid)
    rows: list[dict[str, Any]] = []
    for cell_id in range(grid.n_cells):
        cell = work.loc[work["cell_id"] == cell_id]
        raw_n = int(len(cell))
        classifiable = cell["morph"].astype(str).isin(BIOLOGICAL_MORPHS)
        class_n = int(classifiable.sum())
        mixed_n = raw_n - class_n
        species_counts = cell["species"].astype(str).value_counts() if raw_n else pd.Series(dtype=int)
        top_species = str(species_counts.index[0]) if len(species_counts) else ""
        top_species_n = int(species_counts.iloc[0]) if len(species_counts) else 0
        rows.append(
            {
                "cell_id": cell_id,
                "raw_photos": raw_n,
                "classifiable_photos": class_n,
                "mixed_uncertain_photos": mixed_n,
                "mixed_uncertain_fraction": mixed_n / raw_n if raw_n else np.nan,
                "species": int(cell["species"].astype(str).nunique()) if raw_n else 0,
                "top_species": top_species,
                "top_species_photos": top_species_n,
                "top_species_fraction": top_species_n / raw_n if raw_n else np.nan,
            }
        )
    qc = pd.DataFrame(rows)
    H1_QC_CELLS.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(H1_QC_CELLS, index=False, lineterminator="\n")
    overall_species = photos["species"].astype(str).value_counts()
    return {
        "measured_rows": int(len(photos)),
        "classifiable_rows": int(photos["morph"].astype(str).isin(BIOLOGICAL_MORPHS).sum()),
        "mixed_uncertain_rows": int(photos["morph"].astype(str).eq("mixed_uncertain").sum()),
        "occupied_h1_cells": int((qc["raw_photos"] > 0).sum()),
        "cells_with_at_least_5_classifiable": int((qc["classifiable_photos"] >= 5).sum()),
        "maximum_cell_mixed_uncertain_fraction": float(qc.loc[qc["raw_photos"] > 0, "mixed_uncertain_fraction"].max()),
        "maximum_cell_top_species_fraction": float(qc.loc[qc["raw_photos"] > 0, "top_species_fraction"].max()),
        "global_top_species": str(overall_species.index[0]),
        "global_top_species_photos": int(overall_species.iloc[0]),
        "global_top_species_fraction": float(overall_species.iloc[0] / len(photos)),
        "qc_cell_table": str(H1_QC_CELLS.relative_to(ROOT)),
    }


def run_sensitivities(photos: pd.DataFrame, execution: dict[str, Any]) -> pd.DataFrame:
    primary = execution["h1_primary"]
    spec = execution["h1_predeclared_sensitivities"]
    rows: list[dict[str, Any]] = []
    primary_key = (int(primary["n_lon"]), int(primary["n_sinlat"]), int(primary["species_cap_per_cell_per_replicate"]))
    for grid_spec in spec["factorial_grid_and_species_cap"]["grids"]:
        for cap in spec["factorial_grid_and_species_cap"]["species_caps"]:
            key = (int(grid_spec["n_lon"]), int(grid_spec["n_sinlat"]), int(cap))
            if key == primary_key:
                continue
            grid = EqualAreaGrid(n_lon=key[0], n_sinlat=key[1])
            prepared = prepare_photo_grid(photos, grid=grid)
            capacity = species_capped_sampling_capacity(prepared, species_cap_per_cell=key[2])
            row: dict[str, Any] = {
                "n_lon": key[0],
                "n_sinlat": key[1],
                "species_cap_per_cell": key[2],
                "species_capped_capacity": int(capacity),
                "status": "not_evaluable_fixed_replicate_size" if capacity < int(spec["same_photos_per_replicate"]) else "pending",
                "observed_concentration": np.nan,
                "p_upper": np.nan,
                "supported_edges": 0,
            }
            if capacity >= int(spec["same_photos_per_replicate"]):
                try:
                    observed, _, p_upper = persistence_null_test_cached(
                        photos,
                        grid=grid,
                        target_n=int(spec["same_photos_per_replicate"]),
                        n_replicates=int(spec["same_replicates"]),
                        species_cap_per_cell=key[2],
                        min_photos_per_cell=int(spec["same_minimum_classifiable_photos_per_cell"]),
                        transition_quantile=float(spec["same_transition_quantile"]),
                        n_permutations=int(spec["same_permutations"]),
                        sampling_seed=int(spec["same_sampling_seed"]),
                        permutation_seed=int(spec["same_permutation_seed"]),
                    )
                except ValueError as exc:
                    row["status"] = f"not_evaluable:{str(exc)[:300]}"
                else:
                    row.update(
                        {
                            "status": "complete_evaluable",
                            "observed_concentration": float(observed.concentration),
                            "p_upper": float(p_upper),
                            "supported_edges": int((observed.edge_table["opportunities"] > 0).sum()),
                        }
                    )
            rows.append(row)
    table = pd.DataFrame(rows)
    H1_SENSITIVITIES.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(H1_SENSITIVITIES, index=False, lineterminator="\n")
    return table


def verified_climate_source(path_text: str, expected_sha: str) -> pd.DataFrame:
    path = ROOT / path_text
    if sha256(path) != expected_sha:
        raise RuntimeError(f"frozen climate source hash changed: {path_text}")
    return pd.read_csv(path)


def run_h2(
    h1_result: dict[str, Any],
    primary_maps: Any,
    grid: EqualAreaGrid,
    h2_contract: dict[str, Any],
) -> dict[str, Any]:
    source_spec = h2_contract["environment_source"]["primary_grid_source"]
    climate_source = verified_climate_source(source_spec["path"], source_spec["content_sha256"])
    climate_cells = aggregate_climate_to_h1_grid(climate_source, grid=grid)
    climate_edges = build_edge_climate_contrasts(climate_cells, grid=grid)
    H2_CLIMATE_CELLS.parent.mkdir(parents=True, exist_ok=True)
    climate_cells.to_csv(H2_CLIMATE_CELLS, index=False, lineterminator="\n")
    climate_edges.to_csv(H2_CLIMATE_EDGES, index=False, lineterminator="\n")

    minimum_edges = int(h2_contract["primary_test"]["minimum_supported_edges"])
    primary = climate_concordance_test(
        primary_maps.observed.edge_table,
        primary_maps.null_persistence,
        primary_maps.edge_ids,
        climate_edges,
        predictor="multivariate_climate_distance",
        subset="global",
        minimum_supported_edges=minimum_edges,
    )
    drivers = climate_driver_decomposition(
        primary_maps.observed.edge_table,
        primary_maps.null_persistence,
        primary_maps.edge_ids,
        climate_edges,
        minimum_supported_edges=minimum_edges,
    )
    drivers.to_csv(H2_DRIVERS, index=False, lineterminator="\n")

    sensitivity_rows: list[dict[str, Any]] = []
    for subset in ("within_biome", "within_realm"):
        try:
            result = climate_concordance_test(
                primary_maps.observed.edge_table,
                primary_maps.null_persistence,
                primary_maps.edge_ids,
                climate_edges,
                predictor="multivariate_climate_distance",
                subset=subset,
                minimum_supported_edges=minimum_edges,
            )
        except ValueError as exc:
            sensitivity_rows.append({"type": "subset", "value": subset, "status": f"not_evaluable:{str(exc)[:300]}", "weighted_r": np.nan, "p_upper": np.nan, "supported_edges": 0})
        else:
            sensitivity_rows.append({"type": "subset", "value": subset, "status": "complete_evaluable", "weighted_r": result.statistic, "p_upper": result.p_upper, "supported_edges": result.supported_edges})

    for scale_spec in h2_contract["environment_source"]["scale_sensitivities"]:
        source = verified_climate_source(scale_spec["path"], scale_spec["content_sha256"])
        cells = aggregate_climate_to_h1_grid(source, grid=grid)
        edges = build_edge_climate_contrasts(cells, grid=grid)
        result = climate_concordance_test(
            primary_maps.observed.edge_table,
            primary_maps.null_persistence,
            primary_maps.edge_ids,
            edges,
            predictor="multivariate_climate_distance",
            subset="global",
            minimum_supported_edges=minimum_edges,
        )
        sensitivity_rows.append({"type": "climate_scale_km", "value": int(scale_spec["scale_km"]), "status": "complete_evaluable", "weighted_r": result.statistic, "p_upper": result.p_upper, "supported_edges": result.supported_edges})
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(H2_SENSITIVITIES, index=False, lineterminator="\n")

    h1_pass = float(h1_result["primary"]["p_upper"]) < float(h1_result["alpha"])
    h2_pass = float(primary.p_upper) < float(h2_contract["primary_test"]["alpha"])
    if h1_pass and h2_pass:
        hierarchy = "support_hierarchical_recurrent_boundary_macroclimate_concordance"
    elif not h1_pass:
        hierarchy = "diagnostic_only_h1_not_supported_no_climate_mechanism_claim"
    else:
        hierarchy = "h1_supported_but_h2_macroclimate_concordance_not_supported"

    payload = {
        "protocol": h2_contract["protocol"],
        "status": "complete_h2_evaluable",
        "hierarchical_decision": hierarchy,
        "h1_primary_p_upper": float(h1_result["primary"]["p_upper"]),
        "h1_primary_supported": bool(h1_pass),
        "primary": {
            "climate_source_scale_km": int(source_spec["scale_km"]),
            "weighted_r": float(primary.statistic),
            "p_upper": float(primary.p_upper),
            "alpha": float(h2_contract["primary_test"]["alpha"]),
            "supported_edges": int(primary.supported_edges),
            "support": bool(h2_pass),
        },
        "driver_decomposition": drivers.to_dict(orient="records"),
        "sensitivities": sensitivity.to_dict(orient="records"),
        "files": {
            "climate_cells": str(H2_CLIMATE_CELLS.relative_to(ROOT)),
            "climate_edges": str(H2_CLIMATE_EDGES.relative_to(ROOT)),
            "drivers": str(H2_DRIVERS.relative_to(ROOT)),
            "sensitivities": str(H2_SENSITIVITIES.relative_to(ROOT)),
        },
        "claim_ceiling": h2_contract["claim_ceiling"],
    }
    write_json(H2_RESULT, payload)
    return payload


def write_not_evaluable(h1_reason: str, measurement: dict[str, Any], qc: dict[str, Any]) -> None:
    h1_payload = {
        "protocol": "random-photo-first-boundary-persistence-v1",
        "status": "not_evaluable_h1_after_complete_measurement",
        "reason": h1_reason,
        "measurement_table_sha256": measurement.get("measurement_table_sha256"),
        "measurement_qc": qc,
        "h2_opened": False,
    }
    write_json(H1_RESULT, h1_payload)
    write_json(
        H2_RESULT,
        {
            "protocol": "random-photo-first-h2-climate-concordance-v1",
            "status": "not_evaluable_h2_because_h1_not_evaluable",
            "h1_reason": h1_reason,
            "climate_colour_join_opened": False,
        },
    )


def main() -> int:
    photos, execution, h1_contract, h2_contract, measurement = validate_inputs()
    primary_spec = execution["h1_primary"]
    grid = EqualAreaGrid(n_lon=int(primary_spec["n_lon"]), n_sinlat=int(primary_spec["n_sinlat"]))
    qc = build_qc(photos, grid)

    try:
        maps = persistence_null_maps_cached(
            photos,
            grid=grid,
            target_n=int(primary_spec["photos_per_replicate"]),
            n_replicates=int(primary_spec["replicates"]),
            species_cap_per_cell=int(primary_spec["species_cap_per_cell_per_replicate"]),
            min_photos_per_cell=int(primary_spec["minimum_classifiable_photos_per_cell"]),
            transition_quantile=float(primary_spec["transition_quantile"]),
            n_permutations=int(primary_spec["permutations"]),
            sampling_seed=int(primary_spec["sampling_seed"]),
            permutation_seed=int(primary_spec["permutation_seed"]),
        )
    except ValueError as exc:
        reason = str(exc)
        if not (
            reason.startswith("not_evaluable_")
            or "not estimable" in reason
            or "null" in reason and "not estimable" in reason
        ):
            raise
        write_not_evaluable(reason, measurement, qc)
        print(json.dumps(load_json(H1_RESULT), indent=2))
        return 0

    H1_EDGES.parent.mkdir(parents=True, exist_ok=True)
    maps.observed.edge_table.to_csv(H1_EDGES, index=False, lineterminator="\n")
    pd.DataFrame({"permutation": np.arange(1, len(maps.null_concentrations) + 1), "concentration": maps.null_concentrations}).to_csv(H1_NULL_CONCENTRATIONS, index=False, lineterminator="\n")
    np.savez_compressed(
        H1_NULL_MAPS,
        edge_ids=np.asarray(maps.edge_ids, dtype="U"),
        null_concentrations=maps.null_concentrations,
        null_persistence=maps.null_persistence,
    )
    sensitivities = run_sensitivities(photos, execution)
    alpha = float(primary_spec["alpha"])
    supported = maps.observed.edge_table.loc[maps.observed.edge_table["opportunities"] > 0]
    decision = (
        "support_excess_recurrent_boundary_concentration"
        if maps.p_upper < alpha
        else "no_support_excess_recurrent_boundary_concentration"
    )
    h1_payload = {
        "protocol": h1_contract["protocol"],
        "status": "complete_h1_evaluable",
        "decision": decision,
        "alpha": alpha,
        "primary": {
            "observed_concentration": float(maps.observed.concentration),
            "realized_transition_rate": float(maps.observed.transition_rate),
            "p_upper": float(maps.p_upper),
            "null_mean": float(np.mean(maps.null_concentrations)),
            "null_median": float(np.median(maps.null_concentrations)),
            "null_q025": float(np.quantile(maps.null_concentrations, 0.025)),
            "null_q975": float(np.quantile(maps.null_concentrations, 0.975)),
            "supported_edges": int(len(supported)),
            "maximum_edge_persistence": float(supported["persistence"].max()) if len(supported) else None,
            "mean_sampled_photos": float(maps.observed.mean_sampled_photos),
            "permutations": int(len(maps.null_concentrations)),
            "sampling_seed": int(primary_spec["sampling_seed"]),
            "permutation_seed": int(primary_spec["permutation_seed"]),
        },
        "measurement_table_sha256": measurement["measurement_table_sha256"],
        "candidate_table_sha256": execution["immutable_input_frame"]["candidate_table_sha256"],
        "measurement_qc": qc,
        "sensitivity_summary": sensitivities.to_dict(orient="records"),
        "files": {
            "edge_table": str(H1_EDGES.relative_to(ROOT)),
            "null_concentrations": str(H1_NULL_CONCENTRATIONS.relative_to(ROOT)),
            "null_persistence_maps": str(H1_NULL_MAPS.relative_to(ROOT)),
            "measurement_qc_cells": str(H1_QC_CELLS.relative_to(ROOT)),
            "sensitivities": str(H1_SENSITIVITIES.relative_to(ROOT)),
        },
        "file_sha256": {
            "edge_table": sha256(H1_EDGES),
            "null_concentrations": sha256(H1_NULL_CONCENTRATIONS),
            "null_persistence_maps": sha256(H1_NULL_MAPS),
            "measurement_qc_cells": sha256(H1_QC_CELLS),
            "sensitivities": sha256(H1_SENSITIVITIES),
        },
    }
    write_json(H1_RESULT, h1_payload)
    h2_payload = run_h2(h1_payload, maps, grid, h2_contract)

    print(json.dumps({"h1": h1_payload, "h2": h2_payload}, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
