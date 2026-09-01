#!/usr/bin/env python3
"""Run the frozen pre-image shared-transition signal-recovery qualification.

This script never reads image pixels or measured colour. It verifies the exact
200-species metadata/geometry artifact, reconstructs the frozen 100-km
species-by-cell opportunity matrix, and qualifies only the cross-species scan layer.
The conditional rank-placement null used here is a pre-image method gate; it is not
the primary null for later real-colour inference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from fcp_pipeline.atlas_shared_transition_v5 import (
    build_coexceedance_reference,
    build_detectability_matrix,
    equal_area_cell_xyz,
    signal_recovery_rates,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = (
    REPO_ROOT
    / "docs"
    / "supporting"
    / "jbi_atlas_shared_transition_v5_signal_recovery_contract.json"
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_panels(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "taxon_id",
        "species",
        "genus",
        "selection_hash",
        "cohort_id",
        "cohort_species_index",
        "target_observations",
    }
    require(rows, "species panel is empty")
    require(required.issubset(rows[0]), "species panel columns changed")
    rows.sort(key=lambda r: (r["cohort_id"], int(r["cohort_species_index"])))
    return rows


def exact_opportunity_summary(D: np.ndarray) -> dict[str, Any]:
    opportunity = D.sum(axis=0).astype(int)
    detectable_per_species = D.sum(axis=1).astype(int)
    thresholds = (1, 2, 3, 4, 5, 10, 20)
    return {
        "detectable_cells_at_least": {
            str(t): int(np.count_nonzero(opportunity >= t)) for t in thresholds
        },
        "maximum_detectable_species_in_cell": int(opportunity.max()),
        "minimum_detectable_cells_per_species": int(detectable_per_species.min()),
        "median_detectable_cells_per_species": float(np.median(detectable_per_species)),
    }


def verify_geometry(
    metadata_dir: Path, contract: dict[str, Any]
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, Any]]:
    artifact = contract["terminal_geometry_artifact"]
    expected_files = artifact["files"]

    manifest_path = metadata_dir / artifact["metadata_manifest_required"]
    require(manifest_path.exists(), f"missing {manifest_path.name}")
    manifest = load_json(manifest_path)
    require(
        manifest.get("candidate_image_pixels_opened") is False,
        "metadata manifest says candidate pixels were opened",
    )

    for filename, expected_sha in expected_files.items():
        path = metadata_dir / filename
        require(path.exists(), f"missing frozen artifact file: {filename}")
        actual_sha = sha256_file(path)
        require(
            actual_sha == expected_sha,
            f"SHA256 mismatch for {filename}: {actual_sha} != {expected_sha}",
        )
        manifest_sha = manifest.get("files", {}).get(filename)
        require(
            manifest_sha == expected_sha,
            f"metadata manifest SHA mismatch for {filename}",
        )

    feasibility = load_json(metadata_dir / "scaleout_metadata_feasibility.json")
    require(
        feasibility.get("candidate_image_pixels_opened") is False,
        "feasibility receipt says candidate pixels were opened",
    )
    require(
        feasibility.get("continuous_colour_used") is False,
        "feasibility receipt used continuous colour",
    )
    require(
        feasibility.get("flower_roi_used") is False,
        "feasibility receipt used flower ROI",
    )
    require(
        int(feasibility.get("frozen_species", -1)) == int(artifact["frozen_species"]),
        "frozen species count changed",
    )
    require(
        int(feasibility.get("frozen_observations", -1))
        == int(artifact["frozen_observations"]),
        "frozen observation count changed",
    )

    panels = read_panels(metadata_dir / "scaleout_species_panels.csv")
    require(len(panels) == int(artifact["frozen_species"]), "panel row count changed")
    require(
        len({row["taxon_id"] for row in panels}) == len(panels),
        "panel taxon IDs are not unique",
    )
    require(
        len({row["species"] for row in panels}) == len(panels),
        "panel species names are not unique",
    )
    require(
        all(
            int(row["target_observations"])
            == int(artifact["observations_per_species"])
            for row in panels
        ),
        "target observations per species changed",
    )

    expected_cohorts = [f"C{i:02d}" for i in range(1, int(artifact["cohorts"]) + 1)]
    cohort_counts = {cohort: 0 for cohort in expected_cohorts}
    for row in panels:
        require(row["cohort_id"] in cohort_counts, f"unexpected cohort: {row['cohort_id']}")
        cohort_counts[row["cohort_id"]] += 1
    require(
        all(v == int(artifact["species_per_cohort"]) for v in cohort_counts.values()),
        f"cohort sizes changed: {cohort_counts}",
    )

    species_results = feasibility.get("species_results")
    require(isinstance(species_results, list), "species_results missing from feasibility receipt")
    by_taxon = {str(row["taxon_id"]): row for row in species_results}
    require(len(by_taxon) == len(species_results), "duplicate taxon_id in feasibility species_results")

    scale_km = int(contract["scan_layer"]["primary_scale_km"])
    species_cells: list[list[int]] = []
    cohort_species_cells: dict[str, list[list[int]]] = {c: [] for c in expected_cohorts}

    for panel in panels:
        result = by_taxon.get(panel["taxon_id"])
        require(result is not None, f"selected taxon missing from feasibility: {panel['taxon_id']}")
        require(result.get("species") == panel["species"], f"species name mismatch for {panel['taxon_id']}")
        require(bool(result.get("gate_pass")), f"selected species did not pass gate: {panel['species']}")
        scales = result.get("geometry_scale_results")
        require(isinstance(scales, list), f"missing geometry scales: {panel['species']}")
        matches = [s for s in scales if int(s.get("scale_km", -1)) == scale_km]
        require(len(matches) == 1, f"expected one {scale_km}-km geometry row: {panel['species']}")
        geometry = matches[0]
        require(bool(geometry.get("geometry_evaluable")), f"geometry not evaluable: {panel['species']}")
        cells = [int(x) for x in geometry.get("detectable_cell_ids", [])]
        require(
            len(cells) == int(geometry.get("detectable_cells", -1)),
            f"detectable-cell count mismatch: {panel['species']}",
        )
        species_cells.append(cells)
        cohort_species_cells[panel["cohort_id"]].append(cells)

    grid = contract["scan_layer"]["equal_area_grid"]
    n_cells = int(grid["n_cells"])
    D = build_detectability_matrix(species_cells, n_cells=n_cells)
    cohort_matrices = {
        cohort: build_detectability_matrix(cells, n_cells=n_cells)
        for cohort, cells in cohort_species_cells.items()
    }

    pooled_summary = exact_opportunity_summary(D)
    expected_pooled = contract["exact_geometry_expectations"]["pooled_200_species"]
    require(pooled_summary == expected_pooled, f"pooled opportunity changed: {pooled_summary}")

    per_cohort_summary: dict[str, Any] = {}
    expected_per_cohort = contract["exact_geometry_expectations"]["per_cohort"]
    for cohort, matrix in cohort_matrices.items():
        opportunity = matrix.sum(axis=0).astype(int)
        summary = {
            "cells_at_least_4_species": int(np.count_nonzero(opportunity >= 4)),
            "maximum_detectable_species_in_cell": int(opportunity.max()),
        }
        require(summary == expected_per_cohort[cohort], f"{cohort} opportunity changed: {summary}")
        per_cohort_summary[cohort] = summary

    return D, cohort_matrices, {
        "pooled_200_species": pooled_summary,
        "per_cohort": per_cohort_summary,
    }


def flatten_repetitions(result: dict[str, Any], alpha: float) -> list[dict[str, Any]]:
    series = [result["null_result"], result["heterogeneous_result"], *result["shared_results"]]
    rows: list[dict[str, Any]] = []
    for item in series:
        p_values = np.asarray(item["p_values"], dtype=float)
        statistics = np.asarray(item["statistics"], dtype=float)
        for i, (p_value, statistic) in enumerate(
            zip(p_values, statistics, strict=True), start=1
        ):
            rows.append(
                {
                    "scenario": item["scenario"],
                    "effect_size": float(item["effect_size"]),
                    "repetition": i,
                    "statistic": float(statistic),
                    "p_value": float(p_value),
                    "reject": bool(p_value <= alpha),
                }
            )
    return rows


def summarize_result(
    result: dict[str, Any], contract: dict[str, Any], opportunity: dict[str, Any]
) -> dict[str, Any]:
    qual = contract["preimage_signal_recovery"]
    gates = qual["gates"]
    null_rate = float(result["null_result"]["rate"])
    heterogeneous_rate = float(result["heterogeneous_result"]["rate"])
    shared_power = {
        str(float(item["effect_size"])): float(item["rate"])
        for item in result["shared_results"]
    }
    ordered_power = [shared_power[str(x)] for x in (0.5, 1.0, 2.0)]
    gate_results = {
        "null_false_positive_rate": null_rate
        <= float(gates["maximum_null_false_positive_rate"]),
        "heterogeneous_false_sharing": heterogeneous_rate
        <= float(gates["maximum_heterogeneous_false_sharing_rate_at_effect_2"]),
        "shared_effect_2_power": shared_power["2.0"]
        >= float(gates["minimum_shared_boundary_power_at_effect_2"]),
        "monotone_shared_power": ordered_power == sorted(ordered_power),
    }
    if not bool(gates["require_monotone_shared_power_across_effects_0_5_1_2"]):
        gate_results["monotone_shared_power"] = True

    null_distribution = np.asarray(result["null_distribution"], dtype=float)
    return {
        "contract": contract["contract"],
        "status": "pass" if all(gate_results.values()) else "fail",
        "pixel_status": "not_revealed",
        "method_gate_only": True,
        "biological_support_claimed": False,
        "geometry_opportunity": opportunity,
        "qualification": {
            "seed": int(qual["seed"]),
            "repetitions_per_scenario": int(qual["repetitions_per_scenario"]),
            "null_permutations": int(qual["null_permutations"]),
            "alpha": float(qual["alpha"]),
            "boundary_sigma_radians": float(qual["boundary_sigma_radians"]),
            "null_false_positive_rate": null_rate,
            "heterogeneous_false_sharing_rate_effect_2": heterogeneous_rate,
            "shared_boundary_power": shared_power,
            "null_scan_distribution": {
                "median": float(np.median(null_distribution)),
                "q95": float(np.quantile(null_distribution, 0.95)),
                "q99": float(np.quantile(null_distribution, 0.99)),
            },
            "gates": gate_results,
        },
        "firewall": {
            "qualification_null_only": True,
            "real_colour_primary_null": contract["real_colour_inference_firewall"][
                "real_colour_primary_null"
            ],
            "not_evaluable_never_means_unsupported": True,
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["scenario", "effect_size", "repetition", "statistic", "p_value", "reject"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    contract = load_json(args.contract)
    require(contract.get("pixel_status_at_freeze") == "not_revealed", "contract pixel firewall changed")
    require(
        contract["preimage_signal_recovery"].get("qualification_null_only") is True,
        "pre-image qualification null must remain qualification-only",
    )
    require(
        contract["real_colour_inference_firewall"].get(
            "qualification_null_must_not_be_used_as_real_colour_primary_null"
        )
        is True,
        "real-colour null firewall changed",
    )

    D, _cohort_matrices, opportunity = verify_geometry(args.metadata_dir, contract)
    scan = contract["scan_layer"]
    reference = build_coexceedance_reference(
        D,
        high_transition_quantile=float(scan["high_transition_quantile"]),
        min_detectable_species=int(scan["minimum_detectable_species_per_tested_cell"]),
    )
    grid = scan["equal_area_grid"]
    xyz = equal_area_cell_xyz(
        n_lon=int(grid["n_lon"]),
        n_sinlat=int(grid["n_sinlat"]),
    )
    qual = contract["preimage_signal_recovery"]
    result = signal_recovery_rates(
        reference,
        xyz,
        n_repetitions=int(qual["repetitions_per_scenario"]),
        n_permutations=int(qual["null_permutations"]),
        alpha=float(qual["alpha"]),
        effect_sizes=tuple(float(x) for x in qual["effect_sizes"]),
        boundary_sigma_radians=float(qual["boundary_sigma_radians"]),
        seed=int(qual["seed"]),
    )
    rows = flatten_repetitions(result, float(qual["alpha"]))
    summary = summarize_result(result, contract, opportunity)

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_csv(args.output_csv, rows)

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
