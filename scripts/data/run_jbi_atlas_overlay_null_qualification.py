#!/usr/bin/env python3
"""Qualify the spatially constrained atlas overlay concordance null."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_expansion import validate_expansion_contract
from fcp_pipeline.atlas_overlay_null import (
    MoranSignBasis,
    build_moran_sign_basis,
    equal_area_rook_adjacency,
    geographic_design,
    spectral_family_test,
)
from fcp_pipeline.shared_transition_surface import equal_area_cell_centers
from scripts.data.run_jbi_atlas_signal_recovery import build_geometries, read_csv


def numeric_seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


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
    parser.add_argument("--randomizations", type=int)
    return parser.parse_args()


def smooth_coefficients(basis: MoranSignBasis, rng: np.random.Generator) -> np.ndarray:
    eigenvalues = basis.eigenvalues
    span = float(np.ptp(eigenvalues))
    scaled = (eigenvalues - eigenvalues.max()) / span if span > 0 else eigenvalues * 0
    return rng.normal(size=eigenvalues.size) * np.exp(2.0 * scaled)


def unit(values: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(values))
    if norm <= 0:
        raise ValueError("synthetic coefficient vector has zero norm")
    return values / norm


def surface(
    coefficients: np.ndarray,
    basis: MoranSignBasis,
    *,
    trend: np.ndarray | None = None,
) -> np.ndarray:
    values = (basis.vectors @ coefficients) / basis.sqrt_weights
    return values if trend is None else values + trend


def run_repetition(
    basis: MoranSignBasis,
    design: np.ndarray,
    *,
    scenario: str,
    effect_size: float,
    family_size: int,
    randomizations: int,
    seed_parts: tuple[object, ...],
) -> tuple[float, float]:
    generator = np.random.default_rng(numeric_seed("overlay-fields", *seed_parts))
    randomizer = np.random.default_rng(numeric_seed("overlay-signs", *seed_parts))
    overlay_coefficients = [smooth_coefficients(basis, generator) for _ in range(family_size)]
    flower_coefficients = smooth_coefficients(basis, generator)
    shared_trend = 3.0 * design[:, 1] - 2.0 * design[:, 3] + design[:, 6]

    if scenario == "injected_local_boundary_alignment":
        flower_coefficients = unit(flower_coefficients) + effect_size * unit(
            overlay_coefficients[0]
        )
        trend = None
    elif scenario == "shared_broad_geographic_trend_only":
        trend = shared_trend
    elif scenario == "independent_smooth_fields":
        trend = None
    else:
        raise ValueError(f"unknown overlay-null scenario: {scenario}")

    flower = surface(flower_coefficients, basis, trend=trend)
    overlays = {
        f"overlay_{index + 1}": surface(coefficients, basis, trend=trend)
        for index, coefficients in enumerate(overlay_coefficients)
    }
    result = spectral_family_test(
        flower,
        overlays,
        basis,
        randomizations=randomizations,
        rng=randomizer,
    )
    return float(result["family_statistic"]), float(result["p_value"])


def main() -> None:
    args = parse_args()
    expansion = json.loads(args.expansion_contract.read_text(encoding="utf-8"))
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    validate_expansion_contract(expansion)
    qualification = expansion["overlay_null_qualification"]
    repetitions = args.repetitions or int(qualification["simulation_repetitions"])
    randomizations = args.randomizations or int(
        qualification["randomizations_per_repetition"]
    )
    formal = (
        repetitions == int(qualification["simulation_repetitions"])
        and randomizations == int(qualification["randomizations_per_repetition"])
    )
    observation_rows = read_csv(args.observations)
    rows: list[dict[str, Any]] = []
    jobs = [
        ("independent_smooth_fields", 0.0),
        ("shared_broad_geographic_trend_only", 0.0),
        *[
            ("injected_local_boundary_alignment", float(effect))
            for effect in qualification["effect_sizes"]
        ],
    ]
    for scale_km in qualification["scales_km"]:
        _geometries, grid, opportunity = build_geometries(
            observation_rows, atlas, int(scale_km)
        )
        minimum_species = int(
            atlas["geometry_only_scale_selection"]["passing_criteria"][
                "minimum_detectable_species_per_shared_cell"
            ]
        )
        valid = opportunity >= minimum_species
        cell_ids, latitude, longitude = equal_area_cell_centers(grid)
        admitted_ids = cell_ids[valid]
        design = geographic_design(latitude[valid], longitude[valid])
        adjacency = equal_area_rook_adjacency(
            admitted_ids, n_lon=grid.n_lon, n_sinlat=grid.n_sinlat
        )
        basis = build_moran_sign_basis(adjacency, opportunity[valid], design)
        for scenario, effect in jobs:
            for repetition in range(repetitions):
                statistic, p_value = run_repetition(
                    basis,
                    design,
                    scenario=scenario,
                    effect_size=effect,
                    family_size=int(qualification["overlay_family_size"]),
                    randomizations=randomizations,
                    seed_parts=(scale_km, scenario, effect, repetition),
                )
                rows.append(
                    {
                        "scale_km": int(scale_km),
                        "scenario": scenario,
                        "effect_size": effect,
                        "repetition": repetition + 1,
                        "statistic": statistic,
                        "p_value": p_value,
                        "detected": p_value <= float(qualification["alpha"]),
                    }
                )
            print(
                f"scale={scale_km} scenario={scenario} effect={effect:g} "
                f"repetitions={repetitions}",
                flush=True,
            )

    rates: dict[str, float] = {}
    checks: dict[str, bool] = {}
    for scale_km in qualification["scales_km"]:
        for scenario, effect in jobs:
            subset = [
                row
                for row in rows
                if int(row["scale_km"]) == int(scale_km)
                and row["scenario"] == scenario
                and float(row["effect_size"]) == effect
            ]
            key = f"{scale_km}|{scenario}|{effect:g}"
            rates[key] = sum(bool(row["detected"]) for row in subset) / len(subset)
        for scenario in (
            "independent_smooth_fields",
            "shared_broad_geographic_trend_only",
        ):
            key = f"{scale_km}|{scenario}|0"
            checks[f"false_positive_{key}"] = rates[key] <= float(
                qualification["maximum_false_positive_rate"]
            )
        power = [
            rates[f"{scale_km}|injected_local_boundary_alignment|{effect:g}"]
            for effect in qualification["effect_sizes"]
        ]
        checks[f"power_effect_2_{scale_km}"] = power[-1] >= float(
            qualification["minimum_alignment_power_at_effect_2"]
        )
        checks[f"monotone_power_{scale_km}"] = power == sorted(power)

    passed = formal and all(checks.values())
    result = {
        "status": (
            "pass_spatially_constrained_overlay_null"
            if passed
            else "smoke_only_not_qualification"
            if not formal
            else "stop_overlay_null_failed"
        ),
        "formal_contract_counts": formal,
        "repetitions": repetitions,
        "randomizations_per_repetition": randomizations,
        "detection_rates": rates,
        "checks": checks,
        "environmental_and_pollinator_colour_join_permitted": passed,
    }
    write_csv(args.output_dir / "overlay_null_repetitions.csv", rows)
    write_json(args.output_dir / "overlay_null_result.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if formal and not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
