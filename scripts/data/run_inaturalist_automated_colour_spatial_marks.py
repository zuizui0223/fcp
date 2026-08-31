#!/usr/bin/env python3
"""Run the frozen locked-partition iNaturalist spatial random-mark tests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scripts.data.extract_inaturalist_automated_colour_states import (
        DEFAULT_SEED,
        PROTOCOL_VERSION,
        sha256,
    )
except ModuleNotFoundError:
    from extract_inaturalist_automated_colour_states import (  # type: ignore[no-redef]
        DEFAULT_SEED,
        PROTOCOL_VERSION,
        sha256,
    )


LOCKED_ENCOUNTERS_PER_SPECIES = 120
MIN_LOCKED_ADMITTED_SHARE = 0.70
MIN_LOCKED_BACKGROUND_SHARE = 0.70
PERMUTATIONS = 9_999
FDR_Q = 0.05
DISTANCE_BINS = 5
EARTH_RADIUS_KM = 6_371.0088


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def average_ranks(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0 + 1.0
        start = end
    return ranks


def spearman(first: np.ndarray, second: np.ndarray) -> float:
    first_ranks = average_ranks(first)
    second_ranks = average_ranks(second)
    first_centered = first_ranks - first_ranks.mean()
    second_centered = second_ranks - second_ranks.mean()
    denominator = float(
        np.linalg.norm(first_centered) * np.linalg.norm(second_centered)
    )
    if denominator <= 0 or not math.isfinite(denominator):
        return math.nan
    return float(np.dot(first_centered, second_centered) / denominator)


def robust_standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    median = np.median(values, axis=0)
    q25, q75 = np.quantile(values, [0.25, 0.75], axis=0)
    scale = q75 - q25
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError("non-finite values or zero IQR in locked marks")
    return (values - median) / scale


def pair_indices(n_rows: int) -> tuple[np.ndarray, np.ndarray]:
    return np.triu_indices(n_rows, k=1)


def geographic_pair_distances(
    latitude: np.ndarray, longitude: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    latitude_radians = np.radians(np.asarray(latitude, dtype=float))
    longitude_radians = np.radians(np.asarray(longitude, dtype=float))
    left, right = pair_indices(len(latitude_radians))
    delta_latitude = latitude_radians[right] - latitude_radians[left]
    delta_longitude = longitude_radians[right] - longitude_radians[left]
    a = (
        np.sin(delta_latitude / 2.0) ** 2
        + np.cos(latitude_radians[left])
        * np.cos(latitude_radians[right])
        * np.sin(delta_longitude / 2.0) ** 2
    )
    distance = 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))
    return distance, left, right


def mark_pair_distances(
    marks: np.ndarray, left: np.ndarray, right: np.ndarray, order: np.ndarray | None = None
) -> np.ndarray:
    indices = np.arange(len(marks)) if order is None else order
    return np.linalg.norm(marks[indices[left]] - marks[indices[right]], axis=1)


def numeric_seed(*tokens: str) -> int:
    token = "\x1f".join((DEFAULT_SEED, *tokens))
    return int(hashlib.sha256(token.encode("utf-8")).hexdigest()[:16], 16)


def random_mark_test(
    latitude: np.ndarray,
    longitude: np.ndarray,
    marks: np.ndarray,
    *,
    species: str,
    permutations: int = PERMUTATIONS,
) -> dict[str, Any]:
    standardized = robust_standardize(marks)
    geographic, left, right = geographic_pair_distances(latitude, longitude)
    observed = spearman(geographic, mark_pair_distances(standardized, left, right))
    if not math.isfinite(observed):
        return {"rho": math.nan, "p_greater": math.nan, "permutations": permutations}
    rng = np.random.default_rng(numeric_seed(species, "primary"))
    greater_or_equal = 0
    for _ in range(permutations):
        order = rng.permutation(len(standardized))
        null = spearman(
            geographic, mark_pair_distances(standardized, left, right, order)
        )
        greater_or_equal += null >= observed
    return {
        "rho": observed,
        "p_greater": (1 + greater_or_equal) / (1 + permutations),
        "permutations": permutations,
    }


def flower_background_contrast_test(
    latitude: np.ndarray,
    longitude: np.ndarray,
    flower_marks: np.ndarray,
    background_marks: np.ndarray,
    *,
    species: str,
    permutations: int = PERMUTATIONS,
) -> dict[str, Any]:
    flower = robust_standardize(flower_marks)
    background = robust_standardize(background_marks)
    geographic, left, right = geographic_pair_distances(latitude, longitude)
    flower_rho = spearman(geographic, mark_pair_distances(flower, left, right))
    background_rho = spearman(geographic, mark_pair_distances(background, left, right))
    observed = flower_rho - background_rho
    if not all(math.isfinite(value) for value in (flower_rho, background_rho, observed)):
        return {
            "flower_rho_common_subset": flower_rho,
            "background_rho": background_rho,
            "flower_minus_background_rho": observed,
            "p_greater": math.nan,
            "permutations": permutations,
        }
    rng = np.random.default_rng(numeric_seed(species, "flower-background"))
    greater_or_equal = 0
    for _ in range(permutations):
        order = rng.permutation(len(flower))
        flower_null = spearman(
            geographic, mark_pair_distances(flower, left, right, order)
        )
        background_null = spearman(
            geographic, mark_pair_distances(background, left, right, order)
        )
        greater_or_equal += flower_null - background_null >= observed
    return {
        "flower_rho_common_subset": flower_rho,
        "background_rho": background_rho,
        "flower_minus_background_rho": observed,
        "p_greater": (1 + greater_or_equal) / (1 + permutations),
        "permutations": permutations,
    }


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_index in range(len(values) - 1, -1, -1):
        original_index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, values[original_index] * len(values) / rank)
        adjusted[original_index] = min(running, 1.0)
    return adjusted.tolist()


def colour_fields(row: dict[str, str], prefix: str) -> np.ndarray:
    return np.array(
        [float(row[f"{prefix}_{channel}_mean"]) for channel in ("L", "a", "b")],
        dtype=float,
    )


def variogram_bins(
    species: str,
    latitude: np.ndarray,
    longitude: np.ndarray,
    marks: np.ndarray,
) -> list[dict[str, Any]]:
    standardized = robust_standardize(marks)
    geographic, left, right = geographic_pair_distances(latitude, longitude)
    colour = mark_pair_distances(standardized, left, right)
    order = np.argsort(geographic, kind="mergesort")
    rows: list[dict[str, Any]] = []
    for index, positions in enumerate(np.array_split(order, DISTANCE_BINS), start=1):
        rows.append(
            {
                "canonical_name": species,
                "equal_pair_count_bin": index,
                "n_pairs": len(positions),
                "distance_km_min": float(np.min(geographic[positions])),
                "distance_km_median": float(np.median(geographic[positions])),
                "distance_km_max": float(np.max(geographic[positions])),
                "standardized_colour_distance_median": float(np.median(colour[positions])),
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development-gate", type=Path, required=True)
    parser.add_argument("--locked-table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--public-manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    development_gate = json.loads(args.development_gate.read_text(encoding="utf-8"))
    if development_gate.get("protocol") != PROTOCOL_VERSION:
        raise RuntimeError("development gate protocol mismatch")
    if development_gate.get("status") != "complete_location_free_automated_colour_development_gate":
        raise RuntimeError("development gate is not complete")
    passed_species = sorted(
        row["canonical_name"]
        for row in development_gate["species_results"]
        if row["development_gate_status"] == "pass"
    )
    if not passed_species:
        raise RuntimeError("no species passed development; locked coordinates must remain unopened")

    locked_rows = read_csv(args.locked_table)
    results: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    evaluable_indices: list[int] = []
    primary_p: list[float] = []
    contrast_p: list[float] = []
    for species in passed_species:
        rows = [row for row in locked_rows if row["canonical_name"] == species]
        if len(rows) != LOCKED_ENCOUNTERS_PER_SPECIES:
            raise RuntimeError(f"{species}: locked table must contain exactly 120 encounters")
        if any(row["annotation_partition"] != "locked_60" for row in rows):
            raise RuntimeError(f"{species}: non-locked row in locked table")
        if len({row["encounter_blind_id"] for row in rows}) != len(rows):
            raise RuntimeError(f"{species}: duplicate locked encounter")
        admitted = [
            row
            for row in rows
            if row["encounter_status"] == "automated_colour_state_admitted"
        ]
        background = [
            row
            for row in admitted
            if row["background_control_status"] == "background_control_available"
        ]
        admission_share = len(admitted) / len(rows)
        background_share = len(background) / len(admitted) if admitted else 0.0
        result: dict[str, Any] = {
            "canonical_name": species,
            "locked_encounters": len(rows),
            "admitted_encounters": len(admitted),
            "admitted_encounter_share": admission_share,
            "background_control_encounters": len(background),
            "background_control_encounter_share": background_share,
            "primary_rho": math.nan,
            "primary_p_greater": math.nan,
            "primary_bh_q": math.nan,
            "leave_top_observer_rho": math.nan,
            "flower_minus_background_rho": math.nan,
            "contrast_p_greater": math.nan,
            "contrast_bh_q": math.nan,
            "spatial_claim_status": "not_evaluable",
        }
        if (
            admission_share < MIN_LOCKED_ADMITTED_SHARE
            or background_share < MIN_LOCKED_BACKGROUND_SHARE
        ):
            results.append(result)
            continue
        latitude = np.asarray([float(row["latitude"]) for row in admitted], dtype=float)
        longitude = np.asarray([float(row["longitude"]) for row in admitted], dtype=float)
        marks = np.vstack([colour_fields(row, "flower") for row in admitted])
        if (
            not np.all(np.isfinite(latitude))
            or not np.all(np.isfinite(longitude))
            or np.any((latitude < -90) | (latitude > 90))
            or np.any((longitude < -180) | (longitude > 180))
        ):
            result["spatial_claim_status"] = "not_evaluable_invalid_coordinate"
            results.append(result)
            continue
        try:
            primary = random_mark_test(latitude, longitude, marks, species=species)
        except ValueError:
            result["spatial_claim_status"] = "not_evaluable_nonfinite_or_zero_iqr"
            results.append(result)
            continue
        if not all(math.isfinite(primary[key]) for key in ("rho", "p_greater")):
            result["spatial_claim_status"] = "not_evaluable_constant_distance_or_mark"
            results.append(result)
            continue
        observer_counts = Counter(row["observer_id"] for row in admitted)
        top_observer = min(
            observer_counts,
            key=lambda observer: (-observer_counts[observer], str(observer)),
        )
        leave_rows = [row for row in admitted if row["observer_id"] != top_observer]
        try:
            leave = random_mark_test(
                np.asarray([float(row["latitude"]) for row in leave_rows]),
                np.asarray([float(row["longitude"]) for row in leave_rows]),
                np.vstack([colour_fields(row, "flower") for row in leave_rows]),
                species=f"{species}-leave-top-observer",
                permutations=0,
            )
        except ValueError:
            result["spatial_claim_status"] = "not_evaluable_leave_observer_zero_iqr"
            results.append(result)
            continue
        if not math.isfinite(leave["rho"]):
            result["spatial_claim_status"] = "not_evaluable_leave_observer_constant_distance"
            results.append(result)
            continue
        common_latitude = np.asarray([float(row["latitude"]) for row in background])
        common_longitude = np.asarray([float(row["longitude"]) for row in background])
        try:
            contrast = flower_background_contrast_test(
                common_latitude,
                common_longitude,
                np.vstack([colour_fields(row, "flower") for row in background]),
                np.vstack([colour_fields(row, "background") for row in background]),
                species=species,
            )
        except ValueError:
            result["spatial_claim_status"] = "not_evaluable_background_zero_iqr"
            results.append(result)
            continue
        if not all(
            math.isfinite(contrast[key])
            for key in ("flower_minus_background_rho", "background_rho", "p_greater")
        ):
            result["spatial_claim_status"] = "not_evaluable_constant_background_distance"
            results.append(result)
            continue
        result.update(
            {
                "primary_rho": primary["rho"],
                "primary_p_greater": primary["p_greater"],
                "leave_top_observer_n": observer_counts[top_observer],
                "leave_top_observer_rho": leave["rho"],
                "flower_minus_background_rho": contrast["flower_minus_background_rho"],
                "background_rho": contrast["background_rho"],
                "contrast_p_greater": contrast["p_greater"],
            }
        )
        evaluable_indices.append(len(results))
        primary_p.append(primary["p_greater"])
        contrast_p.append(contrast["p_greater"])
        curves.extend(variogram_bins(species, latitude, longitude, marks))
        results.append(result)

    primary_q = benjamini_hochberg(primary_p)
    contrast_q = benjamini_hochberg(contrast_p)
    for index, q_primary, q_contrast in zip(evaluable_indices, primary_q, contrast_q):
        result = results[index]
        result["primary_bh_q"] = q_primary
        result["contrast_bh_q"] = q_contrast
        supported = (
            result["primary_rho"] > 0
            and q_primary < FDR_Q
            and result["leave_top_observer_rho"] > 0
            and result["flower_minus_background_rho"] > 0
            and q_contrast < FDR_Q
        )
        result["spatial_claim_status"] = (
            "spatial_organization_supported" if supported else "spatial_organization_not_detected"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.output_dir / "species_spatial_mark_results.csv"
    curve_path = args.output_dir / "descriptive_equal_pair_variogram.csv"
    write_csv(result_path, results)
    if curves:
        write_csv(curve_path, curves)
    public = {
        "status": "complete_locked_species_conditioned_random_mark_test",
        "protocol": PROTOCOL_VERSION,
        "seed": DEFAULT_SEED,
        "permutations": PERMUTATIONS,
        "fdr_q": FDR_Q,
        "species_results": results,
        "universality_claim_allowed": False,
        "mechanism_claim_allowed": False,
        "claim_ceiling": (
            "Within selected development-passing species, model-consensus flower-candidate colour "
            "states did or did not reject a frozen locked species-conditioned random-mark null. "
            "The result does not verify flower tissue, botanical morphs, mechanisms or universality."
        ),
        "source_sha256": {
            "development_gate": sha256(args.development_gate),
            "locked_table": sha256(args.locked_table),
            "runner": sha256(Path(__file__).resolve()),
        },
        "output_sha256": {
            result_path.name: sha256(result_path),
            **({curve_path.name: sha256(curve_path)} if curve_path.exists() else {}),
        },
    }
    write_json(args.public_manifest, public)
    print(json.dumps(json_safe(public), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
