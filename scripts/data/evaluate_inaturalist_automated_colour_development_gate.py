#!/usr/bin/env python3
"""Evaluate the frozen, location-free automated colour development gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    from scripts.data.extract_inaturalist_automated_colour_states import (
        DEFAULT_SEED,
        PROTOCOL_VERSION,
        sha256,
        stable_rank,
    )
except ModuleNotFoundError:
    from extract_inaturalist_automated_colour_states import (  # type: ignore[no-redef]
        DEFAULT_SEED,
        PROTOCOL_VERSION,
        sha256,
        stable_rank,
    )


EXPECTED_SPECIES = 6
EXPECTED_ENCOUNTERS_PER_SPECIES = 80
EXPECTED_TOTAL_ENCOUNTERS = 480
MIN_ADMITTED_ENCOUNTER_SHARE = 0.70
MIN_REPEATABILITY_ENCOUNTERS = 10
REPEATABILITY_PERMUTATIONS = 9_999
REPEATABILITY_ALPHA = 0.05
TECHNICAL_FOLDS = 5
TECHNICAL_RIDGE_ALPHA = 1.0
MAX_TECHNICAL_CV_R2 = 0.80
MIN_BACKGROUND_CONTROL_SHARE = 0.70

TECHNICAL_FIELDS = (
    "luminance_mean",
    "fraction_luminance_le_5",
    "fraction_luminance_ge_250",
    "log_pixel_count",
    "log_aspect_ratio",
    "log1p_edge_variance",
)
FORBIDDEN_LOCATION_FIELDS = {
    "latitude",
    "longitude",
    "observed_on",
    "observation_date",
    "observer_id",
    "user_id",
    "place_guess",
}


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


def lab(row: dict[str, str], prefix: str = "flower") -> np.ndarray:
    return np.array(
        [float(row[f"{prefix}_{channel}_mean"]) for channel in ("L", "a", "b")],
        dtype=float,
    )


def repeatability_permutation(
    photo_rows: Iterable[dict[str, str]],
    species: str,
    *,
    permutations: int = REPEATABILITY_PERMUTATIONS,
    seed: str = DEFAULT_SEED,
) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in photo_rows:
        if (
            row["canonical_name"] == species
            and row["automated_colour_state_status"] == "automated_colour_state_admitted"
        ):
            grouped[row["encounter_blind_id"]].append(row)
    eligible = [
        (encounter_id, rows)
        for encounter_id, rows in sorted(grouped.items())
        if len(rows) >= 2
    ]
    if len(eligible) < MIN_REPEATABILITY_ENCOUNTERS:
        return {
            "repeatability_status": "not_evaluable_insufficient_multi_photo_encounters",
            "repeatability_encounters": len(eligible),
            "repeatability_pairs": 0,
            "within_encounter_delta_e_median": math.nan,
            "permutation_p_lower": math.nan,
            "repeatability_pass": False,
        }

    features: list[np.ndarray] = []
    pair_left: list[int] = []
    pair_right: list[int] = []
    offset = 0
    for _encounter_id, rows in eligible:
        vectors = [lab(row) for row in rows]
        features.extend(vectors)
        for left, right in itertools.combinations(range(offset, offset + len(vectors)), 2):
            pair_left.append(left)
            pair_right.append(right)
        offset += len(vectors)
    values = np.vstack(features)
    distances = np.linalg.norm(values[:, None, :] - values[None, :, :], axis=2)
    pair_left_array = np.asarray(pair_left, dtype=int)
    pair_right_array = np.asarray(pair_right, dtype=int)
    observed = float(np.median(distances[pair_left_array, pair_right_array]))
    numeric_seed = int(hashlib.sha256(f"{seed}\x1f{species}".encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(numeric_seed)
    lower_or_equal = 0
    for _ in range(permutations):
        order = rng.permutation(len(values))
        null_value = float(
            np.median(distances[order[pair_left_array], order[pair_right_array]])
        )
        lower_or_equal += null_value <= observed
    p_value = (1 + lower_or_equal) / (1 + permutations)
    passed = p_value < REPEATABILITY_ALPHA
    return {
        "repeatability_status": "pass" if passed else "fail",
        "repeatability_encounters": len(eligible),
        "repeatability_pairs": len(pair_left),
        "within_encounter_delta_e_median": observed,
        "permutation_p_lower": p_value,
        "repeatability_pass": passed,
    }


def fixed_fold_ids(encounter_ids: list[str], species: str) -> np.ndarray:
    ranked = sorted(encounter_ids, key=lambda value: stable_rank(DEFAULT_SEED, species, value))
    assignment = {encounter_id: index % TECHNICAL_FOLDS for index, encounter_id in enumerate(ranked)}
    return np.asarray([assignment[value] for value in encounter_ids], dtype=int)


def ridge_cv_r2(
    predictors: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    *,
    alpha: float = TECHNICAL_RIDGE_ALPHA,
) -> float:
    if predictors.ndim != 2 or outcomes.ndim != 2 or len(predictors) != len(outcomes):
        raise ValueError("predictor/outcome shape mismatch")
    predictions = np.full_like(outcomes, np.nan, dtype=float)
    baselines = np.full_like(outcomes, np.nan, dtype=float)
    for fold in range(TECHNICAL_FOLDS):
        test = folds == fold
        train = ~test
        if np.sum(test) == 0 or np.sum(train) <= predictors.shape[1]:
            return math.nan
        x_mean = predictors[train].mean(axis=0)
        x_sd = predictors[train].std(axis=0)
        x_sd[x_sd == 0] = 1.0
        y_mean = outcomes[train].mean(axis=0)
        y_sd = outcomes[train].std(axis=0)
        y_sd[y_sd == 0] = 1.0
        x_train = (predictors[train] - x_mean) / x_sd
        y_train = (outcomes[train] - y_mean) / y_sd
        penalty = alpha * np.eye(x_train.shape[1])
        coefficients = np.linalg.solve(x_train.T @ x_train + penalty, x_train.T @ y_train)
        standardized_prediction = ((predictors[test] - x_mean) / x_sd) @ coefficients
        predictions[test] = standardized_prediction * y_sd + y_mean
        baselines[test] = y_mean
    if not np.all(np.isfinite(predictions)):
        return math.nan
    residual = float(np.sum((outcomes - predictions) ** 2))
    baseline_error = float(np.sum((outcomes - baselines) ** 2))
    if baseline_error <= 0 or not math.isfinite(baseline_error):
        return math.nan
    return 1.0 - residual / baseline_error


def technical_vector(row: dict[str, str]) -> np.ndarray:
    width = float(row["image_width"])
    height = float(row["image_height"])
    return np.array(
        [
            float(row["luminance_mean"]),
            float(row["fraction_luminance_le_5"]),
            float(row["fraction_luminance_ge_250"]),
            math.log(max(width * height, 1.0)),
            math.log(max(float(row["aspect_ratio"]), 1e-12)),
            math.log1p(max(float(row["edge_variance_descriptive"]), 0.0)),
        ],
        dtype=float,
    )


def technical_dependence(
    photo_rows: list[dict[str, str]],
    encounter_rows: list[dict[str, str]],
    technical_rows: list[dict[str, str]],
    species: str,
) -> dict[str, Any]:
    profile = {row["photo_blind_id"]: row for row in technical_rows}
    admitted_photos: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in photo_rows:
        if (
            row["canonical_name"] == species
            and row["automated_colour_state_status"] == "automated_colour_state_admitted"
        ):
            admitted_photos[row["encounter_blind_id"]].append(row)
    admitted_encounters = [
        row
        for row in encounter_rows
        if row["canonical_name"] == species
        and row["encounter_status"] == "automated_colour_state_admitted"
    ]
    encounter_ids: list[str] = []
    predictors: list[np.ndarray] = []
    outcomes: list[np.ndarray] = []
    for row in admitted_encounters:
        encounter_id = row["encounter_blind_id"]
        attached = admitted_photos[encounter_id]
        if not attached or any(photo["photo_blind_id"] not in profile for photo in attached):
            return {
                "technical_status": "not_evaluable_missing_profile",
                "technical_cv_r2": math.nan,
                "technical_dependence_pass": False,
            }
        encounter_ids.append(encounter_id)
        predictors.append(
            np.median(
                np.vstack([technical_vector(profile[photo["photo_blind_id"]]) for photo in attached]),
                axis=0,
            )
        )
        outcomes.append(lab(row))
    if len(encounter_ids) < TECHNICAL_FOLDS * 5:
        return {
            "technical_status": "not_evaluable_insufficient_encounters",
            "technical_cv_r2": math.nan,
            "technical_dependence_pass": False,
        }
    x = np.vstack(predictors)
    y = np.vstack(outcomes)
    folds = fixed_fold_ids(encounter_ids, species)
    score = ridge_cv_r2(x, y, folds)
    passed = math.isfinite(score) and score < MAX_TECHNICAL_CV_R2
    return {
        "technical_status": "pass" if passed else "fail",
        "technical_predictors": ";".join(TECHNICAL_FIELDS),
        "technical_encounters": len(encounter_ids),
        "technical_cv_r2": score,
        "technical_dependence_pass": passed,
    }


def validate_complete_inputs(
    extraction_dir: Path,
    review_artifact: Path,
) -> tuple[dict[str, Any], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    run_manifest_path = extraction_dir / "run_manifest.json"
    photo_path = extraction_dir / "photo_features.csv"
    encounter_path = extraction_dir / "encounter_features.csv"
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    photo_rows = read_csv(photo_path)
    encounter_rows = read_csv(encounter_path)
    technical_rows = read_csv(review_artifact / "technical_image_profile.csv")
    if run_manifest.get("protocol") != PROTOCOL_VERSION:
        raise RuntimeError("extraction protocol does not match the frozen development protocol")
    if run_manifest.get("status") != "complete_automated_colour_state_development_feasibility_not_spatial":
        raise RuntimeError("extraction manifest is not complete")
    if run_manifest.get("limit_encounters_per_species") != 0:
        raise RuntimeError("development gate requires the complete development partition")
    if run_manifest.get("selected_encounters") != EXPECTED_TOTAL_ENCOUNTERS:
        raise RuntimeError("development extraction does not contain exactly 480 encounters")
    if run_manifest.get("selected_photos") != len(technical_rows) or len(photo_rows) != len(technical_rows):
        raise RuntimeError("development extraction photo count does not match the sealed packet")
    if len(encounter_rows) != EXPECTED_TOTAL_ENCOUNTERS:
        raise RuntimeError("encounter feature table does not contain exactly 480 rows")
    counts = Counter(row["canonical_name"] for row in encounter_rows)
    if len(counts) != EXPECTED_SPECIES or set(counts.values()) != {EXPECTED_ENCOUNTERS_PER_SPECIES}:
        raise RuntimeError("development extraction is not six species by 80 encounters")
    if len({row["photo_blind_id"] for row in photo_rows}) != len(photo_rows):
        raise RuntimeError("duplicate photo blind IDs")
    if len({row["photo_blind_id"] for row in technical_rows}) != len(technical_rows):
        raise RuntimeError("duplicate technical-profile photo blind IDs")
    if {row["photo_blind_id"] for row in photo_rows} != {
        row["photo_blind_id"] for row in technical_rows
    }:
        raise RuntimeError("feature and technical-profile photo ID sets differ")
    if len({row["encounter_blind_id"] for row in encounter_rows}) != len(encounter_rows):
        raise RuntimeError("duplicate encounter blind IDs")
    if FORBIDDEN_LOCATION_FIELDS.intersection(photo_rows[0]) or FORBIDDEN_LOCATION_FIELDS.intersection(encounter_rows[0]):
        raise RuntimeError("location or observer field leaked into development features")
    for row in photo_rows:
        if row["automated_colour_state_status"] == "automated_colour_state_admitted":
            if not np.all(np.isfinite(lab(row))):
                raise RuntimeError("admitted photo has non-finite flower mean Lab")
    for row in encounter_rows:
        if row["encounter_status"] == "automated_colour_state_admitted":
            if not np.all(np.isfinite(lab(row))):
                raise RuntimeError("admitted encounter has non-finite flower mean Lab")
        if row["background_control_status"] == "background_control_available":
            if not np.all(np.isfinite(lab(row, "background"))):
                raise RuntimeError("available encounter background control is non-finite")
    declared_hashes = run_manifest.get("private_output_sha256", {})
    for path in (photo_path, encounter_path):
        if declared_hashes.get(path.name) != sha256(path):
            raise RuntimeError(f"extraction hash mismatch: {path.name}")
    return run_manifest, photo_rows, encounter_rows, technical_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-dir", type=Path, required=True)
    parser.add_argument("--review-artifact", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--public-manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    extraction_dir = args.extraction_dir.resolve()
    review_artifact = args.review_artifact.resolve()
    output_dir = args.output_dir.resolve()
    run_manifest, photo_rows, encounter_rows, technical_rows = validate_complete_inputs(
        extraction_dir, review_artifact
    )
    species_results: list[dict[str, Any]] = []
    for species in sorted({row["canonical_name"] for row in encounter_rows}):
        species_encounters = [row for row in encounter_rows if row["canonical_name"] == species]
        admitted = [
            row
            for row in species_encounters
            if row["encounter_status"] == "automated_colour_state_admitted"
        ]
        background = [
            row
            for row in admitted
            if row["background_control_status"] == "background_control_available"
        ]
        admission_share = len(admitted) / len(species_encounters)
        background_share = len(background) / len(admitted) if admitted else 0.0
        repeatability = repeatability_permutation(photo_rows, species)
        technical = technical_dependence(photo_rows, encounter_rows, technical_rows, species)
        gates = {
            "admission_pass": admission_share >= MIN_ADMITTED_ENCOUNTER_SHARE,
            "repeatability_pass": repeatability["repeatability_pass"],
            "technical_dependence_pass": technical["technical_dependence_pass"],
            "background_control_pass": background_share >= MIN_BACKGROUND_CONTROL_SHARE,
        }
        passed = all(gates.values())
        species_results.append(
            {
                "canonical_name": species,
                "development_gate_status": "pass" if passed else "not_evaluable",
                "development_encounters": len(species_encounters),
                "admitted_encounters": len(admitted),
                "admitted_encounter_share": admission_share,
                "background_control_encounters": len(background),
                "background_control_encounter_share": background_share,
                **repeatability,
                **technical,
                **gates,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    species_path = output_dir / "species_development_gate.csv"
    write_csv(species_path, species_results)
    public_core = {
        "status": "complete_location_free_automated_colour_development_gate",
        "protocol": PROTOCOL_VERSION,
        "source_extraction_contract_sha256": run_manifest["contract_sha256"],
        "fixed_gates": {
            "min_admitted_encounter_share": MIN_ADMITTED_ENCOUNTER_SHARE,
            "min_repeatability_encounters": MIN_REPEATABILITY_ENCOUNTERS,
            "repeatability_permutations": REPEATABILITY_PERMUTATIONS,
            "repeatability_alpha": REPEATABILITY_ALPHA,
            "technical_folds": TECHNICAL_FOLDS,
            "technical_ridge_alpha": TECHNICAL_RIDGE_ALPHA,
            "max_technical_cv_r2": MAX_TECHNICAL_CV_R2,
            "min_background_control_share": MIN_BACKGROUND_CONTROL_SHARE,
        },
        "seed": DEFAULT_SEED,
        "species_results": species_results,
        "species_passed": sum(row["development_gate_status"] == "pass" for row in species_results),
        "species_not_evaluable": sum(
            row["development_gate_status"] != "pass" for row in species_results
        ),
        "spatial_colour_outcome_opened": False,
        "claim_ceiling": (
            "Location-free image-measurement development diagnostics only. A pass permits a later "
            "predeclared species-conditioned spatial random-mark test; it is not evidence of spatial "
            "organization, botanical flower tissue, named colours, discrete morphs or mechanism."
        ),
    }
    private_report = {
        **public_core,
        "source_paths": {
            "extraction_dir": str(extraction_dir),
            "review_artifact": str(review_artifact),
        },
        "source_sha256": {
            "run_manifest.json": sha256(extraction_dir / "run_manifest.json"),
            "photo_features.csv": sha256(extraction_dir / "photo_features.csv"),
            "encounter_features.csv": sha256(extraction_dir / "encounter_features.csv"),
            "technical_image_profile.csv": sha256(review_artifact / "technical_image_profile.csv"),
        },
        "output_sha256": {species_path.name: sha256(species_path)},
    }
    private_path = output_dir / "development_gate_report.json"
    write_json(private_path, private_report)
    public = {
        **public_core,
        "private_report_sha256": sha256(private_path),
        "species_table_sha256": sha256(species_path),
    }
    write_json(args.public_manifest, public)
    print(json.dumps(json_safe(public), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
