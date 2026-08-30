#!/usr/bin/env python3
"""Finalize continuous CIELAB features for the frozen Chapter 1 split.

All 480 calibration photographs and all 720 evaluation photographs must be present and
successfully measured.  Species-specific component means and sample standard deviations
are estimated from the 80 calibration photographs only.  Those frozen scalers are then
applied unchanged to the 120 evaluation photographs of the same species.

The script fails closed: no row is silently removed, imputed, winsorized, relabelled, or
rescued with evaluation-derived statistics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

EXPECTED_SPECIES = 6
CALIBRATION_PER_SPECIES = 80
EVALUATION_PER_SPECIES = 120
EXPECTED_CALIBRATION = EXPECTED_SPECIES * CALIBRATION_PER_SPECIES
EXPECTED_EVALUATION = EXPECTED_SPECIES * EVALUATION_PER_SPECIES
EXPECTED_CALIBRATION_SHARDS = EXPECTED_SPECIES * 4
EXPECTED_EVALUATION_SHARDS = EXPECTED_SPECIES * 6
COMPONENT_NAMES = ("L_star", "a_star", "b_star")


def scalar_id(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            pass
    return text


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def read_jsonl_tree(root: Path) -> tuple[list[dict[str, Any]], list[Path]]:
    paths = sorted(path for path in root.rglob("*.jsonl") if path.is_file())
    rows: list[dict[str, Any]] = []
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(f"{path}:{line_number}: invalid JSON") from exc
                if not isinstance(row, dict):
                    raise ValueError(f"{path}:{line_number}: expected a JSON object")
                item = dict(row)
                item["_colour_source_artifact_file"] = str(path.relative_to(root))
                rows.append(item)
    return rows, paths


def load_frozen_split(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    mapping: dict[str, dict[str, str]] = {}
    species: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("frozen split has no header")
        required = {"species", "photo_id", "split", "latitude", "longitude"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"frozen split is missing {sorted(missing)}")
        for row in reader:
            photo_id = scalar_id(row.get("photo_id"))
            if not photo_id:
                raise ValueError("frozen split contains an empty photo_id")
            if photo_id in mapping:
                raise ValueError(f"frozen split contains duplicate photo_id {photo_id}")
            item = dict(row)
            item["photo_id"] = photo_id
            item["species"] = str(row.get("species", "")).strip()
            item["split"] = str(row.get("split", "")).strip().lower()
            if item["split"] not in {"calibration", "evaluation"}:
                raise ValueError(f"photo_id {photo_id}: invalid split {item['split']!r}")
            mapping[photo_id] = item
            species.add(item["species"])
    if len(species) != EXPECTED_SPECIES:
        raise ValueError(f"expected {EXPECTED_SPECIES} species, found {len(species)}")
    return mapping, sorted(species)


def parse_colour_vector(row: Mapping[str, Any]) -> np.ndarray:
    value = row.get("continuous_colour_vector")
    array = np.asarray(value, dtype=float)
    if array.shape != (len(COMPONENT_NAMES),):
        raise ValueError(
            f"continuous_colour_vector must have {len(COMPONENT_NAMES)} components, "
            f"found shape {array.shape}"
        )
    if not np.isfinite(array).all():
        raise ValueError("continuous_colour_vector contains a non-finite value")
    names = row.get("continuous_colour_component_names")
    if names is not None and tuple(str(name) for name in names) != COMPONENT_NAMES:
        raise ValueError(f"unexpected component names {names!r}")
    return array


def compute_species_scalers(
    species: Sequence[str],
    values: np.ndarray,
) -> dict[str, dict[str, Any]]:
    """Return component-wise calibration mean/SD for each species.

    This pure function accepts calibration data only.  Its signature deliberately has no
    evaluation argument, making accidental evaluation-derived scaling structurally
    impossible.
    """

    labels = np.asarray(species, dtype=object)
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != len(labels) or matrix.shape[1] < 2:
        raise ValueError("calibration species and value matrix have incompatible shapes")
    if not np.isfinite(matrix).all():
        raise ValueError("calibration matrix contains non-finite values")
    scalers: dict[str, dict[str, Any]] = {}
    for label in sorted(set(str(value) for value in labels)):
        subset = matrix[labels == label]
        if len(subset) < 2:
            raise ValueError(f"{label}: at least two calibration rows are required")
        mean = subset.mean(axis=0)
        sd = subset.std(axis=0, ddof=1)
        if not np.isfinite(mean).all() or not np.isfinite(sd).all() or np.any(sd <= 0):
            raise ValueError(f"{label}: calibration mean/SD is non-finite or degenerate")
        scalers[label] = {
            "n_calibration": int(len(subset)),
            "mean": mean.tolist(),
            "sd": sd.tolist(),
        }
    return scalers


def apply_species_scalers(
    species: Sequence[str],
    values: np.ndarray,
    scalers: Mapping[str, Mapping[str, Any]],
) -> np.ndarray:
    labels = np.asarray(species, dtype=object)
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != len(labels):
        raise ValueError("species and value matrix have incompatible shapes")
    result = np.empty_like(matrix, dtype=float)
    for label in sorted(set(str(value) for value in labels)):
        if label not in scalers:
            raise ValueError(f"no calibration scaler exists for {label}")
        mean = np.asarray(scalers[label]["mean"], dtype=float)
        sd = np.asarray(scalers[label]["sd"], dtype=float)
        if mean.shape != (matrix.shape[1],) or sd.shape != mean.shape or np.any(sd <= 0):
            raise ValueError(f"{label}: scaler shape/SD is invalid")
        idx = labels == label
        result[idx] = (matrix[idx] - mean) / sd
    if not np.isfinite(result).all():
        raise ValueError("calibration standardization produced non-finite values")
    return result


def validate_rows(
    rows: list[dict[str, Any]],
    *,
    expected_split: str,
    expected_total: int,
    expected_per_species: int,
    frozen_by_photo: Mapping[str, Mapping[str, str]],
    species_order: Sequence[str],
) -> tuple[list[dict[str, Any]], np.ndarray]:
    if len(rows) != expected_total:
        raise ValueError(f"{expected_split}: expected {expected_total} rows, found {len(rows)}")
    seen: set[str] = set()
    normalized: list[dict[str, Any]] = []
    vectors: list[np.ndarray] = []
    failures: list[str] = []
    counts: Counter[str] = Counter()
    for source in rows:
        photo_id = scalar_id(source.get("photo_id"))
        species = str(source.get("species", "")).strip()
        if not photo_id or not species:
            raise ValueError(f"{expected_split}: every row requires photo_id and species")
        if photo_id in seen:
            raise ValueError(f"{expected_split}: duplicate photo_id {photo_id}")
        seen.add(photo_id)
        frozen = frozen_by_photo.get(photo_id)
        if frozen is None:
            raise ValueError(f"{expected_split}: photo_id {photo_id} is absent from frozen split")
        if str(frozen["split"]) != expected_split:
            raise ValueError(
                f"photo_id {photo_id}: measured as {expected_split} but frozen split is {frozen['split']}"
            )
        if str(frozen["species"]) != species:
            raise ValueError(
                f"photo_id {photo_id}: measured species {species!r} != frozen {frozen['species']!r}"
            )
        status = str(source.get("colour_feature_status", "")).strip().lower()
        if status != "ok":
            failures.append(
                f"{photo_id}:{status or 'missing'}:{source.get('colour_feature_error', '')}"
            )
            continue
        vector = parse_colour_vector(source)
        item = dict(source)
        item["photo_id"] = photo_id
        item["species"] = species
        item["split"] = expected_split
        item["evaluation_row"] = expected_split == "evaluation"
        item["calibration_only"] = expected_split == "calibration"
        item["final_label"] = False
        for key, value in frozen.items():
            if key not in item or key in {
                "species",
                "photo_id",
                "split",
                "latitude",
                "longitude",
                "observation_id",
                "photo_url",
                "photo_license",
                "attribution",
                "observed_on",
                "observed_month",
                "spatial_cell",
                "selection_hash",
                "split_rank_hash",
            }:
                item[key] = value
        # Coordinates are materialized as finite numbers for the Stage-A graph builder.
        try:
            item["latitude"] = float(frozen["latitude"])
            item["longitude"] = float(frozen["longitude"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"photo_id {photo_id}: invalid frozen coordinates") from exc
        if not math.isfinite(item["latitude"]) or not math.isfinite(item["longitude"]):
            raise ValueError(f"photo_id {photo_id}: non-finite frozen coordinates")
        normalized.append(item)
        vectors.append(vector)
        counts[species] += 1
    if failures:
        preview = "; ".join(failures[:10])
        raise ValueError(
            f"{expected_split}: {len(failures)} colour measurements failed; no row was dropped: {preview}"
        )
    if sorted(counts) != list(species_order) or set(counts.values()) != {expected_per_species}:
        raise ValueError(f"{expected_split}: invalid per-species counts {dict(counts)}")
    frozen_ids = {
        photo_id for photo_id, row in frozen_by_photo.items() if row["split"] == expected_split
    }
    if seen != frozen_ids:
        raise ValueError(
            f"{expected_split}: measured photo IDs do not exactly equal frozen split IDs "
            f"(missing={len(frozen_ids - seen)}, unexpected={len(seen - frozen_ids)})"
        )
    return normalized, np.vstack(vectors)


def vector_summary(values: np.ndarray) -> dict[str, list[float]]:
    matrix = np.asarray(values, dtype=float)
    return {
        "min": np.min(matrix, axis=0).tolist(),
        "q025": np.quantile(matrix, 0.025, axis=0).tolist(),
        "median": np.quantile(matrix, 0.5, axis=0).tolist(),
        "mean": np.mean(matrix, axis=0).tolist(),
        "q975": np.quantile(matrix, 0.975, axis=0).tolist(),
        "max": np.max(matrix, axis=0).tolist(),
    }


def add_standardized_vectors(
    rows: list[dict[str, Any]],
    standardized: np.ndarray,
    scalers: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if len(rows) != len(standardized):
        raise ValueError("row and standardized matrices differ in length")
    output: list[dict[str, Any]] = []
    for source, vector in zip(rows, standardized, strict=True):
        item = dict(source)
        item["standardized_colour_vector"] = np.asarray(vector, dtype=float).tolist()
        item["standardized_colour_component_names"] = list(COMPONENT_NAMES)
        item["standardization_basis"] = "species_specific_calibration_mean_sample_sd"
        item["standardization_uses_evaluation_rows"] = False
        item["calibration_scaler_sha256"] = canonical_sha256(scalers[str(item["species"])])
        output.append(item)
    output.sort(key=lambda row: (str(row["species"]), scalar_id(row["photo_id"])))
    return output


def write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n")


def fail(message: str, manifest_path: Path, diagnostics: dict[str, Any]) -> None:
    diagnostics.update(
        {
            "validation_status": "failed",
            "validation_error": message,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    raise SystemExit(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration-artifacts", type=Path, required=True)
    parser.add_argument("--evaluation-artifacts", type=Path, required=True)
    parser.add_argument("--frozen-split", type=Path, required=True)
    parser.add_argument("--output-calibration-jsonl", type=Path, required=True)
    parser.add_argument("--output-evaluation-jsonl", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path, required=True)
    args = parser.parse_args()

    diagnostics: dict[str, Any] = {
        "protocol": "jbi_ch1_continuous_colour_features_v1",
        "expected_calibration_rows": EXPECTED_CALIBRATION,
        "expected_evaluation_rows": EXPECTED_EVALUATION,
        "expected_calibration_per_species": CALIBRATION_PER_SPECIES,
        "expected_evaluation_per_species": EVALUATION_PER_SPECIES,
        "component_names": list(COMPONENT_NAMES),
        "measurement_rule": "fixed_inner_ellipse_componentwise_10pct_trimmed_mean_cielab_d65",
        "standardization_rule": "species_specific_calibration_mean_sample_sd",
        "standardization_uses_evaluation_rows": False,
        "complete_case_required": True,
        "imputation": False,
        "outcome_blind": True,
        "final_labels_created": False,
    }
    try:
        frozen_by_photo, species_order = load_frozen_split(args.frozen_split)
        calibration_rows_raw, calibration_files = read_jsonl_tree(args.calibration_artifacts)
        evaluation_rows_raw, evaluation_files = read_jsonl_tree(args.evaluation_artifacts)
        diagnostics["calibration_artifact_files"] = len(calibration_files)
        diagnostics["evaluation_artifact_files"] = len(evaluation_files)
        if len(calibration_files) != EXPECTED_CALIBRATION_SHARDS:
            raise ValueError(
                f"expected {EXPECTED_CALIBRATION_SHARDS} calibration shard files, "
                f"found {len(calibration_files)}"
            )
        if len(evaluation_files) != EXPECTED_EVALUATION_SHARDS:
            raise ValueError(
                f"expected {EXPECTED_EVALUATION_SHARDS} evaluation shard files, "
                f"found {len(evaluation_files)}"
            )
        calibration_rows, calibration_values = validate_rows(
            calibration_rows_raw,
            expected_split="calibration",
            expected_total=EXPECTED_CALIBRATION,
            expected_per_species=CALIBRATION_PER_SPECIES,
            frozen_by_photo=frozen_by_photo,
            species_order=species_order,
        )
        evaluation_rows, evaluation_values = validate_rows(
            evaluation_rows_raw,
            expected_split="evaluation",
            expected_total=EXPECTED_EVALUATION,
            expected_per_species=EVALUATION_PER_SPECIES,
            frozen_by_photo=frozen_by_photo,
            species_order=species_order,
        )
        calibration_species = [str(row["species"]) for row in calibration_rows]
        evaluation_species = [str(row["species"]) for row in evaluation_rows]
        scalers = compute_species_scalers(calibration_species, calibration_values)
        if sorted(scalers) != species_order:
            raise ValueError("calibration scalers do not cover all frozen species")
        for species in species_order:
            if int(scalers[species]["n_calibration"]) != CALIBRATION_PER_SPECIES:
                raise ValueError(f"{species}: scaler did not use exactly 80 calibration rows")
        calibration_z = apply_species_scalers(calibration_species, calibration_values, scalers)
        evaluation_z = apply_species_scalers(evaluation_species, evaluation_values, scalers)
        output_calibration = add_standardized_vectors(calibration_rows, calibration_z, scalers)
        output_evaluation = add_standardized_vectors(evaluation_rows, evaluation_z, scalers)
        write_jsonl(args.output_calibration_jsonl, output_calibration)
        write_jsonl(args.output_evaluation_jsonl, output_evaluation)

        raw_qc: dict[str, Any] = {}
        z_qc: dict[str, Any] = {}
        for species in species_order:
            cal_idx = np.asarray(calibration_species, dtype=object) == species
            eval_idx = np.asarray(evaluation_species, dtype=object) == species
            raw_qc[species] = {
                "calibration": vector_summary(calibration_values[cal_idx]),
                "evaluation": vector_summary(evaluation_values[eval_idx]),
            }
            z_qc[species] = {
                "calibration": vector_summary(calibration_z[cal_idx]),
                "evaluation": vector_summary(evaluation_z[eval_idx]),
            }
        diagnostics.update(
            {
                "validation_status": "success",
                "calibration_rows": len(output_calibration),
                "evaluation_rows": len(output_evaluation),
                "species": species_order,
                "per_species_calibration": dict(Counter(calibration_species)),
                "per_species_evaluation": dict(Counter(evaluation_species)),
                "calibration_scalers": scalers,
                "calibration_scalers_sha256": canonical_sha256(scalers),
                "raw_feature_qc": raw_qc,
                "standardized_feature_qc": z_qc,
                "frozen_split_sha256": sha256_file(args.frozen_split),
                "calibration_jsonl_sha256": sha256_file(args.output_calibration_jsonl),
                "evaluation_jsonl_sha256": sha256_file(args.output_evaluation_jsonl),
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
        args.manifest_json.write_text(
            json.dumps(diagnostics, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"finalized {len(output_calibration)} calibration and {len(output_evaluation)} "
            "evaluation continuous-colour rows; evaluation never entered scaler estimation"
        )
    except Exception as exc:
        fail(str(exc), args.manifest_json, diagnostics)


if __name__ == "__main__":
    main()
