#!/usr/bin/env python3
"""Prepare frozen Florence ROI rows for continuous colour measurement.

Two immutable sources are combined only by ``photo_id``:

* the original 720-photo evaluation Florence run;
* a temporary 720-photo pseudo-split Florence run used solely to recover boxes for the
  480 original calibration photographs.

Original split membership is restored from the committed split CSV.  The 240 padding
rows in the pseudo-split are discarded by frozen photo ID, never by image/model output.
No colour value, coordinate, boundary statistic, or biological result is used to retain
or remove a photograph.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

EXPECTED_SPECIES = 6
EVALUATION_PER_SPECIES = 120
CALIBRATION_PER_SPECIES = 80
EXPECTED_EVALUATION = EXPECTED_SPECIES * EVALUATION_PER_SPECIES
EXPECTED_CALIBRATION = EXPECTED_SPECIES * CALIBRATION_PER_SPECIES
EXPECTED_PSEUDO_TOTAL = EXPECTED_SPECIES * EVALUATION_PER_SPECIES
EXPECTED_FILES_PER_SOURCE = EXPECTED_SPECIES * 6


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_id(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            pass
    return text


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
                    raise ValueError(f"{path}:{line_number}: row is not a JSON object")
                item = dict(row)
                item["_source_artifact_file"] = str(path.relative_to(root))
                rows.append(item)
    return rows, paths


def load_frozen_split(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    mapping: dict[str, dict[str, str]] = {}
    species: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("frozen split has no CSV header")
        required = {"species", "photo_id", "split", "photo_url", "latitude", "longitude"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"frozen split is missing {sorted(missing)}")
        for row in reader:
            photo_id = scalar_id(row.get("photo_id"))
            if not photo_id:
                raise ValueError("frozen split contains an empty photo_id")
            if photo_id in mapping:
                raise ValueError(f"frozen split contains duplicate photo_id {photo_id}")
            split = str(row.get("split", "")).strip().lower()
            if split not in {"calibration", "evaluation"}:
                raise ValueError(f"photo_id {photo_id}: invalid frozen split {split!r}")
            item = dict(row)
            item["photo_id"] = photo_id
            item["split"] = split
            mapping[photo_id] = item
            species.add(str(row.get("species", "")).strip())
    if len(species) != EXPECTED_SPECIES:
        raise ValueError(f"expected {EXPECTED_SPECIES} frozen species, found {len(species)}")
    return mapping, sorted(species)


def validate_source_rows(
    rows: list[dict[str, Any]],
    *,
    source_label: str,
    expected_total: int,
) -> dict[str, dict[str, Any]]:
    if len(rows) != expected_total:
        raise ValueError(f"{source_label}: expected {expected_total} rows, found {len(rows)}")
    by_photo: dict[str, dict[str, Any]] = {}
    by_blind: set[str] = set()
    failure_rows: list[tuple[str, str, str]] = []
    for row in rows:
        photo_id = scalar_id(row.get("photo_id"))
        blind_id = scalar_id(row.get("blind_id"))
        if not photo_id or not blind_id:
            raise ValueError(f"{source_label}: every row requires photo_id and blind_id")
        if photo_id in by_photo:
            raise ValueError(f"{source_label}: duplicate photo_id {photo_id}")
        if blind_id in by_blind:
            raise ValueError(f"{source_label}: duplicate blind_id {blind_id}")
        by_photo[photo_id] = row
        by_blind.add(blind_id)
        status = str(row.get("feature_status", "")).strip().lower()
        if status != "ok":
            failure_rows.append((photo_id, status or "missing", str(row.get("feature_error", ""))))
    if failure_rows:
        preview = "; ".join(f"{p}:{s}:{e}" for p, s, e in failure_rows[:10])
        raise ValueError(f"{source_label}: {len(failure_rows)} non-ok Florence rows: {preview}")
    return by_photo


def enrich(source: dict[str, Any], frozen: dict[str, str], source_role: str) -> dict[str, Any]:
    result = dict(source)
    result["photo_id"] = scalar_id(source.get("photo_id"))
    result["species"] = str(frozen["species"]).strip()
    result["split"] = str(frozen["split"]).strip().lower()
    result["frozen_split"] = result["split"]
    result["florence_source_role"] = source_role
    # Frozen metadata wins for identity, provenance, image URL and geometry.  Florence
    # model output remains untouched in all other fields.
    for key, value in frozen.items():
        if key in {
            "species",
            "photo_id",
            "split",
            "photo_url",
            "photo_url_api",
            "photo_license",
            "attribution",
            "observation_id",
            "latitude",
            "longitude",
            "positional_accuracy_m",
            "observed_on",
            "observed_month",
            "spatial_cell",
            "selection_hash",
            "split_rank_hash",
        }:
            result[key] = value
    result["evaluation_row"] = result["split"] == "evaluation"
    result["calibration_only"] = result["split"] == "calibration"
    result["final_label"] = False
    return result


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation-artifacts", type=Path, required=True)
    parser.add_argument("--calibration-pseudo-artifacts", type=Path, required=True)
    parser.add_argument("--frozen-split", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--evaluation-run-id", required=True)
    parser.add_argument("--calibration-run-id", required=True)
    args = parser.parse_args()

    frozen_by_photo, species_order = load_frozen_split(args.frozen_split)
    frozen_eval_ids = {
        photo_id for photo_id, row in frozen_by_photo.items() if row["split"] == "evaluation"
    }
    frozen_calibration_ids = {
        photo_id for photo_id, row in frozen_by_photo.items() if row["split"] == "calibration"
    }
    if len(frozen_eval_ids) != EXPECTED_EVALUATION:
        raise ValueError(
            f"frozen split has {len(frozen_eval_ids)} evaluation rows, expected {EXPECTED_EVALUATION}"
        )
    if len(frozen_calibration_ids) != EXPECTED_CALIBRATION:
        raise ValueError(
            f"frozen split has {len(frozen_calibration_ids)} calibration rows, expected {EXPECTED_CALIBRATION}"
        )

    evaluation_rows, evaluation_files = read_jsonl_tree(args.evaluation_artifacts)
    pseudo_rows, pseudo_files = read_jsonl_tree(args.calibration_pseudo_artifacts)
    if len(evaluation_files) != EXPECTED_FILES_PER_SOURCE:
        raise ValueError(
            f"evaluation source: expected {EXPECTED_FILES_PER_SOURCE} JSONL files, "
            f"found {len(evaluation_files)}"
        )
    if len(pseudo_files) != EXPECTED_FILES_PER_SOURCE:
        raise ValueError(
            f"calibration pseudo source: expected {EXPECTED_FILES_PER_SOURCE} JSONL files, "
            f"found {len(pseudo_files)}"
        )
    evaluation_by_photo = validate_source_rows(
        evaluation_rows,
        source_label="original evaluation Florence source",
        expected_total=EXPECTED_EVALUATION,
    )
    pseudo_by_photo = validate_source_rows(
        pseudo_rows,
        source_label="calibration pseudo-split Florence source",
        expected_total=EXPECTED_PSEUDO_TOTAL,
    )

    if set(evaluation_by_photo) != frozen_eval_ids:
        missing = sorted(frozen_eval_ids - set(evaluation_by_photo))
        unexpected = sorted(set(evaluation_by_photo) - frozen_eval_ids)
        raise ValueError(
            "original evaluation Florence photo IDs do not reproduce the frozen evaluation set: "
            f"missing={len(missing)}, unexpected={len(unexpected)}"
        )
    recovered_calibration_ids = set(pseudo_by_photo) & frozen_calibration_ids
    if recovered_calibration_ids != frozen_calibration_ids:
        missing = sorted(frozen_calibration_ids - recovered_calibration_ids)
        raise ValueError(
            f"calibration pseudo-split is missing {len(missing)} original calibration photo IDs"
        )
    ignored_padding_ids = set(pseudo_by_photo) - frozen_calibration_ids
    if len(ignored_padding_ids) != EXPECTED_SPECIES * 40:
        raise ValueError(
            f"expected {EXPECTED_SPECIES * 40} pseudo-split padding rows, found {len(ignored_padding_ids)}"
        )

    prepared_evaluation = [
        enrich(evaluation_by_photo[photo_id], frozen_by_photo[photo_id], "original_evaluation_run")
        for photo_id in frozen_eval_ids
    ]
    prepared_calibration = [
        enrich(pseudo_by_photo[photo_id], frozen_by_photo[photo_id], "temporary_calibration_pseudosplit_run")
        for photo_id in frozen_calibration_ids
    ]
    prepared_evaluation.sort(key=lambda row: (str(row["species"]), scalar_id(row["photo_id"])))
    prepared_calibration.sort(key=lambda row: (str(row["species"]), scalar_id(row["photo_id"])))

    eval_counts = Counter(str(row["species"]) for row in prepared_evaluation)
    cal_counts = Counter(str(row["species"]) for row in prepared_calibration)
    if sorted(eval_counts) != species_order or set(eval_counts.values()) != {EVALUATION_PER_SPECIES}:
        raise ValueError(f"evaluation species counts invalid: {dict(eval_counts)}")
    if sorted(cal_counts) != species_order or set(cal_counts.values()) != {CALIBRATION_PER_SPECIES}:
        raise ValueError(f"calibration species counts invalid: {dict(cal_counts)}")

    evaluation_path = args.output_dir / "evaluation_florence_boxes.jsonl"
    calibration_path = args.output_dir / "calibration_florence_boxes.jsonl"
    write_jsonl(evaluation_path, prepared_evaluation)
    write_jsonl(calibration_path, prepared_calibration)
    manifest = {
        "protocol": "jbi_ch1_prepare_frozen_florence_box_rows_v1",
        "validation_status": "success",
        "evaluation_run_id": str(args.evaluation_run_id),
        "calibration_pseudo_run_id": str(args.calibration_run_id),
        "evaluation_artifact_files": len(evaluation_files),
        "calibration_pseudo_artifact_files": len(pseudo_files),
        "evaluation_rows": len(prepared_evaluation),
        "calibration_rows": len(prepared_calibration),
        "ignored_padding_rows": len(ignored_padding_ids),
        "per_species_evaluation": dict(sorted(eval_counts.items())),
        "per_species_calibration": dict(sorted(cal_counts.items())),
        "species": species_order,
        "selection_of_calibration_rows": "frozen_original_split_membership_by_photo_id",
        "uses_florence_output_to_select_rows": False,
        "uses_colour_to_select_rows": False,
        "uses_coordinates_to_select_rows": False,
        "uses_outcomes_to_select_rows": False,
        "final_labels_created": False,
        "frozen_split_sha256": sha256_file(args.frozen_split),
        "evaluation_boxes_sha256": sha256_file(evaluation_path),
        "calibration_boxes_sha256": sha256_file(calibration_path),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = args.output_dir / "prepare_boxes_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"prepared {len(prepared_calibration)} calibration and {len(prepared_evaluation)} "
        "evaluation Florence box rows"
    )


if __name__ == "__main__":
    main()
