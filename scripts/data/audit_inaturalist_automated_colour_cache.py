#!/usr/bin/env python3
"""Audit resumable automated-colour cache records without opening geography."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

try:
    from scripts.data.extract_inaturalist_automated_colour_states import (
        MODEL_ID,
        MODEL_REVISION,
        PROTOCOL_VERSION,
        contract_sha256,
        sha256,
    )
except ModuleNotFoundError:
    from extract_inaturalist_automated_colour_states import (  # type: ignore[no-redef]
        MODEL_ID,
        MODEL_REVISION,
        PROTOCOL_VERSION,
        contract_sha256,
        sha256,
    )


EXPECTED_TOTAL_PHOTOS = 886
ALLOWED_STATUSES = {
    "automated_colour_state_admitted",
    "automated_colour_state_not_evaluable",
}
FORBIDDEN_FIELDS = {
    "latitude",
    "longitude",
    "observed_on",
    "observer_id",
    "user_id",
    "place_guess",
    "annotation_partition",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def expected_photos(review_artifact: Path) -> dict[str, dict[str, str]]:
    expected: dict[str, dict[str, str]] = {}
    for encounter in read_csv(review_artifact / "reviewer_A_annotation_sheet.csv"):
        for image_file in encounter["image_files"].split("|"):
            photo_id = Path(image_file).stem
            if photo_id in expected:
                raise RuntimeError(f"duplicate expected photo blind ID: {photo_id}")
            expected[photo_id] = {
                "canonical_name": encounter["canonical_name"],
                "encounter_blind_id": encounter["encounter_blind_id"],
                "image_file": image_file,
            }
    return expected


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def audit_cache(review_artifact: Path, cache_dir: Path, model_dir: Path) -> dict[str, Any]:
    expected = expected_photos(review_artifact)
    model_hash = sha256(model_dir / "model.safetensors")
    expected_contract = contract_sha256(model_hash)
    cache_files = {path.stem: path for path in cache_dir.glob("*.json")}
    errors: list[str] = []
    if len(expected) != EXPECTED_TOTAL_PHOTOS:
        errors.append(f"expected_packet_photo_count:{len(expected)}")
    unexpected = sorted(set(cache_files) - set(expected))
    if unexpected:
        errors.append(f"unexpected_cache_records:{len(unexpected)}")
    status_counts: Counter[str] = Counter()
    species_counts: Counter[str] = Counter()
    valid_records = 0
    for photo_id in sorted(set(cache_files) & set(expected)):
        path = cache_files[photo_id]
        expected_row = expected[photo_id]
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            errors.append(f"unreadable_json:{photo_id}")
            continue
        record_errors: list[str] = []
        if FORBIDDEN_FIELDS.intersection(record):
            record_errors.append("forbidden_field")
        for field in ("canonical_name", "encounter_blind_id", "image_file"):
            if record.get(field) != expected_row[field]:
                record_errors.append(f"{field}_mismatch")
        if record.get("photo_blind_id") != photo_id:
            record_errors.append("photo_blind_id_mismatch")
        if record.get("model_id") != MODEL_ID or record.get("model_revision") != MODEL_REVISION:
            record_errors.append("model_identity_mismatch")
        if record.get("contract_sha256") != expected_contract:
            record_errors.append("contract_mismatch")
        image_path = review_artifact / expected_row["image_file"]
        if record.get("image_sha256") != sha256(image_path):
            record_errors.append("image_hash_mismatch")
        status = record.get("automated_colour_state_status")
        if status not in ALLOWED_STATUSES:
            record_errors.append("invalid_status")
        if status == "automated_colour_state_admitted":
            for channel in ("L", "a", "b"):
                if not finite_number(record.get(f"flower_{channel}_mean")):
                    record_errors.append(f"nonfinite_flower_{channel}")
        background_available = record.get("background_features_available")
        if not isinstance(background_available, bool):
            record_errors.append("invalid_background_availability")
        elif background_available:
            for channel in ("L", "a", "b"):
                if not finite_number(record.get(f"background_{channel}_mean")):
                    record_errors.append(f"nonfinite_background_{channel}")
        if record_errors:
            errors.append(f"{photo_id}:{','.join(sorted(set(record_errors)))}")
            continue
        valid_records += 1
        status_counts[status] += 1
        species_counts[record["canonical_name"]] += 1
    partial_files = [
        path for path in cache_dir.parent.rglob("*") if path.is_file() and path.suffix in {".part", ".partial", ".tmp"}
    ]
    if partial_files:
        errors.append(f"partial_files:{len(partial_files)}")
    complete = valid_records == len(expected)
    status = (
        "invalid_cache"
        if errors
        else "complete_validated_cache"
        if complete
        else "running_validated_partial_cache"
    )
    return {
        "status": status,
        "protocol": PROTOCOL_VERSION,
        "contract_sha256": expected_contract,
        "expected_photos": len(expected),
        "cache_files_present": len(cache_files),
        "valid_cache_records": valid_records,
        "missing_cache_records": len(expected) - valid_records,
        "unexpected_cache_records": len(unexpected),
        "partial_files": len(partial_files),
        "errors": errors,
        "status_counts": dict(sorted(status_counts.items())),
        "valid_records_by_species": dict(sorted(species_counts.items())),
        "spatial_fields_opened": False,
        "biological_conclusion_allowed": False,
        "claim_ceiling": "Cache integrity/progress only; partial admission counts are not biological or spatial results.",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review-artifact", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = audit_cache(
        args.review_artifact.resolve(), args.cache_dir.resolve(), args.model_dir.resolve()
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] == "invalid_cache":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
