#!/usr/bin/env python3
"""Validate the location-free locked image packet before colour extraction."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from PIL import Image

if __package__:
    from .build_inaturalist_development_review_packet import (
        PROHIBITED_REVIEW_FIELDS,
        REVIEW_OUTCOME_FIELDS,
        read_csv,
    )
    from .build_inaturalist_photo_development_sample import sha256
else:
    from build_inaturalist_development_review_packet import (  # type: ignore[no-redef]
        PROHIBITED_REVIEW_FIELDS,
        REVIEW_OUTCOME_FIELDS,
        read_csv,
    )
    from build_inaturalist_photo_development_sample import sha256  # type: ignore[no-redef]


EXPECTED_SPECIES = {
    "Erythranthe lewisii",
    "Hesperis matronalis",
    "Orchis mascula",
}
EXPECTED_ENCOUNTERS_PER_SPECIES = 120
EXPECTED_TOTAL_ENCOUNTERS = 360
EXPECTED_TOTAL_PHOTOS = 717


def validate_packet(artifact: Path) -> dict[str, object]:
    manifest_path = artifact / "artifact_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    ledger_path = artifact / "private_photo_provenance.csv"
    technical_path = artifact / "technical_image_profile.csv"
    input_path = artifact / "reviewer_A_annotation_sheet.csv"
    errors: list[str] = []
    if manifest.get("status") != "complete_location_free_locked_image_packet_pending_automated_measurement":
        errors.append("manifest_status")
    if manifest.get("locked_encounters") != EXPECTED_TOTAL_ENCOUNTERS:
        errors.append("manifest_encounter_count")
    if manifest.get("locked_photos") != EXPECTED_TOTAL_PHOTOS:
        errors.append("manifest_photo_count")
    for path in (ledger_path, technical_path, input_path):
        if manifest.get("private_files_sha256", {}).get(path.name) != sha256(path):
            errors.append(f"private_hash:{path.name}")
    ledger = read_csv(ledger_path)
    technical = read_csv(technical_path)
    inputs = read_csv(input_path)
    if len(ledger) != EXPECTED_TOTAL_PHOTOS or len(technical) != EXPECTED_TOTAL_PHOTOS:
        errors.append("photo_table_count")
    if len(inputs) != EXPECTED_TOTAL_ENCOUNTERS:
        errors.append("input_encounter_count")
    input_counts = Counter(row["canonical_name"] for row in inputs)
    if set(input_counts) != EXPECTED_SPECIES or set(input_counts.values()) != {
        EXPECTED_ENCOUNTERS_PER_SPECIES
    }:
        errors.append("input_species_balance")
    if len({row["encounter_blind_id"] for row in inputs}) != len(inputs):
        errors.append("duplicate_encounter_blind_id")
    ledger_ids = {row["photo_blind_id"] for row in ledger}
    technical_ids = {row["photo_blind_id"] for row in technical}
    if len(ledger_ids) != len(ledger) or len(technical_ids) != len(technical):
        errors.append("duplicate_photo_blind_id")
    if ledger_ids != technical_ids:
        errors.append("photo_id_set_mismatch")
    if PROHIBITED_REVIEW_FIELDS.intersection(inputs[0]):
        errors.append("prohibited_input_column")
    if any(str(row.get(field, "")).strip() for row in inputs for field in REVIEW_OUTCOME_FIELDS):
        errors.append("nonblank_human_outcome")
    bad_hashes = 0
    bad_decodes = 0
    missing_images = 0
    image_paths: set[str] = set()
    for row in ledger:
        relative = row["image_file"]
        image_paths.add(relative)
        path = artifact / relative
        if not path.exists():
            missing_images += 1
            continue
        if sha256(path) != row["image_sha256"]:
            bad_hashes += 1
        try:
            with Image.open(path) as image:
                image.verify()
        except Exception:
            bad_decodes += 1
    referenced = {
        image_file
        for row in inputs
        for image_file in row["image_files"].split("|")
    }
    if referenced != image_paths:
        errors.append("input_image_reference_set_mismatch")
    if missing_images:
        errors.append(f"missing_images:{missing_images}")
    if bad_hashes:
        errors.append(f"bad_image_hashes:{bad_hashes}")
    if bad_decodes:
        errors.append(f"bad_image_decodes:{bad_decodes}")
    partials = [
        path
        for path in artifact.rglob("*")
        if path.is_file() and path.suffix in {".part", ".partial", ".tmp"}
    ]
    if partials:
        errors.append(f"partial_files:{len(partials)}")
    return {
        "status": "pass" if not errors else "fail",
        "protocol": "fcp-inaturalist-automated-colour-locked-packet-validation-v2",
        "artifact_manifest_sha256": sha256(manifest_path),
        "locked_encounters": len(inputs),
        "locked_photos": len(ledger),
        "species": sorted(input_counts),
        "human_outcome_cells_nonblank": sum(
            bool(str(row.get(field, "")).strip())
            for row in inputs
            for field in REVIEW_OUTCOME_FIELDS
        ),
        "prohibited_input_columns": sorted(PROHIBITED_REVIEW_FIELDS.intersection(inputs[0])),
        "missing_images": missing_images,
        "bad_image_hashes": bad_hashes,
        "bad_image_decodes": bad_decodes,
        "partial_files": len(partials),
        "errors": errors,
        "coordinates_joined": False,
        "biological_conclusion_allowed": False,
        "claim_ceiling": "Locked packet integrity only; no colour or spatial result.",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = validate_packet(args.artifact.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
