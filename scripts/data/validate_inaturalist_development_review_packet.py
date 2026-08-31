#!/usr/bin/env python3
"""Fail-closed validator for the private iNaturalist development review packet."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from PIL import Image

if __package__:
    from .build_inaturalist_development_review_packet import (
        PROHIBITED_REVIEW_FIELDS,
        REVIEW_FIELDS,
    )
    from .build_inaturalist_photo_development_sample import sha256
else:
    from build_inaturalist_development_review_packet import (  # type: ignore[no-redef]
        PROHIBITED_REVIEW_FIELDS,
        REVIEW_FIELDS,
    )
    from build_inaturalist_photo_development_sample import sha256  # type: ignore[no-redef]


EXPECTED_SPECIES = {
    "Digitalis purpurea",
    "Erythranthe lewisii",
    "Hepatica nobilis",
    "Hesperis matronalis",
    "Orchis mascula",
    "Protea repens",
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def validate_packet(artifact: Path, public_manifest_path: Path) -> dict[str, Any]:
    errors: list[str] = []
    private_manifest_path = artifact / "artifact_manifest.json"
    private_manifest = json.loads(private_manifest_path.read_text(encoding="utf-8"))
    public_manifest = json.loads(public_manifest_path.read_text(encoding="utf-8"))

    expected_hash = public_manifest.get("artifact_manifest_sha256")
    if expected_hash != sha256(private_manifest_path):
        errors.append("public/private artifact manifest hash mismatch")
    public_serialized = public_manifest_path.read_text(encoding="utf-8")
    if "C:\\\\Users\\\\" in public_serialized or "C:/Users/" in public_serialized:
        errors.append("public manifest exposes a local absolute user path")

    for filename, expected in private_manifest.get("private_files_sha256", {}).items():
        path = artifact / filename
        if not path.exists():
            errors.append(f"missing private file: {filename}")
        elif sha256(path) != expected:
            errors.append(f"private file hash mismatch: {filename}")

    ledger_fields, ledger = read_csv(artifact / "private_photo_provenance.csv")
    profile_fields, profiles = read_csv(artifact / "technical_image_profile.csv")
    a_fields, reviewer_a = read_csv(artifact / "reviewer_A_annotation_sheet.csv")
    b_fields, reviewer_b = read_csv(artifact / "reviewer_B_annotation_sheet.csv")
    codebook_fields, codebook = read_csv(artifact / "species_codebook_template.csv")

    required_schemas = {
        "private_photo_provenance.csv": (
            set(ledger_fields),
            {"encounter_blind_id", "photo_blind_id", "canonical_name", "image_file", "image_sha256"},
        ),
        "technical_image_profile.csv": (
            set(profile_fields),
            {"encounter_blind_id", "photo_blind_id", "canonical_name", "technical_status"},
        ),
        "reviewer_A_annotation_sheet.csv": (set(a_fields), set(REVIEW_FIELDS)),
        "reviewer_B_annotation_sheet.csv": (set(b_fields), set(REVIEW_FIELDS)),
        "species_codebook_template.csv": (set(codebook_fields), {"canonical_name"}),
    }
    schema_failed = False
    for filename, (actual, required) in required_schemas.items():
        missing = sorted(required - actual)
        if missing:
            schema_failed = True
            errors.append(f"missing required columns in {filename}: {missing}")
    if schema_failed:
        raise RuntimeError(json.dumps({"status": "fail", "errors": errors}, indent=2, sort_keys=True))

    if len(ledger) != 886 or len(profiles) != 886:
        errors.append(f"expected 886 photo rows, got ledger={len(ledger)} profiles={len(profiles)}")
    if len(reviewer_a) != 480 or len(reviewer_b) != 480:
        errors.append(
            f"expected 480 encounter rows per reviewer, got A={len(reviewer_a)} B={len(reviewer_b)}"
        )
    if a_fields != REVIEW_FIELDS or b_fields != REVIEW_FIELDS:
        errors.append("reviewer sheet schema differs from the frozen review schema")
    if PROHIBITED_REVIEW_FIELDS.intersection(a_fields) or PROHIBITED_REVIEW_FIELDS.intersection(
        b_fields
    ):
        errors.append("reviewer sheet contains prohibited provenance columns")
    if set(ledger_fields).intersection({"flower_visible", "anonymous_morph_code"}):
        errors.append("private provenance unexpectedly contains outcome labels")
    if "technical_status" not in profile_fields:
        errors.append("technical profile lacks its non-biological status field")

    ledger_photo_ids = [row["photo_blind_id"] for row in ledger]
    profile_photo_ids = [row["photo_blind_id"] for row in profiles]
    if len(set(ledger_photo_ids)) != len(ledger_photo_ids):
        errors.append("duplicate photo blind ID in private provenance")
    if set(ledger_photo_ids) != set(profile_photo_ids):
        errors.append("technical profile and provenance photo keys differ")

    a_ids = [row["encounter_blind_id"] for row in reviewer_a]
    b_ids = [row["encounter_blind_id"] for row in reviewer_b]
    if len(set(a_ids)) != 480 or set(a_ids) != set(b_ids):
        errors.append("reviewer encounter keys are duplicate or disagree")
    if a_ids == b_ids:
        errors.append("reviewer encounter orderings are not independent")
    if any(row["reviewer_id"] or row["anonymous_morph_code"] for row in reviewer_a + reviewer_b):
        errors.append("review sheet already contains reviewer or morph outcomes")

    encounter_photo_counts = Counter(row["encounter_blind_id"] for row in ledger)
    for row in reviewer_a:
        files = row["image_files"].split("|") if row["image_files"] else []
        if int(row["n_photos"]) != len(files):
            errors.append(f"reviewer image list count mismatch: {row['encounter_blind_id']}")
        if int(row["n_photos"]) != encounter_photo_counts[row["encounter_blind_id"]]:
            errors.append(f"reviewer/provenance photo count mismatch: {row['encounter_blind_id']}")

    missing_images = 0
    bad_hashes = 0
    bad_decodes = 0
    for row in ledger:
        image_path = artifact / row["image_file"]
        if not image_path.exists():
            missing_images += 1
            continue
        if sha256(image_path) != row["image_sha256"]:
            bad_hashes += 1
        try:
            with Image.open(image_path) as image:
                image.verify()
        except Exception:
            bad_decodes += 1
    if missing_images or bad_hashes or bad_decodes:
        errors.append(
            f"image integrity failed: missing={missing_images} hash={bad_hashes} decode={bad_decodes}"
        )

    species_counts = Counter(row["canonical_name"] for row in reviewer_a)
    if set(species_counts) != EXPECTED_SPECIES or set(species_counts.values()) != {80}:
        errors.append(f"unexpected development species counts: {dict(species_counts)}")
    if {row["canonical_name"] for row in codebook} != EXPECTED_SPECIES:
        errors.append("codebook species universe differs from reviewer universe")
    if any(path.suffix == ".partial" for path in artifact.rglob("*.partial")):
        errors.append("partial image file remains in artifact")

    result = {
        "status": "pass" if not errors else "fail",
        "development_encounters": len(reviewer_a),
        "development_photos": len(ledger),
        "unique_photo_blind_ids": len(set(ledger_photo_ids)),
        "multi_photo_encounters": sum(count > 1 for count in encounter_photo_counts.values()),
        "species_counts": dict(sorted(species_counts.items())),
        "missing_images": missing_images,
        "bad_image_hashes": bad_hashes,
        "bad_image_decodes": bad_decodes,
        "reviewer_prohibited_columns": sorted(
            PROHIBITED_REVIEW_FIELDS.intersection(a_fields) | PROHIBITED_REVIEW_FIELDS.intersection(b_fields)
        ),
        "errors": errors,
        "claim_ceiling": (
            "This validates review readiness and data integrity only; it does not validate "
            "flower visibility, morph labels, reviewer agreement, or spatial organization."
        ),
    }
    if errors:
        raise RuntimeError(json.dumps(result, indent=2, sort_keys=True))
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--public-manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = validate_packet(args.artifact.resolve(), args.public_manifest.resolve())
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
