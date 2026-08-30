#!/usr/bin/env python3
"""Build a temporary 120/80 pseudo-split for calibration ROI recovery.

The committed 80/120 calibration/evaluation assignment is never altered.  This script
creates a workspace-only view in which all 80 original calibration photographs plus a
fixed 40-photo padding subset of the original evaluation photographs are marked as
``evaluation``.  The unchanged Florence evaluation extractor can then process exactly
120 photographs per species with its original 6 x 20 sharding contract.  Downstream,
only original calibration photo IDs are retained; padding rows are discarded before
colour measurement.

Padding is selected by the already-frozen split rank (then photo_id as a deterministic
tie-breaker).  No image, colour, coordinate, model output, or biological outcome is
consulted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ORIGINAL_CALIBRATION_PER_SPECIES = 80
ORIGINAL_EVALUATION_PER_SPECIES = 120
PSEUDO_EVALUATION_PER_SPECIES = 120
PSEUDO_CALIBRATION_PER_SPECIES = 80
PADDING_PER_SPECIES = 40


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


def assignment_digest(rows: list[dict[str, str]]) -> str:
    payload = "".join(
        f"{row['species']}\t{scalar_id(row['photo_id'])}\t{row['split']}\n"
        for row in sorted(rows, key=lambda r: (r["species"], scalar_id(r["photo_id"])))
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_pseudosplit(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    if not rows:
        raise ValueError("source split is empty")
    by_species: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    photo_ids: set[str] = set()
    for source in rows:
        row = dict(source)
        species = str(row.get("species", "")).strip()
        photo_id = scalar_id(row.get("photo_id"))
        split = str(row.get("split", "")).strip().lower()
        if not species or not photo_id or split not in {"calibration", "evaluation"}:
            raise ValueError("every row requires species, photo_id, and calibration/evaluation split")
        if photo_id in photo_ids:
            raise ValueError(f"duplicate photo_id {photo_id}")
        photo_ids.add(photo_id)
        row["photo_id"] = photo_id
        row["original_split"] = split
        by_species[species].append(row)

    output: list[dict[str, str]] = []
    per_species: dict[str, Any] = {}
    for species in sorted(by_species):
        species_rows = by_species[species]
        calibration = [row for row in species_rows if row["original_split"] == "calibration"]
        evaluation = [row for row in species_rows if row["original_split"] == "evaluation"]
        if len(calibration) != ORIGINAL_CALIBRATION_PER_SPECIES:
            raise ValueError(
                f"{species}: expected {ORIGINAL_CALIBRATION_PER_SPECIES} original calibration rows, "
                f"found {len(calibration)}"
            )
        if len(evaluation) != ORIGINAL_EVALUATION_PER_SPECIES:
            raise ValueError(
                f"{species}: expected {ORIGINAL_EVALUATION_PER_SPECIES} original evaluation rows, "
                f"found {len(evaluation)}"
            )
        evaluation.sort(
            key=lambda row: (
                str(row.get("split_rank_hash", "")),
                scalar_id(row.get("photo_id")),
            )
        )
        padding_ids = {scalar_id(row["photo_id"]) for row in evaluation[:PADDING_PER_SPECIES]}
        role_counts: Counter[str] = Counter()
        split_counts: Counter[str] = Counter()
        for source in species_rows:
            row = dict(source)
            photo_id = scalar_id(row["photo_id"])
            if row["original_split"] == "calibration":
                row["split"] = "evaluation"
                role = "target_original_calibration"
            elif photo_id in padding_ids:
                row["split"] = "evaluation"
                role = "padding_original_evaluation"
            else:
                row["split"] = "calibration"
                role = "unused_original_evaluation"
            row["calibration_recovery_role"] = role
            role_counts[role] += 1
            split_counts[row["split"]] += 1
            output.append(row)
        expected_splits = {
            "calibration": PSEUDO_CALIBRATION_PER_SPECIES,
            "evaluation": PSEUDO_EVALUATION_PER_SPECIES,
        }
        if dict(split_counts) != expected_splits:
            raise RuntimeError(f"{species}: pseudo split counts {dict(split_counts)} != {expected_splits}")
        per_species[species] = {
            "pseudo_split_counts": dict(sorted(split_counts.items())),
            "role_counts": dict(sorted(role_counts.items())),
            "padding_photo_ids_sha256": hashlib.sha256(
                "\n".join(sorted(padding_ids)).encode("utf-8")
            ).hexdigest(),
        }

    if len(output) != len(rows):
        raise RuntimeError("pseudo-split changed the number of photographs")
    return output, per_species


def write_csv(path: Path, rows: list[dict[str, str]], source_fields: list[str]) -> None:
    extra = ["original_split", "calibration_recovery_role"]
    fields = list(source_fields)
    for field in extra:
        if field not in fields:
            fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-split", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--output-split", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    args = parser.parse_args()

    with args.source_split.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("source split has no header")
        rows = [dict(row) for row in reader]
        source_fields = list(reader.fieldnames)

    output, per_species = build_pseudosplit(rows)
    write_csv(args.output_split, output, source_fields)

    original_manifest: dict[str, Any] = {}
    if args.source_manifest and args.source_manifest.exists():
        original_manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    species = sorted(per_species)
    manifest = {
        **original_manifest,
        "protocol": "jbi_ch1_calibration_recovery_pseudosplit_v1",
        "temporary_workspace_only": True,
        "committed_split_unchanged": True,
        "purpose": "apply the unchanged frozen Florence ROI extractor to original calibration photos",
        "selection_basis": "original split plus frozen split_rank_hash/photo_id padding order",
        "uses_image_content_for_selection": False,
        "uses_colour_for_selection": False,
        "uses_coordinates_for_selection": False,
        "uses_outcomes_for_selection": False,
        "source_split_sha256": sha256_file(args.source_split),
        "output_csv_sha256": sha256_file(args.output_split),
        "assignment_sha256": assignment_digest(output),
        "expected_species": species,
        "total_rows": len(output),
        "calibration_per_species": PSEUDO_CALIBRATION_PER_SPECIES,
        "evaluation_per_species": PSEUDO_EVALUATION_PER_SPECIES,
        "original_calibration_per_species": ORIGINAL_CALIBRATION_PER_SPECIES,
        "original_evaluation_per_species": ORIGINAL_EVALUATION_PER_SPECIES,
        "padding_original_evaluation_per_species": PADDING_PER_SPECIES,
        "per_species_split_counts": {
            species_name: values["pseudo_split_counts"]
            for species_name, values in per_species.items()
        },
        "per_species_recovery_counts": per_species,
        "outcome_blind": True,
        "evaluation_opened_for_rule_tuning": False,
        "final_labels_created": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"built temporary pseudo-split: {len(output)} rows, {len(species)} species, "
        f"{PSEUDO_EVALUATION_PER_SPECIES} Florence-target rows/species"
    )


if __name__ == "__main__":
    main()
