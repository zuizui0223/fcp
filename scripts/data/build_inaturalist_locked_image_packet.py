#!/usr/bin/env python3
"""Build a location-free all-photo locked packet for development-passing species."""

from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

if __package__:
    from .build_inaturalist_development_review_packet import (
        REVIEW_FIELDS,
        file_extension,
        load_observation_details,
        materialize_photo,
        quantile,
        read_csv,
        reviewer_rows,
        reusable_photos,
        technical_profile,
        write_csv,
    )
    from .build_inaturalist_photo_development_sample import medium_url, sha256
    from .extract_inaturalist_automated_colour_states import PROTOCOL_VERSION
else:
    from build_inaturalist_development_review_packet import (  # type: ignore[no-redef]
        REVIEW_FIELDS,
        file_extension,
        load_observation_details,
        materialize_photo,
        quantile,
        read_csv,
        reviewer_rows,
        reusable_photos,
        technical_profile,
        write_csv,
    )
    from build_inaturalist_photo_development_sample import medium_url, sha256  # type: ignore[no-redef]
    from extract_inaturalist_automated_colour_states import PROTOCOL_VERSION  # type: ignore[no-redef]


DEFAULT_SEED = "fcp-inaturalist-automated-colour-locked-v2"
EXPECTED_LOCKED_PER_SPECIES = 120


def development_passed_species(gate: dict[str, Any]) -> list[str]:
    if gate.get("protocol") != PROTOCOL_VERSION:
        raise RuntimeError("development-gate protocol mismatch")
    if gate.get("status") != "complete_location_free_automated_colour_development_gate":
        raise RuntimeError("development gate is not complete")
    species = sorted(
        row["canonical_name"]
        for row in gate.get("species_results", [])
        if row.get("development_gate_status") == "pass"
    )
    if len(species) != int(gate.get("species_passed", -1)):
        raise RuntimeError("development-gate species count mismatch")
    if not species:
        raise RuntimeError("no development-passing species; locked packet must remain unopened")
    return species


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-artifact", type=Path, required=True)
    parser.add_argument("--development-gate", type=Path, required=True)
    parser.add_argument("--output-artifact", type=Path, required=True)
    parser.add_argument("--public-manifest", type=Path, required=True)
    parser.add_argument("--seed", default=DEFAULT_SEED)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.source_artifact.resolve()
    output = args.output_artifact.resolve()
    gate = json.loads(args.development_gate.read_text(encoding="utf-8"))
    passed_species = development_passed_species(gate)
    provenance = read_csv(source / "private_provenance.csv")
    locked = [
        row
        for row in provenance
        if row["annotation_partition"] == "locked_60"
        and row["canonical_name"] in passed_species
    ]
    counts = Counter(row["canonical_name"] for row in locked)
    if set(counts) != set(passed_species) or set(counts.values()) != {EXPECTED_LOCKED_PER_SPECIES}:
        raise RuntimeError(f"locked partition mismatch: {dict(sorted(counts.items()))}")
    if len({row["blind_id"] for row in locked}) != len(locked):
        raise RuntimeError("locked encounter blind IDs are not unique")

    output.mkdir(parents=True, exist_ok=True)
    details = load_observation_details(source / "raw_api")
    encounters: list[dict[str, Any]] = []
    photo_ledger: list[dict[str, Any]] = []
    technical_rows: list[dict[str, Any]] = []
    species_stats: dict[str, Counter[str]] = defaultdict(Counter)
    acquired_at = datetime.now(timezone.utc).isoformat()
    for row in sorted(locked, key=lambda value: value["blind_id"]):
        observation_id = int(row["observation_id"])
        observation = details.get(observation_id)
        if observation is None:
            raise RuntimeError(f"missing cached detail for locked observation {observation_id}")
        photos = reusable_photos(observation)
        if not photos:
            raise RuntimeError(f"locked observation {observation_id} has no reusable photo")
        representative_id = int(row["photo_id"])
        if representative_id not in {int(photo["id"]) for photo in photos}:
            raise RuntimeError(f"representative photo missing from reusable locked set: {observation_id}")
        image_files: list[str] = []
        for index, photo in enumerate(photos, start=1):
            photo_blind_id = f"{row['blind_id']}-P{index:02d}"
            relative = Path("images") / row["blind_id"] / (
                photo_blind_id + file_extension(photo)
            )
            target = output / relative
            existing = (
                source / "images" / row["image_file"]
                if int(photo["id"]) == representative_id
                else None
            )
            url = medium_url(photo)
            materialize_photo(url, target, existing)
            profile = technical_profile(target)
            image_files.append(relative.as_posix())
            photo_ledger.append(
                {
                    "encounter_blind_id": row["blind_id"],
                    "photo_blind_id": photo_blind_id,
                    "canonical_name": row["canonical_name"],
                    "observation_id": observation_id,
                    "photo_id": int(photo["id"]),
                    "photo_license": str(photo.get("license_code", "")).casefold(),
                    "photo_attribution": photo.get("attribution", ""),
                    "photo_medium_url": url,
                    "image_file": relative.as_posix(),
                    "image_sha256": sha256(target),
                    "acquired_at_utc": acquired_at,
                }
            )
            technical_rows.append(
                {
                    "encounter_blind_id": row["blind_id"],
                    "photo_blind_id": photo_blind_id,
                    "canonical_name": row["canonical_name"],
                    "image_file": relative.as_posix(),
                    **profile,
                    "technical_status": "profile_complete_no_biological_admission",
                }
            )
        encounters.append(
            {
                "canonical_name": row["canonical_name"],
                "encounter_blind_id": row["blind_id"],
                "image_files": image_files,
            }
        )
        species_stats[row["canonical_name"]]["locked_encounters"] += 1
        species_stats[row["canonical_name"]]["reusable_photos"] += len(photos)
        species_stats[row["canonical_name"]]["multi_photo_encounters"] += int(len(photos) > 1)

    compatibility_rows = reviewer_rows(encounters, "automated_locked_input", args.seed)
    photo_path = output / "private_photo_provenance.csv"
    technical_path = output / "technical_image_profile.csv"
    input_path = output / "reviewer_A_annotation_sheet.csv"
    write_csv(photo_path, photo_ledger)
    write_csv(technical_path, technical_rows)
    write_csv(input_path, compatibility_rows, REVIEW_FIELDS)
    if any(
        field in compatibility_rows[0]
        for field in ("latitude", "longitude", "observer_id", "observed_on", "annotation_partition")
    ):
        raise RuntimeError("locked image input sheet leaked spatial or observer fields")

    public_species: list[dict[str, Any]] = []
    for species in passed_species:
        counter = species_stats[species]
        profiles = [row for row in technical_rows if row["canonical_name"] == species]
        public_species.append(
            {
                "canonical_name": species,
                "locked_encounters": counter["locked_encounters"],
                "reusable_photos": counter["reusable_photos"],
                "multi_photo_encounters": counter["multi_photo_encounters"],
                "photos_per_encounter_mean": counter["reusable_photos"]
                / counter["locked_encounters"],
                "image_width_min": min(row["image_width"] for row in profiles),
                "image_height_min": min(row["image_height"] for row in profiles),
                "luminance_p01_median": quantile(
                    [float(row["luminance_p01"]) for row in profiles], 0.5
                ),
                "luminance_p99_median": quantile(
                    [float(row["luminance_p99"]) for row in profiles], 0.5
                ),
            }
        )
    private_manifest = {
        "status": "complete_location_free_locked_image_packet_pending_automated_measurement",
        "protocol": "fcp-inaturalist-automated-colour-locked-packet-v2",
        "automated_colour_protocol": PROTOCOL_VERSION,
        "created_at_utc": acquired_at,
        "source_artifact": str(source),
        "source_artifact_manifest_sha256": sha256(source / "artifact_manifest.json"),
        "development_gate_sha256": sha256(args.development_gate),
        "annotation_partition": "locked_60_only",
        "passed_species": passed_species,
        "locked_encounters": len(encounters),
        "locked_photos": len(photo_ledger),
        "multi_photo_encounters": sum(
            counter["multi_photo_encounters"] for counter in species_stats.values()
        ),
        "human_judgement_used": False,
        "input_sheet_role": "location-free extractor compatibility only; blank human fields are ignored",
        "coordinates_joined": False,
        "species_summary": public_species,
        "private_files_sha256": {
            path.name: sha256(path) for path in (photo_path, technical_path, input_path)
        },
        "claim_ceiling": (
            "Locked image acquisition and integrity only. No candidate colour, spatial organization, "
            "botanical morph, population frequency, mechanism or universality is established."
        ),
    }
    private_path = output / "artifact_manifest.json"
    private_path.write_text(
        json.dumps(private_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    public = {
        key: value
        for key, value in private_manifest.items()
        if key not in {"source_artifact", "private_files_sha256"}
    }
    public["source_artifact_reference"] = "private moving-API development sample"
    public["artifact_manifest_sha256"] = sha256(private_path)
    args.public_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.public_manifest.write_text(
        json.dumps(public, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(public, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
