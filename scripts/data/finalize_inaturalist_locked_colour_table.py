#!/usr/bin/env python3
"""Apply locked completeness gates and join coordinates only for passing species."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

if __package__:
    from .extract_inaturalist_automated_colour_states import PROTOCOL_VERSION, sha256
else:
    from extract_inaturalist_automated_colour_states import PROTOCOL_VERSION, sha256  # type: ignore[no-redef]


EXPECTED_ENCOUNTERS_PER_SPECIES = 120
MIN_ADMITTED_SHARE = 0.70
MIN_BACKGROUND_SHARE = 0.70


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write empty locked table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def locked_species_gate(rows: list[dict[str, str]]) -> dict[str, Any]:
    if len(rows) != EXPECTED_ENCOUNTERS_PER_SPECIES:
        raise RuntimeError("locked species does not contain exactly 120 encounters")
    admitted = [
        row for row in rows if row["encounter_status"] == "automated_colour_state_admitted"
    ]
    background = [
        row
        for row in admitted
        if row["background_control_status"] == "background_control_available"
    ]
    admission_share = len(admitted) / len(rows)
    background_share = len(background) / len(admitted) if admitted else 0.0
    passed = admission_share >= MIN_ADMITTED_SHARE and background_share >= MIN_BACKGROUND_SHARE
    return {
        "locked_gate_status": "pass" if passed else "not_evaluable",
        "locked_encounters": len(rows),
        "admitted_encounters": len(admitted),
        "admitted_encounter_share": admission_share,
        "background_control_encounters": len(background),
        "background_control_encounter_share": background_share,
        "coordinates_joined": passed,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-dir", type=Path, required=True)
    parser.add_argument("--locked-packet", type=Path, required=True)
    parser.add_argument("--source-artifact", type=Path, required=True)
    parser.add_argument("--development-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--public-manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    extraction = args.extraction_dir.resolve()
    packet = args.locked_packet.resolve()
    source = args.source_artifact.resolve()
    output = args.output_dir.resolve()
    development = json.loads(args.development_gate.read_text(encoding="utf-8"))
    development_species = sorted(
        row["canonical_name"]
        for row in development["species_results"]
        if row["development_gate_status"] == "pass"
    )
    if development.get("protocol") != PROTOCOL_VERSION or not development_species:
        raise RuntimeError("invalid or empty development gate")
    run_manifest_path = extraction / "run_manifest.json"
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    encounter_path = extraction / "encounter_features.csv"
    photo_path = extraction / "photo_features.csv"
    encounters = read_csv(encounter_path)
    photos = read_csv(photo_path)
    packet_manifest_path = packet / "artifact_manifest.json"
    packet_manifest = json.loads(packet_manifest_path.read_text(encoding="utf-8"))
    if run_manifest.get("protocol") != PROTOCOL_VERSION:
        raise RuntimeError("locked extraction protocol mismatch")
    if run_manifest.get("status") != "complete_automated_colour_state_development_feasibility_not_spatial":
        raise RuntimeError("locked extraction is incomplete")
    if run_manifest.get("selected_encounters") != len(development_species) * EXPECTED_ENCOUNTERS_PER_SPECIES:
        raise RuntimeError("locked extraction encounter count mismatch")
    if run_manifest.get("selected_photos") != packet_manifest.get("locked_photos"):
        raise RuntimeError("locked extraction photo count mismatch")
    if len(encounters) != run_manifest["selected_encounters"] or len(photos) != run_manifest["selected_photos"]:
        raise RuntimeError("locked feature table count mismatch")
    for path in (encounter_path, photo_path):
        if run_manifest.get("private_output_sha256", {}).get(path.name) != sha256(path):
            raise RuntimeError(f"locked feature hash mismatch: {path.name}")
    counts = Counter(row["canonical_name"] for row in encounters)
    if set(counts) != set(development_species) or set(counts.values()) != {
        EXPECTED_ENCOUNTERS_PER_SPECIES
    }:
        raise RuntimeError("locked feature species balance mismatch")
    if len({row["encounter_blind_id"] for row in encounters}) != len(encounters):
        raise RuntimeError("duplicate locked encounter feature ID")

    species_results: list[dict[str, Any]] = []
    gate_by_species: dict[str, dict[str, Any]] = {}
    for species in development_species:
        gate = locked_species_gate(
            [row for row in encounters if row["canonical_name"] == species]
        )
        gate_by_species[species] = gate
        species_results.append({"canonical_name": species, **gate})

    provenance = read_csv(source / "private_provenance.csv")
    provenance = [
        row
        for row in provenance
        if row["annotation_partition"] == "locked_60"
        and row["canonical_name"] in development_species
    ]
    provenance_by_id = {row["blind_id"]: row for row in provenance}
    if len(provenance_by_id) != len(development_species) * EXPECTED_ENCOUNTERS_PER_SPECIES:
        raise RuntimeError("locked private provenance mismatch")
    joined: list[dict[str, Any]] = []
    for row in sorted(encounters, key=lambda value: (value["canonical_name"], value["encounter_blind_id"])):
        private = provenance_by_id.get(row["encounter_blind_id"])
        if private is None or private["canonical_name"] != row["canonical_name"]:
            raise RuntimeError("locked feature/provenance identity mismatch")
        coordinate_allowed = gate_by_species[row["canonical_name"]]["coordinates_joined"]
        joined.append(
            {
                **row,
                "annotation_partition": "locked_60",
                "latitude": private["latitude"] if coordinate_allowed else "",
                "longitude": private["longitude"] if coordinate_allowed else "",
                "observer_id": private["observer_id"] if coordinate_allowed else "",
            }
        )
    for row in joined:
        allowed = gate_by_species[row["canonical_name"]]["coordinates_joined"]
        has_sensitive = all(str(row[field]).strip() for field in ("latitude", "longitude", "observer_id"))
        if has_sensitive != allowed:
            raise RuntimeError("coordinate firewall mismatch")

    output.mkdir(parents=True, exist_ok=True)
    joined_path = output / "locked_spatial_input.csv"
    write_csv(joined_path, joined)
    public = {
        "status": "complete_locked_automatic_colour_gate_and_coordinate_firewall",
        "protocol": PROTOCOL_VERSION,
        "fixed_gates": {
            "expected_encounters_per_species": EXPECTED_ENCOUNTERS_PER_SPECIES,
            "min_admitted_encounter_share": MIN_ADMITTED_SHARE,
            "min_background_control_share": MIN_BACKGROUND_SHARE,
        },
        "species_results": species_results,
        "species_passing_locked_gate": sum(
            row["locked_gate_status"] == "pass" for row in species_results
        ),
        "coordinate_joined_species": [
            row["canonical_name"] for row in species_results if row["coordinates_joined"]
        ],
        "coordinate_withheld_species": [
            row["canonical_name"] for row in species_results if not row["coordinates_joined"]
        ],
        "locked_spatial_input_sha256": sha256(joined_path),
        "source_sha256": {
            "development_gate": sha256(args.development_gate),
            "locked_extraction_manifest": sha256(run_manifest_path),
            "locked_packet_manifest": sha256(packet_manifest_path),
            "finalizer": sha256(Path(__file__).resolve()),
        },
        "claim_ceiling": (
            "Locked image-measurement completeness and coordinate firewall only. A pass permits the "
            "predeclared random-mark test; it is not spatial, botanical-morph or mechanism evidence."
        ),
    }
    private_report = {
        **public,
        "private_locked_spatial_input": str(joined_path),
        "source_artifact": str(source),
    }
    private_path = output / "locked_gate_report.json"
    private_path.write_text(
        json.dumps(private_report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    public["private_report_sha256"] = sha256(private_path)
    args.public_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.public_manifest.write_text(
        json.dumps(public, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(public, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
