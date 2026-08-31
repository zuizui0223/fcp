#!/usr/bin/env python3
"""Create the shareable locked analysis table without exposing observer IDs.

The coordinate firewall has already been evaluated by
``finalize_inaturalist_locked_colour_table.py``.  This exporter does not decide
which species or encounters are evaluable.  It only removes coordinates from
non-admitted encounters and replaces observer identifiers with deterministic
protocol-scoped group labels so the leave-top-observer sensitivity remains
exactly reproducible.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


PROTOCOL = "fcp-inaturalist-automated-colour-state-v2"
EXPECTED_SPECIES = 3
EXPECTED_ROWS_PER_SPECIES = 120
ADMITTED = "automated_colour_state_admitted"
OBSERVER_DOMAIN = f"{PROTOCOL}:public-observer-group-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def observer_group(species: str, observer_id: str, rank: int) -> str:
    if not observer_id:
        raise ValueError("admitted encounter has no observer identifier")
    token = "\x1f".join((OBSERVER_DOMAIN, species, observer_id))
    digest = hashlib.sha256(token.encode("utf-8")).hexdigest()[:8]
    return f"observer_{rank:04d}_{digest}"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def export_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts = Counter(row["canonical_name"] for row in rows)
    if len(counts) != EXPECTED_SPECIES or set(counts.values()) != {
        EXPECTED_ROWS_PER_SPECIES
    }:
        raise RuntimeError("locked table must contain three species x 120 encounters")
    if len({row["encounter_blind_id"] for row in rows}) != len(rows):
        raise RuntimeError("duplicate encounter blind ID")

    observer_maps: dict[str, dict[str, str]] = {}
    for species in sorted(counts):
        identifiers = sorted(
            {
                row["observer_id"]
                for row in rows
                if row["canonical_name"] == species
                and row["encounter_status"] == ADMITTED
            }
        )
        if not identifiers or any(not value for value in identifiers):
            raise RuntimeError(f"{species}: admitted observer identifier is missing")
        observer_maps[species] = {
            value: observer_group(species, value, rank)
            for rank, value in enumerate(identifiers, start=1)
        }

    exported: list[dict[str, str]] = []
    for row in rows:
        item = dict(row)
        if item["encounter_status"] == ADMITTED:
            for field in ("latitude", "longitude", "observer_id"):
                if not item.get(field, "").strip():
                    raise RuntimeError(f"admitted encounter missing {field}")
            item["observer_id"] = observer_maps[item["canonical_name"]][
                item["observer_id"]
            ]
        else:
            item["latitude"] = ""
            item["longitude"] = ""
            item["observer_id"] = ""
        exported.append(item)
    return exported


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_rows(args.input)
    exported = export_rows(rows)
    write_csv(args.output, exported)
    admitted = [row for row in exported if row["encounter_status"] == ADMITTED]
    manifest: dict[str, Any] = {
        "status": "complete_public_locked_analysis_table",
        "protocol": PROTOCOL,
        "rows": len(exported),
        "admitted_rows": len(admitted),
        "species": sorted({row["canonical_name"] for row in exported}),
        "privacy_transform": {
            "observer_id": (
                "species-scoped rank-preserving pseudonym with an 8-hex "
                "protocol-scoped SHA-256 suffix"
            ),
            "non_admitted_coordinates_removed": True,
            "observation_and_photo_identifiers_present": False,
        },
        "source_private_sha256": sha256(args.input),
        "public_table_sha256": sha256(args.output),
        "exporter_sha256": sha256(Path(__file__).resolve()),
        "claim_ceiling": (
            "Reproducible input for the frozen three-species random-mark test only. "
            "Rows are model-consensus flower-candidate colour measurements, not "
            "verified flower tissue, botanical morphs or population frequencies."
        ),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
