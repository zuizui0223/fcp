#!/usr/bin/env python3
"""Build the location-blind 60k-image worker packet and sealed join keys."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_expansion import validate_expansion_contract
from fcp_pipeline.atlas_measurement import (
    build_measurement_firewall,
    validate_inference_contract,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write an empty firewall table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation-manifest", type=Path, required=True)
    parser.add_argument(
        "--expansion-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json"),
    )
    parser.add_argument(
        "--inference-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.expansion_contract.read_text(encoding="utf-8"))
    validate_expansion_contract(contract)
    inference = json.loads(args.inference_contract.read_text(encoding="utf-8"))
    validate_inference_contract(inference)
    rows = read_csv(args.observation_manifest)
    expected = int(contract["random_cohort_scaleout"]["total_observations"])
    if len(rows) != expected:
        raise RuntimeError(
            f"firewall input must equal the frozen {expected}-observation denominator"
        )
    split = build_measurement_firewall(
        rows, salt=str(inference["scaleout_measurement_gate"]["blinding_salt"])
    )
    worker_path = args.output_dir / "worker_packet" / "measurement_manifest.csv"
    species_path = args.output_dir / "sealed_keys" / "species_key.csv"
    coordinate_path = args.output_dir / "sealed_keys" / "acquisition_coordinate_key.csv"
    write_csv(worker_path, split["measurement_manifest"])
    write_csv(species_path, split["sealed_species_key"])
    write_csv(coordinate_path, split["sealed_coordinate_key"])
    manifest = {
        "status": "pass_scaleout_measurement_firewall",
        "protocol": inference["protocol"],
        "frozen_measurements": len(rows),
        "candidate_image_pixels_opened": False,
        "coordinate_key_opened_by_measurement_worker": False,
        "worker_packet": {
            "path": worker_path.relative_to(args.output_dir).as_posix(),
            "sha256": sha256(worker_path),
            "allowed_fields": list(split["measurement_manifest"][0]),
        },
        "sealed_keys": {
            species_path.name: sha256(species_path),
            coordinate_path.name: sha256(coordinate_path),
        },
        "source_observation_manifest_sha256": sha256(args.observation_manifest),
        "claim_ceiling": (
            "Pre-measurement blinding and denominator evidence only; no image, colour, "
            "spatial, environmental or pollinator conclusion is allowed."
        ),
    }
    manifest_path = args.output_dir / "measurement_firewall_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
