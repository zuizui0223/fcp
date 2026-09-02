#!/usr/bin/env python3
"""Build the v5 location-blind 60k measurement firewall after all pre-image gates."""

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

from fcp_pipeline.atlas_inference_cascade import validate_contract as validate_v5_inference
from fcp_pipeline.atlas_measurement_v5 import (
    build_v5_measurement_firewall,
    validate_measurement_execution_contract,
    validate_preimage_gates,
)


DEFAULT_CONTRACT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
DEFAULT_INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"
DEFAULT_ROI = ROOT / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"
DEFAULT_SHARED = ROOT / "docs/supporting/jbi_atlas_shared_transition_v5_signal_recovery_result.json"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty v5 firewall table: {path}")
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


def git_blob_sha(path: Path) -> str:
    payload = path.read_bytes()
    return hashlib.sha1(f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def verify_repo_parent_blobs(contract: Mapping[str, Any]) -> None:
    parents = contract["immutable_parents"]
    for key in (
        "inference_v5",
        "real_colour_inference_v5",
        "inference_v3_technical_measurement_rules_only",
        "scaleout_worker_execution_v1",
        "roi_v4_contract",
        "roi_v4_locked_result",
        "shared_transition_preimage_qualification",
        "source_role_v5",
        "source_v3_stop_result",
        "source_v4_stop_result",
        "dated_source_uuid_v4_historical_parent",
    ):
        row = parents[key]
        path = ROOT / row["path"]
        require(path.is_file(), f"missing immutable v5 measurement parent: {row['path']}")
        actual = git_blob_sha(path)
        require(
            actual == row["git_blob_sha"],
            f"immutable v5 measurement parent changed: {row['path']} {actual}",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation-manifest", type=Path, required=True)
    parser.add_argument("--source-v5-manifest", type=Path, required=True)
    parser.add_argument("--source-v5-result", type=Path, required=True)
    parser.add_argument("--environmental-coverage-result", type=Path, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--inference-v5", type=Path, default=DEFAULT_INFERENCE)
    parser.add_argument("--roi-locked-result", type=Path, default=DEFAULT_ROI)
    parser.add_argument("--shared-qualification-result", type=Path, default=DEFAULT_SHARED)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = load_json(args.measurement_contract)
    inference = load_json(args.inference_v5)
    validate_v5_inference(inference)
    validate_measurement_execution_contract(contract, inference)
    verify_repo_parent_blobs(contract)

    observation_hash = sha256(args.observation_manifest)
    source_manifest = load_json(args.source_v5_manifest)
    source_result = load_json(args.source_v5_result)
    environmental = load_json(args.environmental_coverage_result)
    roi = load_json(args.roi_locked_result)
    shared = load_json(args.shared_qualification_result)
    validate_preimage_gates(
        dated_reconciliation=source_result,
        dated_manifest=source_manifest,
        observation_manifest_name=args.observation_manifest.name,
        observation_manifest_sha256=observation_hash,
        environmental_coverage=environmental,
        roi_locked_result=roi,
        shared_qualification_result=shared,
        contract=contract,
    )

    rows = read_csv(args.observation_manifest)
    require(len(rows) == 60000, "v5 firewall input must preserve all 60,000 rows")
    require(
        len({row.get("photo_id", "") for row in rows}) == 60000,
        "v5 firewall input photo IDs are not unique",
    )
    require(
        len({row.get("observation_id", "") for row in rows}) == 60000,
        "v5 firewall input live observation IDs are not unique",
    )
    split = build_v5_measurement_firewall(rows, contract)
    require(
        len(split["measurement_manifest"]) == 60000,
        "v5 blind measurement denominator changed",
    )
    require(
        len(split["sealed_coordinate_key"]) == 60000,
        "v5 sealed acquisition denominator changed",
    )

    worker_path = args.output_dir / "worker_packet" / "measurement_manifest.csv"
    species_path = args.output_dir / "sealed_keys" / "species_key.csv"
    coordinate_path = args.output_dir / "sealed_keys" / "acquisition_coordinate_key.csv"
    write_csv(worker_path, split["measurement_manifest"])
    write_csv(species_path, split["sealed_species_key"])
    write_csv(coordinate_path, split["sealed_coordinate_key"])

    manifest = {
        "status": "pass_scaleout_measurement_firewall_v5",
        "protocol": contract["protocol"],
        "inference_version": inference["version"],
        "frozen_measurements": 60000,
        "candidate_image_pixels_opened": False,
        "terminal_scaleout_colour_measured": False,
        "coordinate_key_opened_by_measurement_worker": False,
        "superseded_v3_ordered_inference_used": False,
        "preimage_gate_sha256": {
            "source_v5_manifest": sha256(args.source_v5_manifest),
            "source_v5_result": sha256(args.source_v5_result),
            "environmental_coverage": sha256(args.environmental_coverage_result),
            "roi_locked_result": sha256(args.roi_locked_result),
            "shared_transition_qualification": sha256(args.shared_qualification_result),
        },
        "contract_git_blob_sha": git_blob_sha(args.measurement_contract),
        "inference_v5_git_blob_sha": git_blob_sha(args.inference_v5),
        "real_colour_inference_v5_git_blob_sha": contract["immutable_parents"]["real_colour_inference_v5"]["git_blob_sha"],
        "source_observation_manifest_sha256": observation_hash,
        "worker_packet": {
            "path": worker_path.relative_to(args.output_dir).as_posix(),
            "sha256": sha256(worker_path),
            "allowed_fields": list(split["measurement_manifest"][0]),
        },
        "sealed_keys": {
            species_path.name: sha256(species_path),
            coordinate_path.name: sha256(coordinate_path),
        },
        "claim_ceiling": contract["claim_ceiling"],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "measurement_firewall_v5_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
