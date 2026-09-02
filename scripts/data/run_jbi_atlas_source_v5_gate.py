#!/usr/bin/env python3
"""Validate immutable terminal selection plus already-proven dated snapshot identity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_source_v5 import (
    PASS_LABEL,
    audit_terminal_source_rows,
    build_source_v5_result,
    validate_snapshot_identity_evidence,
    validate_source_v5_contract,
)

DEFAULT_CONTRACT = ROOT / "docs/supporting/jbi_atlas_source_role_amendment_v5.json"
DEFAULT_V3_STOP = ROOT / "docs/supporting/jbi_atlas_dated_source_streaming_v3_stop_result.json"
DEFAULT_V4_STOP = ROOT / "docs/supporting/jbi_atlas_dated_source_uuid_v4_stop_result.json"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--v3-stop", type=Path, default=DEFAULT_V3_STOP)
    parser.add_argument("--v4-stop", type=Path, default=DEFAULT_V4_STOP)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = load_json(args.contract)
    validate_source_v5_contract(contract)
    require(
        git_blob_sha(args.v3_stop) == contract["trigger"]["v3_stop"]["git_blob_sha"],
        "preserved v3 STOP blob changed",
    )
    require(
        git_blob_sha(args.v4_stop) == contract["trigger"]["v4_stop"]["git_blob_sha"],
        "preserved v4 STOP blob changed",
    )
    v3_stop = load_json(args.v3_stop)
    validate_snapshot_identity_evidence(v3_stop, contract)

    selection = contract["immutable_scientific_selection"]
    for filename, expected_sha in selection["files"].items():
        path = args.metadata_dir / filename
        require(path.is_file(), f"missing immutable terminal selection file: {filename}")
        require(sha256(path) == expected_sha, f"immutable terminal selection SHA changed: {filename}")
    manifest_path = args.metadata_dir / "scaleout_metadata_manifest.json"
    require(manifest_path.is_file(), "terminal metadata artifact manifest missing")
    artifact_manifest = load_json(manifest_path)
    require(
        artifact_manifest.get("candidate_image_pixels_opened") is False,
        "terminal metadata artifact says candidate pixels were opened",
    )
    feasibility = load_json(args.metadata_dir / "scaleout_metadata_feasibility.json")
    require(
        feasibility.get("candidate_image_pixels_opened") is False
        and feasibility.get("continuous_colour_used") is False,
        "terminal feasibility audit contains image outcomes",
    )

    rows = read_csv(args.metadata_dir / "scaleout_observation_manifest.csv")
    panels = read_csv(args.metadata_dir / "scaleout_species_panels.csv")
    source_audit = audit_terminal_source_rows(rows, panels, contract)
    result = build_source_v5_result(source_audit=source_audit, contract=contract)
    result["parents"] = {
        "contract_git_blob_sha": git_blob_sha(args.contract),
        "v3_stop_git_blob_sha": git_blob_sha(args.v3_stop),
        "v4_stop_git_blob_sha": git_blob_sha(args.v4_stop),
        "terminal_metadata_artifact_digest": selection["artifact_digest"],
        "scaleout_metadata_feasibility_sha256": sha256(args.metadata_dir / "scaleout_metadata_feasibility.json"),
        "scaleout_observation_manifest_sha256": sha256(args.metadata_dir / "scaleout_observation_manifest.csv"),
        "scaleout_species_panels_sha256": sha256(args.metadata_dir / "scaleout_species_panels.csv"),
    }
    result["dated_snapshot_identity"] = {
        "snapshot_date": contract["dated_snapshot_provenance"]["snapshot_date"],
        "bytes": v3_stop["snapshot_identity"]["bytes_read"],
        "sha256": v3_stop["snapshot_identity"]["computed_sha256"],
        "identity_passed": v3_stop["snapshot_identity"]["identity_passed"],
        "reused_existing_full_stream_proof": True,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.output_dir / "source_v5_result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "protocol": contract["protocol"],
        "status": result["status"],
        "candidate_image_pixels_opened": False,
        "replacement_permitted": False,
        "files": {result_path.name: sha256(result_path)},
    }
    manifest_path = args.output_dir / "source_v5_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == PASS_LABEL else 2


if __name__ == "__main__":
    raise SystemExit(main())
