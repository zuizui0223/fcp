#!/usr/bin/env python3
"""Validate one committed ROI-v4 development or locked-test evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4_evidence import sha256, validate_gate_artifacts
from fcp_pipeline.flower_roi_v4_runtime import validate_scaleout_authorization


CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"


def git_blob_sha256(commit: str, path: str) -> str:
    payload = subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return hashlib.sha256(payload).hexdigest()


def validate_committed_gate(evidence_dir: Path, *, phase: str) -> dict[str, Any]:
    evidence_dir = evidence_dir.resolve()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    manifest_path = evidence_dir / "gate_evidence_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("protocol") != f"jbi-atlas-roi-v4-{phase}-evidence-v1"
        or manifest.get("phase") != phase
        or manifest.get("firewall", {}).get("scaleout_candidate_pixels_opened")
        is not False
        or manifest.get("firewall", {}).get(
            "jrc_locked_test_images_decoded_or_scored"
        )
        is not (phase == "locked_test")
    ):
        raise RuntimeError("ROI v4 gate evidence manifest changed")

    training_item = manifest["training_evidence"]
    training_path = ROOT / training_item["path"]
    if sha256(training_path) != training_item["sha256_exact"]:
        raise RuntimeError("ROI v4 training evidence parent changed")
    training = json.loads(training_path.read_text(encoding="utf-8"))
    weight_item = training.get("evidence", {}).get("trained_weight", {})
    weight_path = ROOT / weight_item["path"]
    weight_sha = sha256(weight_path)
    if (
        training.get("status")
        != "complete_roi_v4_training_frozen_before_prediction"
        or weight_sha != weight_item.get("sha256_exact")
        or weight_sha != manifest.get("trained_weight_sha256")
    ):
        raise RuntimeError("ROI v4 trained weight parent changed")

    files: dict[str, Path] = {}
    for name, item in manifest.get("evidence", {}).items():
        path = ROOT / item["path"]
        if not path.is_file() or sha256(path) != item["sha256_exact"]:
            raise RuntimeError(f"ROI v4 gate artifact changed: {name}")
        files[name] = path
    expected_rows = 400 if phase == "development" else 100
    if manifest["evidence"]["rows"].get("rows") != expected_rows:
        raise RuntimeError("ROI v4 gate row denominator changed")
    execution = manifest["execution"]
    if (
        sha256(files["runner"]) != execution.get("runner_sha256_exact")
        or git_blob_sha256(execution["git_commit"], execution["runner_path"])
        != execution.get("runner_sha256_exact")
    ):
        raise RuntimeError("ROI v4 gate executable identity changed")

    result = validate_gate_artifacts(
        files["rows"],
        files["result"],
        contract,
        phase=phase,
        trained_weight_sha256=weight_sha,
    )
    if manifest.get("status") != result["status"]:
        raise RuntimeError("ROI v4 gate decision changed")
    expected_next = (
        "run the locked JRC test exactly once"
        if phase == "development" and result["jrc_locked_test_permitted"]
        else "open the scaleout acquisition gate"
        if phase == "locked_test" and result["scaleout_candidate_pixels_permitted"]
        else "STOP this estimator version"
    )
    if manifest.get("next_gate") != expected_next:
        raise RuntimeError("ROI v4 next-gate decision changed")
    return {
        "status": "pass_committed_roi_v4_gate_evidence",
        "phase": phase,
        "decision": result["status"],
        "images": result["metrics"]["images"],
        "next_gate": expected_next,
        "scaleout_candidate_pixels_opened": False,
    }


def load_committed_locked_scaleout_result(evidence_dir: Path) -> dict[str, Any]:
    """Load the locked result only after its complete committed bundle validates."""

    evidence_dir = evidence_dir.resolve()
    decision = validate_committed_gate(evidence_dir, phase="locked_test")
    if (
        decision.get("decision") != "pass_roi_v4_locked_test"
        or decision.get("next_gate") != "open the scaleout acquisition gate"
    ):
        raise RuntimeError("committed ROI v4 evidence does not authorize scale-out")
    manifest = json.loads(
        (evidence_dir / "gate_evidence_manifest.json").read_text(encoding="utf-8")
    )
    result_path = ROOT / manifest["evidence"]["result"]["path"]
    result = json.loads(result_path.read_text(encoding="utf-8"))
    validate_scaleout_authorization(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("development", "locked_test"), required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate_committed_gate(args.evidence_dir, phase=args.phase), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
