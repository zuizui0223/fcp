#!/usr/bin/env python3
"""Freeze one completed ROI-v4 JRC gate without changing its decision."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4_evidence import sha256, validate_gate_artifacts


CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
RUNNER_PATH = "scripts/data/run_jbi_atlas_roi_v4_jrc_gate.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("development", "locked_test"), required=True)
    parser.add_argument("--gate-run-dir", type=Path, required=True)
    parser.add_argument("--trained-weight", type=Path, required=True)
    parser.add_argument("--training-evidence-manifest", type=Path, required=True)
    parser.add_argument("--execution-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gate_run_dir = args.gate_run_dir.resolve()
    trained_weight = args.trained_weight.resolve()
    training_evidence_manifest = args.training_evidence_manifest.resolve()
    output_dir = args.output_dir.resolve()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    weight_sha = sha256(trained_weight)
    training = json.loads(
        training_evidence_manifest.read_text(encoding="utf-8")
    )
    training_weight = training.get("evidence", {}).get("trained_weight", {})
    if (
        training.get("status")
        != "complete_roi_v4_training_frozen_before_prediction"
        or training_weight.get("sha256_exact") != weight_sha
        or any(value is not False for value in training.get("firewall", {}).values())
    ):
        raise RuntimeError("committed training evidence does not authorize a JRC gate")
    prefix = f"jrc_roi_v4_{args.phase}"
    rows_source = gate_run_dir / f"{prefix}_rows.csv"
    result_source = gate_run_dir / f"{prefix}_result.json"
    result = validate_gate_artifacts(
        rows_source,
        result_source,
        contract,
        phase=args.phase,
        trained_weight_sha256=weight_sha,
    )
    if output_dir.exists():
        raise RuntimeError("refusing to replace frozen ROI v4 gate evidence")
    output_dir.mkdir(parents=True)
    rows_path = output_dir / rows_source.name
    result_path = output_dir / result_source.name
    shutil.copy2(rows_source, rows_path)
    shutil.copy2(result_source, result_path)
    runner_payload = subprocess.run(
        ["git", "show", f"{args.execution_commit}:{RUNNER_PATH}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    runner_path = output_dir / "gate_executable.py"
    runner_path.write_bytes(runner_payload)
    evidence = {
        "protocol": f"jbi-atlas-roi-v4-{args.phase}-evidence-v1",
        "status": result["status"],
        "phase": args.phase,
        "trained_weight_sha256": weight_sha,
        "execution": {
            "git_commit": args.execution_commit,
            "runner_path": RUNNER_PATH,
            "runner_sha256_exact": sha256(runner_path),
        },
        "training_evidence": {
            "path": str(training_evidence_manifest.relative_to(ROOT)).replace(
                "\\", "/"
            ),
            "sha256_exact": sha256(training_evidence_manifest),
        },
        "evidence": {
            "rows": {
                "path": str(rows_path.relative_to(ROOT)).replace("\\", "/"),
                "sha256_exact": sha256(rows_path),
                "rows": 400 if args.phase == "development" else 100,
            },
            "result": {
                "path": str(result_path.relative_to(ROOT)).replace("\\", "/"),
                "sha256_exact": sha256(result_path),
            },
            "runner": {
                "path": str(runner_path.relative_to(ROOT)).replace("\\", "/"),
                "sha256_exact": sha256(runner_path),
            },
        },
        "firewall": {
            "jrc_locked_test_images_decoded_or_scored": args.phase == "locked_test",
            "scaleout_candidate_pixels_opened": False,
        },
        "next_gate": (
            "run the locked JRC test exactly once"
            if args.phase == "development" and result["jrc_locked_test_permitted"]
            else "open the scaleout acquisition gate"
            if args.phase == "locked_test"
            and result["scaleout_candidate_pixels_permitted"]
            else "STOP this estimator version"
        ),
        "claim_ceiling": "Estimator gate evidence only; no atlas colour or ecological outcome.",
    }
    manifest_path = output_dir / "gate_evidence_manifest.json"
    manifest_path.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(evidence, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
