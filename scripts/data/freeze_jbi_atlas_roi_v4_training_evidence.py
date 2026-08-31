#!/usr/bin/env python3
"""Freeze the completed ROI v4 training run before any JRC prediction."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract


CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
TRAINER_PATH = "scripts/data/train_jbi_atlas_roi_v4_detector.py"
FROZEN_FILENAMES = {
    "trained_weight": "jrc_yolo11n_last_v4.pt",
    "training_results": "training_results.csv",
    "training_arguments": "training_args.yaml",
    "training_result": "training_result_manifest.json",
    "materialization": "training_materialization_manifest.json",
    "trainer": "training_executable.py",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256_bytes(payload)


def git_blob(commit: str, path: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return completed.stdout


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-run-dir", type=Path, required=True)
    parser.add_argument("--materialization-manifest", type=Path, required=True)
    parser.add_argument("--training-code-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def validate_source_run(
    run_dir: Path, materialization_path: Path, contract: dict[str, object]
) -> tuple[dict[str, object], dict[str, object], dict[str, Path]]:
    paths = {
        "trained_weight": run_dir / "jrc_yolo11n_last_v4.pt",
        "training_results": run_dir / "frozen_train/results.csv",
        "training_arguments": run_dir / "frozen_train/args.yaml",
        "training_result": run_dir / "training_result_manifest.json",
        "materialization": materialization_path,
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise RuntimeError(f"ROI v4 training evidence incomplete: {missing}")
    result = json.loads(paths["training_result"].read_text(encoding="utf-8"))
    materialization = json.loads(
        materialization_path.read_text(encoding="utf-8")
    )
    if (
        result.get("protocol") != contract["protocol"]
        or result.get("status")
        != "complete_roi_v4_detector_training_not_yet_qualified"
        or result.get("epochs") != 50
        or result.get("weight_selection") != "last epoch only; never best epoch"
        or result.get("trained_weight_sha256") != sha256(paths["trained_weight"])
        or result.get("training_results_sha256")
        != sha256(paths["training_results"])
        or result.get("jrc_test_images_decoded_or_scored") is not False
        or result.get("scaleout_candidate_pixels_opened") is not False
    ):
        raise RuntimeError("completed ROI v4 training result changed or is incomplete")
    if (
        materialization.get("protocol") != contract["protocol"]
        or materialization.get("status")
        != "pass_roi_v4_training_materialization"
        or materialization.get("images") != 400
        or materialization.get("evaluable_training_boxes") != 6991
        or materialization.get("training_validation_enabled") is not False
        or materialization.get("jrc_test_directory_read") is not False
        or materialization.get("jrc_test_images_decoded_or_scored") is not False
        or materialization.get("scaleout_candidate_pixels_opened") is not False
    ):
        raise RuntimeError("ROI v4 training materialization changed")
    return result, materialization, paths


def main() -> None:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_roi_v4_contract(contract)
    result, materialization, source_paths = validate_source_run(
        args.training_run_dir, args.materialization_manifest, contract
    )
    trainer_payload = git_blob(args.training_code_commit, TRAINER_PATH)
    if not trainer_payload.startswith(b"#!/usr/bin/env python3"):
        raise RuntimeError("training executable identity could not be resolved")
    if args.output_dir.exists():
        raise RuntimeError("refusing to replace frozen ROI v4 training evidence")
    args.output_dir.mkdir(parents=True)
    frozen_paths: dict[str, Path] = {}
    for key, source in source_paths.items():
        destination = args.output_dir / FROZEN_FILENAMES[key]
        shutil.copy2(source, destination)
        frozen_paths[key] = destination
    trainer_path = args.output_dir / FROZEN_FILENAMES["trainer"]
    trainer_path.write_bytes(trainer_payload)
    frozen_paths["trainer"] = trainer_path
    evidence = {
        "protocol": "jbi-atlas-roi-v4-training-evidence-v1",
        "status": "complete_roi_v4_training_frozen_before_prediction",
        "parent_contract": {
            "path": str(CONTRACT.relative_to(ROOT)).replace("\\", "/"),
            "sha256_lf_canonical_v1": canonical_sha256(CONTRACT),
        },
        "training_execution": {
            "git_commit": args.training_code_commit,
            "trainer_path": TRAINER_PATH,
            "trainer_sha256_exact": sha256_bytes(trainer_payload),
            "weight_selection": result["weight_selection"],
            "epochs": result["epochs"],
            "environment": result["environment"],
        },
        "evidence": {
            key: {
                "path": str(path.relative_to(ROOT)).replace("\\", "/"),
                "sha256_exact": sha256(path),
                "bytes": path.stat().st_size,
            }
            for key, path in frozen_paths.items()
        },
        "denominators": {
            "training_images": materialization["images"],
            "source_annotation_boxes": materialization["source_annotation_boxes"],
            "evaluable_training_boxes": materialization[
                "evaluable_training_boxes"
            ],
            "source_not_evaluable_boxes": materialization[
                "source_not_evaluable_boxes"
            ],
        },
        "firewall": {
            "validation_during_training": False,
            "jrc_test_directory_read": False,
            "jrc_test_images_decoded_or_scored": False,
            "oxford102_locked_proxy_scored_by_v4": False,
            "scaleout_candidate_pixels_opened": False,
        },
        "next_gate": "run all 400 JRC development images once with this exact trained weight",
        "claim_ceiling": "Training completion is provenance, not estimator validity or flower-colour evidence.",
    }
    evidence_path = args.output_dir / "training_evidence_manifest.json"
    evidence_path.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(evidence, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
