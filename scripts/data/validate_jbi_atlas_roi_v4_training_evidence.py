#!/usr/bin/env python3
"""Validate committed ROI v4 training evidence before JRC development."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract


CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
EVIDENCE_DIR = ROOT / "data/atlas/qualification/roi_v4_training"
MANIFEST = EVIDENCE_DIR / "training_evidence_manifest.json"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256_bytes(payload)


def main() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_roi_v4_contract(contract)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if (
        manifest.get("protocol") != "jbi-atlas-roi-v4-training-evidence-v1"
        or manifest.get("status")
        != "complete_roi_v4_training_frozen_before_prediction"
        or manifest.get("parent_contract", {}).get("sha256_lf_canonical_v1")
        != canonical_sha256(CONTRACT)
        or any(value is not False for value in manifest.get("firewall", {}).values())
    ):
        raise RuntimeError("ROI v4 training evidence manifest changed")
    execution = manifest["training_execution"]
    automatic_summary = manifest.get("automatic_post_training_summary", {})
    if (
        automatic_summary.get("trainer_argument_val") is not False
        or automatic_summary.get("automatic_post_training_training_set_summary_emitted")
        is not True
        or automatic_summary.get("summary_epoch") != 50
        or automatic_summary.get("summary_metrics_used_for_weight_selection") is not False
        or automatic_summary.get("best_pt_frozen_or_used") is not False
    ):
        raise RuntimeError("ROI v4 post-training summary provenance changed")
    evidence = manifest["evidence"]
    for item in evidence.values():
        path = ROOT / item["path"]
        if (
            not path.is_file()
            or sha256(path) != item["sha256_exact"]
            or path.stat().st_size != item["bytes"]
        ):
            raise RuntimeError(f"ROI v4 training artifact changed: {path.name}")
    trainer_path = ROOT / evidence["trainer"]["path"]
    if sha256(trainer_path) != execution["trainer_sha256_exact"]:
        raise RuntimeError("ROI v4 training executable identity changed")
    result_path = ROOT / evidence["training_result"]["path"]
    result = json.loads(result_path.read_text(encoding="utf-8"))
    weight_path = ROOT / evidence["trained_weight"]["path"]
    results_path = ROOT / evidence["training_results"]["path"]
    materialization_path = ROOT / evidence["materialization"]["path"]
    materialization = json.loads(materialization_path.read_text(encoding="utf-8"))
    if (
        result.get("status")
        != "complete_roi_v4_detector_training_not_yet_qualified"
        or result.get("trained_weight_sha256") != sha256(weight_path)
        or result.get("training_results_sha256") != sha256(results_path)
        or result.get("jrc_test_images_decoded_or_scored") is not False
        or result.get("scaleout_candidate_pixels_opened") is not False
        or materialization.get("jrc_test_directory_read") is not False
        or materialization.get("jrc_test_images_decoded_or_scored") is not False
        or materialization.get("scaleout_candidate_pixels_opened") is not False
    ):
        raise RuntimeError("ROI v4 training provenance no longer closes")
    print(
        json.dumps(
            {
                "status": "pass_committed_roi_v4_training_evidence",
                "trained_weight_sha256": sha256(weight_path),
                "training_images": manifest["denominators"]["training_images"],
                "jrc_test_images_decoded_or_scored": False,
                "scaleout_candidate_pixels_opened": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
