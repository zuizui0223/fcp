#!/usr/bin/env python3
"""Validate the prospective flower-specific ROI v4 contract and parents."""

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
V3_MANIFEST = ROOT / "docs/supporting/jbi_atlas_roi_v3_development_evidence_manifest.json"
V3_RESULT = ROOT / "data/atlas/qualification/roi_v3_development/jrc_segformer_development_result_v1.json"


def exact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_roi_v4_contract(contract)
    v3_manifest = json.loads(V3_MANIFEST.read_text(encoding="utf-8"))
    v3_result = json.loads(V3_RESULT.read_text(encoding="utf-8"))
    parent = contract["immutable_v3_stop"]
    if (
        v3_manifest.get("status") != "stop_jrc_development_failed"
        or v3_result.get("status") != "stop_jrc_development_failed"
        or exact_sha256(V3_RESULT) != parent["result_sha256_exact"]
        or v3_result.get("jrc_locked_test_permitted") is not False
        or v3_result.get("jrc_test_images_decoded_or_scored") is not False
        or v3_result.get("scaleout_candidate_pixels_opened") is not False
    ):
        raise RuntimeError("immutable ROI v3 STOP changed")
    print(
        json.dumps(
            {
                "status": "pass_prospective_roi_v4_contract",
                "jrc_train_images": 400,
                "jrc_locked_test_images": 100,
                "jrc_locked_test_images_decoded_or_scored": False,
                "scaleout_candidate_pixels_opened": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
