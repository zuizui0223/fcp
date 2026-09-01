#!/usr/bin/env python3
"""Validate the terminal atlas v5 pre-image measurement freeze and immutable parents."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_inference_cascade import validate_contract as validate_v5_inference
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs


CONTRACT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"
ROI = ROOT / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"
SHARED = ROOT / "docs/supporting/jbi_atlas_shared_transition_v5_signal_recovery_result.json"


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    contract = load(CONTRACT)
    inference = load(INFERENCE)
    validate_v5_inference(inference)
    validate_measurement_execution_contract(contract, inference)
    verify_repo_parent_blobs(contract)

    roi = load(ROI)
    if (
        roi.get("status") != "pass_roi_v4_locked_test"
        or roi.get("scaleout_candidate_pixels_opened") is not False
        or roi.get("scaleout_candidate_pixels_permitted") is not True
    ):
        raise RuntimeError("committed ROI v4 locked result no longer authorizes scaleout")

    shared = load(SHARED)
    scope = shared.get("scope", {})
    if (
        shared.get("status") != "pass"
        or scope.get("method_gate_only") is not True
        or scope.get("biological_support_claimed") is not False
        or scope.get("candidate_image_pixels_opened") is not False
    ):
        raise RuntimeError("committed shared-transition result is not a preimage method pass")

    print(
        json.dumps(
            {
                "status": "pass_v5_measurement_execution_freeze",
                "candidate_image_pixels_opened": False,
                "terminal_inference_contract": inference["version"],
                "shared_transition_gate_scope": "method_only",
                "roi_v4_locked": "pass",
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
