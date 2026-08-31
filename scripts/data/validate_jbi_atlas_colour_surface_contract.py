#!/usr/bin/env python3
"""Validate the pre-pixel atlas colour-surface and inference contract."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_colour_inference import validate_colour_inference_contract
from fcp_pipeline.atlas_environment import validate_environment_contract
from fcp_pipeline.atlas_measurement import validate_inference_contract
from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract


CONTRACT = ROOT / "docs/supporting/jbi_atlas_colour_surface_contract_v1.json"


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_colour_inference_contract(contract)
    parents = contract["parents"]
    for parent in parents.values():
        path = ROOT / parent["path"]
        if canonical_sha256(path) != parent["sha256_lf_canonical_v1"]:
            raise RuntimeError(f"atlas colour-inference parent changed: {path.name}")
    inference = json.loads(
        (ROOT / parents["inference_v3"]["path"]).read_text(encoding="utf-8")
    )
    environment = json.loads(
        (ROOT / parents["environment_v1"]["path"]).read_text(encoding="utf-8")
    )
    roi = json.loads(
        (ROOT / parents["roi_v4"]["path"]).read_text(encoding="utf-8")
    )
    overlay = json.loads(
        (ROOT / parents["qualified_overlay_null"]["path"]).read_text(
            encoding="utf-8"
        )
    )
    validate_inference_contract(inference)
    validate_environment_contract(environment)
    validate_roi_v4_contract(roi)
    if (
        overlay.get("status") != "pass_spatially_constrained_overlay_null"
        or overlay.get("environmental_and_pollinator_colour_join_permitted")
        is not True
    ):
        raise RuntimeError("qualified overlay null no longer authorizes inference")
    print(
        json.dumps(
            {
                "status": "pass_preimage_colour_surface_contract",
                "scaleout_candidate_pixels_opened": False,
                "coordinate_colour_key_joined": False,
                "randomizations": contract["joint_inference"]["randomizations"],
                "primary_scales_km": contract["transition_surface"]["scales_km"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
