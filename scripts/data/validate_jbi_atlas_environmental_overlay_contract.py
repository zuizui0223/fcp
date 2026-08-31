#!/usr/bin/env python3
"""Validate the colour-blind environmental overlay freeze contract."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import validate_environment_contract
from fcp_pipeline.atlas_measurement import validate_inference_contract


def main() -> None:
    path = Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json")
    contract = json.loads(path.read_text(encoding="utf-8"))
    validate_environment_contract(contract)
    parent = json.loads(
        Path(contract["parent_inference_contract"]).read_text(encoding="utf-8")
    )
    validate_inference_contract(parent)
    print(
        json.dumps(
            {
                "status": "pass_precolour_environmental_overlay_contract",
                "protocol": contract["protocol"],
                "scaleout_colour_opened": False,
                "families": sorted(contract["sources"]),
                "scales_km": contract["grid"]["scales_km"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
