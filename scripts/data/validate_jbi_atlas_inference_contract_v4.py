#!/usr/bin/env python3
"""Validate the strict v4 pre-image atlas decision contract."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_decision_v4 import validate_v4_contract

DEFAULT_CONTRACT = Path("docs/supporting/jbi_image_first_atlas_inference_contract_v4.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_v4_contract(contract, root=ROOT)
    print(json.dumps({
        "status": "pass_preimage_inference_contract_v4",
        "candidate_pixels_opened": False,
        "sampling_repetitions": 8,
        "legacy_6_34_3": "immutable",
        "shared_transition_preimage_state": "not_evaluable",
        "confirmatory_cascade_can_advance": False,
        "advance_only_on": "unsupported"
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
