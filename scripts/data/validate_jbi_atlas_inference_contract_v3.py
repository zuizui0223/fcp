#!/usr/bin/env python3
"""Validate the prospective post-simulation atlas decision tree and evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import validate_inference_contract


DEFAULT_CONTRACT = Path(
    "docs/supporting/jbi_image_first_atlas_inference_contract_v3.json"
)


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def exact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_inference_contract(contract)
    parent = contract["parent_expansion_contract"]
    parent_path = Path(parent["path"])
    if canonical_sha256(parent_path) != parent["sha256_lf_canonical_v1"]:
        raise RuntimeError("v2 expansion contract identity mismatch")
    checked = []
    for branch in ("geographic_shared_boundary", "environmental_pollinator_overlay_null"):
        evidence = contract["qualification_evidence"][branch]
        path = Path(evidence["committed_summary"])
        observed = exact_sha256(path)
        if observed != evidence["committed_summary_sha256"]:
            raise RuntimeError(f"qualification evidence hash mismatch: {branch}")
        checked.append({"branch": branch, "path": path.as_posix(), "sha256": observed})
    print(
        json.dumps(
            {
                "status": "pass_preimage_inference_contract_v3",
                "protocol": contract["protocol"],
                "scaleout_candidate_pixels_opened": False,
                "geographic_branch": "not_evaluable",
                "next_evaluable_branch": "environmental_concordance",
                "qualification_evidence_verified": checked,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
