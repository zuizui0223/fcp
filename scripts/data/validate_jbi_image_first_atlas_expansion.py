#!/usr/bin/env python3
"""Validate the prospective FCP atlas expansion before any candidate pixels open."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from fcp_pipeline.atlas_expansion import validate_expansion_contract


DEFAULT_CONTRACT = Path(
    "docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json"
)


def sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_expansion_contract(contract)

    checked: list[dict[str, str]] = []
    for name, item in contract["parent_artifacts"].items():
        path = Path(item["path"])
        if not path.is_file():
            raise RuntimeError(f"missing frozen parent artifact: {path}")
        observed = sha256(path)
        expected = str(item["sha256_lf_canonical_v1"])
        if observed != expected:
            raise RuntimeError(
                f"frozen parent hash mismatch for {name}: expected {expected}, observed {observed}"
            )
        checked.append({"name": name, "path": path.as_posix(), "sha256": observed})

    print(
        json.dumps(
            {
                "status": "pass_preimage_expansion_contract",
                "protocol": contract["protocol"],
                "parent_artifacts_verified": checked,
                "atlas_candidate_pixels_opened": False,
                "random_cohorts": contract["random_cohort_scaleout"]["cohorts"],
                "random_species": contract["random_cohort_scaleout"]["total_species"],
                "random_observations": contract["random_cohort_scaleout"][
                    "total_observations"
                ],
                "next_gate": "independent ROI and signal-recovery qualification",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
