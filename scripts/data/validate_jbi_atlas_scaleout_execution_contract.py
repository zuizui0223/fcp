#!/usr/bin/env python3
"""Validate frozen scale-out worker partition and retry rules."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import validate_execution_contract


CONTRACT = ROOT / "docs/supporting/jbi_atlas_scaleout_execution_contract_v1.json"


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    value = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_execution_contract(value)
    for parent in value["parents"].values():
        path = ROOT / parent["path"]
        if canonical_sha256(path) != parent["sha256_lf_canonical_v1"]:
            raise RuntimeError(f"scale-out execution parent changed: {path}")
    print(
        json.dumps(
            {
                "status": "pass_scaleout_execution_contract",
                "shards": value["measurement"]["shard_count"],
                "maximum_concurrent_workers": value["measurement"][
                    "maximum_concurrent_workers"
                ],
                "candidate_image_pixels_opened": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
