#!/usr/bin/env python3
"""Validate the prospective geometry-admission correction and frozen parents."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_scaleout import validate_geometry_admission_amendment


AMENDMENT = Path(
    "docs/supporting/jbi_atlas_scaleout_geometry_admission_amendment_v1.json"
)


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    amendment = json.loads(AMENDMENT.read_text(encoding="utf-8"))
    validate_geometry_admission_amendment(amendment)
    verified = []
    for name, parent in amendment["parent_contracts"].items():
        path = Path(parent["path"])
        observed = canonical_sha256(path)
        if observed != parent["sha256_lf_canonical_v1"]:
            raise RuntimeError(f"geometry amendment parent changed: {name}")
        verified.append({"name": name, "path": path.as_posix(), "sha256": observed})
    print(
        json.dumps(
            {
                "status": "pass_preimage_scaleout_geometry_admission_amendment",
                "protocol": amendment["protocol"],
                "candidate_image_pixels_opened": False,
                "parent_contracts_verified": verified,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
