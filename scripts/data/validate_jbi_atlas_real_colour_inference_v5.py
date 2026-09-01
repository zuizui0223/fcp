#!/usr/bin/env python3
"""Validate the exact terminal v5 real-colour inference amendment and parent blobs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_real_inference_v5 import validate_real_inference_amendment


CONTRACT = ROOT / "docs/supporting/jbi_atlas_real_colour_inference_amendment_v5.json"


def git_blob_sha(path: Path) -> str:
    payload = path.read_bytes()
    return hashlib.sha1(f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def main() -> int:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_real_inference_amendment(contract)
    for name, row in contract["immutable_parents"].items():
        path = ROOT / row["path"]
        if not path.is_file():
            raise RuntimeError(f"missing real-colour inference parent: {name}: {row['path']}")
        actual = git_blob_sha(path)
        if actual != row["git_blob_sha"]:
            raise RuntimeError(
                f"real-colour inference parent changed: {row['path']} {actual}"
            )
    print(
        json.dumps(
            {
                "status": "pass_real_colour_inference_v5_freeze",
                "pixel_status": contract["pixel_status_at_freeze"],
                "spatial_randomizations": contract["species_conditioned_spatial_organization"]["randomizations"],
                "shared_randomizations": contract["shared_transition"]["randomizations"],
                "environment_randomizations": contract["environmental_concordance"]["randomizations"],
                "not_evaluable_advances": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
