#!/usr/bin/env python3
"""Validate the original dated-source STOP and its pre-colour M:M successor."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_dated_source import validate_dated_source_amendment
from fcp_pipeline.atlas_dated_source_m2m import validate_m2m_amendment


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    v1_path = ROOT / "docs/supporting/jbi_atlas_dated_source_amendment_v1.json"
    v1 = json.loads(v1_path.read_text(encoding="utf-8"))
    validate_dated_source_amendment(v1)
    v1_sha = sha256(v1_path)
    stop_path = ROOT / "docs/supporting/jbi_atlas_dated_source_v1_stop_result.json"
    stop = json.loads(stop_path.read_text(encoding="utf-8"))
    if (
        stop.get("status") != "not_evaluable_dated_source_reconciliation"
        or stop.get("candidate_image_pixels_opened") is not False
        or stop.get("continuous_colour_used") is not False
        or stop.get("selected_photo_association_rows_inspected") is not False
        or stop.get("replacement_permitted") is not False
        or stop.get("v1_retry_permitted") is not False
        or stop.get("parent_amendment_sha256_exact") != v1_sha
    ):
        raise RuntimeError("dated-source v1 STOP evidence changed")
    v2_path = ROOT / "docs/supporting/jbi_atlas_dated_source_m2m_amendment_v2.json"
    v2 = json.loads(v2_path.read_text(encoding="utf-8"))
    validate_m2m_amendment(v2)
    if v2["trigger"]["v1_stop_result"] != stop_path.relative_to(ROOT).as_posix():
        raise RuntimeError("dated-source v2 no longer points to the v1 STOP")
    print(
        json.dumps(
            {
                "status": "pass",
                "v1_protocol": v1["protocol"],
                "v1_outcome": stop["status"],
                "active_protocol": v2["protocol"],
            }
        )
    )


if __name__ == "__main__":
    main()
