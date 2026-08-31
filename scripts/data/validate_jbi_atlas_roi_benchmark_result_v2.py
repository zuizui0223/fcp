#!/usr/bin/env python3
"""Validate the immutable Oxford-17 ROI STOP evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


RESULT = Path("docs/supporting/jbi_atlas_oxford17_roi_benchmark_result_v2.json")
ROWS = Path("data/atlas/qualification/oxford17_roi_benchmark_rows_v2.csv")
MANIFEST = Path("docs/supporting/jbi_atlas_roi_benchmark_evidence_manifest_v2.json")
EXPECTED = {
    ROWS.as_posix(): "ed043d92428fe3a282eca5a291ae75c212f6bb48e6531203275b8dbd1796ab3b",
    RESULT.as_posix(): "7db59dd8936f2b72bf95630be59cac90141049c6678ddb0117d011f0c460940a",
}


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    observed_hashes = {path: canonical_sha256(Path(path)) for path in EXPECTED}
    if observed_hashes != EXPECTED or manifest.get("files_sha256_lf_canonical_v1") != EXPECTED:
        raise RuntimeError("Oxford-17 ROI evidence identity changed")
    checks = result.get("checks", {})
    metrics = result.get("metrics", {})
    if (
        result.get("status") != "stop_roi_benchmark_failed"
        or result.get("atlas_pixels_permitted_by_roi_gate") is not False
        or checks.get("minimum_admitted_fraction") is not False
        or not all(value is True for key, value in checks.items() if key != "minimum_admitted_fraction")
        or metrics.get("admitted_images") != 676
        or metrics.get("scored_images") != 848
        or not math.isclose(metrics.get("admitted_fraction", -1), 676 / 848)
        or not metrics.get("admitted_fraction", 1) < manifest.get(
            "frozen_minimum_admitted_fraction", 0
        )
        or manifest.get("decision")
        != "the v2 estimator remains stopped; do not lower the gate or retune on Oxford-17"
    ):
        raise RuntimeError("Oxford-17 ROI STOP decision changed")
    with ROWS.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 848 or len({row["image_id"] for row in rows}) != 848:
        raise RuntimeError("Oxford-17 ROI row denominator changed")
    admitted = sum(row["estimator_admitted"] == "True" for row in rows)
    if admitted != 676 or any(not row["image_sha256"] or not row["trimap_sha256"] for row in rows):
        raise RuntimeError("Oxford-17 ROI row evidence changed")
    print(
        json.dumps(
            {
                "status": result["status"],
                "scored_images": len(rows),
                "admitted_images": admitted,
                "admitted_fraction": metrics["admitted_fraction"],
                "failed_check": "minimum_admitted_fraction",
                "all_other_checks_passed": True,
                "atlas_pixels_permitted": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
