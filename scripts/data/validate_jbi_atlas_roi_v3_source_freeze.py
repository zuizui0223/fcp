#!/usr/bin/env python3
"""Validate the committed pre-prediction ROI v3 source inventories."""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.segformer_roi import validate_roi_v3_contract


CONTRACT = Path("docs/supporting/jbi_atlas_roi_estimator_contract_v3.json")
SOURCE_ROOT = Path("data/atlas/qualification/roi_v3_sources")
MANIFEST = SOURCE_ROOT / "roi_v3_source_inventory_manifest.json"
EXPECTED = {
    "jrc_flower_detection_source_inventory_v1.csv": (
        "c30d40a5333327c1c6d8ea183d274d8a784a6990bf3cc98e69e158db6a8e9e4c"
    ),
    "oxford102_roi_proxy_source_inventory_v1.csv": (
        "032a5a79ddee53aad15a7f1b88f848f480706bf4aacb135917ff1324fba3b6cf"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_roi_v3_contract(contract)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if (
        manifest.get("status") != "pass_roi_v3_source_inventory_freeze"
        or manifest.get("contract_sha256_lf_canonical_v1") != canonical_sha256(CONTRACT)
        or manifest.get("scaleout_candidate_pixels_opened") is not False
        or manifest.get("estimator_predictions_run") is not False
        or manifest.get("locked_images_decoded") is not False
        or manifest.get("cross_source_exact_duplicates") != 0
        or manifest.get("files") != EXPECTED
    ):
        raise RuntimeError("ROI v3 source manifest identity changed")
    for name, expected in EXPECTED.items():
        if sha256(SOURCE_ROOT / name) != expected:
            raise RuntimeError(f"ROI v3 source inventory hash changed: {name}")

    jrc = rows(SOURCE_ROOT / "jrc_flower_detection_source_inventory_v1.csv")
    if (
        len(jrc) != 500
        or Counter(row["split"] for row in jrc) != {"train": 400, "test": 100}
        or sum(int(row["annotation_boxes"]) for row in jrc) != 9516
        or len({row["image_sha256"] for row in jrc}) != 500
        or any(
            row["image_pixels_decoded"] != "False"
            or row["estimator_prediction_run"] != "False"
            for row in jrc
        )
    ):
        raise RuntimeError("JRC pre-prediction inventory changed")

    oxford = rows(SOURCE_ROOT / "oxford102_roi_proxy_source_inventory_v1.csv")
    role_counts = Counter(row["role"] for row in oxford)
    locked = [row for row in oxford if row["role"] == "locked_proxy"]
    overlap_ids = sorted(
        int(row["image_id"])
        for row in oxford
        if row["full_oxford17_exact_overlap"] == "True"
    )
    internal_ids = sorted(
        int(row["image_id"])
        for row in oxford
        if int(row["internal_exact_duplicate_group_size"]) > 1
    )
    if (
        len(oxford) != 8189
        or len({row["image_sha256"] for row in oxford}) != 8185
        or role_counts
        != {
            "development": 2036,
            "excluded_exact_duplicate": 13,
            "locked_proxy": 2040,
            "quarantine": 4100,
        }
        or Counter(int(row["class_id"]) for row in locked)
        != {class_id: 20 for class_id in range(1, 103)}
        or overlap_ids != [3448, 3456, 4657, 4691, 6241]
        or internal_ids != [7307, 7309, 7323, 7328, 7885, 7886, 8067, 8077]
        or any(
            row["image_pixels_decoded_for_estimator"] != "False"
            or row["estimator_prediction_run"] != "False"
            for row in oxford
        )
    ):
        raise RuntimeError("Oxford-102 pre-prediction inventory changed")
    print(
        json.dumps(
            {
                "status": "pass_committed_roi_v3_source_inventory",
                "jrc_images": len(jrc),
                "jrc_boxes": 9516,
                "oxford102_images": len(oxford),
                "oxford102_locked_proxy": len(locked),
                "scaleout_candidate_pixels_opened": False,
                "estimator_predictions_run": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
