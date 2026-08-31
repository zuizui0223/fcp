#!/usr/bin/env python3
"""Validate the immutable failed JRC development result for ROI v3."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.segformer_roi import (
    summarize_jrc_gate,
    validate_jrc_box_edge_amendment,
    validate_jrc_box_edge_amendment_v2,
    validate_roi_v3_contract,
)


CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v3.json"
AMENDMENT_V1 = ROOT / "docs/supporting/jbi_atlas_roi_v3_jrc_box_edge_amendment_v1.json"
AMENDMENT_V2 = ROOT / "docs/supporting/jbi_atlas_roi_v3_jrc_box_edge_amendment_v2.json"
MANIFEST = ROOT / "docs/supporting/jbi_atlas_roi_v3_development_evidence_manifest.json"
SOURCE = ROOT / "data/atlas/qualification/roi_v3_sources/jrc_flower_detection_source_inventory_v1.csv"
RESULT = ROOT / "data/atlas/qualification/roi_v3_development/jrc_segformer_development_result_v1.json"
ROWS = ROOT / "data/atlas/qualification/roi_v3_development/jrc_segformer_development_rows_v1.csv"

EXPECTED_HASHES = {
    CONTRACT: "99b66870767462a5d8e2de74b365fffa1b23f2480f4befe3b07f77f51633f9ad",
    AMENDMENT_V1: "828540c36969386822c7656bcaa3534a820da520670d379e50304192ec66f484",
    AMENDMENT_V2: "2e5cabd3599c8030288909d64f64e9ce46244883dd3f8ae9bfacf033e82413c2",
    SOURCE: "c30d40a5333327c1c6d8ea183d274d8a784a6990bf3cc98e69e158db6a8e9e4c",
    RESULT: "c29ec31e105b8054dab4766e252d8295ed4235893dc7d87f4ed402770036662e",
    ROWS: "40e3ebe75c21c779c8690f3073b7ec4c542b7aaf814860e4a5bcab157ee0dafe",
}


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def exact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def normalized_row(row: dict[str, str]) -> dict[str, Any]:
    value: dict[str, Any] = dict(row)
    if row["estimator_admitted"] not in {"True", "False"}:
        raise RuntimeError("invalid estimator admission value")
    value["estimator_admitted"] = row["estimator_admitted"] == "True"
    integer_fields = (
        "flower_pixels",
        "plant_background_control_pixels",
        "predicted_flower_pixels",
        "predicted_flower_pixels_inside_box_union",
        "reference_boxes",
        "source_annotation_boxes",
        "source_not_evaluable_boxes",
        "hit_boxes",
        "small_reference_boxes",
        "medium_reference_boxes",
        "large_reference_boxes",
        "small_hit_boxes",
        "medium_hit_boxes",
        "large_hit_boxes",
    )
    float_fields = (
        "horizontal_flip_mask_iou",
        "horizontal_flip_colour_delta_e",
        "image_predicted_pixel_precision_inside_box_union",
    )
    for field in integer_fields:
        value[field] = int(row[field])
    for field in float_fields:
        value[field] = float(row[field])
    return value


def same_value(observed: Any, expected: Any) -> bool:
    if isinstance(expected, float):
        return math.isclose(float(observed), expected, rel_tol=1e-12, abs_tol=1e-12)
    return observed == expected


def main() -> None:
    for path, expected in EXPECTED_HASHES.items():
        if canonical_sha256(path) != expected:
            raise RuntimeError(f"committed ROI v3 evidence identity changed: {path.name}")

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    amendment_v1 = json.loads(AMENDMENT_V1.read_text(encoding="utf-8"))
    amendment_v2 = json.loads(AMENDMENT_V2.read_text(encoding="utf-8"))
    validate_roi_v3_contract(contract)
    validate_jrc_box_edge_amendment(amendment_v1)
    validate_jrc_box_edge_amendment_v2(amendment_v2)

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    evidence = manifest.get("evidence", {})
    firewall = manifest.get("firewall", {})
    if (
        manifest.get("protocol") != "jbi-atlas-roi-v3-development-evidence-v1"
        or manifest.get("status") != "stop_jrc_development_failed"
        or manifest.get("parent_contract", {}).get("sha256_lf_canonical_v1")
        != EXPECTED_HASHES[CONTRACT]
        or [item.get("sha256_lf_canonical_v1") for item in manifest.get("annotation_amendments", [])]
        != [EXPECTED_HASHES[AMENDMENT_V1], EXPECTED_HASHES[AMENDMENT_V2]]
        or evidence.get("result_sha256_lf_canonical_v1") != EXPECTED_HASHES[RESULT]
        or evidence.get("result_sha256_exact")
        != "6be2b006aead4e27b4b8bb86444af5cb1d71652fafeb9b8682ed9902b547674c"
        or evidence.get("rows_sha256_lf_canonical_v1") != EXPECTED_HASHES[ROWS]
        or any(value is not False for value in firewall.values())
    ):
        raise RuntimeError("ROI v3 development evidence manifest changed")
    if exact_sha256(RESULT) != evidence["result_sha256_exact"]:
        raise RuntimeError("ROI v3 development exact result bytes changed")

    source_rows = [row for row in read_csv(SOURCE) if row["split"] == "train"]
    result_rows_raw = read_csv(ROWS)
    source_identity = {(row["image_id"], row["file_name"], row["image_sha256"]) for row in source_rows}
    result_identity = {
        (row["image_id"], row["file_name"], row["image_sha256"])
        for row in result_rows_raw
    }
    if (
        len(source_rows) != 400
        or len(result_rows_raw) != 400
        or len(result_identity) != 400
        or result_identity != source_identity
        or any(row["split"] != "train" for row in result_rows_raw)
    ):
        raise RuntimeError("JRC development row denominator or identity changed")

    recomputed = summarize_jrc_gate(
        [normalized_row(row) for row in result_rows_raw], contract, phase="development"
    )
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    if recomputed["status"] != "stop_jrc_development_failed":
        raise RuntimeError("failed JRC development decision no longer reproduces")
    for section in ("metrics", "checks"):
        if set(recomputed[section]) != set(result[section]):
            raise RuntimeError(f"JRC development {section} fields changed")
        for key, expected in result[section].items():
            if not same_value(recomputed[section][key], expected):
                raise RuntimeError(f"JRC development {section} changed: {key}")

    if (
        result.get("status") != "stop_jrc_development_failed"
        or result.get("contract_sha256_lf_canonical_v1") != EXPECTED_HASHES[CONTRACT]
        or result.get("box_edge_amendment_sha256_lf_canonical_v1")
        != EXPECTED_HASHES[AMENDMENT_V2]
        or result.get("rows_sha256") != EXPECTED_HASHES[ROWS]
        or result.get("source_annotation_sha256")
        != "7e1d7b45b5720fcf8463aab9c8b154bc03870aa7036a350d95fc40b8c36f35ac"
        or result.get("jrc_locked_test_permitted") is not False
        or result.get("jrc_test_images_decoded_or_scored") is not False
        or result.get("scaleout_candidate_pixels_opened") is not False
        or result.get("atlas_pixels_permitted_by_roi_v3") is not False
        or result["metrics"].get("source_annotation_boxes") != 6992
        or result["metrics"].get("reference_boxes") != 6991
        or result["metrics"].get("source_not_evaluable_boxes") != 1
    ):
        raise RuntimeError("ROI v3 failure or firewall identity changed")

    print(
        json.dumps(
            {
                "status": "pass_committed_failed_roi_v3_development_evidence",
                "development_images": 400,
                "admitted_images": result["metrics"]["admitted_images"],
                "pooled_object_recall": result["metrics"]["pooled_object_recall"],
                "jrc_locked_test_permitted": False,
                "scaleout_candidate_pixels_opened": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
