#!/usr/bin/env python3
"""Run the frozen 110-image Monarda flower-region agreement diagnostic.

This script is intentionally narrow. It verifies the exact user-supplied archive,
checks every exported JPEG/COCO alignment before model execution, and then runs the
unchanged FCP ROI-v4 runtime once on all 110 images. It never trains, tunes, measures
colour, joins geography, or publishes raw images/polygon vertices.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

from scripts.analysis.audit_rgfca_monarda_archive import inspect_archive

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "docs/supporting/rgfca_monarda_region_agreement_contract_v1.json"
DEFAULT_INTAKE = ROOT / "docs/supporting/rgfca_monarda_archive_intake_v1.json"
DEFAULT_MAPPING = ROOT / "docs/supporting/rgfca_monarda_source_mapping_v1.json"
DEFAULT_ROI_CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
DEFAULT_LOCKED = ROOT / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"

SPLIT_ORDER = ("train", "valid", "test")
EXIF_ORIENTATION_TAG = 274


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _round_clip(value: float, maximum: int) -> int:
    if not math.isfinite(float(value)):
        raise ValueError("polygon coordinate is non-finite")
    return max(0, min(int(maximum) - 1, int(round(float(value)))))


def rasterize_polygon_union(
    width: int,
    height: int,
    segmentations: list[list[float]],
) -> np.ndarray:
    """Rasterize the contract-frozen generic-flower polygon union."""
    if int(width) <= 0 or int(height) <= 0:
        raise ValueError("reference dimensions must be positive")
    canvas = Image.new("1", (int(width), int(height)), 0)
    draw = ImageDraw.Draw(canvas)
    for polygon in segmentations:
        if not isinstance(polygon, list) or len(polygon) < 6 or len(polygon) % 2:
            raise ValueError("unsupported COCO polygon component")
        points = [
            (_round_clip(x, width), _round_clip(y, height))
            for x, y in zip(polygon[::2], polygon[1::2])
        ]
        draw.polygon(points, fill=1)
    return np.asarray(canvas, dtype=bool)


def score_masks(reference: np.ndarray, predicted: np.ndarray) -> dict[str, float | int]:
    ref = np.asarray(reference, dtype=bool)
    pred = np.asarray(predicted, dtype=bool)
    if ref.shape != pred.shape or ref.ndim != 2:
        raise ValueError("reference/predicted masks must be equal 2D arrays")
    reference_pixels = int(np.count_nonzero(ref))
    predicted_pixels = int(np.count_nonzero(pred))
    intersection = int(np.count_nonzero(ref & pred))
    union = int(np.count_nonzero(ref | pred))
    if reference_pixels <= 0:
        raise ValueError("positive-reference scoring requires nonempty reference mask")
    precision = intersection / predicted_pixels if predicted_pixels else 0.0
    recall = intersection / reference_pixels
    iou = intersection / union if union else 0.0
    denom = reference_pixels + predicted_pixels
    dice = 2.0 * intersection / denom if denom else 0.0
    return {
        "reference_pixels": reference_pixels,
        "predicted_pixels": predicted_pixels,
        "intersection_pixels": intersection,
        "union_pixels": union,
        "prediction_precision": float(precision),
        "reference_recall": float(recall),
        "iou": float(iou),
        "dice": float(dice),
    }


def aggregate_positive_rows(frame: pd.DataFrame) -> dict[str, object]:
    required = {
        "reference_pixels",
        "predicted_pixels",
        "intersection_pixels",
        "prediction_precision",
        "reference_recall",
        "iou",
        "dice",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"positive scoring frame missing columns: {missing}")
    if len(frame) != 109:
        raise ValueError(f"expected 109 positive-reference rows, found {len(frame)}")
    reference_total = int(pd.to_numeric(frame["reference_pixels"], errors="raise").sum())
    predicted_total = int(pd.to_numeric(frame["predicted_pixels"], errors="raise").sum())
    intersection_total = int(pd.to_numeric(frame["intersection_pixels"], errors="raise").sum())
    pooled_precision = intersection_total / predicted_total if predicted_total else 0.0
    pooled_recall = intersection_total / reference_total if reference_total else 0.0
    median_precision = float(pd.to_numeric(frame["prediction_precision"], errors="raise").median())
    summary = {
        "positive_reference_images": 109,
        "reference_pixels_total": reference_total,
        "predicted_pixels_total": predicted_total,
        "intersection_pixels_total": intersection_total,
        "pooled_prediction_precision": float(pooled_precision),
        "pooled_reference_recall": float(pooled_recall),
        "median_positive_image_prediction_precision": median_precision,
        "median_positive_image_reference_recall": float(pd.to_numeric(frame["reference_recall"], errors="raise").median()),
        "median_positive_image_iou": float(pd.to_numeric(frame["iou"], errors="raise").median()),
        "median_positive_image_dice": float(pd.to_numeric(frame["dice"], errors="raise").median()),
        "fraction_positive_images_iou_ge_0_25": float((pd.to_numeric(frame["iou"], errors="raise") >= 0.25).mean()),
        "fraction_positive_images_with_nonempty_prediction": float((pd.to_numeric(frame["predicted_pixels"], errors="raise") > 0).mean()),
    }
    summary["gate_components"] = {
        "pooled_prediction_precision_ge_0_70": bool(pooled_precision >= 0.70),
        "pooled_reference_recall_ge_0_35": bool(pooled_recall >= 0.35),
        "median_positive_image_prediction_precision_ge_0_70": bool(median_precision >= 0.70),
    }
    summary["limited_gate_pass"] = bool(all(summary["gate_components"].values()))
    return summary


def _load_contract(path: Path) -> dict:
    contract = json.loads(path.read_text(encoding="utf-8"))
    expected = "prospectively_frozen_after_complete_metadata_intake_before_any_monarda_pixel_decode_or_model_execution"
    if contract.get("status") != expected:
        raise RuntimeError("Monarda region-agreement contract is not frozen pre-outcome")
    if int(contract["fixed_denominator"]["all_images"]) != 110:
        raise RuntimeError("Monarda denominator drifted")
    if contract["limited_gate"]["requirements"] != {
        "complete_reference_alignment": True,
        "all_110_images_accounted_for": True,
        "pooled_prediction_precision_minimum": 0.70,
        "pooled_reference_recall_minimum": 0.35,
        "median_positive_image_prediction_precision_minimum": 0.70,
    }:
        raise RuntimeError("Monarda limited gate drifted")
    return contract


def _reference_census(archive: zipfile.ZipFile) -> list[dict]:
    rows: list[dict] = []
    for split in SPLIT_ORDER:
        coco = json.loads(archive.read(f"{split}/_annotations.coco.json"))
        images = {int(row["id"]): row for row in coco["images"]}
        annotations: dict[int, list[dict]] = defaultdict(list)
        for ann in coco["annotations"]:
            if int(ann["category_id"]) != 1:
                raise RuntimeError("unexpected Monarda annotation category")
            annotations[int(ann["image_id"])].append(ann)
        for image_id in sorted(images):
            image = images[image_id]
            segmentations: list[list[float]] = []
            for ann in annotations[image_id]:
                seg = ann.get("segmentation")
                if not isinstance(seg, list):
                    raise RuntimeError("Monarda diagnostic supports polygon segmentation only")
                segmentations.extend(seg)
            rows.append(
                {
                    "split": split,
                    "coco_image_id": image_id,
                    "member_path": f"{split}/{image['file_name']}",
                    "width_declared": int(image["width"]),
                    "height_declared": int(image["height"]),
                    "annotation_count": len(annotations[image_id]),
                    "segmentations": segmentations,
                }
            )
    if len(rows) != 110:
        raise RuntimeError(f"Monarda reference census drifted: {len(rows)}")
    return rows


def _alignment_audit(archive: zipfile.ZipFile, census: list[dict]) -> pd.DataFrame:
    rows = []
    for row in census:
        status = "alignment_pass"
        decoded_width = decoded_height = None
        exif_orientation = None
        reference_pixels = None
        try:
            raw = archive.read(row["member_path"])
            with Image.open(BytesIO(raw)) as image:
                image.load()
                decoded_width, decoded_height = image.size
                exif_orientation = image.getexif().get(EXIF_ORIENTATION_TAG, 1)
            if (decoded_width, decoded_height) != (row["width_declared"], row["height_declared"]):
                status = "declared_dimension_mismatch"
            elif exif_orientation not in (None, 1):
                status = "unsupported_exif_orientation"
            else:
                mask = rasterize_polygon_union(
                    row["width_declared"], row["height_declared"], row["segmentations"]
                )
                reference_pixels = int(np.count_nonzero(mask))
                if row["annotation_count"] > 0 and reference_pixels <= 0:
                    status = "positive_annotations_rasterized_empty"
                if row["annotation_count"] == 0 and reference_pixels != 0:
                    status = "zero_annotation_rasterized_nonempty"
        except Exception as exc:
            status = f"alignment_runtime_failure:{type(exc).__name__}"
        rows.append(
            {
                "split": row["split"],
                "coco_image_id": row["coco_image_id"],
                "member_path": row["member_path"],
                "annotation_count": row["annotation_count"],
                "width_declared": row["width_declared"],
                "height_declared": row["height_declared"],
                "decoded_width": decoded_width,
                "decoded_height": decoded_height,
                "exif_orientation": exif_orientation,
                "reference_pixels": reference_pixels,
                "alignment_status": status,
            }
        )
    return pd.DataFrame(rows)


def _write_not_evaluable(output_dir: Path, audit: pd.DataFrame, contract: dict, reason: str, archive_sha: str) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / "monarda_region_agreement_rows_v1.csv", index=False, lineterminator="\n")
    summary = {
        "protocol": contract["protocol"],
        "status": "not_evaluable_reference_alignment",
        "reason": reason,
        "all_110_images_accounted_for": bool(len(audit) == 110),
        "model_executed": False,
        "ecological_results_changed": False,
        "archive_sha256": archive_sha,
    }
    (output_dir / "monarda_region_agreement_result_v1.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--detector-weight", type=Path, required=True)
    ap.add_argument("--efficient-sam-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    ap.add_argument("--intake", type=Path, default=DEFAULT_INTAKE)
    ap.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    ap.add_argument("--roi-contract", type=Path, default=DEFAULT_ROI_CONTRACT)
    ap.add_argument("--locked-result", type=Path, default=DEFAULT_LOCKED)
    ap.add_argument("--torch-threads", type=int, default=2)
    args = ap.parse_args()

    contract = _load_contract(args.contract)
    if sha256_file(args.archive) != contract["parent_checkpoint"]["archive_sha256"]:
        raise RuntimeError("Monarda archive SHA differs from frozen contract")
    if sha256_file(args.intake) != contract["parent_checkpoint"]["archive_intake_receipt_sha256"]:
        raise RuntimeError("Monarda intake receipt SHA differs from frozen contract")
    if sha256_file(args.mapping) != contract["parent_checkpoint"]["source_mapping_receipt_sha256"]:
        raise RuntimeError("Monarda source-mapping receipt SHA differs from frozen contract")
    # Re-run the complete metadata/CRC intake before any pixel decode.
    intake = inspect_archive(args.archive)
    if intake["total_images"] != 110 or intake["total_annotations"] != 788:
        raise RuntimeError("Monarda archive intake census drifted")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.archive) as archive:
        census = _reference_census(archive)
        alignment = _alignment_audit(archive, census)
        if len(alignment) != 110 or not alignment["alignment_status"].eq("alignment_pass").all():
            return _write_not_evaluable(
                output_dir,
                alignment,
                contract,
                "one_or_more_reference_images_failed_exact_decode_dimension_orientation_or_rasterization_gate",
                contract["parent_checkpoint"]["archive_sha256"],
            )
        positive = alignment["annotation_count"].astype(int).gt(0)
        if int(positive.sum()) != 109 or int((~positive).sum()) != 1:
            raise RuntimeError("Monarda positive/unknown reference census drifted")

        # Heavy runtime imports only after the entire reference/alignment gate passes.
        from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
        from fcp_pipeline.flower_roi_v4_runtime import (
            FrozenFlowerColourEstimator,
            file_sha256,
            validate_scaleout_authorization,
        )

        roi_contract = json.loads(args.roi_contract.read_text(encoding="utf-8"))
        locked = json.loads(args.locked_result.read_text(encoding="utf-8"))
        validate_roi_v4_contract(roi_contract)
        detector_sha = file_sha256(args.detector_weight)
        if detector_sha != contract["frozen_fcp_runtime"]["detector_weight_sha256"]:
            raise RuntimeError("detector weight SHA differs from frozen Monarda contract")
        validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)
        for name, path, expected_sha in (
            (
                "encoder",
                args.efficient_sam_dir / Path(roi_contract["mask_generator"]["encoder_path"]).name,
                contract["frozen_fcp_runtime"]["efficient_sam_encoder_sha256"],
            ),
            (
                "decoder",
                args.efficient_sam_dir / Path(roi_contract["mask_generator"]["decoder_path"]).name,
                contract["frozen_fcp_runtime"]["efficient_sam_decoder_sha256"],
            ),
        ):
            if file_sha256(path) != expected_sha:
                raise RuntimeError(f"EfficientSAM {name} SHA differs from frozen Monarda contract")
        estimator = FrozenFlowerColourEstimator(
            args.detector_weight,
            args.efficient_sam_dir,
            roi_contract,
            torch_threads=int(args.torch_threads),
        )

        scored_rows = []
        census_key = {(r["split"], r["coco_image_id"]): r for r in census}
        for record in alignment.to_dict("records"):
            key = (str(record["split"]), int(record["coco_image_id"]))
            source_row = census_key[key]
            reference_mask = rasterize_polygon_union(
                source_row["width_declared"], source_row["height_declared"], source_row["segmentations"]
            )
            raw = archive.read(source_row["member_path"])
            model_status = "model_success"
            retained_instances = 0
            automated_status = ""
            predicted_mask = np.zeros_like(reference_mask, dtype=bool)
            try:
                with Image.open(BytesIO(raw)) as image:
                    image.load()
                    measured = estimator.measure(image)
                predicted_mask = np.asarray(measured["flower_mask"], dtype=bool)
                if predicted_mask.shape != reference_mask.shape:
                    raise RuntimeError("predicted mask shape differs from frozen reference frame")
                retained_instances = int(measured.get("retained_instances") or 0)
                automated_status = str(measured.get("automated_colour_state_status") or "")
            except Exception as exc:
                model_status = f"model_runtime_failure:{type(exc).__name__}"
                predicted_mask = np.zeros_like(reference_mask, dtype=bool)

            row = {
                "split": key[0],
                "coco_image_id": key[1],
                "member_path": source_row["member_path"],
                "annotation_count": int(source_row["annotation_count"]),
                "reference_role": "positive_generic_flower_region" if source_row["annotation_count"] > 0 else "reference_unknown_empty",
                "model_status": model_status,
                "retained_instances": retained_instances,
                "automated_colour_state_status_descriptive_only": automated_status,
            }
            if source_row["annotation_count"] > 0:
                row.update(score_masks(reference_mask, predicted_mask))
            else:
                row.update(
                    {
                        "reference_pixels": 0,
                        "predicted_pixels": int(np.count_nonzero(predicted_mask)),
                        "intersection_pixels": None,
                        "union_pixels": None,
                        "prediction_precision": None,
                        "reference_recall": None,
                        "iou": None,
                        "dice": None,
                    }
                )
            scored_rows.append(row)
            print(f"monarda_region_scored={len(scored_rows)}/110", flush=True)

    scored = pd.DataFrame(scored_rows)
    if len(scored) != 110:
        raise RuntimeError("Monarda scored row census is not 110")
    positive_frame = scored.loc[scored["reference_role"].eq("positive_generic_flower_region")].copy()
    aggregate = aggregate_positive_rows(positive_frame)
    split_descriptive = {}
    for split in SPLIT_ORDER:
        part = positive_frame.loc[positive_frame["split"].eq(split)].copy()
        split_descriptive[split] = {
            "positive_reference_images": int(len(part)),
            "median_iou": float(pd.to_numeric(part["iou"], errors="raise").median()),
            "median_precision": float(pd.to_numeric(part["prediction_precision"], errors="raise").median()),
            "median_recall": float(pd.to_numeric(part["reference_recall"], errors="raise").median()),
        }

    scored_path = output_dir / "monarda_region_agreement_rows_v1.csv"
    scored.to_csv(scored_path, index=False, lineterminator="\n")
    unknown = scored.loc[scored["reference_role"].eq("reference_unknown_empty")].iloc[0]
    summary = {
        "protocol": contract["protocol"],
        "status": "complete_monarda_limited_region_agreement_diagnostic",
        "inferential_role": contract["role"],
        "all_110_images_accounted_for": True,
        "reference_alignment_complete": True,
        "positive_reference_images": 109,
        "reference_unknown_images": 1,
        "unknown_zero_annotation_image": {
            "split": str(unknown["split"]),
            "coco_image_id": int(unknown["coco_image_id"]),
            "predicted_pixels_descriptive_only": int(unknown["predicted_pixels"]),
            "model_status": str(unknown["model_status"]),
            "counted_as_verified_negative": False,
        },
        "aggregate": aggregate,
        "limited_gate_pass": bool(aggregate["limited_gate_pass"]),
        "provider_split_descriptive_only": split_descriptive,
        "model_runtime_failure_images": int(scored["model_status"].astype(str).str.startswith("model_runtime_failure").sum()),
        "ecological_results_changed": False,
        "reserve_flower_specific_gate_reclassified": False,
        "claim_ceiling": contract["claim_ceiling"],
        "limits": [
            "Generic flowers polygons are not verified focal-taxon petal truth.",
            "The released 110-image selection is not reproducibly explained by the author's unseeded notebook sampling step.",
            "Provider train/valid/test labels are not treated as FCP holdouts.",
            "Zero photo-ID/observation-ID overlap with the two FCP frames does not establish event, observer, JRC or foundation-model independence.",
            "A pass is measurement-localization evidence only and cannot rescue the reserve flower-minus-background p=0.087 ecological gate."
        ],
        "lineage": {
            "archive_sha256": sha256_file(args.archive),
            "contract_sha256": sha256_file(args.contract),
            "intake_receipt_sha256": sha256_file(args.intake),
            "source_mapping_receipt_sha256": sha256_file(args.mapping),
            "detector_weight_sha256": detector_sha,
            "efficient_sam_encoder_sha256": contract["frozen_fcp_runtime"]["efficient_sam_encoder_sha256"],
            "efficient_sam_decoder_sha256": contract["frozen_fcp_runtime"]["efficient_sam_decoder_sha256"],
            "rows_sha256": sha256_file(scored_path),
        },
    }
    result_path = output_dir / "monarda_region_agreement_result_v1.json"
    result_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
