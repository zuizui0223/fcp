#!/usr/bin/env python3
"""Measure one location-blind atlas shard with the locked ROI-v4 estimator."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any

from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import (
    select_measurement_shard,
    validate_inference_contract,
    validate_measurement_result_rows,
)
from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
from fcp_pipeline.flower_roi_v4_runtime import (
    FrozenFlowerColourEstimator,
    file_sha256,
    validate_scaleout_authorization,
)


MODEL_ID = "jbi-atlas-roi-estimator-v4"
SUMMARY_FIELDS = [
    *[
        f"{prefix}_{channel}_{metric}"
        for prefix in ("flower", "background")
        for channel in ("L", "a", "b")
        for metric in ("mean", "sd", "q10", "q50", "q90")
    ],
    "flower_effective_pixels",
    "background_effective_pixels",
    "flip_effective_pixels",
    "flip_background_pixels",
    "horizontal_flip_mask_iou",
    "horizontal_flip_colour_delta_e",
    "detector_predictions",
    "retained_instances",
    "flip_detector_predictions",
    "flip_retained_instances",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def failed_record(
    row: dict[str, str], reason: str, *, trained_weight_sha256: str
) -> dict[str, Any]:
    return {
        "measurement_id": row["measurement_id"],
        "species_blind_id": row["species_blind_id"],
        "image_sha256": "",
        "model_id": MODEL_ID,
        "model_revision": trained_weight_sha256,
        **{field: None for field in SUMMARY_FIELDS},
        "background_features_available": False,
        "automated_colour_state_status": "image_acquisition_failed",
        "failure_reasons": reason[:500],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-manifest", type=Path, required=True)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--trained-weight", type=Path, required=True)
    parser.add_argument("--roi-result", type=Path, required=True)
    parser.add_argument("--efficient-sam-weights-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--roi-contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"),
    )
    parser.add_argument(
        "--inference-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json"),
    )
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--torch-threads", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inference = json.loads(args.inference_contract.read_text(encoding="utf-8"))
    validate_inference_contract(inference)
    contract = json.loads(args.roi_contract.read_text(encoding="utf-8"))
    validate_roi_v4_contract(contract)
    trained_weight_sha = file_sha256(args.trained_weight)
    locked = json.loads(args.roi_result.read_text(encoding="utf-8"))
    validate_scaleout_authorization(
        locked, trained_weight_sha256=trained_weight_sha
    )
    selected = select_measurement_shard(
        read_csv(args.measurement_manifest),
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )
    if not selected:
        raise RuntimeError("measurement shard is empty; reduce shard_count")
    estimator = FrozenFlowerColourEstimator(
        args.trained_weight,
        args.efficient_sam_weights_dir,
        contract,
        torch_threads=args.torch_threads,
    )
    contract_hash = canonical_sha256(args.roi_contract)
    cache_dir = args.output_dir / "photo_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict[str, Any]] = {}
    started = time.time()
    for index, row in enumerate(selected, start=1):
        image_path = args.images_dir / row["image_filename"]
        cache_path = cache_dir / f"{row['measurement_id']}.json"
        if not image_path.is_file():
            records[row["measurement_id"]] = failed_record(
                row,
                "image_file_missing",
                trained_weight_sha256=trained_weight_sha,
            )
            continue
        try:
            image_hash = file_sha256(image_path)
            if cache_path.is_file():
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                if (
                    cached.get("image_sha256") == image_hash
                    and cached.get("model_revision") == trained_weight_sha
                    and cached.get("contract_sha256_lf_canonical_v1")
                    == contract_hash
                ):
                    records[row["measurement_id"]] = cached
                    continue
            with Image.open(image_path) as image:
                measured = estimator.measure(image)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            records[row["measurement_id"]] = failed_record(
                row,
                f"image_measurement_failed:{type(exc).__name__}",
                trained_weight_sha256=trained_weight_sha,
            )
            continue
        record = {
            "measurement_id": row["measurement_id"],
            "species_blind_id": row["species_blind_id"],
            "image_sha256": image_hash,
            "model_id": MODEL_ID,
            "model_revision": trained_weight_sha,
            **{field: measured[field] for field in SUMMARY_FIELDS if field in measured},
            "detector_predictions": len(measured["predictions"]),
            "retained_instances": int(measured["retained_instances"]),
            "flip_detector_predictions": int(measured["flip_detector_predictions"]),
            "flip_retained_instances": int(measured["flip_retained_instances"]),
            "background_features_available": bool(
                measured["background_features_available"]
            ),
            "automated_colour_state_status": measured[
                "automated_colour_state_status"
            ],
            "failure_reasons": measured["failure_reasons"],
            "contract_sha256_lf_canonical_v1": contract_hash,
        }
        records[row["measurement_id"]] = record
        cache_path.write_text(
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        if index % 10 == 0 or index == len(selected):
            print(f"measured_or_failed={index}/{len(selected)}", flush=True)

    ordered = [records[row["measurement_id"]] for row in selected]
    validate_measurement_result_rows(ordered)
    result_path = args.output_dir / f"measurement_shard_{args.shard_index:04d}.csv"
    fields = list(ordered[0])
    union = {key for row in ordered for key in row}
    fields.extend(sorted(union - set(fields)))
    with result_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(ordered)

    import onnxruntime
    import torch
    import ultralytics

    manifest = {
        "status": "complete_location_blind_roi_v4_measurement_shard",
        "protocol": inference["protocol"],
        "roi_protocol": contract["protocol"],
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "frozen_shard_denominator": len(selected),
        "terminal_records": len(ordered),
        "coordinates_opened": False,
        "taxon_names_opened": False,
        "result_sha256": file_sha256(result_path),
        "model_id": MODEL_ID,
        "trained_weight_sha256": trained_weight_sha,
        "roi_contract_sha256_lf_canonical_v1": contract_hash,
        "elapsed_seconds_new_inference": time.time() - started,
        "environment": {
            "python": platform.python_version(),
            "torch": torch.__version__,
            "ultralytics": ultralytics.__version__,
            "onnxruntime": onnxruntime.__version__,
        },
    }
    (args.output_dir / f"measurement_shard_{args.shard_index:04d}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
