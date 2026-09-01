#!/usr/bin/env python3
"""Measure one terminal atlas shard under the v5 inference firewall and locked ROI v4."""

from __future__ import annotations

import argparse
import csv
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

from fcp_pipeline.atlas_measurement import select_measurement_shard, validate_measurement_result_rows
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
from fcp_pipeline.flower_roi_v4_runtime import (
    FrozenFlowerColourEstimator,
    file_sha256,
    validate_scaleout_authorization,
)
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs
from scripts.data.measure_jbi_atlas_blinded_images_v4 import (
    MODEL_ID,
    SUMMARY_FIELDS,
    canonical_sha256,
    failed_record,
    read_csv,
)
from scripts.data.validate_jbi_atlas_roi_v4_gate_evidence import load_committed_locked_scaleout_result


DEFAULT_MEASUREMENT_CONTRACT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
DEFAULT_INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"
DEFAULT_ROI_CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-manifest", type=Path, required=True)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--trained-weight", type=Path, required=True)
    parser.add_argument("--roi-evidence-dir", type=Path, required=True)
    parser.add_argument("--efficient-sam-weights-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--roi-contract", type=Path, default=DEFAULT_ROI_CONTRACT)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT_CONTRACT)
    parser.add_argument("--inference-v5", type=Path, default=DEFAULT_INFERENCE)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--torch-threads", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inference = load_json(args.inference_v5)
    execution = load_json(args.measurement_contract)
    validate_measurement_execution_contract(execution, inference)
    verify_repo_parent_blobs(execution)
    technical = execution["technical_execution"]
    if (
        args.shard_count != technical["measurement_shard_count"]
        or args.torch_threads != technical["measurement_torch_threads_per_worker"]
    ):
        raise ValueError("v5 measurement worker settings differ from the frozen contract")

    roi_contract = load_json(args.roi_contract)
    validate_roi_v4_contract(roi_contract)
    trained_weight_sha = file_sha256(args.trained_weight)
    locked = load_committed_locked_scaleout_result(args.roi_evidence_dir)
    validate_scaleout_authorization(locked, trained_weight_sha256=trained_weight_sha)

    selected = select_measurement_shard(
        read_csv(args.measurement_manifest),
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )
    if not selected:
        raise RuntimeError("v5 measurement shard is empty")

    estimator = FrozenFlowerColourEstimator(
        args.trained_weight,
        args.efficient_sam_weights_dir,
        roi_contract,
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
                row, "image_file_missing", trained_weight_sha256=trained_weight_sha
            )
            continue
        try:
            image_hash = file_sha256(image_path)
            if cache_path.is_file():
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                if (
                    cached.get("image_sha256") == image_hash
                    and cached.get("model_revision") == trained_weight_sha
                    and cached.get("contract_sha256_lf_canonical_v1") == contract_hash
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
            "background_features_available": bool(measured["background_features_available"]),
            "automated_colour_state_status": measured["automated_colour_state_status"],
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
    result_path = args.output_dir / f"measurement_v5_shard_{args.shard_index:04d}.csv"
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
        "status": "complete_location_blind_roi_v4_measurement_v5_shard",
        "protocol": execution["protocol"],
        "inference_version": inference["version"],
        "superseded_v3_ordered_inference_used": False,
        "roi_protocol": roi_contract["protocol"],
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
    (args.output_dir / f"measurement_v5_shard_{args.shard_index:04d}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
