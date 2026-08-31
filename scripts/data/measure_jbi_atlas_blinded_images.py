#!/usr/bin/env python3
"""Measure one location-blind atlas image shard with the pinned estimator."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
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
from scripts.data.extract_inaturalist_automated_colour_states import (
    MODEL_ID,
    MODEL_REVISION,
    contract_sha256,
    image_array,
    inference_batch,
    json_safe,
    summarize_photo,
)


MODEL_SHA256 = "d00ca85d6b859f9d07b7cfb8ef26fe9771cb275b34c9368f2ecf603139307f55"
FEATURE_FIELDS = [
    *[
        f"{prefix}_{channel}_{metric}"
        for prefix in ("flower", "background")
        for channel in ("L", "a", "b")
        for metric in ("mean", "sd", "q10", "q50", "q90")
    ],
    "flower_effective_pixels",
    "background_effective_pixels",
    "valid_positive_prompts",
    "prompt_max_delta_e",
    "flip_delta_e",
    "flip_soft_iou",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def failed_record(row: dict[str, str], reason: str) -> dict[str, Any]:
    return {
        "measurement_id": row["measurement_id"],
        "species_blind_id": row["species_blind_id"],
        "image_sha256": "",
        "model_id": MODEL_ID,
        "model_revision": MODEL_REVISION,
        **{field: None for field in FEATURE_FIELDS},
        "background_features_available": False,
        "automated_colour_state_status": "image_acquisition_failed",
        "failure_reasons": reason[:500],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-manifest", type=Path, required=True)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--inference-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json"),
    )
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--torch-threads", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inference = json.loads(args.inference_contract.read_text(encoding="utf-8"))
    validate_inference_contract(inference)
    selected = select_measurement_shard(
        read_csv(args.measurement_manifest),
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )
    if not selected:
        raise RuntimeError("measurement shard is empty; reduce shard_count")
    model_path = args.model_dir / "model.safetensors"
    if sha256(model_path) != MODEL_SHA256:
        raise RuntimeError("pinned CLIPSeg model SHA-256 mismatch")
    estimator_contract = contract_sha256(MODEL_SHA256)
    cache_dir = args.output_dir / "photo_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict[str, Any]] = {}
    pending: list[tuple[dict[str, str], Image.Image, str]] = []
    for row in selected:
        image_path = args.images_dir / row["image_filename"]
        cache_path = cache_dir / f"{row['measurement_id']}.json"
        if not image_path.is_file():
            records[row["measurement_id"]] = failed_record(row, "image_file_missing")
            continue
        try:
            image_hash = sha256(image_path)
            if cache_path.is_file():
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                if (
                    cached.get("image_sha256") == image_hash
                    and cached.get("contract_sha256") == estimator_contract
                ):
                    records[row["measurement_id"]] = cached
                    continue
            image = Image.open(image_path).convert("RGB")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            records[row["measurement_id"]] = failed_record(
                row, f"image_decode_failed:{type(exc).__name__}"
            )
            continue
        pending.append((row, image, image_hash))

    os.environ.setdefault("USE_TF", "0")
    os.environ.setdefault("TRANSFORMERS_NO_TF", "1")
    import torch
    import transformers
    from transformers import CLIPSegForImageSegmentation, CLIPSegProcessor

    torch.set_num_threads(args.torch_threads)
    processor = CLIPSegProcessor.from_pretrained(
        args.model_dir, local_files_only=True, use_fast=False
    )
    model = CLIPSegForImageSegmentation.from_pretrained(
        args.model_dir, local_files_only=True
    ).eval()
    started = time.time()
    for offset in range(0, len(pending), args.batch_size):
        batch = pending[offset : offset + args.batch_size]
        images = [item[1] for item in batch]
        flipped = [image.transpose(Image.Transpose.FLIP_LEFT_RIGHT) for image in images]
        original_logits = inference_batch(model, processor, images)
        flipped_logits = inference_batch(model, processor, flipped)
        for index, (row, image, image_hash) in enumerate(batch):
            summary = summarize_photo(
                image_array(image),
                original_logits[index],
                image_array(flipped[index]),
                flipped_logits[index],
            )
            record = {
                "measurement_id": row["measurement_id"],
                "species_blind_id": row["species_blind_id"],
                "image_sha256": image_hash,
                "model_id": MODEL_ID,
                "model_revision": MODEL_REVISION,
                **summary,
                "contract_sha256": estimator_contract,
            }
            record = json_safe(record)
            records[row["measurement_id"]] = record
            (cache_dir / f"{row['measurement_id']}.json").write_text(
                json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
        print(
            f"measured={min(offset + len(batch), len(pending))}/{len(pending)}",
            flush=True,
        )

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
    manifest = {
        "status": "complete_location_blind_measurement_shard",
        "protocol": inference["protocol"],
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "frozen_shard_denominator": len(selected),
        "terminal_records": len(ordered),
        "coordinates_opened": False,
        "taxon_names_opened": False,
        "result_sha256": sha256(result_path),
        "model_id": MODEL_ID,
        "model_revision": MODEL_REVISION,
        "model_sha256": MODEL_SHA256,
        "elapsed_seconds_new_inference": time.time() - started,
        "environment": {
            "python": platform.python_version(),
            "torch": torch.__version__,
            "transformers": transformers.__version__,
        },
    }
    (args.output_dir / f"measurement_shard_{args.shard_index:04d}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
