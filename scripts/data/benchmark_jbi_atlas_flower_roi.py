#!/usr/bin/env python3
"""Run the frozen Oxford-17 flower-tissue benchmark for the atlas estimator."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np
from PIL import Image

# Direct execution sets sys.path[0] to scripts/data rather than the repository root.
# Keep the CLI behaviour identical on local Windows and GitHub's Linux runner.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fcp_pipeline.atlas_expansion import validate_expansion_contract
from fcp_pipeline.roi_benchmark import score_flower_weights, summarize_roi_benchmark
from scripts.data.extract_inaturalist_automated_colour_states import (
    MODEL_SIZE,
    candidate_weights,
    image_array,
    inference_batch,
    summarize_photo,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json"),
    )
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--trimaps-dir", type=Path, required=True)
    parser.add_argument("--images-archive", type=Path, required=True)
    parser.add_argument("--trimaps-archive", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--torch-threads", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_expansion_contract(contract)
    benchmark = contract["estimator_qualification"]["independent_roi_benchmark"]
    estimator = contract["estimator_qualification"]["estimator"]

    if sha256(args.images_archive) != benchmark["images_sha256"]:
        raise RuntimeError("Oxford-17 image archive SHA-256 mismatch")
    if sha256(args.trimaps_archive) != benchmark["trimaps_sha256"]:
        raise RuntimeError("Oxford-17 trimap archive SHA-256 mismatch")
    model_path = args.model_dir / "model.safetensors"
    if sha256(model_path) != estimator["model_safetensors_sha256"]:
        raise RuntimeError("pinned CLIPSeg model SHA-256 mismatch")

    trimaps = sorted(args.trimaps_dir.glob("image_*.png"))
    image_pairs: list[tuple[Path, Path]] = []
    for trimap in trimaps:
        image_path = args.images_dir / f"{trimap.stem}.jpg"
        if not image_path.is_file():
            raise RuntimeError(f"missing Oxford image for {trimap.name}")
        image_pairs.append((image_path, trimap))
    if len(image_pairs) < int(benchmark["minimum_scored_images"]):
        raise RuntimeError(
            f"official trimap set is incomplete: found {len(image_pairs)} matching images"
        )

    os.environ.setdefault("USE_TF", "0")
    os.environ.setdefault("TRANSFORMERS_NO_TF", "1")
    import torch
    from transformers import CLIPSegForImageSegmentation, CLIPSegProcessor

    torch.set_num_threads(args.torch_threads)
    processor = CLIPSegProcessor.from_pretrained(
        args.model_dir, local_files_only=True, use_fast=False
    )
    model = CLIPSegForImageSegmentation.from_pretrained(
        args.model_dir, local_files_only=True
    ).eval()

    rows: list[dict[str, Any]] = []
    for offset in range(0, len(image_pairs), args.batch_size):
        batch = image_pairs[offset : offset + args.batch_size]
        images = [Image.open(image_path).convert("RGB") for image_path, _ in batch]
        flipped = [image.transpose(Image.Transpose.FLIP_LEFT_RIGHT) for image in images]
        original_logits = inference_batch(model, processor, images)
        flipped_logits = inference_batch(model, processor, flipped)
        for index, (image_path, trimap_path) in enumerate(batch):
            rgb = image_array(images[index])
            _prompt_weights, ensemble, _negative = candidate_weights(original_logits[index])
            trimap = np.asarray(
                Image.open(trimap_path).resize(
                    (MODEL_SIZE, MODEL_SIZE), Image.Resampling.NEAREST
                ),
                dtype=np.uint8,
            )
            metrics = score_flower_weights(
                rgb,
                ensemble,
                trimap,
                foreground_label=int(benchmark["foreground_label"]),
                background_labels=tuple(
                    int(label) for label in benchmark["background_labels"]
                ),
            )
            admission = summarize_photo(
                rgb,
                original_logits[index],
                image_array(flipped[index]),
                flipped_logits[index],
            )
            rows.append(
                {
                    "image_id": image_path.stem,
                    "image_sha256": sha256(image_path),
                    "trimap_sha256": sha256(trimap_path),
                    "estimator_admitted": admission["automated_colour_state_status"]
                    == "automated_colour_state_admitted",
                    "failure_reasons": admission["failure_reasons"],
                    **metrics,
                }
            )
        print(f"scored={min(offset + len(batch), len(image_pairs))}/{len(image_pairs)}", flush=True)

    summary = summarize_roi_benchmark(rows, benchmark)
    result = {
        "protocol": contract["protocol"],
        "benchmark_dataset": benchmark["dataset"],
        "selection": benchmark["selection"],
        "model_id": estimator["model_id"],
        "model_revision": estimator["model_revision"],
        "model_safetensors_sha256": estimator["model_safetensors_sha256"],
        "images_archive_sha256": sha256(args.images_archive),
        "trimaps_archive_sha256": sha256(args.trimaps_archive),
        **summary,
    }
    write_csv(args.output_dir / "oxford17_roi_benchmark_rows.csv", rows)
    write_json(args.output_dir / "oxford17_roi_benchmark_result.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "pass_independent_roi_benchmark":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
