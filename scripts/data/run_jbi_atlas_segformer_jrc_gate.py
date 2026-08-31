#!/usr/bin/env python3
"""Run the frozen SegFormer v3 smoke, development, or locked JRC field gate."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any, Mapping, Sequence

import numpy as np
from PIL import Image, ImageOps


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.segformer_roi import (
    CANVAS_SIZE,
    class_masks_from_labels,
    evaluate_flip_stable_admission,
    score_jrc_boxes,
    summarize_jrc_gate,
    validate_jrc_box_edge_amendment,
    validate_jrc_box_edge_amendment_v2,
    validate_roi_v3_contract,
)


SOURCE_INVENTORY = Path(
    "data/atlas/qualification/roi_v3_sources/jrc_flower_detection_source_inventory_v1.csv"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write empty JRC gate rows")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def load_model(model_dir: Path, contract: Mapping[str, Any], torch_threads: int):
    estimator = contract["estimator"]
    expected_files = {
        "model.safetensors": estimator["model_safetensors_sha256"],
        "config.json": estimator["config_sha256"],
        "preprocessor_config.json": estimator["preprocessor_config_sha256"],
    }
    for name, expected in expected_files.items():
        path = model_dir / name
        if sha256(path) != expected:
            raise RuntimeError(f"pinned SegFormer file hash mismatch: {name}")
    if (model_dir / "model.safetensors").stat().st_size != int(
        estimator["model_safetensors_bytes"]
    ):
        raise RuntimeError("pinned SegFormer weight byte count changed")
    os.environ.setdefault("USE_TF", "0")
    os.environ.setdefault("TRANSFORMERS_NO_TF", "1")
    import torch
    from transformers import SegformerForSemanticSegmentation, SegformerImageProcessor

    torch.set_num_threads(torch_threads)
    processor = SegformerImageProcessor.from_pretrained(model_dir, local_files_only=True)
    model = SegformerForSemanticSegmentation.from_pretrained(
        model_dir, local_files_only=True
    ).eval()
    return torch, processor, model


def infer_labels(image: Image.Image, torch, processor, model) -> np.ndarray:
    inputs = processor(images=image, return_tensors="pt")
    with torch.inference_mode():
        logits = model(**inputs).logits
        resized = torch.nn.functional.interpolate(
            logits,
            size=(CANVAS_SIZE, CANVAS_SIZE),
            mode="bilinear",
            align_corners=False,
        )
        labels = resized.argmax(dim=1)[0]
    return labels.cpu().numpy().astype(np.int16, copy=False)


def load_annotations(path: Path, expected_sha256: str) -> dict[str, Any]:
    if sha256(path) != expected_sha256:
        raise RuntimeError(f"JRC annotation SHA-256 mismatch: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("categories") != [{"id": 1, "name": "flower", "supercategory": ""}]:
        raise RuntimeError("JRC category identity changed")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_roi_estimator_contract_v3.json"),
    )
    parser.add_argument(
        "--box-edge-amendment",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_roi_v3_jrc_box_edge_amendment_v2.json"),
    )
    parser.add_argument("--phase", choices=("smoke", "development", "locked_test"), required=True)
    parser.add_argument("--jrc-root", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--development-result", type=Path)
    parser.add_argument("--torch-threads", type=int, default=8)
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_roi_v3_contract(contract)
    amendment = json.loads(args.box_edge_amendment.read_text(encoding="utf-8"))
    validate_jrc_box_edge_amendment_v2(amendment)
    contract_hash = canonical_sha256(args.contract)
    amendment_hash = canonical_sha256(args.box_edge_amendment)
    if args.phase == "locked_test":
        if args.limit is not None:
            raise RuntimeError("locked JRC test cannot be limited")
        if args.development_result is None:
            raise RuntimeError("locked JRC test requires the committed development result")
        development = json.loads(args.development_result.read_text(encoding="utf-8"))
        if (
            development.get("status") != "pass_jrc_development"
            or development.get("contract_sha256_lf_canonical_v1") != contract_hash
            or development.get("box_edge_amendment_sha256_lf_canonical_v1")
            != amendment_hash
        ):
            raise RuntimeError("JRC development gate is not passed for this contract")

    split = "test" if args.phase == "locked_test" else "train"
    expected_annotation = contract["jrc_field_gate"][f"{split}_annotation_sha256"]
    annotation_path = args.jrc_root / "annotations" / f"{split}.json"
    data = load_annotations(annotation_path, expected_annotation)
    inventory = {
        (row["split"], int(row["image_id"])): row for row in read_csv(SOURCE_INVENTORY)
    }
    annotations: dict[int, list[list[float]]] = defaultdict(list)
    for row in data["annotations"]:
        if int(row["category_id"]) != 1:
            raise RuntimeError("non-flower object reached the JRC gate")
        annotations[int(row["image_id"])].append([float(value) for value in row["bbox"]])
    images = sorted(data["images"], key=lambda row: int(row["id"]))
    if args.phase == "smoke":
        images = images[:1]
    elif args.limit is not None:
        if args.limit < 1:
            raise ValueError("development limit must be positive")
        images = images[: args.limit]

    torch, processor, model = load_model(args.model_dir, contract, args.torch_threads)
    if args.phase == "smoke":
        metadata = images[0]
        image_path = args.jrc_root / "images" / split / str(metadata["file_name"])
        image = ImageOps.exif_transpose(Image.open(image_path)).convert("RGB")
        started = time.perf_counter()
        first = infer_labels(image, torch, processor, model)
        first_seconds = time.perf_counter() - started
        started = time.perf_counter()
        second = infer_labels(image, torch, processor, model)
        second_seconds = time.perf_counter() - started
        result = {
            "protocol": contract["protocol"],
            "phase": "smoke",
            "status": "pass_deterministic_cpu_smoke" if np.array_equal(first, second) else "stop_nondeterministic_cpu_smoke",
            "contract_sha256_lf_canonical_v1": contract_hash,
            "box_edge_amendment_sha256_lf_canonical_v1": amendment_hash,
            "image_id": int(metadata["id"]),
            "image_sha256": sha256(image_path),
            "first_seconds": first_seconds,
            "second_seconds": second_seconds,
            "identical_repeated_labels": bool(np.array_equal(first, second)),
            "flower_pixels": int(np.count_nonzero(first == int(contract["estimator"]["flower_label"]))),
            "scaleout_candidate_pixels_opened": False,
            "jrc_test_images_decoded_or_scored": False,
        }
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_json(args.output_dir / "jrc_segformer_smoke_result.json", result)
        print(json.dumps(result, indent=2, sort_keys=True))
        if not result["identical_repeated_labels"]:
            raise SystemExit(2)
        return

    output_rows: list[dict[str, Any]] = []
    for index, metadata in enumerate(images, start=1):
        image_id = int(metadata["id"])
        source = inventory.get((split, image_id))
        if source is None:
            raise RuntimeError(f"JRC image {split}/{image_id} is absent from source freeze")
        image_path = args.jrc_root / "images" / split / str(metadata["file_name"])
        if sha256(image_path) != source["image_sha256"]:
            raise RuntimeError(f"JRC image hash changed: {split}/{image_id}")
        image = ImageOps.exif_transpose(Image.open(image_path)).convert("RGB")
        if image.size != (int(metadata["width"]), int(metadata["height"])):
            raise RuntimeError(f"JRC oriented dimensions changed: {split}/{image_id}")
        rgb = np.asarray(
            image.resize((CANVAS_SIZE, CANVAS_SIZE), Image.Resampling.BILINEAR),
            dtype=np.uint8,
        )
        labels = infer_labels(image, torch, processor, model)
        flipped_image = image.transpose(Image.Transpose.FLIP_LEFT_RIGHT)
        flipped_labels = infer_labels(flipped_image, torch, processor, model)
        flower, plant = class_masks_from_labels(
            labels,
            flower_label=int(contract["estimator"]["flower_label"]),
            plant_label=int(contract["estimator"]["plant_background_control_label"]),
        )
        flipped_flower, _ = class_masks_from_labels(
            flipped_labels,
            flower_label=int(contract["estimator"]["flower_label"]),
            plant_label=int(contract["estimator"]["plant_background_control_label"]),
        )
        admission = evaluate_flip_stable_admission(
            rgb,
            flower,
            plant,
            np.fliplr(flipped_flower),
            contract,
        )
        box_score = score_jrc_boxes(
            flower,
            annotations[image_id],
            source_width=int(metadata["width"]),
            source_height=int(metadata["height"]),
        )
        output_rows.append(
            {
                "split": split,
                "image_id": image_id,
                "file_name": str(metadata["file_name"]),
                "image_sha256": source["image_sha256"],
                **admission,
                **box_score,
            }
        )
        if index % 10 == 0 or index == len(images):
            print(f"scored={index}/{len(images)} phase={args.phase}", flush=True)

    summary = summarize_jrc_gate(output_rows, contract, phase=args.phase)
    rows_path = args.output_dir / f"jrc_segformer_{args.phase}_rows.csv"
    write_csv(rows_path, output_rows)
    result = {
        **summary,
        "contract_sha256_lf_canonical_v1": contract_hash,
        "box_edge_amendment_sha256_lf_canonical_v1": amendment_hash,
        "model_id": contract["estimator"]["model_id"],
        "model_revision": contract["estimator"]["revision"],
        "model_safetensors_sha256": contract["estimator"]["model_safetensors_sha256"],
        "source_annotation_sha256": expected_annotation,
        "rows_sha256": sha256(rows_path),
        "scaleout_candidate_pixels_opened": False,
        "jrc_test_images_decoded_or_scored": args.phase == "locked_test",
    }
    result_path = args.output_dir / f"jrc_segformer_{args.phase}_result.json"
    write_json(result_path, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["status"].startswith("pass_"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
