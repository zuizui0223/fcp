#!/usr/bin/env python3
"""Run the frozen composite flower detector and box-prompt mask JRC gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
from PIL import Image, ImageOps
from skimage.color import rgb2lab


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4 import (
    CANVAS_SIZE,
    box_to_canvas,
    greedy_detection_matches,
    letterbox_geometry,
    select_prompt_mask,
    summarize_composite_gate,
    validate_reference_size_amendment,
    validate_roi_v4_contract,
)


CONTRACT = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
AMENDMENT = ROOT / "docs/supporting/jbi_atlas_roi_v4_reference_size_amendment_v1.json"
SOURCE_INVENTORY = ROOT / "data/atlas/qualification/roi_v3_sources/jrc_flower_detection_source_inventory_v1.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("development", "locked_test"), required=True)
    parser.add_argument("--jrc-root", type=Path, required=True)
    parser.add_argument("--trained-weight", type=Path, required=True)
    parser.add_argument("--training-result-manifest", type=Path, required=True)
    parser.add_argument("--efficient-sam-weights-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--development-result", type=Path)
    parser.add_argument("--torch-threads", type=int, default=8)
    return parser.parse_args()


def clipped_references(
    annotation: dict[str, Any], image: dict[str, Any]
) -> tuple[list[dict[str, Any]], int]:
    references = []
    source_not_evaluable = 0
    width = int(image["width"])
    height = int(image["height"])
    for row in annotation["annotations"]:
        if int(row["image_id"]) != int(image["id"]):
            continue
        x, y, box_width, box_height = (float(value) for value in row["bbox"])
        x0 = max(0.0, min(float(width), x))
        y0 = max(0.0, min(float(height), y))
        x1 = max(0.0, min(float(width), x + box_width))
        y1 = max(0.0, min(float(height), y + box_height))
        if x1 <= x0 or y1 <= y0:
            source_not_evaluable += 1
            continue
        area_512 = ((x1 - x0) * 512 / width) * ((y1 - y0) * 512 / height)
        size = "small" if area_512 < 1024 else "medium" if area_512 < 9216 else "large"
        references.append(
            {
                "annotation_id": int(row["id"]),
                "box_xyxy": [x0, y0, x1, y1],
                "size_bin": size,
            }
        )
    return references, source_not_evaluable


def make_reference_union(
    references: list[dict[str, Any]], *, width: int, height: int
) -> np.ndarray:
    union = np.zeros((height, width), dtype=bool)
    for reference in references:
        x0, y0, x1, y1 = reference["box_xyxy"]
        ix0 = max(0, min(width, math.floor(x0)))
        iy0 = max(0, min(height, math.floor(y0)))
        ix1 = max(0, min(width, math.ceil(x1)))
        iy1 = max(0, min(height, math.ceil(y1)))
        union[iy0:iy1, ix0:ix1] = True
    return union


def letterboxed_rgb(image: Image.Image) -> tuple[np.ndarray, dict[str, int | float]]:
    width, height = image.size
    geometry = letterbox_geometry(width, height)
    resized = image.resize(
        (int(geometry["resized_width"]), int(geometry["resized_height"])),
        Image.Resampling.BILINEAR,
    )
    canvas = Image.new("RGB", (CANVAS_SIZE, CANVAS_SIZE), (114, 114, 114))
    canvas.paste(resized, (int(geometry["pad_left"]), int(geometry["pad_top"])))
    return np.asarray(canvas, dtype=np.uint8), geometry


def canvas_mask_to_original(
    mask: np.ndarray,
    geometry: dict[str, int | float],
    *,
    width: int,
    height: int,
) -> np.ndarray:
    left = int(geometry["pad_left"])
    top = int(geometry["pad_top"])
    right = left + int(geometry["resized_width"])
    bottom = top + int(geometry["resized_height"])
    cropped = Image.fromarray(np.asarray(mask, dtype=np.uint8)[top:bottom, left:right] * 255)
    restored = cropped.resize((width, height), Image.Resampling.NEAREST)
    return np.asarray(restored, dtype=np.uint8) > 0


def background_annulus(
    boxes: list[list[float]], flower_mask: np.ndarray, *, width: int, height: int
) -> np.ndarray:
    expanded = np.zeros((height, width), dtype=bool)
    core = np.zeros((height, width), dtype=bool)
    for raw in boxes:
        x0, y0, x1, y1 = (float(value) for value in raw)
        box_width = x1 - x0
        box_height = y1 - y0
        ex0 = max(0, min(width, math.floor(x0 - 0.25 * box_width)))
        ey0 = max(0, min(height, math.floor(y0 - 0.25 * box_height)))
        ex1 = max(0, min(width, math.ceil(x1 + 0.25 * box_width)))
        ey1 = max(0, min(height, math.ceil(y1 + 0.25 * box_height)))
        ix0 = max(0, min(width, math.floor(x0)))
        iy0 = max(0, min(height, math.floor(y0)))
        ix1 = max(0, min(width, math.ceil(x1)))
        iy1 = max(0, min(height, math.ceil(y1)))
        expanded[ey0:ey1, ex0:ex1] = True
        core[iy0:iy1, ix0:ix1] = True
    return expanded & ~core & ~flower_mask


def masked_mean_lab(rgb: np.ndarray, mask: np.ndarray) -> np.ndarray | None:
    pixels = rgb[mask]
    if not len(pixels):
        return None
    lab = rgb2lab(pixels.reshape(-1, 1, 3).astype(np.float32) / 255.0)
    return np.mean(lab[:, 0, :], axis=0)


class CompositeEstimator:
    def __init__(
        self,
        detector_weight: Path,
        segmenter_dir: Path,
        contract: dict[str, Any],
        *,
        torch_threads: int,
    ) -> None:
        import onnxruntime as ort
        import torch
        from ultralytics import YOLO

        torch.set_num_threads(torch_threads)
        self.detector = YOLO(str(detector_weight))
        options = ort.SessionOptions()
        options.intra_op_num_threads = torch_threads
        options.inter_op_num_threads = 1
        options.execution_mode = ort.ExecutionMode.ORT_SEQUENTIAL
        self.encoder = ort.InferenceSession(
            str(segmenter_dir / "efficient_sam_vitt_encoder.onnx"),
            sess_options=options,
            providers=["CPUExecutionProvider"],
        )
        self.decoder = ort.InferenceSession(
            str(segmenter_dir / "efficient_sam_vitt_decoder.onnx"),
            sess_options=options,
            providers=["CPUExecutionProvider"],
        )
        self.contract = contract

    def analyze(self, image: Image.Image) -> dict[str, Any]:
        inference = self.contract["detector"]["inference"]
        prediction = self.detector.predict(
            source=image,
            imgsz=int(inference["image_size"]),
            conf=float(inference["confidence_minimum"]),
            iou=float(inference["nms_iou"]),
            max_det=int(inference["maximum_detections"]),
            augment=False,
            device="cpu",
            verbose=False,
        )[0]
        if prediction.boxes is None:
            boxes_xyxy = np.empty((0, 4), dtype=float)
            confidences = np.empty(0, dtype=float)
        else:
            boxes_xyxy = prediction.boxes.xyxy.detach().cpu().numpy().astype(float)
            confidences = prediction.boxes.conf.detach().cpu().numpy().astype(float)
        width, height = image.size
        canvas, geometry = letterboxed_rgb(image)
        image_embeddings = self.encoder.run(
            None,
            {"batched_images": canvas.transpose(2, 0, 1)[None].astype(np.float32) / 255.0},
        )[0]
        union_canvas = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=bool)
        retained = 0
        for box in boxes_xyxy:
            canvas_box = box_to_canvas(box, width=width, height=height)
            points = np.array(
                [[[[canvas_box[0], canvas_box[1]], [canvas_box[2], canvas_box[3]]]]],
                dtype=np.float32,
            )
            labels = np.array([[[2.0, 3.0]]], dtype=np.float32)
            output_masks, predicted_iou, _ = self.decoder.run(
                None,
                {
                    "image_embeddings": image_embeddings,
                    "batched_point_coords": points,
                    "batched_point_labels": labels,
                    "orig_im_size": np.array([CANVAS_SIZE, CANVAS_SIZE], dtype=np.int64),
                },
            )
            selected = select_prompt_mask(
                output_masks[0, 0], predicted_iou[0, 0], canvas_box
            )
            if selected.any():
                retained += 1
                union_canvas |= selected
        flower_mask = canvas_mask_to_original(
            union_canvas, geometry, width=width, height=height
        )
        boxes = [list(map(float, box)) for box in boxes_xyxy]
        return {
            "flower_mask": flower_mask,
            "background_mask": background_annulus(
                boxes, flower_mask, width=width, height=height
            ),
            "boxes": boxes,
            "predictions": [
                {
                    "prediction_id": index,
                    "confidence": float(confidences[index]),
                    "box_xyxy": boxes[index],
                }
                for index in range(len(boxes))
            ],
            "retained_instances": retained,
        }


def score_image(
    image_path: Path,
    image_metadata: dict[str, Any],
    annotation: dict[str, Any],
    estimator: CompositeEstimator,
    contract: dict[str, Any],
) -> dict[str, Any]:
    image = ImageOps.exif_transpose(Image.open(image_path)).convert("RGB")
    width, height = image.size
    if (width, height) != (int(image_metadata["width"]), int(image_metadata["height"])):
        raise RuntimeError("JRC image dimensions changed after orientation")
    original = estimator.analyze(image)
    flipped = estimator.analyze(image.transpose(Image.Transpose.FLIP_LEFT_RIGHT))
    flip_flower = np.fliplr(flipped["flower_mask"])
    flip_background = np.fliplr(flipped["background_mask"])
    flower = original["flower_mask"]
    background = original["background_mask"]
    intersection = int(np.count_nonzero(flower & flip_flower))
    union = int(np.count_nonzero(flower | flip_flower))
    flip_iou = intersection / union if union else 0.0
    rgb = np.asarray(image, dtype=np.uint8)
    original_lab = masked_mean_lab(rgb, flower)
    flipped_lab = masked_mean_lab(rgb, flip_flower)
    delta_e = (
        float(np.linalg.norm(original_lab - flipped_lab))
        if original_lab is not None and flipped_lab is not None
        else math.inf
    )
    references, source_not_evaluable = clipped_references(annotation, image_metadata)
    matching = greedy_detection_matches(original["predictions"], references)
    hits = {"small": 0, "medium": 0, "large": 0}
    totals = {"small": 0, "medium": 0, "large": 0}
    for reference in references:
        totals[reference["size_bin"]] += 1
    for match in matching["matches"]:
        hits[references[int(match["reference_index"])]["size_bin"]] += 1
    reference_union = make_reference_union(references, width=width, height=height)
    mask_pixels = int(flower.sum())
    inside = int(np.count_nonzero(flower & reference_union))
    measurement = contract["image_measurement"]
    failures = []
    if original["retained_instances"] < 1:
        failures.append("no_retained_flower_instance")
    if mask_pixels < int(measurement["minimum_union_flower_pixels_on_original"]):
        failures.append("insufficient_flower_pixels")
    if int(background.sum()) < int(measurement["minimum_background_pixels_on_original"]):
        failures.append("insufficient_background_pixels")
    if flip_iou < float(measurement["horizontal_flip_mask_iou_minimum"]):
        failures.append("horizontal_flip_mask_instability")
    if delta_e > float(measurement["horizontal_flip_colour_delta_e_maximum"]):
        failures.append("horizontal_flip_colour_instability")
    return {
        "image_id": int(image_metadata["id"]),
        "file_name": str(image_metadata["file_name"]),
        "image_sha256": sha256(image_path),
        "detector_predictions": len(original["predictions"]),
        "retained_instances": int(original["retained_instances"]),
        "flip_detector_predictions": len(flipped["predictions"]),
        "flip_retained_instances": int(flipped["retained_instances"]),
        "true_positive": int(matching["true_positive"]),
        "false_positive": int(matching["false_positive"]),
        "false_negative": int(matching["false_negative"]),
        "mask_pixels": mask_pixels,
        "mask_pixels_inside_reference_box_union": inside,
        "image_mask_pixels_inside_reference_box_union": inside / mask_pixels if mask_pixels else 0.0,
        "background_pixels": int(background.sum()),
        "flip_background_pixels": int(flip_background.sum()),
        "horizontal_flip_mask_iou": float(flip_iou),
        "horizontal_flip_colour_delta_e": None if not math.isfinite(delta_e) else delta_e,
        "source_annotation_boxes": len(references) + source_not_evaluable,
        "reference_boxes": len(references),
        "source_not_evaluable_boxes": source_not_evaluable,
        **{f"{size}_reference_boxes": totals[size] for size in totals},
        **{f"{size}_hit_boxes": hits[size] for size in hits},
        "background_features_available": int(background.sum())
        >= int(measurement["minimum_background_pixels_on_original"]),
        "estimator_admitted": not failures,
        "failure_reasons": ";".join(failures),
    }


def main() -> None:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    amendment = json.loads(AMENDMENT.read_text(encoding="utf-8"))
    validate_roi_v4_contract(contract)
    validate_reference_size_amendment(amendment)
    if amendment["parent_contract"]["sha256_lf_canonical_v1"] != canonical_sha256(CONTRACT):
        raise RuntimeError("ROI v4 amendment parent changed")
    training = json.loads(args.training_result_manifest.read_text(encoding="utf-8"))
    trained_weight_sha = sha256(args.trained_weight)
    if (
        training.get("status") != "complete_roi_v4_detector_training_not_yet_qualified"
        or training.get("epochs") != 50
        or training.get("weight_selection") != "last epoch only; never best epoch"
        or training.get("trained_weight_sha256") != trained_weight_sha
        or training.get("jrc_test_images_decoded_or_scored") is not False
    ):
        raise RuntimeError("ROI v4 trained-weight evidence changed")
    if args.phase == "locked_test":
        if args.development_result is None:
            raise RuntimeError("locked test requires committed development evidence")
        development = json.loads(args.development_result.read_text(encoding="utf-8"))
        if (
            development.get("status") != "pass_roi_v4_development"
            or development.get("trained_weight_sha256") != trained_weight_sha
            or development.get("jrc_locked_test_permitted") is not True
            or development.get("jrc_test_images_decoded_or_scored") is not False
        ):
            raise RuntimeError("ROI v4 development did not authorize the locked test")
    elif args.development_result is not None:
        raise RuntimeError("development must not read a prior development outcome")

    segmenter = contract["mask_generator"]
    encoder = args.efficient_sam_weights_dir / Path(segmenter["encoder_path"]).name
    decoder = args.efficient_sam_weights_dir / Path(segmenter["decoder_path"]).name
    if sha256(encoder) != segmenter["encoder_sha256"] or sha256(decoder) != segmenter["decoder_sha256"]:
        raise RuntimeError("EfficientSAM v4 weight identity changed")
    split = "train" if args.phase == "development" else "test"
    annotation_path = args.jrc_root / f"annotations/{split}.json"
    expected_annotation = contract["jrc_source"][f"{split}_annotation_sha256"]
    if sha256(annotation_path) != expected_annotation:
        raise RuntimeError(f"JRC {split} annotation identity changed")
    annotation = json.loads(annotation_path.read_text(encoding="utf-8"))
    images = sorted(annotation["images"], key=lambda row: int(row["id"]))
    expected_images = 400 if split == "train" else 100
    if len(images) != expected_images:
        raise RuntimeError("JRC gate image denominator changed")
    inventory = {
        (row["split"], row["file_name"]): row for row in read_csv(SOURCE_INVENTORY)
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output_dir / "per_image_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    estimator = CompositeEstimator(
        args.trained_weight,
        args.efficient_sam_weights_dir,
        contract,
        torch_threads=args.torch_threads,
    )
    identity = {
        "protocol": contract["protocol"],
        "phase": args.phase,
        "trained_weight_sha256": trained_weight_sha,
        "encoder_sha256": segmenter["encoder_sha256"],
        "decoder_sha256": segmenter["decoder_sha256"],
    }
    rows = []
    started = time.time()
    for index, image_metadata in enumerate(images, start=1):
        name = str(image_metadata["file_name"])
        image_path = args.jrc_root / f"images/{split}" / name
        source = inventory.get((split, name))
        if source is None or sha256(image_path) != source["image_sha256"]:
            raise RuntimeError(f"JRC {split} image identity changed: {name}")
        cache_path = cache_dir / f"{int(image_metadata['id']):04d}.json"
        if cache_path.is_file():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("identity") != identity:
                raise RuntimeError("incompatible ROI v4 per-image cache")
            row = cached["row"]
        else:
            row = score_image(image_path, image_metadata, annotation, estimator, contract)
            cache_path.write_text(
                json.dumps({"identity": identity, "row": row}, indent=2, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
        rows.append(row)
        print(f"{args.phase}={index}/{len(images)}", flush=True)
    rows_path = args.output_dir / f"jrc_roi_v4_{args.phase}_rows.csv"
    with rows_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    decision = summarize_composite_gate(rows, contract, phase=args.phase)
    decision.update(
        {
            "trained_weight_sha256": trained_weight_sha,
            "contract_sha256_lf_canonical_v1": canonical_sha256(CONTRACT),
            "reference_size_amendment_sha256_lf_canonical_v1": canonical_sha256(AMENDMENT),
            "source_annotation_sha256": expected_annotation,
            "rows_sha256": sha256(rows_path),
            "source_annotation_boxes": sum(int(row["source_annotation_boxes"]) for row in rows),
            "reference_boxes": sum(int(row["reference_boxes"]) for row in rows),
            "source_not_evaluable_boxes": sum(int(row["source_not_evaluable_boxes"]) for row in rows),
            "elapsed_seconds_new_inference": time.time() - started,
            "jrc_test_images_decoded_or_scored": args.phase == "locked_test",
            "scaleout_candidate_pixels_opened": False,
        }
    )
    result_path = args.output_dir / f"jrc_roi_v4_{args.phase}_result.json"
    result_path.write_text(
        json.dumps(decision, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
