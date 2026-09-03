"""Shared frozen ROI-v4 image measurement used by qualification and the atlas.

The public interface is one estimator with one ``measure`` operation.  Model
loading, letterboxing, box prompting, hard-mask colour summaries and reflection
checks stay behind that seam so JRC qualification and scale-out cannot drift.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from PIL import Image, ImageOps
from skimage.color import rgb2lab

from .flower_roi_v4 import (
    CANVAS_SIZE,
    PROTOCOL,
    box_to_canvas,
    letterbox_geometry,
    select_prompt_mask,
    validate_roi_v4_contract,
)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_scaleout_authorization(
    locked_result: Mapping[str, Any], *, trained_weight_sha256: str | None = None
) -> str:
    """Require the one permitted locked-JRC outcome before atlas image access."""

    observed_weight = str(locked_result.get("trained_weight_sha256") or "")
    if (
        locked_result.get("protocol") != PROTOCOL
        or locked_result.get("phase") != "locked_test"
        or locked_result.get("status") != "pass_roi_v4_locked_test"
        or locked_result.get("jrc_test_images_decoded_or_scored") is not True
        or locked_result.get("scaleout_candidate_pixels_permitted") is not True
        or locked_result.get("scaleout_candidate_pixels_opened") is not False
        or not observed_weight
    ):
        raise RuntimeError("ROI v4 locked test has not authorized scale-out pixels")
    if trained_weight_sha256 is not None and observed_weight != trained_weight_sha256:
        raise RuntimeError("ROI v4 locked result and trained weight do not match")
    return observed_weight


def _letterboxed_rgb(image: Image.Image) -> tuple[np.ndarray, dict[str, int | float]]:
    width, height = image.size
    geometry = letterbox_geometry(width, height)
    resized = image.resize(
        (int(geometry["resized_width"]), int(geometry["resized_height"])),
        Image.Resampling.BILINEAR,
    )
    canvas = Image.new("RGB", (CANVAS_SIZE, CANVAS_SIZE), (114, 114, 114))
    canvas.paste(resized, (int(geometry["pad_left"]), int(geometry["pad_top"])))
    return np.asarray(canvas, dtype=np.uint8), geometry


def _canvas_mask_to_original(
    mask: np.ndarray,
    geometry: Mapping[str, int | float],
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


def _background_annulus(
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


def hard_mask_lab_summary(
    rgb: np.ndarray, mask: np.ndarray, prefix: str
) -> dict[str, float | int | None]:
    """Summarize a hard original-resolution mask without calibrated-colour claims."""

    if rgb.ndim != 3 or rgb.shape[2] != 3 or mask.shape != rgb.shape[:2]:
        raise ValueError("RGB image and hard mask shapes are incompatible")
    pixels = int(np.count_nonzero(mask))
    result: dict[str, float | int | None] = {f"{prefix}_effective_pixels": pixels}
    if pixels == 0:
        for channel in ("L", "a", "b"):
            for metric in ("mean", "sd", "q10", "q50", "q90"):
                result[f"{prefix}_{channel}_{metric}"] = None
        return result
    lab = rgb2lab(rgb.astype(np.float32) / 255.0)
    selected = lab[mask].astype(float)
    for index, channel in enumerate(("L", "a", "b")):
        values = selected[:, index]
        result[f"{prefix}_{channel}_mean"] = float(np.mean(values))
        result[f"{prefix}_{channel}_sd"] = float(np.std(values, ddof=0))
        for metric, probability in (("q10", 0.1), ("q50", 0.5), ("q90", 0.9)):
            result[f"{prefix}_{channel}_{metric}"] = float(
                np.quantile(values, probability, method="linear")
            )
    return result


def _mean_lab(summary: Mapping[str, Any], prefix: str) -> np.ndarray | None:
    values = [summary.get(f"{prefix}_{channel}_mean") for channel in ("L", "a", "b")]
    if any(value is None or not math.isfinite(float(value)) for value in values):
        return None
    return np.asarray(values, dtype=float)


def summarize_hard_mask_measurement(
    rgb: np.ndarray,
    flower_mask: np.ndarray,
    background_mask: np.ndarray,
    flip_flower_mask: np.ndarray,
    flip_background_mask: np.ndarray,
    contract: Mapping[str, Any],
    *,
    retained_instances: int,
) -> dict[str, Any]:
    """Apply the frozen v4 colour and admission rules to two orientation runs."""

    validate_roi_v4_contract(contract)
    masks = (flower_mask, background_mask, flip_flower_mask, flip_background_mask)
    if any(mask.shape != rgb.shape[:2] for mask in masks):
        raise ValueError("ROI v4 masks must use original image dimensions")
    flower = hard_mask_lab_summary(rgb, flower_mask, "flower")
    background = hard_mask_lab_summary(rgb, background_mask, "background")
    flip = hard_mask_lab_summary(rgb, flip_flower_mask, "flip")
    intersection = int(np.count_nonzero(flower_mask & flip_flower_mask))
    union = int(np.count_nonzero(flower_mask | flip_flower_mask))
    flip_iou = intersection / union if union else 0.0
    flower_mean = _mean_lab(flower, "flower")
    flip_mean = _mean_lab(flip, "flip")
    flip_delta = (
        float(np.linalg.norm(flower_mean - flip_mean))
        if flower_mean is not None and flip_mean is not None
        else None
    )
    measurement = contract["image_measurement"]
    failures: list[str] = []
    if retained_instances < 1:
        failures.append("no_retained_flower_instance")
    if int(flower["flower_effective_pixels"]) < int(
        measurement["minimum_union_flower_pixels_on_original"]
    ):
        failures.append("insufficient_flower_pixels")
    if int(background["background_effective_pixels"]) < int(
        measurement["minimum_background_pixels_on_original"]
    ):
        failures.append("insufficient_background_pixels")
    if flip_iou < float(measurement["horizontal_flip_mask_iou_minimum"]):
        failures.append("horizontal_flip_mask_instability")
    if flip_delta is None or flip_delta > float(
        measurement["horizontal_flip_colour_delta_e_maximum"]
    ):
        failures.append("horizontal_flip_colour_instability")
    return {
        **flower,
        **background,
        "flip_effective_pixels": int(flip["flip_effective_pixels"]),
        "flip_background_pixels": int(np.count_nonzero(flip_background_mask)),
        "horizontal_flip_mask_iou": float(flip_iou),
        "horizontal_flip_colour_delta_e": flip_delta,
        "background_features_available": int(
            background["background_effective_pixels"]
        )
        >= int(measurement["minimum_background_pixels_on_original"]),
        "automated_colour_state_status": (
            "automated_colour_state_admitted"
            if not failures
            else "automated_colour_state_not_evaluable"
        ),
        "failure_reasons": ";".join(failures),
    }


class FrozenFlowerColourEstimator:
    """Load and apply the exact detector-plus-segmenter frozen by ROI v4."""

    def __init__(
        self,
        detector_weight: Path,
        segmenter_dir: Path,
        contract: Mapping[str, Any],
        *,
        torch_threads: int,
    ) -> None:
        import onnxruntime as ort
        import torch
        from ultralytics import YOLO

        validate_roi_v4_contract(contract)
        segmenter = contract["mask_generator"]
        encoder = segmenter_dir / Path(segmenter["encoder_path"]).name
        decoder = segmenter_dir / Path(segmenter["decoder_path"]).name
        if (
            file_sha256(encoder) != segmenter["encoder_sha256"]
            or file_sha256(decoder) != segmenter["decoder_sha256"]
        ):
            raise RuntimeError("EfficientSAM v4 weight identity changed")
        torch.set_num_threads(torch_threads)
        self.detector = YOLO(str(detector_weight))
        options = ort.SessionOptions()
        options.intra_op_num_threads = torch_threads
        options.inter_op_num_threads = 1
        options.execution_mode = ort.ExecutionMode.ORT_SEQUENTIAL
        self.encoder = ort.InferenceSession(
            str(encoder), sess_options=options, providers=["CPUExecutionProvider"]
        )
        self.decoder = ort.InferenceSession(
            str(decoder), sess_options=options, providers=["CPUExecutionProvider"]
        )
        self.contract = dict(contract)

    def _analyze_orientation(self, image: Image.Image) -> dict[str, Any]:
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
        union_canvas = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=bool)
        retained = 0
        if len(boxes_xyxy):
            canvas, geometry = _letterboxed_rgb(image)
            embeddings = self.encoder.run(
                None,
                {"batched_images": canvas.transpose(2, 0, 1)[None].astype(np.float32) / 255.0},
            )[0]
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
                        "image_embeddings": embeddings,
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
        else:
            geometry = letterbox_geometry(width, height)
        flower_mask = _canvas_mask_to_original(
            union_canvas, geometry, width=width, height=height
        )
        boxes = [list(map(float, box)) for box in boxes_xyxy]
        return {
            "flower_mask": flower_mask,
            "background_mask": _background_annulus(
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

    def measure(self, image: Image.Image) -> dict[str, Any]:
        """Return one complete, location-free v4 measurement record."""

        oriented = ImageOps.exif_transpose(image).convert("RGB")
        original = self._analyze_orientation(oriented)
        flipped = self._analyze_orientation(
            oriented.transpose(Image.Transpose.FLIP_LEFT_RIGHT)
        )
        flip_flower = np.fliplr(flipped["flower_mask"])
        flip_background = np.fliplr(flipped["background_mask"])
        rgb = np.asarray(oriented, dtype=np.uint8)
        summary = summarize_hard_mask_measurement(
            rgb,
            original["flower_mask"],
            original["background_mask"],
            flip_flower,
            flip_background,
            self.contract,
            retained_instances=int(original["retained_instances"]),
        )
        return {
            **summary,
            "image_size": oriented.size,
            "flower_mask": original["flower_mask"],
            "background_mask": original["background_mask"],
            "flip_flower_mask": flip_flower,
            "flip_background_mask": flip_background,
            "boxes": original["boxes"],
            "predictions": original["predictions"],
            "retained_instances": int(original["retained_instances"]),
            "flip_detector_predictions": len(flipped["predictions"]),
            "flip_retained_instances": int(flipped["retained_instances"]),
        }
