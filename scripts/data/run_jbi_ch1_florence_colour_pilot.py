#!/usr/bin/env python3
"""Quota-independent six-image flower ROI pilot using Florence-2 + numeric colour features.

This pilot is calibration-only and diagnostic. Florence-2 is used only to localize the
visible flower region. Biological colour-state suggestions are then derived from pixels
inside the localized ROI using a fixed reference palette declared in this file. No
candidate state produced here is a final label.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import torch
from transformers import AutoModelForMultimodalLM, AutoProcessor

MODEL_ID = "florence-community/Florence-2-base-ft"
PROTOCOL = "jbi-ch1-florence-colour-pilot-v1"

# Fixed before running this pilot. These are generic sRGB visual anchors, not fitted to
# the six pilot images. Green/brown/black are nuisance/background anchors.
REFERENCE_RGB = {
    "white": (245, 245, 245),
    "yellow": (245, 210, 45),
    "orange": (235, 130, 35),
    "red": (195, 40, 45),
    "pink": (235, 115, 165),
    "magenta": (180, 45, 145),
    "purple": (115, 65, 155),
    "blue": (65, 95, 185),
    "bronze": (150, 100, 55),
    "green": (75, 130, 65),
    "brown": (100, 70, 45),
    "black": (25, 25, 25),
}
NUISANCE = {"green", "brown", "black"}


def srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    x = np.asarray(rgb, dtype=np.float64) / 255.0
    x = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
    m = np.array([
        [0.4124564, 0.3575761, 0.1804375],
        [0.2126729, 0.7151522, 0.0721750],
        [0.0193339, 0.1191920, 0.9503041],
    ])
    xyz = x @ m.T
    xyz = xyz / np.array([0.95047, 1.0, 1.08883])
    eps = 216 / 24389
    kappa = 24389 / 27
    f = np.where(xyz > eps, np.cbrt(xyz), (kappa * xyz + 16) / 116)
    L = 116 * f[..., 1] - 16
    a = 500 * (f[..., 0] - f[..., 1])
    b = 200 * (f[..., 1] - f[..., 2])
    return np.stack([L, a, b], axis=-1)


REFERENCE_LAB = {
    name: srgb_to_lab(np.array(rgb, dtype=np.uint8).reshape(1, 3))[0]
    for name, rgb in REFERENCE_RGB.items()
}


def nearest_palette_counts(image: Image.Image, bbox: list[float]) -> dict[str, int]:
    w, h = image.size
    x0, y0, x1, y1 = bbox
    x0 = max(0, min(w - 1, int(math.floor(x0))))
    y0 = max(0, min(h - 1, int(math.floor(y0))))
    x1 = max(x0 + 1, min(w, int(math.ceil(x1))))
    y1 = max(y0 + 1, min(h, int(math.ceil(y1))))
    crop = image.crop((x0, y0, x1, y1)).convert("RGB")
    crop.thumbnail((256, 256), Image.Resampling.LANCZOS)
    arr = np.asarray(crop, dtype=np.uint8).reshape(-1, 3)
    lab = srgb_to_lab(arr)
    names = list(REFERENCE_LAB)
    refs = np.stack([REFERENCE_LAB[n] for n in names])
    distances = ((lab[:, None, :] - refs[None, :, :]) ** 2).sum(axis=2)
    nearest = distances.argmin(axis=1)
    counts = Counter(names[int(i)] for i in nearest)
    return {name: int(counts.get(name, 0)) for name in names}


def flower_only_fractions(counts: dict[str, int]) -> dict[str, float]:
    usable = {k: v for k, v in counts.items() if k not in NUISANCE}
    total = sum(usable.values())
    if total <= 0:
        return {k: 0.0 for k in usable}
    return {k: v / total for k, v in usable.items()}


def deterministic_candidate(species: str, fractions: dict[str, float]) -> tuple[str, dict[str, float]]:
    f = lambda *names: sum(fractions.get(name, 0.0) for name in names)
    if species == "Ipomoea purpurea":
        scores = {"white": f("white"), "pink": f("pink", "magenta"), "blue_purple": f("blue", "purple")}
    elif species == "Raphanus sativus":
        scores = {"white": f("white"), "yellow": f("yellow"), "pink": f("pink", "magenta", "purple"), "bronze": f("bronze", "orange")}
    elif species == "Gentiana lutea":
        scores = {"yellow": f("yellow"), "orange": f("orange", "bronze")}
    elif species == "Dactylorhiza sambucina":
        scores = {"yellow": f("yellow", "white"), "purple": f("purple", "magenta", "pink")}
    elif species == "Lysimachia arvensis":
        scores = {"blue": f("blue", "purple"), "red": f("red", "orange", "magenta", "pink")}
    elif species == "Antirrhinum majus":
        magenta = f("magenta", "pink", "purple", "red")
        yellow = f("yellow", "orange")
        scores = {
            "magenta_pseudomajus_like": magenta,
            "yellow_striatum_like": yellow,
            "intermediate_or_other": min(magenta, yellow),
        }
    else:
        raise ValueError(species)
    best = max(scores, key=scores.get)
    return best, scores


def extract_boxes(parsed: Any, task: str) -> list[list[float]]:
    if not isinstance(parsed, dict):
        return []
    payload = parsed.get(task, parsed)
    if not isinstance(payload, dict):
        return []
    boxes = payload.get("bboxes", [])
    out = []
    for box in boxes:
        if isinstance(box, (list, tuple)) and len(box) == 4:
            vals = [float(v) for v in box]
            if vals[2] > vals[0] and vals[3] > vals[1]:
                out.append(vals)
    return out


def choose_box(boxes: list[list[float]], image_size: tuple[int, int]) -> list[float] | None:
    if not boxes:
        return None
    w, h = image_size
    cx, cy = w / 2, h / 2
    def key(box):
        x0, y0, x1, y1 = box
        area = (x1 - x0) * (y1 - y0)
        bx, by = (x0 + x1) / 2, (y0 + y1) / 2
        center_penalty = ((bx - cx) / max(w, 1)) ** 2 + ((by - cy) / max(h, 1)) ** 2
        # Prefer substantial, centrally located localized flowers while avoiding a
        # full-frame box if a more specific flower is available.
        return area / (1.0 + 3.0 * center_penalty)
    candidates = [b for b in boxes if ((b[2]-b[0])*(b[3]-b[1])) < 0.95*w*h] or boxes
    return max(candidates, key=key)


def run_task(model, processor, image: Image.Image, task: str, text_input: str = "") -> Any:
    prompt = task + text_input
    inputs = processor(text=prompt, images=image, return_tensors="pt")
    with torch.inference_mode():
        generated = model.generate(
            input_ids=inputs["input_ids"],
            pixel_values=inputs["pixel_values"],
            max_new_tokens=512,
            num_beams=3,
            do_sample=False,
        )
    text = processor.batch_decode(generated, skip_special_tokens=False)[0]
    return processor.post_process_generation(text, task=task, image_size=image.size)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-manifest", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--overlay-dir", type=Path, required=True)
    args = parser.parse_args()

    frame = pd.read_csv(args.pilot_manifest)
    if len(frame) != 6 or frame["species"].nunique() != 6:
        raise RuntimeError("Florence pilot requires six calibration-only images")
    if frame.get("evaluation_row", pd.Series([False] * len(frame))).astype(bool).any():
        raise RuntimeError("evaluation leakage")

    processor = AutoProcessor.from_pretrained(MODEL_ID)
    model = AutoModelForMultimodalLM.from_pretrained(MODEL_ID)
    model.eval()
    args.overlay_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for _, row in frame.sort_values("species", kind="mergesort").iterrows():
        image = Image.open(str(row["image_path"])).convert("RGB")
        detection = run_task(model, processor, image, "<OPEN_VOCABULARY_DETECTION>", "flower")
        boxes = extract_boxes(detection, "<OPEN_VOCABULARY_DETECTION>")
        if not boxes:
            detection = run_task(model, processor, image, "<OPEN_VOCABULARY_DETECTION>", "flower petals")
            boxes = extract_boxes(detection, "<OPEN_VOCABULARY_DETECTION>")
        box = choose_box(boxes, image.size)
        if box is None:
            # Diagnostic fallback only. The pilot records this explicitly; it is not a
            # valid segmentation result for later inference.
            w, h = image.size
            box = [0.15*w, 0.15*h, 0.85*w, 0.85*h]
            localization_status = "fallback_central_crop"
        else:
            localization_status = "florence_open_vocab_box"

        counts = nearest_palette_counts(image, box)
        fractions = flower_only_fractions(counts)
        candidate, scores = deterministic_candidate(str(row["species"]), fractions)

        overlay = image.copy()
        draw = ImageDraw.Draw(overlay)
        draw.rectangle(box, outline="red", width=max(2, min(image.size)//250))
        overlay_path = args.overlay_dir / f"{row['blind_id']}.jpg"
        overlay.save(overlay_path, quality=90)

        results.append({
            "protocol": PROTOCOL,
            "model": MODEL_ID,
            "species": str(row["species"]),
            "blind_id": str(row["blind_id"]),
            "evaluation_row": False,
            "pilot_only": True,
            "final_label": False,
            "localization_status": localization_status,
            "n_detected_boxes": len(boxes),
            "selected_bbox": [round(float(x), 2) for x in box],
            "palette_counts": counts,
            "flower_only_fractions": {k: round(float(v), 6) for k, v in fractions.items()},
            "deterministic_candidate_state": candidate,
            "candidate_scores": {k: round(float(v), 6) for k, v in scores.items()},
            "overlay_path": str(overlay_path),
        })
        print(json.dumps(results[-1], ensure_ascii=False), flush=True)

    payload = {
        "protocol": PROTOCOL,
        "status": "pilot_complete_not_final",
        "model": MODEL_ID,
        "n_images": len(results),
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "reference_palette_fitted_to_pilot": False,
        "results": results,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
