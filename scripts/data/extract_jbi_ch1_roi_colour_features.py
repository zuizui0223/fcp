#!/usr/bin/env python3
"""Extract deterministic continuous CIELAB vectors from frozen Florence ROIs.

This is a measurement conversion, not a classifier.  It does not use coordinates,
species outcomes, graph structure, thresholds learned from evaluation data, or any
biological labels.  The Florence-selected flower box is kept fixed.  Within that box we
use a fixed inner ellipse and a 10% component-wise trimmed mean in CIE L*a*b* (D65).
The same conversion is applied to calibration and evaluation photographs.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from PIL import Image, ImageOps

USER_AGENT = "fcp-jbi-ch1-colour-measurement/1.0"
MAX_SIDE = 512
INNER_MARGIN_FRACTION = 0.05
ELLIPSE_RADIUS_FRACTION = 0.48
TRIM_FRACTION = 0.10

EXACT_BOX_PATHS = (
    "selected_box_xyxy",
    "selected_bbox_xyxy",
    "flower_box_xyxy",
    "flower_bbox_xyxy",
    "floral_box_xyxy",
    "floral_bbox_xyxy",
    "roi_box_xyxy",
    "roi_bbox_xyxy",
    "bbox_xyxy",
    "box_xyxy",
    "features.selected_box_xyxy",
    "features.flower_box_xyxy",
    "features.roi_box_xyxy",
    "florence.selected_box_xyxy",
    "florence.flower_box_xyxy",
    "florence.roi_box_xyxy",
    "florence_result.selected_box_xyxy",
    "florence_result.flower_box_xyxy",
    "measurement.selected_box_xyxy",
    "measurement.flower_box_xyxy",
    "selected_box",
    "selected_bbox",
    "flower_box",
    "flower_bbox",
    "roi_box",
    "roi_bbox",
)
LIST_BOX_PATHS = (
    "bboxes",
    "boxes",
    "features.bboxes",
    "features.boxes",
    "florence.bboxes",
    "florence.boxes",
    "florence_result.bboxes",
    "florence_result.boxes",
)


def get_path(row: dict[str, Any], path: str) -> Any:
    value: Any = row
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(path)
        value = value[part]
    return value


def iter_leaf_paths(value: Any, prefix: str = "") -> Iterable[tuple[str, Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            yield from iter_leaf_paths(child, path)
    else:
        yield prefix, value


def _numeric_box(value: Any) -> np.ndarray:
    if isinstance(value, dict):
        key_orders = (
            ("x1", "y1", "x2", "y2"),
            ("xmin", "ymin", "xmax", "ymax"),
            ("left", "top", "right", "bottom"),
        )
        for keys in key_orders:
            if all(key in value for key in keys):
                value = [value[key] for key in keys]
                break
        else:
            raise ValueError("unrecognized box object")
    arr = np.asarray(value, dtype=float)
    if arr.shape != (4,) or not np.isfinite(arr).all():
        raise ValueError("box is not a finite four-vector")
    return arr


def _largest_box(values: Any) -> np.ndarray:
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError("box list is empty")
    candidates = []
    for value in values:
        try:
            box = _numeric_box(value)
        except (TypeError, ValueError):
            continue
        area = max(0.0, float(box[2] - box[0])) * max(0.0, float(box[3] - box[1]))
        candidates.append((area, tuple(float(x) for x in box), box))
    if not candidates:
        raise ValueError("box list contains no finite four-vector")
    # Geometry-only and deterministic.  The tuple breaks equal-area ties reproducibly.
    return max(candidates, key=lambda item: (item[0], item[1]))[2]


def find_box(row: dict[str, Any]) -> tuple[np.ndarray, str]:
    for path in EXACT_BOX_PATHS:
        try:
            value = get_path(row, path)
        except KeyError:
            continue
        try:
            return _numeric_box(value), path
        except (TypeError, ValueError):
            pass
    for path in LIST_BOX_PATHS:
        try:
            value = get_path(row, path)
        except KeyError:
            continue
        try:
            return _largest_box(value), path + "[largest_geometry_only]"
        except (TypeError, ValueError):
            pass

    # Last structural fallback: only semantically explicit box/bbox leaves are eligible.
    candidates: list[tuple[int, str, np.ndarray]] = []
    for path, value in iter_leaf_paths(row):
        lower = path.lower()
        if "box" not in lower and "bbox" not in lower:
            continue
        if any(token in lower for token in ("count", "status", "area", "fraction", "label")):
            continue
        try:
            box = _numeric_box(value)
        except (TypeError, ValueError):
            try:
                box = _largest_box(value)
            except (TypeError, ValueError):
                continue
        priority = 0
        priority += 8 if "selected" in lower or "primary" in lower else 0
        priority += 4 if any(token in lower for token in ("flower", "floral", "corolla", "roi")) else 0
        priority += 2 if "xyxy" in lower else 0
        candidates.append((priority, path, box))
    if not candidates:
        raise ValueError("no frozen Florence flower box was found")
    candidates.sort(key=lambda item: (-item[0], item[1]))
    return candidates[0][2], candidates[0][1]


def normalize_box(box: np.ndarray, width: int, height: int, source_path: str) -> tuple[int, int, int, int]:
    arr = np.asarray(box, dtype=float).copy()
    lower = source_path.lower()
    if "xywh" in lower:
        arr[2] = arr[0] + arr[2]
        arr[3] = arr[1] + arr[3]
    elif not (arr[2] > arr[0] and arr[3] > arr[1]):
        # Only use xywh as a structural rescue when xyxy is impossible.
        arr[2] = arr[0] + arr[2]
        arr[3] = arr[1] + arr[3]
    if np.min(arr) >= 0.0 and np.max(arr) <= 1.000001:
        arr[[0, 2]] *= width
        arr[[1, 3]] *= height
    x1, y1, x2, y2 = arr
    x1 = int(max(0, min(width - 1, math.floor(x1))))
    y1 = int(max(0, min(height - 1, math.floor(y1))))
    x2 = int(max(x1 + 1, min(width, math.ceil(x2))))
    y2 = int(max(y1 + 1, min(height, math.ceil(y2))))
    if x2 - x1 < 4 or y2 - y1 < 4:
        raise ValueError(f"frozen box is too small after clamping: {(x1, y1, x2, y2)}")
    return x1, y1, x2, y2


def download_image(url: str, retries: int = 4) -> tuple[Image.Image, str]:
    if not url:
        raise ValueError("photo_url is empty")
    error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=45) as response:
                payload = response.read()
            if not payload:
                raise ValueError("empty image response")
            image = Image.open(io.BytesIO(payload))
            image = ImageOps.exif_transpose(image).convert("RGB")
            return image, hashlib.sha256(payload).hexdigest()
        except (OSError, ValueError, urllib.error.URLError, TimeoutError) as exc:
            error = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"image download/decode failed after {retries} attempts: {error}")


def srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    x = np.asarray(rgb, dtype=np.float64) / 255.0
    x = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
    matrix = np.array(
        [
            [0.4124564, 0.3575761, 0.1804375],
            [0.2126729, 0.7151522, 0.0721750],
            [0.0193339, 0.1191920, 0.9503041],
        ],
        dtype=np.float64,
    )
    xyz = x @ matrix.T
    white = np.array([0.95047, 1.0, 1.08883], dtype=np.float64)
    ratio = xyz / white
    delta = 6.0 / 29.0
    f = np.where(ratio > delta**3, np.cbrt(ratio), ratio / (3 * delta**2) + 4.0 / 29.0)
    lstar = 116.0 * f[..., 1] - 16.0
    astar = 500.0 * (f[..., 0] - f[..., 1])
    bstar = 200.0 * (f[..., 1] - f[..., 2])
    return np.stack((lstar, astar, bstar), axis=-1)


def trimmed_mean(values: np.ndarray, trim_fraction: float) -> np.ndarray:
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("expected an n x 3 Lab matrix")
    n = len(values)
    trim = int(math.floor(n * trim_fraction))
    if n - 2 * trim < 20:
        trim = max(0, (n - 20) // 2)
    result = []
    for column in range(3):
        ordered = np.sort(values[:, column])
        kept = ordered[trim : n - trim] if trim else ordered
        result.append(float(np.mean(kept)))
    return np.asarray(result, dtype=float)


def measure(image: Image.Image, box: tuple[int, int, int, int]) -> tuple[np.ndarray, dict[str, Any]]:
    crop = image.crop(box)
    width, height = crop.size
    mx = int(round(width * INNER_MARGIN_FRACTION))
    my = int(round(height * INNER_MARGIN_FRACTION))
    if width - 2 * mx >= 4 and height - 2 * my >= 4:
        crop = crop.crop((mx, my, width - mx, height - my))
    width, height = crop.size
    scale = min(1.0, MAX_SIDE / max(width, height))
    if scale < 1.0:
        crop = crop.resize(
            (max(4, int(round(width * scale))), max(4, int(round(height * scale)))),
            Image.Resampling.LANCZOS,
        )
    rgb = np.asarray(crop, dtype=np.uint8)
    height, width = rgb.shape[:2]
    yy, xx = np.ogrid[:height, :width]
    cx = (width - 1) / 2.0
    cy = (height - 1) / 2.0
    rx = max(1.0, width * ELLIPSE_RADIUS_FRACTION)
    ry = max(1.0, height * ELLIPSE_RADIUS_FRACTION)
    mask = ((xx - cx) / rx) ** 2 + ((yy - cy) / ry) ** 2 <= 1.0
    pixels = rgb[mask]
    if len(pixels) < 50:
        raise ValueError(f"fewer than 50 pixels remain in the fixed inner ellipse ({len(pixels)})")
    lab = srgb_to_lab(pixels)
    finite = np.isfinite(lab).all(axis=1)
    lab = lab[finite]
    # Remove only unusable sensor/encoding extremes, symmetrically and with fixed limits.
    lab = lab[(lab[:, 0] >= 1.0) & (lab[:, 0] <= 99.0)]
    if len(lab) < 50:
        raise ValueError(f"fewer than 50 finite non-extreme Lab pixels remain ({len(lab)})")
    vector = trimmed_mean(lab, TRIM_FRACTION)
    if not np.isfinite(vector).all():
        raise ValueError("continuous Lab vector is non-finite")
    diagnostics = {
        "crop_width_px": int(width),
        "crop_height_px": int(height),
        "ellipse_pixels": int(mask.sum()),
        "measured_pixels": int(len(lab)),
        "trim_fraction": TRIM_FRACTION,
    }
    return vector, diagnostics


def read_rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            if not isinstance(row, dict):
                raise ValueError(f"{path}:{line_number}: expected a JSON object")
            rows.append(row)
    return rows


def stable_shard(rows: list[dict[str, Any]], species: str, index: int, count: int) -> list[dict[str, Any]]:
    selected = [row for row in rows if str(row.get("species", "")).strip() == species]
    selected.sort(key=lambda row: str(row.get("photo_id", "")))
    if len(selected) % count:
        raise ValueError(f"{species}: {len(selected)} rows are not divisible by shard_count={count}")
    size = len(selected) // count
    start = index * size
    return selected[start : start + size]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--boxes-jsonl", type=Path, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--split", choices=("calibration", "evaluation"), required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output-jsonl", type=Path, required=True)
    args = parser.parse_args()

    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index/count")
    rows = stable_shard(read_rows(args.boxes_jsonl), args.species, args.shard_index, args.shard_count)
    if not rows:
        raise ValueError("selected colour-measurement shard is empty")

    output: list[dict[str, Any]] = []
    failures = 0
    for position, source in enumerate(rows, start=1):
        photo_id = str(source.get("photo_id", "")).strip()
        row = {
            "protocol": "jbi_ch1_frozen_florence_roi_cielab_v1",
            "species": args.species,
            "split": args.split,
            "photo_id": photo_id,
            "blind_id": source.get("blind_id"),
            "source_box_feature_status": source.get("feature_status"),
            "continuous_colour_representation": "fixed_inner_ellipse_trimmed_mean_cielab_d65",
            "continuous_colour_dimension": 3,
            "continuous_colour_component_names": ["L_star", "a_star", "b_star"],
            "measurement_uses_species_outcome": False,
            "measurement_uses_coordinates": False,
            "measurement_uses_graph": False,
            "measurement_rule_tuned_on_evaluation": False,
        }
        try:
            if not photo_id:
                raise ValueError("photo_id is empty")
            box, box_path = find_box(source)
            image, image_sha = download_image(str(source.get("photo_url", "")))
            normalized_box = normalize_box(box, image.width, image.height, box_path)
            vector, diagnostics = measure(image, normalized_box)
            row.update(
                {
                    "colour_feature_status": "ok",
                    "continuous_colour_vector": vector.tolist(),
                    "source_box_path": box_path,
                    "source_box_xyxy": list(normalized_box),
                    "source_image_width_px": image.width,
                    "source_image_height_px": image.height,
                    "source_image_sha256": image_sha,
                    **diagnostics,
                }
            )
        except Exception as exc:  # row-level diagnostics are retained; workflow later fails closed.
            failures += 1
            row.update(
                {
                    "colour_feature_status": "failed",
                    "colour_feature_error_type": type(exc).__name__,
                    "colour_feature_error": str(exc),
                }
            )
        output.append(row)
        if position % 5 == 0 or position == len(rows):
            print(
                f"{args.species} {args.split} shard {args.shard_index}/{args.shard_count}: "
                f"{position}/{len(rows)} measured; failures={failures}",
                flush=True,
            )

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.output_jsonl.open("w", encoding="utf-8") as handle:
        for row in output:
            handle.write(json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n")
    if len(output) != len(rows):
        raise RuntimeError("output row count changed")
    if failures:
        raise SystemExit(f"{failures}/{len(rows)} colour measurements failed; diagnostics retained")


if __name__ == "__main__":
    main()
