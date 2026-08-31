"""Independent flower-tissue localization metrics for atlas qualification."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
from typing import Any

import numpy as np
from skimage.color import rgb2lab


def score_flower_weights(
    rgb: np.ndarray,
    weights: np.ndarray,
    trimap: np.ndarray,
    *,
    foreground_label: int = 1,
    background_label: int = 2,
) -> dict[str, float]:
    """Score soft flower weights against labelled Oxford trimap pixels.

    Labels other than the frozen flower foreground and background labels are excluded
    from every localization metric.  This follows the benchmark's unlabelled-region
    convention and prevents ambiguous boundary pixels from being counted as errors.
    """

    rgb = np.asarray(rgb)
    weights = np.asarray(weights, dtype=float)
    trimap = np.asarray(trimap)
    if rgb.ndim != 3 or rgb.shape[2] != 3:
        raise ValueError("rgb must have shape (height, width, 3)")
    if weights.ndim != 2 or trimap.ndim != 2:
        raise ValueError("weights and trimap must be two-dimensional")
    if rgb.shape[:2] != weights.shape or weights.shape != trimap.shape:
        raise ValueError("rgb, weights, and trimap dimensions must match")
    if not np.isfinite(weights).all() or np.any((weights < 0.0) | (weights > 1.0)):
        raise ValueError("weights must be finite and lie in [0, 1]")

    foreground = trimap == int(foreground_label)
    background = trimap == int(background_label)
    scored = foreground | background
    n_foreground = int(foreground.sum())
    n_background = int(background.sum())
    if n_foreground < 1 or n_background < 1:
        raise ValueError("trimap must contain scored foreground and background pixels")

    scored_weights = np.where(scored, weights, 0.0)
    foreground_weight = float(scored_weights[foreground].sum())
    total_weight = float(scored_weights.sum())
    soft_precision = foreground_weight / total_weight if total_weight > 0 else 0.0
    soft_recall = foreground_weight / n_foreground
    denominator = float(np.where(scored, np.maximum(weights, foreground.astype(float)), 0.0).sum())
    soft_iou = foreground_weight / denominator if denominator > 0 else 0.0

    lab = rgb2lab(rgb.astype(np.float32) / 255.0)
    gold_lab = np.mean(lab[foreground], axis=0)
    if total_weight > 0:
        estimated_lab = np.sum(lab * scored_weights[..., None], axis=(0, 1))
        estimated_lab /= total_weight
        colour_delta_e = float(np.linalg.norm(estimated_lab - gold_lab))
    else:
        colour_delta_e = math.inf

    return {
        "foreground_pixels": float(n_foreground),
        "background_pixels": float(n_background),
        "scored_weight": total_weight,
        "soft_precision": float(soft_precision),
        "soft_recall": float(soft_recall),
        "soft_iou": float(soft_iou),
        "colour_delta_e": colour_delta_e,
    }


def summarize_roi_benchmark(
    rows: Sequence[Mapping[str, Any]],
    gates: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply the predeclared aggregate ROI gates without threshold searching."""

    if not rows:
        raise ValueError("ROI benchmark has no scored rows")
    admitted = [row for row in rows if bool(row.get("estimator_admitted", False))]
    finite_rows = [
        row
        for row in rows
        if all(
            math.isfinite(float(row[key]))
            for key in ("soft_precision", "soft_recall", "soft_iou", "colour_delta_e")
        )
    ]
    if not finite_rows:
        raise ValueError("ROI benchmark has no finite localization rows")

    def percentile(key: str, q: float) -> float:
        return float(np.percentile([float(row[key]) for row in finite_rows], q))

    metrics = {
        "scored_images": len(rows),
        "finite_images": len(finite_rows),
        "admitted_images": len(admitted),
        "admitted_fraction": len(admitted) / len(rows),
        "median_soft_precision": percentile("soft_precision", 50),
        "median_soft_recall": percentile("soft_recall", 50),
        "median_soft_iou": percentile("soft_iou", 50),
        "median_colour_delta_e": percentile("colour_delta_e", 50),
        "p90_colour_delta_e": percentile("colour_delta_e", 90),
    }
    checks = {
        "minimum_scored_images": metrics["scored_images"]
        >= int(gates["minimum_scored_images"]),
        "minimum_admitted_fraction": metrics["admitted_fraction"]
        >= float(gates["minimum_admitted_fraction"]),
        "minimum_median_soft_precision": metrics["median_soft_precision"]
        >= float(gates["minimum_median_soft_precision"]),
        "minimum_median_soft_recall": metrics["median_soft_recall"]
        >= float(gates["minimum_median_soft_recall"]),
        "minimum_median_soft_iou": metrics["median_soft_iou"]
        >= float(gates["minimum_median_soft_iou"]),
        "maximum_median_colour_delta_e": metrics["median_colour_delta_e"]
        <= float(gates["maximum_median_colour_delta_e"]),
        "maximum_p90_colour_delta_e": metrics["p90_colour_delta_e"]
        <= float(gates["maximum_p90_colour_delta_e"]),
    }
    passed = all(checks.values())
    return {
        "status": "pass_independent_roi_benchmark" if passed else "stop_roi_benchmark_failed",
        "metrics": metrics,
        "checks": checks,
        "atlas_pixels_permitted_by_roi_gate": passed,
    }
