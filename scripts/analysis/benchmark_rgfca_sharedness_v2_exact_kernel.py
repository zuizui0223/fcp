#!/usr/bin/env python3
"""Performance-only qualification of the frozen RGFCA-v2 exact all-pairs axis score.

No observed flower/background colour is read. This benchmark uses synthetic arrays only
and is not a scientific/power outcome. It verifies numerical equivalence of a vectorized
implementation to an explicit pair loop before the full synthetic qualification runner
is allowed to exist.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[2]
MAPPING = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_technical_mapping_v1.json"
AXIS_SEED = 202609070801
THRESHOLDS = np.array([-0.5, -0.25, 0.0, 0.25, 0.5], dtype=float)


def unit_vectors(rng: np.random.Generator, shape: tuple[int, ...]) -> np.ndarray:
    x = rng.normal(size=shape)
    n = np.linalg.norm(x, axis=-1, keepdims=True)
    return x / n


def canonical_axes() -> np.ndarray:
    rng = np.random.default_rng(AXIS_SEED)
    x = unit_vectors(rng, (96, 3))
    for i in range(len(x)):
        for j in range(3):
            if not np.isclose(x[i, j], 0.0):
                if x[i, j] < 0:
                    x[i] *= -1
                break
    return x


def exact_axis_scores_vectorized(xyz: np.ndarray, flower: np.ndarray, background: np.ndarray) -> np.ndarray:
    """Return 96x5 exact raw scores, using every retained unordered pair."""
    n = len(xyz)
    if xyz.shape != (n, 3) or flower.shape != (n, 3) or background.shape != (n, 3):
        raise ValueError("expected n x 3 arrays")
    c = cdist(flower, flower, metric="euclidean") - cdist(background, background, metric="euclidean")
    p = xyz @ canonical_axes().T
    out = np.full((96, 5), np.nan, dtype=float)
    for ti, threshold in enumerate(THRESHOLDS):
        g = (p > threshold).astype(float)
        ng = g.sum(axis=0)
        nh = n - ng
        valid = (ng > 0) & (nh > 0)
        if not np.any(valid):
            continue
        cg = c @ g
        total_row = c.sum(axis=1, keepdims=True)
        ch = total_row - cg
        cross_sum = np.sum(g * ch, axis=0)
        same_sum = 0.5 * (np.sum(g * cg, axis=0) + np.sum((1.0 - g) * ch, axis=0))
        cross_n = ng * nh
        same_n = ng * (ng - 1.0) / 2.0 + nh * (nh - 1.0) / 2.0
        valid &= same_n > 0
        out[valid, ti] = cross_sum[valid] / cross_n[valid] - same_sum[valid] / same_n[valid]
    return out


def exact_axis_scores_reference(xyz: np.ndarray, flower: np.ndarray, background: np.ndarray) -> np.ndarray:
    axes = canonical_axes()
    n = len(xyz)
    out = np.full((96, 5), np.nan, dtype=float)
    for k, axis in enumerate(axes):
        projection = xyz @ axis
        for ti, threshold in enumerate(THRESHOLDS):
            side = projection > threshold
            cross, same = [], []
            for i in range(n - 1):
                for j in range(i + 1, n):
                    d = float(np.linalg.norm(flower[i] - flower[j]) - np.linalg.norm(background[i] - background[j]))
                    (cross if bool(side[i]) != bool(side[j]) else same).append(d)
            if cross and same:
                out[k, ti] = float(np.mean(cross) - np.mean(same))
    return out


def main() -> int:
    mapping = json.loads(MAPPING.read_text())
    if mapping["status"] != "frozen_after_complete_150x300_geometry_before_any_v2_synthetic_score":
        raise RuntimeError("technical mapping not frozen")
    if mapping["pair_computation"]["all_retained_pairs_required"] is not True:
        raise RuntimeError("all-pairs requirement drift")
    rng = np.random.default_rng(202609081751)
    # Numerical check on a deliberately small case.
    n0 = 18
    xyz0 = unit_vectors(rng, (n0, 3))
    f0 = rng.normal(size=(n0, 3))
    b0 = rng.normal(size=(n0, 3))
    a = exact_axis_scores_vectorized(xyz0, f0, b0)
    r = exact_axis_scores_reference(xyz0, f0, b0)
    keep = np.isfinite(a) & np.isfinite(r)
    if not np.array_equal(np.isfinite(a), np.isfinite(r)) or not np.allclose(a[keep], r[keep], atol=1e-11, rtol=1e-11):
        raise RuntimeError("vectorized exact score disagrees with explicit all-pairs reference")
    timings = {}
    finite_cells = {}
    for n in (120, 180, 240):
        xyz = unit_vectors(rng, (n, 3))
        flower = rng.normal(size=(n, 3))
        background = rng.normal(size=(n, 3))
        t0 = time.perf_counter()
        score = exact_axis_scores_vectorized(xyz, flower, background)
        timings[str(n)] = float(time.perf_counter() - t0)
        finite_cells[str(n)] = int(np.isfinite(score).sum())
        if not np.isfinite(score).any():
            raise RuntimeError(f"no finite axis scores at n={n}")
    payload = {
        "status": "performance_only_exact_kernel_qualified",
        "numerical_reference_max_abs_error": float(np.max(np.abs(a[keep] - r[keep]))),
        "retained_photo_counts": [120, 180, 240],
        "seconds_single_species_single_world": timings,
        "finite_axis_threshold_cells": finite_cells,
        "observed_colour_opened": False,
        "image_pixels_opened": False,
        "scientific_outcome_generated": False,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
