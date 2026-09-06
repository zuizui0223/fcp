"""External-only ecological/geographic profiles for sympatric and allopatric pairs.

These descriptors are frozen before flower-colour outcomes. They are deliberately
non-causal and nonexclusive. A pair can receive multiple labels, or none.

The same external contrast variables used by G3/G4 are summarized relative to a
pre-frozen reference distribution of all eligible plant-pair opportunities.
Flower colour is never used to create thresholds or labels.

For sympatric pairs, low external distance/turnover can describe shared ecological
context (e.g. climate-similar, edaphic-similar, pollinator-similar).
For allopatric pairs, high external distance or positive barrier intensity can
describe the kind of separation (e.g. climate-separated, marine-gap-separated).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class PairProfile:
    pair_id: str
    geography: str
    labels: tuple[str, ...]
    n_evaluable_axes: int


def _finite_quantile(values: Sequence[float], q: float, *, minimum_reference: int) -> float:
    x = np.asarray(values, dtype=float)
    finite = x[np.isfinite(x)]
    if len(finite) < int(minimum_reference):
        return float("nan")
    return float(np.quantile(finite, float(q)))


def build_pair_context_profiles(
    *,
    pair_ids: Sequence[str],
    geography: Sequence[str],
    pair_scores: Mapping[str, Sequence[float]],
    reference_scores: Mapping[str, Sequence[float]],
    low_similarity_axes: Sequence[str] = (
        "pollinator_turnover",
        "climate_turnover",
        "edaphic_turnover",
        "terrain_turnover",
    ),
    high_separation_axes: Sequence[str] = (
        "pollinator_turnover",
        "climate_turnover",
        "edaphic_turnover",
        "terrain_turnover",
        "marine_gap",
        "major_river_crossing",
        "biogeographic_boundary_crossing",
        "mountain_boundary_crossing",
    ),
    direct_barrier_axes: Sequence[str] = (
        "marine_gap",
        "major_river_crossing",
        "biogeographic_boundary_crossing",
        "mountain_boundary_crossing",
    ),
    low_quantile: float = 0.25,
    high_quantile: float = 0.75,
    minimum_reference: int = 100,
) -> dict[str, object]:
    """Assign external-only context labels to sympatric/allopatric plant pairs.

    ``pair_scores`` and ``reference_scores`` must be outcome-blind. All scores
    increase with ecological turnover or geographic-barrier intensity.

    Sympatric turnover/distance axes receive ``<axis>-similar`` when at or below
    the external-only 25th percentile. Allopatric turnover/distance axes receive
    ``<axis>-separated`` at or above the 75th percentile.

    Zero-inflated explicit barrier axes need one extra guard. If their 75th
    percentile is zero, a zero score cannot be called separated merely because it
    equals that percentile; only a strictly positive barrier score is labelled.
    If Q75 is positive, the usual >=Q75 rule applies. Constant/sparse references
    remain not evaluable.
    """
    ids = np.asarray(pair_ids, dtype=object)
    geo = np.asarray(geography, dtype=object)
    if ids.ndim != 1 or geo.ndim != 1 or len(ids) != len(geo):
        raise ValueError("pair_ids and geography must be matching 1D sequences")
    valid_geo = {"sympatric", "allopatric"}
    if any(str(value) not in valid_geo for value in geo):
        raise ValueError("geography values must be sympatric or allopatric")
    if len(set(str(x) for x in ids)) != len(ids):
        raise ValueError("pair_ids must be unique")
    if not 0 < float(low_quantile) < float(high_quantile) < 1:
        raise ValueError("require 0 < low_quantile < high_quantile < 1")
    minimum_reference = int(minimum_reference)
    if minimum_reference < 10:
        raise ValueError("minimum_reference must be at least 10")

    low_axes = tuple(str(x) for x in low_similarity_axes)
    high_axes = tuple(str(x) for x in high_separation_axes)
    barrier_axes = set(str(x) for x in direct_barrier_axes)
    if not barrier_axes.issubset(set(high_axes)):
        raise ValueError("direct_barrier_axes must be a subset of high_separation_axes")
    requested_axes = tuple(dict.fromkeys((*low_axes, *high_axes)))
    score_arrays: dict[str, np.ndarray] = {}
    thresholds: dict[str, dict[str, float | str]] = {}

    for axis in requested_axes:
        if axis not in pair_scores or axis not in reference_scores:
            thresholds[axis] = {
                "status": "not_evaluable_missing_axis",
                "low": float("nan"),
                "high": float("nan"),
            }
            continue
        scores = np.asarray(pair_scores[axis], dtype=float)
        if scores.ndim != 1 or len(scores) != len(ids):
            raise ValueError(f"pair score for {axis} must align with pair_ids")
        if np.any(scores[np.isfinite(scores)] < 0):
            raise ValueError(f"pair scores for {axis} cannot be negative")
        score_arrays[axis] = scores
        reference = np.asarray(reference_scores[axis], dtype=float)
        finite_reference = reference[np.isfinite(reference)]
        if np.any(finite_reference < 0):
            raise ValueError(f"reference scores for {axis} cannot be negative")
        if len(finite_reference) < minimum_reference:
            thresholds[axis] = {
                "status": "not_evaluable_reference_coverage",
                "low": float("nan"),
                "high": float("nan"),
            }
            continue
        if float(np.ptp(finite_reference)) <= 0:
            thresholds[axis] = {
                "status": "not_evaluable_constant_reference",
                "low": float("nan"),
                "high": float("nan"),
            }
            continue
        thresholds[axis] = {
            "status": "evaluated",
            "low": _finite_quantile(reference, low_quantile, minimum_reference=minimum_reference),
            "high": _finite_quantile(reference, high_quantile, minimum_reference=minimum_reference),
        }

    profiles: dict[str, PairProfile] = {}
    label_counts: dict[str, int] = {}
    for i, pair_id in enumerate(ids):
        labels: list[str] = []
        evaluable_axes = 0
        if str(geo[i]) == "sympatric":
            for axis in low_axes:
                info = thresholds.get(axis, {})
                if info.get("status") != "evaluated" or axis not in score_arrays:
                    continue
                value = float(score_arrays[axis][i])
                if not np.isfinite(value):
                    continue
                evaluable_axes += 1
                if value <= float(info["low"]):
                    label = f"{axis.replace('_turnover', '').replace('_', '-')}-similar"
                    labels.append(label)
                    label_counts[label] = label_counts.get(label, 0) + 1
        else:
            for axis in high_axes:
                info = thresholds.get(axis, {})
                if info.get("status") != "evaluated" or axis not in score_arrays:
                    continue
                value = float(score_arrays[axis][i])
                if not np.isfinite(value):
                    continue
                evaluable_axes += 1
                high = float(info["high"])
                if axis in barrier_axes and high <= 0:
                    separated = value > 0
                else:
                    separated = value >= high
                if separated:
                    base = axis.replace("_turnover", "").replace("_crossing", "").replace("_", "-")
                    label = f"{base}-separated"
                    labels.append(label)
                    label_counts[label] = label_counts.get(label, 0) + 1
        profiles[str(pair_id)] = PairProfile(
            pair_id=str(pair_id),
            geography=str(geo[i]),
            labels=tuple(labels),
            n_evaluable_axes=evaluable_axes,
        )

    return {
        "status": "evaluated",
        "threshold_source": "external-only eligible-pair reference distribution",
        "low_quantile": float(low_quantile),
        "high_quantile": float(high_quantile),
        "direct_barrier_zero_guard": True,
        "thresholds": thresholds,
        "profiles": profiles,
        "label_counts": label_counts,
        "labels_are_nonexclusive": True,
        "labels_are_descriptive_not_causal": True,
    }


__all__ = ["PairProfile", "build_pair_context_profiles"]
