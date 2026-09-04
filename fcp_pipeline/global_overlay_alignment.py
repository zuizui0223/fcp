"""Prospective alignment tests between a supported colour-barrier field and external geography.

The external predictor surfaces are frozen independently of flower-colour outcomes.
Observed and species-conditioned null flower-colour fields are evaluated against the
same predictor surface and the same opportunity weights. This lets the original
within-species colour permutation provide the spatial null for overlay alignment.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import rankdata


@dataclass(frozen=True)
class OverlayAlignmentResult:
    predictor: str
    status: str
    observed_rho: float
    null_mean: float
    null_q025: float
    null_q975: float
    p_upper: float
    n_cells: int
    permutations: int


def _weighted_corr(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    if x.shape != y.shape or x.shape != w.shape:
        raise ValueError("x, y and weights must share a shape")
    if x.ndim != 1:
        raise ValueError("weighted correlation inputs must be one-dimensional")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)) or np.any(~np.isfinite(w)):
        raise ValueError("weighted correlation inputs must be finite")
    if np.any(w <= 0):
        raise ValueError("weights must be strictly positive")
    total = float(w.sum())
    if total <= 0:
        raise ValueError("weights sum to zero")
    mx = float(np.sum(w * x) / total)
    my = float(np.sum(w * y) / total)
    dx = x - mx
    dy = y - my
    vx = float(np.sum(w * dx * dx) / total)
    vy = float(np.sum(w * dy * dy) / total)
    if vx <= 0 or vy <= 0:
        raise ValueError("weighted correlation is undefined for a constant surface")
    cov = float(np.sum(w * dx * dy) / total)
    return cov / float(np.sqrt(vx * vy))


def weighted_spearman(
    first: Sequence[float],
    second: Sequence[float],
    weights: Sequence[float],
) -> float:
    """Opportunity-weighted Spearman correlation using average ranks for ties."""
    x = np.asarray(first, dtype=float)
    y = np.asarray(second, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.shape != y.shape or x.shape != w.shape:
        raise ValueError("first, second and weights must share a shape")
    if x.ndim != 1:
        raise ValueError("inputs must be one-dimensional")
    return float(_weighted_corr(rankdata(x, method="average"), rankdata(y, method="average"), w))


def overlay_alignment_permutation_test(
    *,
    predictor_name: str,
    observed_field: Sequence[float],
    null_fields: Sequence[Sequence[float]],
    predictor_surface: Sequence[float],
    opportunity: Sequence[float],
    minimum_cells: int = 50,
) -> OverlayAlignmentResult:
    """Test positive spatial alignment against the original colour-permutation null.

    High predictor values must always mean stronger turnover/barrier intensity before
    calling this function. Predictor direction is therefore frozen upstream rather
    than selected after the observed colour field is seen.
    """
    observed = np.asarray(observed_field, dtype=float)
    null = np.asarray(null_fields, dtype=float)
    predictor = np.asarray(predictor_surface, dtype=float)
    weights = np.asarray(opportunity, dtype=float)
    if observed.ndim != 1 or predictor.shape != observed.shape or weights.shape != observed.shape:
        raise ValueError("observed, predictor and opportunity must be matching 1D arrays")
    if null.ndim != 2 or null.shape[1] != observed.shape[0]:
        raise ValueError("null_fields must be permutations x cells")
    minimum_cells = int(minimum_cells)
    if minimum_cells < 3:
        raise ValueError("minimum_cells must be at least three")

    mask = (
        np.isfinite(observed)
        & np.isfinite(predictor)
        & np.isfinite(weights)
        & (weights > 0)
    )
    n_cells = int(mask.sum())
    if n_cells < minimum_cells:
        return OverlayAlignmentResult(
            predictor=str(predictor_name),
            status="not_evaluable_external_surface_coverage",
            observed_rho=float("nan"),
            null_mean=float("nan"),
            null_q025=float("nan"),
            null_q975=float("nan"),
            p_upper=float("nan"),
            n_cells=n_cells,
            permutations=int(null.shape[0]),
        )
    if not np.all(np.isfinite(null[:, mask])):
        raise ValueError("null fields must be finite wherever the observed overlay is evaluable")

    observed_rho = weighted_spearman(observed[mask], predictor[mask], weights[mask])
    null_rho = np.empty(null.shape[0], dtype=float)
    for i in range(null.shape[0]):
        null_rho[i] = weighted_spearman(null[i, mask], predictor[mask], weights[mask])
    p_upper = float((1 + np.count_nonzero(null_rho >= observed_rho)) / (len(null_rho) + 1))
    return OverlayAlignmentResult(
        predictor=str(predictor_name),
        status="evaluated",
        observed_rho=float(observed_rho),
        null_mean=float(np.mean(null_rho)),
        null_q025=float(np.quantile(null_rho, 0.025)),
        null_q975=float(np.quantile(null_rho, 0.975)),
        p_upper=p_upper,
        n_cells=n_cells,
        permutations=int(len(null_rho)),
    )


def holm_adjust(p_values: Mapping[str, float]) -> dict[str, float]:
    """Holm family-wise adjustment for a fixed named predictor family."""
    clean = {str(k): float(v) for k, v in p_values.items()}
    if not clean:
        return {}
    for value in clean.values():
        if not np.isfinite(value) or value < 0 or value > 1:
            raise ValueError("p-values must be finite and lie in [0, 1]")
    ordered = sorted(clean.items(), key=lambda item: item[1])
    m = len(ordered)
    adjusted_sorted: list[tuple[str, float]] = []
    running = 0.0
    for i, (name, value) in enumerate(ordered):
        candidate = min(1.0, (m - i) * value)
        running = max(running, candidate)
        adjusted_sorted.append((name, min(1.0, running)))
    return dict(adjusted_sorted)


def evaluate_primary_overlays(
    *,
    g1_supported: bool,
    g2_stable: bool,
    observed_field: Sequence[float],
    null_fields: Sequence[Sequence[float]],
    predictor_surfaces: Mapping[str, Sequence[float]],
    opportunity: Sequence[float],
    minimum_cells: int = 50,
    alpha: float = 0.05,
) -> dict[str, object]:
    """Run the fixed overlay family only after the G1+G2 hierarchical gate."""
    if not bool(g1_supported):
        return {"status": "not_run_g1_hierarchical_gate", "results": {}}
    if not bool(g2_stable):
        return {"status": "not_run_g2_stability_gate", "results": {}}
    alpha = float(alpha)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    results: dict[str, OverlayAlignmentResult] = {}
    raw_p: dict[str, float] = {}
    for name, surface in predictor_surfaces.items():
        result = overlay_alignment_permutation_test(
            predictor_name=name,
            observed_field=observed_field,
            null_fields=null_fields,
            predictor_surface=surface,
            opportunity=opportunity,
            minimum_cells=minimum_cells,
        )
        results[str(name)] = result
        if result.status == "evaluated":
            raw_p[str(name)] = result.p_upper
    adjusted = holm_adjust(raw_p)
    payload: dict[str, object] = {}
    for name, result in results.items():
        p_holm = adjusted.get(name, float("nan"))
        payload[name] = {
            "status": result.status,
            "observed_rho": result.observed_rho,
            "null_mean": result.null_mean,
            "null_q025": result.null_q025,
            "null_q975": result.null_q975,
            "p_upper": result.p_upper,
            "p_holm": p_holm,
            "supported": bool(result.status == "evaluated" and result.observed_rho > 0 and p_holm < alpha),
            "n_cells": result.n_cells,
            "permutations": result.permutations,
        }
    return {
        "status": "evaluated",
        "alpha": alpha,
        "multiplicity": "Holm across the fixed primary overlay family",
        "results": payload,
    }


__all__ = [
    "OverlayAlignmentResult",
    "evaluate_primary_overlays",
    "holm_adjust",
    "overlay_alignment_permutation_test",
    "weighted_spearman",
]
