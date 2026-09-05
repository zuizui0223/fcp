"""Prospective alignment tests between flower-colour geography and external geography.

Two inferential levels are deliberately separated.

1. Edge mechanism alignment is evaluable whenever adequately replicated species
   exist. It asks whether within-species flower-colour discontinuity is stronger
   on photo-pair graph edges that also show stronger pollinator, climate,
   topographic or explicit-barrier contrast. This analysis does not require a
   shared global colour zone and therefore cannot be used to rescue a null G1.
2. Global field overlay is only evaluable after the recurrent colour field (G1)
   is supported and its geography passes the frozen G2 stability gate.

In both levels, external predictors are frozen independently of flower-colour
outcomes and the same species-conditioned colour permutation null is reused.
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


@dataclass(frozen=True)
class EdgeAlignmentResult:
    predictor: str
    status: str
    mean_species_rho: float
    median_species_rho: float
    positive_species_fraction: float
    null_mean: float
    null_q025: float
    null_q975: float
    p_upper: float
    n_species: int
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


def _species_rho(
    colour_score: np.ndarray,
    external_score: np.ndarray,
    species_index: np.ndarray,
    *,
    minimum_edges_per_species: int,
) -> np.ndarray:
    if colour_score.shape != external_score.shape or colour_score.shape != species_index.shape:
        raise ValueError("colour, external and species_index must have identical shapes")
    minimum_edges_per_species = int(minimum_edges_per_species)
    if minimum_edges_per_species < 3:
        raise ValueError("minimum_edges_per_species must be at least three")
    values: list[float] = []
    for species_id in np.unique(species_index):
        idx = np.flatnonzero(species_index == species_id)
        keep = np.isfinite(colour_score[idx]) & np.isfinite(external_score[idx])
        idx = idx[keep]
        if len(idx) < minimum_edges_per_species:
            continue
        colour = colour_score[idx]
        external = external_score[idx]
        if np.ptp(colour) <= 0 or np.ptp(external) <= 0:
            continue
        rho = float(np.corrcoef(rankdata(colour, method="average"), rankdata(external, method="average"))[0, 1])
        if np.isfinite(rho):
            values.append(rho)
    return np.asarray(values, dtype=float)


def edge_alignment_permutation_test(
    *,
    predictor_name: str,
    colour_scores: Sequence[float],
    null_colour_scores: Sequence[Sequence[float]],
    external_edge_scores: Sequence[float],
    species_index: Sequence[int],
    minimum_edges_per_species: int = 5,
    minimum_species: int = 30,
) -> EdgeAlignmentResult:
    """Equal-species test of colour-edge versus external-edge contrast.

    Each species contributes one Spearman rho, irrespective of its number of graph
    edges. The primary statistic is the arithmetic mean of those species rhos.
    Null rows must contain the precomputed within-species colour permutations used
    by the global analysis, preserving geometry and external predictors exactly.
    """
    colour = np.asarray(colour_scores, dtype=float)
    null = np.asarray(null_colour_scores, dtype=float)
    external = np.asarray(external_edge_scores, dtype=float)
    species = np.asarray(species_index, dtype=np.int64)
    if colour.ndim != 1 or external.shape != colour.shape or species.shape != colour.shape:
        raise ValueError("edge inputs must be matching one-dimensional arrays")
    if null.ndim != 2 or null.shape[1] != colour.shape[0]:
        raise ValueError("null_colour_scores must be permutations x edges")
    minimum_species = int(minimum_species)
    if minimum_species < 2:
        raise ValueError("minimum_species must be at least two")
    observed_species_rho = _species_rho(
        colour,
        external,
        species,
        minimum_edges_per_species=minimum_edges_per_species,
    )
    n_species = int(len(observed_species_rho))
    if n_species < minimum_species:
        return EdgeAlignmentResult(
            predictor=str(predictor_name),
            status="not_evaluable_species_coverage",
            mean_species_rho=float("nan"),
            median_species_rho=float("nan"),
            positive_species_fraction=float("nan"),
            null_mean=float("nan"),
            null_q025=float("nan"),
            null_q975=float("nan"),
            p_upper=float("nan"),
            n_species=n_species,
            permutations=int(null.shape[0]),
        )
    null_stat = np.empty(null.shape[0], dtype=float)
    for i in range(null.shape[0]):
        permuted_rho = _species_rho(
            null[i],
            external,
            species,
            minimum_edges_per_species=minimum_edges_per_species,
        )
        if len(permuted_rho) != n_species:
            raise ValueError("null permutation changed the set of evaluable species")
        null_stat[i] = float(np.mean(permuted_rho))
    observed_mean = float(np.mean(observed_species_rho))
    p_upper = float((1 + np.count_nonzero(null_stat >= observed_mean)) / (len(null_stat) + 1))
    return EdgeAlignmentResult(
        predictor=str(predictor_name),
        status="evaluated",
        mean_species_rho=observed_mean,
        median_species_rho=float(np.median(observed_species_rho)),
        positive_species_fraction=float(np.mean(observed_species_rho > 0)),
        null_mean=float(np.mean(null_stat)),
        null_q025=float(np.quantile(null_stat, 0.025)),
        null_q975=float(np.quantile(null_stat, 0.975)),
        p_upper=p_upper,
        n_species=n_species,
        permutations=int(len(null_stat)),
    )


def overlay_alignment_permutation_test(
    *,
    predictor_name: str,
    observed_field: Sequence[float],
    null_fields: Sequence[Sequence[float]],
    predictor_surface: Sequence[float],
    opportunity: Sequence[float],
    minimum_cells: int = 50,
) -> OverlayAlignmentResult:
    """Test positive field alignment against the original colour-permutation null.

    High predictor values must always mean stronger turnover/barrier intensity before
    calling this function. Predictor direction is frozen upstream rather than
    selected after the observed colour field is seen.
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


def evaluate_primary_edge_mechanisms(
    *,
    colour_scores: Sequence[float],
    null_colour_scores: Sequence[Sequence[float]],
    external_edge_scores: Mapping[str, Sequence[float]],
    species_index: Sequence[int],
    minimum_edges_per_species: int = 5,
    minimum_species: int = 30,
    alpha: float = 0.05,
) -> dict[str, object]:
    """Evaluate the fixed mechanism family without requiring shared-zone support."""
    alpha = float(alpha)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    results: dict[str, EdgeAlignmentResult] = {}
    raw_p: dict[str, float] = {}
    for name, external in external_edge_scores.items():
        result = edge_alignment_permutation_test(
            predictor_name=name,
            colour_scores=colour_scores,
            null_colour_scores=null_colour_scores,
            external_edge_scores=external,
            species_index=species_index,
            minimum_edges_per_species=minimum_edges_per_species,
            minimum_species=minimum_species,
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
            "mean_species_rho": result.mean_species_rho,
            "median_species_rho": result.median_species_rho,
            "positive_species_fraction": result.positive_species_fraction,
            "null_mean": result.null_mean,
            "null_q025": result.null_q025,
            "null_q975": result.null_q975,
            "p_upper": result.p_upper,
            "p_holm": p_holm,
            "supported": bool(result.status == "evaluated" and result.mean_species_rho > 0 and p_holm < alpha),
            "n_species": result.n_species,
            "permutations": result.permutations,
        }
    return {
        "status": "evaluated",
        "alpha": alpha,
        "multiplicity": "Holm across the fixed primary edge-mechanism family",
        "does_not_require_G1": True,
        "cannot_rescue_null_G1": True,
        "results": payload,
    }


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
    """Run the fixed shared-zone overlay family only after the G1+G2 gate."""
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
        "multiplicity": "Holm across the fixed primary shared-zone overlay family",
        "results": payload,
    }


__all__ = [
    "EdgeAlignmentResult",
    "OverlayAlignmentResult",
    "edge_alignment_permutation_test",
    "evaluate_primary_edge_mechanisms",
    "evaluate_primary_overlays",
    "holm_adjust",
    "overlay_alignment_permutation_test",
    "weighted_spearman",
]
