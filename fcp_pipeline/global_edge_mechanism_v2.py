"""Distance-controlled within-species mechanism concordance for the global atlas.

This module supersedes the raw edge-Spearman statistic in the first prospective
overlay draft. The correction is frozen before any global Monte Carlo colour
field is available. It removes the generic tendency for colour, climate and
community dissimilarity all to increase with geographic separation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import rankdata


@dataclass(frozen=True)
class EdgeMechanismResult:
    predictor: str
    status: str
    mean_species_partial_rho: float
    median_species_partial_rho: float
    positive_species_fraction: float
    null_mean: float
    null_q025: float
    null_q975: float
    p_upper: float
    n_species: int
    permutations: int


def _residualize_on_distance(values: np.ndarray, distance_rank: np.ndarray) -> np.ndarray:
    """Residualize ranked values on ranked distance with an intercept."""
    if values.shape != distance_rank.shape or values.ndim != 1:
        raise ValueError("values and distance_rank must be matching 1D arrays")
    if np.ptp(distance_rank) <= 0:
        return values - float(np.mean(values))
    design = np.column_stack([np.ones(len(values), dtype=float), distance_rank])
    coef, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coef


def partial_spearman_distance(
    colour: Sequence[float],
    external: Sequence[float],
    geographic_distance: Sequence[float],
) -> float:
    """Spearman partial correlation of colour and external contrast given distance."""
    x = np.asarray(colour, dtype=float)
    y = np.asarray(external, dtype=float)
    d = np.asarray(geographic_distance, dtype=float)
    if x.shape != y.shape or x.shape != d.shape or x.ndim != 1:
        raise ValueError("colour, external and distance must be matching 1D arrays")
    if len(x) < 3 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(d)):
        raise ValueError("partial Spearman requires at least three finite observations")
    if np.ptp(x) <= 0 or np.ptp(y) <= 0:
        raise ValueError("colour and external contrast must vary")
    rx = rankdata(x, method="average").astype(float)
    ry = rankdata(y, method="average").astype(float)
    rd = rankdata(d, method="average").astype(float)
    ex = _residualize_on_distance(rx, rd)
    ey = _residualize_on_distance(ry, rd)
    sx = float(np.sqrt(np.sum(ex * ex)))
    sy = float(np.sqrt(np.sum(ey * ey)))
    if sx <= 1e-12 or sy <= 1e-12:
        raise ValueError("partial Spearman undefined after distance residualization")
    return float(np.sum(ex * ey) / (sx * sy))


def _species_partial_rho(
    colour_score: np.ndarray,
    external_score: np.ndarray,
    geographic_distance: np.ndarray,
    species_index: np.ndarray,
    *,
    minimum_edges_per_species: int,
) -> np.ndarray:
    if not (
        colour_score.shape
        == external_score.shape
        == geographic_distance.shape
        == species_index.shape
    ):
        raise ValueError("all edge arrays must have identical shapes")
    minimum_edges_per_species = int(minimum_edges_per_species)
    if minimum_edges_per_species < 5:
        raise ValueError("minimum_edges_per_species must be at least five")
    values: list[float] = []
    for species_id in np.unique(species_index):
        idx = np.flatnonzero(species_index == species_id)
        keep = (
            np.isfinite(colour_score[idx])
            & np.isfinite(external_score[idx])
            & np.isfinite(geographic_distance[idx])
        )
        idx = idx[keep]
        if len(idx) < minimum_edges_per_species:
            continue
        try:
            rho = partial_spearman_distance(
                colour_score[idx], external_score[idx], geographic_distance[idx]
            )
        except ValueError:
            continue
        if np.isfinite(rho):
            values.append(float(rho))
    return np.asarray(values, dtype=float)


def edge_mechanism_permutation_test(
    *,
    predictor_name: str,
    colour_scores: Sequence[float],
    null_colour_scores: Sequence[Sequence[float]],
    external_edge_scores: Sequence[float],
    geographic_edge_distance_km: Sequence[float],
    species_index: Sequence[int],
    minimum_edges_per_species: int = 5,
    minimum_species: int = 30,
) -> EdgeMechanismResult:
    """Equal-species partial-Spearman test controlling great-circle edge distance."""
    colour = np.asarray(colour_scores, dtype=float)
    null = np.asarray(null_colour_scores, dtype=float)
    external = np.asarray(external_edge_scores, dtype=float)
    distance = np.asarray(geographic_edge_distance_km, dtype=float)
    species = np.asarray(species_index, dtype=np.int64)
    if colour.ndim != 1 or not (
        colour.shape == external.shape == distance.shape == species.shape
    ):
        raise ValueError("edge inputs must be matching one-dimensional arrays")
    if np.any(distance < 0):
        raise ValueError("geographic distance cannot be negative")
    if null.ndim != 2 or null.shape[1] != len(colour):
        raise ValueError("null_colour_scores must be permutations x edges")
    minimum_species = int(minimum_species)
    if minimum_species < 2:
        raise ValueError("minimum_species must be at least two")

    observed = _species_partial_rho(
        colour,
        external,
        distance,
        species,
        minimum_edges_per_species=minimum_edges_per_species,
    )
    n_species = int(len(observed))
    if n_species < minimum_species:
        return EdgeMechanismResult(
            predictor=str(predictor_name),
            status="not_evaluable_species_coverage",
            mean_species_partial_rho=float("nan"),
            median_species_partial_rho=float("nan"),
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
        permuted = _species_partial_rho(
            null[i],
            external,
            distance,
            species,
            minimum_edges_per_species=minimum_edges_per_species,
        )
        if len(permuted) != n_species:
            raise ValueError("null permutation changed the set of evaluable species")
        null_stat[i] = float(np.mean(permuted))

    observed_mean = float(np.mean(observed))
    p_upper = float((1 + np.count_nonzero(null_stat >= observed_mean)) / (len(null_stat) + 1))
    return EdgeMechanismResult(
        predictor=str(predictor_name),
        status="evaluated",
        mean_species_partial_rho=observed_mean,
        median_species_partial_rho=float(np.median(observed)),
        positive_species_fraction=float(np.mean(observed > 0)),
        null_mean=float(np.mean(null_stat)),
        null_q025=float(np.quantile(null_stat, 0.025)),
        null_q975=float(np.quantile(null_stat, 0.975)),
        p_upper=p_upper,
        n_species=n_species,
        permutations=int(len(null_stat)),
    )


def holm_adjust(p_values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(k): float(v) for k, v in p_values.items()}
    if not clean:
        return {}
    for value in clean.values():
        if not np.isfinite(value) or value < 0 or value > 1:
            raise ValueError("p-values must be finite and lie in [0, 1]")
    ordered = sorted(clean.items(), key=lambda item: item[1])
    m = len(ordered)
    running = 0.0
    out: dict[str, float] = {}
    for i, (name, value) in enumerate(ordered):
        running = max(running, min(1.0, (m - i) * value))
        out[name] = min(1.0, running)
    return out


def evaluate_edge_mechanism_family(
    *,
    colour_scores: Sequence[float],
    null_colour_scores: Sequence[Sequence[float]],
    external_edge_scores: Mapping[str, Sequence[float]],
    geographic_edge_distance_km: Sequence[float],
    species_index: Sequence[int],
    minimum_edges_per_species: int = 5,
    minimum_species: int = 30,
    alpha: float = 0.05,
) -> dict[str, object]:
    """Evaluate all frozen mechanisms with Holm correction after distance control."""
    alpha = float(alpha)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    results: dict[str, EdgeMechanismResult] = {}
    raw_p: dict[str, float] = {}
    for name, external in external_edge_scores.items():
        result = edge_mechanism_permutation_test(
            predictor_name=name,
            colour_scores=colour_scores,
            null_colour_scores=null_colour_scores,
            external_edge_scores=external,
            geographic_edge_distance_km=geographic_edge_distance_km,
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
            "mean_species_partial_rho": result.mean_species_partial_rho,
            "median_species_partial_rho": result.median_species_partial_rho,
            "positive_species_fraction": result.positive_species_fraction,
            "null_mean": result.null_mean,
            "null_q025": result.null_q025,
            "null_q975": result.null_q975,
            "p_upper": result.p_upper,
            "p_holm": p_holm,
            "supported": bool(
                result.status == "evaluated"
                and result.mean_species_partial_rho > 0
                and p_holm < alpha
            ),
            "n_species": result.n_species,
            "permutations": result.permutations,
        }
    return {
        "status": "evaluated",
        "statistic": "equal-species mean within-species partial Spearman controlling ranked great-circle edge distance",
        "alpha": alpha,
        "multiplicity": "Holm across the fixed primary edge-mechanism family",
        "does_not_require_G1": True,
        "cannot_rescue_null_G1": True,
        "results": payload,
    }


__all__ = [
    "EdgeMechanismResult",
    "edge_mechanism_permutation_test",
    "evaluate_edge_mechanism_family",
    "holm_adjust",
    "partial_spearman_distance",
]
