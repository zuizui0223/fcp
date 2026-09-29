"""Generic distance-turnover estimands for trait and interaction data.

This module is deliberately agnostic to flower colour. It provides two primitives
for the cross-program spatiotemporal turnover development lane:

1. a common rank-scale estimand for how biological dissimilarity changes with
   separation in evolutionary time or geographic space; and
2. an additive binary interaction-beta decomposition separating species
   turnover from rewiring among shared species.

These are building blocks, not a new standalone statistical method.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import numpy as np
from scipy.stats import rankdata


@dataclass(frozen=True)
class TurnoverNullResult:
    observed_rho: float
    null_mean_rho: float
    delta_rho: float
    p_upper: float
    null_rho: np.ndarray


@dataclass(frozen=True)
class InteractionBetaResult:
    beta_wn: float
    beta_st: float
    beta_os: float
    shared_interactions: int
    unique_a: int
    unique_b: int
    shared_species_a: int
    shared_species_b: int


def _rank_pearson(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(np.asarray(x, dtype=float), method="average")
    ry = rankdata(np.asarray(y, dtype=float), method="average")
    xc = rx - rx.mean()
    yc = ry - ry.mean()
    den = float(np.linalg.norm(xc) * np.linalg.norm(yc))
    return float(np.dot(xc, yc) / den) if den > 1e-15 else float("nan")


def turnover_rho(
    separation: Sequence[float],
    dissimilarity: Sequence[float],
) -> float:
    """Spearman association between separation and biological dissimilarity.

    Higher values mean faster turnover with separation. Separation can be
    phylogenetic distance, geographic distance, temporal lag, or another
    predeclared distance axis. Dissimilarity can be trait distance,
    partner-profile distance, or a network dissimilarity.

    The function is intentionally scale-free. Cross-domain comparisons should
    compare direction and standardized effect structure rather than assume that
    raw rho values from different biological representations are metrically
    interchangeable.
    """
    x = np.asarray(separation, dtype=float)
    y = np.asarray(dissimilarity, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape or len(x) < 3:
        raise ValueError("separation and dissimilarity must be equal 1D vectors with n>=3")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError("inputs must be finite")
    if np.ptp(x) <= 1e-15:
        raise ValueError("separation has no variation")
    if np.ptp(y) <= 1e-15:
        return 0.0
    return _rank_pearson(x, y)


def null_centered_turnover_rho(
    separation: Sequence[float],
    dissimilarity: Sequence[float],
    null_dissimilarities: np.ndarray,
) -> TurnoverNullResult:
    """Return null-centered distance-turnover effect on a common rho scale.

    Each row of null_dissimilarities must be one null world on the same
    pairwise geometry as dissimilarity. The biological effect is delta_rho:
    observed Spearman rho minus the mean null rho. Null SD is intentionally
    not used as the cross-system biological effect size.
    """
    observed = turnover_rho(separation, dissimilarity)
    null = np.asarray(null_dissimilarities, dtype=float)
    y = np.asarray(dissimilarity, dtype=float)
    if null.ndim != 2 or null.shape[1] != len(y) or null.shape[0] < 1:
        raise ValueError("null_dissimilarities must have shape (n_null>=1, n_pairs)")
    if np.any(~np.isfinite(null)):
        raise ValueError("null dissimilarities must be finite")
    null_rho = np.asarray([turnover_rho(separation, row) for row in null], dtype=float)
    mean = float(np.mean(null_rho))
    delta = float(observed - mean)
    p_upper = float((1 + np.sum(null_rho >= observed)) / (len(null_rho) + 1))
    return TurnoverNullResult(
        observed_rho=float(observed),
        null_mean_rho=mean,
        delta_rho=delta,
        p_upper=p_upper,
        null_rho=null_rho,
    )


def _sorensen_binary(a: np.ndarray, b: np.ndarray) -> tuple[float, int, int, int]:
    va = np.asarray(a, dtype=bool).ravel()
    vb = np.asarray(b, dtype=bool).ravel()
    if va.shape != vb.shape:
        raise ValueError("binary arrays must have identical shape")
    shared = int(np.sum(va & vb))
    only_a = int(np.sum(va & ~vb))
    only_b = int(np.sum(~va & vb))
    den = 2 * shared + only_a + only_b
    beta = 0.0 if den == 0 else float((only_a + only_b) / den)
    return beta, shared, only_a, only_b


def interaction_beta_partition(
    network_a: np.ndarray,
    network_b: np.ndarray,
) -> InteractionBetaResult:
    """Partition binary interaction turnover into species-turnover and rewiring.

    The two matrices must use the same global row and column species universe.
    Presence at a site is inferred from having at least one observed interaction
    in that network.

    beta_wn is whole-network Sorensen dissimilarity. beta_os is the
    common-denominator contribution from interaction changes among species
    present in both networks (rewiring). beta_st is the complementary
    contribution attributable to species turnover. By construction,
    beta_wn equals beta_st plus beta_os up to floating-point precision.

    This follows the additive common-denominator partition described by
    Frund (2021) as an alternative implementation of the Poisot et al. (2012)
    interaction-beta logic.

    The function does not infer species presence independently of observed
    interactions. Applications with separate occurrence or detection data
    should build occurrence-aware aligned matrices upstream.
    """
    a = np.asarray(network_a)
    b = np.asarray(network_b)
    if a.ndim != 2 or b.ndim != 2 or a.shape != b.shape:
        raise ValueError("networks must be equal-shape 2D matrices")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("network matrices must be finite")
    aa = a > 0
    bb = b > 0

    beta_wn, shared, only_a, only_b = _sorensen_binary(aa, bb)
    denom = 2 * shared + only_a + only_b
    if denom == 0:
        return InteractionBetaResult(
            beta_wn=0.0,
            beta_st=0.0,
            beta_os=0.0,
            shared_interactions=0,
            unique_a=0,
            unique_b=0,
            shared_species_a=0,
            shared_species_b=0,
        )

    row_a = aa.any(axis=1)
    row_b = bb.any(axis=1)
    col_a = aa.any(axis=0)
    col_b = bb.any(axis=0)
    shared_rows = row_a & row_b
    shared_cols = col_a & col_b

    shared_mask = np.outer(shared_rows, shared_cols)
    b_os = int(np.sum(aa & ~bb & shared_mask))
    c_os = int(np.sum(~aa & bb & shared_mask))

    beta_os = float((b_os + c_os) / denom)
    beta_st = float(beta_wn - beta_os)

    if beta_st < -1e-12 or beta_os < -1e-12:
        raise RuntimeError("invalid negative beta component")
    beta_st = max(0.0, beta_st)
    beta_os = max(0.0, beta_os)

    return InteractionBetaResult(
        beta_wn=float(beta_wn),
        beta_st=beta_st,
        beta_os=beta_os,
        shared_interactions=shared,
        unique_a=only_a,
        unique_b=only_b,
        shared_species_a=int(shared_rows.sum()),
        shared_species_b=int(shared_cols.sum()),
    )
