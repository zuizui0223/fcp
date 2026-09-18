"""Species-specific signed slopes and random-effects meta-analysis."""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import numpy as np
from scipy.stats import chi2, norm


@dataclass(frozen=True)
class SlopeEstimate:
    slope: float
    variance: float
    standard_error: float
    n: int


@dataclass(frozen=True)
class RandomEffectsSlopeResult:
    n_species: int
    fixed_mean: float
    random_mean: float
    random_mean_se: float
    p_mean_two_sided: float
    q: float
    q_df: int
    p_heterogeneity: float
    tau2: float
    p_omnibus: float


def species_slope_estimate(
    position: Sequence[float],
    trait: Sequence[float],
    *,
    variance_floor: float = 1e-12,
) -> SlopeEstimate:
    """OLS slope and sampling variance for one species."""
    x = np.asarray(position, dtype=float)
    y = np.asarray(trait, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 3:
        raise ValueError("position and trait must be equal vectors with at least 3 rows")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError("position and trait must be finite")
    xc = x - x.mean()
    yc = y - y.mean()
    denom = float(xc @ xc)
    if denom <= 1e-15:
        raise ValueError("position has no within-species variation")
    slope = float((xc @ yc) / denom)
    residual = yc - slope * xc
    df = len(x) - 2
    mse = float((residual @ residual) / df)
    variance = max(float(mse / denom), float(variance_floor))
    return SlopeEstimate(
        slope=slope,
        variance=variance,
        standard_error=float(np.sqrt(variance)),
        n=int(len(x)),
    )


def random_effects_slope_summary(
    slopes: Sequence[float],
    variances: Sequence[float],
) -> RandomEffectsSlopeResult:
    """DerSimonian-Laird random-effects summary plus a mean/heterogeneity omnibus.

    The omnibus combines two distinct signed-slope questions:
    (i) whether the average signed slope is nonzero and
    (ii) whether among-species slope heterogeneity exceeds within-species
    sampling uncertainty. The final probability is a Bonferroni combination
    of those two tests and is deliberately conservative.
    """
    b = np.asarray(slopes, dtype=float)
    v = np.asarray(variances, dtype=float)
    if b.ndim != 1 or v.shape != b.shape or len(b) < 2:
        raise ValueError("slopes and variances must be equal vectors for >=2 species")
    if np.any(~np.isfinite(b)) or np.any(~np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("slopes must be finite and variances finite positive")

    w = 1.0 / v
    sw = float(w.sum())
    fixed_mean = float(np.sum(w * b) / sw)
    q = float(np.sum(w * np.square(b - fixed_mean)))
    q_df = int(len(b) - 1)
    p_heterogeneity = float(chi2.sf(q, q_df))
    c = float(sw - np.sum(np.square(w)) / sw)
    tau2 = max(0.0, (q - q_df) / c) if c > 0 else 0.0

    wr = 1.0 / (v + tau2)
    swr = float(wr.sum())
    random_mean = float(np.sum(wr * b) / swr)
    random_mean_se = float(np.sqrt(1.0 / swr))
    z = random_mean / random_mean_se
    p_mean = float(2.0 * norm.sf(abs(z)))
    p_omnibus = float(min(1.0, 2.0 * min(p_mean, p_heterogeneity)))

    return RandomEffectsSlopeResult(
        n_species=int(len(b)),
        fixed_mean=fixed_mean,
        random_mean=random_mean,
        random_mean_se=random_mean_se,
        p_mean_two_sided=p_mean,
        q=q,
        q_df=q_df,
        p_heterogeneity=p_heterogeneity,
        tau2=float(tau2),
        p_omnibus=p_omnibus,
    )


def random_effects_from_groups(
    groups: Sequence[tuple[Sequence[float], Sequence[float]]],
) -> tuple[RandomEffectsSlopeResult, list[SlopeEstimate]]:
    """Estimate one signed slope per species and summarize them."""
    estimates = [species_slope_estimate(x, y) for x, y in groups]
    result = random_effects_slope_summary(
        [x.slope for x in estimates],
        [x.variance for x in estimates],
    )
    return result, estimates
