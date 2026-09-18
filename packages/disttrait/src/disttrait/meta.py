"""Species-specific signed slopes and random-effects meta-analysis."""
from __future__ import annotations

import hashlib
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
class PermutationCalibratedSlopeResult:
    observed: "RandomEffectsSlopeResult"
    p_mean_permutation: float
    p_heterogeneity_permutation: float
    p_omnibus_permutation: float
    null_random_mean: np.ndarray
    null_q: np.ndarray


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
    """DerSimonian-Laird random-effects summary with asymptotic component tests.

    The asymptotic p-values are descriptive/model-based and can be imperfect in
    small or heterogeneous species samples. Use random_effects_slope_permutation_test
    when a within-species exchangeability null is scientifically appropriate.

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


def _permutation_seed(seed: int, key: str, species_index: int, permutation_index: int) -> int:
    raw = f"{seed}|{key}|{species_index}|{permutation_index}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "little")


def random_effects_slope_permutation_test(
    groups: Sequence[tuple[Sequence[float], Sequence[float]]],
    *,
    n_permutations: int = 99,
    seed: int = 0,
    key: str = "",
) -> PermutationCalibratedSlopeResult:
    """Calibrate signed-slope mean and heterogeneity by within-species permutation.

    Positions and the complete observed trait values remain fixed within each
    species; only the assignment of trait values to positions is permuted. Each
    null world re-estimates species slopes, sampling variances, the random-effects
    mean and Cochran Q. The two component permutation p-values are combined with
    a Bonferroni correction.

    At least 39 null worlds are recommended; 99 are used by the v0.8 benchmark
    so the Bonferroni omnibus can attain values below 0.05.
    """
    if int(n_permutations) < 1:
        raise ValueError("n_permutations must be positive")
    prepared: list[tuple[np.ndarray, np.ndarray]] = []
    for position, trait in groups:
        x = np.asarray(position, dtype=float)
        y = np.asarray(trait, dtype=float)
        if x.ndim != 1 or y.shape != x.shape or len(x) < 3:
            raise ValueError("every group must contain equal vectors with at least 3 rows")
        if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
            raise ValueError("group values must be finite")
        prepared.append((x, y))

    observed, _ = random_effects_from_groups(prepared)
    null_mean = np.empty(int(n_permutations), dtype=float)
    null_q = np.empty(int(n_permutations), dtype=float)

    for permutation_index in range(int(n_permutations)):
        slopes: list[float] = []
        variances: list[float] = []
        for species_index, (x, y) in enumerate(prepared):
            rng = np.random.default_rng(
                _permutation_seed(seed, key, species_index, permutation_index)
            )
            permuted = y[rng.permutation(len(y))]
            estimate = species_slope_estimate(x, permuted)
            slopes.append(estimate.slope)
            variances.append(estimate.variance)
        summary = random_effects_slope_summary(slopes, variances)
        null_mean[permutation_index] = summary.random_mean
        null_q[permutation_index] = summary.q

    p_mean = float(
        (1 + np.sum(np.abs(null_mean) >= abs(observed.random_mean)))
        / (len(null_mean) + 1)
    )
    p_heterogeneity = float(
        (1 + np.sum(null_q >= observed.q))
        / (len(null_q) + 1)
    )
    p_omnibus = float(min(1.0, 2.0 * min(p_mean, p_heterogeneity)))

    return PermutationCalibratedSlopeResult(
        observed=observed,
        p_mean_permutation=p_mean,
        p_heterogeneity_permutation=p_heterogeneity,
        p_omnibus_permutation=p_omnibus,
        null_random_mean=null_mean,
        null_q=null_q,
    )

