#!/usr/bin/env python3
"""Model-based comparator surface for species-conditioned spatial inference.

Adds a species fixed-intercept logistic model with a common nonnegative
within-species slope to the v0.4 comparator design. Species intercepts are
profiled out separately for each candidate slope. The one-sided test uses the
0.5*chi-square_1 likelihood-ratio tail for a slope constrained to beta >= 0.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit
from scipy.stats import chi2, spearmanr, ttest_1samp

from disttrait import spatial_permutation_null, species_equal_spatial_omnibus


EFFECT_SIZES = (0.0, 0.8, 1.6, 2.4)
IMBALANCE_RATIOS = (1.0, 4.0)
MISSING_FRACTIONS = (0.0, 0.25, 0.5)
WORLDS_PER_CELL = 40
N_SPECIES = 20
BASE_OBSERVATIONS = 16
N_PERMUTATIONS = 39
POOLED_PAIR_SAMPLE = 8_000
ALPHA = 0.05


def _one_hot_binary(values: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=int)
    return np.column_stack([1 - y, y]).astype(float)


def _cell_seed(effect_size: float, imbalance_ratio: float, missing_fraction: float, world: int) -> int:
    return (
        100_000
        + int(effect_size * 1000) * 10_000
        + int(imbalance_ratio * 10) * 100
        + int(missing_fraction * 100)
        + int(world)
    )


def _weighted_matched_omnibus(
    observed: np.ndarray,
    null: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    obs = float(np.sum(w * np.asarray(observed, dtype=float)))
    null_values = np.sum(w[:, None] * np.asarray(null, dtype=float), axis=0)
    p = float((1 + np.sum(null_values >= obs)) / (len(null_values) + 1))
    return obs, p


def _profile_intercept(
    y: np.ndarray,
    z: np.ndarray,
    beta: float,
    *,
    bound: float = 20.0,
) -> float:
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if np.all(y == 0):
        return -float(bound)
    if np.all(y == 1):
        return float(bound)

    def score(alpha: float) -> float:
        return float(np.sum(y - expit(alpha + float(beta) * z)))

    lo, hi = -float(bound), float(bound)
    if score(lo) < 0:
        return lo
    if score(hi) > 0:
        return hi
    return float(brentq(score, lo, hi))


def _profile_loglikelihood(
    groups: list[tuple[np.ndarray, np.ndarray]],
    beta: float,
) -> float:
    value = 0.0
    for y, z in groups:
        alpha = _profile_intercept(y, z, beta)
        eta = alpha + float(beta) * z
        value += float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    return value


def _fixed_effect_logit_test(
    groups: list[tuple[np.ndarray, np.ndarray]],
    *,
    max_beta: float = 8.0,
) -> tuple[float, float, float]:
    """One-sided profile-LRT for a common beta >= 0 with species intercepts."""
    ll0 = _profile_loglikelihood(groups, 0.0)
    fit = minimize_scalar(
        lambda beta: -_profile_loglikelihood(groups, float(beta)),
        bounds=(0.0, float(max_beta)),
        method="bounded",
        options={"xatol": 1e-6},
    )
    beta_hat = float(fit.x)
    ll1 = -float(fit.fun)
    lrt = max(0.0, 2.0 * (ll1 - ll0))
    p = 1.0 if lrt <= 1e-12 else float(0.5 * chi2.sf(lrt, 1))
    return beta_hat, p, lrt


def run_world(
    *,
    seed: int,
    effect_size: float,
    imbalance_ratio: float,
    missing_fraction: float,
) -> dict[str, float]:
    rng = np.random.default_rng(int(seed))
    all_x: list[np.ndarray] = []
    all_state: list[np.ndarray] = []
    observed: list[float] = []
    nulls: list[np.ndarray] = []
    retained_n: list[int] = []
    model_groups: list[tuple[np.ndarray, np.ndarray]] = []

    centres = np.linspace(-5.0, 5.0, N_SPECIES)
    weights = np.linspace(1.0, float(imbalance_ratio), N_SPECIES)
    counts = np.maximum(8, np.rint(BASE_OBSERVATIONS * weights).astype(int))

    for i, (centre, n0) in enumerate(zip(centres, counts, strict=True)):
        local = rng.normal(0.0, 0.2, int(n0))
        x = centre + local
        z = local / 0.2
        baseline = 1.0 / (1.0 + np.exp(-centre))
        logit = np.log(baseline / (1.0 - baseline))
        probability = 1.0 / (
            1.0 + np.exp(-(logit + float(effect_size) * z))
        )
        state = rng.binomial(1, probability)

        if float(missing_fraction) > 0:
            keep = rng.random(len(x)) >= float(missing_fraction)
            if int(keep.sum()) < 8:
                chosen = rng.choice(len(x), size=8, replace=False)
                keep = np.zeros(len(x), dtype=bool)
                keep[chosen] = True
            x = x[keep]
            z = z[keep]
            state = state[keep]

        rho, null = spatial_permutation_null(
            latitude=np.zeros_like(x),
            longitude=x,
            traits=_one_hot_binary(state),
            n_permutations=N_PERMUTATIONS,
            seed=99123,
            key=f"sp{i}|{seed}",
        )
        observed.append(rho)
        nulls.append(null)
        all_x.append(x)
        all_state.append(state)
        retained_n.append(len(x))
        model_groups.append((state.astype(float), z.astype(float)))

    observed_array = np.asarray(observed, dtype=float)
    null_array = np.vstack(nulls)

    equal = species_equal_spatial_omnibus(observed_array, null_array)

    pair_weights = np.asarray(
        [n * (n - 1) / 2 for n in retained_n],
        dtype=float,
    )
    weighted_rho, weighted_p = _weighted_matched_omnibus(
        observed_array,
        null_array,
        pair_weights,
    )

    species_test = ttest_1samp(
        observed_array,
        popmean=0.0,
        alternative="greater",
    )

    beta_hat, fixed_effect_p, fixed_effect_lrt = _fixed_effect_logit_test(
        model_groups
    )

    x = np.concatenate(all_x)
    state = np.concatenate(all_state)
    n = len(x)
    left = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    right = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    keep = left != right
    left, right = left[keep], right[keep]
    naive = spearmanr(
        np.abs(x[left] - x[right]),
        (state[left] != state[right]).astype(float),
    )

    return {
        "naive_p": float(naive.pvalue),
        "equal_matched_p": float(equal.p_upper),
        "pair_weighted_matched_p": float(weighted_p),
        "species_rho_ttest_p": float(species_test.pvalue),
        "fixed_effect_logit_p": float(fixed_effect_p),
        "fixed_effect_logit_beta": float(beta_hat),
        "fixed_effect_logit_lrt": float(fixed_effect_lrt),
        "naive_rho": float(naive.statistic),
        "equal_rho": float(equal.observed),
        "pair_weighted_rho": float(weighted_rho),
        "species_rho_mean": float(np.mean(observed_array)),
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []
    for effect_size in EFFECT_SIZES:
        for imbalance_ratio in IMBALANCE_RATIOS:
            for missing_fraction in MISSING_FRACTIONS:
                rows = [
                    run_world(
                        seed=_cell_seed(effect_size, imbalance_ratio, missing_fraction, world),
                        effect_size=effect_size,
                        imbalance_ratio=imbalance_ratio,
                        missing_fraction=missing_fraction,
                    )
                    for world in range(int(worlds_per_cell))
                ]
                frame = pd.DataFrame(rows)
                cells.append(
                    {
                        "effect_size": float(effect_size),
                        "imbalance_ratio": float(imbalance_ratio),
                        "missing_fraction": float(missing_fraction),
                        "worlds": int(worlds_per_cell),
                        "naive_detection": float(np.mean(frame["naive_p"] < ALPHA)),
                        "equal_matched_detection": float(
                            np.mean(frame["equal_matched_p"] < ALPHA)
                        ),
                        "pair_weighted_matched_detection": float(
                            np.mean(frame["pair_weighted_matched_p"] < ALPHA)
                        ),
                        "species_rho_ttest_detection": float(
                            np.mean(frame["species_rho_ttest_p"] < ALPHA)
                        ),
                        "fixed_effect_logit_detection": float(
                            np.mean(frame["fixed_effect_logit_p"] < ALPHA)
                        ),
                        "fixed_effect_beta_median": float(
                            frame["fixed_effect_logit_beta"].median()
                        ),
                        "fixed_effect_p_median": float(
                            frame["fixed_effect_logit_p"].median()
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                    }
                )

    table = pd.DataFrame(cells)
    null = table[table["effect_size"] == 0.0]
    e08 = table[table["effect_size"] == 0.8]
    e16 = table[table["effect_size"] == 1.6]
    return {
        "cells": cells,
        "summary": {
            "max_equal_matched_null_fpr": float(
                null["equal_matched_detection"].max()
            ),
            "max_pair_weighted_matched_null_fpr": float(
                null["pair_weighted_matched_detection"].max()
            ),
            "max_species_rho_ttest_null_fpr": float(
                null["species_rho_ttest_detection"].max()
            ),
            "max_fixed_effect_logit_null_fpr": float(
                null["fixed_effect_logit_detection"].max()
            ),
            "mean_fixed_effect_logit_null_fpr": float(
                null["fixed_effect_logit_detection"].mean()
            ),
            "effect_0_8_fixed_effect_detection_range": [
                float(e08["fixed_effect_logit_detection"].min()),
                float(e08["fixed_effect_logit_detection"].max()),
            ],
            "effect_0_8_equal_detection_range": [
                float(e08["equal_matched_detection"].min()),
                float(e08["equal_matched_detection"].max()),
            ],
            "effect_1_6_fixed_effect_detection_range": [
                float(e16["fixed_effect_logit_detection"].min()),
                float(e16["fixed_effect_logit_detection"].max()),
            ],
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
