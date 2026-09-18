#!/usr/bin/env python3
"""Comparator benchmark for species-conditioned spatial inference.

Compares:
1. naive pooled pair analysis;
2. equal-species matched-null omnibus;
3. pair-count-weighted matched-null omnibus;
4. one-sample t test on species-specific rho values.

All methods receive the same synthetic worlds. The benchmark is designed to
separate the effect of species-conditioning from the effect of equal species
weighting.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_1samp

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

    centres = np.linspace(-5.0, 5.0, N_SPECIES)
    weights = np.linspace(1.0, float(imbalance_ratio), N_SPECIES)
    counts = np.maximum(8, np.rint(BASE_OBSERVATIONS * weights).astype(int))

    for i, (centre, n0) in enumerate(zip(centres, counts, strict=True)):
        local = rng.normal(0.0, 0.2, int(n0))
        x = centre + local
        baseline = 1.0 / (1.0 + np.exp(-centre))
        logit = np.log(baseline / (1.0 - baseline))
        probability = 1.0 / (
            1.0 + np.exp(-(logit + float(effect_size) * local / 0.2))
        )
        state = rng.binomial(1, probability)

        if float(missing_fraction) > 0:
            keep = rng.random(len(x)) >= float(missing_fraction)
            if int(keep.sum()) < 8:
                chosen = rng.choice(len(x), size=8, replace=False)
                keep = np.zeros(len(x), dtype=bool)
                keep[chosen] = True
            x = x[keep]
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
        "naive_rho": float(naive.statistic),
        "equal_rho": float(equal.observed),
        "pair_weighted_rho": float(weighted_rho),
        "species_rho_mean": float(np.mean(observed_array)),
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float]] = []
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
                        "naive_detection": float(np.mean(frame["naive_p"] < ALPHA)),
                        "equal_matched_detection": float(np.mean(frame["equal_matched_p"] < ALPHA)),
                        "pair_weighted_matched_detection": float(
                            np.mean(frame["pair_weighted_matched_p"] < ALPHA)
                        ),
                        "species_rho_ttest_detection": float(
                            np.mean(frame["species_rho_ttest_p"] < ALPHA)
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                        "pair_weighted_rho_median": float(
                            frame["pair_weighted_rho"].median()
                        ),
                        "ttest_mean_rho_median": float(
                            frame["species_rho_mean"].median()
                        ),
                    }
                )

    table = pd.DataFrame(cells)
    null = table[table["effect_size"] == 0.0]
    weak_imbalanced = table[
        (table["effect_size"] == 0.8) & (table["imbalance_ratio"] == 4.0)
    ]
    return {
        "cells": cells,
        "summary": {
            "max_equal_matched_null_fpr": float(null["equal_matched_detection"].max()),
            "max_pair_weighted_matched_null_fpr": float(
                null["pair_weighted_matched_detection"].max()
            ),
            "max_species_rho_ttest_null_fpr": float(
                null["species_rho_ttest_detection"].max()
            ),
            "min_naive_null_fpr": float(null["naive_detection"].min()),
            "effect_0_8_imbalance4_equal_detection_range": [
                float(weak_imbalanced["equal_matched_detection"].min()),
                float(weak_imbalanced["equal_matched_detection"].max()),
            ],
            "effect_0_8_imbalance4_pair_weighted_detection_range": [
                float(weak_imbalanced["pair_weighted_matched_detection"].min()),
                float(weak_imbalanced["pair_weighted_matched_detection"].max()),
            ],
            "effect_0_8_imbalance4_ttest_detection_range": [
                float(weak_imbalanced["species_rho_ttest_detection"].min()),
                float(weak_imbalanced["species_rho_ttest_detection"].max()),
            ],
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
