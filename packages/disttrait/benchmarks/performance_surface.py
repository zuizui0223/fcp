#!/usr/bin/env python3
"""Performance surface for species-conditioned spatial inference.

The benchmark varies within-species effect size, between-species observation
imbalance and MCAR observation loss under a fixed strong between-species
geographic turnover process.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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


def run_world(
    *,
    seed: int,
    effect_size: float,
    imbalance_ratio: float,
    missing_fraction: float,
    n_species: int = N_SPECIES,
    base_observations: int = BASE_OBSERVATIONS,
    n_permutations: int = N_PERMUTATIONS,
    pooled_pair_sample: int = POOLED_PAIR_SAMPLE,
) -> dict[str, float]:
    rng = np.random.default_rng(int(seed))
    all_x: list[np.ndarray] = []
    all_state: list[np.ndarray] = []
    observed: list[float] = []
    nulls: list[np.ndarray] = []
    n_retained: list[int] = []

    centres = np.linspace(-5.0, 5.0, int(n_species))
    weights = np.linspace(1.0, float(imbalance_ratio), int(n_species))
    counts = np.maximum(8, np.rint(int(base_observations) * weights).astype(int))

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
            n_permutations=int(n_permutations),
            seed=99123,
            key=f"sp{i}|{seed}",
        )
        observed.append(rho)
        nulls.append(null)
        all_x.append(x)
        all_state.append(state)
        n_retained.append(len(x))

    conditioned = species_equal_spatial_omnibus(observed, np.vstack(nulls))

    x = np.concatenate(all_x)
    state = np.concatenate(all_state)
    n = len(x)
    left = rng.integers(0, n, int(pooled_pair_sample))
    right = rng.integers(0, n, int(pooled_pair_sample))
    keep = left != right
    left, right = left[keep], right[keep]
    pooled_distance = np.abs(x[left] - x[right])
    pooled_difference = (state[left] != state[right]).astype(float)
    naive = spearmanr(pooled_distance, pooled_difference)

    return {
        "naive_rho": float(naive.statistic),
        "naive_p": float(naive.pvalue),
        "conditioned_rho": float(conditioned.observed),
        "conditioned_p": float(conditioned.p_upper),
        "mean_observations": float(np.mean(n_retained)),
        "min_observations": int(np.min(n_retained)),
        "max_observations": int(np.max(n_retained)),
    }


def run_surface(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
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
                        "conditioned_detection": float(np.mean(frame["conditioned_p"] < ALPHA)),
                        "naive_rho_median": float(frame["naive_rho"].median()),
                        "conditioned_rho_median": float(frame["conditioned_rho"].median()),
                        "conditioned_p_median": float(frame["conditioned_p"].median()),
                        "mean_observations": float(frame["mean_observations"].mean()),
                        "min_observations": int(frame["min_observations"].min()),
                        "max_observations": int(frame["max_observations"].max()),
                    }
                )

    cell_frame = pd.DataFrame(cells)
    null = cell_frame[cell_frame["effect_size"] == 0.0]
    e08 = cell_frame[cell_frame["effect_size"] == 0.8]
    e16 = cell_frame[cell_frame["effect_size"] == 1.6]
    e24 = cell_frame[cell_frame["effect_size"] == 2.4]

    return {
        "cells": cells,
        "summary": {
            "null_cells": int(len(null)),
            "max_species_conditioned_null_false_positive_fraction": float(
                null["conditioned_detection"].max()
            ),
            "min_naive_null_false_positive_fraction": float(
                null["naive_detection"].min()
            ),
            "effect_0_8_conditioned_detection_range": [
                float(e08["conditioned_detection"].min()),
                float(e08["conditioned_detection"].max()),
            ],
            "effect_1_6_conditioned_detection_range": [
                float(e16["conditioned_detection"].min()),
                float(e16["conditioned_detection"].max()),
            ],
            "effect_2_4_conditioned_detection_range": [
                float(e24["conditioned_detection"].min()),
                float(e24["conditioned_detection"].max()),
            ],
        },
    }


def main() -> int:
    print(json.dumps(run_surface(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
