#!/usr/bin/env python3
"""Synthetic benchmark for species-conditioning under geographic turnover.

This benchmark deliberately targets one failure mode: pooling observations across
species can confound between-species geographic turnover with within-species
spatial trait organization.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import spearmanr

from disttrait import (
    spatial_permutation_null,
    species_equal_spatial_omnibus,
)


def _one_hot_binary(values: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=int)
    return np.column_stack([1 - y, y]).astype(float)


def run_world(
    *,
    seed: int,
    signal: bool,
    n_species: int = 24,
    observations_per_species: int = 18,
    n_permutations: int = 39,
    pooled_pair_sample: int = 10_000,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    all_x: list[np.ndarray] = []
    all_state: list[np.ndarray] = []
    observed: list[float] = []
    nulls: list[np.ndarray] = []

    for i, centre in enumerate(np.linspace(-5.0, 5.0, int(n_species))):
        local = rng.normal(0.0, 0.2, int(observations_per_species))
        x = centre + local
        baseline = 1.0 / (1.0 + np.exp(-centre))

        if signal:
            logit = np.log(baseline / (1.0 - baseline))
            probability = 1.0 / (1.0 + np.exp(-(logit + 2.2 * local / 0.2)))
        else:
            # No within-species spatial dependence. Only the species baseline
            # changes along the between-species geographic gradient.
            probability = np.full(int(observations_per_species), baseline)

        state = rng.binomial(1, probability)
        rho, null = spatial_permutation_null(
            latitude=np.zeros_like(x),
            longitude=x,
            traits=_one_hot_binary(state),
            n_permutations=n_permutations,
            seed=99123,
            key=f"sp{i}|{seed}",
        )
        observed.append(rho)
        nulls.append(null)
        all_x.append(x)
        all_state.append(state)

    conditioned = species_equal_spatial_omnibus(
        observed,
        np.vstack(nulls),
    )

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
        "naive_pooled_rho": float(naive.statistic),
        "naive_pooled_p": float(naive.pvalue),
        "species_conditioned_mean_rho": float(conditioned.observed),
        "species_conditioned_p": float(conditioned.p_upper),
    }


def run_benchmark(
    *,
    worlds: int = 40,
    null_seed_base: int = 1000,
    signal_seed_base: int = 2000,
) -> dict:
    null = [run_world(seed=null_seed_base + i, signal=False) for i in range(int(worlds))]
    signal = [run_world(seed=signal_seed_base + i, signal=True) for i in range(int(worlds))]

    def matrix(rows: list[dict[str, float]]) -> np.ndarray:
        return np.asarray(
            [
                [
                    row["naive_pooled_rho"],
                    row["naive_pooled_p"],
                    row["species_conditioned_mean_rho"],
                    row["species_conditioned_p"],
                ]
                for row in rows
            ],
            dtype=float,
        )

    null_m = matrix(null)
    signal_m = matrix(signal)
    return {
        "worlds_per_condition": int(worlds),
        "null": {
            "naive_pooled_false_positive_fraction": float(np.mean(null_m[:, 1] < 0.05)),
            "species_conditioned_false_positive_fraction": float(np.mean(null_m[:, 3] < 0.05)),
            "naive_pooled_rho_median": float(np.median(null_m[:, 0])),
            "species_conditioned_mean_rho_median": float(np.median(null_m[:, 2])),
            "species_conditioned_p_median": float(np.median(null_m[:, 3])),
        },
        "signal": {
            "naive_pooled_detection_fraction": float(np.mean(signal_m[:, 1] < 0.05)),
            "species_conditioned_detection_fraction": float(np.mean(signal_m[:, 3] < 0.05)),
            "species_conditioned_mean_rho_median": float(np.median(signal_m[:, 2])),
            "species_conditioned_p_median": float(np.median(signal_m[:, 3])),
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
