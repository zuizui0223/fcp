#!/usr/bin/env python3
"""Multivariate continuous-trait benchmark for disttrait.

Each species has a two-dimensional continuous trait vector. A local spatial
effect moves individuals through trait space along a species-specific direction.
The orientation can be shared across species or spread around trait space.

Compared estimands:
1. naive pooled pairwise geographic distance vs Euclidean trait distance;
2. equal-species matched-null multivariate spatial organization;
3. species fixed-intercept common multivariate signed slope vector.

The benchmark asks whether direction-invariant trait-space organization remains
detectable when species differ in the orientation of their multivariate response.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr

from disttrait import (
    multivariate_spatial_permutation_null,
    species_equal_spatial_omnibus,
)


EFFECT_SIZES = (0.0, 0.4, 0.8, 1.2)
ORIENTATION_SPREADS = (0.0, 0.5, 1.0)
MISSING_FRACTIONS = (0.0, 0.5)
WORLDS_PER_CELL = 40
N_SPECIES = 20
BASE_OBSERVATIONS = 20
NOISE_SD = 0.55
N_PERMUTATIONS = 39
POOLED_PAIR_SAMPLE = 8000
ALPHA = 0.05


def _cell_seed(
    effect_size: float,
    orientation_spread: float,
    missing_fraction: float,
    world: int,
) -> int:
    return (
        12_000_000
        + int(effect_size * 1000) * 10_000
        + int(orientation_spread * 100) * 100
        + int(missing_fraction * 100)
        + int(world)
    )


def _orientation_vectors(
    n_species: int,
    spread: float,
    *,
    phase: float,
) -> np.ndarray:
    if float(spread) == 0.0:
        angle = np.full(int(n_species), float(phase), dtype=float)
    else:
        # spread=1 spans a complete circle with balanced response directions;
        # spread=0.5 spans half a circle.
        offset = (
            np.arange(int(n_species), dtype=float) / float(n_species) - 0.5
        )
        angle = float(phase) + 2.0 * np.pi * float(spread) * offset
    return np.column_stack([np.cos(angle), np.sin(angle)])


def _common_multivariate_slope_test(
    groups: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, float]:
    z_blocks: list[np.ndarray] = []
    y_blocks: list[np.ndarray] = []
    n = 0

    for z, y in groups:
        z = np.asarray(z, dtype=float)
        y = np.asarray(y, dtype=float)
        if z.ndim != 1 or y.ndim != 2 or y.shape[0] != len(z) or y.shape[1] != 2:
            raise ValueError("invalid multivariate model group")
        z_blocks.append(z - z.mean())
        y_blocks.append(y - y.mean(axis=0))
        n += len(z)

    zc = np.concatenate(z_blocks)
    yc = np.vstack(y_blocks)
    denom = float(zc @ zc)
    if denom <= 1e-15:
        return np.full(2, np.nan), float("nan")

    beta = np.sum(zc[:, None] * yc, axis=0) / denom
    residual = yc - zc[:, None] * beta[None, :]
    df = int(n - len(groups) - 1)
    sigma = (residual.T @ residual) / float(df)
    covariance = sigma / denom
    inverse = np.linalg.pinv(covariance)
    wald = float(beta @ inverse @ beta)
    p = float(chi2.sf(wald, 2))
    return beta.astype(float), p


def run_world(
    *,
    seed: int,
    effect_size: float,
    orientation_spread: float,
    missing_fraction: float,
) -> dict[str, float]:
    rng = np.random.default_rng(int(seed))
    centres = np.linspace(-5.0, 5.0, N_SPECIES)
    phase = float(rng.uniform(0.0, 2.0 * np.pi))
    directions = _orientation_vectors(
        N_SPECIES,
        orientation_spread,
        phase=phase,
    )

    observed_rho: list[float] = []
    spatial_nulls: list[np.ndarray] = []
    model_groups: list[tuple[np.ndarray, np.ndarray]] = []
    all_position: list[np.ndarray] = []
    all_trait: list[np.ndarray] = []

    for i, centre in enumerate(centres):
        z = rng.normal(0.0, 1.0, BASE_OBSERVATIONS)
        local = 0.35 * z
        position = centre + local

        baseline = np.array([0.8 * centre, -0.6 * centre], dtype=float)
        trait = (
            baseline[None, :]
            + float(effect_size) * z[:, None] * directions[i][None, :]
            + rng.normal(0.0, NOISE_SD, size=(BASE_OBSERVATIONS, 2))
        )

        if float(missing_fraction) > 0:
            keep = rng.random(BASE_OBSERVATIONS) >= float(missing_fraction)
            if int(keep.sum()) < 8:
                chosen = rng.choice(BASE_OBSERVATIONS, size=8, replace=False)
                keep = np.zeros(BASE_OBSERVATIONS, dtype=bool)
                keep[chosen] = True
            z = z[keep]
            position = position[keep]
            trait = trait[keep]

        rho, null = multivariate_spatial_permutation_null(
            latitude=np.zeros_like(position),
            longitude=position,
            values=trait,
            n_permutations=N_PERMUTATIONS,
            seed=20260919,
            key=f"sp{i}|{seed}",
        )
        observed_rho.append(float(rho))
        spatial_nulls.append(null)
        model_groups.append((z, trait))
        all_position.append(position)
        all_trait.append(trait)

    equal = species_equal_spatial_omnibus(
        np.asarray(observed_rho, dtype=float),
        np.vstack(spatial_nulls),
    )

    common_beta, common_p = _common_multivariate_slope_test(model_groups)

    x = np.concatenate(all_position)
    y = np.vstack(all_trait)
    n = len(x)
    left = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    right = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    keep = left != right
    left = left[keep]
    right = right[keep]
    naive = spearmanr(
        np.abs(x[left] - x[right]),
        np.linalg.norm(y[left] - y[right], axis=1),
    )

    return {
        "naive_p": float(naive.pvalue),
        "naive_rho": float(naive.statistic),
        "equal_matched_p": float(equal.p_upper),
        "equal_rho": float(equal.observed),
        "common_vector_p": float(common_p),
        "common_vector_norm": float(np.linalg.norm(common_beta)),
        "common_beta_1": float(common_beta[0]),
        "common_beta_2": float(common_beta[1]),
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []

    for effect_size in EFFECT_SIZES:
        for orientation_spread in ORIENTATION_SPREADS:
            for missing_fraction in MISSING_FRACTIONS:
                rows = [
                    run_world(
                        seed=_cell_seed(
                            effect_size,
                            orientation_spread,
                            missing_fraction,
                            world,
                        ),
                        effect_size=effect_size,
                        orientation_spread=orientation_spread,
                        missing_fraction=missing_fraction,
                    )
                    for world in range(int(worlds_per_cell))
                ]
                frame = pd.DataFrame(rows)
                cells.append(
                    {
                        "effect_size": float(effect_size),
                        "orientation_spread": float(orientation_spread),
                        "missing_fraction": float(missing_fraction),
                        "worlds": int(worlds_per_cell),
                        "naive_detection": float(np.mean(frame["naive_p"] < ALPHA)),
                        "equal_matched_detection": float(
                            np.mean(frame["equal_matched_p"] < ALPHA)
                        ),
                        "common_vector_detection": float(
                            np.mean(frame["common_vector_p"] < ALPHA)
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                        "common_vector_norm_median": float(
                            frame["common_vector_norm"].median()
                        ),
                    }
                )

    table = pd.DataFrame(cells)
    null = table[table["effect_size"] == 0.0]
    weak = table[table["effect_size"] == 0.4]
    strong = table[table["effect_size"] == 0.8]

    def detection_range(
        subset: pd.DataFrame,
        column: str,
        spread: float,
    ) -> list[float]:
        x = subset[subset["orientation_spread"] == spread][column]
        return [float(x.min()), float(x.max())]

    return {
        "cells": cells,
        "summary": {
            "max_equal_matched_null_fpr": float(
                null["equal_matched_detection"].max()
            ),
            "max_common_vector_null_fpr": float(
                null["common_vector_detection"].max()
            ),
            "min_naive_null_fpr": float(
                null["naive_detection"].min()
            ),
            "effect_0_4_spread0_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.0
            ),
            "effect_0_4_spread0_common_detection_range": detection_range(
                weak, "common_vector_detection", 0.0
            ),
            "effect_0_4_spread1_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 1.0
            ),
            "effect_0_4_spread1_common_detection_range": detection_range(
                weak, "common_vector_detection", 1.0
            ),
            "effect_0_8_spread1_equal_detection_range": detection_range(
                strong, "equal_matched_detection", 1.0
            ),
            "effect_0_8_spread1_common_detection_range": detection_range(
                strong, "common_vector_detection", 1.0
            ),
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
