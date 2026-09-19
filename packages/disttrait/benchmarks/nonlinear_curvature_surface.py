#!/usr/bin/env python3
"""Nonlinear within-species benchmark for disttrait.

Each species has a strong between-species baseline gradient and a local
quadratic trait response. The curvature sign may reverse among species.

Compared estimands:
1. naive pooled pairwise distance-dissimilarity;
2. equal-species matched-null distance-dissimilarity;
3. species fixed-intercept common linear slope;
4. species fixed-intercept common quadratic curvature.

The benchmark asks whether method performance follows the scientific estimand:
linear signed response, shared nonlinear curvature, or direction-invariant
spatial organization.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t

from disttrait import (
    continuous_spatial_permutation_null,
    species_equal_spatial_omnibus,
)


EFFECT_SIZES = (0.0, 0.4, 0.8, 1.2)
CURVATURE_REVERSAL_FRACTIONS = (0.0, 0.25, 0.5)
MISSING_FRACTIONS = (0.0, 0.5)
WORLDS_PER_CELL = 40
N_SPECIES = 20
BASE_OBSERVATIONS = 24
NOISE_SD = 0.35
N_PERMUTATIONS = 39
POOLED_PAIR_SAMPLE = 8000
ALPHA = 0.05


def _cell_seed(
    effect_size: float,
    reversal_fraction: float,
    missing_fraction: float,
    world: int,
) -> int:
    return (
        11_000_000
        + int(effect_size * 1000) * 10_000
        + int(reversal_fraction * 100) * 100
        + int(missing_fraction * 100)
        + int(world)
    )


def _fixed_effect_linear_test(
    groups: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[float, float]:
    z_blocks: list[np.ndarray] = []
    y_blocks: list[np.ndarray] = []
    n = 0
    for z, y in groups:
        z = np.asarray(z, dtype=float)
        y = np.asarray(y, dtype=float)
        z_blocks.append(z - z.mean())
        y_blocks.append(y - y.mean())
        n += len(z)

    zc = np.concatenate(z_blocks)
    yc = np.concatenate(y_blocks)
    denom = float(zc @ zc)
    if denom <= 1e-15:
        return float("nan"), float("nan")
    beta = float((zc @ yc) / denom)
    residual = yc - beta * zc
    df = int(n - len(groups) - 1)
    mse = float((residual @ residual) / df)
    if mse <= 0:
        return beta, 0.0
    se = math.sqrt(mse / denom)
    p = float(2.0 * student_t.sf(abs(beta / se), df))
    return beta, p


def _fixed_effect_quadratic_test(
    groups: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[float, float, float]:
    x_blocks: list[np.ndarray] = []
    y_blocks: list[np.ndarray] = []
    n = 0
    for z, y in groups:
        z = np.asarray(z, dtype=float)
        y = np.asarray(y, dtype=float)
        q = np.square(z)
        x_blocks.append(
            np.column_stack([
                z - z.mean(),
                q - q.mean(),
            ])
        )
        y_blocks.append(y - y.mean())
        n += len(z)

    x = np.vstack(x_blocks)
    yc = np.concatenate(y_blocks)
    xtx = x.T @ x
    beta = np.linalg.solve(xtx, x.T @ yc)
    residual = yc - x @ beta
    df = int(n - len(groups) - 2)
    mse = float((residual @ residual) / df)
    covariance = mse * np.linalg.inv(xtx)
    se_quadratic = math.sqrt(float(covariance[1, 1]))
    t_quadratic = float(beta[1] / se_quadratic)
    p_quadratic = float(2.0 * student_t.sf(abs(t_quadratic), df))
    return float(beta[0]), float(beta[1]), p_quadratic


def run_world(
    *,
    seed: int,
    effect_size: float,
    reversal_fraction: float,
    missing_fraction: float,
) -> dict[str, float]:
    rng = np.random.default_rng(int(seed))
    centres = np.linspace(-5.0, 5.0, N_SPECIES)

    signs = np.ones(N_SPECIES, dtype=float)
    n_reversed = int(round(float(reversal_fraction) * N_SPECIES))
    if n_reversed:
        order = rng.permutation(N_SPECIES)
        signs[order[:n_reversed]] = -1.0

    observed_rho: list[float] = []
    spatial_nulls: list[np.ndarray] = []
    model_groups: list[tuple[np.ndarray, np.ndarray]] = []
    all_position: list[np.ndarray] = []
    all_trait: list[np.ndarray] = []

    for i, centre in enumerate(centres):
        z = rng.uniform(-1.0, 1.0, BASE_OBSERVATIONS)
        local = 0.35 * z
        position = centre + local

        # Center z^2 so effect_size controls curvature without shifting the
        # species-specific mean trait level in expectation.
        curvature = np.square(z) - (1.0 / 3.0)
        baseline = 0.8 * centre
        trait = (
            baseline
            + signs[i] * float(effect_size) * curvature
            + rng.normal(0.0, NOISE_SD, BASE_OBSERVATIONS)
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

        rho, null = continuous_spatial_permutation_null(
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

    linear_beta, linear_p = _fixed_effect_linear_test(model_groups)
    _, quadratic_gamma, quadratic_p = _fixed_effect_quadratic_test(model_groups)

    x = np.concatenate(all_position)
    y = np.concatenate(all_trait)
    n = len(x)
    left = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    right = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    keep = left != right
    left = left[keep]
    right = right[keep]
    naive = spearmanr(
        np.abs(x[left] - x[right]),
        np.abs(y[left] - y[right]),
    )

    return {
        "naive_p": float(naive.pvalue),
        "naive_rho": float(naive.statistic),
        "equal_matched_p": float(equal.p_upper),
        "equal_rho": float(equal.observed),
        "linear_beta": float(linear_beta),
        "linear_p": float(linear_p),
        "quadratic_gamma": float(quadratic_gamma),
        "quadratic_p": float(quadratic_p),
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []

    for effect_size in EFFECT_SIZES:
        for reversal_fraction in CURVATURE_REVERSAL_FRACTIONS:
            for missing_fraction in MISSING_FRACTIONS:
                rows = [
                    run_world(
                        seed=_cell_seed(
                            effect_size,
                            reversal_fraction,
                            missing_fraction,
                            world,
                        ),
                        effect_size=effect_size,
                        reversal_fraction=reversal_fraction,
                        missing_fraction=missing_fraction,
                    )
                    for world in range(int(worlds_per_cell))
                ]
                frame = pd.DataFrame(rows)
                cells.append(
                    {
                        "effect_size": float(effect_size),
                        "reversal_fraction": float(reversal_fraction),
                        "missing_fraction": float(missing_fraction),
                        "worlds": int(worlds_per_cell),
                        "naive_detection": float(np.mean(frame["naive_p"] < ALPHA)),
                        "equal_matched_detection": float(
                            np.mean(frame["equal_matched_p"] < ALPHA)
                        ),
                        "linear_detection": float(
                            np.mean(frame["linear_p"] < ALPHA)
                        ),
                        "quadratic_detection": float(
                            np.mean(frame["quadratic_p"] < ALPHA)
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                        "linear_beta_median": float(frame["linear_beta"].median()),
                        "quadratic_gamma_median": float(
                            frame["quadratic_gamma"].median()
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
        reversal_fraction: float,
    ) -> list[float]:
        x = subset[subset["reversal_fraction"] == reversal_fraction][column]
        return [float(x.min()), float(x.max())]

    return {
        "cells": cells,
        "summary": {
            "max_equal_matched_null_fpr": float(
                null["equal_matched_detection"].max()
            ),
            "max_linear_null_fpr": float(
                null["linear_detection"].max()
            ),
            "max_quadratic_null_fpr": float(
                null["quadratic_detection"].max()
            ),
            "min_naive_null_fpr": float(
                null["naive_detection"].min()
            ),
            "effect_0_4_reversal0_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.0
            ),
            "effect_0_4_reversal0_linear_detection_range": detection_range(
                weak, "linear_detection", 0.0
            ),
            "effect_0_4_reversal0_quadratic_detection_range": detection_range(
                weak, "quadratic_detection", 0.0
            ),
            "effect_0_4_reversal0_5_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.5
            ),
            "effect_0_4_reversal0_5_quadratic_detection_range": detection_range(
                weak, "quadratic_detection", 0.5
            ),
            "effect_0_8_reversal0_5_equal_detection_range": detection_range(
                strong, "equal_matched_detection", 0.5
            ),
            "effect_0_8_reversal0_5_quadratic_detection_range": detection_range(
                strong, "quadratic_detection", 0.5
            ),
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
