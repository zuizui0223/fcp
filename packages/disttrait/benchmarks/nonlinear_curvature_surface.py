#!/usr/bin/env python3
"""Nonlinear within-species benchmark for disttrait.

Each species has strong between-species geographic turnover plus a symmetric
quadratic within-species trait pattern. The sign of the curvature can reverse
among species.

Compared estimands:
1. naive pooled distance-difference correlation;
2. direction-invariant species-conditioned distance-dissimilarity matched null;
3. one common signed linear within-species slope;
4. one common signed quadratic coefficient;
5. species-specific quadratic coefficients summarized by the permutation-
   calibrated random-effects slope/meta layer.

The benchmark separates nonlinear model specification from species-specific
direction heterogeneity.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t

from disttrait import (
    continuous_spatial_permutation_null,
    random_effects_slope_permutation_test,
    species_equal_spatial_omnibus,
)


EFFECT_SIZES = (0.0, 0.4, 0.8, 1.2)
REVERSAL_FRACTIONS = (0.0, 0.25, 0.5)
MISSING_FRACTIONS = (0.0, 0.5)
WORLDS_PER_CELL = 40
N_SPECIES = 20
BASE_OBSERVATIONS = 20
N_SPATIAL_PERMUTATIONS = 39
N_META_PERMUTATIONS = 99
POOLED_PAIR_SAMPLE = 8_000
ALPHA = 0.05


def _cell_seed(
    effect_size: float,
    reversal_fraction: float,
    missing_fraction: float,
    world: int,
) -> int:
    return (
        10_000_000
        + int(effect_size * 1000) * 10_000
        + int(reversal_fraction * 100) * 100
        + int(missing_fraction * 100)
        + int(world)
    )


def _common_predictor_test(
    groups: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[float, float]:
    x_blocks: list[np.ndarray] = []
    y_blocks: list[np.ndarray] = []
    n = 0
    for x, y in groups:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.shape != y.shape or x.ndim != 1 or len(x) < 3:
            raise ValueError("invalid species group")
        x_blocks.append(x - x.mean())
        y_blocks.append(y - y.mean())
        n += len(x)

    xc = np.concatenate(x_blocks)
    yc = np.concatenate(y_blocks)
    denom = float(xc @ xc)
    if denom <= 1e-15:
        return float("nan"), float("nan")
    beta = float((xc @ yc) / denom)
    residual = yc - beta * xc
    df = int(n - len(groups) - 1)
    mse = float((residual @ residual) / df)
    if mse <= 0:
        return beta, 0.0 if beta > 0 else 1.0
    se = math.sqrt(mse / denom)
    return beta, float(student_t.sf(beta / se, df))


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
        reversed_order = rng.permutation(N_SPECIES)
        signs[reversed_order[:n_reversed]] = -1.0

    all_position: list[np.ndarray] = []
    all_trait: list[np.ndarray] = []
    observed_rho: list[float] = []
    spatial_nulls: list[np.ndarray] = []
    linear_groups: list[tuple[np.ndarray, np.ndarray]] = []
    quadratic_groups: list[tuple[np.ndarray, np.ndarray]] = []

    for i, centre in enumerate(centres):
        z = rng.normal(0.0, 1.0, BASE_OBSERVATIONS)
        local = 0.35 * z
        quadratic = np.square(z)

        baseline = 0.8 * centre
        trait = (
            baseline
            + signs[i] * float(effect_size) * (quadratic - 1.0)
            + rng.normal(0.0, 1.0, BASE_OBSERVATIONS)
        )

        if float(missing_fraction) > 0:
            keep = rng.random(BASE_OBSERVATIONS) >= float(missing_fraction)
            if int(keep.sum()) < 8:
                chosen = rng.choice(BASE_OBSERVATIONS, size=8, replace=False)
                keep = np.zeros(BASE_OBSERVATIONS, dtype=bool)
                keep[chosen] = True
            z = z[keep]
            local = local[keep]
            quadratic = quadratic[keep]
            trait = trait[keep]

        position = centre + local
        rho, null = continuous_spatial_permutation_null(
            latitude=np.zeros_like(position),
            longitude=position,
            values=trait,
            n_permutations=N_SPATIAL_PERMUTATIONS,
            seed=20260919,
            key=f"nonlinear|sp{i}|{seed}",
        )
        observed_rho.append(float(rho))
        spatial_nulls.append(null)
        linear_groups.append((z, trait))
        quadratic_groups.append((quadratic, trait))
        all_position.append(position)
        all_trait.append(trait)

    observed_array = np.asarray(observed_rho, dtype=float)
    spatial_null_array = np.vstack(spatial_nulls)
    equal = species_equal_spatial_omnibus(observed_array, spatial_null_array)

    linear_beta, linear_p = _common_predictor_test(linear_groups)
    quadratic_beta, quadratic_p = _common_predictor_test(quadratic_groups)

    meta = random_effects_slope_permutation_test(
        quadratic_groups,
        n_permutations=N_META_PERMUTATIONS,
        seed=20260919,
        key=f"nonlinear_quadratic|{seed}",
    )

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
        "common_linear_beta": float(linear_beta),
        "common_linear_p": float(linear_p),
        "common_quadratic_beta": float(quadratic_beta),
        "common_quadratic_p": float(quadratic_p),
        "meta_mean_p": float(meta.p_mean_permutation),
        "meta_heterogeneity_p": float(meta.p_heterogeneity_permutation),
        "meta_omnibus_p": float(meta.p_omnibus_permutation),
        "meta_random_mean": float(meta.observed.random_mean),
        "meta_tau2": float(meta.observed.tau2),
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []

    for effect_size in EFFECT_SIZES:
        for reversal_fraction in REVERSAL_FRACTIONS:
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
                        "common_linear_detection": float(
                            np.mean(frame["common_linear_p"] < ALPHA)
                        ),
                        "common_quadratic_detection": float(
                            np.mean(frame["common_quadratic_p"] < ALPHA)
                        ),
                        "meta_mean_detection": float(
                            np.mean(frame["meta_mean_p"] < ALPHA)
                        ),
                        "meta_heterogeneity_detection": float(
                            np.mean(frame["meta_heterogeneity_p"] < ALPHA)
                        ),
                        "meta_omnibus_detection": float(
                            np.mean(frame["meta_omnibus_p"] < ALPHA)
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                        "common_linear_beta_median": float(
                            frame["common_linear_beta"].median()
                        ),
                        "common_quadratic_beta_median": float(
                            frame["common_quadratic_beta"].median()
                        ),
                        "meta_random_mean_median": float(
                            frame["meta_random_mean"].median()
                        ),
                        "meta_tau2_median": float(frame["meta_tau2"].median()),
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
            "max_common_linear_null_fpr": float(
                null["common_linear_detection"].max()
            ),
            "max_common_quadratic_null_fpr": float(
                null["common_quadratic_detection"].max()
            ),
            "max_meta_omnibus_null_fpr": float(
                null["meta_omnibus_detection"].max()
            ),
            "effect_0_4_reversal0_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.0
            ),
            "effect_0_4_reversal0_linear_detection_range": detection_range(
                weak, "common_linear_detection", 0.0
            ),
            "effect_0_4_reversal0_quadratic_detection_range": detection_range(
                weak, "common_quadratic_detection", 0.0
            ),
            "effect_0_4_reversal0_meta_detection_range": detection_range(
                weak, "meta_omnibus_detection", 0.0
            ),
            "effect_0_4_reversal0_5_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.5
            ),
            "effect_0_4_reversal0_5_quadratic_detection_range": detection_range(
                weak, "common_quadratic_detection", 0.5
            ),
            "effect_0_4_reversal0_5_meta_heterogeneity_detection_range": detection_range(
                weak, "meta_heterogeneity_detection", 0.5
            ),
            "effect_0_4_reversal0_5_meta_omnibus_detection_range": detection_range(
                weak, "meta_omnibus_detection", 0.5
            ),
            "effect_0_8_reversal0_5_equal_detection_range": detection_range(
                strong, "equal_matched_detection", 0.5
            ),
            "effect_0_8_reversal0_5_meta_omnibus_detection_range": detection_range(
                strong, "meta_omnibus_detection", 0.5
            ),
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
