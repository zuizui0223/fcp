#!/usr/bin/env python3
"""Random-effects species-slope benchmark under direction heterogeneity.

Extends the v0.7 continuous-trait benchmark with a model-based summary that
estimates one signed slope per species and then separates:
1. an average signed slope;
2. between-species slope heterogeneity;
3. a conservative omnibus combining those two tests.

The benchmark compares estimands rather than declaring a universal winner.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t, ttest_1samp

from disttrait import (
    continuous_spatial_permutation_null,
    random_effects_from_groups,
    random_effects_slope_permutation_test,
    species_equal_spatial_omnibus,
)


EFFECT_SIZES = (0.0, 0.4, 0.8, 1.2)
REVERSAL_FRACTIONS = (0.0, 0.25, 0.5)
MISSING_FRACTIONS = (0.0, 0.5)
WORLDS_PER_CELL = 40
N_SPECIES = 20
BASE_OBSERVATIONS = 20
N_PERMUTATIONS = 39
META_PERMUTATIONS = 99
POOLED_PAIR_SAMPLE = 8_000
ALPHA = 0.05


def _cell_seed(
    effect_size: float,
    reversal_fraction: float,
    missing_fraction: float,
    world: int,
) -> int:
    return (
        8_000_000
        + int(effect_size * 1000) * 10_000
        + int(reversal_fraction * 100) * 100
        + int(missing_fraction * 100)
        + int(world)
    )


def _fixed_effect_common_slope_test(
    groups: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[float, float]:
    """One-sided common-slope OLS after demeaning within species."""
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
    beta = float((zc @ yc) / denom)
    residual = yc - beta * zc
    df = int(n - len(groups) - 1)
    mse = float((residual @ residual) / df)
    if mse <= 0:
        return beta, 0.0 if beta > 0 else 1.0
    se = math.sqrt(mse / denom)
    p = float(student_t.sf(beta / se, df))
    return beta, p


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
        reverse_order = rng.permutation(N_SPECIES)
        signs[reverse_order[:n_reversed]] = -1.0

    all_position: list[np.ndarray] = []
    all_trait: list[np.ndarray] = []
    observed_rho: list[float] = []
    spatial_nulls: list[np.ndarray] = []
    model_groups: list[tuple[np.ndarray, np.ndarray]] = []

    for i, centre in enumerate(centres):
        local = rng.normal(0.0, 0.35, BASE_OBSERVATIONS)
        z = local / 0.35
        baseline = 0.8 * centre
        trait = (
            baseline
            + signs[i] * float(effect_size) * z
            + rng.normal(0.0, 1.0, BASE_OBSERVATIONS)
        )

        if float(missing_fraction) > 0:
            keep = rng.random(BASE_OBSERVATIONS) >= float(missing_fraction)
            if int(keep.sum()) < 8:
                chosen = rng.choice(BASE_OBSERVATIONS, size=8, replace=False)
                keep = np.zeros(BASE_OBSERVATIONS, dtype=bool)
                keep[chosen] = True
            local = local[keep]
            z = z[keep]
            trait = trait[keep]

        position = centre + local
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

    observed_array = np.asarray(observed_rho, dtype=float)
    null_array = np.vstack(spatial_nulls)
    equal = species_equal_spatial_omnibus(observed_array, null_array)

    species_test = ttest_1samp(
        observed_array,
        popmean=0.0,
        alternative="greater",
    )
    common_beta, common_p = _fixed_effect_common_slope_test(model_groups)

    meta, slopes = random_effects_from_groups(model_groups)
    calibrated_meta = random_effects_slope_permutation_test(
        model_groups,
        n_permutations=META_PERMUTATIONS,
        seed=20260919,
        key=f"world|{seed}",
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

    slope_values = np.asarray([s.slope for s in slopes], dtype=float)
    return {
        "naive_p": float(naive.pvalue),
        "naive_rho": float(naive.statistic),
        "equal_matched_p": float(equal.p_upper),
        "equal_rho": float(equal.observed),
        "species_rho_ttest_p": float(species_test.pvalue),
        "common_slope_beta": float(common_beta),
        "common_slope_p": float(common_p),
        "meta_random_mean": float(meta.random_mean),
        "meta_mean_p": float(meta.p_mean_two_sided),
        "meta_heterogeneity_p": float(meta.p_heterogeneity),
        "meta_tau2": float(meta.tau2),
        "meta_omnibus_p": float(meta.p_omnibus),
        "calibrated_meta_mean_p": float(calibrated_meta.p_mean_permutation),
        "calibrated_meta_heterogeneity_p": float(
            calibrated_meta.p_heterogeneity_permutation
        ),
        "calibrated_meta_omnibus_p": float(
            calibrated_meta.p_omnibus_permutation
        ),
        "species_slope_sd": float(np.std(slope_values, ddof=1)),
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
                        "species_rho_ttest_detection": float(
                            np.mean(frame["species_rho_ttest_p"] < ALPHA)
                        ),
                        "common_slope_detection": float(
                            np.mean(frame["common_slope_p"] < ALPHA)
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
                        "calibrated_meta_mean_detection": float(
                            np.mean(frame["calibrated_meta_mean_p"] < ALPHA)
                        ),
                        "calibrated_meta_heterogeneity_detection": float(
                            np.mean(frame["calibrated_meta_heterogeneity_p"] < ALPHA)
                        ),
                        "calibrated_meta_omnibus_detection": float(
                            np.mean(frame["calibrated_meta_omnibus_p"] < ALPHA)
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                        "common_slope_beta_median": float(
                            frame["common_slope_beta"].median()
                        ),
                        "meta_random_mean_median": float(
                            frame["meta_random_mean"].median()
                        ),
                        "meta_tau2_median": float(frame["meta_tau2"].median()),
                        "species_slope_sd_median": float(
                            frame["species_slope_sd"].median()
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
            "max_common_slope_null_fpr": float(
                null["common_slope_detection"].max()
            ),
            "max_meta_mean_null_fpr": float(
                null["meta_mean_detection"].max()
            ),
            "max_meta_heterogeneity_null_fpr": float(
                null["meta_heterogeneity_detection"].max()
            ),
            "max_meta_omnibus_null_fpr": float(
                null["meta_omnibus_detection"].max()
            ),
            "max_calibrated_meta_mean_null_fpr": float(
                null["calibrated_meta_mean_detection"].max()
            ),
            "max_calibrated_meta_heterogeneity_null_fpr": float(
                null["calibrated_meta_heterogeneity_detection"].max()
            ),
            "max_calibrated_meta_omnibus_null_fpr": float(
                null["calibrated_meta_omnibus_detection"].max()
            ),
            "effect_0_4_reversal0_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.0
            ),
            "effect_0_4_reversal0_common_slope_detection_range": detection_range(
                weak, "common_slope_detection", 0.0
            ),
            "effect_0_4_reversal0_meta_omnibus_detection_range": detection_range(
                weak, "meta_omnibus_detection", 0.0
            ),
            "effect_0_4_reversal0_calibrated_meta_omnibus_detection_range": detection_range(
                weak, "calibrated_meta_omnibus_detection", 0.0
            ),
            "effect_0_4_reversal0_5_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.5
            ),
            "effect_0_4_reversal0_5_common_slope_detection_range": detection_range(
                weak, "common_slope_detection", 0.5
            ),
            "effect_0_4_reversal0_5_meta_heterogeneity_detection_range": detection_range(
                weak, "meta_heterogeneity_detection", 0.5
            ),
            "effect_0_4_reversal0_5_meta_omnibus_detection_range": detection_range(
                weak, "meta_omnibus_detection", 0.5
            ),
            "effect_0_4_reversal0_5_calibrated_meta_heterogeneity_detection_range": detection_range(
                weak, "calibrated_meta_heterogeneity_detection", 0.5
            ),
            "effect_0_4_reversal0_5_calibrated_meta_omnibus_detection_range": detection_range(
                weak, "calibrated_meta_omnibus_detection", 0.5
            ),
            "effect_0_8_reversal0_5_equal_detection_range": detection_range(
                strong, "equal_matched_detection", 0.5
            ),
            "effect_0_8_reversal0_5_meta_omnibus_detection_range": detection_range(
                strong, "meta_omnibus_detection", 0.5
            ),
            "effect_0_8_reversal0_5_calibrated_meta_omnibus_detection_range": detection_range(
                strong, "calibrated_meta_omnibus_detection", 0.5
            ),
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
