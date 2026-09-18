#!/usr/bin/env python3
"""Hierarchical random-slope comparator under direction heterogeneity.

Extends the v0.7 benchmark with an observation-level Gaussian mixed model:

- fixed common within-species slope;
- species random intercept;
- optional species random slope and intercept-slope covariance.

The random-slope comparator uses an ML likelihood-ratio statistic between
random-intercept-only and random-intercept-plus-slope models. The reference
p-value uses a chi-square distribution with 2 degrees of freedom for the two
added covariance parameters. Because variance-component boundaries make this
an approximation rather than an exact finite-sample test, the benchmark
explicitly evaluates its null rejection rate.
"""
from __future__ import annotations

import json
import math
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr, t as student_t
import statsmodels.formula.api as smf

from disttrait import (
    continuous_spatial_permutation_null,
    species_equal_spatial_omnibus,
)


EFFECT_SIZES = (0.0, 0.4, 0.8, 1.2)
REVERSAL_FRACTIONS = (0.0, 0.25, 0.5)
MISSING_FRACTIONS = (0.0, 0.5)
WORLDS_PER_CELL = 40
N_SPECIES = 20
BASE_OBSERVATIONS = 20
N_PERMUTATIONS = 39
POOLED_PAIR_SAMPLE = 8_000
ALPHA = 0.05


def _cell_seed(
    effect_size: float,
    reversal_fraction: float,
    missing_fraction: float,
    world: int,
) -> int:
    # Deliberately identical to v0.7 so all comparators see the same worlds.
    return (
        7_000_000
        + int(effect_size * 1000) * 10_000
        + int(reversal_fraction * 100) * 100
        + int(missing_fraction * 100)
        + int(world)
    )


def _fixed_effect_common_slope_test(
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
    beta = float((zc @ yc) / denom)
    residual = yc - beta * zc
    df = int(n - len(groups) - 1)
    mse = float((residual @ residual) / df)
    if mse <= 0:
        return beta, 0.0 if beta > 0 else 1.0
    se = math.sqrt(mse / denom)
    return beta, float(student_t.sf(beta / se, df))


def _fit_mixed(
    frame: pd.DataFrame,
    *,
    re_formula: str,
):
    last_error: Exception | None = None
    for method in ("lbfgs", "powell"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = smf.mixedlm(
                    "trait ~ z",
                    frame,
                    groups=frame["species"],
                    re_formula=re_formula,
                ).fit(
                    reml=False,
                    method=method,
                    disp=False,
                    maxiter=500,
                )
            if np.isfinite(float(fit.llf)):
                return fit, method
        except Exception as exc:
            last_error = exc
    if last_error is not None:
        raise RuntimeError("mixed-model fit failed") from last_error
    raise RuntimeError("mixed-model fit failed without exception")


def _random_slope_lrt(
    groups: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[float, float, float, float, bool, str, str]:
    rows: list[tuple[str, float, float]] = []
    for i, (z, y) in enumerate(groups):
        for zz, yy in zip(np.asarray(z, float), np.asarray(y, float), strict=True):
            rows.append((f"sp{i}", float(zz), float(yy)))
    frame = pd.DataFrame(rows, columns=["species", "z", "trait"])

    null_fit, null_method = _fit_mixed(frame, re_formula="1")
    full_fit, full_method = _fit_mixed(frame, re_formula="~z")

    lrt = max(0.0, 2.0 * (float(full_fit.llf) - float(null_fit.llf)))
    # Approximate: full model adds random-slope variance and intercept-slope
    # covariance. Boundary behavior is assessed empirically in this benchmark.
    p = float(chi2.sf(lrt, 2))
    fixed_beta = float(full_fit.fe_params["z"])
    slope_variance = float(full_fit.cov_re.loc["z", "z"])
    converged = bool(null_fit.converged and full_fit.converged)
    return (
        lrt,
        p,
        fixed_beta,
        slope_variance,
        converged,
        str(null_method),
        str(full_method),
    )


def run_world(
    *,
    seed: int,
    effect_size: float,
    reversal_fraction: float,
    missing_fraction: float,
) -> dict[str, float | bool | str]:
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
    common_beta, common_p = _fixed_effect_common_slope_test(model_groups)
    (
        mixed_lrt,
        mixed_p,
        mixed_fixed_beta,
        mixed_slope_variance,
        mixed_converged,
        null_method,
        full_method,
    ) = _random_slope_lrt(model_groups)

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
        "common_slope_beta": float(common_beta),
        "common_slope_p": float(common_p),
        "mixed_lrt": float(mixed_lrt),
        "mixed_random_slope_p": float(mixed_p),
        "mixed_fixed_beta": float(mixed_fixed_beta),
        "mixed_slope_variance": float(mixed_slope_variance),
        "mixed_converged": bool(mixed_converged),
        "mixed_null_method": null_method,
        "mixed_full_method": full_method,
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []
    method_counts: dict[str, int] = {}

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
                for col in ("mixed_null_method", "mixed_full_method"):
                    for name, n in frame[col].value_counts().items():
                        key = f"{col}:{name}"
                        method_counts[key] = method_counts.get(key, 0) + int(n)

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
                        "common_slope_detection": float(
                            np.mean(frame["common_slope_p"] < ALPHA)
                        ),
                        "mixed_random_slope_detection": float(
                            np.mean(frame["mixed_random_slope_p"] < ALPHA)
                        ),
                        "mixed_convergence_fraction": float(
                            np.mean(frame["mixed_converged"].astype(bool))
                        ),
                        "equal_rho_median": float(frame["equal_rho"].median()),
                        "common_slope_beta_median": float(
                            frame["common_slope_beta"].median()
                        ),
                        "mixed_fixed_beta_median": float(
                            frame["mixed_fixed_beta"].median()
                        ),
                        "mixed_slope_variance_median": float(
                            frame["mixed_slope_variance"].median()
                        ),
                        "mixed_lrt_median": float(frame["mixed_lrt"].median()),
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
        "fit_method_counts": method_counts,
        "summary": {
            "max_equal_matched_null_fpr": float(
                null["equal_matched_detection"].max()
            ),
            "max_common_slope_null_fpr": float(
                null["common_slope_detection"].max()
            ),
            "max_mixed_random_slope_null_fpr": float(
                null["mixed_random_slope_detection"].max()
            ),
            "min_mixed_convergence_fraction": float(
                table["mixed_convergence_fraction"].min()
            ),
            "effect_0_4_reversal0_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.0
            ),
            "effect_0_4_reversal0_common_slope_detection_range": detection_range(
                weak, "common_slope_detection", 0.0
            ),
            "effect_0_4_reversal0_mixed_detection_range": detection_range(
                weak, "mixed_random_slope_detection", 0.0
            ),
            "effect_0_4_reversal0_5_equal_detection_range": detection_range(
                weak, "equal_matched_detection", 0.5
            ),
            "effect_0_4_reversal0_5_common_slope_detection_range": detection_range(
                weak, "common_slope_detection", 0.5
            ),
            "effect_0_4_reversal0_5_mixed_detection_range": detection_range(
                weak, "mixed_random_slope_detection", 0.5
            ),
            "effect_0_8_reversal0_5_equal_detection_range": detection_range(
                strong, "equal_matched_detection", 0.5
            ),
            "effect_0_8_reversal0_5_mixed_detection_range": detection_range(
                strong, "mixed_random_slope_detection", 0.5
            ),
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
