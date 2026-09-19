#!/usr/bin/env python3
"""MNAR observation-process stress test for species-conditioned inference.

The latent biological null has no within-species spatial trait dependence.
Selection is then applied before analysis under four observation mechanisms:

1. MCAR;
2. trait-only selection;
3. position-only selection;
4. joint trait-by-position selection.

The first three preserve trait-position independence within species under the
latent null. The joint mechanism can induce trait-position association in the
observed sample even when the latent biological process is null.

This benchmark therefore distinguishes conditional calibration on the observed
sample from protection against outcome-dependent observation processes.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit
from scipy.stats import chi2, spearmanr

from disttrait import spatial_permutation_null, species_equal_spatial_omnibus


MECHANISMS = ("mcar", "trait_only", "position_only", "joint_trait_position")
STRENGTHS = (0.8, 1.6)
WORLDS_PER_CELL = 40
N_SPECIES = 20
LATENT_OBSERVATIONS = 80
N_PERMUTATIONS = 39
POOLED_PAIR_SAMPLE = 8000
ALPHA = 0.05
SELECTION_INTERCEPT = -0.35


def _one_hot_binary(values: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=int)
    return np.column_stack([1 - y, y]).astype(float)


def _cell_seed(mechanism: str, strength: float, world: int) -> int:
    code = {
        "mcar": 1,
        "trait_only": 2,
        "position_only": 3,
        "joint_trait_position": 4,
    }[mechanism]
    return 10_000_000 + code * 1_000_000 + int(strength * 1000) * 100 + int(world)


def _selection_probability(
    state: np.ndarray,
    z: np.ndarray,
    *,
    mechanism: str,
    strength: float,
) -> np.ndarray:
    signed_state = 2.0 * np.asarray(state, dtype=float) - 1.0
    z = np.asarray(z, dtype=float)
    if mechanism == "mcar":
        score = np.zeros_like(z)
    elif mechanism == "trait_only":
        score = float(strength) * signed_state
    elif mechanism == "position_only":
        score = float(strength) * z
    elif mechanism == "joint_trait_position":
        score = float(strength) * signed_state * z
    else:
        raise ValueError(f"unknown mechanism: {mechanism}")
    return expit(SELECTION_INTERCEPT + score)


def _profile_intercept(y: np.ndarray, z: np.ndarray, beta: float, bound: float = 20.0) -> float:
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


def _profile_loglikelihood(groups: list[tuple[np.ndarray, np.ndarray]], beta: float) -> float:
    total = 0.0
    for y, z in groups:
        alpha = _profile_intercept(y, z, beta)
        eta = alpha + float(beta) * z
        total += float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    return total


def _fixed_effect_logit_test(groups: list[tuple[np.ndarray, np.ndarray]]) -> tuple[float, float]:
    ll0 = _profile_loglikelihood(groups, 0.0)
    fit = minimize_scalar(
        lambda beta: -_profile_loglikelihood(groups, float(beta)),
        bounds=(0.0, 8.0),
        method="bounded",
        options={"xatol": 1e-6},
    )
    beta_hat = float(fit.x)
    ll1 = -float(fit.fun)
    lrt = max(0.0, 2.0 * (ll1 - ll0))
    p = 1.0 if lrt <= 1e-12 else float(0.5 * chi2.sf(lrt, 1))
    return beta_hat, p


def run_world(
    *,
    seed: int,
    mechanism: str,
    strength: float,
) -> dict[str, float]:
    rng = np.random.default_rng(int(seed))
    centres = np.linspace(-4.0, 4.0, N_SPECIES)

    observed_rho: list[float] = []
    nulls: list[np.ndarray] = []
    model_groups: list[tuple[np.ndarray, np.ndarray]] = []
    all_position: list[np.ndarray] = []
    all_state: list[np.ndarray] = []
    retained: list[int] = []
    within_species_state_position_rho: list[float] = []

    for i, centre in enumerate(centres):
        z = rng.normal(0.0, 1.0, LATENT_OBSERVATIONS)
        local = 0.25 * z
        position = centre + local

        # Latent biological null: state probability changes among species but
        # is constant across local positions within every species.
        baseline = float(expit(centre))
        state = rng.binomial(1, baseline, LATENT_OBSERVATIONS)

        keep_probability = _selection_probability(
            state,
            z,
            mechanism=mechanism,
            strength=strength,
        )
        keep = rng.random(LATENT_OBSERVATIONS) < keep_probability

        # This minimum is a data-availability gate, not outcome-adaptive rescue:
        # worlds/species failing it are underidentified and excluded from the
        # cell rather than topped up with favorable rows.
        if int(keep.sum()) < 8:
            return {
                "evaluable": 0.0,
                "equal_p": float("nan"),
                "equal_rho": float("nan"),
                "fixed_logit_p": float("nan"),
                "fixed_logit_beta": float("nan"),
                "naive_p": float("nan"),
                "naive_rho": float("nan"),
                "retained_mean": float("nan"),
                "retained_min": float("nan"),
                "observed_state_position_abs_rho": float("nan"),
            }

        z_obs = z[keep]
        position_obs = position[keep]
        state_obs = state[keep]

        rho, null = spatial_permutation_null(
            latitude=np.zeros_like(position_obs),
            longitude=position_obs,
            traits=_one_hot_binary(state_obs),
            n_permutations=N_PERMUTATIONS,
            seed=20260919,
            key=f"sp{i}|{seed}|{mechanism}|{strength}",
        )
        observed_rho.append(float(rho))
        nulls.append(null)
        model_groups.append((state_obs.astype(float), z_obs.astype(float)))
        all_position.append(position_obs)
        all_state.append(state_obs)
        retained.append(int(len(state_obs)))

        if np.ptp(state_obs) <= 0:
            within_species_state_position_rho.append(0.0)
        else:
            within_species_state_position_rho.append(
                abs(float(spearmanr(z_obs, state_obs).statistic))
            )

    equal = species_equal_spatial_omnibus(
        np.asarray(observed_rho, dtype=float),
        np.vstack(nulls),
    )
    beta_hat, fixed_p = _fixed_effect_logit_test(model_groups)

    x = np.concatenate(all_position)
    y = np.concatenate(all_state)
    n = len(x)
    left = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    right = rng.integers(0, n, POOLED_PAIR_SAMPLE)
    pair_keep = left != right
    left = left[pair_keep]
    right = right[pair_keep]
    naive = spearmanr(
        np.abs(x[left] - x[right]),
        (y[left] != y[right]).astype(float),
    )

    return {
        "evaluable": 1.0,
        "equal_p": float(equal.p_upper),
        "equal_rho": float(equal.observed),
        "fixed_logit_p": float(fixed_p),
        "fixed_logit_beta": float(beta_hat),
        "naive_p": float(naive.pvalue),
        "naive_rho": float(naive.statistic),
        "retained_mean": float(np.mean(retained)),
        "retained_min": float(np.min(retained)),
        "observed_state_position_abs_rho": float(
            np.mean(within_species_state_position_rho)
        ),
    }


def run_benchmark(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []

    for mechanism in MECHANISMS:
        for strength in STRENGTHS:
            rows = [
                run_world(
                    seed=_cell_seed(mechanism, strength, world),
                    mechanism=mechanism,
                    strength=strength,
                )
                for world in range(int(worlds_per_cell))
            ]
            frame = pd.DataFrame(rows)
            evaluable = frame[frame["evaluable"] == 1.0].copy()
            cells.append(
                {
                    "mechanism": mechanism,
                    "strength": float(strength),
                    "worlds": int(worlds_per_cell),
                    "evaluable_worlds": int(len(evaluable)),
                    "equal_rejection_fraction": float(
                        np.mean(evaluable["equal_p"] < ALPHA)
                    ) if len(evaluable) else float("nan"),
                    "fixed_logit_rejection_fraction": float(
                        np.mean(evaluable["fixed_logit_p"] < ALPHA)
                    ) if len(evaluable) else float("nan"),
                    "naive_rejection_fraction": float(
                        np.mean(evaluable["naive_p"] < ALPHA)
                    ) if len(evaluable) else float("nan"),
                    "equal_rho_median": float(evaluable["equal_rho"].median())
                    if len(evaluable) else float("nan"),
                    "fixed_logit_beta_median": float(
                        evaluable["fixed_logit_beta"].median()
                    ) if len(evaluable) else float("nan"),
                    "retained_mean": float(evaluable["retained_mean"].mean())
                    if len(evaluable) else float("nan"),
                    "minimum_species_rows": int(evaluable["retained_min"].min())
                    if len(evaluable) else 0,
                    "observed_state_position_abs_rho_mean": float(
                        evaluable["observed_state_position_abs_rho"].mean()
                    ) if len(evaluable) else float("nan"),
                }
            )

    table = pd.DataFrame(cells)
    safe = table[table["mechanism"].isin(["mcar", "trait_only", "position_only"])]
    joint = table[table["mechanism"] == "joint_trait_position"]

    return {
        "cells": cells,
        "summary": {
            "minimum_evaluable_worlds": int(table["evaluable_worlds"].min()),
            "max_equal_rejection_nonjoint": float(
                safe["equal_rejection_fraction"].max()
            ),
            "max_fixed_logit_rejection_nonjoint": float(
                safe["fixed_logit_rejection_fraction"].max()
            ),
            "joint_equal_rejection_range": [
                float(joint["equal_rejection_fraction"].min()),
                float(joint["equal_rejection_fraction"].max()),
            ],
            "joint_fixed_logit_rejection_range": [
                float(joint["fixed_logit_rejection_fraction"].min()),
                float(joint["fixed_logit_rejection_fraction"].max()),
            ],
            "joint_observed_abs_state_position_rho_range": [
                float(joint["observed_state_position_abs_rho_mean"].min()),
                float(joint["observed_state_position_abs_rho_mean"].max()),
            ],
        },
    }


def main() -> int:
    print(json.dumps(run_benchmark(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
