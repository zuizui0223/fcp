#!/usr/bin/env python3
"""Non-flower synthetic demonstration for disttrait.

The trait is an abstract binary morph. No flower-colour assumptions are used.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import spearmanr

from disttrait import distribution_spatial_association, gini_simpson, spatial_rho


def _one_hot_binary(values: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=int)
    return np.column_stack([1 - y, y]).astype(float)


def run_signal_demo(
    *,
    seed: int = 20260918,
    n_species: int = 80,
    observations_per_species: int = 60,
) -> dict[str, float]:
    """Recover an induced diversity-spatial gradient across synthetic species."""
    rng = np.random.default_rng(seed)
    diversities: list[float] = []
    spatial: list[float] = []
    latent_strength: list[float] = []

    for baseline in np.linspace(0.05, 0.48, int(n_species)):
        strength = 6.0 * float(baseline)
        intercept = float(np.log(baseline / (1.0 - baseline)))
        x = rng.uniform(-1.0, 1.0, int(observations_per_species))
        probability = 1.0 / (1.0 + np.exp(-(intercept + strength * x)))
        state = rng.binomial(1, probability)

        counts = np.bincount(state, minlength=2)
        diversities.append(gini_simpson(counts))
        spatial.append(
            spatial_rho(
                latitude=np.zeros_like(x),
                longitude=x,
                traits=_one_hot_binary(state),
            )
        )
        latent_strength.append(strength)

    association = distribution_spatial_association(diversities, spatial)
    strength_recovery = float(spearmanr(latent_strength, spatial).statistic)
    return {
        "n_species": int(n_species),
        "rho_diversity_spatial": float(association.observed),
        "rho_latent_strength_spatial": strength_recovery,
    }


def run_pooled_confounding_demo(
    *,
    seed: int = 123,
    n_species: int = 40,
    observations_per_species: int = 25,
    pooled_pairs: int = 50_000,
) -> dict[str, float]:
    """Show why species-conditioning matters under geographic species turnover.

    Species differ strongly in baseline morph frequency across a geographic
    gradient, but there is no within-species spatial dependence by construction.
    A naive pooled analysis therefore confounds between-species turnover with
    within-species organization.
    """
    rng = np.random.default_rng(seed)
    all_x: list[np.ndarray] = []
    all_state: list[np.ndarray] = []
    species_rhos: list[float] = []

    for centre in np.linspace(-5.0, 5.0, int(n_species)):
        x = centre + rng.normal(0.0, 0.15, int(observations_per_species))
        probability = 1.0 / (1.0 + np.exp(-1.2 * centre))
        state = rng.binomial(1, probability, int(observations_per_species))

        all_x.append(x)
        all_state.append(state)
        species_rhos.append(
            spatial_rho(
                latitude=np.zeros_like(x),
                longitude=x,
                traits=_one_hot_binary(state),
            )
        )

    x = np.concatenate(all_x)
    state = np.concatenate(all_state)
    n = len(x)
    left = rng.integers(0, n, int(pooled_pairs))
    right = rng.integers(0, n, int(pooled_pairs))
    keep = left != right
    left, right = left[keep], right[keep]
    pooled_distance = np.abs(x[left] - x[right])
    pooled_difference = (state[left] != state[right]).astype(float)
    naive_pooled_rho = float(spearmanr(pooled_distance, pooled_difference).statistic)

    finite = np.asarray([v for v in species_rhos if np.isfinite(v)], dtype=float)
    return {
        "n_species": int(n_species),
        "naive_pooled_rho": naive_pooled_rho,
        "species_conditioned_mean_rho": float(np.mean(finite)),
        "species_conditioned_median_rho": float(np.median(finite)),
    }


def main() -> int:
    out = {
        "signal_demo": run_signal_demo(),
        "pooled_confounding_demo": run_pooled_confounding_demo(),
    }
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
