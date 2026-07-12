"""Finite-population bridge from deterministic phase theory to observed polymorphism.

The deterministic invasion criterion answers whether a rare morph grows when
rare in an infinite population. Observed polymorphism additionally depends on:

1. origin: a new morph must arise (mutation/recombination/introgression),
2. establishment: it must escape early stochastic loss,
3. persistence: drift can still remove a protected morph in finite N,
4. observation: sampling can miss a morph that is truly present.

This module provides small dependency-free utilities for that bridge. It is
not intended to replace diffusion theory; it makes the observation pipeline
explicit and testable.
"""
from __future__ import annotations

from math import exp
from random import Random
from typing import Tuple


def branching_establishment_probability(r: float) -> float:
    """Approximate establishment probability 1-exp(-2r) for a rare beneficial type.

    Returns zero for non-positive invasion exponent. This is a weak-selection
    branching approximation and should be treated as a bridge, not an exact
    universal formula.
    """
    if r <= 0.0:
        return 0.0
    return 1.0 - exp(-2.0 * r)


def detection_probability(true_frequency: float, sample_size: int) -> float:
    """Probability of detecting at least one copy of a morph in a random sample."""
    if not 0.0 <= true_frequency <= 1.0:
        raise ValueError("true_frequency must lie in [0, 1]")
    if sample_size < 0:
        raise ValueError("sample_size must be non-negative")
    return 1.0 - (1.0 - true_frequency) ** sample_size


def polymorphism_detection_probability(p: float, sample_size: int) -> float:
    """Probability that both morphs are observed at least once in n samples."""
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must lie in [0, 1]")
    if sample_size < 0:
        raise ValueError("sample_size must be non-negative")
    if sample_size == 0:
        return 0.0
    return 1.0 - p ** sample_size - (1.0 - p) ** sample_size


def expected_origins(population_size: int, mutation_rate: float, generations: float) -> float:
    """Poisson mean for new copies arising over time in a haploid convention."""
    if population_size < 0 or mutation_rate < 0 or generations < 0:
        raise ValueError("arguments must be non-negative")
    return population_size * mutation_rate * generations


def probability_at_least_one_established_origin(
    population_size: int,
    mutation_rate: float,
    generations: float,
    invasion_exponent: float,
) -> float:
    """Poisson approximation for at least one origin that establishes."""
    mean_origins = expected_origins(population_size, mutation_rate, generations)
    p_est = branching_establishment_probability(invasion_exponent)
    return 1.0 - exp(-mean_origins * p_est)


def wright_fisher_step(p: float, population_size: int, selection_delta: float, rng: Random) -> float:
    """One haploid Wright-Fisher step with viability selection.

    A has relative fitness exp(selection_delta), B has fitness 1.
    """
    if population_size <= 0:
        raise ValueError("population_size must be positive")
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must lie in [0, 1]")
    w = exp(selection_delta)
    q = (p * w) / (p * w + (1.0 - p))
    count = sum(1 for _ in range(population_size) if rng.random() < q)
    return count / population_size


def observed_polymorphism_pipeline(
    p_true: float,
    sample_size: int,
) -> Tuple[float, float]:
    """Return probabilities of observed polymorphism and false monomorphism."""
    p_obs_poly = polymorphism_detection_probability(p_true, sample_size)
    return p_obs_poly, 1.0 - p_obs_poly
