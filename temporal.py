"""Temporal extension of the FCP phase-boundary model.

Key result for a single well-mixed population:

    dp/dt = p(1-p) [d(t) + b(1-2p)]

At the B boundary (p ~ 0), the rare A morph obeys

    d log p / dt = b + d(t).

Hence its long-run invasion exponent is

    lambda_A|B = b + mean_t[d(t)],

and similarly

    lambda_B|A = b - mean_t[d(t)].

Therefore purely additive temporal fluctuation in directional selection does
not generate an extra balancing term. Its variance cancels from the invasion
criterion; only the temporal mean matters. Temporal variability can promote
coexistence only when it interacts with additional structure/nonlinearity
(e.g. spatial structure, storage, stage structure, density dependence).
"""
from __future__ import annotations

from dataclasses import dataclass
from math import fsum
from typing import Iterable, Tuple


@dataclass(frozen=True)
class TemporalSummary:
    b: float
    mean_d: float


def time_mean(values: Iterable[float]) -> float:
    vals = tuple(float(x) for x in values)
    if not vals:
        raise ValueError("at least one temporal value is required")
    return fsum(vals) / len(vals)


def invasion_rates_from_series(b: float, d_series: Iterable[float]) -> Tuple[float, float]:
    """Return long-run rare-invasion exponents for additive temporal forcing."""
    mean_d = time_mean(d_series)
    return b + mean_d, b - mean_d


def classify_temporal_phase(b: float, d_series: Iterable[float], tol: float = 1e-12) -> str:
    a_into_b, b_into_a = invasion_rates_from_series(b, d_series)
    a_pos = a_into_b > tol
    b_pos = b_into_a > tol
    if a_pos and b_pos:
        return "protected_polymorphism"
    if a_pos and not b_pos:
        return "A_fixation"
    if not a_pos and b_pos:
        return "B_fixation"
    return "bistability"
