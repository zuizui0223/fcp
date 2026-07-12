"""Minimal monomorphism–polymorphism phase model for floral colour polymorphism.

Log relative fitness advantage of morph A over B:
    Delta(p) = d + b * (1 - 2p)

Replicator dynamics:
    dp/dt = p * (1 - p) * Delta(p)

Protected polymorphism exists iff both boundary invasion rates are positive:
    r_A|B = d + b > 0
    r_B|A = b - d > 0

Therefore:
    b > |d|
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class PhaseResult:
    balancing_strength: float
    directional_bias: float
    invasion_A_into_B: float
    invasion_B_into_A: float
    phase: str
    equilibrium_p_A: float | None


def delta_log_fitness(p: float, balancing_strength: float, directional_bias: float) -> float:
    """Return log-fitness difference log(W_A) - log(W_B)."""
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must be in [0, 1]")
    if not all(isfinite(x) for x in (p, balancing_strength, directional_bias)):
        raise ValueError("all inputs must be finite")
    return directional_bias + balancing_strength * (1.0 - 2.0 * p)


def derivative(p: float, balancing_strength: float, directional_bias: float) -> float:
    """Replicator derivative dp/dt."""
    return p * (1.0 - p) * delta_log_fitness(p, balancing_strength, directional_bias)


def classify_phase(balancing_strength: float, directional_bias: float) -> PhaseResult:
    """Classify the evolutionary regime from the two boundary invasion rates."""
    if not all(isfinite(x) for x in (balancing_strength, directional_bias)):
        raise ValueError("parameters must be finite")

    r_a_into_b = directional_bias + balancing_strength
    r_b_into_a = balancing_strength - directional_bias

    if r_a_into_b > 0.0 and r_b_into_a > 0.0:
        phase = "protected_polymorphism"
        equilibrium = 0.5 * (1.0 + directional_bias / balancing_strength)
    elif r_a_into_b <= 0.0 and r_b_into_a > 0.0:
        phase = "B_fixation"
        equilibrium = 0.0
    elif r_a_into_b > 0.0 and r_b_into_a <= 0.0:
        phase = "A_fixation"
        equilibrium = 1.0
    else:
        phase = "bistability"
        equilibrium = None

    return PhaseResult(
        balancing_strength=balancing_strength,
        directional_bias=directional_bias,
        invasion_A_into_B=r_a_into_b,
        invasion_B_into_A=r_b_into_a,
        phase=phase,
        equilibrium_p_A=equilibrium,
    )


def simulate(
    p0: float,
    balancing_strength: float,
    directional_bias: float,
    dt: float = 0.01,
    steps: int = 10_000,
) -> list[float]:
    """Euler-integrate the replicator equation and return the frequency trajectory."""
    if not 0.0 <= p0 <= 1.0:
        raise ValueError("p0 must be in [0, 1]")
    if dt <= 0.0:
        raise ValueError("dt must be > 0")
    if steps < 1:
        raise ValueError("steps must be >= 1")

    p = p0
    trajectory = [p]
    for _ in range(steps):
        p += dt * derivative(p, balancing_strength, directional_bias)
        p = min(1.0, max(0.0, p))
        trajectory.append(p)
    return trajectory
