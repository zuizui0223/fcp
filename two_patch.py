"""Two-patch extension of the FCP phase-boundary model.

Model
-----
For morph-A frequencies p1 and p2 in two equal patches,

    dp1/dt = p1(1-p1) * [d + h + b(1-2p1)] + m(p2-p1)
    dp2/dt = p2(1-p2) * [d - h + b(1-2p2)] + m(p1-p2)

Parameters
----------
b : within-patch balancing strength (negative frequency dependence if b > 0)
d : net directional bias shared by both patches
h : spatial contrast; local directional biases are d+h and d-h
m : symmetric migration rate (m >= 0)

The dominant rare-invasion exponents are

    lambda_A|B = d + b - m + sqrt(h^2 + m^2)
    lambda_B|A = b - d - m + sqrt(h^2 + m^2)

Therefore protected polymorphism exists iff

    b + sqrt(h^2 + m^2) - m > |d|.

Define H_eff = sqrt(h^2 + m^2) - m. Then the boundary is

    b + H_eff > |d|.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Tuple


@dataclass(frozen=True)
class TwoPatchParameters:
    b: float
    d: float
    h: float
    m: float

    def __post_init__(self) -> None:
        if self.m < 0:
            raise ValueError("m must be non-negative")


def effective_heterogeneity(h: float, m: float) -> float:
    """Return the migration-attenuated contribution of spatial heterogeneity."""
    if m < 0:
        raise ValueError("m must be non-negative")
    return sqrt(h * h + m * m) - m


def invasion_rates(params: TwoPatchParameters) -> Tuple[float, float]:
    """Return (lambda_A_given_B, lambda_B_given_A)."""
    bonus = effective_heterogeneity(params.h, params.m)
    return params.d + params.b + bonus, params.b - params.d + bonus


def classify_phase(params: TwoPatchParameters, tol: float = 1e-12) -> str:
    """Classify the two boundary states from the two rare-invasion exponents."""
    a_into_b, b_into_a = invasion_rates(params)
    a_pos = a_into_b > tol
    b_pos = b_into_a > tol

    if a_pos and b_pos:
        return "protected_polymorphism"
    if a_pos and not b_pos:
        return "A_fixation"
    if not a_pos and b_pos:
        return "B_fixation"
    return "bistability"


def rhs(p1: float, p2: float, params: TwoPatchParameters) -> Tuple[float, float]:
    """Continuous-time frequency dynamics."""
    delta1 = params.d + params.h + params.b * (1.0 - 2.0 * p1)
    delta2 = params.d - params.h + params.b * (1.0 - 2.0 * p2)
    dp1 = p1 * (1.0 - p1) * delta1 + params.m * (p2 - p1)
    dp2 = p2 * (1.0 - p2) * delta2 + params.m * (p1 - p2)
    return dp1, dp2


def simulate(
    params: TwoPatchParameters,
    p1_0: float = 0.2,
    p2_0: float = 0.8,
    dt: float = 0.01,
    steps: int = 100_000,
) -> Tuple[float, float]:
    """Integrate by forward Euler and return the final patch frequencies.

    The integrator is intentionally dependency-free and is used for coarse
    consistency checks rather than precision bifurcation analysis.
    """
    if not (0.0 <= p1_0 <= 1.0 and 0.0 <= p2_0 <= 1.0):
        raise ValueError("initial frequencies must lie in [0, 1]")
    if dt <= 0 or steps <= 0:
        raise ValueError("dt and steps must be positive")

    p1, p2 = p1_0, p2_0
    for _ in range(steps):
        dp1, dp2 = rhs(p1, p2, params)
        p1 = min(1.0, max(0.0, p1 + dt * dp1))
        p2 = min(1.0, max(0.0, p2 + dt * dp2))
    return p1, p2
