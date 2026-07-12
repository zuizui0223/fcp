"""General spatiotemporal invasion framework for FCP.

For n patches, let x be the abundance vector of a rare morph and let

    dx/dt = A(t) x,

where A(t) combines local rare-morph growth and migration. The long-run
invasion exponent is the top Lyapunov exponent of the linear cocycle.
For a periodic environment this is the Floquet exponent

    gamma = (1/T) log rho(Phi(T)),

where Phi(T) is the period map and rho is its spectral radius.

For the two reciprocal morph invasions, define gamma_A and gamma_B. Then
protected polymorphism is equivalent to

    gamma_A > 0 and gamma_B > 0.

Compression:

    B_ST = (gamma_A + gamma_B)/2
    D_ST = (gamma_A - gamma_B)/2

so the phase condition becomes

    B_ST > |D_ST|.

This is exact by definition and generalizes the static spectral criterion.
It also makes clear that temporal variance alone is not a universal predictor:
ordering, autocorrelation, covariance with migration, and noncommuting A(t)
matrices can change the Floquet/Lyapunov exponent even when marginal temporal
means and variances are identical.
"""
from __future__ import annotations

from math import log, sqrt
from typing import Callable, List, Sequence, Tuple

Vector = List[float]
Matrix = List[List[float]]


def matvec(a: Sequence[Sequence[float]], x: Sequence[float]) -> Vector:
    return [sum(aij * xj for aij, xj in zip(row, x)) for row in a]


def add_scaled(x: Sequence[float], y: Sequence[float], scale: float) -> Vector:
    return [xi + scale * yi for xi, yi in zip(x, y)]


def norm2(x: Sequence[float]) -> float:
    return sqrt(sum(xi * xi for xi in x))


def rk4_step(x: Sequence[float], t: float, dt: float, a_of_t: Callable[[float], Matrix]) -> Vector:
    k1 = matvec(a_of_t(t), x)
    k2 = matvec(a_of_t(t + 0.5 * dt), add_scaled(x, k1, 0.5 * dt))
    k3 = matvec(a_of_t(t + 0.5 * dt), add_scaled(x, k2, 0.5 * dt))
    k4 = matvec(a_of_t(t + dt), add_scaled(x, k3, dt))
    return [
        xi + dt * (a + 2.0 * b + 2.0 * c + d) / 6.0
        for xi, a, b, c, d in zip(x, k1, k2, k3, k4)
    ]


def top_lyapunov_exponent(
    a_of_t: Callable[[float], Matrix],
    dimension: int,
    total_time: float,
    dt: float = 1e-3,
    renormalize_every: int = 50,
) -> float:
    """Estimate the top Lyapunov exponent of dx/dt=A(t)x.

    Uses RK4 with periodic renormalization. Suitable for deterministic or
    externally generated time-dependent environments.
    """
    if dimension <= 0 or total_time <= 0 or dt <= 0:
        raise ValueError("dimension, total_time, and dt must be positive")
    x = [1.0 / sqrt(dimension)] * dimension
    t = 0.0
    log_growth = 0.0
    steps = int(total_time / dt)
    for step in range(steps):
        x = rk4_step(x, t, dt, a_of_t)
        t += dt
        if (step + 1) % renormalize_every == 0:
            nrm = norm2(x)
            if nrm == 0.0:
                return float("-inf")
            log_growth += log(nrm)
            x = [xi / nrm for xi in x]
    nrm = norm2(x)
    if nrm == 0.0:
        return float("-inf")
    log_growth += log(nrm)
    return log_growth / t


def compress_invasion_exponents(gamma_a: float, gamma_b: float) -> Tuple[float, float]:
    """Return (effective balancing, effective directional asymmetry)."""
    return 0.5 * (gamma_a + gamma_b), 0.5 * (gamma_a - gamma_b)


def protected_polymorphism(gamma_a: float, gamma_b: float, tol: float = 1e-12) -> bool:
    return gamma_a > tol and gamma_b > tol
