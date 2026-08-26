#!/usr/bin/env python3
"""Execution-safe wrapper for detection-aware latent C/S models.

The v1 model specification is retained unchanged. This wrapper only fixes the SciPy
optimizer call by binding keyword-only objective arguments through a closure; scipy.optimize.minimize
does not accept a `kwargs=` argument.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from scipy.optimize import minimize

import run_full_fcp_latent_models as base


def fit_model(
    obs: list[base.SpeciesObs],
    *,
    epsilon: float,
    beta_fixed_zero: bool,
    start: np.ndarray | None = None,
    weights: dict[str, float] | None = None,
) -> Any:
    k = 3 if beta_fixed_zero else 4
    if start is None or len(start) != k:
        start = np.zeros(k, dtype=float)
        start[-1] = 0.0

    def objective(theta: np.ndarray) -> float:
        return base.neg_log_posterior(
            theta,
            obs,
            epsilon=epsilon,
            beta_fixed_zero=beta_fixed_zero,
            weights=weights,
        )

    return minimize(
        objective,
        np.asarray(start, dtype=float),
        method="L-BFGS-B",
        bounds=[(-12, 12)] * k,
        options={"maxiter": 600, "ftol": 1e-10},
    )


base.fit_model = fit_model


if __name__ == "__main__":
    base.main()
