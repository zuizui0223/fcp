"""Reusable statistical helpers for the frozen 34-species paper analysis.

Only production logic shared by the manuscript analyses belongs here. Exploratory
range/control models intentionally remain outside this module.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm

from .constants import MODEL_FORMULA

MIN_MODEL_SPECIES = 20
CI_METHOD = "Wald 95% CI from family-clustered sandwich standard errors"


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index, dtype=float)
    return (x - x.mean()) / sd


def prepare_model_data(data: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = data.copy()
    d["among"] = (d["spatial_scale"] == "among_population").astype(int)
    d["metric_z"] = zscore(d[metric])
    d["effort_z"] = zscore(np.log1p(pd.to_numeric(d["n_climate_cells"], errors="coerce")))
    return d.dropna(subset=["among", "metric_z", "effort_z", "family"]).copy()


def fit_model(data: pd.DataFrame, metric: str, clustered: bool = True):
    d = prepare_model_data(data, metric)
    if len(d) < MIN_MODEL_SPECIES or d["among"].nunique() < 2 or d["family"].nunique() < 2:
        return None, d
    X = sm.add_constant(d[["metric_z", "effort_z"]], has_constant="add")
    model = sm.GLM(d["among"], X, family=sm.families.Binomial())
    if clustered:
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": d["family"]})
    else:
        fit = model.fit()
    return fit, d


def fit_beta(data: pd.DataFrame, metric: str) -> float:
    fit, _ = fit_model(data, metric, clustered=False)
    return float(fit.params["metric_z"]) if fit is not None else np.nan


def analyse_metrics(
    data: pd.DataFrame,
    analysis_set: str,
    metrics: list[str],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[list[dict], list[dict]]:
    """Fit manuscript GLMs, label permutations and leave-one-family-out refits."""
    rows: list[dict] = []
    loo_rows: list[dict] = []
    for metric in metrics:
        fit, model_data = fit_model(data, metric, clustered=True)
        estimable = fit is not None
        observed = float(fit.params["metric_z"]) if estimable else np.nan
        se = float(fit.bse["metric_z"]) if estimable else np.nan

        permuted: list[float] = []
        if estimable:
            for _ in range(permutations):
                xp = model_data.copy()
                xp["spatial_scale"] = rng.permutation(xp["spatial_scale"].to_numpy())
                permuted.append(fit_beta(xp, metric))
        valid = np.asarray(permuted, dtype=float)
        valid = valid[np.isfinite(valid)]
        p_two = (
            float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))
            if estimable else np.nan
        )

        loo: list[float] = []
        for family in sorted(model_data["family"].dropna().astype(str).unique()):
            if not estimable:
                break
            beta = fit_beta(model_data.loc[model_data["family"].astype(str) != family], metric)
            loo.append(beta)
            loo_rows.append({
                "analysis_set": analysis_set,
                "metric": metric,
                "omitted_family": family,
                "estimate": beta,
                "odds_ratio": float(np.exp(beta)) if np.isfinite(beta) else np.nan,
            })
        valid_loo = np.asarray(loo, dtype=float)
        valid_loo = valid_loo[np.isfinite(valid_loo)]

        lo = observed - 1.96 * se if estimable else np.nan
        hi = observed + 1.96 * se if estimable else np.nan
        fit_history = getattr(fit, "fit_history", {}) if estimable else {}
        rows.append({
            "analysis_set": analysis_set,
            "metric": metric,
            "analysis_status": "complete" if estimable else "not_estimable",
            "analysis_reason": "" if estimable else "insufficient species, families, or response classes",
            "model_formula": MODEL_FORMULA,
            "estimator": "statsmodels GLM Binomial(logit)",
            "covariance": "family-clustered sandwich",
            "ci_method": CI_METHOD,
            "statsmodels_version": statsmodels.__version__,
            "n_species": int(len(model_data)),
            "n_families": int(model_data["family"].nunique()),
            "n_within": int((model_data["spatial_scale"] == "within_population").sum()),
            "n_among": int((model_data["spatial_scale"] == "among_population").sum()),
            "estimate": observed,
            "std_error_clustered": se,
            "estimate_ci_low": lo,
            "estimate_ci_high": hi,
            "odds_ratio": float(np.exp(observed)) if estimable else np.nan,
            "odds_ratio_ci_low": float(np.exp(lo)) if estimable else np.nan,
            "odds_ratio_ci_high": float(np.exp(hi)) if estimable else np.nan,
            "wald_p_value_clustered": float(fit.pvalues["metric_z"]) if estimable else np.nan,
            "permutation_p_two_sided": p_two,
            "permutations_requested": int(permutations),
            "permutations_valid": int(len(valid)),
            "converged": bool(getattr(fit, "converged", False)) if estimable else False,
            "iterations": int(fit_history.get("iteration", -1)) if estimable else -1,
            "predicted_probability_min": float(np.min(fit.fittedvalues)) if estimable else np.nan,
            "predicted_probability_max": float(np.max(fit.fittedvalues)) if estimable else np.nan,
            "loo_min_odds_ratio": float(np.exp(np.nanmin(valid_loo))) if len(valid_loo) else np.nan,
            "loo_max_odds_ratio": float(np.exp(np.nanmax(valid_loo))) if len(valid_loo) else np.nan,
            "loo_same_direction_fraction": (
                float(np.mean(np.sign(valid_loo) == np.sign(observed)))
                if len(valid_loo) and np.isfinite(observed) else np.nan
            ),
        })
    return rows, loo_rows
