#!/usr/bin/env python3
"""Design-based power/precision simulation for the frozen 34-species JBI comparison.

The simulation preserves the observed 34-species predictor values, family labels,
class prevalence target, and two-predictor model structure. It asks how often a
sample of this exact size/design recovers negative climatic effects of specified
magnitudes. This is a precision/power diagnostic, not a post-hoc proof of adequacy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import brentq

METRICS = [
    "temperature_breadth",
    "moisture_breadth",
    "climatic_heterogeneity",
    "pca_dispersion",
    "pca_hull_area",
]
TARGET_ORS = [1.0, 0.8, 0.7, 0.6, 0.5, 0.4]


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = float(x.std(ddof=0))
    return (x - x.mean()) / sd


def prepare(data: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = data.copy()
    d["among"] = (d["spatial_scale"] == "among_population").astype(int)
    d["metric_z"] = zscore(d[metric])
    d["effort_z"] = zscore(np.log1p(pd.to_numeric(d["n_climate_cells"], errors="coerce")))
    return d.dropna(subset=["among", "metric_z", "effort_z", "family"]).copy()


def fit_glm(d: pd.DataFrame, y: np.ndarray, clustered: bool = True):
    X = sm.add_constant(d[["metric_z", "effort_z"]], has_constant="add")
    model = sm.GLM(y, X, family=sm.families.Binomial())
    if clustered:
        return model.fit(cov_type="cluster", cov_kwds={"groups": d["family"]})
    return model.fit()


def intercept_for_prevalence(xb: np.ndarray, target: float) -> float:
    def f(a: float) -> float:
        p = 1.0 / (1.0 + np.exp(-(a + xb)))
        return float(p.mean() - target)
    return float(brentq(f, -30.0, 30.0))


def simulate_design(d: pd.DataFrame, beta_metric: float, beta_effort: float,
                    target_prev: float, reps: int, rng: np.random.Generator) -> dict:
    xb = beta_metric * d["metric_z"].to_numpy() + beta_effort * d["effort_z"].to_numpy()
    intercept = intercept_for_prevalence(xb, target_prev)
    p = 1.0 / (1.0 + np.exp(-(intercept + xb)))
    estimates, ses, pvals, covers = [], [], [], []
    failures = 0
    for _ in range(reps):
        y = rng.binomial(1, p)
        if np.unique(y).size < 2:
            failures += 1
            continue
        try:
            fit = fit_glm(d, y, clustered=True)
            est = float(fit.params["metric_z"])
            se = float(fit.bse["metric_z"])
            pv = float(fit.pvalues["metric_z"])
            estimates.append(est)
            ses.append(se)
            pvals.append(pv)
            covers.append((est - 1.96 * se <= beta_metric) and (beta_metric <= est + 1.96 * se))
        except Exception:
            failures += 1
    e = np.asarray(estimates)
    s = np.asarray(ses)
    pv = np.asarray(pvals)
    return {
        "beta_true": beta_metric,
        "odds_ratio_true": float(np.exp(beta_metric)),
        "beta_effort_fixed": beta_effort,
        "target_among_prevalence": target_prev,
        "simulations_requested": reps,
        "simulations_valid": int(len(e)),
        "fit_failures": int(failures),
        "prob_estimate_negative": float(np.mean(e < 0)) if len(e) else np.nan,
        "prob_p_lt_0_05": float(np.mean(pv < 0.05)) if len(e) else np.nan,
        "prob_negative_and_p_lt_0_05": float(np.mean((e < 0) & (pv < 0.05))) if len(e) else np.nan,
        "median_estimated_or": float(np.exp(np.median(e))) if len(e) else np.nan,
        "median_ci_width_logodds": float(np.median(2 * 1.96 * s)) if len(e) else np.nan,
        "wald_95_coverage": float(np.mean(covers)) if len(e) else np.nan,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--reps", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=20260812)
    args = ap.parse_args()

    data = pd.read_csv(args.dataset)
    if len(data) != 34:
        raise SystemExit(f"Expected frozen 34-species dataset, found {len(data)}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    target_prev = float((data["spatial_scale"] == "among_population").mean())

    rows = []
    observed = []
    for metric in METRICS:
        d = prepare(data, metric)
        obs = fit_glm(d, d["among"].to_numpy(), clustered=False)
        beta_obs = float(obs.params["metric_z"])
        effort_obs = float(obs.params["effort_z"])
        observed.append({
            "metric": metric,
            "n_species": len(d),
            "n_families": int(d["family"].nunique()),
            "observed_beta_unclustered": beta_obs,
            "observed_or_unclustered": float(np.exp(beta_obs)),
            "observed_effort_beta_unclustered": effort_obs,
        })
        # Generic OR grid using the metric-specific observed effort coefficient.
        for target_or in TARGET_ORS:
            res = simulate_design(d, np.log(target_or), effort_obs, target_prev, args.reps, rng)
            rows.append({"metric": metric, "scenario": f"OR_{target_or:g}", **res})
        # Exact observed-effect scenario.
        res = simulate_design(d, beta_obs, effort_obs, target_prev, args.reps, rng)
        rows.append({"metric": metric, "scenario": "observed_effect", **res})

    table = pd.DataFrame(rows)
    obs_table = pd.DataFrame(observed)
    table.to_csv(outdir / "jbi_34species_power_precision_simulation.csv", index=False)
    obs_table.to_csv(outdir / "jbi_34species_power_precision_observed_design.csv", index=False)

    moisture = table[(table.metric == "moisture_breadth") & (table.scenario == "observed_effect")].iloc[0]
    manifest = {
        "status": "complete",
        "role": "design-based power/precision diagnostic, not a criterion for declaring the observed result valid",
        "n_species": 34,
        "n_families": int(data["family"].nunique()),
        "n_within": int((data["spatial_scale"] == "within_population").sum()),
        "n_among": int((data["spatial_scale"] == "among_population").sum()),
        "repetitions_per_scenario": args.reps,
        "seed": args.seed,
        "target_or_grid": TARGET_ORS,
        "simulation_design": (
            "Observed standardized predictor values and family labels are fixed; the effort coefficient is fixed "
            "at its observed unclustered estimate for each metric; the intercept is calibrated so expected among "
            "prevalence matches 14/34; binary responses are simulated and refit with the manuscript family-clustered GLM."
        ),
        "moisture_observed_effect_diagnostic": moisture.to_dict(),
        "interpretation_guard": (
            "Low simulated rejection probability indicates limited inferential precision at n=34; it does not invalidate "
            "the effect estimate. High simulated power under a large assumed effect does not prove the biological effect."
        ),
    }
    (outdir / "jbi_34species_power_precision_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
