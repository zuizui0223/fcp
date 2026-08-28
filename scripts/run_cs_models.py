#!/usr/bin/env python3
"""Primary mixed-preserving climate models for local coexistence (C) and spatial segregation (S).

C and S are fitted as separate positive documented-evidence outcomes. A species with
C=1,S=1 contributes a positive case to both models. Neither outcome is the complement
of the other.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

METRICS = [
    "temperature_breadth",
    "moisture_breadth",
    "climatic_heterogeneity",
    "pca_dispersion",
    "pca_hull_area",
]
OUTCOMES = {
    "C_local_coexistence_documented": "C",
    "S_spatial_segregation_documented": "S",
}
MIN_SPECIES = 20


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index, dtype=float)
    return (x - x.mean()) / sd


def prepare(data: pd.DataFrame, outcome: str, metric: str) -> pd.DataFrame:
    d = data.sort_values("canonical_name", kind="stable").reset_index(drop=True).copy()
    d["outcome"] = pd.to_numeric(d[outcome], errors="coerce")
    d["metric_z"] = zscore(d[metric])
    d["effort_z"] = zscore(np.log1p(pd.to_numeric(d["n_climate_cells"], errors="coerce")))
    d["source_effort_z"] = zscore(np.log1p(pd.to_numeric(d.get("n_resolved_sources", 0), errors="coerce")))
    return d.dropna(subset=["outcome", "metric_z", "effort_z", "family"]).copy()


def fit_glm(d: pd.DataFrame, clustered: bool = True, source_sensitivity: bool = False):
    predictors = ["metric_z", "effort_z"] + (["source_effort_z"] if source_sensitivity else [])
    required = ["outcome", "family", *predictors]
    x = d.dropna(subset=required).copy()
    if len(x) < MIN_SPECIES or x.outcome.nunique() < 2 or x.family.nunique() < 2:
        return None, x
    X = sm.add_constant(x[predictors], has_constant="add")
    model = sm.GLM(x.outcome, X, family=sm.families.Binomial())
    fit = model.fit(
        cov_type="cluster",
        cov_kwds={"groups": x.family},
    ) if clustered else model.fit()
    return fit, x


def unclustered_beta(d: pd.DataFrame) -> float:
    fit, _ = fit_glm(d, clustered=False, source_sensitivity=False)
    if fit is None:
        return np.nan
    return float(fit.params["metric_z"])


def fit_one(data: pd.DataFrame, outcome: str, metric: str, permutations: int, rng: np.random.Generator):
    d = prepare(data, outcome, metric)
    fit, model_data = fit_glm(d, clustered=True, source_sensitivity=False)
    if fit is None:
        return {
            "outcome": outcome,
            "metric": metric,
            "analysis_status": "not_estimable",
            "n_species": int(len(model_data)),
            "n_families": int(model_data.family.nunique()),
            "n_positive": int(model_data.outcome.sum()) if len(model_data) else 0,
            "n_negative": int(len(model_data) - model_data.outcome.sum()) if len(model_data) else 0,
        }, []

    beta = float(fit.params["metric_z"])
    se = float(fit.bse["metric_z"])
    lo = beta - 1.96 * se
    hi = beta + 1.96 * se

    perm = []
    labels = model_data.outcome.to_numpy().copy()
    for _ in range(permutations):
        xp = model_data.copy()
        xp["outcome"] = rng.permutation(labels)
        b = unclustered_beta(xp)
        if np.isfinite(b):
            perm.append(b)
    perm = np.asarray(perm, dtype=float)
    p_perm = float((1 + np.sum(np.abs(perm) >= abs(beta))) / (1 + len(perm))) if len(perm) else np.nan

    loo_rows = []
    loo_betas = []
    for family in sorted(model_data.family.dropna().astype(str).unique()):
        sub = model_data.loc[model_data.family.astype(str) != family].copy()
        b = unclustered_beta(sub)
        loo_betas.append(b)
        loo_rows.append({
            "outcome": outcome,
            "metric": metric,
            "omitted_family": family,
            "estimate": b,
            "odds_ratio": float(np.exp(b)) if np.isfinite(b) else np.nan,
        })
    loo = np.asarray([x for x in loo_betas if np.isfinite(x)], dtype=float)

    source_fit, source_data = fit_glm(d, clustered=True, source_sensitivity=True)
    source_beta = float(source_fit.params["metric_z"]) if source_fit is not None else np.nan
    source_se = float(source_fit.bse["metric_z"]) if source_fit is not None else np.nan

    row = {
        "outcome": outcome,
        "outcome_short": OUTCOMES[outcome],
        "metric": metric,
        "analysis_status": "complete",
        "model_formula": f"{OUTCOMES[outcome]} ~ metric_z + effort_z",
        "estimator": "statsmodels GLM Binomial(logit)",
        "covariance": "family-clustered sandwich",
        "statsmodels_version": statsmodels.__version__,
        "n_species": int(len(model_data)),
        "n_families": int(model_data.family.nunique()),
        "n_positive": int(model_data.outcome.sum()),
        "n_negative": int(len(model_data) - model_data.outcome.sum()),
        "estimate": beta,
        "std_error_clustered": se,
        "estimate_ci_low": lo,
        "estimate_ci_high": hi,
        "odds_ratio": float(np.exp(beta)),
        "odds_ratio_ci_low": float(np.exp(lo)),
        "odds_ratio_ci_high": float(np.exp(hi)),
        "wald_p_value_clustered": float(fit.pvalues["metric_z"]),
        "permutation_p_two_sided": p_perm,
        "permutations_requested": int(permutations),
        "permutations_valid": int(len(perm)),
        "loo_min_odds_ratio": float(np.exp(np.min(loo))) if len(loo) else np.nan,
        "loo_max_odds_ratio": float(np.exp(np.max(loo))) if len(loo) else np.nan,
        "loo_same_direction_fraction": float(np.mean(np.sign(loo) == np.sign(beta))) if len(loo) else np.nan,
        "source_effort_sensitivity_n": int(len(source_data)),
        "source_effort_sensitivity_estimate": source_beta,
        "source_effort_sensitivity_se": source_se,
        "source_effort_sensitivity_odds_ratio": float(np.exp(source_beta)) if np.isfinite(source_beta) else np.nan,
        "semantic_guard": "Outcome is documented positive evidence; zero does not imply biological absence.",
    }
    return row, loo_rows


def add_holm(results: pd.DataFrame) -> pd.DataFrame:
    out = results.copy()
    out["wald_p_holm_within_outcome"] = np.nan
    out["permutation_p_holm_within_outcome"] = np.nan
    for outcome, idx in out.groupby("outcome").groups.items():
        idx = list(idx)
        wald = pd.to_numeric(out.loc[idx, "wald_p_value_clustered"], errors="coerce")
        perm = pd.to_numeric(out.loc[idx, "permutation_p_two_sided"], errors="coerce")
        wmask = wald.notna()
        pmask = perm.notna()
        if wmask.any():
            out.loc[np.asarray(idx)[wmask.to_numpy()], "wald_p_holm_within_outcome"] = multipletests(wald[wmask], method="holm")[1]
        if pmask.any():
            out.loc[np.asarray(idx)[pmask.to_numpy()], "permutation_p_holm_within_outcome"] = multipletests(perm[pmask], method="holm")[1]
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260826)
    args = p.parse_args()

    data = pd.read_csv(args.dataset)
    required = {"canonical_name", "family", "n_climate_cells", *METRICS, *OUTCOMES}
    missing = required - set(data.columns)
    if missing:
        raise SystemExit(f"Dataset missing required columns: {sorted(missing)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    rows = []
    loo = []
    for outcome in OUTCOMES:
        for metric in METRICS:
            row, family_rows = fit_one(data, outcome, metric, args.permutations, rng)
            rows.append(row)
            loo.extend(family_rows)

    results = add_holm(pd.DataFrame(rows))
    results.to_csv(outdir / "cs_five_metric_models.csv", index=False)
    pd.DataFrame(loo).to_csv(outdir / "cs_leave_one_family_out.csv", index=False)

    summary = {
        "dataset": str(args.dataset),
        "n_species": int(len(data)),
        "n_families": int(data.family.nunique()),
        "C_positive": int(pd.to_numeric(data.C_local_coexistence_documented).sum()),
        "S_positive": int(pd.to_numeric(data.S_spatial_segregation_documented).sum()),
        "C_and_S_positive": int(((data.C_local_coexistence_documented == 1) & (data.S_spatial_segregation_documented == 1)).sum()),
        "metrics": METRICS,
        "outcomes": list(OUTCOMES),
        "permutations": int(args.permutations),
        "seed": int(args.seed),
        "primary_formulae": ["C ~ metric_z + effort_z", "S ~ metric_z + effort_z"],
        "source_effort_sensitivity": "adds z(log1p(n_resolved_sources)) when available",
        "semantic_guard": "C and S are separate positive documented-evidence outcomes; neither is the complement of the other.",
    }
    (outdir / "cs_model_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))

    if len(results) != 10:
        raise SystemExit(f"Expected 10 primary rows, found {len(results)}")
    if not set(results.analysis_status).issubset({"complete", "not_estimable"}):
        raise SystemExit("Unexpected analysis status")


if __name__ == "__main__":
    main()
