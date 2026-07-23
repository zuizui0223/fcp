#!/usr/bin/env python3
"""Robustness checks for within- versus among-population climatic niche contrasts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm

METRICS = ["pca_hull_area", "moisture_breadth"]
MIN_MODEL_SPECIES = 20
MODEL_FORMULA = "among ~ metric_z + effort_z"
CI_METHOD = "Wald 95% CI from family-clustered sandwich standard errors"


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def prepare_model_data(d: pd.DataFrame, metric: str) -> pd.DataFrame:
    x = d.copy()
    x["among"] = (x["spatial_scale"] == "among_population").astype(int)
    x["metric_z"] = zscore(x[metric])
    x["effort_z"] = zscore(np.log1p(pd.to_numeric(x["n_climate_cells"], errors="coerce")))
    return x.dropna(subset=["among", "metric_z", "effort_z", "family"])


def fit_model(d: pd.DataFrame, metric: str, clustered: bool = True):
    x = prepare_model_data(d, metric)
    if (
        len(x) < MIN_MODEL_SPECIES
        or x["among"].nunique() < 2
        or x["family"].nunique() < 2
    ):
        return None, x
    model = sm.GLM(
        x["among"],
        sm.add_constant(x[["metric_z", "effort_z"]], has_constant="add"),
        family=sm.families.Binomial(),
    )
    if clustered:
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": x["family"]})
    else:
        fit = model.fit()
    return fit, x


def fit_beta(d: pd.DataFrame, metric: str) -> float:
    fit, _ = fit_model(d, metric, clustered=False)
    return float(fit.params["metric_z"]) if fit is not None else np.nan


def analyse_set(
    data: pd.DataFrame,
    analysis_set: str,
    metrics: list[str],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    loo_rows: list[dict] = []

    for metric in metrics:
        fit, model_data = fit_model(data, metric, clustered=True)
        estimable = fit is not None
        observed = float(fit.params["metric_z"]) if estimable else np.nan
        se = float(fit.bse["metric_z"]) if estimable else np.nan

        if estimable:
            permuted = []
            for _ in range(permutations):
                xp = model_data.copy()
                xp["spatial_scale"] = rng.permutation(xp["spatial_scale"].to_numpy())
                permuted.append(fit_beta(xp, metric))
            valid = np.asarray(permuted, dtype=float)
            valid = valid[np.isfinite(valid)]
        else:
            valid = np.asarray([], dtype=float)

        p_two = (
            float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))
            if np.isfinite(observed)
            else np.nan
        )

        families = sorted(model_data["family"].dropna().astype(str).unique())
        loo = []
        if estimable:
            for fam in families:
                beta = fit_beta(model_data.loc[model_data["family"].astype(str) != fam], metric)
                loo.append(beta)
                loo_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "metric": metric,
                        "omitted_family": fam,
                        "estimate": beta,
                        "odds_ratio": float(np.exp(beta)) if np.isfinite(beta) else np.nan,
                    }
                )
        valid_loo = np.asarray(loo, dtype=float)
        valid_loo = valid_loo[np.isfinite(valid_loo)]

        ci_low_beta = observed - 1.96 * se if estimable else np.nan
        ci_high_beta = observed + 1.96 * se if estimable else np.nan
        fit_history = getattr(fit, "fit_history", {}) if estimable else {}
        rows.append(
            {
                "analysis_set": analysis_set,
                "metric": metric,
                "analysis_status": "complete" if estimable else "not_estimable",
                "analysis_reason": "" if estimable else (
                    f"fewer than {MIN_MODEL_SPECIES} species, fewer than two families, "
                    "or only one spatial-scale class"
                ),
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
                "estimate_ci_low": ci_low_beta,
                "estimate_ci_high": ci_high_beta,
                "odds_ratio": float(np.exp(observed)) if np.isfinite(observed) else np.nan,
                "odds_ratio_ci_low": float(np.exp(ci_low_beta)) if np.isfinite(ci_low_beta) else np.nan,
                "odds_ratio_ci_high": float(np.exp(ci_high_beta)) if np.isfinite(ci_high_beta) else np.nan,
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
                    if len(valid_loo) and np.isfinite(observed)
                    else np.nan
                ),
            }
        )
    return rows, loo_rows


def write_frozen_manifest(data: pd.DataFrame, outdir: Path) -> dict:
    strict = data.loc[data["classification_source"] == "baseline_unambiguous"].copy()
    columns = [
        "canonical_name",
        "family",
        "spatial_scale",
        "classification_source",
        "n_climate_cells",
    ]
    optional = [c for c in ("evidence_source", "source_id", "decision_note") if c in strict.columns]
    columns.extend(optional)
    strict = strict[columns].sort_values("canonical_name", kind="stable").reset_index(drop=True)
    csv_bytes = strict.to_csv(index=False).encode("utf-8")
    digest = hashlib.sha256(csv_bytes).hexdigest()
    strict.to_csv(outdir / "baseline_unambiguous_frozen_manifest.csv", index=False)
    correction_log = pd.DataFrame(
        columns=["date", "canonical_name", "field", "old_value", "new_value", "reason", "commit"]
    )
    correction_log.to_csv(outdir / "baseline_unambiguous_correction_log.csv", index=False)
    return {
        "freeze_rule": "classification fixed without reference to climatic model results",
        "n_species": int(len(strict)),
        "n_families": int(strict["family"].nunique()),
        "n_within": int((strict["spatial_scale"] == "within_population").sum()),
        "n_among": int((strict["spatial_scale"] == "among_population").sum()),
        "sha256": digest,
        "columns": columns,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260719)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.dataset)
    required = {
        "canonical_name",
        "family",
        "spatial_scale",
        "classification_source",
        "n_climate_cells",
        *METRICS,
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    eligible = data.loc[data["n_climate_cells"] >= 20].copy()
    analysis_sets = {
        "all_evidence_classified": eligible,
        "baseline_unambiguous_only": eligible.loc[
            eligible["classification_source"] == "baseline_unambiguous"
        ].copy(),
        "high_confidence_enrichment_only": eligible.loc[
            eligible["classification_source"] == "high_confidence_enrichment"
        ].copy(),
    }

    enrichment_audit = analysis_sets["high_confidence_enrichment_only"].copy()
    enrichment_audit["scale_priority"] = enrichment_audit["spatial_scale"].map(
        {"among_population": 0, "within_population": 1}
    ).fillna(2)
    enrichment_audit = enrichment_audit.sort_values(
        ["scale_priority", "n_climate_cells", "canonical_name"],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    enrichment_audit["audit_priority_rank"] = range(1, len(enrichment_audit) + 1)
    enrichment_audit["audit_reason"] = (
        "Classification derives from high-confidence enrichment and materially affects the "
        "source-stratified moisture-breadth comparison; verify explicit within- versus "
        "among-population evidence in the retained source."
    )
    audit_columns = [
        "audit_priority_rank", "canonical_name", "family", "spatial_scale",
        "classification_source", "n_climate_cells", "pca_hull_area",
        "moisture_breadth", "audit_reason",
    ]
    enrichment_audit[audit_columns].to_csv(
        outdir / "climatic_niche_spatial_scale_enrichment_audit.csv", index=False
    )

    frozen_manifest = write_frozen_manifest(eligible, outdir)
    rng = np.random.default_rng(args.seed)
    rows: list[dict] = []
    loo_rows: list[dict] = []
    for name, frame in analysis_sets.items():
        set_rows, set_loo = analyse_set(frame, name, METRICS, args.permutations, rng)
        rows.extend(set_rows)
        loo_rows.extend(set_loo)

    summary = pd.DataFrame(rows)
    loo_table = pd.DataFrame(loo_rows)
    summary.to_csv(outdir / "climatic_niche_spatial_scale_robustness.csv", index=False)
    summary.to_csv(outdir / "manuscript_spatial_scale_model_table.csv", index=False)
    loo_table.to_csv(outdir / "climatic_niche_spatial_scale_leave_one_family_out.csv", index=False)

    manifest = {
        "min_cells": 20,
        "minimum_model_species": MIN_MODEL_SPECIES,
        "model_formula": MODEL_FORMULA,
        "estimator": "statsmodels GLM Binomial(logit)",
        "covariance": "family-clustered sandwich",
        "confidence_interval_method": CI_METHOD,
        "statsmodels_version": statsmodels.__version__,
        "metrics": METRICS,
        "frozen_baseline_manifest": frozen_manifest,
        "analysis_sets": {
            name: {
                "n_species": int(len(frame)),
                "n_families": int(frame["family"].nunique()),
                "n_within": int((frame["spatial_scale"] == "within_population").sum()),
                "n_among": int((frame["spatial_scale"] == "among_population").sum()),
                "estimable": bool(
                    len(frame) >= MIN_MODEL_SPECIES
                    and frame["spatial_scale"].nunique() >= 2
                    and frame["family"].nunique() >= 2
                ),
            }
            for name, frame in analysis_sets.items()
        },
        "enrichment_audit_rows": int(len(enrichment_audit)),
        "permutations": args.permutations,
        "seed": args.seed,
        "results": summary.to_dict("records"),
        "interpretation_guard": (
            "Exploratory robustness analysis of two focal reported metrics; source-stratified "
            "estimates diagnose classification sensitivity rather than establish effect heterogeneity. "
            "Confidence intervals use family-clustered sandwich uncertainty. Permutation p-values and "
            "leave-one-family-out stability do not establish causation, local adaptation, or morph-specific tolerance."
        ),
    }
    (outdir / "climatic_niche_spatial_scale_robustness_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
