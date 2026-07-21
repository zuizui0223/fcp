#!/usr/bin/env python3
"""Robustness checks for within- versus among-population climatic niche contrasts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

METRICS = ["pca_hull_area", "moisture_breadth"]
MIN_MODEL_SPECIES = 20


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def fit_beta(d: pd.DataFrame, metric: str) -> float:
    x = d.copy()
    x["among"] = (x["spatial_scale"] == "among_population").astype(int)
    x["metric_z"] = zscore(x[metric])
    x["effort_z"] = zscore(np.log1p(pd.to_numeric(x["n_climate_cells"], errors="coerce")))
    x = x.dropna(subset=["among", "metric_z", "effort_z"])
    if len(x) < MIN_MODEL_SPECIES or x["among"].nunique() < 2:
        return np.nan
    model = sm.GLM(
        x["among"],
        sm.add_constant(x[["metric_z", "effort_z"]], has_constant="add"),
        family=sm.families.Binomial(),
    )
    return float(model.fit().params["metric_z"])


def analyse_set(
    data: pd.DataFrame,
    analysis_set: str,
    metrics: list[str],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    loo_rows: list[dict] = []
    estimable = len(data) >= MIN_MODEL_SPECIES and data["spatial_scale"].nunique() >= 2

    for metric in metrics:
        observed = fit_beta(data, metric) if estimable else np.nan
        if estimable:
            permuted = []
            for _ in range(permutations):
                xp = data.copy()
                xp["spatial_scale"] = rng.permutation(xp["spatial_scale"].to_numpy())
                permuted.append(fit_beta(xp, metric))
            permuted = np.asarray(permuted, dtype=float)
            valid = permuted[np.isfinite(permuted)]
        else:
            valid = np.asarray([], dtype=float)

        p_two = (
            float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))
            if np.isfinite(observed)
            else np.nan
        )

        families = sorted(data["family"].dropna().astype(str).unique())
        loo = []
        if estimable:
            for fam in families:
                beta = fit_beta(data.loc[data["family"].astype(str) != fam], metric)
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
        loo_arr = np.asarray(loo, dtype=float)
        valid_loo = loo_arr[np.isfinite(loo_arr)]
        rows.append(
            {
                "analysis_set": analysis_set,
                "metric": metric,
                "analysis_status": "complete" if estimable else "not_estimable",
                "analysis_reason": "" if estimable else f"fewer than {MIN_MODEL_SPECIES} species or only one spatial-scale class",
                "n_species": int(len(data)),
                "n_within": int((data["spatial_scale"] == "within_population").sum()),
                "n_among": int((data["spatial_scale"] == "among_population").sum()),
                "estimate": observed,
                "odds_ratio": float(np.exp(observed)) if np.isfinite(observed) else np.nan,
                "permutation_p_two_sided": p_two,
                "permutations_valid": int(len(valid)),
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
        "audit_priority_rank",
        "canonical_name",
        "family",
        "spatial_scale",
        "classification_source",
        "n_climate_cells",
        "pca_hull_area",
        "moisture_breadth",
        "audit_reason",
    ]
    enrichment_audit[audit_columns].to_csv(
        outdir / "climatic_niche_spatial_scale_enrichment_audit.csv", index=False
    )

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
    loo_table.to_csv(outdir / "climatic_niche_spatial_scale_leave_one_family_out.csv", index=False)

    manifest = {
        "min_cells": 20,
        "minimum_model_species": MIN_MODEL_SPECIES,
        "metrics": METRICS,
        "analysis_sets": {
            name: {
                "n_species": int(len(frame)),
                "n_within": int((frame["spatial_scale"] == "within_population").sum()),
                "n_among": int((frame["spatial_scale"] == "among_population").sum()),
                "estimable": bool(len(frame) >= MIN_MODEL_SPECIES and frame["spatial_scale"].nunique() >= 2),
            }
            for name, frame in analysis_sets.items()
        },
        "enrichment_audit_rows": int(len(enrichment_audit)),
        "permutations": args.permutations,
        "results": summary.to_dict("records"),
        "interpretation_guard": (
            "Exploratory robustness analysis of two preselected near-threshold metrics; "
            "source-stratified estimates diagnose classification sensitivity rather than establish "
            "effect heterogeneity. Subsets below the prespecified 20-species minimum are retained "
            "for audit but not assigned inferential estimates. Permutation p-values and leave-one-family-out "
            "stability do not establish causation or physiological niche breadth."
        ),
    }
    (outdir / "climatic_niche_spatial_scale_robustness_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
