#!/usr/bin/env python3
"""Design-based power/precision diagnostic for the frozen 34-species paper."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import brentq

from fcp_pipeline.constants import METRICS
from fcp_pipeline.models import prepare_model_data
from fcp_pipeline.validation import validate_frozen_dataset

TARGET_ORS = [1.0, 0.8, 0.7, 0.6, 0.5, 0.4]


def fit_glm(d: pd.DataFrame, y: np.ndarray, clustered: bool = True):
    X = sm.add_constant(d[["metric_z", "effort_z"]], has_constant="add")
    model = sm.GLM(y, X, family=sm.families.Binomial())
    return model.fit(cov_type="cluster", cov_kwds={"groups": d["family"]}) if clustered else model.fit()


def intercept_for_prevalence(xb: np.ndarray, target: float) -> float:
    def f(a: float) -> float:
        p = 1.0 / (1.0 + np.exp(-(a + xb)))
        return float(p.mean() - target)
    return float(brentq(f, -30.0, 30.0))


def simulate(d, beta_metric, beta_effort, target_prev, reps, rng):
    xb = beta_metric * d["metric_z"].to_numpy() + beta_effort * d["effort_z"].to_numpy()
    intercept = intercept_for_prevalence(xb, target_prev)
    p = 1.0 / (1.0 + np.exp(-(intercept + xb)))
    est, se, pv, cover = [], [], [], []
    failures = 0
    for _ in range(reps):
        y = rng.binomial(1, p)
        if np.unique(y).size < 2:
            failures += 1
            continue
        try:
            fit = fit_glm(d, y, clustered=True)
            b, s = float(fit.params["metric_z"]), float(fit.bse["metric_z"])
            est.append(b); se.append(s); pv.append(float(fit.pvalues["metric_z"]))
            cover.append(b - 1.96*s <= beta_metric <= b + 1.96*s)
        except Exception:
            failures += 1
    e, s, pvals = np.asarray(est), np.asarray(se), np.asarray(pv)
    return {
        "beta_true": beta_metric,
        "odds_ratio_true": float(np.exp(beta_metric)),
        "beta_effort_fixed": beta_effort,
        "target_among_prevalence": target_prev,
        "simulations_requested": reps,
        "simulations_valid": int(len(e)),
        "fit_failures": failures,
        "prob_estimate_negative": float(np.mean(e < 0)),
        "prob_p_lt_0_05": float(np.mean(pvals < .05)),
        "prob_negative_and_p_lt_0_05": float(np.mean((e < 0) & (pvals < .05))),
        "median_estimated_or": float(np.exp(np.median(e))),
        "median_ci_width_logodds": float(np.median(2*1.96*s)),
        "wald_95_coverage": float(np.mean(cover)),
    }


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--dataset",required=True); ap.add_argument("--outdir",required=True); ap.add_argument("--reps",type=int,default=3000); ap.add_argument("--seed",type=int,default=20260812); args=ap.parse_args()
    data=validate_frozen_dataset(pd.read_csv(args.dataset)); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed); target=float((data["spatial_scale"]=="among_population").mean())
    rows=[]; observed=[]
    for metric in METRICS:
        d=prepare_model_data(data,metric); obs=fit_glm(d,d["among"].to_numpy(),clustered=False)
        beta=float(obs.params["metric_z"]); effort=float(obs.params["effort_z"])
        observed.append({"metric":metric,"n_species":len(d),"n_families":d["family"].nunique(),"observed_beta_unclustered":beta,"observed_or_unclustered":float(np.exp(beta)),"observed_effort_beta_unclustered":effort})
        for target_or in TARGET_ORS:
            rows.append({"metric":metric,"scenario":f"OR_{target_or:g}",**simulate(d,np.log(target_or),effort,target,args.reps,rng)})
        rows.append({"metric":metric,"scenario":"observed_effect",**simulate(d,beta,effort,target,args.reps,rng)})
    table=pd.DataFrame(rows); pd.DataFrame(observed).to_csv(out/"jbi_34species_power_precision_observed_design.csv",index=False); table.to_csv(out/"jbi_34species_power_precision_simulation.csv",index=False)
    moisture=table[(table.metric=="moisture_breadth")&(table.scenario=="observed_effect")].iloc[0]
    manifest={"status":"complete","role":"design-based power/precision diagnostic, not an adequacy criterion","n_species":34,"n_families":25,"n_within":20,"n_among":14,"repetitions_per_scenario":args.reps,"seed":args.seed,"target_or_grid":TARGET_ORS,"moisture_observed_effect_diagnostic":moisture.to_dict()}
    (out/"jbi_34species_power_precision_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8"); print(json.dumps(manifest,indent=2))

if __name__=="__main__": main()
