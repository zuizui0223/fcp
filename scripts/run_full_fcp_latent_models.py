#!/usr/bin/env python3
"""Detection-aware latent C/S models for the rebuilt all-species FCP universe.

Strict C/S source evidence is high-specificity but incomplete. A zero at species level is
therefore a non-detection, not a biological absence. This script models that explicitly.

For each species i and axis A (C or S):
  z_i ~ Bernoulli(psi_i)
  logit(psi_i) = alpha + beta_climate * climate_i + beta_geo * geography_i

For each eligible source j belonging to species i:
  y_ij | z_i=1 ~ Bernoulli(p_A)
  y_ij | z_i=0 ~ Bernoulli(epsilon)

The latent z_i is analytically marginalized. epsilon defaults to zero because the strict
positive rules were designed for high specificity; fixed nonzero epsilon values can be
used as false-positive sensitivities. Source count enters as repeated detection
opportunities rather than as an outcome-derived regression covariate.

Primary uncertainty is calibrated by a parametric bootstrap of the climate likelihood-
ratio statistic under beta_climate=0. Family-cluster bootstrap intervals and leave-one-
family-out fits assess taxonomic dependence. Core and expanded universes are fit with the
same climate/geography dataset, but their source replicates are kept distinct: the core
uses only display-core-eligible source rows, whereas expanded uses all expanded-eligible
source rows.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from statsmodels.stats.multitest import multipletests

METRICS = [
    "temperature_breadth",
    "moisture_breadth",
    "climatic_heterogeneity",
    "pca_dispersion",
    "pca_hull_area",
]
AXES = {
    "C": "C_local_coexistence_documented_strict",
    "S": "S_spatial_segregation_documented_strict",
}
PRIOR_SD = 5.0


@dataclass
class SpeciesObs:
    name: str
    family: str
    x_metric: float
    x_geo: float
    y: np.ndarray


def zscore(x: pd.Series) -> pd.Series:
    v = pd.to_numeric(x, errors="coerce")
    sd = float(v.std(ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=v.index, dtype=float)
    return (v - float(v.mean())) / sd


def safe_log_bern(y: np.ndarray, p: float) -> float:
    p = min(max(float(p), 1e-12), 1 - 1e-12)
    return float(np.sum(y * math.log(p) + (1 - y) * math.log1p(-p)))


def build_species_obs(
    climate: pd.DataFrame,
    sources: pd.DataFrame,
    species_names: set[str],
    axis_col: str,
    metric: str,
) -> list[SpeciesObs]:
    d = climate.loc[climate.canonical_name.astype(str).isin(species_names)].copy()
    d = d.loc[d.climate_eligible.astype(bool)].copy() if "climate_eligible" in d.columns else d
    d["metric_z"] = zscore(d[metric])
    d["geo_z"] = zscore(np.log1p(pd.to_numeric(d["geographic_radius_95_km"], errors="coerce")))
    d = d.dropna(subset=["canonical_name", "family", "metric_z", "geo_z"])

    s = sources.copy()
    # Core membership-source artifacts and expanded source audits both retain these
    # columns. Filtering again is harmless and protects against accidental extra rows.
    if "FCP_eligible_source" in s.columns:
        s = s.loc[pd.to_numeric(s["FCP_eligible_source"], errors="coerce").fillna(0).eq(1)]
    if "taxon_resolution_status" in s.columns:
        s = s.loc[s.taxon_resolution_status.astype(str).eq("resolved_unique")]
    s = s.loc[s.accepted_name.astype(str).isin(set(d.canonical_name.astype(str)))]
    s["y"] = pd.to_numeric(s[axis_col], errors="coerce").fillna(0).astype(int)

    y_by_species = {
        name: grp.y.to_numpy(dtype=int)
        for name, grp in s.groupby("accepted_name", sort=False)
    }
    out: list[SpeciesObs] = []
    for _, row in d.iterrows():
        name = str(row.canonical_name)
        y = y_by_species.get(name)
        if y is None or len(y) == 0:
            continue
        out.append(
            SpeciesObs(
                name=name,
                family=str(row.family),
                x_metric=float(row.metric_z),
                x_geo=float(row.geo_z),
                y=y,
            )
        )
    return out


def neg_log_posterior(
    theta: np.ndarray,
    obs: list[SpeciesObs],
    *,
    epsilon: float,
    beta_fixed_zero: bool,
    weights: dict[str, float] | None = None,
) -> float:
    if beta_fixed_zero:
        alpha, beta_geo, det_intercept = theta
        beta = 0.0
    else:
        alpha, beta, beta_geo, det_intercept = theta
    p_det = float(expit(det_intercept))
    total = 0.0
    for item in obs:
        psi = float(expit(alpha + beta * item.x_metric + beta_geo * item.x_geo))
        log_l1 = safe_log_bern(item.y, p_det)
        log_l0 = safe_log_bern(item.y, epsilon)
        ll = float(logsumexp([math.log(max(psi, 1e-15)) + log_l1,
                              math.log(max(1.0 - psi, 1e-15)) + log_l0]))
        w = 1.0 if weights is None else float(weights.get(item.family, 0.0))
        total -= w * ll
    total += 0.5 * float(np.sum((theta / PRIOR_SD) ** 2))
    return total


def fit_model(
    obs: list[SpeciesObs],
    *,
    epsilon: float,
    beta_fixed_zero: bool,
    start: np.ndarray | None = None,
    weights: dict[str, float] | None = None,
) -> Any:
    # Kept here for backward provenance; execution workflows use
    # run_full_fcp_latent_models_v2.py, which replaces this optimizer call with an
    # equivalent closure because scipy.optimize.minimize has no kwargs= parameter.
    k = 3 if beta_fixed_zero else 4
    if start is None or len(start) != k:
        start = np.zeros(k, dtype=float)
        start[-1] = 0.0
    fit = minimize(
        neg_log_posterior,
        np.asarray(start, dtype=float),
        args=(obs,),
        method="L-BFGS-B",
        bounds=[(-12, 12)] * k,
        options={"maxiter": 600, "ftol": 1e-10},
    )
    return fit


def simulate_under_null(
    obs: list[SpeciesObs],
    null_theta: np.ndarray,
    *,
    epsilon: float,
    rng: np.random.Generator,
) -> list[SpeciesObs]:
    alpha, beta_geo, det_intercept = null_theta
    p_det = float(expit(det_intercept))
    sim: list[SpeciesObs] = []
    for item in obs:
        psi = float(expit(alpha + beta_geo * item.x_geo))
        z = int(rng.random() < psi)
        p = p_det if z else epsilon
        y = (rng.random(len(item.y)) < p).astype(int)
        sim.append(SpeciesObs(item.name, item.family, item.x_metric, item.x_geo, y))
    return sim


def fit_one(
    obs: list[SpeciesObs],
    *,
    epsilon: float,
    parametric_bootstraps: int,
    family_bootstraps: int,
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    alt = fit_model(obs, epsilon=epsilon, beta_fixed_zero=False)
    null = fit_model(obs, epsilon=epsilon, beta_fixed_zero=True)
    if not alt.success or not null.success:
        return {"analysis_status":"fit_failed","n_species":len(obs),"alt_message":str(alt.message),"null_message":str(null.message)}, [], []
    beta=float(alt.x[1]); beta_geo=float(alt.x[2]); detection=float(expit(alt.x[3])); lr_obs=max(0.0,2.0*(float(null.fun)-float(alt.fun)))
    lr_boot=[]
    for _ in range(parametric_bootstraps):
        sim=simulate_under_null(obs,null.x,epsilon=epsilon,rng=rng)
        a=fit_model(sim,epsilon=epsilon,beta_fixed_zero=False,start=alt.x); n=fit_model(sim,epsilon=epsilon,beta_fixed_zero=True,start=null.x)
        if a.success and n.success and np.isfinite(a.fun) and np.isfinite(n.fun): lr_boot.append(max(0.0,2.0*(float(n.fun)-float(a.fun))))
    lr_arr=np.asarray(lr_boot,float); p_boot=float((1+np.sum(lr_arr>=lr_obs))/(1+len(lr_arr))) if len(lr_arr) else np.nan
    families=sorted({o.family for o in obs}); fam_beta=[]
    for _ in range(family_bootstraps):
        draws=rng.choice(families,size=len(families),replace=True); counts={}
        for f in draws: counts[str(f)]=counts.get(str(f),0.0)+1.0
        f=fit_model(obs,epsilon=epsilon,beta_fixed_zero=False,start=alt.x,weights=counts)
        if f.success and np.isfinite(f.x[1]): fam_beta.append(float(f.x[1]))
    fb=np.asarray(fam_beta,float)
    if len(fb): ci_low,ci_high=np.quantile(fb,[.025,.975]); sign_fraction=float(np.mean(np.sign(fb)==np.sign(beta)))
    else: ci_low=ci_high=sign_fraction=np.nan
    loo_rows=[]; loo_beta=[]
    for family in families:
        sub=[o for o in obs if o.family!=family]; f=fit_model(sub,epsilon=epsilon,beta_fixed_zero=False,start=alt.x); b=float(f.x[1]) if f.success else np.nan
        if np.isfinite(b): loo_beta.append(b)
        loo_rows.append({'omitted_family':family,'estimate':b,'odds_ratio':float(np.exp(b)) if np.isfinite(b) else np.nan,'fit_success':bool(f.success)})
    loo_arr=np.asarray(loo_beta,float); posterior_rows=[]; alpha=float(alt.x[0]); p_det=detection
    for item in obs:
        psi=float(expit(alpha+beta*item.x_metric+beta_geo*item.x_geo)); l1=math.log(max(psi,1e-15))+safe_log_bern(item.y,p_det); l0=math.log(max(1-psi,1e-15))+safe_log_bern(item.y,epsilon); denom=float(logsumexp([l1,l0])); posterior=float(math.exp(l1-denom))
        posterior_rows.append({'canonical_name':item.name,'family':item.family,'n_eligible_sources':int(len(item.y)),'n_positive_sources':int(np.sum(item.y)),'latent_state_prior_probability':psi,'latent_state_posterior_probability':posterior})
    row={'analysis_status':'complete','n_species':int(len(obs)),'n_families':int(len(families)),'n_sources':int(sum(len(o.y) for o in obs)),'n_species_with_positive_source':int(sum(np.any(o.y==1) for o in obs)),'estimate':beta,'odds_ratio':float(np.exp(beta)),'geographic_extent_estimate':beta_geo,'geographic_extent_odds_ratio':float(np.exp(beta_geo)),'estimated_source_detection_probability':detection,'fixed_false_positive_probability':float(epsilon),'penalized_lr_statistic':lr_obs,'parametric_bootstrap_p':p_boot,'parametric_bootstraps_requested':int(parametric_bootstraps),'parametric_bootstraps_valid':int(len(lr_arr)),'family_bootstrap_ci_low':float(np.exp(ci_low)) if np.isfinite(ci_low) else np.nan,'family_bootstrap_ci_high':float(np.exp(ci_high)) if np.isfinite(ci_high) else np.nan,'family_bootstraps_requested':int(family_bootstraps),'family_bootstraps_valid':int(len(fb)),'family_bootstrap_same_direction_fraction':sign_fraction,'loo_min_odds_ratio':float(np.exp(np.min(loo_arr))) if len(loo_arr) else np.nan,'loo_max_odds_ratio':float(np.exp(np.max(loo_arr))) if len(loo_arr) else np.nan,'loo_same_direction_fraction':float(np.mean(np.sign(loo_arr)==np.sign(beta))) if len(loo_arr) else np.nan}
    return row,loo_rows,posterior_rows


def main() -> None:
    p=argparse.ArgumentParser(); p.add_argument('--analysis',required=True); p.add_argument('--source-audit',required=True); p.add_argument('--core-source-audit',required=True); p.add_argument('--core-species',required=True); p.add_argument('--outdir',required=True); p.add_argument('--epsilon',type=float,default=0.0); p.add_argument('--parametric-bootstraps',type=int,default=999); p.add_argument('--family-bootstraps',type=int,default=499); p.add_argument('--seed',type=int,default=20260826); args=p.parse_args()
    climate=pd.read_csv(args.analysis); expanded_sources=pd.read_csv(args.source_audit); core_sources=pd.read_csv(args.core_source_audit); core=pd.read_csv(args.core_species); core_names=set(core.canonical_name.astype(str)); expanded_names=set(climate.canonical_name.astype(str)); rng=np.random.default_rng(args.seed)
    required_climate={'canonical_name','family','geographic_radius_95_km',*METRICS}; missing=required_climate-set(climate.columns)
    if missing: raise SystemExit(f'Analysis dataset missing columns: {sorted(missing)}')
    required_sources={'accepted_name','family',*AXES.values()}
    for label,sdf in [('expanded',expanded_sources),('core',core_sources)]:
        m=required_sources-set(sdf.columns)
        if m: raise SystemExit(f'{label} source audit missing columns: {sorted(m)}')
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True); rows=[]; loo_all=[]; post_all=[]
    scopes={'core':(core_names,core_sources),'expanded':(expanded_names,expanded_sources)}
    for scope,(names,source_df) in scopes.items():
        for axis_short,axis_col in AXES.items():
            for metric in METRICS:
                obs=build_species_obs(climate,source_df,names,axis_col,metric)
                row,loo,posterior=fit_one(obs,epsilon=args.epsilon,parametric_bootstraps=args.parametric_bootstraps,family_bootstraps=args.family_bootstraps,rng=rng); row.update({'scope':scope,'axis':axis_short,'axis_column':axis_col,'metric':metric}); rows.append(row)
                for x in loo: x.update({'scope':scope,'axis':axis_short,'metric':metric}); loo_all.extend(loo)
                for x in posterior: x.update({'scope':scope,'axis':axis_short,'metric':metric}); post_all.extend(posterior)
                print({k:row.get(k) for k in ['scope','axis','metric','n_species','n_sources','odds_ratio','estimated_source_detection_probability','parametric_bootstrap_p']},flush=True)
    results=pd.DataFrame(rows); results['parametric_bootstrap_p_holm_within_scope_axis']=np.nan
    for (scope,axis),idx in results.groupby(['scope','axis']).groups.items():
        idx=list(idx); vals=pd.to_numeric(results.loc[idx,'parametric_bootstrap_p'],errors='coerce'); mask=vals.notna()
        if mask.any(): results.loc[np.asarray(idx)[mask.to_numpy()],'parametric_bootstrap_p_holm_within_scope_axis']=multipletests(vals[mask],method='holm')[1]
    results.to_csv(outdir/'full_fcp_latent_cs_models.csv',index=False); pd.DataFrame(loo_all).to_csv(outdir/'full_fcp_latent_cs_leave_one_family_out.csv',index=False); pd.DataFrame(post_all).to_csv(outdir/'full_fcp_latent_cs_species_posteriors.csv',index=False)
    manifest={'status':'complete','model':'marginalized latent-state repeated-source detection model','scopes':['core','expanded'],'axes':AXES,'metrics':METRICS,'occupancy_formula':'latent C or S ~ metric_z + z(log1p(geographic_radius_95_km))','detection_model':'constant source-level detection probability within each axis/metric/scope fit','core_source_boundary':'core uses only display-core membership sources; expanded uses all expanded FCP-eligible sources','false_positive_probability':float(args.epsilon),'weak_regularization_sd':PRIOR_SD,'parametric_bootstraps':int(args.parametric_bootstraps),'family_cluster_bootstraps':int(args.family_bootstraps),'seed':int(args.seed),'zero_semantics':'source-level zero is a non-detection; species-level biological absence is latent and never assigned directly','historical_34_role':'none in model fitting'}
    (outdir/'full_fcp_latent_cs_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    if len(results)!=20: raise SystemExit(f'Expected 20 model rows, found {len(results)}')
    if not set(results.analysis_status).issubset({'complete','fit_failed'}): raise SystemExit('Unexpected latent model status')

if __name__=='__main__': main()
