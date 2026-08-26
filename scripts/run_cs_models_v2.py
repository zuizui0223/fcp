#!/usr/bin/env python3
"""Canonical primary C/S climatic models.

C and S are separate positive documented-evidence outcomes. Mixed species are positive
for both. The only covariate in primary models besides the focal climatic metric is
occupied-climate sampling effort. Outcome-path-derived source counts are deliberately
not used as observation-effort covariates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

METRICS = ['temperature_breadth','moisture_breadth','climatic_heterogeneity','pca_dispersion','pca_hull_area']
OUTCOMES = {'C_local_coexistence_documented':'C','S_spatial_segregation_documented':'S'}


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()


def zscore(x: pd.Series) -> pd.Series:
    x=pd.to_numeric(x,errors='coerce'); sd=x.std(ddof=0)
    return (x-x.mean())/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index,dtype=float)


def prepare(data: pd.DataFrame,outcome:str,metric:str)->pd.DataFrame:
    d=data.sort_values('canonical_name',kind='stable').reset_index(drop=True).copy()
    d['outcome']=pd.to_numeric(d[outcome],errors='coerce')
    d['metric_z']=zscore(d[metric])
    d['effort_z']=zscore(np.log1p(pd.to_numeric(d['n_climate_cells'],errors='coerce')))
    return d.dropna(subset=['outcome','metric_z','effort_z','family']).copy()


def fit(d:pd.DataFrame,clustered:bool):
    if len(d)<20 or d.outcome.nunique()<2 or d.family.nunique()<2:return None
    X=sm.add_constant(d[['metric_z','effort_z']],has_constant='add')
    model=sm.GLM(d.outcome,X,family=sm.families.Binomial())
    return model.fit(cov_type='cluster',cov_kwds={'groups':d.family}) if clustered else model.fit()


def beta_unclustered(d:pd.DataFrame)->float:
    f=fit(d,False)
    return float(f.params['metric_z']) if f is not None else np.nan


def one(data,outcome,metric,permutations,rng):
    d=prepare(data,outcome,metric); f=fit(d,True)
    if f is None: raise RuntimeError(f'not estimable: {outcome} {metric}')
    b=float(f.params['metric_z']); se=float(f.bse['metric_z']); lo=b-1.96*se; hi=b+1.96*se
    labels=d.outcome.to_numpy().copy(); perm=[]
    for _ in range(permutations):
        x=d.copy(); x['outcome']=rng.permutation(labels); bp=beta_unclustered(x)
        if np.isfinite(bp):perm.append(bp)
    perm=np.asarray(perm,float)
    pp=float((1+np.sum(np.abs(perm)>=abs(b)))/(1+len(perm)))
    loo=[]; loo_rows=[]
    for fam in sorted(d.family.astype(str).unique()):
        sub=d.loc[d.family.astype(str)!=fam].copy(); bl=beta_unclustered(sub); loo.append(bl)
        loo_rows.append({'outcome':outcome,'outcome_short':OUTCOMES[outcome],'metric':metric,'omitted_family':fam,'estimate':bl,'odds_ratio':float(np.exp(bl)) if np.isfinite(bl) else np.nan})
    lv=np.asarray([x for x in loo if np.isfinite(x)],float)
    return {
        'outcome':outcome,'outcome_short':OUTCOMES[outcome],'metric':metric,'analysis_status':'complete',
        'model_formula':f"{OUTCOMES[outcome]} ~ metric_z + effort_z",'estimator':'statsmodels GLM Binomial(logit)',
        'covariance':'family-clustered sandwich','statsmodels_version':statsmodels.__version__,
        'n_species':int(len(d)),'n_families':int(d.family.nunique()),'n_positive':int(d.outcome.sum()),'n_negative':int(len(d)-d.outcome.sum()),
        'estimate':b,'std_error_clustered':se,'estimate_ci_low':lo,'estimate_ci_high':hi,
        'odds_ratio':float(np.exp(b)),'odds_ratio_ci_low':float(np.exp(lo)),'odds_ratio_ci_high':float(np.exp(hi)),
        'wald_p_value_clustered':float(f.pvalues['metric_z']),'permutation_p_two_sided':pp,
        'permutations_requested':int(permutations),'permutations_valid':int(len(perm)),
        'loo_min_odds_ratio':float(np.exp(np.min(lv))) if len(lv) else np.nan,
        'loo_max_odds_ratio':float(np.exp(np.max(lv))) if len(lv) else np.nan,
        'loo_same_direction_fraction':float(np.mean(np.sign(lv)==np.sign(b))) if len(lv) else np.nan,
        'semantic_guard':'Outcome is positive documented evidence; zero means no positive evidence in the frozen audit, not demonstrated biological absence.'
    },loo_rows


def holm(df):
    out=df.copy(); out['wald_p_holm_within_outcome']=np.nan; out['permutation_p_holm_within_outcome']=np.nan
    for outcome,idx in out.groupby('outcome').groups.items():
        ids=list(idx); w=out.loc[ids,'wald_p_value_clustered'].astype(float); p=out.loc[ids,'permutation_p_two_sided'].astype(float)
        out.loc[ids,'wald_p_holm_within_outcome']=multipletests(w,method='holm')[1]
        out.loc[ids,'permutation_p_holm_within_outcome']=multipletests(p,method='holm')[1]
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--dataset',required=True); ap.add_argument('--outdir',required=True); ap.add_argument('--permutations',type=int,default=9999); ap.add_argument('--seed',type=int,default=20260826); a=ap.parse_args()
    path=Path(a.dataset); data=pd.read_csv(path)
    required={'canonical_name','family','n_climate_cells',*METRICS,*OUTCOMES}; missing=required-set(data.columns)
    if missing:raise SystemExit(f'missing: {sorted(missing)}')
    rng=np.random.default_rng(a.seed); rows=[]; loos=[]
    for outcome in OUTCOMES:
        for metric in METRICS:
            r,l=one(data,outcome,metric,a.permutations,rng); rows.append(r); loos.extend(l)
    results=holm(pd.DataFrame(rows)); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    results.to_csv(out/'cs_five_metric_models.csv',index=False); pd.DataFrame(loos).to_csv(out/'cs_leave_one_family_out.csv',index=False)
    manifest={'dataset':str(path),'dataset_sha256':sha256(path),'n_species':int(len(data)),'n_families':int(data.family.nunique()),'C_positive':int(data.C_local_coexistence_documented.sum()),'S_positive':int(data.S_spatial_segregation_documented.sum()),'C_and_S_positive':int(((data.C_local_coexistence_documented==1)&(data.S_spatial_segregation_documented==1)).sum()),'metrics':METRICS,'permutations':a.permutations,'seed':a.seed,'primary_formulae':['C ~ metric_z + effort_z','S ~ metric_z + effort_z'],'observation_process_guard':'n_resolved_sources is outcome-path-derived and is not used as an effort covariate. Independent literature-attention sensitivity is run separately.','semantic_guard':'C and S are separate positive documented-evidence outcomes; neither is the complement of the other.'}
    (out/'cs_model_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(results[['outcome_short','metric','odds_ratio','permutation_p_two_sided','permutation_p_holm_within_outcome','loo_same_direction_fraction']].to_string(index=False))
    if len(results)!=10 or not (results.permutations_valid==a.permutations).all():raise SystemExit('model output invariant failed')

if __name__=='__main__':main()
