#!/usr/bin/env python3
"""Two-stage observation/organization analysis for the rebuilt full FCP universe.

Stage A diagnoses whether spatial organization (C or S) is documented at all. Stage B
conditions on documented organization and compares C-only, S-only and C+S without
pretending unresolved species are biological negatives. Documentation propensity weights
are estimated from outcome-independent effort variables and used as an IPW sensitivity.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression

METRICS=["temperature_breadth","moisture_breadth","climatic_heterogeneity","pca_dispersion","pca_hull_area"]
CLASSES=["local_coexistence_only","spatial_segregation_only","coexistence_and_segregation"]


def z(v: pd.Series) -> pd.Series:
    x=pd.to_numeric(v,errors='coerce').astype(float)
    sd=float(x.std(ddof=0))
    return (x-float(x.mean()))/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)


def scope_data(analysis: pd.DataFrame, core: pd.DataFrame, attention: pd.DataFrame, scope: str) -> pd.DataFrame:
    d=analysis.copy()
    if scope=='core':
        c=core.copy()
        keep=['canonical_name','organization_state','C_local_coexistence_documented','S_spatial_segregation_documented','n_FCP_eligible_sources']
        c=c[keep].rename(columns={x:f'core_{x}' for x in keep if x!='canonical_name'})
        d=d.merge(c,on='canonical_name',how='inner',validate='one_to_one')
        for x in ['organization_state','C_local_coexistence_documented','S_spatial_segregation_documented','n_FCP_eligible_sources']:
            d[x]=d[f'core_{x}']
    d=d.merge(attention[['canonical_name','n_v22_exact_name_records']],on='canonical_name',how='left',validate='one_to_one')
    d['n_v22_exact_name_records']=d.n_v22_exact_name_records.fillna(0)
    d['D_documented']=((pd.to_numeric(d.C_local_coexistence_documented,errors='coerce').fillna(0)>0) |
                       (pd.to_numeric(d.S_spatial_segregation_documented,errors='coerce').fillna(0)>0)).astype(int)
    d['log_sources_z']=z(np.log1p(pd.to_numeric(d.n_FCP_eligible_sources,errors='coerce')))
    d['attention_z']=z(np.log1p(pd.to_numeric(d.n_v22_exact_name_records,errors='coerce')))
    d['geo_z']=z(np.log1p(pd.to_numeric(d.geographic_radius_95_km,errors='coerce')))
    return d


def fit_documentation(d: pd.DataFrame, metric: str) -> dict:
    x=d.copy(); x['metric_z']=z(x[metric])
    x=x.dropna(subset=['D_documented','metric_z','geo_z','log_sources_z','attention_z','family'])
    X=sm.add_constant(x[['metric_z','geo_z','log_sources_z','attention_z']],has_constant='add')
    try:
        fit=sm.GLM(x.D_documented.astype(int),X,family=sm.families.Binomial()).fit(
            cov_type='cluster',cov_kwds={'groups':x.family.astype(str)})
        return {'status':'complete','n_species':len(x),'n_families':x.family.nunique(),
                'climate_estimate':float(fit.params['metric_z']),'climate_or':float(np.exp(fit.params['metric_z'])),
                'climate_ci_low':float(np.exp(fit.conf_int().loc['metric_z',0])),
                'climate_ci_high':float(np.exp(fit.conf_int().loc['metric_z',1])),
                'climate_p':float(fit.pvalues['metric_z']),
                'geo_or':float(np.exp(fit.params['geo_z'])),'source_count_or':float(np.exp(fit.params['log_sources_z'])),
                'attention_or':float(np.exp(fit.params['attention_z']))}
    except Exception as e:
        return {'status':'fit_failed','n_species':len(x),'error':f'{type(e).__name__}:{e}'[:300]}


def propensity_weights(d: pd.DataFrame) -> pd.DataFrame:
    x=d.dropna(subset=['D_documented','log_sources_z','attention_z']).copy()
    X=sm.add_constant(x[['log_sources_z','attention_z']],has_constant='add')
    fit=sm.GLM(x.D_documented.astype(int),X,family=sm.families.Binomial()).fit()
    p=np.clip(np.asarray(fit.predict(X),dtype=float),0.05,0.95)
    prev=float(x.D_documented.mean())
    x['documentation_probability']=p
    x['stabilized_ipw']=np.where(x.D_documented.eq(1),prev/p,(1-prev)/(1-p))
    return x[['canonical_name','D_documented','documentation_probability','stabilized_ipw']]


def pairwise_from_model(model: LogisticRegression, classes: list[str], metric_index: int=0) -> dict[str,float]:
    coef={c:float(model.coef_[list(model.classes_).index(c),metric_index]) for c in classes}
    return {
        'S_vs_C_OR':float(np.exp(coef['spatial_segregation_only']-coef['local_coexistence_only'])),
        'mixed_vs_C_OR':float(np.exp(coef['coexistence_and_segregation']-coef['local_coexistence_only'])),
        'mixed_vs_S_OR':float(np.exp(coef['coexistence_and_segregation']-coef['spatial_segregation_only'])),
    }


def fit_multinomial(d: pd.DataFrame, metric: str, rng: np.random.Generator, nboot: int) -> tuple[dict,list[dict]]:
    x=d.loc[d.D_documented.eq(1) & d.organization_state.isin(CLASSES)].copy()
    x['metric_z']=z(x[metric])
    x=x.dropna(subset=['metric_z','geo_z','family','stabilized_ipw'])
    if set(x.organization_state)!=set(CLASSES):
        return {'status':'missing_class','n_species':len(x)},[]
    X=x[['metric_z','geo_z']].to_numpy(float); y=x.organization_state.astype(str).to_numpy(); w=x.stabilized_ipw.to_numpy(float)
    model=LogisticRegression(C=1.0,solver='lbfgs',max_iter=2000)
    model.fit(X,y,sample_weight=w)
    point=pairwise_from_model(model,CLASSES)
    families=sorted(x.family.astype(str).unique())
    boots=[]
    for b in range(nboot):
        draws=rng.choice(families,size=len(families),replace=True)
        parts=[]
        for j,f in enumerate(draws):
            q=x.loc[x.family.astype(str).eq(str(f))].copy(); q['_cluster_draw']=j; parts.append(q)
        if not parts: continue
        q=pd.concat(parts,ignore_index=True)
        if set(q.organization_state)!=set(CLASSES): continue
        try:
            m=LogisticRegression(C=1.0,solver='lbfgs',max_iter=2000)
            m.fit(q[['metric_z','geo_z']].to_numpy(float),q.organization_state.astype(str).to_numpy(),sample_weight=q.stabilized_ipw.to_numpy(float))
            vals=pairwise_from_model(m,CLASSES); vals['bootstrap']=b; boots.append(vals)
        except Exception:
            pass
    out={'status':'complete','n_species':len(x),'n_families':len(families),**point,'family_bootstraps_valid':len(boots)}
    for key in ['S_vs_C_OR','mixed_vs_C_OR','mixed_vs_S_OR']:
        vals=np.asarray([r[key] for r in boots if np.isfinite(r[key])],float)
        if len(vals):
            out[key+'_ci_low']=float(np.quantile(vals,.025)); out[key+'_ci_high']=float(np.quantile(vals,.975))
    return out,boots


def main():
    p=argparse.ArgumentParser(); p.add_argument('--analysis',required=True); p.add_argument('--core-species',required=True); p.add_argument('--attention',required=True); p.add_argument('--outdir',required=True); p.add_argument('--family-bootstraps',type=int,default=499); p.add_argument('--seed',type=int,default=20260826); a=p.parse_args()
    analysis=pd.read_csv(a.analysis); core=pd.read_csv(a.core_species); attention=pd.read_csv(a.attention); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(a.seed)
    doc_rows=[]; multi_rows=[]; boot_rows=[]; prop_rows=[]
    for scope in ['core','expanded']:
        d=scope_data(analysis,core,attention,scope)
        prop=propensity_weights(d); d=d.merge(prop[['canonical_name','documentation_probability','stabilized_ipw']],on='canonical_name',how='left',validate='one_to_one')
        q=prop.copy(); q['scope']=scope; prop_rows.append(q)
        for metric in METRICS:
            r=fit_documentation(d,metric); r.update({'scope':scope,'metric':metric}); doc_rows.append(r)
            m,b=fit_multinomial(d,metric,rng,a.family_bootstraps); m.update({'scope':scope,'metric':metric}); multi_rows.append(m)
            for x in b: x.update({'scope':scope,'metric':metric}); boot_rows.extend(b)
    pd.DataFrame(doc_rows).to_csv(out/'full_fcp_documentation_models.csv',index=False)
    pd.DataFrame(multi_rows).to_csv(out/'full_fcp_conditional_multinomial.csv',index=False)
    pd.DataFrame(boot_rows).to_csv(out/'full_fcp_conditional_multinomial_bootstrap.csv',index=False)
    pd.concat(prop_rows,ignore_index=True).to_csv(out/'full_fcp_documentation_propensity.csv',index=False)
    manifest={'status':'complete','design':'two-stage hurdle: documentation process then conditional three-state organization','scopes':['core','expanded'],'metrics':METRICS,'documentation_predictors':'metric + geographic extent + FCP source count + outcome-independent literature attention','conditional_model':'IPW L2 multinomial on metric + geographic extent; C-only/S-only/C+S; family bootstrap','unresolved_semantics':'not treated as biological absence','family_bootstraps':a.family_bootstraps,'seed':a.seed}
    (out/'full_fcp_hurdle_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if len(doc_rows)!=10 or len(multi_rows)!=10: raise SystemExit('Expected 10 documentation and 10 conditional rows')
    print(json.dumps({'documentation_rows':len(doc_rows),'conditional_rows':len(multi_rows),'bootstrap_rows':len(boot_rows)},indent=2))

if __name__=='__main__': main()
