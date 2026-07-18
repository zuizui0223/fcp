#!/usr/bin/env python3
"""Research analysis for natural flower-colour polymorphism."""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def zscore(s):
    x=pd.to_numeric(s,errors='coerce'); sd=x.std(ddof=0)
    return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0

def load(path): return pd.read_csv(Path(path))

def build_labels(species,review):
    d=species.copy(); accepted=set(review.loc[review.review_priority.isin(['P0','P1','P2','P3']),'canonical_name'])
    strict=set(review.loc[review.review_priority.isin(['P0','P1']),'canonical_name'])
    d['validated_case']=d.canonical_name.isin(accepted).astype(int)
    d['strict_case']=d.canonical_name.isin(strict).astype(int)
    for c in ['n_works','n_title_matches','n_context_matches','max_score','total_score']:
        d[c]=pd.to_numeric(d.get(c,0),errors='coerce').fillna(0)
    d['log_works']=np.log1p(d.n_works); d['family']=d.family.fillna('Unknown').astype(str)
    return d

def fit_effort(d,out):
    linear=smf.glm('validated_case ~ log_works',d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':d.family})
    quad=smf.glm('validated_case ~ log_works + I(log_works**2)',d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':d.family})
    strict=smf.glm('strict_case ~ log_works + I(log_works**2)',d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':d.family})
    ci=quad.conf_int()
    pd.DataFrame({'term':quad.params.index,'estimate':quad.params.values,'std_error':quad.bse.values,'p_value':quad.pvalues.values,'odds_ratio':np.exp(quad.params.values),'ci_low':np.exp(ci[0].values),'ci_high':np.exp(ci[1].values)}).to_csv(out/'discovery_effort_nonlinear_model.csv',index=False)
    pred=[]
    for n in [1,2,3,5,10,20]:
        lw=np.log1p(n); pred.append({'n_works':n,'predicted_validation_probability':float(quad.predict(pd.DataFrame({'log_works':[lw]}))[0])})
    pd.DataFrame(pred).to_csv(out/'effort_standardized_predictions.csv',index=False)
    d=d.copy(); d['predicted']=quad.predict(d); d['var_term']=d.predicted*(1-d.predicted)
    fam=d.groupby('family').agg(candidate_species=('validated_case','size'),observed=('validated_case','sum'),expected=('predicted','sum'),variance=('var_term','sum'),mean_works=('n_works','mean')).reset_index()
    fam['standardized_residual']=(fam.observed-fam.expected)/np.sqrt(np.maximum(fam.variance,1e-9))
    fam.sort_values('standardized_residual',ascending=False).to_csv(out/'family_effort_adjusted_signal.csv',index=False)
    loo=[]
    for family,n in d.family.value_counts().items():
        if n<10: continue
        x=d[d.family!=family]; fit=smf.glm('validated_case ~ log_works + I(log_works**2)',x,family=sm.families.Binomial()).fit()
        loo.append({'held_out_family':family,'n_held_out':int(n),'linear_term':float(fit.params['log_works']),'quadratic_term':float(fit.params['I(log_works ** 2)'])})
    pd.DataFrame(loo).to_csv(out/'leave_one_family_out.csv',index=False)
    return {'n_species':int(len(d)),'validated_cases':int(d.validated_case.sum()),'strict_cases':int(d.strict_case.sum()),'candidate_conditional_fraction':float(d.validated_case.mean()),'linear_aic':float(linear.aic),'quadratic_aic':float(quad.aic),'delta_aic_linear_minus_quadratic':float(linear.aic-quad.aic),'strict_quadratic_p':float(strict.pvalues['I(log_works ** 2)'])}

def fit_macro(d,covariates,out):
    p=Path(covariates)
    if not p.exists(): return {'status':'not_run','reason':'covariate table absent'}
    cov=pd.read_csv(p); x=d.merge(cov,on='canonical_name',how='inner',validate='one_to_one')
    req=['island','absolute_latitude','life_history','pollination_system']; missing=[c for c in req if c not in x]
    if missing:return {'status':'not_run','missing':missing,'matched_species':int(len(x))}
    x['absolute_latitude_z']=zscore(x.absolute_latitude); x['log_works_z']=zscore(x.log_works); x['island']=pd.to_numeric(x.island,errors='coerce')
    x=x.dropna(subset=['validated_case','island','absolute_latitude_z','life_history','pollination_system','log_works_z'])
    formula='validated_case ~ island + absolute_latitude_z + C(life_history) + C(pollination_system) + log_works_z'
    m=smf.glm(formula,x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.family}); ci=m.conf_int()
    pd.DataFrame({'term':m.params.index,'estimate':m.params.values,'odds_ratio':np.exp(m.params.values),'ci_low':np.exp(ci[0].values),'ci_high':np.exp(ci[1].values),'p_value':m.pvalues.values}).to_csv(out/'macroecology_model.csv',index=False)
    return {'status':'complete','matched_species':int(len(x)),'formula':formula,'aic':float(m.aic)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--species',default='data/global_flower_colour_species_ranked.csv'); ap.add_argument('--review',default='data/global_flower_colour_review_queue.csv'); ap.add_argument('--covariates',default='data/species_macroecology_covariates.csv'); ap.add_argument('--outdir',default='analysis_outputs/research'); a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); d=build_labels(load(a.species),load(a.review))
    manifest={'research_question':'Which ecological and geographic conditions predict documented natural flower-colour polymorphism after controlling for ascertainment?','estimand_current':'Validation probability conditional on candidature','effort_model':fit_effort(d,out),'macroecology_model':fit_macro(d,a.covariates,out),'guardrail':'Candidate-only data cannot estimate global prevalence.'}
    (out/'research_analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n'); print(json.dumps(manifest,ensure_ascii=False))
if __name__=='__main__': main()
