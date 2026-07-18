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


def coefficient_table(model):
    ci=model.conf_int()
    return pd.DataFrame({'term':model.params.index,'estimate':model.params.values,'std_error':model.bse.values,
                         'odds_ratio':np.exp(model.params.values),'ci_low':np.exp(ci[0].values),
                         'ci_high':np.exp(ci[1].values),'p_value':model.pvalues.values})


def fit_effort(d,out):
    linear=smf.glm('validated_case ~ log_works',d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':d.family})
    quad=smf.glm('validated_case ~ log_works + I(log_works**2)',d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':d.family})
    strict=smf.glm('strict_case ~ log_works + I(log_works**2)',d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':d.family})
    coefficient_table(quad).to_csv(out/'discovery_effort_nonlinear_model.csv',index=False)
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


def prepare_geography(d,cov):
    x=d.merge(cov,on='canonical_name',how='inner',suffixes=('','_cov'))
    for c in ['absolute_latitude','latitudinal_range','gbif_occurrences','gbif_coordinate_records_sampled']:
        if c in x: x[c]=pd.to_numeric(x[c],errors='coerce')
    x=x.dropna(subset=['validated_case','absolute_latitude','latitudinal_range','gbif_occurrences','log_works','family']).copy()
    x['absolute_latitude_z']=zscore(x.absolute_latitude)
    x['latitudinal_range_z']=zscore(np.log1p(x.latitudinal_range.clip(lower=0)))
    x['latitudinal_range_raw_z']=zscore(x.latitudinal_range.clip(lower=0))
    x['latitudinal_range_sqrt_z']=zscore(np.sqrt(x.latitudinal_range.clip(lower=0)))
    x['log_gbif_occurrences_z']=zscore(np.log1p(x.gbif_occurrences.clip(lower=0)))
    x['log_works_z']=zscore(x.log_works)
    return x


def fit_geo_model(x,response='validated_case',range_term='latitudinal_range_z',cluster=True):
    formula=f'{response} ~ absolute_latitude_z + {range_term} + log_gbif_occurrences_z + log_works_z + I(log_works_z**2)'
    kwargs={}
    if cluster: kwargs={'cov_type':'cluster','cov_kwds':{'groups':x.family}}
    return smf.glm(formula,x,family=sm.families.Binomial()).fit(**kwargs),formula


def fit_geography(d,covariates,out):
    p=Path(covariates)
    if not p.exists(): return {'status':'not_run','reason':'GBIF geography covariate table absent'}
    cov=pd.read_csv(p)
    required=['canonical_name','absolute_latitude','latitudinal_range','gbif_occurrences']
    missing=[c for c in required if c not in cov]
    if missing: return {'status':'not_run','reason':'missing GBIF geography columns','missing':missing}
    x=prepare_geography(d,cov)
    if len(x)<100 or x.validated_case.sum()<15:
        return {'status':'not_run','reason':'insufficient matched geography data','matched_species':int(len(x)),'validated_cases':int(x.validated_case.sum())}
    model,formula=fit_geo_model(x)
    coefficient_table(model).to_csv(out/'geography_macroecology_model.csv',index=False)
    strict,_=fit_geo_model(x,response='strict_case')
    coefficient_table(strict).to_csv(out/'geography_macroecology_strict_model.csv',index=False)
    coverage=pd.DataFrame([{
        'candidate_species_total':len(d),'matched_species':len(x),'matched_fraction':len(x)/len(d),
        'validated_cases_matched':int(x.validated_case.sum()),'strict_cases_matched':int(x.strict_case.sum()),
        'families_matched':int(x.family.nunique()),'median_gbif_occurrences':float(x.gbif_occurrences.median()),
        'median_coordinate_sample':float(x.gbif_coordinate_records_sampled.median()) if 'gbif_coordinate_records_sampled' in x else np.nan
    }])
    coverage.to_csv(out/'geography_covariate_coverage.csv',index=False)

    robust=[]
    specs=[('log_range_all',x,'latitudinal_range_z'),('raw_range_all',x,'latitudinal_range_raw_z'),('sqrt_range_all',x,'latitudinal_range_sqrt_z')]
    if 'gbif_coordinate_records_sampled' in x:
        for threshold in [20,50,100,300]:
            sub=x[x.gbif_coordinate_records_sampled>=threshold].copy()
            if len(sub)>=100 and sub.validated_case.sum()>=15:
                for term in ['latitudinal_range_z','latitudinal_range_raw_z','latitudinal_range_sqrt_z']:
                    specs.append((f'{term}_mincoord_{threshold}',sub,term))
    for label,sub,term in specs:
        fit,_=fit_geo_model(sub,range_term=term)
        robust.append({'specification':label,'n_species':len(sub),'validated_cases':int(sub.validated_case.sum()),
                       'range_term':term,'estimate':float(fit.params[term]),'odds_ratio':float(np.exp(fit.params[term])),
                       'std_error':float(fit.bse[term]),'p_value':float(fit.pvalues[term]),'aic':float(fit.aic)})
    robustness=pd.DataFrame(robust)
    robustness.to_csv(out/'geography_range_robustness.csv',index=False)

    loo=[]
    for family,n in x.family.value_counts().items():
        if n<10: continue
        sub=x[x.family!=family]
        fit,_=fit_geo_model(sub,cluster=False)
        term='latitudinal_range_z'
        loo.append({'held_out_family':family,'n_held_out':int(n),'n_remaining':int(len(sub)),
                    'range_estimate':float(fit.params[term]),'range_odds_ratio':float(np.exp(fit.params[term])),
                    'range_p_value':float(fit.pvalues[term])})
    loo_df=pd.DataFrame(loo)
    loo_df.to_csv(out/'geography_range_leave_one_family_out.csv',index=False)

    range_or=float(np.exp(model.params['latitudinal_range_z']))
    return {'status':'complete','matched_species':int(len(x)),'matched_fraction':float(len(x)/len(d)),
            'validated_cases':int(x.validated_case.sum()),'families':int(x.family.nunique()),
            'formula':formula,'aic':float(model.aic),'strict_aic':float(strict.aic),
            'range_odds_ratio':range_or,'range_p_value':float(model.pvalues['latitudinal_range_z']),
            'strict_range_odds_ratio':float(np.exp(strict.params['latitudinal_range_z'])),
            'strict_range_p_value':float(strict.pvalues['latitudinal_range_z']),
            'robustness_specifications':int(len(robustness)),
            'negative_range_estimate_fraction':float((robustness.estimate<0).mean()),
            'loo_families':int(len(loo_df)),
            'loo_negative_fraction':float((loo_df.range_estimate<0).mean()) if len(loo_df) else None}


def fit_full_macro(d,covariates,out):
    p=Path(covariates)
    if not p.exists(): return {'status':'not_run','reason':'full ecological covariate table absent'}
    cov=pd.read_csv(p); x=d.merge(cov,on='canonical_name',how='inner',validate='one_to_one')
    req=['island','absolute_latitude','life_history','pollination_system']; missing=[c for c in req if c not in x]
    if missing:return {'status':'not_run','missing':missing,'matched_species':int(len(x))}
    x['absolute_latitude_z']=zscore(x.absolute_latitude); x['log_works_z']=zscore(x.log_works); x['island']=pd.to_numeric(x.island,errors='coerce')
    x=x.dropna(subset=['validated_case','island','absolute_latitude_z','life_history','pollination_system','log_works_z'])
    formula='validated_case ~ island + absolute_latitude_z + C(life_history) + C(pollination_system) + log_works_z'
    m=smf.glm(formula,x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.family})
    coefficient_table(m).to_csv(out/'macroecology_model.csv',index=False)
    return {'status':'complete','matched_species':int(len(x)),'formula':formula,'aic':float(m.aic)}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--species',default='data/global_flower_colour_species_ranked.csv'); ap.add_argument('--review',default='data/global_flower_colour_review_queue.csv'); ap.add_argument('--covariates',default='data/species_macroecology_covariates.csv'); ap.add_argument('--geography-covariates',default='data/species_gbif_geography_covariates.csv'); ap.add_argument('--outdir',default='analysis_outputs/research'); a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); d=build_labels(load(a.species),load(a.review))
    manifest={'research_question':'Which ecological and geographic conditions predict documented natural flower-colour polymorphism after controlling for ascertainment?','estimand_current':'Validation probability conditional on candidature','effort_model':fit_effort(d,out),'geography_model':fit_geography(d,a.geography_covariates,out),'full_macroecology_model':fit_full_macro(d,a.covariates,out),'guardrail':'Candidate-only data cannot estimate global prevalence; GBIF coordinate summaries describe documented distributions and are observation-effort sensitive.'}
    (out/'research_analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n'); print(json.dumps(manifest,ensure_ascii=False))
if __name__=='__main__': main()
