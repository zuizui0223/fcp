#!/usr/bin/env python3
"""Research analysis for natural flower-colour polymorphism.

The pipeline has two operational layers:
1. a discovery/ascertainment model that is runnable on the current evidence table;
2. a biological macroecology model that activates when species covariates are supplied.

The current repository contains candidate species rather than a random sample of all
angiosperms, so the discovery layer estimates validation conditional on candidature.
It must not be interpreted as global prevalence.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def zscore(s: pd.Series) -> pd.Series:
    x=pd.to_numeric(s,errors='coerce')
    sd=x.std(ddof=0)
    return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0


def load_table(path: str) -> pd.DataFrame:
    p=Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(p)


def build_case_labels(species: pd.DataFrame, review: pd.DataFrame) -> pd.DataFrame:
    d=species.copy()
    accepted=set(review.loc[review['review_priority'].isin(['P0','P1','P2','P3']),'canonical_name'])
    d['validated_case']=d['canonical_name'].isin(accepted).astype(int)
    for c in ['n_works','n_title_matches','n_context_matches','max_score','total_score']:
        d[c]=pd.to_numeric(d.get(c,0),errors='coerce').fillna(0)
    d['log_works']=np.log1p(d['n_works'])
    d['family']=d['family'].fillna('Unknown').astype(str)
    return d


def fit_discovery_model(d: pd.DataFrame, out: Path) -> dict:
    # Effort-only model: tests how strongly validation is driven by literature effort.
    m1=smf.glm('validated_case ~ log_works',data=d,family=sm.families.Binomial()).fit()
    coef=pd.DataFrame({'term':m1.params.index,'estimate':m1.params.values,'std_error':m1.bse.values,
                       'z':m1.tvalues.values,'p_value':m1.pvalues.values})
    coef.to_csv(out/'discovery_effort_model.csv',index=False)

    # Family deviations after standardising to the same literature effort.
    fam=d.groupby('family').agg(candidate_species=('validated_case','size'),validated=('validated_case','sum'),
                                mean_log_works=('log_works','mean')).reset_index()
    fam['expected_probability']=m1.predict(pd.DataFrame({'log_works':fam['mean_log_works']}))
    fam['expected_validated']=fam['candidate_species']*fam['expected_probability']
    fam['standardized_residual']=(fam['validated']-fam['expected_validated'])/np.sqrt(
        np.maximum(fam['candidate_species']*fam['expected_probability']*(1-fam['expected_probability']),1e-9))
    fam.sort_values('standardized_residual',ascending=False).to_csv(out/'family_effort_adjusted_signal.csv',index=False)

    # Leave-one-family-out stability: refit excluding each sufficiently represented family.
    rows=[]
    for family,n in d['family'].value_counts().items():
        if n<10: continue
        sub=d[d['family']!=family]
        fit=smf.glm('validated_case ~ log_works',data=sub,family=sm.families.Binomial()).fit()
        rows.append({'held_out_family':family,'n_held_out':int(n),'effort_log_odds':float(fit.params['log_works']),
                     'effort_p_value':float(fit.pvalues['log_works'])})
    pd.DataFrame(rows).to_csv(out/'leave_one_family_out.csv',index=False)

    return {'n_species':int(len(d)),'validated_cases':int(d.validated_case.sum()),
            'case_fraction_conditional_on_candidate':float(d.validated_case.mean()),
            'effort_log_odds':float(m1.params['log_works']),'effort_odds_ratio':float(math.exp(m1.params['log_works'])),
            'effort_p_value':float(m1.pvalues['log_works']),'aic':float(m1.aic)}


def fit_macroecology(d: pd.DataFrame, covariates_path: str, out: Path) -> dict:
    p=Path(covariates_path)
    if not p.exists():
        return {'status':'not_run','reason':'species covariate table not present'}
    cov=pd.read_csv(p)
    x=d.merge(cov,on='canonical_name',how='inner',validate='one_to_one')
    required=['island','absolute_latitude','life_history','pollination_system']
    missing=[c for c in required if c not in x]
    if missing:
        return {'status':'not_run','reason':'missing required covariates','missing':missing,'matched_species':int(len(x))}
    x['absolute_latitude_z']=zscore(x['absolute_latitude'])
    x['log_works_z']=zscore(x['log_works'])
    x['island']=pd.to_numeric(x['island'],errors='coerce')
    x=x.dropna(subset=['validated_case','island','absolute_latitude_z','life_history','pollination_system','log_works_z'])
    formula='validated_case ~ island + absolute_latitude_z + C(life_history) + C(pollination_system) + log_works_z'
    model=smf.glm(formula,data=x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x['family']})
    ci=model.conf_int()
    table=pd.DataFrame({'term':model.params.index,'estimate':model.params.values,'odds_ratio':np.exp(model.params.values),
                        'ci_low':np.exp(ci[0].values),'ci_high':np.exp(ci[1].values),'p_value':model.pvalues.values})
    table.to_csv(out/'macroecology_model.csv',index=False)
    return {'status':'complete','matched_species':int(len(x)),'formula':formula,'aic':float(model.aic)}


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--species',default='data/global_flower_colour_species_ranked.csv')
    ap.add_argument('--review',default='data/global_flower_colour_review_queue.csv')
    ap.add_argument('--covariates',default='data/species_macroecology_covariates.csv')
    ap.add_argument('--outdir',default='analysis_outputs/research')
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=build_case_labels(load_table(a.species),load_table(a.review))
    discovery=fit_discovery_model(d,out)
    macro=fit_macroecology(d,a.covariates,out)
    manifest={'research_question':'Which ecological and geographic conditions predict documented natural flower-colour polymorphism after controlling for publication effort?',
              'estimand_current':'Probability of validation conditional on entering the candidate set',
              'estimand_target':'Probability of natural flower-colour polymorphism among sampled angiosperm species',
              'discovery_model':discovery,'macroecology_model':macro,
              'guardrail':'Candidate-only analyses cannot estimate global prevalence or causal ecological effects.'}
    (out/'research_analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(manifest,ensure_ascii=False))

if __name__=='__main__':
    main()
