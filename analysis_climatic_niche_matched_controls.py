#!/usr/bin/env python3
"""Test whether verified FCP species occupy broader climatic niches than matched controls."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit

METRICS=['pca_dispersion','climatic_heterogeneity','pca_hull_area','temperature_breadth','moisture_breadth']

def z(s):
    x=pd.to_numeric(s,errors='coerce'); sd=x.std(ddof=0); return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--review',required=True);ap.add_argument('--controls',required=True);ap.add_argument('--metrics',required=True);ap.add_argument('--candidate-species',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args()
    review=pd.read_csv(a.review); controls=pd.read_csv(a.controls); metrics=pd.read_csv(a.metrics); candidates=set(pd.read_csv(a.candidate_species).canonical_name.astype(str))
    controls=controls[controls.focal_species.isin(set(review.canonical_name)) & ~controls.control_candidate.isin(candidates)].copy()
    fm=metrics[metrics.role.eq('focal')].drop_duplicates('canonical_name').set_index('canonical_name')
    cm=metrics[metrics.role.eq('control')].drop_duplicates('canonical_name').set_index('canonical_name')
    rows=[]
    for _,r in controls.iterrows():
        f=r.focal_species;c=r.control_candidate
        if f not in fm.index or c not in cm.index: continue
        for outcome,name,m in [(1,f,fm.loc[f]),(0,c,cm.loc[c])]:
            row={'stratum':f,'canonical_name':name,'outcome':outcome,'match_level':r.match_level,'family':r.focal_family,'n_climate_cells':m.n_climate_cells,'metric_status':m.metric_status}
            for x in METRICS:row[x]=m.get(x,np.nan)
            rows.append(row)
    d=pd.DataFrame(rows).drop_duplicates(['stratum','canonical_name','outcome']);Path(a.outdir).mkdir(parents=True,exist_ok=True);d.to_csv(Path(a.outdir)/'climatic_niche_matched_dataset.csv',index=False)
    results=[]
    for level in ['all','same_genus']:
      base=d if level=='all' else d[d.match_level.eq('same_genus')]
      for minimum in [10,20,30,50]:
       for metric in METRICS:
        x=base[(base.metric_status=='complete')&(base.n_climate_cells>=minimum)].dropna(subset=[metric]).copy()
        valid=x.groupby('stratum').outcome.agg(['min','max','count']); valid=valid[(valid['min']==0)&(valid['max']==1)&(valid['count']>=2)].index;x=x[x.stratum.isin(valid)]
        if len(valid)<20: continue
        x['metric_z']=z(np.log1p(x[metric].clip(lower=0)) if metric=='pca_hull_area' else x[metric]);x['effort_z']=z(np.log1p(x.n_climate_cells))
        try:
            model=ConditionalLogit(x.outcome,x[['metric_z','effort_z']],groups=x.stratum).fit(disp=False)
            ci=model.conf_int();results.append({'match_level':level,'min_cells':minimum,'metric':metric,'n_strata':len(valid),'n_rows':len(x),'estimate':model.params.metric_z,'odds_ratio':float(np.exp(model.params.metric_z)),'ci_low':float(np.exp(ci.loc['metric_z',0])),'ci_high':float(np.exp(ci.loc['metric_z',1])),'p_value':model.pvalues.metric_z,'effort_or':float(np.exp(model.params.effort_z))})
        except Exception as e:results.append({'match_level':level,'min_cells':minimum,'metric':metric,'n_strata':len(valid),'error':str(e)[:200]})
    out=pd.DataFrame(results);out.to_csv(Path(a.outdir)/'climatic_niche_conditional_logit.csv',index=False)
    primary=out[(out.match_level=='same_genus')&(out.min_cells==20)] if len(out) else out
    manifest={'dataset_rows':len(d),'focal_species':int(d[d.outcome==1].canonical_name.nunique()) if len(d) else 0,'control_species':int(d[d.outcome==0].canonical_name.nunique()) if len(d) else 0,'specifications':len(out),'primary_results':primary.to_dict('records'),'semantic_guard':'associations concern verified FCP versus taxonomically matched species outside the candidate set; controls are not asserted monomorphic; occupied niche breadth is not physiological tolerance or causal expansion'}
    (Path(a.outdir)/'climatic_niche_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))
if __name__=='__main__':main()
