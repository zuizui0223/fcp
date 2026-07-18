#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit


def z(s):
    x=pd.to_numeric(s,errors='coerce'); sd=x.std(ddof=0)
    return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0


def valid(d):
    g=d.groupby('focal_species').is_focal.agg(['sum','count'])
    keep=set(g[(g['sum']==1)&(g['count']>=2)].index)
    return d[d.focal_species.isin(keep)].copy()


def nearest(d,caliper):
    out=[]
    for _,g in d.groupby('focal_species',sort=False):
        f=g[g.is_focal==1]
        c=g[(g.is_focal==0)&(g.match_level=='same_genus')].copy()
        if f.empty or c.empty: continue
        lf=float(np.log1p(f.gbif_occurrences.iloc[0]))
        c['effort_distance']=(np.log1p(c.gbif_occurrences)-lf).abs()
        c=c.sort_values(['effort_distance','species']).head(1)
        if caliper is not None and float(c.effort_distance.iloc[0])>caliper: continue
        f=f.copy(); f['effort_distance']=0.0
        out.extend([f,c])
    return pd.concat(out,ignore_index=True) if out else d.iloc[:0].copy()


def fit(d,label,include_gbif):
    d=valid(d)
    d['range_z']=z(np.log1p(d.latitudinal_range.clip(lower=0)))
    d['latitude_z']=z(d.absolute_latitude)
    d['gbif_z']=z(np.log1p(d.gbif_occurrences.clip(lower=0)))
    terms=['range_z','latitude_z']+(['gbif_z'] if include_gbif else [])
    m=ConditionalLogit(d.is_focal,d[terms],groups=d.focal_species).fit(disp=False,maxiter=400)
    ci=m.conf_int(); t=pd.DataFrame({'specification':label,'term':m.params.index,'estimate':m.params.values,'std_error':m.bse.values,'odds_ratio':np.exp(m.params.values),'ci_low':np.exp(ci[0].values),'ci_high':np.exp(ci[1].values),'p_value':m.pvalues.values,'n_rows':len(d),'n_strata':d.focal_species.nunique()})
    r=t.set_index('term').loc['range_z']
    s={'specification':label,'status':'complete','n_rows':len(d),'n_strata':int(d.focal_species.nunique()),'range_odds_ratio':float(r.odds_ratio),'range_ci_low':float(r.ci_low),'range_ci_high':float(r.ci_high),'range_p_value':float(r.p_value),'gbif_adjusted_in_model':include_gbif}
    return t,s


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--dataset',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.dataset)
    for c in ['absolute_latitude','latitudinal_range','gbif_occurrences']:
        d[c]=pd.to_numeric(d[c],errors='coerce')
    d=d.dropna(subset=['absolute_latitude','latitudinal_range','gbif_occurrences'])
    specs=[]
    for label,caliper,adj in [('nearest_same_genus',None,True),('nearest_same_genus_caliper_1',1.0,False),('nearest_same_genus_caliper_0_5',0.5,False)]:
        x=nearest(d,caliper); x.to_csv(out/f'{label}_dataset.csv',index=False)
        if valid(x).focal_species.nunique()>=20: specs.append((label,x,adj))
    tabs=[]; sums=[]
    for label,x,adj in specs:
        t,s=fit(x,label,adj); tabs.append(t); sums.append(s)
    if tabs: pd.concat(tabs,ignore_index=True).to_csv(out/'effort_balanced_conditional_logit.csv',index=False)
    pd.DataFrame(sums).to_csv(out/'effort_balanced_summary.csv',index=False)
    manifest={'semantic_guard':'screening_controls_are_not_asserted_monomorphic','observation_effort_guard':'nearest same-genus control matched on log1p GBIF occurrences','specifications':sums}
    (out/'effort_balanced_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest))

if __name__=='__main__': main()
