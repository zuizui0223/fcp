#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import wilcoxon, binomtest
from statsmodels.discrete.conditional_models import ConditionalLogit
import rasterio

MORPHS=['white','yellow_orange','red_pink','blue_purple']
THRESHOLDS=[25,50,100,250]

def bool_series(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({'true','1','yes','y'})

def haversine_matrix(lat1,lon1,lat2,lon2):
    R=6371.0088
    a1=np.radians(np.asarray(lat1,float))[:,None]; o1=np.radians(np.asarray(lon1,float))[:,None]
    a2=np.radians(np.asarray(lat2,float))[None,:]; o2=np.radians(np.asarray(lon2,float))[None,:]
    da=a2-a1; do=o2-o1
    z=np.sin(da/2)**2+np.cos(a1)*np.cos(a2)*np.sin(do/2)**2
    return 2*R*np.arcsin(np.minimum(1,np.sqrt(z)))

def sample_raster(path,lon,lat):
    with rasterio.open(path) as src:
        vals=np.array([v[0] for v in src.sample(list(zip(lon.astype(float),lat.astype(float))))],float)
        if src.nodata is not None: vals[np.isclose(vals,src.nodata)]=np.nan
        return vals

def build_pairs(d,max_km):
    rows=[]; pair_id=0
    for sp,g in d.groupby('species'):
        w=g[g.white.eq(1)].copy(); n=g[g.white.eq(0)].copy()
        if w.empty or n.empty: continue
        dist=haversine_matrix(w.latitude,w.longitude,n.latitude,n.longitude)
        rr,cc=linear_sum_assignment(dist)
        for i,j in zip(rr,cc):
            km=float(dist[i,j])
            if km>max_km: continue
            rw=w.iloc[i]; rn=n.iloc[j]
            pair_id+=1
            rows.append({'pair_id':pair_id,'species':sp,'distance_km':km,'bio5_diff_white_minus_nonwhite':float(rw.bio5-rn.bio5),'near_clip_diff_white_minus_nonwhite':float(rw.near_clip_fraction-rn.near_clip_fraction)})
    return pd.DataFrame(rows)

def pair_conditional_model(pairs):
    if pairs.empty: return {'estimable':False}
    long=[]
    for r in pairs.itertuples(index=False):
        long.append({'pair_id':r.pair_id,'white':1,'bio5':r.bio5_diff_white_minus_nonwhite/2,'near_clip':r.near_clip_diff_white_minus_nonwhite/2})
        long.append({'pair_id':r.pair_id,'white':0,'bio5':-r.bio5_diff_white_minus_nonwhite/2,'near_clip':-r.near_clip_diff_white_minus_nonwhite/2})
    x=pd.DataFrame(long)
    X=x[['bio5','near_clip']].to_numpy(float)
    groups=x.pair_id.to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit=ConditionalLogit(x.white.to_numpy(float),X,groups=groups).fit(method='bfgs',maxiter=400,disp=False)
    b=float(fit.params[0]); se=float(fit.bse[0])
    return {'estimable':True,'beta_bio5':b,'se_bio5':se,'OR_bio5_native_unit':float(np.exp(b)),'ci_low':float(np.exp(b-1.96*se)),'ci_high':float(np.exp(b+1.96*se)),'p_bio5':float(fit.pvalues[0]),'beta_near_clip':float(fit.params[1]),'p_near_clip':float(fit.pvalues[1])}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True); p.add_argument('--technical-table',required=True); p.add_argument('--join-key',required=True); p.add_argument('--high-clip-ids',required=True); p.add_argument('--bio5',required=True); p.add_argument('--outdir',required=True)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    use=['photo_id','species','morph','global_classifiable','latitude','longitude']
    d=pd.read_csv(a.measured,usecols=use)
    d=d[bool_series(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    d['white']=(d.morph=='white').astype(int)
    tech=pd.read_csv(a.technical_table,compression='gzip',dtype={'measurement_id':str})
    join=pd.read_csv(a.join_key,dtype={'measurement_id':str})
    high=pd.read_csv(a.high_clip_ids,dtype={'measurement_id':str})
    tech=tech.merge(join,on='measurement_id',how='left',validate='one_to_one')
    tech['high_clip']=tech.measurement_id.astype(str).isin(set(high.measurement_id.astype(str)))
    d=d.merge(tech[['photo_id','near_clip_fraction','high_clip']],on='photo_id',how='left',validate='one_to_one')
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    d['bio5']=sample_raster(Path(a.bio5),d.longitude,d.latitude)
    d=d.dropna(subset=['bio5','latitude','longitude'])
    all_pairs=[]; summaries=[]
    for t in THRESHOLDS:
        pairs=build_pairs(d,t)
        if pairs.empty:
            summaries.append({'threshold_km':t,'n_pairs':0,'n_species':0,'estimable':False})
            continue
        sp=pairs.groupby('species').bio5_diff_white_minus_nonwhite.mean()
        vals=sp.to_numpy(float)
        w=wilcoxon(vals,zero_method='wilcox',alternative='two-sided') if np.any(vals!=0) else None
        pos=int((vals>0).sum()); neg=int((vals<0).sum())
        sign=binomtest(pos,pos+neg,0.5,alternative='greater') if pos+neg else None
        cond=pair_conditional_model(pairs)
        summaries.append({'threshold_km':t,'n_pairs':int(len(pairs)),'n_species':int(len(sp)),'median_species_bio5_diff':float(np.median(vals)),'mean_species_bio5_diff':float(np.mean(vals)),'positive_species_fraction':float(pos/(pos+neg)) if pos+neg else np.nan,'wilcoxon_two_sided_p':float(w.pvalue) if w else np.nan,'directional_sign_p':float(sign.pvalue) if sign else np.nan,**cond})
        q=pairs.copy(); q['threshold_km']=t; all_pairs.append(q)
    result={'schema':'fcp_white_heat_local_match_robustness_v1','status':'complete','role':'post_result_geographic_matching_robustness','prospective_gate_changed':False,'hard_nonclaims':['post-result robustness only','does not establish evolutionary transition direction','does not establish causal heat selection','does not alter the prospective BIO5 gate'],'results':summaries}
    if all_pairs: pd.concat(all_pairs,ignore_index=True).to_csv(out/'matched_pairs.csv',index=False)
    pd.DataFrame(summaries).to_csv(out/'local_match_summary.csv',index=False)
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
