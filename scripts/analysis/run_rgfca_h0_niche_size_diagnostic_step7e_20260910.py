#!/usr/bin/env python3
from __future__ import annotations

import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import statsmodels.api as sm

ROOT=Path(__file__).resolve().parents[2]
MET=ROOT/'results/rgfca_worldclim_old34_transfer_step7d_20260910/species_environment_metrics.csv'
DISC=ROOT/'data/derived/global_monte_carlo_measured_photos_v1.csv'
RES=ROOT/'data/derived/rgfca_reserve_replication_measured_photos_v1.csv'
OUT=ROOT/'results/rgfca_h0_niche_size_diagnostic_step7e_20260910'; OUT.mkdir(parents=True,exist_ok=True)
N=20000; SEED=20260910; R=6371.0088


def centroid_radius95(df):
    rows=[]
    for sp,g in df.groupby('species',sort=True):
        lat=np.radians(g.latitude.to_numpy(float)); lon=np.radians(g.longitude.to_numpy(float))
        x=np.cos(lat)*np.cos(lon); y=np.cos(lat)*np.sin(lon); z=np.sin(lat)
        cx,cy,cz=x.mean(),y.mean(),z.mean(); clon=math.atan2(cy,cx); clat=math.atan2(cz,math.hypot(cx,cy))
        a=np.sin((lat-clat)/2)**2 + np.cos(lat)*math.cos(clat)*np.sin((lon-clon)/2)**2
        d=2*R*np.arcsin(np.sqrt(np.clip(a,0,1)))
        rows.append({'species':sp,'geographic_radius_95_km':float(np.quantile(d,.95))})
    return pd.DataFrame(rows)


def perm_delta(v,s,seed):
    v=np.asarray(v,float); s=np.asarray(s,bool); obs=float(np.median(v[s])-np.median(v[~s])); rng=np.random.default_rng(seed); null=np.empty(N)
    for i in range(N):
        q=s[rng.permutation(len(s))]; null[i]=np.median(v[q])-np.median(v[~q])
    p=float((1+np.sum(np.abs(null)>=abs(obs)-1e-15))/(N+1))
    return {'delta_S_minus_C':obs,'p':p}


def hc3_fit(y,s,x):
    X=sm.add_constant(np.c_[s.astype(float),x]); fit=sm.OLS(y,X).fit(cov_type='HC3'); return fit


def freedman_lane(y,s,x,seed):
    Xr=sm.add_constant(x); red=sm.OLS(y,Xr).fit(); fitted=np.asarray(red.fittedvalues); resid=np.asarray(red.resid)
    obs=float(hc3_fit(y,s,x).params[1]); rng=np.random.default_rng(seed); null=np.empty(N)
    for i in range(N):
        yp=fitted+resid[rng.permutation(len(resid))]; null[i]=sm.OLS(yp,sm.add_constant(np.c_[s.astype(float),x])).fit().params[1]
    p=float((1+np.sum(np.abs(null)>=abs(obs)-1e-15))/(N+1)); return obs,p


def strat_perm(y,s,x,seed):
    ranks=pd.Series(x).rank(method='first').to_numpy(); strata=np.asarray(pd.qcut(ranks,4,labels=False)); obs=float(np.median(y[s])-np.median(y[~s])); rng=np.random.default_rng(seed); null=np.empty(N)
    idxs=[np.where(strata==k)[0] for k in sorted(set(strata))]
    for b in range(N):
        q=s.copy()
        for idx in idxs: q[idx]=s[idx][rng.permutation(len(idx))]
        null[b]=np.median(y[q])-np.median(y[~q])
    p=float((1+np.sum(np.abs(null)>=abs(obs)-1e-15))/(N+1)); return obs,p


def run_one(tranche,coords_path,seed):
    m=pd.read_csv(MET); m=m[(m.tranche==tranche)&m.organization_state.isin(['local_cooccurrence_only','spatial_segregation_only'])&(m.metric_status=='complete')].copy()
    c=pd.read_csv(coords_path,usecols=['species','latitude','longitude']).dropna(); geo=centroid_radius95(c[c.species.isin(m.species)])
    d=m.merge(geo,on='species',how='inner',validate='one_to_one'); d['S']=d.organization_state.eq('spatial_segregation_only'); d['x']=np.log1p(d.geographic_radius_95_km); d['y']=np.log1p(d.hv3d_rarefied20)
    s=d.S.to_numpy(bool); x=d.x.to_numpy(float); y=d.y.to_numpy(float); xz=(x-x.mean())/x.std(ddof=0)
    e1=perm_delta(x,s,seed)
    pr=pearsonr(x,y); sr=spearmanr(x,y)
    fit=hc3_fit(y,s,xz); coef=float(fit.params[1]); ci=fit.conf_int(alpha=.05)[1]
    flcoef,flp=freedman_lane(y,s,xz,seed+1); sd,sp=strat_perm(y,s,x,seed+2)
    out={'tranche':tranche,'n':len(d),'n_C':int((~s).sum()),'n_S':int(s.sum()),'E1_extent_state':e1,'E2_hv_extent':{'pearson_r':float(pr.statistic),'pearson_p':float(pr.pvalue),'spearman_rho':float(sr.statistic),'spearman_p':float(sr.pvalue)},'E3_extent_adjusted':{'S_coef_log1p_hv':coef,'hc3_ci_low':float(ci[0]),'hc3_ci_high':float(ci[1]),'freedman_lane_coef':float(flcoef),'freedman_lane_p':float(flp)},'E4_radius_stratified':{'delta_S_minus_C_log1p_hv':float(sd),'p':float(sp)}}
    d.to_csv(OUT/f'{tranche}_diagnostic_table.csv',index=False)
    return out


def main():
    ds=run_one('discovery',DISC,SEED+100); rs=run_one('reserve',RES,SEED+200)
    retain=lambda z: z['E3_extent_adjusted']['freedman_lane_coef']>0 and z['E3_extent_adjusted']['freedman_lane_p']<=.05 and z['E4_radius_stratified']['delta_S_minus_C_log1p_hv']>0 and z['E4_radius_stratified']['p']<=.05
    extent_insufficient=bool(retain(ds) and retain(rs))
    extent_explains=bool(ds['E3_extent_adjusted']['freedman_lane_p']>.05 and ds['E4_radius_stratified']['p']>.05 and rs['E3_extent_adjusted']['freedman_lane_p']>.05 and rs['E4_radius_stratified']['p']>.05 and (ds['E1_extent_state']['p']<=.05 or ds['E2_hv_extent']['spearman_p']<=.05) and (rs['E1_extent_state']['p']<=.05 or rs['E2_hv_extent']['spearman_p']<=.05))
    cls='extent_insufficient' if extent_insufficient else ('extent_explains_signal' if extent_explains else 'mixed_or_unresolved')
    result={'analysis':'rgfca_h0_niche_size_diagnostic_step7e','discovery':ds,'reserve':rs,'diagnostic_class':cls,'boundary':'Post-result diagnostic only; no Step-7D hypothesis is reopened.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# RGFCA Step 7E — H0 niche-size diagnostic','']
    for name,z in [('Discovery',ds),('Reserve',rs)]:
        lines += [f'## {name}',f"- n: **{z['n']}** ({z['n_C']} C*, {z['n_S']} S*)",f"- E1 S-C log1p(radius95): **{z['E1_extent_state']['delta_S_minus_C']:.4f}**, p **{z['E1_extent_state']['p']:.6g}**",f"- E2 rho(hypervolume, radius95): **{z['E2_hv_extent']['spearman_rho']:.4f}**, p **{z['E2_hv_extent']['spearman_p']:.6g}**",f"- E3 adjusted S* coefficient: **{z['E3_extent_adjusted']['S_coef_log1p_hv']:.4f}**, HC3 CI [{z['E3_extent_adjusted']['hc3_ci_low']:.4f}, {z['E3_extent_adjusted']['hc3_ci_high']:.4f}], Freedman-Lane p **{z['E3_extent_adjusted']['freedman_lane_p']:.6g}**",f"- E4 radius-stratified delta: **{z['E4_radius_stratified']['delta_S_minus_C_log1p_hv']:.4f}**, p **{z['E4_radius_stratified']['p']:.6g}**",'']
    lines += [f"## Diagnostic class\n\n**{cls}**",'', 'This is a post-result diagnostic and does not reopen Step-7D directional hypotheses.']
    (OUT/'RESULT.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__': main()
