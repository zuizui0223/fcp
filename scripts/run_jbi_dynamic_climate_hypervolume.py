#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from sklearn.covariance import LedoitWolf
from sklearn.linear_model import LogisticRegression

CLASSES=['local_coexistence_only','spatial_segregation_only','coexistence_and_segregation']
FEATURES=['tmean','ppt','def','vpd']
METRICS=[
 ('seasonal_centroid_cycle_mean','C','seasonal environmental movement'),
 ('interannual_centroid_drift_mean','C','year-to-year environmental centroid movement'),
 ('interannual_overlap_loss_mean','C','year-to-year Gaussian hypervolume turnover'),
 ('annual_hypervolume_log_sd','C','interannual variability in environmental hypervolume size'),
 ('temporal_variance_component','C','within-location interannual environmental variance'),
 ('spatial_variance_component','S','persistent among-location environmental variance'),
 ('space_time_variance_ratio','S','spatial relative to temporal environmental variance'),
]

def z(s):
    x=pd.to_numeric(s,errors='coerce').astype(float); sd=float(x.std(ddof=0))
    return (x-float(x.mean()))/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)

def gaussian_descriptor(X):
    X=np.asarray(X,float); cov=LedoitWolf().fit(X).covariance_; mu=X.mean(0)
    sign,logdet=np.linalg.slogdet(cov)
    if sign<=0: raise ValueError('non-positive covariance')
    d=X.shape[1]; unit_ball=math.pi**(d/2)/math.gamma(d/2+1); r2=float(chi2.ppf(.95,d))
    volume=unit_ball*(r2**(d/2))*math.exp(.5*logdet)
    return mu,cov,float(volume),float(.5*logdet)

def bhattacharyya(mu1,c1,mu2,c2):
    c=(c1+c2)/2; delta=mu1-mu2
    inv=np.linalg.pinv(c); term1=.125*float(delta@inv@delta)
    s,ld=np.linalg.slogdet(c); s1,ld1=np.linalg.slogdet(c1); s2,ld2=np.linalg.slogdet(c2)
    if min(s,s1,s2)<=0: return np.nan
    db=term1+.5*(ld-.5*(ld1+ld2))
    return float(np.clip(np.exp(-max(db,0)),0,1))

def scale_from_baseline(d,cols,baseline_start,baseline_end,prefix):
    b=d.loc[d.year.between(baseline_start,baseline_end),cols]
    pars={}
    for c in cols:
        mu=float(b[c].mean()); sd=float(b[c].std(ddof=0))
        if not np.isfinite(sd) or sd<=0: raise ValueError(f'bad scaler for {c}')
        d[prefix+c]=(d[c]-mu)/sd; pars[c]={'mean':mu,'sd':sd}
    return pars

def derive_metrics(raw,baseline_start,baseline_end,min_points):
    d=raw.copy(); d['time']=pd.to_datetime(d.time,utc=True,errors='coerce'); d['year']=d.time.dt.year; d['month']=d.time.dt.month
    for c in ['tmin','tmax','ppt','def','vpd']: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['tmean']=(d.tmin+d.tmax)/2
    d=d.dropna(subset=['canonical_name','selection_rank','year','month']+FEATURES)
    monthly_scaler=scale_from_baseline(d,FEATURES,baseline_start,baseline_end,'mz_')
    cent=d.loc[d.year.between(baseline_start,baseline_end)].groupby(['canonical_name','month'])[[f'mz_{c}' for c in FEATURES]].mean().reset_index()
    seasonal={}
    for name,g in cent.groupby('canonical_name'):
        g=g.set_index('month').reindex(range(1,13))
        if g.isna().any().any(): continue
        X=g[[f'mz_{c}' for c in FEATURES]].to_numpy(float); Xn=np.roll(X,-1,axis=0)
        seasonal[name]=float(np.linalg.norm(Xn-X,axis=1).mean())
    group=['canonical_name','family','organization_state','selection_rank','year']
    count=d.groupby(group).month.nunique().rename('n_months')
    annual=d.groupby(group).agg(tmean=('tmean','mean'),ppt=('ppt','sum'),def_=('def','sum'),vpd=('vpd','mean')).reset_index().rename(columns={'def_':'def'})
    annual=annual.merge(count.reset_index(),on=group,validate='one_to_one'); annual=annual.loc[annual.n_months.eq(12)].copy()
    annual_scaler=scale_from_baseline(annual,FEATURES,baseline_start,baseline_end,'az_')
    az=[f'az_{c}' for c in FEATURES]
    rows=[]; hrows=[]
    for (name,year),g in annual.groupby(['canonical_name','year']):
        if g.selection_rank.nunique()<min_points: continue
        try: mu,cov,vol,logscale=gaussian_descriptor(g[az].to_numpy(float))
        except Exception: continue
        hrows.append({'canonical_name':name,'year':int(year),'n_points':len(g),'mu':mu,'cov':cov,'hypervolume_95':vol,'logscale':logscale})
    h=pd.DataFrame(hrows)
    for name,g in annual.groupby('canonical_name'):
        meta=g.iloc[0]; r={'canonical_name':name,'family':meta.family,'organization_state':meta.organization_state,'seasonal_centroid_cycle_mean':seasonal.get(name,np.nan)}
        point_means=g.groupby('selection_rank')[az].mean()
        spatial=float(point_means.var(axis=0,ddof=1).sum()) if len(point_means)>1 else np.nan
        temporal_by_point=g.groupby('selection_rank')[az].var(ddof=1).sum(axis=1)
        temporal=float(temporal_by_point.mean()) if len(temporal_by_point) else np.nan
        r['spatial_variance_component']=spatial; r['temporal_variance_component']=temporal; r['space_time_variance_ratio']=spatial/temporal if temporal>0 else np.nan
        hg=h.loc[h.canonical_name.eq(name)].sort_values('year') if not h.empty else pd.DataFrame()
        drifts=[]; losses=[]
        if len(hg)>=2:
            vals=list(hg.itertuples(index=False))
            for a,b in zip(vals[:-1],vals[1:]):
                if b.year-a.year!=1: continue
                drifts.append(float(np.linalg.norm(b.mu-a.mu)))
                bc=bhattacharyya(a.mu,a.cov,b.mu,b.cov)
                if np.isfinite(bc): losses.append(1-bc)
        r['interannual_centroid_drift_mean']=float(np.mean(drifts)) if drifts else np.nan
        r['interannual_overlap_loss_mean']=float(np.mean(losses)) if losses else np.nan
        r['annual_hypervolume_log_sd']=float(hg.logscale.std(ddof=1)) if len(hg)>=2 else np.nan
        r['n_annual_hypervolumes']=len(hg); r['first_year']=int(g.year.min()); r['last_year']=int(g.year.max())
        rows.append(r)
    return pd.DataFrame(rows),annual,h,{'monthly_scaler':monthly_scaler,'annual_scaler':annual_scaler}

def propensity(core,analysis,attention):
    m=core.merge(analysis[['canonical_name','geographic_radius_95_km']],on='canonical_name',how='inner',validate='one_to_one').merge(attention[['canonical_name','n_v22_exact_name_records']],on='canonical_name',how='left',validate='one_to_one')
    m['n_v22_exact_name_records']=m.n_v22_exact_name_records.fillna(0); m['D_documented']=((m.C_local_coexistence_documented>0)|(m.S_spatial_segregation_documented>0)).astype(int)
    m['log_sources_z']=z(np.log1p(m.n_FCP_eligible_sources)); m['attention_z']=z(np.log1p(m.n_v22_exact_name_records)); m['geo_z']=z(np.log1p(m.geographic_radius_95_km))
    X=sm.add_constant(m[['log_sources_z','attention_z']],has_constant='add'); fit=sm.GLM(m.D_documented,X,family=sm.families.Binomial()).fit(); p=np.clip(np.asarray(fit.predict(X),float),.05,.95); prev=float(m.D_documented.mean())
    m['stabilized_ipw']=np.where(m.D_documented.eq(1),prev/p,(1-prev)/(1-p)); return m

def pairwise(model):
    coef={c:float(model.coef_[list(model.classes_).index(c),0]) for c in CLASSES}
    return {'S_vs_C_OR':float(np.exp(coef['spatial_segregation_only']-coef['local_coexistence_only'])),'mixed_vs_C_OR':float(np.exp(coef['coexistence_and_segregation']-coef['local_coexistence_only'])),'mixed_vs_S_OR':float(np.exp(coef['coexistence_and_segregation']-coef['spatial_segregation_only']))}

def model_metric(d,metric,rng,nboot):
    x=d.loc[d.organization_state.isin(CLASSES)].copy(); x['metric_z']=z(x[metric]); x=x.dropna(subset=['metric_z','geo_z','family','stabilized_ipw'])
    X=x[['metric_z','geo_z']].to_numpy(float); y=x.organization_state.astype(str).to_numpy(); w=x.stabilized_ipw.to_numpy(float)
    model=LogisticRegression(C=1,solver='lbfgs',max_iter=3000).fit(X,y,sample_weight=w); out={'status':'complete','n_species':len(x),'n_families':x.family.nunique(),**pairwise(model)}
    fam=sorted(x.family.astype(str).unique()); boots=[]
    for b in range(nboot):
        draws=rng.choice(fam,size=len(fam),replace=True); parts=[]
        for j,f in enumerate(draws):
            q=x.loc[x.family.astype(str).eq(f)].copy(); q['_draw']=j; parts.append(q)
        q=pd.concat(parts,ignore_index=True)
        if set(q.organization_state)!=set(CLASSES): continue
        try:
            mm=LogisticRegression(C=1,solver='lbfgs',max_iter=3000).fit(q[['metric_z','geo_z']].to_numpy(float),q.organization_state.astype(str).to_numpy(),sample_weight=q.stabilized_ipw.to_numpy(float)); v=pairwise(mm); v['bootstrap']=b; boots.append(v)
        except Exception: pass
    out['family_bootstraps_valid']=len(boots)
    for k in ['S_vs_C_OR','mixed_vs_C_OR','mixed_vs_S_OR']:
        vals=np.asarray([r[k] for r in boots if np.isfinite(r[k])])
        if len(vals): out[k+'_ci_low']=float(np.quantile(vals,.025)); out[k+'_ci_high']=float(np.quantile(vals,.975))
    pval=np.nan; por=np.nan
    q=x.loc[x.organization_state.isin(['local_coexistence_only','spatial_segregation_only'])].copy(); q['S']=(q.organization_state=='spatial_segregation_only').astype(int)
    if q.S.nunique()==2:
        try:
            XX=sm.add_constant(q[['metric_z','geo_z']],has_constant='add'); ff=sm.GLM(q.S,XX,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':q.family.astype(str)})
            por=float(np.exp(ff.params.metric_z)); pval=float(ff.pvalues.metric_z)
        except Exception: pass
    out['pure_S_vs_C_OR']=por; out['pure_S_vs_C_cluster_p']=pval
    return out,boots

def holm(pvals):
    vals=[(i,p) for i,p in enumerate(pvals) if np.isfinite(p)]; vals.sort(key=lambda z:z[1]); m=len(vals); adj=[np.nan]*len(pvals); last=0
    for rank,(i,p) in enumerate(vals):
        v=min(1,(m-rank)*p); v=max(v,last); adj[i]=v; last=v
    return adj

def main():
    p=argparse.ArgumentParser(); p.add_argument('--timeseries',required=True); p.add_argument('--core-species',required=True); p.add_argument('--analysis',required=True); p.add_argument('--attention',required=True); p.add_argument('--outdir',required=True)
    p.add_argument('--baseline-start',type=int,default=1991); p.add_argument('--baseline-end',type=int,default=2020); p.add_argument('--min-points',type=int,default=18); p.add_argument('--family-bootstraps',type=int,default=499); p.add_argument('--seed',type=int,default=20260827); a=p.parse_args()
    raw=pd.read_csv(a.timeseries); core=pd.read_csv(a.core_species); analysis=pd.read_csv(a.analysis); attention=pd.read_csv(a.attention); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    metrics,annual,h,scalers=derive_metrics(raw,a.baseline_start,a.baseline_end,a.min_points)
    base=propensity(core,analysis,attention); d=base.merge(metrics,on=['canonical_name','family','organization_state'],how='left',validate='one_to_one')
    rng=np.random.default_rng(a.seed); rows=[]; brows=[]
    for metric,pred,meaning in METRICS:
        r,b=model_metric(d,metric,rng,a.family_bootstraps); r.update({'metric':metric,'predicted_state_with_higher_metric':pred,'meaning':meaning}); rows.append(r)
        for x in b: x.update({'metric':metric}); brows.extend(b)
    adj=holm([r.get('pure_S_vs_C_cluster_p',np.nan) for r in rows])
    for r,v in zip(rows,adj): r['pure_S_vs_C_cluster_p_holm']=v
    metrics.to_csv(out/'jbi_dynamic_climate_species_metrics.csv',index=False); pd.DataFrame(rows).to_csv(out/'jbi_dynamic_climate_models.csv',index=False); pd.DataFrame(brows).to_csv(out/'jbi_dynamic_climate_family_bootstrap.csv',index=False)
    h.drop(columns=['mu','cov'],errors='ignore').to_csv(out/'jbi_dynamic_annual_hypervolumes.csv',index=False)
    manifest={'status':'complete','protocol':'dynamic-climate-hypervolume-v1','hypothesis_status':'prospectively specified before dynamic TerraClimate result extraction, but not externally preregistered','baseline':[a.baseline_start,a.baseline_end],'features':FEATURES,'metrics':[m[0] for m in METRICS],'hypervolume':'95% shrinkage-Gaussian ellipsoid over 4 globally standardized annual climate dimensions at fixed spatial points','overlap':'Bhattacharyya affinity between consecutive annual Gaussian hypervolumes','space_time_decomposition':'annual standardized environment partitioned into variance among fixed point means (space) and mean within-point variance across years (time)','family_bootstraps':a.family_bootstraps,'seed':a.seed,'scalers':scalers,'TerraClimate_guard':'use one current TerraClimate version only; geometry/variability analysis does not attribute long-term trends independently of parent reanalysis'}
    (out/'jbi_dynamic_climate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__': main()
