#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, binomtest
from statsmodels.discrete.conditional_models import ConditionalLogit
import rasterio

MORPHS=['white','yellow_orange','red_pink','blue_purple']

def bs(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({'true','1','yes','y'})

def sample(path,lon,lat):
    with rasterio.open(path) as src:
        x=np.array([v[0] for v in src.sample(list(zip(lon.astype(float),lat.astype(float))))],float)
        if src.nodata is not None: x[np.isclose(x,src.nodata)]=np.nan
        return x

def signed_summary(vals):
    vals=np.asarray(vals,float)
    if not len(vals): return {}
    w=wilcoxon(vals,zero_method='wilcox',alternative='two-sided')
    sign=binomtest(int((vals>0).sum()),len(vals),0.5,alternative='two-sided')
    return {
        'n_species':int(len(vals)),
        'median_delta':float(np.median(vals)),
        'mean_delta':float(np.mean(vals)),
        'positive_fraction':float((vals>0).mean()),
        'wilcoxon_two_sided_p':float(w.pvalue),
        'sign_two_sided_p':float(sign.pvalue),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True); p.add_argument('--technical-table',required=True); p.add_argument('--join-key',required=True); p.add_argument('--high-clip-ids',required=True)
    p.add_argument('--bio-dir',required=True); p.add_argument('--srad-dir',required=True); p.add_argument('--outdir',required=True)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.measured,usecols=['photo_id','species','observer_id','morph','global_classifiable','latitude','longitude'])
    d=d[bs(d.global_classifiable)&d.morph.isin(MORPHS)].copy(); d['white']=(d.morph=='white').astype(int)
    t=pd.read_csv(a.technical_table,compression='gzip',dtype={'measurement_id':str}); j=pd.read_csv(a.join_key,dtype={'measurement_id':str}); h=set(pd.read_csv(a.high_clip_ids,dtype={'measurement_id':str}).measurement_id.astype(str))
    t=t.merge(j,on='measurement_id',validate='one_to_one'); t['high_clip']=t.measurement_id.astype(str).isin(h)
    d=d.merge(t[['photo_id','near_clip_fraction','high_clip']],on='photo_id',how='left',validate='one_to_one')
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    for var,num in [('bio5',5),('bio14',14)]:
        d[var]=sample(Path(a.bio_dir)/f'wc2.1_10m_bio_{num}.tif',d.longitude,d.latitude)
    sr=[sample(Path(a.srad_dir)/f'wc2.1_10m_srad_{m:02d}.tif',d.longitude,d.latitude) for m in range(1,13)]
    d['srad_mean']=np.nanmean(np.vstack(sr),axis=0)
    for c in ['bio5','bio14','srad_mean','near_clip_fraction']:
        g=d.groupby('species')[c]; d['wz_'+c]=(d[c]-g.transform('mean'))/g.transform(lambda x:x.std(ddof=0)).replace(0,np.nan)

    paired_rows=[]; paired_species=[]; balanced_rows=[]; balanced_species=[]
    for var in ['bio5','bio14','srad_mean']:
        zc='wz_'+var
        q=d.dropna(subset=[zc,'wz_near_clip_fraction']).copy()

        # Sensitivity 1: strict within-observer comparison. This changes the spatial
        # estimand because only observers who photographed both outcomes contribute.
        ss=q.groupby(['species','observer_id']).white.agg(['sum','count']).reset_index()
        good=ss[(ss['sum']>0)&(ss['sum']<ss['count'])][['species','observer_id']]
        x=q.merge(good,on=['species','observer_id'],how='inner',validate='many_to_one')
        obs=[]
        for (sp,ob),g in x.groupby(['species','observer_id']):
            delta=float(g.loc[g.white.eq(1),zc].mean()-g.loc[g.white.eq(0),zc].mean())
            obs.append({'species':sp,'observer_id':str(ob),'delta':delta,'n':len(g)})
        od=pd.DataFrame(obs)
        sd=od.groupby('species').delta.mean() if len(od) else pd.Series(dtype=float)
        vals=sd.to_numpy(float)
        fitout={'estimable':False}
        if len(x) and x.groupby(['species','observer_id']).ngroups>=30:
            x['group']=x.species.astype(str)+'|'+x.observer_id.astype(str)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    fit=ConditionalLogit(x.white.to_numpy(float),x[[zc,'wz_near_clip_fraction']].to_numpy(float),groups=x.group.astype('category').cat.codes.to_numpy()).fit(method='bfgs',maxiter=400,disp=False)
                b=float(fit.params[0]); se=float(fit.bse[0])
                fitout={'estimable':True,'beta':b,'OR_per_within_species_SD':float(np.exp(b)),'ci_low':float(np.exp(b-1.96*se)),'ci_high':float(np.exp(b+1.96*se)),'p':float(fit.pvalues[0]),'near_clip_beta':float(fit.params[1]),'near_clip_p':float(fit.pvalues[1])}
            except Exception as e:
                fitout={'estimable':False,'error':str(e)[:300]}
        pair_summary=signed_summary(vals)
        paired_rows.append({'variable':var,'role':'posthoc_observer_paired_sensitivity','paired_observer_strata':int(len(od)),**pair_summary,**fitout})
        for sp,v in sd.items(): paired_species.append({'variable':var,'species':sp,'mean_paired_observer_delta':float(v),'n_paired_observers':int((od.species==sp).sum())})

        # Sensitivity 2: observer-balanced broad-range comparison. Each observer can
        # contribute at most one mean value to each colour outcome within a species,
        # preserving between-observer geography while removing unequal observer weight.
        om=q.groupby(['species','observer_id','white'],as_index=False)[zc].mean()
        ocount=om.groupby(['species','white']).observer_id.nunique().unstack(fill_value=0)
        good_species=set(ocount.index[(ocount.get(0,0)>=3)&(ocount.get(1,0)>=3)])
        bvals=[]; brows=[]
        for sp,g in om[om.species.isin(good_species)].groupby('species'):
            wv=g.loc[g.white.eq(1),zc].to_numpy(float)
            nv=g.loc[g.white.eq(0),zc].to_numpy(float)
            delta=float(wv.mean()-nv.mean())
            bvals.append(delta)
            brows.append({'variable':var,'species':sp,'observer_balanced_delta':delta,'n_white_observers':int(len(wv)),'n_nonwhite_observers':int(len(nv))})
        bal_summary=signed_summary(bvals)
        balanced_rows.append({'variable':var,'role':'posthoc_observer_balanced_broad_range_sensitivity','minimum_observers_per_colour':3,**bal_summary})
        balanced_species.extend(brows)

    result={
        'schema':'fcp_white_environment_observer_sensitivity_v2',
        'status':'complete',
        'role':'posthoc sensitivity after the primary WorldClim mechanism result was opened; cannot upgrade or redefine the prospective gate',
        'paired_design':'only species-observer strata containing both white and nonwhite outcomes; strongest observer control but greatly contracts geographic/environmental range',
        'observer_balanced_design':'average environment first within species x observer x colour, then compare colour means across observers; retains broad between-observer geography while equalizing observer weight',
        'paired_results':paired_rows,
        'observer_balanced_results':balanced_rows
    }
    pd.DataFrame(paired_rows).to_csv(out/'observer_paired_sensitivity.csv',index=False)
    pd.DataFrame(paired_species).to_csv(out/'observer_paired_species_deltas.csv',index=False)
    pd.DataFrame(balanced_rows).to_csv(out/'observer_balanced_sensitivity.csv',index=False)
    pd.DataFrame(balanced_species).to_csv(out/'observer_balanced_species_deltas.csv',index=False)
    (out/'observer_paired_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
