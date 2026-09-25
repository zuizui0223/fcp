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
PREDICTIONS={
    'bio5': {'direction': 1, 'meaning': 'warmer maximum temperature -> more white (heat/pigment-suppression hypothesis)'},
    'bio14': {'direction': 1, 'meaning': 'wetter driest month -> more white (reduced drought-stress hypothesis)'},
    'srad_mean': {'direction': -1, 'meaning': 'greater long-term solar radiation -> less white (pigment photoprotection hypothesis)'},
}

def bool_series(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({'true','1','yes','y'})

def holm(pvals):
    vals=[(i,float(p)) for i,p in enumerate(pvals) if np.isfinite(p)]
    vals.sort(key=lambda z:z[1]); m=len(vals); out=[np.nan]*len(pvals); prev=0.0
    for rank,(i,p) in enumerate(vals):
        v=min(1.0,(m-rank)*p); v=max(v,prev); out[i]=v; prev=v
    return out

def sample_raster(path, lon, lat):
    with rasterio.open(path) as src:
        pts=list(zip(lon.astype(float),lat.astype(float)))
        vals=np.array([v[0] for v in src.sample(pts)],dtype=float)
        if src.nodata is not None: vals[np.isclose(vals,src.nodata)]=np.nan
        return vals

def within_z(df,col):
    g=df.groupby('species')[col]
    mu=g.transform('mean'); sd=g.transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (df[col]-mu)/sd

def conditional_fit(df,predictor):
    x=df.dropna(subset=['white','species',predictor,'wz_near_clip_fraction']).copy()
    st=x.groupby('species').white.agg(['sum','count']); good=set(st.index[(st['sum']>0)&(st['sum']<st['count'])])
    x=x[x.species.isin(good)].copy()
    groups=x.species.astype('category').cat.codes.to_numpy()
    X=x[[predictor,'wz_near_clip_fraction']].to_numpy(float)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit=ConditionalLogit(x.white.to_numpy(float),X,groups=groups).fit(method='bfgs',maxiter=400,disp=False)
    b=float(fit.params[0]); se=float(fit.bse[0])
    return {'n_rows':int(len(x)),'n_species':int(x.species.nunique()),'beta':b,'se':se,'OR_per_within_species_SD':float(np.exp(b)),'ci_low':float(np.exp(b-1.96*se)),'ci_high':float(np.exp(b+1.96*se)),'p':float(fit.pvalues[0]),'near_clip_beta':float(fit.params[1]),'near_clip_p':float(fit.pvalues[1])}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True); p.add_argument('--technical-table',required=True); p.add_argument('--join-key',required=True); p.add_argument('--high-clip-ids',required=True)
    p.add_argument('--bio-dir',required=True); p.add_argument('--srad-dir',required=True); p.add_argument('--outdir',required=True)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    use=['photo_id','species','morph','global_classifiable','latitude','longitude']
    d=pd.read_csv(a.measured,usecols=use); d=d[bool_series(d.global_classifiable)&d.morph.isin(MORPHS)].copy(); d['white']=(d.morph=='white').astype(int)
    tech=pd.read_csv(a.technical_table,compression='gzip',dtype={'measurement_id':str}); join=pd.read_csv(a.join_key,dtype={'measurement_id':str}); high=pd.read_csv(a.high_clip_ids,dtype={'measurement_id':str})
    tech=tech.merge(join,on='measurement_id',how='left',validate='one_to_one'); highset=set(high.measurement_id.astype(str)); tech['high_clip']=tech.measurement_id.astype(str).isin(highset)
    d=d.merge(tech[['photo_id','near_clip_fraction','high_clip']],on='photo_id',how='left',validate='one_to_one')
    if d.near_clip_fraction.notna().mean()<0.85: raise SystemExit('highlight coverage below 85%')
    for var in ['bio5','bio14']:
        rp=Path(a.bio_dir)/f"wc2.1_10m_bio_{var.replace('bio','')}.tif"
        if not rp.exists(): raise SystemExit(f'missing {rp}')
        d[var]=sample_raster(rp,d.longitude,d.latitude)
    sr=[]
    for m in range(1,13):
        rp=Path(a.srad_dir)/f'wc2.1_10m_srad_{m:02d}.tif'
        if not rp.exists(): raise SystemExit(f'missing {rp}')
        sr.append(sample_raster(rp,d.longitude,d.latitude))
    d['srad_mean']=np.nanmean(np.vstack(sr),axis=0)
    for col in ['bio5','bio14','srad_mean','near_clip_fraction']:
        d['wz_'+col]=within_z(d,col)
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    counts=d.groupby('species').white.agg(['sum','count']); counts['nonwhite']=counts['count']-counts['sum']
    eligible=set(counts.index[(counts['sum']>=5)&(counts.nonwhite>=5)])
    d=d[d.species.isin(eligible)].copy()
    if len(eligible)<100: raise SystemExit(f'underpowered mechanism panel: {len(eligible)} eligible species')
    rows=[]; deltas=[]
    for var,spec in PREDICTIONS.items():
        zc='wz_'+var
        per=[]
        for sp,g in d.dropna(subset=[zc]).groupby('species'):
            if g.white.sum()<5 or (len(g)-g.white.sum())<5: continue
            delta=float(g.loc[g.white.eq(1),zc].mean()-g.loc[g.white.eq(0),zc].mean())
            per.append((sp,delta,len(g),int(g.white.sum())))
        vals=np.array([x[1] for x in per],float)
        if len(vals)<100: raise SystemExit(f'{var}: fewer than 100 species with estimable delta')
        w=wilcoxon(vals,zero_method='wilcox',alternative='two-sided')
        s=binomtest(int((vals*spec['direction']>0).sum()),len(vals),0.5,alternative='greater')
        cond=conditional_fit(d,zc)
        r={'variable':var,'prediction_direction':spec['direction'],'meaning':spec['meaning'],'n_species_delta':int(len(vals)),'median_delta_white_minus_nonwhite_SD':float(np.median(vals)),'mean_delta_white_minus_nonwhite_SD':float(np.mean(vals)),'predicted_direction_fraction':float((vals*spec['direction']>0).mean()),'wilcoxon_two_sided_p':float(w.pvalue),'directional_sign_p':float(s.pvalue),**cond}
        rows.append(r)
        for sp,delta,n,nw in per: deltas.append({'species':sp,'variable':var,'delta_white_minus_nonwhite_SD':delta,'n':n,'n_white':nw})
    adj=holm([r['wilcoxon_two_sided_p'] for r in rows])
    for r,padj in zip(rows,adj):
        r['wilcoxon_holm_p']=padj
        predicted=(np.sign(r['median_delta_white_minus_nonwhite_SD'])==r['prediction_direction'])
        cond_pred=(np.sign(r['beta'])==r['prediction_direction'])
        r['mechanism_gate_pass']=bool(predicted and cond_pred and padj<0.05 and r['p']<0.05)
    manifest={
      'schema':'fcp_white_environment_mechanism_v1','status':'complete','role':'prospectively_specified_environmental_filter_test_after_exploratory_solar_geometry_only','primary_outcome':'white vs nonwhite among frozen globally classifiable third-cohort rows','primary_design':'within-species environmental contrasts; response-blind high-clip exclusion; continuous near-clip adjustment in conditional logistic corroboration','eligibility':'>=5 white and >=5 nonwhite classifiable rows per species after technical availability and high-clip exclusion','minimum_species':100,'multiplicity':'Holm across three prespecified environmental mechanism variables using species-level Wilcoxon p-values','hard_nonclaims':['does not establish evolutionary transition direction','does not establish genetic pigment-loss mechanism','does not test pollinator causation','WorldClim long-term environment is not the weather experienced on the photo date','positive association is environmental sorting evidence, not causal proof'],'rows_after_primary_filters':int(len(d)),'eligible_species':int(len(eligible)),'results':rows
    }
    pd.DataFrame(deltas).to_csv(out/'species_environment_deltas.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'environment_mechanism_models.csv',index=False)
    (out/'result.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()
