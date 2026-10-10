#!/usr/bin/env python3
"""All ORIGINAL FCP photo opportunity: rain-associated colour vs observation selection.

Two propensity estimators (climate+location+source panel; plus observable flower-area
proxy) cross-fit by original photographed 10-degree geographical cell. This
identifies only a MAR-conditional sensitivity, NOT unobserved photo hue. Explicit
MNAR tipping curves and nonrandom missingness strata are the primary safeguards.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,zipfile,io
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score,brier_score_loss

SOURCE_SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
MORPHS=('white','yellow_orange','red_pink','blue_purple')

EXPECTED={'photo_union':100543,'species':42111,'classifiable':45953,'missing':54590,'georain':100255,'cov_complete_classified':45836}
NBOOT=1999
SEED=202610102215
MIN_PROPENSITY=.05
MAX_PROPENSITY=.95

def load_original(archive:Path):
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise RuntimeError('original source archive SHA256 mismatch')
    with zipfile.ZipFile(archive) as z:
        cell=pd.read_csv(io.BytesIO(z.read('fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz')),compression='gzip',low_memory=False)
        breadth=pd.read_csv(io.BytesIO(z.read('fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz')),compression='gzip',low_memory=False)
    if len(cell)!=85337 or len(breadth)!=42111:raise RuntimeError('old source table population drift')
    cell['source_panel']='cell'; breadth['source_panel']='breadth'
    src=pd.concat([cell,breadth],ignore_index=True)
    source_duplicate=src[src.photo_id.duplicated(keep=False)]
    if source_duplicate.photo_id.nunique()!=26905 or source_duplicate.groupby('photo_id').inat_taxon_id.nunique().max()!=1:
        raise RuntimeError('original photo identity conflict')
    assigned=source_duplicate[source_duplicate.measurement_status.eq('classified_four_state_morph')]
    if assigned.groupby('photo_id').morph.nunique().max()>1:
        raise RuntimeError('original photo colour labels contradictory')
    u=src.drop_duplicates('photo_id',keep='first').copy().reset_index(drop=True)
    u['classified']=(u.measurement_status.eq('classified_four_state_morph')&u.morph.isin(MORPHS))
    if len(u)!=EXPECTED['photo_union'] or u.inat_taxon_id.nunique()!=EXPECTED['species'] or int(u.classified.sum())!=EXPECTED['classifiable']:
        raise RuntimeError('100543 original photo union or classifiability drift')
    u['genus']=u.species.astype(str).str.split().str[0]
    u['white']=u.morph.eq('white').astype(int)
    u['site_cell_10degree']=np.floor((u.latitude+90)/10).fillna(-1).astype(int)*37+np.floor((u.longitude+180)/10).fillna(-1).astype(int)
    return u

def features(s:pd.DataFrame,model:str)->np.ndarray:
    lat=s.latitude.to_numpy(float);lon=s.longitude.to_numpy(float)
    v=[np.abs(lat),np.sin(np.deg2rad(lon)),np.cos(np.deg2rad(lon)),
       s.wc_bio1.to_numpy(float),s.wc_bio5.to_numpy(float),s.wc_bio12.to_numpy(float),
       s.wc_bio15.to_numpy(float),s.wc_elevation_m.to_numpy(float),
       s.source_panel.eq('cell').to_numpy(float)]
    if model=='geo_climate_plus_flower_area':
        flower=pd.to_numeric(s.flower_effective_pixels,errors='coerce')
        v.extend([np.log1p(flower.fillna(0).clip(lower=0).to_numpy(float)),
                  flower.isna().to_numpy(float)])
    elif model!='geo_climate':raise ValueError('undefined propensity model')
    a=np.column_stack(v)
    if not np.isfinite(a).all():raise RuntimeError('source environmental feature missing')
    return a

def crossfit(s:pd.DataFrame,model:str):
    X=features(s,model)
    y=s.classified.astype(int).to_numpy()
    group=s.site_cell_10degree.to_numpy(int)
    p=np.full(len(s),np.nan)
    folds=[]
    for i,(tr,te) in enumerate(GroupKFold(n_splits=5).split(X,y,groups=group)):
        scaler=StandardScaler().fit(X[tr])
        xx=scaler.transform(X[tr]); tt=scaler.transform(X[te])
        fit=LogisticRegression(C=1.,max_iter=300,solver='lbfgs',tol=1e-5)
        fit.fit(xx,y[tr])
        p[te]=fit.predict_proba(tt)[:,1]
        folds.append({'fold':i,'train':len(tr),'heldout':len(te),'geographic_cells_heldout':int(np.unique(group[te]).size)})
    if not np.isfinite(p).all():raise RuntimeError('out of region classification models incomplete')
    return p,{'out_of_geographic_region_roc_auc':float(roc_auc_score(y,p)),
              'out_of_geographic_region_brier':float(brier_score_loss(y,p)),
              'true_global_source_photo_classifiability':float(y.mean()),
              'min_probability':float(p.min()),'max_probability':float(p.max()),
              'fraction_clipped_0_05_or_0_95':float(np.mean((p<MIN_PROPENSITY)|(p>MAX_PROPENSITY))),
              'folds':folds}

def fractions(s:pd.DataFrame,weights:np.ndarray,rain_bin=9):
    z=s.loc[s.classified & s.rain_decile.isin((0,rain_bin))].copy()
    w=weights[z.index.to_numpy(int)]
    out={}
    for k in (0,rain_bin):
        m=z.rain_decile.to_numpy(int)==k
        wei=w[m];z1=z.loc[m]
        out[str(k)]={'n_classifiable':int(m.sum()),'n_effective':float(wei.sum()**2/np.square(wei).sum()),
                     'white_fraction':float(np.average(z1.white.to_numpy(float),weights=wei)),
                     'rain_mm_median':float(z1.wc_bio12.median()),
                     'yellow_orange_fraction':float(np.average(z1.morph.eq('yellow_orange').to_numpy(float),weights=wei))}
    return out,out[str(rain_bin)]['white_fraction']-out['0']['white_fraction']

def cluster_bootstrap_original(s, weights,cluster,nboot=NBOOT):
    z=s.loc[s.classified&s.rain_decile.isin((0,9))].copy()
    z['weight']=weights[z.index.to_numpy(int)]
    z['w_white']=z.weight*z.white
    z['wet']=z.rain_decile.eq(9)
    t=z.groupby([cluster,'wet']).agg(white=('w_white','sum'),total=('weight','sum')).unstack(fill_value=0)
    a=np.stack([t[('white',False)].to_numpy(float),t[('total',False)].to_numpy(float),
                t[('white',True)].to_numpy(float),t[('total',True)].to_numpy(float)],axis=1)
    if len(a)<30: raise RuntimeError('insufficient source clusters')
    rng=np.random.default_rng(SEED+sum(map(ord,cluster)))
    idx=rng.integers(len(a),size=(nboot,len(a)))
    agg=a[idx].sum(axis=1)
    if np.any((agg[:,1]<=0)|(agg[:,3]<=0)):
        raise RuntimeError('source rain bin missing from clustered bootstrap')
    vals=agg[:,2]/agg[:,3]-agg[:,0]/agg[:,1]
    return {'n_clusters':len(a),'interval_95':[float(x) for x in np.quantile(vals,[.025,.975])],
            'replicate_count':nboot,'p_one_sided_nonpositive':float((1+(vals<=0).sum())/(nboot+1))}

def tipping(s:pd.DataFrame):
    out={}
    for d in [0.,.25,.5,.75,1.]:
        rows=[]
        for k in (0,9):
            frame=s.loc[s.rain_decile.eq(k)]
            n=len(frame);white=frame.loc[frame.classified,'white'].sum();miss=int((~frame.classified).sum())
            rows.append((n,white,miss))
        nd,wd,md=rows[0];nw,ww,mw=rows[1]
        white_dry=(wd+md*d)/nd
        target=(white_dry*nw-ww)/mw
        out[str(d)]={'assumed_unclassified_dry_white_fraction':d,
                     'required_unclassified_wet_white_fraction_for_EQUAL_true_dry_wet_white':float(target),
                     'feasible_value_between_zero_and_one':bool(0<=target<=1),
                     'dry_unconditional_white_if_assumed':float(white_dry)}
    return out

def quality_strata(s:pd.DataFrame):
    val=pd.to_numeric(s.flower_effective_pixels,errors='coerce')
    if val.isna().sum()>5:raise RuntimeError('original photo quality coverage unexpectedly missing')
    cut=np.nanquantile(val,[.25,.5,.75])
    grades=np.searchsorted(cut,val.fillna(-1),side='right')
    out=[]
    for j in range(4):
        sub=s.loc[(grades==j)&s.rain_decile.isin((0,9))]
        dry=sub.loc[sub.rain_decile.eq(0)];wet=sub.loc[sub.rain_decile.eq(9)]
        n=lambda x:int(x.classified.sum())
        if min(n(dry),n(wet))<100:continue
        out.append({'flower_area_quartile':j+1,'min_px':None if j==0 else float(cut[j-1]),
                    'max_px':None if j==3 else float(cut[j]),
                    'driest_all_photos':len(dry),'wettest_all_photos':len(wet),
                    'driest_classifiable':n(dry),'wettest_classifiable':n(wet),
                    'driest_photo_white':float(dry.loc[dry.classified,'white'].mean()),
                    'wettest_photo_white':float(wet.loc[wet.classified,'white'].mean()),
                    'wet_minus_dry_photo_white':float(wet.loc[wet.classified,'white'].mean()-dry.loc[dry.classified,'white'].mean()),
                    'driest_classifiable_rate':float(dry.classified.mean()),
                    'wettest_classifiable_rate':float(wet.classified.mean())})
    return out

def run(path:Path):
    s=load_original(path)
    s=s.loc[s.latitude.between(-90,90)&s.longitude.between(-180,180)&s.wc_bio12.notna()].copy()
    if len(s)!=EXPECTED['georain']:raise RuntimeError('geographic rain original photo opportunity drift')
    s['rain_decile']=pd.qcut(s.wc_bio12,10,labels=False,duplicates='drop')
    if s.rain_decile.nunique()!=10:raise RuntimeError('not 10 original opportunity rainfall deciles')
    # All climate features are available wherever original BIO12 exists, but do not assume:
    req=['wc_bio1','wc_bio5','wc_bio12','wc_bio15','wc_elevation_m']
    complete=s[req].notna().all(axis=1)
    if not complete.all():raise RuntimeError('source climate/photo opportunity lacks predictors')
    s=s.reset_index(drop=True)
    counts=s.groupby('rain_decile').agg(n=('photo_id','size'),n_class=('classified','sum')).reset_index()
    obs,uncorr=fractions(s,np.ones(len(s)))
    results={}
    for model in ('geo_climate','geo_climate_plus_flower_area'):
        pred,fit=crossfit(s,model)
        p=np.clip(pred,MIN_PROPENSITY,MAX_PROPENSITY)
        w=1/p
        fractions_by_bin,gain=fractions(s,w)
        results[model]={'propensity_oof':fit,'original_classified_photo_weighting':'inverse propensity from original all classified and unclassified, clipped to .05-.95',
             'rain_extremes_ipw_fractions':fractions_by_bin,
             'IPW_wet_minus_dry_photo_white_fraction':gain,
             'genus_cluster_95CI':cluster_bootstrap_original(s,w,'genus'),
             'true_site_10degree_geo_cluster_95CI':cluster_bootstrap_original(s,w,'site_cell_10degree')}
    result={'schema':'fcp_original100543_colour_rainfall_photo_selection_mnar_sensitivity_v1',
         'status':'SOURCE_ORIGINAL_PHOTO_PROPENSITY_CROSSFIT_AND_MNAR_TIPPING_COMPLETED',
         'date_jst':'2026-10-10',
         'source_original_archive_sha256':SOURCE_SHA,
         'all_original_unique_photo_ids':EXPECTED['photo_union'],
         'original_geolocated_rain_photo_opportunities':len(s),
         'source_all_classifiable':EXPECTED['classifiable'],
         'source_all_unclassifiable':EXPECTED['missing'],
         'source_georain_classifiable':int(s.classified.sum()),
         'original_photo_wet_dry_unweighted':obs,
         'original_photo_wet_minus_dry_unweighted':uncorr,
         'rain_bin_all_photo_source_counts':counts.to_dict('records'),
         'modelled_selection_missing_at_random_CROSSED':results,
         'original_photo_flower_region_pixel_count_quartile_diagnostic':quality_strata(s),
         'latent_hue_missing_not_at_random_tipping_curve':tipping(s),
         'no_human_expert_flower_colour_truth_labels':True,
         'no_causal_pigment_tradeoff_cost_benefit_estimated':True,
         'classifiability_correction_identifiability':'IPW assumes missing at random given included original environmental/photograph predictors and original photo-site overlap. MNAR outcome-colour missingness is not identified.',
         'hard_nonclaims':[
           'Using all source photo opportunity rows to model classification does not recover the unknown true colour of unclassified photographs.',
           'Small IPW propensity correction is only a measured-photograph, missing-at-random sensitivity, not proof of missing at random.',
           'Flower-effective-pixels is a post-image-derived quality feature; treating it as a predictor is a diagnostic, not a causal adjustment.',
           'Large differences between source dry and wet unknown photo colour compositions can reverse every observed classified-photo trend.',
           'The original 405-source-photo two-expert botanical organ and hue adjudication has not occurred.',
           'White-versus-nonwhite 4-state photographic classes do not identify UV pigment or fitness tradeoff.',
           'Same photo species classes, original source H1/H2, and unopened 2000+730 future source taxa remain frozen.'
         ]}
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source-archive',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    out=run(a.source_archive)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/'result.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k in ['original_photo_wet_minus_dry_unweighted','source_georain_classifiable','modelled_selection_missing_at_random_CROSSED','latent_hue_missing_not_at_random_tipping_curve','original_photo_flower_region_pixel_count_quartile_diagnostic']},indent=2))
if __name__=='__main__':main()