"""Post-outcome multi-climate screen on frozen third-cohort flower photos.

Uses the SAME photo, classification and technical exposure gate as the original
2026-09-25 white-environment model. Tests additional climate variables but
never rewrites or extends the original three-variable prospective hypothesis.
New 10-variable results and observer/site sensitivities are EXPLORATORY.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import binomtest, spearmanr, wilcoxon

BIO_KEYS = (1, 4, 5, 6, 7, 12, 14, 15, 17)
PREDICTORS = tuple(f"bio{i}" for i in BIO_KEYS) + ("srad_mean",)
FAMILY = {"bio1":"temperature", "bio4":"temperature", "bio5":"temperature", "bio6":"temperature", "bio7":"temperature", "bio12":"moisture", "bio14":"moisture", "bio15":"moisture", "bio17":"moisture", "srad_mean":"radiation"}
MORPHS = {"white", "yellow_orange", "red_pink", "blue_purple"}
MEASURED_SHA256 = "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186"
BATCH_SEED = 2026101401


def file_sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna('').astype(str).str.lower().str.strip().isin(('true','1','yes','y'))


def raster_values(path: Path, lon: pd.Series, lat: pd.Series) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Missing source raster: {path}")
    pts=list(zip(lon.to_numpy(float), lat.to_numpy(float)))
    with rasterio.open(path) as src:
        a=np.fromiter((float(p[0]) for p in src.sample(pts)),dtype=float,count=len(pts))
        if src.nodata is not None:
            a[np.isclose(a,src.nodata)]=np.nan
    return a


def within_z(s: pd.Series, groups: pd.Series) -> pd.Series:
    g=s.groupby(groups)
    mu=g.transform('mean')
    sd=g.transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (s-mu)/sd


def holm_adjust(p: np.ndarray) -> np.ndarray:
    p=np.where(np.isfinite(p),p,1.)
    order=np.argsort(p,kind='stable')
    out=np.empty_like(p,dtype=float)
    out[order]=np.minimum(1.,np.maximum.accumulate((len(p)-np.arange(len(p)))*p[order]))
    return out


def make_panel(measured: Path, technical: Path, join_key: Path, high_ids: Path,
               bio_dir: Path, srad_dir: Path, *, verify_sha:bool=True) -> tuple[pd.DataFrame,dict]:
    if verify_sha and file_sha(measured)!=MEASURED_SHA256:
        raise ValueError('immutable original third-cohort photo SHA256 mismatch')
    need=['photo_id','species','observer_id','morph','global_classifiable','latitude','longitude']
    data=pd.read_csv(measured,usecols=need,low_memory=False)
    if data.photo_id.isna().any() or data.photo_id.duplicated().any():
        raise ValueError('duplicate/missing photo_id')
    data=data[bool_series(data.global_classifiable) & data.morph.isin(MORPHS)].copy()
    data['white']=data.morph.eq('white')
    tech=pd.read_csv(technical,compression='gzip',dtype={'measurement_id':str})
    join=pd.read_csv(join_key,dtype={'measurement_id':str})
    high=pd.read_csv(high_ids,dtype={'measurement_id':str})
    if tech.measurement_id.duplicated().any() or join.measurement_id.duplicated().any():
        raise ValueError('ambiguous technical/seal join')
    tech=tech.merge(join,on='measurement_id',how='left',validate='one_to_one')
    if tech.photo_id.isna().any() or tech.photo_id.duplicated().any():
        raise ValueError('ambiguous original photo-to-technical join')
    tech['high_clip']=tech.measurement_id.isin(set(high.measurement_id))
    data=data.merge(tech[['photo_id','near_clip_fraction','high_clip']],on='photo_id',how='left',validate='one_to_one')
    technical_coverage=float(data.near_clip_fraction.notna().mean())
    if technical_coverage<.85:
        raise ValueError('frozen technical availability below 85%')
    if not data.latitude.between(-90,90).all() or not data.longitude.between(-180,180).all():
        raise ValueError('invalid coordinates')
    for k in BIO_KEYS:
        name=f'bio{k}'
        data[name]=raster_values(bio_dir/f'wc2.1_10m_bio_{k}.tif',data.longitude,data.latitude)
    solar=[raster_values(srad_dir/f'wc2.1_10m_srad_{k:02}.tif',data.longitude,data.latitude) for k in range(1,13)]
    with np.errstate(invalid='ignore'):
        counts=np.sum(np.isfinite(solar),axis=0)
        data['srad_mean']=np.where(counts==0,np.nan,np.nansum(solar,axis=0)/np.where(counts==0,1,counts))
    # Exact chronology of the original 2026-09-25 primary environmental test:
    # standardize each variable across all classifiable photos BEFORE clipping exclusion.
    for k in PREDICTORS:
        data[f'z_{k}']=within_z(data[k],data.species)
    technical_rows=len(data)
    data=data.loc[data.near_clip_fraction.notna() & ~data.high_clip.fillna(False)].copy()
    counts=data.groupby('species').white.agg(['sum','count'])
    eligible=counts.index[(counts['sum']>=5)&((counts['count']-counts['sum'])>=5)]
    data=data.loc[data.species.isin(eligible)].copy()
    metadata={'source_photo_sha256':file_sha(measured), 'technical_coverage':technical_coverage,
              'classifiable_rows_before_quality_exclusion':int(technical_rows),
              'after_high_clip_rows':len(data),'eligible_species':int(data.species.nunique()),
              'variable_family':FAMILY,'climate_variables':PREDICTORS,
              'normalization':'within species across classifiable photos, before exposure exclusion',
              'worldclim_temporal_mean_not_photo_weather':True}
    return data,metadata


def _safe_wilcoxon(a: np.ndarray) -> float:
    a=np.asarray(a,float)
    a=a[np.isfinite(a)]
    if not len(a) or np.all(a==0):
        return 1.
    return float(wilcoxon(a,alternative='two-sided',zero_method='wilcox').pvalue)


def _stat(a: np.ndarray) -> dict:
    a=np.asarray(a,float)
    a=a[np.isfinite(a)]
    if not len(a):return {'n':0,'median':None,'mean':None,'p_wilcoxon_two_sided':None,'positive_fraction':None}
    return {'n':int(len(a)),'median':float(np.median(a)),'mean':float(np.mean(a)),
            'p_wilcoxon_two_sided':_safe_wilcoxon(a), 'positive_fraction':float((a>0).mean())}


def contrasts(panel: pd.DataFrame, *, grid_degrees:float=1.) -> tuple[pd.DataFrame,pd.DataFrame,dict]:
    """Frozen baseline effect plus observer-balanced/paired and local-grid diagnostics."""
    rows=[]
    per_species=[]
    for var in PREDICTORS:
        z='z_'+var
        sub=panel[['species','observer_id','latitude','longitude','white',z,'near_clip_fraction']].dropna(subset=[z]).copy()
        per_species_var=[]; balanced=[]; paired=[]; cells=[]; clip_d=[]
        for sp,g in sub.groupby('species',sort=True):
            yes=g.loc[g.white,z]; no=g.loc[~g.white,z]
            if len(yes)<5 or len(no)<5:continue
            delta=float(yes.mean()-no.mean())
            per_species_var.append(delta)
            clip_delta=float(g.loc[g.white,'near_clip_fraction'].mean()-g.loc[~g.white,'near_clip_fraction'].mean())
            clip_d.append(clip_delta)
            # Equal-weight each observer within colour. Preserve coverage gate separately.
            obs=g.dropna(subset=['observer_id']).groupby(['observer_id','white'],sort=True)[z].mean().unstack()
            bval=np.nan
            if True in obs.columns and False in obs.columns and int(obs[True].notna().sum())>=3 and int(obs[False].notna().sum())>=3:
                bval=float(obs[True].mean()-obs[False].mean())
                balanced.append(bval)
            pval=np.nan
            pn=0
            if True in obs.columns and False in obs.columns:
                one=obs.dropna(subset=[True,False])
                pn=int(len(one))
                if pn:
                    pval=float((one[True]-one[False]).mean())
                    paired.append(pval)
            # Sites in same 1-degree grid block containing both photographic outcomes.
            lats=np.floor((g.latitude.to_numpy(float)+90.)/grid_degrees).astype(int)
            lons=np.floor((g.longitude.to_numpy(float)+180.)/grid_degrees).astype(int)
            grid=g.assign(grid_lat=lats,grid_lon=lons).groupby(['grid_lat','grid_lon','white'],sort=True)[z].mean().unstack()
            kval=0; kdelta=np.nan
            if True in grid.columns and False in grid.columns:
                matched=grid.dropna(subset=[True,False]);kval=int(len(matched))
                if kval>=2:
                    kdelta=float((matched[True]-matched[False]).mean())
                    cells.append(kdelta)
            per_species.append({'species':sp,'variable':var,'white_minus_nonwhite_sd':delta,
                'n_white':int(len(yes)),'n_nonwhite':int(len(no)),
                'observer_balanced_delta':bval,'paired_observer_delta':pval,'n_paired_observers':pn,
                'site_cell_delta':kdelta,'n_mixed_grid_cells':kval,
                'near_clip_white_minus_nonwhite':clip_delta})
        vals=np.asarray(per_species_var,float)
        # Confidence interval for median is descriptive, not a confirmatory discovery gate.
        if len(vals):
            rng=np.random.default_rng(BATCH_SEED+list(PREDICTORS).index(var))
            boot=np.median(vals[rng.integers(0,len(vals),size=(2000,len(vals)))],axis=1)
            ci=[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))]
        else:ci=[None,None]
        result={'variable':var,'family':FAMILY[var],**_stat(vals),
                'median_bootstrap_CI_95':ci,'observer_balanced':_stat(balanced),
                'within_observer_paired':_stat(paired),'within_1degree_mixed_grid_cells':_stat(cells),
                'n_source_photo_species':int(sub.species.nunique())}
        rows.append(result)
    table=pd.DataFrame(per_species)
    results=pd.DataFrame(rows)
    results['holm_10_p']=holm_adjust(results['p_wilcoxon_two_sided'].to_numpy(float))
    results['holm_10_supported']=results['holm_10_p']<.05
    return table,results,{'climate_hypotheses_after_results_exposed':True,
        'confirmatory_evidence':False, 'grid_size_degrees':grid_degrees,
        'bootstrap_n':2000, 'holm_family_n':len(PREDICTORS),
        'observer_cell_controls_change_estimand':True,
        'observer_stratification_does_not_identify_genetic_adaptation':True}


def frozen_gate(per_species: pd.DataFrame, frozen_result: dict, panel: pd.DataFrame,
               *, tolerance:float=1e-10) -> dict:
    if len(panel) < 10000 or panel.species.nunique()<100:
        raise ValueError('original frozen environmental panel lost too many records')
    expected=frozen_result['results']
    if len(expected)!=3:
        raise ValueError('frozen three-variable source has wrong length')
    checks={}
    for d in expected:
        x=per_species.loc[per_species.variable.eq(d['variable']),'white_minus_nonwhite_sd'].to_numpy(float)
        ref=int(d['n_species_delta'])
        if len(x)!=ref:
            raise ValueError(f"{d['variable']}: n={len(x)} differs from original n={ref}")
        summary={'mean':float(x.mean()),'median':float(np.median(x)),'n':len(x)}
        for k,refkey in [('mean','mean_delta_white_minus_nonwhite_SD'),('median','median_delta_white_minus_nonwhite_SD')]:
            if not np.isclose(summary[k],d[refkey],atol=tolerance,rtol=0):
                raise ValueError(f"{d['variable']}: frozen 2026-09-25 source numerical mismatch for {k}")
        checks[d['variable']]=summary
    return checks


def run(measured:Path,technical:Path,join:Path,high:Path,bio:Path,srad:Path,frozen:Path,outdir:Path) -> dict:
    panel,metadata=make_panel(measured,technical,join,high,bio,srad)
    effects,summary,contract=contrasts(panel)
    ref=json.loads(frozen.read_text(encoding='utf-8'))
    checks=frozen_gate(effects,ref,panel)
    if metadata['eligible_species']!=281 or len(panel)!=ref['rows_after_primary_filters']:
        raise ValueError('original reference cohort size mismatch')
    outdir.mkdir(parents=True,exist_ok=True)
    effects.to_csv(outdir/'multi_climate_species_effects.csv',index=False,float_format='%.17g')
    summary.to_json(outdir/'multi_climate_summary.json',orient='records',indent=2)
    metadata['frozen_3var_reproduced']=checks
    metadata['frozen_reference_sha256']=file_sha(frozen)
    metadata['raster_source']='WorldClim v2.1 10-minute BIO and monthly solar, release SHA256 locked by workflow'
    metadata['results_sha256']=file_sha(outdir/'multi_climate_species_effects.csv')
    metadata['status']='EXPLORATORY_10VARIABLE_SCREEN'
    metadata['contract']=contract
    (outdir/'multi_climate_receipt.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return metadata


def main():
    p=argparse.ArgumentParser()
    for arg in ('measured','technical','join_key','high_ids','bio_dir','srad_dir','frozen','outdir'):
        p.add_argument('--'+arg.replace('_','-'),type=Path,required=True)
    a=p.parse_args()
    r=run(a.measured,a.technical,a.join_key,a.high_ids,a.bio_dir,a.srad_dir,a.frozen,a.outdir)
    print(json.dumps({'status':r['status'],'n_eligible':r['eligible_species'],
        'frozen_reproduction':r['frozen_3var_reproduced'],'variables':len(PREDICTORS)}))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
