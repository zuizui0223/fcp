#!/usr/bin/env python3
"""Unadjusted descriptive white-minus-nonwhite climate separation in five photo cohorts.

No p-values, causal claims, observer/background adjustment, or individual-photo
outputs. Climate surfaces are long-term normals, not weather on observation dates.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import zipfile

import numpy as np
import pandas as pd
from PIL import Image, __version__ as pillow_version

COHORTS=('discovery','validation','prospective','wave4','new_U100')
MORPHS=('white','yellow_orange','red_pink','blue_purple')
VARIABLES=('BIO1','BIO4','BIO5','BIO12','BIO14','BIO15','SRADmean')
PHOTO_SHA='8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6'
ZIP_SHA={'wc2.1_10m_bio.zip':'00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf',
         'wc2.1_10m_srad.zip':'c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026'}
DEFINITIONS={
 'BIO1':{'name':'Annual Mean Temperature','native_unit':'degree Celsius'},
 'BIO4':{'name':'Temperature Seasonality','native_unit':'temperature standard deviation multiplied by 100; native raster index'},
 'BIO5':{'name':'Max Temperature of Warmest Month','native_unit':'degree Celsius'},
 'BIO12':{'name':'Annual Precipitation','native_unit':'mm'},
 'BIO14':{'name':'Precipitation of Driest Month','native_unit':'mm'},
 'BIO15':{'name':'Precipitation Seasonality','native_unit':'coefficient of variation; native WorldClim index'},
 'SRADmean':{'name':'Arithmetic mean of 12 monthly solar-radiation normals','native_unit':'kJ m^-2 day^-1'},
}


def demand(ok,message):
    if not ok: raise ValueError('FCP_ENVIRONMENT_DESCRIPTION_INVALID: '+message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def validate_tags(tags,size):
    demand(size==(2160,1080),'unexpected raster dimensions')
    demand(np.allclose(tags.get(33550,()),(1/6,1/6,0),rtol=0,atol=1e-14),'pixel scale changed')
    demand(np.allclose(tags.get(33922,()),(0,0,0,-180,90,0),rtol=0,atol=1e-14),'tiepoint changed')
    keys=tuple(tags.get(34735,()))
    demand(len(keys)>=4 and len(keys)==4+4*keys[3],'GeoKeyDirectory invalid')
    lookup={keys[i]:(keys[i+1],keys[i+2],keys[i+3]) for i in range(4,len(keys),4)}
    demand(lookup.get(1024)==(0,1,2) and lookup.get(1025)==(0,1,1)
           and lookup.get(2048)==(0,1,4326),'EPSG4326 geographic PixelIsArea required')
    demand(42113 in tags,'nodata absent')


def sample_array(array,lat,lon,nodata,*,step=1/6):
    """Top-left area-cell floor sampling; out-of-grid coordinates remain missing."""
    lat,lon=np.asarray(lat,float),np.asarray(lon,float)
    good=np.isfinite(lat)&np.isfinite(lon)
    rr=np.full(len(lat),-1,dtype=np.int64); cc=rr.copy()
    rr[good]=np.floor((90-lat[good])*(1/step)).astype(np.int64)
    cc[good]=np.floor((lon[good]+180)*(1/step)).astype(np.int64)
    good &= (rr>=0)&(rr<array.shape[0])&(cc>=0)&(cc<array.shape[1])
    out=np.full(len(lat),np.nan)
    nd=np.asarray(nodata,dtype=array.dtype).item()
    values=array[rr[good],cc[good]].astype(float)
    values[(values==nd)|~np.isfinite(values)]=np.nan
    out[good]=values
    return out


def verify_raster_member(archive,path):
    """Link the bytes sampled to the already SHA-authenticated climate ZIP."""
    h=hashlib.sha256()
    with zipfile.ZipFile(archive) as z, z.open(path.name) as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    demand(sha(path)==h.hexdigest(),path.name+': extracted raster differs from archive member')


def read_and_sample(path,lat,lon,archive):
    verify_raster_member(archive,path)
    with Image.open(path) as im:
        tags=im.tag_v2
        validate_tags(tags,im.size)
        nodata=float(str(tags[42113]).strip('\x00'))
        array=np.array(im)
        metadata={'file':path.name,'sha256':sha(path),'width':im.width,'height':im.height,
                  'dtype':str(array.dtype),'nodata':nodata,'epsg':4326,'pixel_is_area':True,
                  'pixel_scale':list(tags[33550]),'tiepoint':list(tags[33922]),
                  'authenticated_archive_member_match':True,'archive':archive.name}
    return sample_array(array,lat,lon,nodata),metadata


def effect(x,white):
    x=np.asarray(x,float); white=np.asarray(white,bool)
    demand(np.isfinite(x).all() and white.any() and (~white).any(),'effect sample invalid')
    raw=float(x[white].mean()-x[~white].mean())
    sd=float(x.std(ddof=0))
    threshold=32*np.finfo(float).eps*max(1.,float(np.abs(x).max()))
    state='zero' if np.ptp(x)==0 else ('near_zero' if sd<=threshold else 'estimable')
    return raw,sd,(raw/sd if state=='estimable' else np.nan),state


def species_results(photo,env,*,minimum_each):
    demand(env.shape==(len(photo),len(VARIABLES)),'environment matrix shape differs')
    use=np.isfinite(env).all(axis=1)
    complete=photo.loc[use].reset_index(drop=True)
    e=env[use]
    records=[]; panels=[]
    for tid,ix in complete.groupby('inat_taxon_id',sort=True).indices.items():
        g=complete.iloc[ix]; white=g.morph.eq('white').to_numpy()
        nw=int(white.sum()); nc=int((~white).sum())
        if len(g)<40 or min(nw,nc)<minimum_each: continue
        base={'inat_taxon_id':int(tid),'species':str(g.species.iloc[0]),
              'source_cohort':str(g.fcp_source_cohort.iloc[0]),
              'n_complete_classifiable':len(g),'n_white':nw,'n_nonwhite':nc}
        panels.append(base)
        for j,variable in enumerate(VARIABLES):
            raw,sd,standardized,state=effect(e[ix,j],white)
            records.append({**base,'variable':variable,'white_mean':float(e[ix,j][white].mean()),
                            'nonwhite_mean':float(e[ix,j][~white].mean()),'raw_difference':raw,
                            'within_species_sd_ddof0':sd,'standardized_difference':standardized,
                            'sd_status':state,'native_unit':DEFINITIONS[variable]['native_unit']})
    return pd.DataFrame(records),pd.DataFrame(panels)


def summarize(rows,panel):
    summaries=[]
    for group in (*COHORTS,'pooled'):
        sub=rows if group=='pooled' else rows.loc[rows.source_cohort.eq(group)]
        for variable,v in sub.groupby('variable',sort=False):
            for metric in ('raw_difference','standardized_difference'):
                x=v[metric].to_numpy(float); x=x[np.isfinite(x)]
                row={'panel':panel,'group':group,'variable':variable,'metric':metric,
                     'n_panel_species':len(v),'n_estimable':len(x),'n_missing':len(v)-len(x),
                     'positive_count':int((x>0).sum()),'negative_count':int((x<0).sum()),
                     'zero_count':int((x==0).sum()),'species_equal_weight':True}
                for key,value in (('mean',np.mean),('median',np.median),
                                  ('q25',lambda a:np.quantile(a,.25)),('q75',lambda a:np.quantile(a,.75)),
                                  ('positive_fraction',lambda a:np.mean(a>0)),
                                  ('zero_fraction',lambda a:np.mean(a==0))):
                    row[key]=float(value(x)) if len(x) else np.nan
                row['iqr']=row['q75']-row['q25']
                summaries.append(row)
    return pd.DataFrame(summaries)


def json_safe(value):
    if isinstance(value,dict): return {k:json_safe(v) for k,v in value.items()}
    if isinstance(value,list): return [json_safe(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value): return None
    if isinstance(value,np.generic): return value.item()
    return value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--photos',type=Path,required=True)
    p.add_argument('--worldclim-root',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    started=datetime.now(timezone.utc).isoformat()
    demand(sha(a.photos)==PHOTO_SHA,'original measured-photo SHA mismatch')
    for name,digest in ZIP_SHA.items(): demand(sha(a.worldclim_root/name)==digest,name+' SHA mismatch')
    cols=['inat_taxon_id','species','photo_id','fcp_source_cohort','morph','global_classifiable','latitude','longitude']
    full=pd.read_csv(a.photos,usecols=cols,low_memory=False)
    demand(len(full)==447500 and full.inat_taxon_id.nunique()==4475 and full.photo_id.is_unique,'source identities differ')
    demand(set(full.fcp_source_cohort)==set(COHORTS)
           and (full.groupby('inat_taxon_id').fcp_source_cohort.nunique()==1).all(),'cohort partition differs')
    classified=full.global_classifiable.fillna(False).astype(str).str.lower().isin(['true','1','yes','y'])
    demand(int(classified.sum())==227107 and full.loc[classified].morph.isin(MORPHS).all(),'classified source differs')
    photo=full.loc[classified].reset_index(drop=True)
    lat=photo.latitude.to_numpy(float); lon=photo.longitude.to_numpy(float)
    sampled=[]; raster_meta=[]
    for number in (1,4,5,12,14,15):
        values,meta=read_and_sample(a.worldclim_root/'bio'/f'wc2.1_10m_bio_{number}.tif',lat,lon,
                                    a.worldclim_root/'wc2.1_10m_bio.zip')
        sampled.append(values); raster_meta.append(meta)
    monthly=[]
    for month in range(1,13):
        values,meta=read_and_sample(a.worldclim_root/'srad'/f'wc2.1_10m_srad_{month:02d}.tif',lat,lon,
                                    a.worldclim_root/'wc2.1_10m_srad.zip')
        monthly.append(values); raster_meta.append(meta)
    monthly=np.stack(monthly,axis=1)
    # np.mean, not nanmean: every monthly raster is required for SRADmean.
    sampled.append(monthly.mean(axis=1)); env=np.stack(sampled,axis=1)
    no_raster=np.zeros_like(env)
    raw_rows,raw_panel=species_results(photo,no_raster,minimum_each=5)
    demand(len(raw_panel)==1574,'raw white/nonwhite candidate denominator changed')
    cc=np.isfinite(env).all(axis=1)
    profile={'schema':'fcp_real_447500_unadjusted_environment_description_v1',
             'status':'complete_post_hoc_descriptive_no_pvalues','source_photo_rows':len(full),
             'source_measured_species':4475,'classifiable_photo_rows':len(photo),
             'complete_case_classifiable_photo_rows':int(cc.sum()),
             'missing_any_raster_classifiable_photo_rows':int((~cc).sum()),
             'raw_ge40_white_ge5_nonwhite_ge5_species':len(raw_panel),
             'variable_missing_photos':{v:int((~np.isfinite(env[:,j])).sum()) for j,v in enumerate(VARIABLES)},
             'definitions':DEFINITIONS,'climate_period':'1970-2000 long-term normals',
             'source_documentation':['https://www.worldclim.org/data/bioclim.html',
                                     'https://www.worldclim.org/data/worldclim21.html'],
             'method':{'formula':'(mean white - mean nonwhite) / SD(all complete-case classifiable photos within species)',
                       'sd_ddof':0,'sd_near_zero_threshold':'32 * float64 epsilon * max(1, max(abs(x)))',
                       'same_seven_variable_complete_cases':True,'minimum_complete_classifiable':40,
                       'primary_minimum_each':5,'sensitivity_minimum_each':10,'species_equal_weights':True,
                       'sampling':'floor PixelIsArea; no interpolation, no clipping; invalid and nodata -> missing',
                       'srad_mean_requires_all_12_months':True,'p_values_computed':False},
             'limitations':['Unadjusted photographed white versus nonwhite geographic climate separation, not causal effects',
                'No complete five-cohort near_clip, background or observer controls',
                'Positive-sign fractions are not fractions of true biological responding species',
                'One within-species SD is not the same physical environmental change across species',
                'Long-term climate normals are not weather at the photograph date',
                'Spatial sampling, cultivation, classification and geographic observer biases remain',
                'No genetic, pigment biochemical, fitness or balancing-selection inference'], 'panels':{}}
    outputs={}; all_summaries=[]
    for name,minimum in (('primary_ge5',5),('sensitivity_ge10',10)):
        rows,panel=species_results(photo,env,minimum_each=minimum)
        summary=summarize(rows,name); all_summaries.append(summary)
        profile['panels'][name]={'minimum_each':minimum,'n_species':len(panel),
            'n_complete_classifiable_photos':int(panel.n_complete_classifiable.sum()),
            'by_cohort_species':panel.groupby('source_cohort').size().to_dict(),
            'sd_status_counts':rows.groupby(['variable','sd_status']).size().unstack(fill_value=0).to_dict(orient='index')}
        outputs[f'{name}_species_environment.csv']=rows
        print('FCP_ENV_DESCRIPTION_PANEL',name,'species',len(panel),'photos',int(panel.n_complete_classifiable.sum()),flush=True)
        print(summary.loc[(summary.group=='pooled') & (summary.metric=='standardized_difference'),
                          ['variable','n_estimable','median','mean','positive_fraction','zero_fraction']].to_string(index=False),flush=True)
    summary=pd.concat(all_summaries,ignore_index=True)
    outputs['group_environment_summary.csv']=summary
    profile['summary']=summary.to_dict(orient='records')
    a.outdir.mkdir(parents=True,exist_ok=True)
    for name,frame in outputs.items(): frame.to_csv(a.outdir/name,index=False)
    (a.outdir/'profile.json').write_text(json.dumps(json_safe(profile),indent=2,sort_keys=True,allow_nan=False)+'\n')
    receipt={'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),
             'command':sys.argv,'python':sys.version,'platform':platform.platform(),
             'packages':{'numpy':np.__version__,'pandas':pd.__version__,'Pillow':pillow_version},
             'input_photo_sha256':PHOTO_SHA,'worldclim_zip_sha256':ZIP_SHA,'rasters':raster_meta,
             'code_sha256':sha(Path(__file__)),
             'output_sha256':{name:sha(a.outdir/name) for name in (*outputs,'profile.json')}}
    (a.outdir/'execution_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('FCP_447500_UNADJUSTED_ENVIRONMENT_DESCRIPTION_COMPLETE',flush=True)


if __name__=='__main__': main()
