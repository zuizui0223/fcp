#!/usr/bin/env python3
"""Joint environmental transfer across unseen genera AND raw-photo geography.

All36 measured covariates jointly; no significance-based variable selection.
The 4 fixed longitude sectors are not named biogeographical realms. Training
purges whole species based on ALL original photo footprints, including unknown
colours, not merely averaged or classified source locations. This is an
observational source-exposed generalization test, not pigment cost/fitness.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from sklearn.model_selection import GroupKFold
import fcp_audit_36_preflight_20261011 as src
from fcp_all36_global_region_cv_20261011 import aggregate_species
from fcp_joint_environment_unseen_genera_20261011 import (
    predict, loss_components, verify_file, SHA_PHOTO, SHA_ENV, FULL, CLIMATE, GEO, CLASSES)

EARTH_KM=6371.0088
SEED=2026101151
SCHEMES=('genus_only','sector_purged','buffer500')
COMPONENTS=('total','white_axis','remaining_hue')
SCHEMA='fcp_full_environment_joint_genus_geography_transport_20261011_v1'


def sector_ids(longitude):
    lon=np.asarray(longitude,float)
    if not np.isfinite(lon).all() or (abs(lon)>180).any():
        raise ValueError('invalid source longitude')
    return np.minimum(3,np.floor((lon+180.)/90.).astype(int))


def spatial_exclusions(raw,buffer_km=500.):
    if buffer_km<0:raise ValueError('negative buffer')
    valid=raw.latitude.between(-90,90)&raw.longitude.between(-180,180)
    missing=set(raw.loc[~valid,'inat_taxon_id'].astype(int))
    r=raw.loc[valid].copy()
    r['sector']=sector_ids(r.longitude)
    coords=np.deg2rad(r[['latitude','longitude']].to_numpy(float))
    masks={};stats={}
    for s in range(4):
        inside=r.sector.eq(s).to_numpy()
        if not inside.any():
            masks[s]={'sector_purged':set(missing),'buffer500':set(missing)}
            stats[s]={'n_region_source_photos':0,'status':'NO_SOURCE_TARGET_REGION'}
            continue
        tree=BallTree(coords[inside],metric='haversine')
        km=tree.query(coords,k=1,return_distance=True)[0][:,0]*EARTH_KM
        whole_sector=set(r.loc[inside,'inat_taxon_id'].astype(int))|missing
        whole_buffer=whole_sector|set(r.loc[km<=buffer_km+1e-7,'inat_taxon_id'].astype(int))
        masks[s]={'sector_purged':whole_sector,'buffer500':whole_buffer}
        stats[s]={'n_region_source_photos':int(inside.sum()),
           'n_source_species_with_missing_coordinate':len(missing),
           'n_source_species_excluded_by_region_or_missing':len(whole_sector),
           'n_source_species_excluded_by_region_buffer_or_missing':len(whole_buffer),
           'n_outside_region_source_photos_within_buffer':int(((~inside)&(km<=buffer_km+1e-7)).sum()),
           'minimum_retained_raw_photo_distance_km':float(km[~r.inat_taxon_id.isin(whole_buffer)].min())
               if (~r.inat_taxon_id.isin(whole_buffer)).any() else None,
           'all_original_photo_states_used':True}
    return masks,stats


def select_train(sp,heldout_fold,forbidden):
    return sp.loc[sp.fold.ne(heldout_fold)&~sp.inat_taxon_id.isin(forbidden)].copy()


def environmental_support(train_x,target_x):
    x=np.asarray(train_x,float);y=np.asarray(target_x,float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('nonfinite environment')
    mu=x.mean(0);sd=x.std(0);sd[sd<1e-8]=1
    z=(x-mu)/sd;t=(y-mu)/sd
    cov=.9*(z.T@z/len(z))+.1*np.eye(z.shape[1])
    inv=np.linalg.inv(cov)
    tr=np.einsum('ij,jk,ik->i',z,inv,z)
    te=np.einsum('ij,jk,ik->i',t,inv,t)
    cutoff=float(np.quantile(tr,.99))
    return te<=cutoff,te,cutoff


def summarize_rows(records,reps=1999,seed=SEED):
    if records.empty:return {'status':'HOLD_NO_TARGET_RECORDS'}
    cols=[prefix+'_'+c for prefix in ('base','full') for c in COMPONENTS]
    # Each species receives total weight1, regardless of its number of sectors/photos.
    sp=records.groupby('inat_taxon_id',sort=True).agg({**{c:'mean' for c in cols},'genus':'first'})
    genus=sp.groupby('genus',sort=True)
    gs=genus[cols].sum().to_numpy(float);ct=genus.size().to_numpy(float)
    gm=gs/ct[:,None];ng=len(ct)
    if ng<5:return {'status':'HOLD_FEWER_THAN_5_GENERA','n_species':len(sp),'n_genera':ng}
    rng=np.random.default_rng(seed)
    draw={k:[] for k in ('genus_equal','species_equal')}
    for start in range(0,reps,80):
        ix=rng.integers(ng,size=(min(80,reps-start),ng))
        for key,vals,den in [('genus_equal',gm,None),('species_equal',gs,ct)]:
            totals=vals[ix].mean(1) if den is None else vals[ix].sum(1)/den[ix].sum(1)[:,None]
            delta=totals[:,:3]-totals[:,3:]
            relative=np.divide(100*delta,totals[:,:3],out=np.zeros_like(delta),where=totals[:,:3]>1e-15)
            draw[key].append(np.stack((delta,relative),axis=-1))
    out={'status':'OK','n_species':len(sp),'n_genera':ng,'n_species_sector_records':len(records),
        'n_sectors':int(records.sector.nunique()),'bootstrap_replicates':reps,
        'fraction_genera_with_positive_mean_gain':float(np.mean(gm[:,0]-gm[:,3]>0))}
    for key in draw:
        point=gm.mean(0) if key=='genus_equal' else sp[cols].mean().to_numpy(float)
        values=np.concatenate(draw[key]);est={}
        for j,c in enumerate(COMPONENTS):
            gain=float(point[j]-point[j+3])
            est[c]={'baseline_loss':float(point[j]),'environment_loss':float(point[j+3]),'gain':gain,
                 'gain_CI95':np.quantile(values[:,j,0],[.025,.975]).tolist(),
                 'gain_CI97_5':np.quantile(values[:,j,0],[.0125,.9875]).tolist(),
                 'relative_gain_percent':float(100*gain/point[j]) if point[j]>1e-15 else None,
                 'relative_gain_percent_CI95':np.quantile(values[:,j,1],[.025,.975]).tolist()}
        out[key]=est
    return out


def evaluate_case(raw,mode,quality,exclusions,reps,outdir):
    d=raw.copy()
    if quality:
        d=d.loc[(d.flower_effective_pixels>=23826)&(d.nuisance_pixel_fraction<=0.1375717742195135)].copy()
    d['sector']=sector_ids(d.longitude)
    global_sp=aggregate_species(d).sort_values('inat_taxon_id',kind='stable').reset_index(drop=True)
    mapping={}
    for fold,(_,te) in enumerate(GroupKFold(n_splits=5).split(global_sp,groups=global_sp.genus)):
        for g in global_sp.iloc[te].genus.unique():mapping[g]=fold
    global_sp['fold']=global_sp.genus.map(mapping)
    features={'geography':GEO,'climate26':GEO+CLIMATE}
    if mode=='soil':features['environment36']=GEO+FULL
    envname='environment36' if mode=='soil' else 'climate26'
    actual_env=FULL if mode=='soil' else CLIMATE
    case=mode+('_quality' if quality else '_all')
    rows=[];audits=[]
    for sector in range(4):
        target_sp=aggregate_species(d[d.sector.eq(sector)]).sort_values('inat_taxon_id').reset_index(drop=True)
        target_sp['fold']=target_sp.genus.map(mapping)
        for fold in range(5):
            target=target_sp[target_sp.fold.eq(fold)].copy()
            if not len(target):continue
            truth=target[['hue_'+h for h in CLASSES]].to_numpy(float)
            for scheme in SCHEMES:
                forbidden=set() if scheme=='genus_only' else exclusions[sector][scheme]
                train=select_train(global_sp,fold,forbidden)
                guard={'sector':sector,'fold':fold,'scheme':scheme,'n_train_species':len(train),
                    'n_train_genera':int(train.genus.nunique()),'n_target_species':len(target),
                    'n_target_photos':int(target.n_photos.sum()),
                    'shared_train_test_genus':len(set(train.genus)&set(target.genus)),
                    'shared_train_test_species':len(set(train.inat_taxon_id)&set(target.inat_taxon_id)),
                    'training_species_intersect_spatial_blacklist':len(set(train.inat_taxon_id)&forbidden)}
                if any(guard[k] for k in ('shared_train_test_genus','shared_train_test_species','training_species_intersect_spatial_blacklist')):
                    raise RuntimeError('source leakage '+str(guard))
                if len(train)<300 or train.genus.nunique()<30:
                    guard['status']='HOLD_TRAINING_SUPPORT';audits.append(guard);continue
                y=train[['hue_'+h for h in CLASSES]].to_numpy(float)
                prediction={'constant':np.tile(y.mean(0),(len(target),1))}
                for name,feat in features.items():prediction[name]=predict(train,target,feat,'absolute')
                losses={name:loss_components(truth-p) for name,p in prediction.items()}
                support,ds,cut=environmental_support(train[list(actual_env)],target[list(actual_env)])
                guard.update(status='OK',n_target_within_training_environment_envelope=int(support.sum()),
                             train_environment_d2_99_cutoff=cut)
                audits.append(guard)
                r=target[['inat_taxon_id','genus','n_photos','cell']].copy()
                r['sector']=sector;r['fold']=fold;r['scheme']=scheme;r['environment_support']=support
                r['environment_squared_distance']=ds
                for name,lo in losses.items():
                    for j,c in enumerate(COMPONENTS):r[name+'__'+c]=lo[:,j]
                rows.append(r)
        print(json.dumps({'case':case,'sector_completed':sector,'n_species_sector_targets':len(target_sp)}),flush=True)
    if not rows:raise RuntimeError('no test predictions')
    records=pd.concat(rows,ignore_index=True)
    if records.duplicated(['inat_taxon_id','sector','scheme']).any():raise RuntimeError('duplicate species-sector target')
    records.to_csv(outdir/(case+'_oof_losses.csv.gz'),index=False,compression='gzip')
    pd.DataFrame(audits).to_csv(outdir/(case+'_fold_audit.csv'),index=False)
    contrasts=[(envname+'_vs_geography','geography',envname),
               (envname+'_vs_constant','constant',envname)]
    if mode=='soil':contrasts += [('soil_increment_same_source','climate26','environment36')]
    result={'n_input_photos':len(d),'n_input_species':len(global_sp),'n_input_genera':int(global_sp.genus.nunique()),
            'loss_target':'per-species sector-specific source photo colour proportions, squared composition loss',
            'all_eligible_source_photos_accounted':int(d.groupby('sector').size().sum())==len(d),
            'schemes':{},'fold_status_counts':pd.Series([a['status'] for a in audits]).value_counts().to_dict()}
    for scheme in SCHEMES:
        rr=records.loc[records.scheme.eq(scheme)].copy()
        result['schemes'][scheme]={'n_evaluated_species':int(rr.inat_taxon_id.nunique()),
            'n_evaluated_source_photos':int(rr.n_photos.sum()),'n_species_sector_records':len(rr),
            'n_supported_environment_records':int(rr.environment_support.sum()),'contrasts':{}}
        for cname,a,b in contrasts:
            rename={a+'__'+k:'base_'+k for k in COMPONENTS}|{b+'__'+k:'full_'+k for k in COMPONENTS}
            x=rr[['inat_taxon_id','genus','sector','environment_support',*rename]].rename(columns=rename)
            val={'overall':summarize_rows(x,reps,SEED+sum(map(ord,case+scheme+cname))),
                 'regions':{},'within_train_environment_envelope':summarize_rows(x.loc[x.environment_support],reps,
                     SEED+1+sum(map(ord,case+scheme+cname)))}
            for s in range(4):val['regions'][str(s)]=summarize_rows(x.loc[x.sector.eq(s)],reps,
                SEED+11+s+sum(map(ord,case+scheme+cname)))
            result['schemes'][scheme]['contrasts'][cname]=val
        z=result['schemes'][scheme]['contrasts'][envname+'_vs_geography']['overall']
        print(json.dumps({'case':case,'scheme':scheme,'n_species':z.get('n_species'),
            'total':z.get('genus_equal',{}).get('total')}),flush=True)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path,default=Path('/mnt/data'))
    p.add_argument('--outdir',required=True,type=Path);p.add_argument('--bootstrap',type=int,default=1999)
    p.add_argument('--mode',choices=('both','soil','climate'),default='both')
    p.add_argument('--quality',choices=('both','all','strict'),default='both');a=p.parse_args()
    if a.bootstrap<99:raise ValueError('bootstrap below99')
    src.A=a.source_root/'fcp_original_photo_site_climate_soil_85337_20261010.zip'
    src.B=a.source_root/'fcp_100543_all36_expanded_photo_site_environment_20261011.zip'
    verify_file(src.A,SHA_PHOTO);verify_file(src.B,SHA_ENV)
    d,den,masks=src.run()
    exclusion,geo_audit=spatial_exclusions(d,buffer_km=500)
    out={'schema':SCHEMA,'source_counts':den,'source_hashes':{'photo':SHA_PHOTO,'environment':SHA_ENV},
         'source_commit_parent':'6571ce101d3a414a319ff49c103b75ea81ccda1d',
         'all_environment_features':list(FULL),'n_bootstrap':a.bootstrap,'buffer_km':500,
         'longitude_sectors':[[-180,-90],[-90,0],[0,90],[90,180]],
         'source_raw_footprint_exclusion_audit':geo_audit,'cases':{},
         'new_colour_measurements':False,'main_H1H2_changed':False,
         'inference_status':'POSTHOC_OBSERVATIONAL_TAXON_AND_GEOGRAPHY_TRANSPORT',
         'limitations':['Four longitude sectors are artificial spatial domains, not evolutionary replicates.',
          'Named genus separation is not direct dated phylogeny or family independence.',
          'Unknown original photo colours are missing; measurement errors/MNAR remain.',
          'All original photo footprints are used for purging, but unseen real populations are unobserved.',
          'Bootstrap conditions on trained models; only four geographic domains limit geographic generalization.',
          'Environmental covariance support is a training-only ellipsoidal diagnostic, not complete overlap proof.',
          'Test target is sector-specific photographic colour composition, not individual fitness or genotype.',
          'Failure may reflect extrapolation or model misspecification, not biological absence of environmental effect.'],
         'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}}
    a.outdir.mkdir(parents=True,exist_ok=True)
    for mode in (('soil','climate') if a.mode=='both' else (a.mode,)):
        for quality in ((False,True) if a.quality=='both' else (a.quality=='strict',)):
            name=mode+('_quality' if quality else '_all')
            out['cases'][name]=evaluate_case(d.loc[masks[mode]],mode,quality,exclusion,a.bootstrap,a.outdir)
            (a.outdir/'result.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print('SOURCE_REPLAY_COMPLETE',flush=True)
if __name__=='__main__':main()
