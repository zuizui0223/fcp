#!/usr/bin/env python3
"""All measured 36 abiotic predictors, seven blocks, three geographic scopes.

Observational original FCP photograph-only, classifiable outcome species means.
5-fold out-of-geo-cell heldout genus-intercept ridge multiclass Brier scores.
Global fit and original longitude E90-180 / complement; source species equal
and genus-equal evaluation and genus-cluster intervals. No independent true
phylogeny, source photo quality or causal genotype/fitness inference.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from fcp_audit_36_preflight_20261011 import run as source_run,FEATURES,BLOCKS,CLASSES

SEED=20261011
ALPHA=40.0
N_FOLDS=5
N_BOOT=999
GEO=('abs_lat','lon_sin','lon_cos')
MIN_TRAIN_GENUS_SPECIES=3


def cell_index(latitude,longitude):
    lat=np.asarray(latitude,float);lon=np.asarray(longitude,float)
    y=np.minimum(8,np.maximum(0,np.floor((np.sin(np.deg2rad(lat))+1)*4.5).astype(int)))
    x=np.minimum(17,np.maximum(0,np.floor((lon+180)/20).astype(int)))
    return 18*y+x


def aggregate_species(d):
    d=d.copy()
    for c in CLASSES:
        d['hue_'+c]=(d.morph==c).astype(float)
    species=d.groupby('inat_taxon_id',sort=False).agg(
        {**{x:'mean' for x in [*GEO,*FEATURES,*['hue_'+c for c in CLASSES]]},
        'latitude':'mean','longitude':'mean', 'genus':'first','photo_id':'size'}).reset_index()
    species['cell']=cell_index(species.latitude.to_numpy(float),species.longitude.to_numpy(float))
    species['n_photos']=species.photo_id
    species.drop(columns='photo_id',inplace=True)
    assert species.inat_taxon_id.is_unique
    return species


def construct_models(soil):
    blocks={k:v for k,v in BLOCKS.items() if k!='soil' or soil}
    allf=tuple(x for bb in blocks.values() for x in bb)
    m={'GEOGRAPHY_ONLY':GEO,'GEO_ALL_ABIOTIC':GEO+allf}
    for k,b in blocks.items():
        m['GEO_PLUS_'+k.upper()]=GEO+b
        m['FULL_MINUS_'+k.upper()]=GEO+tuple(x for x in allf if x not in b)
    for feature in allf:
        m['FULL_MINUS_SINGLE_'+feature.upper()]=GEO+tuple(x for x in allf if x!=feature)
    return m,blocks


def fit_predict(train,test,cols,train_genus_n=None):
    # TRAIN species only determines both genus intercept and centered slope.
    trgen=train.genus.to_numpy(str)
    tegen=test.genus.to_numpy(str)
    yt=train[['hue_'+c for c in CLASSES]].to_numpy(float)
    cy=['hue_'+c for c in CLASSES]
    ym=pd.DataFrame(yt,columns=cy).assign(genus=trgen).groupby('genus',sort=False)[cy].mean()
    test_genus_y=ym.loc[tegen].to_numpy(float)
    xx=train[list(cols)].to_numpy(float);xt=test[list(cols)].to_numpy(float)
    mu=xx.mean(axis=0);sd=xx.std(axis=0)
    sd=np.where(sd<1e-8,1.,sd)
    xx=(xx-mu)/sd;xt=(xt-mu)/sd
    xs=pd.DataFrame(xx,columns=cols).assign(genus=trgen).groupby('genus',sort=False)[list(cols)].mean()
    group_tr_x=xs.loc[trgen].to_numpy(float)
    group_tr_y=ym.loc[trgen].to_numpy(float)
    beta=np.linalg.solve((xx-group_tr_x).T@(xx-group_tr_x)+ALPHA*np.eye(len(cols)),
                          (xx-group_tr_x).T@(yt-group_tr_y))
    score=test_genus_y+(xt-xs.loc[tegen].to_numpy(float))@beta
    # All models use identical simplex projection to get comparably bounded Brier scores.
    score=np.maximum(0.,score);norm=score.sum(axis=1,keepdims=True)
    score=np.divide(score,norm,out=np.full_like(score,1/len(CLASSES)),where=norm>0)
    return score


def point_and_ci(delta,group,rep=N_BOOT,seed=SEED):
    keys,gid=np.unique(np.asarray(group,str),return_inverse=True)
    n=len(keys)
    vs=np.bincount(gid,weights=delta,minlength=n)
    cs=np.bincount(gid,minlength=n)
    rng=np.random.default_rng(seed)
    ix=rng.integers(0,n,size=(rep,n))
    draws=vs[ix].sum(axis=1)/cs[ix].sum(axis=1)
    pergenus=vs/cs
    draws_equal=pergenus[ix].mean(axis=1)
    return {'photo_species_equal_gain':float(delta.mean()),
        'source_genus_cluster_CI95':[float(x) for x in np.quantile(draws,[.025,.975])],
        'genus_equal_gain':float(pergenus.mean()),
        'genus_equal_cluster_CI95':[float(x) for x in np.quantile(draws_equal,[.025,.975])],
        'n_genus_clusters':int(n)}


def evaluate(species,soil,sector,nboot=N_BOOT):
    models,blocks=construct_models(soil)
    categories=tuple('hue_'+c for c in CLASSES)
    n=len(species)
    if n<400 or species.cell.nunique()<N_FOLDS:
        return {'status':'HOLD_INSUFFICIENT_SPECIES_OR_SPATIAL_CELLS','n_species':n}
    yy=species[list(categories)].to_numpy(float)
    pred={k:np.full((n,4),np.nan,float) for k in models}
    valid=np.zeros(n,bool)
    folds=[]
    for fold,(tr,te) in enumerate(GroupKFold(n_splits=N_FOLDS).split(species,groups=species.cell)):
        train=species.iloc[tr];test=species.iloc[te]
        genera=train.genus.value_counts()
        ix=np.flatnonzero(test.genus.isin(genera[genera>=MIN_TRAIN_GENUS_SPECIES].index).to_numpy())
        target=test.iloc[ix]
        if set(train.inat_taxon_id)&set(target.inat_taxon_id):
            raise RuntimeError('Source species crossed groupfold boundary')
        valid[te[ix]]=True
        folds.append({'fold':fold,'n_test':len(te),'n_supported_genus_test_species':len(target),
                      'n_geographic_cells_heldout':test.cell.nunique()})
        if not len(target):continue
        for key,cols in models.items():
            pred[key][te[ix]]=fit_predict(train,target,cols)
    eligible=species.loc[valid]
    if valid.sum()<200:raise RuntimeError('No evaluation source species in original geo holdout')
    if any(~np.isfinite(a[valid]).all() for a in pred.values()):
        raise RuntimeError('Incomplete original photo comparison population')
    loss={name:np.square(values[valid]-yy[valid]).sum(axis=1) for name,values in pred.items()}
    out={}
    for i,block in enumerate(blocks):
        gains={}
        for label,src,tgt in [('marginal_vs_geography','GEOGRAPHY_ONLY','GEO_PLUS_'+block.upper()),
                             ('unique_beyond_all_others','FULL_MINUS_'+block.upper(),'GEO_ALL_ABIOTIC')]:
            delta=loss[src]-loss[tgt]
            gains[label]=point_and_ci(delta,eligible.genus.to_numpy(str),nboot,seed=SEED+i*97+(0 if label.startswith('marginal') else 23)+len(sector)*17)
            gains[label]['positive_species_and_genus_equal_intervals']=bool(gains[label]['source_genus_cluster_CI95'][0]>0 and gains[label]['genus_equal_cluster_CI95'][0]>0)
        out[block]=gains
    single={}
    for j,f in enumerate(x for blockvars in blocks.values() for x in blockvars):
        d=loss['FULL_MINUS_SINGLE_'+f.upper()]-loss['GEO_ALL_ABIOTIC']
        single[f]=point_and_ci(d,eligible.genus,nboot,seed=SEED+3043+j*13+len(sector)*7)
    delta=loss['GEOGRAPHY_ONLY']-loss['GEO_ALL_ABIOTIC']
    return {'status':'EXPLORATORY_ORIGINAL_ALL_BLOCKS_GEOCELL_OUT_OF_FOLD_GENUS_BASELINE',
      'n_source_photo_species':len(species),
      'n_source_photos':int(species.n_photos.sum()),
      'n_source_test_photo_species':int(valid.sum()),
      'n_test_source_genera':int(eligible.genus.nunique()),
      'n_heldout_original_geo_cells':int(eligible.cell.nunique()),
      'folds':folds,'n_orig_abiotic_variables':sum(map(len,blocks.values())),
      'models_original_photo_colour_response':'species 4 hue PHOTO-average categorical fraction, one equal-weight unit per nominal source taxon',
      'geography_controls':list(GEO),
      'ridge_regularization_alpha':ALPHA,
      'outcome_unobserved_species_never_classified_as_white':True,
      'baseline_heldout_brier':float(loss['GEOGRAPHY_ONLY'].mean()),
      'full_heldout_brier':float(loss['GEO_ALL_ABIOTIC'].mean()),
      'all_environment_beyond_geography':point_and_ci(delta,eligible.genus,nboot,seed=SEED+1),
      'all_seven_blocks':out,
      'all_measured_single_predictor_ablation':single,
      'all_source_feature_hypotheses_are_posthoc_and_multiplicity_exposed':True,
      'nonclaims':['original photo colours are not verified botanical organ/genetic morphs',
        'genus-conditioned heldout prediction is not direct dated phylogeny',
        'geographic source-cell holdout does not fully remove true photo-site spatial confounding',
        'retrospective 7 blocks x regions x missingness cases require multiplicity and independent validation',
        'models do not measure actual anthocyanin production, pollination/defense benefit or reproductive cost']}


def run(nboot=N_BOOT,only_region=None,only_mode=None):
    d,counts,pools=source_run()
    r={'schema':'fcp_original100543_36variables_seven_abiotic_blocks_regional_species_genus_predictive_audit_v1',
       'date_jst':'2026-10-11',
       'source_archive_n_original_unique_photo_ids':100543,
       'source_original_species_opportunity':42111,
       'source_original_classifiable':45953,
       'source_original_photo_unclassifiable':54590,
       'variable_blocks':{k:list(v) for k,v in BLOCKS.items()},
       'source_full_abiotic_coverage':counts,
       'n_source_genus_bootstrap_replications':nboot,
       'region_definition':'east 90 <= true photo lon <= 180; rest other original photo longitudes; global all original sites',
       'results':{},'confirmatory_previous_H1H2_unchanged':True}
    for mode,mask in pools.items():
        if only_mode and mode!=only_mode: continue
        src=d.loc[mask].copy()
        r['results'][mode]={}
        for name,frame in [('world',src),('east',src[src.east]),('elsewhere',src[~src.east])]:
            if only_region and name!=only_region: continue
            print('starting',mode,name,len(frame),flush=True)
            sp=aggregate_species(frame)
            z=evaluate(sp,soil=mode=='soil',sector=name,nboot=nboot)
            r['results'][mode][name]=z
            print('ended',mode,name,z['n_source_test_photo_species'],
                   {k:round(v['unique_beyond_all_others']['genus_equal_gain'],5) for k,v in z['all_seven_blocks'].items()},flush=True)
    return r

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=Path('/mnt/data/fcp_all36_source_results_20261011'))
    p.add_argument('--bootstrap',type=int,default=N_BOOT)
    p.add_argument('--only-region',choices=['world','east','elsewhere'],default=None)
    p.add_argument('--only-mode',choices=['climate','soil'],default=None)
    a=p.parse_args();res=run(a.bootstrap,only_region=a.only_region,only_mode=a.only_mode)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/'result.json').write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n')
    print('DONE',a.outdir,flush=True)