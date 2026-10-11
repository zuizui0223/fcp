#!/usr/bin/env python3
"""Transfer ALL measured environmental information to genera never used in fitting.

Three noninterchangeable ecological targets: absolute species colour composition,
within-genus differences between species, and within-species photo differences.
The latter two are centred-contrast risks, NOT operational absolute predictions.
No measured pigment costs, genotypes or fitness. Source-exposed observational
analysis; genus holdout is not a dated phylogeny or spatially independent test.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
import fcp_audit_36_preflight_20261011 as src
from fcp_all36_global_region_cv_20261011 import aggregate_species

CLASSES=src.CLASSES
GEO=('abs_lat','lon_sin','lon_cos')
CLIMATE=tuple(v for k,vs in src.BLOCKS.items() if k!='soil' for v in vs)
FULL=tuple(src.FEATURES)
ALPHA=40.
SEED=2026101127
SHA_PHOTO='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
SHA_ENV='6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
LEVELS=('absolute','within_genus','within_species')
COMPONENTS=('total','white_axis','remaining_hue')


def verify_file(path,expected):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    if h.hexdigest()!=expected:raise ValueError('original source SHA256 mismatch: '+str(path))


def center(x,groups):
    x=np.asarray(x,float);_,idx=np.unique(np.asarray(groups),return_inverse=True)
    n=np.bincount(idx)
    means=np.column_stack([np.bincount(idx,weights=x[:,j])/n for j in range(x.shape[1])])
    return x-means[idx]


def make_folds(d):
    return list(GroupKFold(n_splits=5).split(d,groups=d.genus.to_numpy(str)))


def predict(train,test,features,level):
    if set(train.genus)&set(test.genus):raise ValueError('genus leakage into training')
    x=train[list(features)].to_numpy(float); xt=test[list(features)].to_numpy(float)
    y=train[['hue_'+h for h in CLASSES]].to_numpy(float)
    if not np.isfinite(x).all() or not np.isfinite(xt).all() or not np.isfinite(y).all():
        raise ValueError('nonfinite training response/environment')
    if level not in LEVELS:raise ValueError('unknown ecological level')
    w=np.ones(len(train))
    if level=='within_species':
        counts=train.inat_taxon_id.value_counts()
        w=1/counts.reindex(train.inat_taxon_id).to_numpy(float)
    mu=np.average(x,axis=0,weights=w)
    sd=np.sqrt(np.average((x-mu)**2,axis=0,weights=w));sd[sd<1e-8]=1.
    x=(x-mu)/sd;xt=(xt-mu)/sd
    if level=='absolute':
        my=np.average(y,axis=0,weights=w)
        yc=y-my
    else:
        key='genus' if level=='within_genus' else 'inat_taxon_id'
        x=center(x,train[key]);yc=center(y,train[key]);xt=center(xt,test[key]);my=0.
    beta=np.linalg.solve(x.T@(w[:,None]*x)+ALPHA*np.eye(len(features)),x.T@(w[:,None]*yc))
    score=my+xt@beta
    if level=='absolute':
        score=np.maximum(score,0.)
        norm=score.sum(1,keepdims=True)
        score=np.divide(score,norm,out=np.full_like(score,.25),where=norm>0)
    return score


def loss_components(residual):
    r=np.asarray(residual,float)
    if not np.allclose(r.sum(1),0,atol=1e-9):raise ValueError('four-colour contrast must sum to zero')
    total=np.square(r).sum(1);axis=(4./3)*r[:,0]**2
    return np.column_stack((total,axis,total-axis))


def genus_error_partition(residual,groups):
    """ANOVA identity of heldout residual risk; exact at genus/global averages."""
    r=np.asarray(residual,float);within=center(r,groups);mean=r-within
    return np.column_stack(((r*r).sum(1),(mean*mean).sum(1),(within*within).sum(1)))


def summarize(base,full,genera,reps=1999,seed=SEED,components=COMPONENTS):
    """Input rows are species means of loss; resample whole genera in pairs."""
    base=np.asarray(base,float);full=np.asarray(full,float)
    names,idx=np.unique(np.asarray(genera,str),return_inverse=True);ng=len(names)
    counts=np.bincount(idx); n=len(base)
    sums_b=np.column_stack([np.bincount(idx,weights=base[:,j]) for j in range(3)])
    sums_f=np.column_stack([np.bincount(idx,weights=full[:,j]) for j in range(3)])
    means_b=sums_b/counts[:,None];means_f=sums_f/counts[:,None]
    rng=np.random.default_rng(seed)
    draws={k:[] for k in ('genus_equal','species_equal')}
    for start in range(0,reps,60):
        ind=rng.integers(0,ng,size=(min(60,reps-start),ng))
        for scheme,b,f,den in [('genus_equal',means_b,means_f,None),('species_equal',sums_b,sums_f,counts)]:
            if den is None:bb=b[ind].mean(1);ff=f[ind].mean(1)
            else:
                divisor=den[ind].sum(1)[:,None];bb=b[ind].sum(1)/divisor;ff=f[ind].sum(1)/divisor
            gain=bb-ff;rel=100*np.divide(gain,bb,out=np.zeros_like(bb),where=bb>1e-15)
            draws[scheme].append(np.stack((gain,rel),axis=2))
    out={'n_species':n,'n_genera':ng,'bootstrap_repetitions':reps,
         'fraction_genera_positive_mean_total_gain':float(np.mean((means_b-means_f)[:,0]>0)),
         'fraction_species_positive_total_gain':float(np.mean((base-full)[:,0]>0))}
    for scheme in draws:
        if scheme=='genus_equal':b=means_b.mean(0);f=means_f.mean(0)
        else:b=base.mean(0);f=full.mean(0)
        a=np.concatenate(draws[scheme]); vals={}
        for j,key in enumerate(components):
            vals[key]={'baseline_loss':float(b[j]),'environment_loss':float(f[j]),'gain':float(b[j]-f[j]),
                'gain_CI95':np.quantile(a[:,j,0],[.025,.975]).tolist(),
                'gain_CI98_333_three_primary_endpoints':np.quantile(a[:,j,0],[.05/6,1-.05/6]).tolist(),
                'relative_gain_percent':float(100*(b[j]-f[j])/b[j]) if b[j]>1e-15 else None,
                'relative_gain_percent_CI95':np.quantile(a[:,j,1],[.025,.975]).tolist()}
        out[scheme]=vals
    return out


def run_case(photo,mode,quality,reps,outdir):
    photo=photo.copy()
    if quality:
        photo=photo[(photo.flower_effective_pixels>=23826)&(photo.nuisance_pixel_fraction<=0.1375717742195135)].copy()
    for h in CLASSES:photo['hue_'+h]=photo.morph.eq(h).astype(float)
    sp=aggregate_species(photo).sort_values('inat_taxon_id',kind='stable').reset_index(drop=True)
    # Build the genus split ONCE from species counts; all photo-level analyses reuse it.
    fold_map={}
    for fold,(_,te) in enumerate(make_folds(sp)):
        for g in sp.iloc[te].genus.unique():fold_map[g]=fold
    sp['fold']=sp.genus.map(fold_map);photo['fold']=photo.genus.map(fold_map)
    models={'geography':GEO,'climate26':GEO+CLIMATE}
    if mode=='soil':models['environment36']=GEO+FULL
    prefix=mode+('_quality' if quality else '_all')
    outcomes={};by_genus=[]
    for level in LEVELS:
        d=photo.copy().reset_index(drop=True) if level=='within_species' else sp.copy()
        pred={key:np.full((len(d),4),np.nan) for key in models}
        records=[]
        for fold in range(5):
            train=d[d.fold.ne(fold)];test=d[d.fold.eq(fold)]
            if not len(train) or not len(test):raise ValueError('empty genus holdout')
            if set(train.inat_taxon_id)&set(test.inat_taxon_id):raise ValueError('source species leakage')
            for key,feat in models.items():pred[key][test.index]=predict(train,test,feat,level)
            records.append({'fold':fold,'n_train_rows':len(train),'n_test_rows':len(test),
                            'n_test_genera':int(test.genus.nunique()),'shared_train_test_genera':0})
        y=d[['hue_'+h for h in CLASSES]].to_numpy(float)
        if level=='within_genus':y=center(y,d.genus)
        elif level=='within_species':y=center(y,d.inat_taxon_id)
        if any(not np.isfinite(x).all() for x in pred.values()):raise ValueError('source predictions incomplete')
        losses={key:loss_components(y-val) for key,val in pred.items()}
        if level=='within_genus':
            counts=d.genus.value_counts();good=d.genus.isin(counts[counts>=3].index).to_numpy()
        elif level=='within_species':
            counts=d.inat_taxon_id.value_counts();good=d.inat_taxon_id.isin(counts[counts>=2].index).to_numpy()
        else:good=np.ones(len(d),bool)
        ids=d.loc[good,'inat_taxon_id'].to_numpy()
        unique,ix=np.unique(ids,return_inverse=True);sizes=np.bincount(ix)
        species_losses={key:np.column_stack([np.bincount(ix,weights=arr[good,j])/sizes for j in range(3)])
                        for key,arr in losses.items()}
        info=d.loc[good].groupby('inat_taxon_id',sort=True).genus.first()
        if not np.array_equal(info.index,unique):raise ValueError('species aggregation alignment')
        genera=info.to_numpy(str)
        nphotos=int(good.sum()) if level=='within_species' else int(d.loc[good,'n_photos'].sum())
        entry={'n_evaluated_original_photos':nphotos,'n_evaluated_species':len(unique),
               'n_evaluated_genera':len(set(genera)),'folds':records,
               'target':'absolute_probabilities' if level=='absolute' else 'heldout_group_centred_contrasts_NOT_absolute_probabilities',
               'results':{}}
        contrasts=[('climate26_vs_geography','geography','climate26')]
        if mode=='soil':contrasts += [('environment36_vs_geography','geography','environment36'),
                                      ('soil_increment_same_photos','climate26','environment36')]
        for name,a,b in contrasts:
            entry['results'][name]=summarize(species_losses[a],species_losses[b],genera,reps,
                seed=SEED+sum(ord(c) for c in prefix+level+name))
            frame=pd.DataFrame({'species_id':unique,'genus':genera})
            for j,component in enumerate(COMPONENTS):
                frame['base_'+component]=species_losses[a][:,j]
                frame['full_'+component]=species_losses[b][:,j]
            frame.groupby('genus',sort=True).agg({**{col:'mean' for col in frame if col.startswith(('base_','full_'))},
                                                'species_id':'size'}).rename(columns={'species_id':'n_species'}).to_csv(
                outdir/f'{prefix}__{level}__{name}__genus_losses.csv')
        if level=='absolute':
            parts={k:genus_error_partition(y-v,d.genus) for k,v in pred.items()}
            entry['same_predictions_genus_mean_vs_species_deviation']={}
            for name,a,b in contrasts:
                dec=summarize(parts[a],parts[b],d.genus,reps,
                    seed=SEED+sum(ord(c) for c in prefix+'absolute'+name),
                    components=('total','genus_mean','within_genus_deviation'))
                for w in ('genus_equal','species_equal'):
                    if abs(dec[w]['total']['gain']-dec[w]['genus_mean']['gain']-
                           dec[w]['within_genus_deviation']['gain'])>1e-12:
                        raise ValueError('genus mean/deviation gain did not reassemble')
                entry['same_predictions_genus_mean_vs_species_deviation'][name]=dec
                frame=pd.DataFrame({'genus':d.genus})
                for j,comp in enumerate(('total','genus_mean','within_genus_deviation')):
                    frame['base_'+comp]=parts[a][:,j];frame['full_'+comp]=parts[b][:,j]
                frame.groupby('genus',sort=True).mean().to_csv(outdir/f'{prefix}__absolute__{name}__partition_genus_losses.csv')
        outcomes[level]=entry
        primary='environment36_vs_geography' if mode=='soil' else 'climate26_vs_geography'
        v=entry['results'][primary]['genus_equal']['total']
        print(json.dumps({'case':prefix,'level':level,'species':len(unique),'genera':len(set(genera)),
              'gain':v['gain'],'CI':v['gain_CI95'],'relative_percent':v['relative_gain_percent']},ensure_ascii=False),flush=True)
    return {'mode':mode,'quality':quality,'n_input_photos':len(photo),'n_input_species':len(sp),
            'n_input_genera':int(sp.genus.nunique()),'all_eligible_genera_tested_once':True,'levels':outcomes}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path,default=Path('/mnt/data'))
    p.add_argument('--outdir',type=Path,required=True);p.add_argument('--bootstrap',type=int,default=1999)
    p.add_argument('--mode',choices=('soil','climate','both'),default='both')
    p.add_argument('--quality',choices=('all','strict','both'),default='both');a=p.parse_args()
    if a.bootstrap<99:raise ValueError('At least 99 bootstrap draws')
    src.A=a.source_root/'fcp_original_photo_site_climate_soil_85337_20261010.zip'
    src.B=a.source_root/'fcp_100543_all36_expanded_photo_site_environment_20261011.zip'
    verify_file(src.A,SHA_PHOTO);verify_file(src.B,SHA_ENV)
    d,den,masks=src.run()
    if not (len(d)==100543 and den['classified']==45953):raise ValueError('source opportunity drift')
    a.outdir.mkdir(parents=True,exist_ok=True)
    result={'schema':'fcp_all36_environment_transfer_to_unseen_genera_20261011_v1','source_counts':den,
       'source_SHA256':{'photos':SHA_PHOTO,'environment':SHA_ENV},'full36_predictors':list(FULL),
       'alpha':ALPHA,'n_genus_folds':5,'bootstrap':a.bootstrap,'case_results':{},
       'scope':'source-exposed exploratory cross-genus transport; no causal selection or fitness identification',
       'original_main_H1H2_changed':False,
       'limitations':['Named-genus disjointness is not complete dated phylogeny control.',
          'Geography coordinates are adjusted but geographic cells are not withheld.',
          'Test contrast outcomes are centred ONLY for scoring, never used in predictions.',
          'Within-genus and within-species contrast losses are not absolute-probability Brier losses.',
          'Species/photo sources and complete cases differ; 26 vs36 soil effect uses same-soil-sample contrast only.',
          'Bootstrap intervals condition on fixed fitted folds and are not full model-refit intervals.',
          'Classification quality and MNAR, plasticity, lineage history and geographic covariance remain.',
          'No biochemical pigments, plant genotypes, pollination or fitness were measured.'],
       'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}}
    for mode in (('soil','climate') if a.mode=='both' else (a.mode,)):
        for qual in ((False,True) if a.quality=='both' else (a.quality=='strict',)):
            k=mode+('_quality' if qual else '_all')
            result['case_results'][k]=run_case(d.loc[masks[mode]],mode,qual,a.bootstrap,a.outdir)
            (a.outdir/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('SOURCE_REPLAY_COMPLETE',flush=True)
if __name__=='__main__':main()
