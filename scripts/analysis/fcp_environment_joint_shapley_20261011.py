#!/usr/bin/env python3
"""All 36 environmental variables together: out-of-source-cell whole-environment gain and grouped Shapley decomposition.

FCP source original 42,111 species / 100,543 photo ID union. Seven mechanistically distinct environmental
blocks are an explanatory grouping, NOT seven causal ecological factors. 5-fold original geography holdout,
source training genus means (>=3 training species), species photo means. Identical heldout species across all
2**7 models. Shapley contribution allocates shared nonlinear normalization/correlated predictor information,
not unique biophysical effects or mechanisms. Fixed source original model and alpha from 2026-10-11.
"""
from __future__ import annotations
import argparse,math,json,sys,hashlib,os
from pathlib import Path
from itertools import combinations
import numpy as np,pandas as pd
from sklearn.model_selection import GroupKFold
ROOT=Path(os.environ.get('FCP_SOURCE_ROOT','/mnt/data'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fcp_audit_36_preflight_20261011 import run as load_source,BLOCKS,FEATURES,CLASSES
from fcp_all36_global_region_cv_20261011 import aggregate_species,cell_index,ALPHA,SEED
GEO=('abs_lat','lon_sin','lon_cos')
SCHEMA='fcp_joint_abiotic_shapley_all_original_photo_species_geocell_v1'
N_FOLDS=5

def score(d,mode,region='world',nboot=999):
    blocks={k:v for k,v in BLOCKS.items() if mode=='soil' or k!='soil'}
    blocknames=list(blocks)
    nblock=len(blocknames)
    allcols=GEO+tuple(f for bb in blocks.values() for f in bb)
    colindex={x:i for i,x in enumerate(allcols)}
    masks=list(range(1<<nblock))
    groups={mask:np.array([colindex[f] for f in GEO+tuple(y for bit,b in enumerate(blocknames) if mask & (1<<bit) for y in blocks[b])],int) for mask in masks}
    d=aggregate_species(d)
    n=len(d);y=d[['hue_'+c for c in CLASSES]].to_numpy(float)
    predloss=np.full((1<<nblock,n),np.nan,dtype=np.float32)
    keep=np.zeros(n,bool)
    folds=[]
    for t,(train_ix,test_ix) in enumerate(GroupKFold(n_splits=N_FOLDS).split(d,groups=d.cell)):
        tr=d.iloc[train_ix];te=d.iloc[test_ix]
        counts=tr.genus.value_counts()
        select=np.flatnonzero(te.genus.isin(counts[counts>=3].index).to_numpy(bool))
        ix=test_ix[select]
        if len(ix)==0:continue
        tgt=te.iloc[select];keep[ix]=True
        # Exactly same source train species and original geographic-heldout tested rows for all 2**block models.
        assert not(set(tr.inat_taxon_id)&set(tgt.inat_taxon_id))
        xtr=tr[list(allcols)].to_numpy(float);xte=tgt[list(allcols)].to_numpy(float)
        mu=xtr.mean(axis=0);sd=xtr.std(axis=0);sd=np.where(sd<1e-8,1.,sd)
        xtr=(xtr-mu)/sd;xte=(xte-mu)/sd
        trgen=tr.genus.to_numpy(str);tegen=tgt.genus.to_numpy(str)
        yr=tr[['hue_'+c for c in CLASSES]].to_numpy(float)
        cy=['hue_'+c for c in CLASSES]
        means_y=pd.DataFrame(yr,columns=cy).assign(genus=trgen).groupby('genus',sort=False)[cy].mean()
        baseline=means_y.loc[tegen].to_numpy(float)
        ycenter=yr-means_y.loc[trgen].to_numpy(float)
        means_x=pd.DataFrame(xtr,columns=allcols).assign(genus=trgen).groupby('genus',sort=False)[list(allcols)].mean()
        xc=xtr-means_x.loc[trgen].to_numpy(float)
        xt=xte-means_x.loc[tegen].to_numpy(float)
        gram=xc.T@xc;rhs=xc.T@ycenter
        yt=y[ix]
        for mask in masks:
            use=groups[mask]
            coeff=np.linalg.solve(gram[np.ix_(use,use)]+ALPHA*np.eye(len(use)),rhs[use,:])
            pp=baseline+xt[:,use]@coeff
            pp=np.maximum(0,pp)
            tot=pp.sum(axis=1,keepdims=True)
            pp=np.divide(pp,tot,out=np.full_like(pp,.25),where=tot>0)
            predloss[mask,ix]=((pp-yt)**2).sum(axis=1)
        folds.append({'fold':t,'n_test_original_source_species':len(te),'n_classifiable_eligible_evaluated':len(ix),'n_original_geocells_heldout':te.cell.nunique()})
        print(mode,region,'fold',t,len(ix),'done',flush=True)
    if not keep.sum() or np.isnan(predloss[:,keep]).any():raise RuntimeError('Incomplete common heldout rows')
    L=predloss[:,keep].astype(float)
    # Brier gain is decrease in loss relative to geography-only model (mask 0).
    gains=L[0,None,:]-L
    gain_full=gains[-1]
    shapley=np.zeros((nblock,keep.sum()),float)
    for j in range(nblock):
        for S in masks:
            if S&(1<<j):continue
            m=S.bit_count()
            w=1/(nblock*math.comb(nblock-1,m))
            shapley[j]+=w*(gains[S|(1<<j)]-gains[S])
    assert np.allclose(shapley.sum(axis=0),gain_full,atol=1e-12)
    tax=d.loc[keep]
    genus=tax.genus.to_numpy(str)
    names,gid=np.unique(genus,return_inverse=True)
    ns=np.bincount(gid)
    data=np.vstack((gain_full,shapley))
    m=data.shape[0]
    # Cluster/bootstrap both species-weight and equal-genus weighting; both same heldout source rows.
    genusmean=np.vstack([np.bincount(gid,weights=row)/ns for row in data])
    gensum=genusmean*ns[None,:]
    rng=np.random.default_rng(SEED+len(region)*50+(17 if mode=='soil' else 0))
    # bootstrap in chunks to bound memory for many original genera
    boots_equal=[];boots_species=[]
    for _ in range(0,nboot,100):
        nb=min(100,nboot-len(boots_equal)*100)
        ids=rng.integers(0,len(names),size=(nb,len(names)))
        boots_equal.append(np.mean(genusmean[:,ids],axis=2).T)
        boots_species.append((np.sum(gensum[:,ids],axis=2)/np.sum(ns[ids],axis=1)[None,:]).T)
    be=np.concatenate(boots_equal);bs=np.concatenate(boots_species)
    assert len(be)==nboot
    base_loss=L[0].mean();full_loss=L[-1].mean()
    total={'species_equal_delta_brier':float(gain_full.mean()),
           'genus_equal_delta_brier':float(genusmean[0].mean()),
           'genus_equal_cluster95':np.quantile(be[:,0],[.025,.975]).tolist(),
           'species_equal_genus_cluster95':np.quantile(bs[:,0],[.025,.975]).tolist(),
           'geography_only_Brier':float(base_loss),'geography_plus_environment_Brier':float(full_loss),
           'relative_species_equal_Brier_reduction_pct':float(100*gain_full.mean()/base_loss),
           'relative_genus_equal_reduction_pct':float(100*genusmean[0].mean()/(np.bincount(gid,weights=L[0])/ns).mean())}
    sh={}
    for j,name in enumerate(blocknames):
        sh[name]={'species_equal_shapley_Brier':float(shapley[j].mean()),
                  'genus_equal_shapley_Brier':float(genusmean[j+1].mean()),
                  'genus_equal_cluster95':np.quantile(be[:,j+1],[.025,.975]).tolist(),
                  'species_equal_genus_cluster95':np.quantile(bs[:,j+1],[.025,.975]).tolist()}
    val={'schema':SCHEMA,'source_photo_species':n,'source_source_photo_count':int(d.n_photos.sum()),
         'n_oof_species':int(keep.sum()),'n_oof_genera':len(names),'n_oof_geocells':int(tax.cell.nunique()),
         'n_environment_blocks':nblock,'n_subsets_scored':1<<nblock,'n_bootstrap':nboot,
         'source_photo_hues_four_classes':list(CLASSES),'folds':folds,
         'full_environment_increment_beyond_geo':total,
         'grouped_shapley_of_joint_out_of_geo_cell_predictive_information':sh,
         'shapley_telescopes_to_full_environment_gain':True,
         'predictor_group_inventory':{name:list(blocks[name]) for name in blocknames},
         'hard_nonclaims':['These Shapley values allocate model-specific predictive information of correlated predictors; not measured causal benefits or costs.',
         'Species means from sampled photographs, not natural population colour phenotype, independent genomes or fitness.',
         'Genus-conditioned only, not full dated phylogeny or spatial Gaussian process.',
         'Source colour classification and actual photographed date/observer may correlate with real environmental exposures.',
         'This is a 5-fold spatial cell CV of per-species photo-location centroids, not a test of environmental selection within species.',
         'Biological pollinators/herbivores/sunscreen pigments and reproduction fitness not observed.']}
    return val

def run(modes=('soil','climate'),regions=('world',),boots=999,quality=False):
    d,counts,pools=load_source()
    result={'schema':SCHEMA,'date_jst':'2026-10-11','source':{'photo_union':100543,'species_opportunity':42111,'photo_original_classified':45953,
       'original_unclassifiable_kept_missing':54590,'all_36_abiotic_predictors':36,
       'original_full_climate_soil_case_coverage':counts},'cohorts':{}}
    all_classified=d[d.measurement_status.eq('classified_four_state_morph') & d.morph.isin(CLASSES)]
    quality_pixels=23826.0 # Frozen earlier technical quality global photo gate
    quality_nuisance=0.1375717742195135 # Frozen earlier technical quality global photo gate
    result['quality_rule']={'median_flower_effective_pixels_threshold':quality_pixels,
       'median_nuisance_pixel_fraction_max':quality_nuisance} if quality else None
    for mode in modes:
        dsrc=d.loc[pools[mode]].copy()
        if quality:
            dsrc=dsrc.loc[(dsrc.flower_effective_pixels>=quality_pixels) &
                          (dsrc.nuisance_pixel_fraction<=quality_nuisance)].copy()
        for region in regions:
            dr=dsrc if region=='world' else dsrc[dsrc.east] if region=='east' else dsrc[~dsrc.east]
            result['cohorts'][mode+'_'+region+('_quality' if quality else '')]=score(dr,mode,region,boots)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=['soil','climate','both'],default='both')
    p.add_argument('--region',choices=['world','east','elsewhere','all'],default='world')
    p.add_argument('--bootstrap',type=int,default=499)
    p.add_argument('--quality',action='store_true')
    p.add_argument('--outdir',type=Path,default=ROOT/'fcp_joint_environment_results_20261011')
    a=p.parse_args()
    modes=('soil','climate') if a.mode=='both' else (a.mode,)
    regions=('world','east','elsewhere') if a.region=='all' else (a.region,)
    r=run(modes,regions,a.bootstrap,quality=a.quality)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/'result.json').write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
    print('RESULT',json.dumps({k:{'n':v['n_oof_species'],'total':v['full_environment_increment_beyond_geo'],'shapley':{x:round(y['genus_equal_shapley_Brier'],7) for x,y in v['grouped_shapley_of_joint_out_of_geo_cell_predictive_information'].items()}} for k,v in r['cohorts'].items()},indent=2))