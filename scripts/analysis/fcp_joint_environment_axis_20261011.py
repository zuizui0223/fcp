#!/usr/bin/env python3
"""Joint-environment prediction: white/nonwhite axis vs remaining hue contrasts.

Retrospective, same-source, same original models as PR166/167. This is squared
prediction error on species' photographic colour-composition means, not fitness,
causal variance partition, or the independently frozen 9-colour H2 statistic.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from fcp_audit_36_preflight_20261011 import run as load_source, BLOCKS, CLASSES
from fcp_all36_global_region_cv_20261011 import aggregate_species, fit_predict

GEO=('abs_lat','lon_sin','lon_cos')
COMPONENTS=('white_axis','remaining_hue','total')
ARCHIVES={
 'fcp_original_photo_site_climate_soil_85337_20261010.zip':'2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0',
 'fcp_100543_all36_expanded_photo_site_environment_20261011.zip':'6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca',
}


def components(truth:np.ndarray,prediction:np.ndarray)->np.ndarray:
    y,p=np.asarray(truth,float),np.asarray(prediction,float)
    if y.shape!=p.shape or y.ndim!=2 or y.shape[1]!=4:
        raise ValueError('Expected identical n by 4 simplex arrays')
    for a in (y,p):
        if not np.isfinite(a).all() or (a < -1e-12).any() or not np.allclose(a.sum(axis=1),1,atol=1e-10,rtol=0):
            raise ValueError('Both truth and predictions must be finite probability compositions')
    r=p-y
    total=np.sum(r*r,axis=1)
    axis=(4./3.)*r[:,0]**2
    hue=total-axis
    if (hue < -1e-11).any():
        raise ValueError('Orthogonal residual loss unexpectedly negative')
    hue=np.maximum(hue,0.)
    out=np.column_stack([axis,hue,total])
    if not np.allclose(out[:,0]+out[:,1],out[:,2],atol=1e-12,rtol=0):
        raise RuntimeError('Orthogonal loss does not reassemble total')
    return out


def ratio(numer,denom):
    return float(numer/denom) if abs(denom)>1e-14 else None


def ci(a):
    a=np.asarray(a,float)
    if not np.isfinite(a).all():return None
    return np.quantile(a,[.025,.975]).tolist()


def paired_summary(before,after,genus,repetitions=1999,seed=2026101142):
    before,after=np.asarray(before,float),np.asarray(after,float)
    if before.shape!=after.shape or before.ndim!=2 or before.shape[1]!=3:
        raise ValueError('Expected n by 3 loss components')
    if len(genus)!=len(before) or not np.isfinite(before).all() or not np.isfinite(after).all():
        raise ValueError('Invalid paired observations')
    names,gid=np.unique(np.asarray(genus,str),return_inverse=True)
    k=len(names)
    if k<2 or repetitions<2:raise ValueError('At least 2 genera and bootstrap draws required')
    size=np.bincount(gid,minlength=k)
    a=np.stack([np.bincount(gid,weights=before[:,j],minlength=k) for j in range(3)],axis=1)
    b=np.stack([np.bincount(gid,weights=after[:,j],minlength=k) for j in range(3)],axis=1)
    am,bm=a/size[:,None],b/size[:,None]
    rng=np.random.default_rng(seed)
    boot_a={w:[] for w in ('equal_genus','equal_species')};boot_b={w:[] for w in boot_a}
    for start in range(0,repetitions,64):
        draw=rng.integers(k,size=(min(64,repetitions-start),k))
        boot_a['equal_genus'].append(am[draw].mean(axis=1))
        boot_b['equal_genus'].append(bm[draw].mean(axis=1))
        denominator=size[draw].sum(axis=1)[:,None]
        boot_a['equal_species'].append(a[draw].sum(axis=1)/denominator)
        boot_b['equal_species'].append(b[draw].sum(axis=1)/denominator)
    out={'n_genera':k,'n_evaluation_species':len(before),'bootstrap_repetitions':repetitions,
         'uncertainty':'paired resampling of held-out genera, conditional on fixed original fitted models'}
    for weight in boot_a:
        b0=np.concatenate(boot_a[weight]);b1=np.concatenate(boot_b[weight]);bd=b0-b1
        m0=am.mean(axis=0) if weight=='equal_genus' else before.mean(axis=0)
        m1=bm.mean(axis=0) if weight=='equal_genus' else after.mean(axis=0)
        gain=m0-m1
        s={}
        for j,name in enumerate(COMPONENTS):
            rel=ratio(gain[j],m0[j]);rel_boot=bd[:,j]/b0[:,j] if np.all(b0[:,j]>1e-14) else None
            s[name]={'baseline_loss':float(m0[j]),'environment_loss':float(m1[j]),
                     'gain':float(gain[j]),'gain_ci95':ci(bd[:,j]),
                     'relative_error_reduction':rel,
                     'relative_error_reduction_ci95':ci(rel_boot) if rel_boot is not None else None}
        s['white_share_of_baseline_loss']=ratio(m0[0],m0[2])
        s['raw_signed_white_gain_ratio']=ratio(gain[0],gain[2])
        unstable=bool(ci(bd[:,2])[0]<=0<=ci(bd[:,2])[1])
        s['gain_share_is_unstable_if_total_ci_includes_zero']=unstable
        s['white_share_of_total_gain']=ratio(gain[0],gain[2]) if not unstable and gain[2]>0 and np.all(gain[:2]>=0) else None
        left,right=s['white_axis']['relative_error_reduction'],s['remaining_hue']['relative_error_reduction']
        diff=None if left is None or right is None else left-right
        diff_b=None
        if np.all(b0[:,:2]>1e-14):diff_b=bd[:,0]/b0[:,0]-bd[:,1]/b0[:,1]
        s['axis_minus_hue_relative_skill']={'estimate':diff,'paired_genus_ci95':ci(diff_b) if diff_b is not None else None}
        # Extra axis gain after allocating gains in proportion to baseline losses.
        base_fraction=ratio(m0[0],m0[2])
        s['excess_axis_gain_vs_baseline_error_share']={
            'estimate':float(gain[0]-base_fraction*gain[2]) if base_fraction is not None else None,
            'ci95':ci(bd[:,0]-(b0[:,0]/b0[:,2])*bd[:,2]) if np.all(b0[:,2]>1e-14) else None}
        out[weight]=s
    return out


def compute(rows,mode,region,quality,repetitions,outdir):
    species=aggregate_species(rows)
    target=species[['hue_'+c for c in CLASSES]].to_numpy(float)
    p0=np.full_like(target,np.nan);p1=np.full_like(target,np.nan)
    use=list(GEO)+[col for block,cols in BLOCKS.items() if mode=='soil' or block!='soil' for col in cols]
    fold_ids=np.full(len(species),-1,dtype=int);ledger=[]
    for fold,(itr,ite) in enumerate(GroupKFold(n_splits=5).split(species,groups=species.cell)):
        train,test=species.iloc[itr],species.iloc[ite]
        count=train.genus.value_counts()
        selected=np.flatnonzero(test.genus.isin(count[count>=3].index).to_numpy())
        focus=test.iloc[selected]
        if set(train.inat_taxon_id)&set(focus.inat_taxon_id):raise RuntimeError('Species training leakage')
        if set(train.cell)&set(focus.cell):raise RuntimeError('Species-centroid cell leakage')
        ledger.append({'fold':fold,'n_test':len(test),'n_evaluable':len(focus)})
        if not len(focus):continue
        p0[ite[selected]]=fit_predict(train,focus,list(GEO))
        p1[ite[selected]]=fit_predict(train,focus,use)
        fold_ids[ite[selected]]=fold
    valid=fold_ids>=0
    if int(valid.sum())<200:raise RuntimeError('Insufficient evaluation population')
    before=components(target[valid],p0[valid]);after=components(target[valid],p1[valid])
    key=mode+'_'+region+('_strict_quality' if quality else '')
    seed=int.from_bytes(hashlib.sha256(('fcp-axis-relative-skill-v1|'+key).encode()).digest()[:8],'little')
    result=paired_summary(before,after,species.loc[valid,'genus'],repetitions,seed)
    result.update({'case':key,'n_source_photo_rows':len(rows),'n_source_species':len(species),
                   'n_evaluation_photo_rows':int(species.loc[valid,'n_photos'].sum()),
                   'n_evaluation_centroid_cells':int(species.loc[valid,'cell'].nunique()),
                   'mode':mode,'region':region,'strict_photo_quality':quality,'folds':ledger})
    original_white_delta=((p0[valid,0]-target[valid,0])**2-(p1[valid,0]-target[valid,0])**2)
    if not np.allclose(before[:,0]-after[:,0],4/3*original_white_delta,atol=1e-12):
        raise RuntimeError('White factor audit failed')
    # Preserve compact per-genus sufficient summaries, no individual-photo labels.
    detail=pd.DataFrame({'genus':species.loc[valid,'genus'].to_numpy()})
    for j,name in enumerate(COMPONENTS):
        detail['baseline_'+name]=before[:,j];detail['environment_'+name]=after[:,j]
    summary=detail.groupby('genus').agg(['size','mean'])
    summary.columns=['_'.join(x) for x in summary.columns]
    summary.to_csv(outdir/(key+'_genus_loss_audit.csv'))
    return result


def run(source_root,outdir,repetitions=1999):
    root=Path(source_root);outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    for name,expected in ARCHIVES.items():
        got=hashlib.sha256((root/name).read_bytes()).hexdigest()
        if got!=expected:raise RuntimeError('Original archive SHA256 mismatch: '+name)
    # Resolve the original portable preflight without changing any file content.
    import fcp_audit_36_preflight_20261011 as pre
    pre.A=root/'fcp_original_photo_site_climate_soil_85337_20261010.zip'
    pre.B=root/'fcp_100543_all36_expanded_photo_site_environment_20261011.zip'
    data,counts,pools=load_source()
    results={}
    expected_n={'climate_world':17183,'climate_east':4555,'climate_elsewhere':12653,
                'soil_world':13961,'soil_east':3307,'soil_elsewhere':10307,
                'climate_world_strict_quality':4512,'soil_world_strict_quality':3399}
    for mode in ('soil','climate'):
        for quality in (False,True):
            rows=data.loc[pools[mode]].copy()
            if quality:rows=rows[(rows.flower_effective_pixels>=23826)&(rows.nuisance_pixel_fraction<=0.1375717742195135)]
            for region in (('world',) if quality else ('world','east','elsewhere')):
                subset=rows if region=='world' else rows[rows.east] if region=='east' else rows[~rows.east]
                res=compute(subset,mode,region,quality,repetitions,outdir)
                key=res['case']
                if res['n_evaluation_species']!=expected_n[key]:raise RuntimeError('Historical evaluation count drift: '+key)
                results[key]=res
                print(key, res['n_evaluation_species'], 'GAIN',res['equal_genus']['total']['gain'],
                      'AXIS-SKILL',res['equal_genus']['white_axis']['relative_error_reduction'],
                      'HUE-SKILL',res['equal_genus']['remaining_hue']['relative_error_reduction'],flush=True)
    import sklearn,scipy
    out={'schema':'fcp_joint_environment_orthogonal_white_hue_relative_skill_20261011_v1',
         'source_sha256':ARCHIVES,'source_counts':counts,'bootstrap_repetitions':repetitions,
         'primary_case':'soil_world','primary_specificity_contrast':'axis minus hue relative error reduction',
         'models':'unchanged genus-centred ridge alpha40; original species-centroid-cell fivefold holdout',
         'environment_blocks':{k:list(v) for k,v in BLOCKS.items()},
         'definition':'L_total=L_white_axis+L_hue; L_white_axis=4/3*(predicted_white-observed_white)^2',
         'cases':results,
         'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
                    'scikit_learn':sklearn.__version__,'scipy':scipy.__version__},
         'boundaries':['Four-state photographic response; not the continuous nine-colour prospective H2 statistic.',
                       'Baseline-error-normalized skills control unequal error opportunity, not all class-dependent measurement bias.',
                       'Species means use all eligible original photographs; 54590 colour-missing photos not imputed.',
                       'Geographic holdout uses species centroids, not all original photograph sites.',
                       'Fixed-model cluster resampling does not include refitting, phylogenetic or spatial GP uncertainty.',
                       'Other regions and quality subsets are exploratory, differing target taxa; no causal cross-region contrast.',
                       'No genotype, flower pigment chemistry, pollination, herbivory, individual fitness or energy costs observed.',
                       'Current main manuscript, H1/H2 and reserved taxa unchanged.']}
    (outdir/'result.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source-root',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--bootstrap',type=int,default=1999)
    a=p.parse_args();run(a.source_root,a.outdir,a.bootstrap)
