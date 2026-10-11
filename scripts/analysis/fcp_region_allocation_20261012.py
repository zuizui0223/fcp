#!/usr/bin/env python3
"""Posthoc species-centered geographical colour allocation: source-frozen predictions.

Only PHOTO sector fractions; biological within-population C / among-population S
cannot be inferred. Does not fit/retrain any colour or environment model.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
import platform
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd

SOURCE_SHA = 'e51547cf9daea2bac68d0ab1fbffe17200cdecd4b810ebb660040df4c40f05a1'
MEMBER_ROOT = 'fcp_nonlinear_transport_20261011/results/'
SOURCE = 'FCP_全36環境_非線形と地域外予測_20261011.zip'
CLASSES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')
MODELS = ('ridge_geo','ridge_full','nonlinear_geo','nonlinear_full')
SCHEMES = ('genus_only', 'sector_purged','buffer500')
COMPONENTS = ('total','white_axis','remaining_hue')
CONTRASTS = {
    'environment_nonlinear': ('nonlinear_geo','nonlinear_full'),
    'environment_ridge': ('ridge_geo','ridge_full'),
    'full_nonlinear_vs_ridge': ('ridge_full','nonlinear_full'),
    'geo_nonlinear_vs_ridge': ('ridge_geo','nonlinear_geo'),
    'zero_vs_nonlinear_full': ('zero','nonlinear_full'),
    'zero_vs_ridge_full': ('zero','ridge_full'),
}
KEY=['inat_taxon_id','sector','scheme']
SEED_BOOT=2026101201
SEED_NULL=2026101202


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for v in iter(lambda:f.read(1<<20),b''):h.update(v)
    return h.hexdigest()


def _member(case): return MEMBER_ROOT+case+'_oof_predictions.csv.gz'


def load_case(z,case):
    d=pd.read_csv(io.BytesIO(z.read(_member(case))),compression='gzip')
    if d.empty or d.duplicated(KEY).any():raise ValueError('empty or duplicate species-sector-scheme records')
    if d[KEY+['genus','fold','n_photos']].isna().any().any():raise ValueError('bad IDs')
    if set(d.scheme.unique())!=set(SCHEMES):raise ValueError('missing schemes')
    if not set(d.sector.unique())=={0,1,2,3}:raise ValueError('wrong sector labels')
    if (d.n_photos<=0).any() or not np.allclose(d.n_photos,np.rint(d.n_photos),rtol=0,atol=1e-10):
        raise ValueError('invalid integer source photo counts')
    base=d[d.scheme.eq('genus_only')].set_index(['inat_taxon_id','sector']).sort_index()
    truthcols=['truth__'+h for h in CLASSES]
    invariant=['genus','fold','n_photos',*truthcols]
    for scheme in SCHEMES[1:]:
        x=d[d.scheme.eq(scheme)].set_index(['inat_taxon_id','sector']).sort_index()
        if not x.index.equals(base.index) or not x[invariant].equals(base[invariant]):
            raise ValueError('evaluation target drift between original source schemes')
    if d.groupby('inat_taxon_id').genus.nunique().max()!=1 or d.groupby('genus').fold.nunique().max()!=1:
        raise ValueError('genus/fold consistency failure')
    for name,cols in [('truth',truthcols),*[(m,[m+'__p__'+h for h in CLASSES]) for m in MODELS]]:
        p=d[cols].to_numpy(float)
        if not np.isfinite(p).all() or (p < -1e-10).any() or (p>1+1e-10).any() or not np.allclose(p.sum(1),1.,rtol=0,atol=1e-10):
            raise ValueError(name+' source probability simplex violation')
    if len(d)!= (70416 if case=='soil_all' else 21006):raise ValueError('original prediction row-count drift')
    return d


def opportunity_gate(d,min_n=2):
    # Colour values are never read to decide which source species qualify.
    slots=d[d.n_photos.ge(min_n)][['inat_taxon_id','sector']]
    ok=slots.groupby('inat_taxon_id').size()
    selected=set(ok.index[ok.ge(2)])
    return d[d.inat_taxon_id.isin(selected)&d.n_photos.ge(min_n)].copy()


def all_model_losses(frame):
    if frame.empty:return frame
    frame=frame.sort_values(['inat_taxon_id','sector']).copy()
    group=frame.groupby('inat_taxon_id',sort=False)
    tr=frame[['truth__'+h for h in CLASSES]].to_numpy(float)
    y=tr-group[['truth__'+h for h in CLASSES]].transform('mean').to_numpy(float)
    result=frame[['inat_taxon_id','genus','sector','scheme','n_photos','fold','environment_support',*['truth__'+h for h in CLASSES]]].copy()
    def components(r):
        total=(r*r).sum(axis=1)
        white=(4./3.)*r[:,0]**2
        hue=total-white
        if (hue < -1e-10).any():raise ValueError('orthogonal component invalid')
        return np.column_stack([total,white,hue])
    for j,c in enumerate(COMPONENTS):result['zero__'+c]=components(y)[:,j]
    for m in MODELS:
        pc=[m+'__p__'+h for h in CLASSES]
        p=frame[pc].to_numpy(float) - group[pc].transform('mean').to_numpy(float)
        scores=components(y-p)
        for j,c in enumerate(COMPONENTS):result[m+'__'+c]=scores[:,j]
        # Fixed photo-sector predictions used to check the exact permutation identity.
        for j,h in enumerate(CLASSES):result[m+'__centred_p__'+h]=p[:,j]
    for name,(a,b) in CONTRASTS.items():
        for c in COMPONENTS:result[name+'__'+c]=result[a+'__'+c]-result[b+'__'+c]
    return result


def species_ledger(reg):
    cols=[c for c in reg if any(c.endswith('__'+x) for x in COMPONENTS)]
    g=reg.groupby('inat_taxon_id',sort=True)
    out=g[cols].mean().reset_index()
    out['genus']=g.genus.first().to_numpy()
    out['n_sectors']=g.sector.nunique().to_numpy()
    out['n_original_photos']=g.n_photos.sum().to_numpy()
    if out.inat_taxon_id.nunique()!=len(out):raise ValueError('duplicated summary species')
    for c in COMPONENTS:
        for model in ('zero',*MODELS):
            if not np.allclose(out[model+'__total'] if c=='total' else out[model+'__total'],
                out[model+'__white_axis']+out[model+'__remaining_hue'],rtol=0,atol=1e-11):
                raise ValueError('decomposition fails')
    return out


def summarize_ledger(s,seed,with_ci=True):
    numerical=[c for c in s if any(c.endswith('__'+x) for x in COMPONENTS)]
    data=s.groupby('genus',sort=True)[numerical].mean()
    if with_ci:
        rng=np.random.default_rng(seed)
        ng=len(data)
        idx=rng.integers(0,ng,size=(1999,ng))
        draws_genus=data.to_numpy(float)[idx].mean(axis=1)
        # cluster bootstrap with species-equal estimator, resampling whole genera
        sums=s.groupby('genus',sort=True)[numerical].sum().loc[data.index].to_numpy(float)
        cnt=s.groupby('genus',sort=True).size().loc[data.index].to_numpy(float)
        draws_species=sums[idx].sum(axis=1)/cnt[idx].sum(axis=1)[:,None]
    out={'n_species':len(s),'n_genera':len(data),
        'n_photos':int(s.n_original_photos.sum()),'n_species_sector_records':int(s.n_sectors.sum()),
        'n_sectors_distribution':{str(k):int(v) for k,v in s.n_sectors.value_counts().items()},
        'bootstrap_reps':1999 if with_ci else 0,'seed':seed,'status':'OK' if with_ci else 'HOLD_LOW_SUPPORT',
        'inference': 'fixed source models; paired full-genus bootstrap' if with_ci else 'descriptive only'}
    out['genus_equal']={}
    out['species_equal']={}
    for name,point,draws in [('genus_equal',data.mean(),draws_genus if with_ci else None),('species_equal',s[numerical].mean(),draws_species if with_ci else None)]:
        for key in numerical:
            item={'estimate':float(point[key])}
            if with_ci:
                j=numerical.index(key)
                item['CI95']=[float(z) for z in np.quantile(draws[:,j],[.025,.975])]
                item['CI97_5']=[float(z) for z in np.quantile(draws[:,j],[.0125,.9875])]
            out[name][key]=item
    return out


def photo_sector_null(d, ledger, reps=1999, seed=SEED_NULL):
    """Analytic linear identity plus exact within-species photo label shuffles."""
    if not d.scheme.eq('buffer500').all():raise ValueError('primary null must use buffer500')
    if not d.groupby('inat_taxon_id').size().ge(2).all():raise ValueError('missing paired sector')
    gsort=ledger.groupby('genus').size()
    ng=len(gsort)
    if len(ledger)<100 or ng<50:raise ValueError('insufficient primary support')
    genus_weight=ledger.genus.map((1/gsort/ng).to_dict()).to_numpy(float)
    obs_g=float(np.dot(genus_weight,ledger['environment_nonlinear__total']))
    obs_organisation=float(np.dot(genus_weight,ledger['zero__total']))
    rng=np.random.default_rng(seed)
    delta_perm=np.zeros(reps)
    org_perm=np.zeros(reps)
    count_audit=[]
    permutation_identity_max_error=0.
    idx_by_sp={k:v for k,v in d.sort_values(['inat_taxon_id','sector']).groupby('inat_taxon_id',sort=True)}
    for i,s in ledger.reset_index(drop=True).iterrows():
        group=idx_by_sp[s.inat_taxon_id]
        K=len(group)
        n=group.n_photos.to_numpy(int)
        tr=group[['truth__'+h for h in CLASSES]].to_numpy(float)
        counts=np.rint(tr*n[:,None]).astype(int)
        if not np.allclose(counts/n[:,None],tr,rtol=0,atol=1e-8):
            raise ValueError('nonintegral reconstructed four-class photo counts')
        if np.any(counts<0) or not np.array_equal(counts.sum(axis=1),n):raise ValueError('invalid photo inventory')
        original_total=counts.sum(axis=0)
        labels=np.repeat(np.arange(4),original_total)
        assert labels.size==n.sum()
        # draw a uniform random bijection of fixed original photo label inventory
        permidx=np.argsort(rng.random(size=(reps,len(labels))),axis=1)
        perm=labels[permidx]  # reps x n_total
        new_counts=np.zeros((reps,K,4),dtype=int)
        start=0
        for k,m in enumerate(n):
            new_counts[:,k,:]=np.eye(4,dtype=int)[perm[:,start:start+m]].sum(axis=1)
            start+=m
        if not np.array_equal(new_counts.sum(axis=(1)),np.tile(original_total,(reps,1))):
            raise ValueError('null violates species colour totals')
        if not np.array_equal(new_counts.sum(axis=2),np.tile(n,(reps,1))):
            raise ValueError('null violates geographic photo slot counts')
        # convert labels back to original sector fractions
        yc=new_counts/n[None,:,None]
        yc=yc-yc.mean(axis=1,keepdims=True)
        pg=group[['nonlinear_geo__centred_p__'+h for h in CLASSES]].to_numpy(float)
        pf=group[['nonlinear_full__centred_p__'+h for h in CLASSES]].to_numpy(float)
        obs_colour=tr-tr.mean(0,keepdims=True)
        # MSE difference is 2 Y . (Pfull - Pgeo) + (||Pgeo||² - ||Pfull||²)
        offset=float(((pg**2).sum(1)-(pf**2).sum(1)).mean())
        calc_obs=float(2*np.einsum('sk,sk->',obs_colour,pf-pg)/K+offset)
        permutation_identity_max_error=max(permutation_identity_max_error,abs(calc_obs-s['environment_nonlinear__total']))
        if abs(calc_obs-s['environment_nonlinear__total'])>1e-10:
            raise ValueError('primary permuted-loss algebra mismatch')
        dlt=2*np.einsum('rsk,sk->r',yc,pf-pg)/K+offset
        o=(yc*yc).sum(axis=(1,2))/K
        delta_perm+=genus_weight[i]*dlt
        org_perm+=genus_weight[i]*o
        count_audit.append((int(s.inat_taxon_id),K,int(sum(n)),tuple(map(int,original_total))))
    if len(count_audit)!=len(ledger):raise ValueError('permutation cohort mismatch')
    return {'null_reps':reps,'null_seed':seed,'n_species':len(ledger),'n_genera':ng,
            'null_eligibility':'original fixed sector slots + exact original within-species four-colour photo inventory',
            'primary_observed_environment_gain':obs_g,
            'primary_null_mean':float(delta_perm.mean()),
            'primary_null_std':float(delta_perm.std(ddof=1)),
            'primary_null_CI95':np.quantile(delta_perm,[.025,.975]).tolist(),
            'primary_plus_one_upper_tail_p':float((1+np.sum(delta_perm>=obs_g))/(reps+1)),
            'descriptive_spatial_colour_organisation':{
                'observed':obs_organisation,
                'null_mean':float(org_perm.mean()),
                'plus_one_upper_tail_p':float((1+np.sum(org_perm>=obs_organisation))/(reps+1))},
            'reconstruction_status':'ALL_SOURCE_COUNTS_EXACT_AND_PERMUTATION_INVENTORY_PRESERVED',
            'max_observed_gain_identity_error':float(permutation_identity_max_error),
            'species_photo_slot_counts_sha256':hashlib.sha256(json.dumps(count_audit).encode()).hexdigest()}


def run(source,outdir):
    source=Path(source);outdir=Path(outdir)
    if sha256_file(source)!=SOURCE_SHA:raise ValueError('source ZIP SHA256 drift')
    outdir.mkdir(parents=True,exist_ok=True)
    report={'schema':'fcp_20261012_species_sector_colouroffsets_v1','status':'RUNNING',
        'source_package_sha256':SOURCE_SHA,'source_package_bytes':source.stat().st_size,
        'outcome_exposure':'retrospective, source-exposed, protocol fixed before the NEW region-centred contrast was examined',
        'models':'FOUR_PR170_FIXED_UNSEEN_GENUS_OUT_OF_SAMPLE_MODELS_NO_REFITTING',
        'source_photo_denominators':{'all_original_photo_opportunities':100543,'classifiable_original':45953,
                                      'unclassifiable_original':54590,'all36_classified':31925},
        'primary_endpoint':{'case':'soil_all','scheme':'buffer500','min_photos_per_sector':2,
            'n_required_species':100,'n_required_genera':50,'weight':'genus_equal','contrast':'environment_nonlinear__total'},
        'bootstrap_seed':SEED_BOOT,'null_seed':SEED_NULL,
        'available_models':list(MODELS),'cases':{},'null':None,
        'original_main_H1H2_changed':False}
    with zipfile.ZipFile(source) as z:
        for case in ('soil_all','soil_quality'):
            raw=load_case(z,case)
            report['cases'][case]={}
            for min_n in (2,3):
                tag='min'+str(min_n)
                selected=opportunity_gate(raw[raw.scheme.eq('buffer500')],min_n)
                species=set(selected.inat_taxon_id)
                subset=raw[raw.inat_taxon_id.isin(species)&raw.n_photos.ge(min_n)]
                report['cases'][case][tag]={}
                print(json.dumps({'case':case,'gate':tag,'n_species':len(species),
                       'n_genera':selected.genus.nunique(),'n_sectors':len(selected)}),flush=True)
                for j,scheme in enumerate(SCHEMES):
                    d=subset[subset.scheme.eq(scheme)].copy()
                    reg=all_model_losses(d)
                    ledger=species_ledger(reg)
                    with_ci=len(ledger)>=100 and ledger.genus.nunique()>=50
                    result=summarize_ledger(ledger,SEED_BOOT+17*j+53*min_n+(0 if case=='soil_all' else 500),with_ci)
                    report['cases'][case][tag][scheme]=result
                    if scheme=='buffer500':
                        ledger.to_csv(outdir/f'{case}_{tag}_{scheme}_species_losses.csv.gz',index=False,compression='gzip')
                        if case=='soil_all' and min_n==2:
                            report['null']=photo_sector_null(reg,ledger)
                # Check response sets and same target sector identities across all three schemes.
                targetsets=[set(map(tuple,subset[subset.scheme.eq(s)][['inat_taxon_id','sector']].to_numpy())) for s in SCHEMES]
                if not all(targetsets[0]==v for v in targetsets[1:]):raise ValueError('target drift')
            # Check all three regimes have identical performance on same y/sector means baseline.
            for tag,res in report['cases'][case].items():
                zeroes=[res[s]['genus_equal']['zero__total']['estimate'] for s in SCHEMES]
                if not np.allclose(zeroes,zeroes[0],rtol=0,atol=1e-12):raise ValueError('observed target baseline drift')
    report['status']='ANALYSIS_COMPLETE'
    report['runtime']={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
    (outdir/'result.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True)
    p.add_argument('--outdir',required=True)
    a=p.parse_args()
    r=run(a.source,a.outdir)
    z=r['cases']['soil_all']['min2']['buffer500']
    print(json.dumps({'status':r['status'],'primary_n_species':z['n_species'],'n_genera':z['n_genera'],
         'primary':z['genus_equal']['environment_nonlinear__total'],
         'null':r['null']['primary_plus_one_upper_tail_p']}),flush=True)

if __name__=='__main__':main()