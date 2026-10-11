#!/usr/bin/env python3
"""Independent finite-outcome validation, imports no analysis functions."""
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np,pandas as pd

CLASSES=('white','yellow_orange','red_pink','blue_purple')
MODELS=('ridge_geo','ridge_full','nonlinear_geo','nonlinear_full')
SCH=('genus_only','sector_purged','buffer500')
COMP=('total','white_axis','remaining_hue')
MAPPING={'environment_nonlinear':('nonlinear_geo','nonlinear_full'),
         'environment_ridge':('ridge_geo','ridge_full'),
         'full_nonlinear_vs_ridge':('ridge_full','nonlinear_full'),
         'geo_nonlinear_vs_ridge':('ridge_geo','nonlinear_geo'),
         'zero_vs_nonlinear_full':('zero','nonlinear_full'),
         'zero_vs_ridge_full':('zero','ridge_full')}
SRC='e51547cf9daea2bac68d0ab1fbffe17200cdecd4b810ebb660040df4c40f05a1'
ROOT='fcp_nonlinear_transport_20261011/results/'


def verify(source,results):
    h=hashlib.sha256(source.read_bytes()).hexdigest()
    assert h==SRC
    out=json.loads((results/'result.json').read_text())
    assert out['source_package_sha256']==SRC and out['status']=='ANALYSIS_COMPLETE'
    checks=0; maxdiff=0.; null_mc_gap=0.;cases={}; analytic_null_mean=None
    def close(x,y,name,tol=1e-11):
        nonlocal checks,maxdiff,null_mc_gap
        a=np.asarray(x,dtype=float);b=np.asarray(y,dtype=float)
        assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),name
        err=float(abs(a-b).max()) if a.size else 0
        
        if tol<=1e-11:maxdiff=max(maxdiff,err)
        else:null_mc_gap=max(null_mc_gap,err)
        checks+=1
        if err>tol:raise ValueError(f'{name} mismatch {err}')
    with zipfile.ZipFile(source) as z:
        for case in ('soil_all','soil_quality'):
            d=pd.read_csv(io.BytesIO(z.read(ROOT+case+'_oof_predictions.csv.gz')),compression='gzip')
            assert set(d.scheme)==set(SCH)
            cases[case]={}
            for n in (2,3):
                gate=d[(d.scheme=='buffer500')&(d.n_photos>=n)].groupby('inat_taxon_id').size()
                eligible=set(gate[gate>=2].index)
                dd=d[(d.inat_taxon_id.isin(eligible))&(d.n_photos>=n)]
                cases[case]['min'+str(n)]={}
                for scheme in SCH:
                    group=dd[dd.scheme==scheme]
                    rows=[]; offsets=[]
                    for species, one in group.groupby('inat_taxon_id'):
                        if len(one)<2:raise ValueError('incomplete sector')
                        truth=one[['truth__'+c for c in CLASSES]].to_numpy(float)
                        truth=truth-truth.mean(axis=0)
                        losses={}
                        for m in ('zero',*MODELS):
                            pred=np.zeros_like(truth) if m=='zero' else one[[m+'__p__'+c for c in CLASSES]].to_numpy(float)
                            if m!='zero':pred-=pred.mean(axis=0)
                            residual=truth-pred
                            total=np.mean(np.einsum('ij,ij->i',residual,residual))
                            white=(4/3)*np.mean(residual[:,0]**2)
                            losses[m]=(total,white,total-white)
                        row={'inat_taxon_id':species,'genus':one.genus.iloc[0],
                             'n_sectors':len(one),'n_photos':int(one.n_photos.sum())}
                        for m,vals in losses.items():
                            for c,v in zip(COMP,vals):row[m+'__'+c]=v
                        for name,(a,b) in MAPPING.items():
                            for j,c in enumerate(COMP):row[name+'__'+c]=losses[a][j]-losses[b][j]
                        rows.append(row)
                        if case=='soil_all' and n==2 and scheme=='buffer500':
                            pg=one[['nonlinear_geo__p__'+c for c in CLASSES]].to_numpy(float)
                            pf=one[['nonlinear_full__p__'+c for c in CLASSES]].to_numpy(float)
                            pg-=pg.mean(0);pf-=pf.mean(0)
                            offsets.append((one.genus.iloc[0],species,np.mean((pg*pg-pf*pf).sum(1))))
                    frame=pd.DataFrame(rows)
                    summary=out['cases'][case]['min'+str(n)][scheme]
                    assert len(frame)==summary['n_species'] and len(frame.genus.unique())==summary['n_genera']
                    assert int(frame.n_photos.sum())==summary['n_photos']
                    gen=frame.groupby('genus',sort=True).mean(numeric_only=True)
                    calc_gen=gen.mean();calc_sp=frame.mean(numeric_only=True)
                    for col in (k for k in frame.columns if k.endswith('__total') or k.endswith('__white_axis') or k.endswith('__remaining_hue')):
                        close(calc_gen[col],summary['genus_equal'][col]['estimate'],case+':'+str(n)+':'+scheme+'::'+col+'::genus')
                        close(calc_sp[col],summary['species_equal'][col]['estimate'],case+':'+str(n)+':'+scheme+'::'+col+'::species')
                    if case=='soil_all' and n==2 and scheme=='buffer500':
                        led=pd.read_csv(results/'soil_all_min2_buffer500_species_losses.csv.gz',compression='gzip').set_index('inat_taxon_id')
                        for col in (k for k in frame.columns if k.endswith('__total') or k.endswith('__white_axis') or k.endswith('__remaining_hue')):
                            compare=frame.set_index('inat_taxon_id')[col].reindex(led.index)
                            close(compare,led[col],col+':per-species')
                        off=pd.DataFrame(offsets,columns=['genus','species','v']).groupby('genus').v.mean().mean()
                        analytic_null_mean=float(off)
                        close(analytic_null_mean,out['null']['primary_null_mean'],'permutation theoretical mean',tol=.0003)
                        # Check independent paired 1,999 genus bootstrap using exactly the frozen resampling seed.
                        rng=np.random.default_rng(summary['seed'])
                        ng=len(gen)
                        ix=rng.integers(ng,size=(1999,ng))
                        m=gen['environment_nonlinear__total'].to_numpy()
                        boot=m[ix].mean(1)
                        close(np.quantile(boot,[.025,.975]),summary['genus_equal']['environment_nonlinear__total']['CI95'],'primary bootstrap CI')
                    cases[case]['min'+str(n)][scheme]={'species':len(frame),'genera':len(gen),'photos':int(frame.n_photos.sum())}
    # The exact null draws come from an independent fixed count-shuffling implementation,
    # but only their probability/support is reported here; the main code is hash-locked.
    assert out['null']['reconstruction_status']=='ALL_SOURCE_COUNTS_EXACT_AND_PERMUTATION_INVENTORY_PRESERVED'
    assert out['null']['null_reps']==1999
    receipt={'status':'PASS','checks':checks,'max_exact_absolute_difference':maxdiff,'null_expectation_monte_carlo_gap':null_mc_gap,
      'source_zip_sha256':h,'cases':cases,'analytic_expected_matched_null_gain':analytic_null_mean,
      'permutation_observed_mean':out['null']['primary_null_mean'],
      'primary_observed_gain':out['null']['primary_observed_environment_gain'],
      'primary_permutation_p':out['null']['primary_plus_one_upper_tail_p'],
      'meaning':'Independent per-species sector demeaning and all point estimates verified; primary fixed bootstrap interval verified; null expectation independently derived.'}
    return receipt

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--results',type=Path,required=True)
    a=ap.parse_args();r=verify(a.source,a.results)
    (a.results/'independent_verification.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ('status','checks','max_exact_absolute_difference','null_expectation_monte_carlo_gap','analytic_expected_matched_null_gain','permutation_observed_mean')}))