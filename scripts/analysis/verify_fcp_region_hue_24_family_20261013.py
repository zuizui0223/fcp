#!/usr/bin/env python3
"""Independent point-estimate reader for original-photo 4-sector x 6-colour pairs.

Uses taxon-level weighted design/lstsq, never imports primary audit or
its genus sufficient-statistics bootstrap; only archived source reader and
original colour classes are reused from the earlier independent verifier.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_fcp_same_genus_east_elsewhere_rain_tradeoff_20261011 as old


def fit_direct(photo, side,mu,sd,equal_genus):
    terms=list(old.COLS)
    cols=list(old.HUES)
    data=photo.loc[photo['inside'].eq(side)]
    if data.empty: raise ValueError('empty side')
    x=data.groupby('inat_taxon_id',as_index=False).agg({**{k:'mean' for k in terms+cols},'genus':'first'})
    x[terms]=(x[terms].to_numpy(dtype=float)-mu)/sd
    x[terms]=x[terms]-x.groupby('genus')[terms].transform('mean')
    x[cols]=x[cols]-x.groupby('genus')[cols].transform('mean')
    wt=x.groupby('genus').inat_taxon_id.transform('size').to_numpy(float)
    w=1./wt if equal_genus else np.ones(len(x))
    # Least squares, not the production summed sufficient statistics.
    XX=x[terms].to_numpy(float)*np.sqrt(w)[:,None]
    YY=x[cols].to_numpy(float)*np.sqrt(w)[:,None]
    return np.linalg.lstsq(XX,YY,rcond=None)[0][0,:]


def verify(source_path,result_path):
    src=old.get_source(source_path)
    _,data,_,mu,sd=old.eligible(src)
    data['sector']=np.floor((data.longitude+180)/90).astype(int).clip(0,3)
    stored=json.loads(result_path.read_text())
    counter=0
    checks=[]
    worst=0.0
    for source_type in ('all_classified_original_source_photos','original_photo_quality_strict'):
        photo=data
        if source_type=='original_photo_quality_strict':
            photo=data[(data.flower_effective_pixels>=data.flower_effective_pixels.median())&
                       (data.nuisance_pixel_fraction<=data.nuisance_pixel_fraction.median())]
        for region in range(4):
            frame=photo.copy()
            frame['inside']=frame.sector.eq(region)
            count=frame.groupby(['genus','inside']).inat_taxon_id.nunique().unstack(fill_value=0)
            names=count.index[(count[True]>=3)&(count[False]>=3)]
            joined=frame[frame.genus.isin(names)]
            for weight_name in ('equal_species','equal_genera'):
                wt=weight_name=='equal_genera'
                a=fit_direct(joined,True,mu,sd,wt)
                b=fit_direct(joined,False,mu,sd,wt)
                parent=stored['cases'][source_type]['weighting'][weight_name]
                raw=parent['sector_metadata'][region]
                for k,h in enumerate(old.HUES):
                    d=max(abs(a[k]-raw['four_hue_inside_slopes_per_sd'][h]),
                          abs(b[k]-raw['four_hue_outside_slopes_per_sd'][h]))
                    worst=max(worst,float(d));counter+=2
                    if d>1e-9:raise AssertionError(f'coefficient mismatch {region} {weight_name} {h} {d}')
                rows=[v for v in parent['all_directional_contrasts'] if v['sector']==region]
                for i,j in itertools.combinations(range(4),2):
                    target=[v for v in rows if v['hue_A']==old.HUES[i] and v['hue_B']==old.HUES[j]]
                    assert len(target)==1
                    got=100*((a[i]-b[i])-(a[j]-b[j]))
                    err=abs(got-target[0]['point_pp_per_BIO12_SD'])
                    worst=max(worst,float(err));counter+=1
                    if err>1e-7:raise AssertionError(f'colour difference mismatch {region}: {i}/{j}')
                checks.append({'case':source_type,'sector':region,'weight':weight_name,
                               'genera':len(names),'species':joined.inat_taxon_id.nunique(),
                               'photos':len(joined)})
    if counter!=224:raise AssertionError(f'Expected 224 checks: got {counter}')
    assert stored['source_SHA256']==old.SHA
    return {'status':'PASS','independent_lstsq_checks':counter,
            'max_abs_original_probability_coeff_or_pp_error':worst,
            'n_confirmed_original_sector_weight_conditions':len(checks),
            'source_archive_sha256':old.SHA,
            'method':'independent original photo -> species averages -> genus centered weighted np.linalg.lstsq',
            'authoritative_bootstrap':'computed by primary shared-global-genus-draw procedure; not independently reproduced here',
            'conditions':checks}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source-archive',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    z=verify(args.source_archive,args.result)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in z.items() if k!='conditions'},ensure_ascii=False,indent=2))
