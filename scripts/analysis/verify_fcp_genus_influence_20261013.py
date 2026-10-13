#!/usr/bin/env python3
"""Independent audit: direct species-row weighted LS (does NOT import production audit)."""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np,pandas as pd

EXPECTED='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
RAW=('wc_bio12','wc_bio1','wc_bio5','wc_bio15','wc_elevation_m')
COV=RAW+('abs_lat','lon_sin','lon_cos')
COLOR=('white','yellow_orange','red_pink','blue_purple')


def read_raw(p):
    if hashlib.sha256(p.read_bytes()).hexdigest()!=EXPECTED:raise AssertionError('source SHA changed')
    with zipfile.ZipFile(p) as z:
        def grab(x):return pd.read_csv(io.BytesIO(z.read('fcp42111_real_site_abiotic/'+x+'.csv.gz')),compression='gzip',low_memory=False)
        allr=pd.concat((grab('taxon_cell_flower_colour_climate_soil'),grab('breadth_flower_colour_climate_soil')),ignore_index=True)
    assert len(allr)==127448 and allr.photo_id.duplicated().sum()==26905
    allr=allr.drop_duplicates('photo_id')
    assert len(allr)==100543
    mask=allr.measurement_status.eq('classified_four_state_morph') & allr.morph.isin(COLOR)
    mask &= allr[list(RAW)+['latitude','longitude']].notna().all(axis=1)
    mask &= allr.latitude.between(-90,90)&allr.longitude.between(-180,180)
    a=allr.loc[mask].copy()
    assert len(a)==45836 and a.inat_taxon_id.nunique()==25028
    a['abs_lat']=a.latitude.abs()
    a['lon_sin']=np.sin(np.deg2rad(a.longitude))
    a['lon_cos']=np.cos(np.deg2rad(a.longitude))
    a['genus']=a.species.astype(str).str.split().str[0]
    a['east']=a.longitude.ge(90.)&a.longitude.le(180.)
    for k in COLOR:a[k]=a.morph.eq(k).astype(float)
    return a,a[list(COV)].mean().to_numpy(),a[list(COV)].std(ddof=0).to_numpy()


def direct_slope(data,mu,sd,equal_genus,remove):
    d=data[~data.genus.isin(remove)]
    betas=[]
    for inside in (True,False):
        t=d[d.east.eq(inside)]
        v=t.groupby(['inat_taxon_id','genus'],as_index=False)[list(COV+COLOR)].mean()
        v[list(COV)]=(v[list(COV)].to_numpy(float)-mu)/sd
        xx=v[list(COV)].to_numpy(float)
        yy=v[list(COLOR)].to_numpy(float)
        gen,idx=np.unique(v.genus.to_numpy(str),return_inverse=True)
        sz=np.bincount(idx)
        for j in range(8):
            xxx=np.bincount(idx,weights=xx[:,j])/sz
            xx[:,j]-=xxx[idx]
        for j in range(4):
            yyy=np.bincount(idx,weights=yy[:,j])/sz
            yy[:,j]-=yyy[idx]
        if equal_genus:
            w=1/np.sqrt(sz[idx])
            xx=xx*w[:,None]
            yy=yy*w[:,None]
        coef=np.linalg.lstsq(xx,yy,rcond=None)[0][0]
        betas.append(coef)
    E,O=betas
    return float(100*((E[0]-E[3])-(O[0]-O[3])))


def verify(source_path,result_json,influence_csv):
    src,mu,sd=read_raw(source_path)
    result=json.loads(Path(result_json).read_text());influence=pd.read_csv(influence_csv)
    checks=0;deviations=[]
    for sample,cond in [('all',src),('strict',src.loc[(src.flower_effective_pixels>=src.flower_effective_pixels.median())&(src.nuisance_pixel_fraction<=src.nuisance_pixel_fraction.median())])]:
        groups=cond.groupby(['genus','east']).inat_taxon_id.nunique().unstack(fill_value=0)
        eligible=set(groups.index[(groups[True]>=3)&(groups[False]>=3)])
        use=cond[cond.genus.isin(eligible)]
        for nm,eq in [('species_equal',False),('genus_equal',True)]:
            key=sample+'_'+nm
            a=result['results'][key]
            if a['n_shared_genera']!=len(eligible):raise AssertionError('Genus support mismatch')
            check=[('full',(),a['full_C_pp_per_sd'])]
            if len(a['top_genus_by_absolute_impact'])>=3:
                check.extend((f'omitted_{q["genus"]}',(q['genus'],),q['omit_C_pp']) for q in a['top_genus_by_absolute_impact'][:3])
            rich=influence[(influence.source_sample==sample)&(influence.weight==nm)].sort_values(['n_total_species_units','genus'],ascending=[False,True])['genus'].tolist()
            for k in (1,5,10,20):
                if str(k) in a['predefined_remove_largest_species_rich_genera_C_pp']:
                    check.append((f'richest_{k}',tuple(rich[:k]),a['predefined_remove_largest_species_rich_genera_C_pp'][str(k)]))
            for label,removed,expected in check:
                value=direct_slope(use,mu,sd,eq,set(removed))
                delta=abs(value-expected)
                deviations.append(delta)
                if delta>1e-9:raise AssertionError((sample,nm,label,value,expected))
                checks+=1
    return {'status':'PASS','source_verified':True,'independent_direct_species_row_least_squares_checks':checks,
            'max_abs_difference_pp':max(deviations),'no_production_import':True}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-archive',type=Path,required=True)
    ap.add_argument('--result',type=Path,required=True)
    ap.add_argument('--genus-ledger',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    v=verify(args.source_archive,args.result,args.genus_ledger)
    args.output.write_text(json.dumps(v,indent=2)+'\n')
    print(v)
if __name__=='__main__':main()
