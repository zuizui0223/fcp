#!/usr/bin/env python3
"""Independent raw source, literal botanical lexicon and direct species-row least-squares check.

This deliberately does not import the primary analysis module; bootstrap
uncertainty is NOT independently regenerated, but fit point estimates,
eligibility counts, duplicate rules, and all original-source provenance are.
"""
from __future__ import annotations
import argparse, hashlib, json, zipfile, io, re
from pathlib import Path
import numpy as np
import pandas as pd

FCP_SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
NHM_SHA='e1cfd2d5c518c2061ec91d1a9e8b64fc7ca69389'
COV=('wc_bio12','wc_bio1','wc_bio5','wc_bio15','wc_elevation_m','abs_lat','lon_sin','lon_cos')
FIELD=('petal colour','corolla colour','flower colour')
WHITE={'white'}; BLUE={'blue','purple','violet'}
OTHER={'red','pink','magenta','yellow','orange','green','brown','black','cream','ivory','scarlet','crimson','maroon','lilac','lavender','grey','gray','bronze','gold','golden','purple-red','blue-violet'}


def independent_trait(raw):
    assigned=[]
    for f in FIELD:
        v=raw.get(f)
        if v is None or pd.isna(v) or not str(v).strip():continue
        t=set(re.findall(r'[a-z]+',str(v).lower()))
        choice=[]
        if t&WHITE:choice.append('white')
        if t&BLUE:choice.append('blue_purple')
        if t&OTHER:choice.append('other')
        if len(choice)!=1:return 'HOLD'
        assigned.append(choice[0])
    return assigned[0] if assigned and len(set(assigned))==1 else 'HOLD'


def raw_description_table(path):
    b=path.read_bytes()
    expected_git_blob=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\x00'+b).hexdigest()
    if expected_git_blob!=NHM_SHA:raise ValueError('NHM Git blob hash changed')
    orig=pd.read_csv(path,low_memory=False)
    records=[]
    for row in orig.loc[:,['taxon',*FIELD]].to_dict(orient='records'):
        taxon=re.sub(r'\s+',' ',str(row['taxon']).strip())
        if not re.fullmatch(r'[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+',taxon):continue
        records.append((taxon,independent_trait(row)))
    df=pd.DataFrame(records,columns=['taxon','label'])
    vals=df.groupby('taxon').label.nunique()
    conflicts=set(vals.index[vals>1])
    df.loc[df.taxon.isin(conflicts),'label']='HOLD'
    return df.drop_duplicates('taxon',keep='first')


def raw_locations(path):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=FCP_SHA:raise ValueError('FCP full-photo SHA mismatch')
    with zipfile.ZipFile(path) as z:
        directory='fcp42111_real_site_abiotic/'
        a=pd.read_csv(io.BytesIO(z.read(directory+'breadth_flower_colour_climate_soil.csv.gz')),compression='gzip',low_memory=False)
        b=pd.read_csv(io.BytesIO(z.read(directory+'taxon_cell_flower_colour_climate_soil.csv.gz')),compression='gzip',low_memory=False)
    if len(a)!=42111 or len(b)!=85337:raise ValueError('FCP photo source counts differ')
    subset=pd.concat([a,b],ignore_index=True).drop_duplicates('photo_id')
    if len(subset)!=100543:raise ValueError('FCP original ID shape unexpected')
    loc=subset[['photo_id','species','latitude','longitude','wc_bio1','wc_bio5','wc_bio12','wc_bio15','wc_elevation_m']].copy()
    mask=loc[['latitude','longitude','wc_bio1','wc_bio5','wc_bio12','wc_bio15','wc_elevation_m']].notna().all(axis=1)
    mask=mask&loc.latitude.between(-90,90)&loc.longitude.between(-180,180)
    loc=loc.loc[mask].copy()
    loc['species']=loc.species.astype(str).str.strip().str.replace(r'\s+',' ',regex=True)
    loc['genus']=loc.species.str.split().str[0]
    loc['abs_lat']=np.abs(loc.latitude)
    loc['lon_sin']=np.sin(np.deg2rad(loc.longitude))
    loc['lon_cos']=np.cos(np.deg2rad(loc.longitude))
    loc['sector']=np.minimum(3,np.floor((loc.longitude+180)/90).astype(int))
    return loc


def fit_direct_species_rows(data,equal_genus):
    # Unlike original block-matrix estimator, explicitly residualise every
    # species predictor+binary response and fit a dense weighted design.
    out=[]
    for genus,g in data.groupby('genus',sort=True):
        x=g[list(COV)].to_numpy(float)
        y=(g.label=='white').astype(float).to_numpy()
        x=x-x.mean(0)
        y=y-y.mean()
        n=len(g)
        w=np.full(n,1.0/n if equal_genus else 1.0)
        out.append((x*np.sqrt(w)[:,None],y*np.sqrt(w)))
    X=np.concatenate([i[0] for i in out],axis=0)
    Y=np.concatenate([i[1] for i in out])
    beta,_,rank,_=np.linalg.lstsq(X,Y,rcond=None)
    return None if rank<8 else float(2*beta[0])


def verify(nhm_csv,fcp_zip,stored):
    nhm=raw_description_table(nhm_csv)
    fcp=raw_locations(fcp_zip)
    unlabelled=fcp[list(COV)].to_numpy(float)
    mean=unlabelled.mean(0);sd=unlabelled.std(0)
    lab=nhm[nhm.label.isin(['white','blue_purple'])]
    d=fcp.merge(lab,left_on='species',right_on='taxon',how='inner',validate='many_to_one')
    assert len(d)==stored['overlap']['original_photo_rows']
    assert d.species.nunique()==stored['overlap']['botanically_labeled_species']
    assert d.genus.nunique()==stored['overlap']['botanically_labeled_genera']
    for k in COV:d[k]=(d[k]-mean[COV.index(k)])/sd[COV.index(k)]
    nfit=0;hold=0;max_diff=0;count=0
    for key,all_weights in stored['spatial_strata'].items():
        min_req=int(key[3])
        selected=d.copy()
        if '2plus_' in key:
            count_photo=selected.groupby('species').photo_id.nunique()
            selected=selected.loc[selected.species.isin(count_photo.index[count_photo>=2])]
        for weight_name,weight_data in all_weights.items():
            for num,case in enumerate(weight_data['sectors']):
                # group source taxa after all photo opportunities selected
                frame=selected.copy();frame['selected_region']=frame.sector.eq(num)
                counts=frame.groupby(['genus','selected_region']).species.nunique().unstack(fill_value=0)
                selected_genera=counts.index[(counts.get(True,0)>=min_req)&(counts.get(False,0)>=min_req)]
                frame=frame.loc[frame.genus.isin(selected_genera)].copy()
                if int(frame.genus.nunique())!=case['support']['n_genus']:
                    raise AssertionError(f'Independent species coverage mismatch: {key} {num}')
                count+=1
                if case['status']!='OK':hold+=1;continue
                e=frame.loc[frame.selected_region];o=frame.loc[~frame.selected_region]
                def agg(df):
                    grouped=df.groupby('species',as_index=False).agg({**{k:'mean' for k in COV},'genus':'first','label':'first'})
                    return grouped
                be=fit_direct_species_rows(agg(e),weight_name=='equal_genera')
                bo=fit_direct_species_rows(agg(o),weight_name=='equal_genera')
                if be is None or bo is None:raise AssertionError('Independent rank failure in accepted source')
                actual=be-bo
                expected=case['estimate']['point']
                max_diff=max(max_diff,abs(actual-expected))
                if abs(actual-expected)>1e-9:
                    raise AssertionError(f'Independent original-photo regression mismatch {key}/{weight_name}/{num}: {actual} vs {expected}')
                nfit+=1
    return {'schema':'fcp_nhm_botanical_independent_direct_source_verification_v1','status':'PASS',
      'independent_recomputed_regression_contrasts':nfit,'held_non_fit_sectors':hold,
      'verified_foursector_weighting_support_cells':count,'max_abs_estimate_difference':max_diff,
      'source_nhm_exact_colour_records':len(nhm),'source_fcp_eligible_sites':len(fcp),'matched_original_photo_rows':len(d),
      'never_uses_fcp_photo_hue_outcomes':True,'independent_bootstrap_CI_replayed':False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--nhm-csv',type=Path,required=True)
    p.add_argument('--fcp-original-source',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    d=verify(a.nhm_csv,a.fcp_original_source,json.loads(a.result.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(d,indent=2))
if __name__=='__main__':main()
