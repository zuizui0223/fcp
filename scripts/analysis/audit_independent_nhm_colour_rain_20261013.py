#!/usr/bin/env python3
"""Independent botanical-description flower colour vs rainfall direction in FCP geographic opportunities.

This source-exposed follow-up to FCP PR190 tests a pre-fixed East 90-180 E
white-minus-blue/purple photographic *site-climate* association with botanical
colour labels supplied by Natural History Museum (NHM) plant-traits, not iNat
images. This is independent colour measurement, NOT independent environment,
plant fitness, within-species regional morph, true phylogeny, or a blinded
new species selection. Failure of coverage gates is HOLD, not zero effect.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import zipfile
import io
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd

FCP_SHA = '2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
NHM_BLOB_SHA = 'e1cfd2d5c518c2061ec91d1a9e8b64fc7ca69389'
NHM_BYTES = 17_279_412
NHM_REPO_COMMIT = '4843c22bb2fc91018d9828dcddebf067c83b1097'
DOI = '10.5519/p3dm31kc'
COLS = ('wc_bio12','wc_bio1','wc_bio5','wc_bio15','wc_elevation_m','abs_lat','lon_sin','lon_cos')
ORGANS = ('petal colour','corolla colour','flower colour')
SEED = 2026101321
B = 999
SECTOR_NAMES = ('W180_W90','W90_E0','E0_E90','E90_E180')
# Recognised source description words, not inferred pigments; ambiguous mixed species HOLD.
WHITE = {'white'}
BLUE = {'blue','purple','violet'}
OTHER_HUE = {'red','pink','magenta','yellow','orange','green','brown','black','cream','ivory','scarlet','crimson','maroon','lilac','lavender','grey','gray','bronze','gold','golden','purple-red','blue-violet'}
WORD_RE = re.compile(r'[a-z]+')
NAME_RE = re.compile(r'^[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+$')


def blob_sha(path: Path) -> str:
    b=path.read_bytes()
    return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\x00'+b).hexdigest()


def classify_botanical_description(rec: dict) -> tuple[str,str]:
    """Require all nonempty floral organs to have exactly ONE shared hue family.

    White+blue, white+other, blue+other, and any unrecognized nonempty field
    are HOLD, not a presumed clean photo white/blue label.
    """
    detected=[]; reasons=[]
    for organ in ORGANS:
        raw=rec.get(organ)
        if raw is None or pd.isna(raw) or not str(raw).strip(): continue
        txt=str(raw).lower()
        tokens=set(WORD_RE.findall(txt))
        found=set()
        if tokens&WHITE:found.add('white')
        if tokens&BLUE:found.add('blue_purple')
        if tokens&OTHER_HUE:found.add('other')
        if len(found)!=1:
            return 'HOLD','ambiguous_or_unrecognized_'+organ.replace(' ','_')
        detected.append(next(iter(found)))
        reasons.append(organ)
    if not detected: return 'HOLD','no_unambiguous_flower_organ_colour'
    if len(set(detected))!=1: return 'HOLD','cross_organ_colour_conflict'
    if detected[0]=='other': return 'OTHER','other_colour_not_used_in_white_blue_target'
    return detected[0], '|'.join(reasons)


def extract_botanical(nhm: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    needed=('taxon',)+ORGANS
    missing=[k for k in needed if k not in nhm.columns]
    if missing: raise ValueError('NHM raw fields missing: '+str(missing))
    x=nhm.loc[:,needed].copy()
    x['taxon']=x['taxon'].astype(str).str.strip().str.replace(r'\s+',' ',regex=True)
    # Reject hybrids, trinomial, uncertain taxa. No fuzzy repair.
    valid=x['taxon'].map(lambda s: bool(NAME_RE.fullmatch(s)))
    x=x.loc[valid].copy()
    results=x.apply(lambda row: classify_botanical_description(row.to_dict()),axis=1)
    x['label']=[v[0] for v in results]
    x['reason']=[v[1] for v in results]
    # Multiple independent source rows for same binomial are allowed only
    # if every row gives the same clear label.
    same=x.groupby('taxon').label.nunique()
    conflicts=set(same[same>1].index)
    if conflicts:x.loc[x.taxon.isin(conflicts),'label']='HOLD'
    summary={'source_rows':int(len(nhm)),'strict_binomial_rows':int(len(x)),
      'strict_binomial_unique_names':int(x.taxon.nunique()),
      'botanical_label_rows':x.label.value_counts().to_dict(),
      'ambiguous_conflicted_duplicate_names':len(conflicts)}
    unique=x.loc[~x.taxon.duplicated(keep='first'),['taxon','label','reason']].copy()
    assert unique.taxon.is_unique
    return unique,summary


def read_fcp_photo_opportunities(path:Path) -> tuple[pd.DataFrame,dict]:
    if hashlib.sha256(path.read_bytes()).hexdigest()!=FCP_SHA:
        raise ValueError('original photo opportunity SHA256 mismatch')
    with zipfile.ZipFile(path) as z:
        d='fcp42111_real_site_abiotic/'
        def read(n):return pd.read_csv(io.BytesIO(z.read(d+n+'.csv.gz')),compression='gzip',low_memory=False)
        a=read('breadth_flower_colour_climate_soil')
        b=read('taxon_cell_flower_colour_climate_soil')
    if len(a)!=42111 or len(b)!=85337:raise ValueError('FCP historical opportunity cohort drift')
    x=pd.concat([a,b],ignore_index=True)
    if x.photo_id.duplicated().sum()!=26905:raise ValueError('overlap photo ID count drift')
    # Full original includes 54,590 photos FCP was unable to classify;
    # these stay admitted because phenotype comes from NHM descriptions.
    t=x.drop_duplicates('photo_id').copy()
    if len(t)!=100543 or t.inat_taxon_id.nunique()!=42111:
        raise ValueError('100543 photo/42111 source taxon support drift')
    subset=t[['photo_id','species','inat_taxon_id','latitude','longitude',
      'wc_bio12','wc_bio1','wc_bio5','wc_bio15','wc_elevation_m']].copy()
    climate=subset[list(COLS[:5])+['latitude','longitude']].notna().all(axis=1)
    geo=subset.latitude.between(-90,90)&subset.longitude.between(-180,180)
    subset=subset.loc[climate&geo].copy()
    subset['species']=subset.species.astype(str).str.strip().str.replace(r'\s+',' ',regex=True)
    subset['genus']=subset.species.str.split().str[0]
    subset['abs_lat']=subset.latitude.abs()
    subset['lon_sin']=np.sin(np.deg2rad(subset.longitude))
    subset['lon_cos']=np.cos(np.deg2rad(subset.longitude))
    subset['sector']=np.minimum(3,np.floor((subset.longitude+180)/90).astype(int))
    return subset,{'source_photo_ids':len(t),'source_taxa':42111,
      'source_full_photo_site_climate_opportunities':len(subset),
      'eligible_species_all_photo_climate':int(subset.species.nunique())}


def species_means(rows:pd.DataFrame,means:np.ndarray,sd:np.ndarray) ->pd.DataFrame:
    agg={v:'mean' for v in COLS}; agg.update({'photo_id':'size','genus':'first','label':'first'})
    d=rows.groupby('species',sort=True,as_index=False).agg(agg).rename(columns={'photo_id':'n_photo'})
    d.loc[:,list(COLS)]=(d[list(COLS)].to_numpy(float)-means)/sd
    d['white']=(d.label=='white').astype(float)
    d['blue_purple']=(d.label=='blue_purple').astype(float)
    assert np.allclose(d.white+d.blue_purple,1)
    return d


def shared_genera(photo:pd.DataFrame,min_species=3):
    valid=photo.copy()
    valid['east']=valid.sector.eq(3)
    by=valid.groupby(['genus','east']).species.nunique().unstack(fill_value=0)
    if not {True,False}.issubset(by.columns):return valid.iloc[:0].copy()
    selected=by.index[(by[True]>=min_species)&(by[False]>=min_species)]
    return valid.loc[valid.genus.isin(selected)].copy()


def regression_matrices(species:pd.DataFrame,equal_genus:bool) -> dict:
    mats={}
    for g,x in species.groupby('genus',sort=True):
        a=x[list(COLS)].to_numpy(float); y=x[['white','blue_purple']].to_numpy(float)
        a=a-a.mean(0);y=y-y.mean(0)
        if not np.isfinite(a).all() or not np.isfinite(y).all():raise ValueError('invalid source N/S')
        n=len(x); weight=1/n if equal_genus else 1.0
        mats[g]=(weight*(a.T@a),weight*(a.T@y),n)
    return mats


def fit_sum(a:dict, names:list[str]):
    xx=np.stack([a[g][0] for g in names]).sum(0)
    yy=np.stack([a[g][1] for g in names]).sum(0)
    if np.linalg.cond(xx)>1e10:return None
    return np.linalg.solve(xx,yy)[0,:]


def paired_bootstrap(groups_e:dict, groups_o:dict,nreps:int=B,seed:int=SEED):
    names=sorted(set(groups_e)&set(groups_o))
    n=len(names)
    a,b=fit_sum(groups_e,names),fit_sum(groups_o,names)
    if a is None or b is None:return {'status':'HOLD_DESIGN_RANK','n_genus':n}
    region_contrast=lambda A,B:(A[0]-A[1])-(B[0]-B[1])
    point=float(region_contrast(a,b))
    AE=np.stack([groups_e[g][0] for g in names]); AO=np.stack([groups_o[g][0] for g in names])
    BE=np.stack([groups_e[g][1] for g in names]); BO=np.stack([groups_o[g][1] for g in names])
    rng=np.random.default_rng(seed)
    weights=rng.multinomial(n,np.repeat(1/n,n),size=nreps)
    re=np.einsum('rn,nij->rij',weights,AE,optimize=True)
    ro=np.einsum('rn,nij->rij',weights,AO,optimize=True)
    be=np.einsum('rn,nij->rij',weights,BE,optimize=True)
    bo=np.einsum('rn,nij->rij',weights,BO,optimize=True)
    try:
        pe=np.linalg.solve(re,be)[:,0,:]
        po=np.linalg.solve(ro,bo)[:,0,:]
    except np.linalg.LinAlgError:return {'status':'HOLD_SINGULAR_BOOTSTRAP','n_genus':n}
    vals=(pe[:,0]-pe[:,1])-(po[:,0]-po[:,1])
    return {'status':'OK','n_genera':n,'point':point,'CI95':np.quantile(vals,[.025,.975]).tolist(),
        'bootstrap_frac_positive':float(np.mean(vals>0)),
        'bootstrap_samples':vals.tolist()}


def single_sector(photo:pd.DataFrame,means:np.ndarray,sd:np.ndarray,sector:int,minimum:int,genus_equal:bool,boot:int):
    photo=photo.copy(); photo['east']=photo.sector.eq(sector)
    a=shared_genera(photo.assign(sector=np.where(photo['east'],3,0)),minimum)
    genus=a.genus.nunique()
    sp=a.groupby('east').species.nunique()
    labels=a.groupby('east').apply(lambda g:g.groupby('species',sort=False).label.first().value_counts().to_dict(),include_groups=False).to_dict() if len(a) else {}
    support={'n_genus':int(genus),'source_photos':int(len(a)),
      'species_by_side':{str(k):int(v) for k,v in sp.items()},
      'named_species':int(a.species.nunique()),'botanical_colour_by_side':{str(k):v for k,v in labels.items()},
      'minimum_species_each_side':minimum}
    if genus<30 or min([sp.get(True,0),sp.get(False,0)])<100 or any(labels.get(side,{}).get(col,0)<20 for side in (True,False) for col in ('white','blue_purple')):
        return {'status':'HOLD_INDEPENDENT_BOTANICAL_COLOUR_SUPPORT','support':support}
    yes=species_means(a.loc[a.east],means,sd)
    no=species_means(a.loc[~a.east],means,sd)
    res=paired_bootstrap(regression_matrices(yes,genus_equal),regression_matrices(no,genus_equal),nreps=boot,seed=SEED+sector+(100 if genus_equal else 0)+10*minimum)
    return {'status':res['status'],'support':support,'estimate':res}


def run(nhm_file:Path,fcp_file:Path,bootstrap_n:int=B):
    if nhm_file.stat().st_size!=NHM_BYTES or blob_sha(nhm_file)!=NHM_BLOB_SHA:
        raise ValueError('NHM published independent description-source commit/blob hash mismatch')
    nhm=pd.read_csv(nhm_file,low_memory=False)
    botan,botanical_report=extract_botanical(nhm)
    source,fcp_report=read_fcp_photo_opportunities(fcp_file)
    means=source[list(COLS)].mean().to_numpy(float)
    sd=source[list(COLS)].std(ddof=0).to_numpy(float)
    if (sd<=0).any():raise ValueError('degenerate exposure scale')
    # No FCP measurement_status/morph/palette used in independent outcome.
    outcome=botan.loc[botan.label.isin(['white','blue_purple']),['taxon','label']].copy()
    x=source.merge(outcome,left_on='species',right_on='taxon',how='inner',validate='many_to_one')
    # Independent site exposure for all separately botanically described species.
    if x.empty:raise ValueError('no exact taxonomic overlap')
    result={'schema':'FCP_NHM_INDEPENDENT_BOTANICAL_DESCRIPTION_COLOUR_REPLICATION_V1',
      'design':'post-FCP-photo-hypothesis frozen independent botanical source; same FCP source site-climate, NOT new environmental independent source',
      'date':'2026-10-13','original_fcp_photo_SHA256':FCP_SHA,'NHM_source_repository':'NaturalHistoryMuseum/plant-traits',
      'NHM_source_commit':NHM_REPO_COMMIT,'NHM_source_git_blob_sha1':NHM_BLOB_SHA,'NHM_dataset_DOI':DOI,
      'original_photo_opportunity':fcp_report,'independent_botanical_coverage':botanical_report,
      'overlap':{'original_photo_rows':len(x),'botanically_labeled_species':int(x.species.nunique()),
         'botanically_labeled_genera':int(x.genus.nunique()),
         'botanical_whites':int((outcome.label=='white').sum()),
         'botanical_blue_purples':int((outcome.label=='blue_purple').sum()),
         'overlap_white_species':int(x.loc[x.label.eq('white'),'species'].nunique()),
         'overlap_blue_species':int(x.loc[x.label.eq('blue_purple'),'species'].nunique())},
      'spatial_strata':{},'hard_boundaries':['Pure description petal/corolla/flower colour categories only, mixed and conflicting organs excluded.',
      'Colour label is species-wide botany from textual natural history records; cannot estimate within-species geographic hue-frequency changes.',
      'The environment is from the same old iNaturalist site opportunities and source could bias geographic availability.',
      'No climatic fitness, plant pigment chemistry, phylogenetic covariance, taxonomy synonym expansion, pollinator or selection mechanism.',
      'No historic FCP numeric records or original manuscript H1/H2 changed.']}
    # Four sectors for an already selected east 90-180 test; all four included.
    # All 30-genus gates fixed before examining independent colour outcomes.
    for minimum in (3,2):
        for strict in (False,True):
            # strict photo means >=2 source photo sites per species: independent
            # of PHOTO measured hue and thus no image-quality selection.
            xx=x
            if strict:
                ids=x.groupby('species').photo_id.nunique()
                xx=x[x.species.isin(ids[ids>=2].index)]
            cond=f'min{minimum}_'+('2plus_photo_opportunity' if strict else 'all_photo_opportunity')
            result['spatial_strata'][cond]={}
            for weighting in ('equal_species','equal_genera'):
                results=[]
                for s in range(4):
                    v=single_sector(xx,means,sd,s,minimum,weighting=='equal_genera',bootstrap_n)
                    v['sector']=SECTOR_NAMES[s]
                    # keep individual bootstrap arrays only internally in result JSON for replicate exactness
                    results.append(v)
                # all four must be OK for family-wise inference (no winner-only family)
                all_pass=all(i['status']=='OK' for i in results)
                if all_pass:
                    diff=np.column_stack([q['estimate']['bootstrap_samples'] for q in results])
                    point=np.array([q['estimate']['point'] for q in results])
                    sdv=diff.std(0,ddof=1)
                    if (sdv<=1e-12).any():
                        all_pass=False;crit=None
                    else:
                        crit=float(np.quantile(np.max(np.abs((diff-point[None,:])/sdv[None,:]),axis=1),.95))
                        for k,q in enumerate(results):
                            q['estimate']['maxT_4sector_CI95']=[float(point[k]-crit*sdv[k]),float(point[k]+crit*sdv[k])]
                            q['estimate']['maxT_family_4_sectors']=4
                result['spatial_strata'][cond][weighting]={
                    'status':'OK_4SECTOR_FAMILY' if all_pass else 'HOLD_COMPLETE_4SECTOR_FAMILY',
                    'maxT_4sector_critical95':crit if all_pass else None,
                    'sectors':results}
    # avoid 8*4*999 samples bloating deliverable; retain exact CI only
    for strata in result['spatial_strata'].values():
        for w in strata.values():
            for s in w['sectors']:
                if 'estimate' in s:s['estimate'].pop('bootstrap_samples',None)
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--nhm-csv',type=Path,required=True)
    ap.add_argument('--fcp-original-source',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--bootstrap',type=int,default=B)
    args=ap.parse_args()
    z=run(args.nhm_csv,args.fcp_original_source,args.bootstrap)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/'result.json').write_text(json.dumps(z,indent=2,ensure_ascii=False)+'\n')
    summary={'original_overlap':z['overlap'],'cases':{cond:{w:{'status':v['status'],
      'n_pass_sectors':sum(s['status']=='OK' for s in v['sectors']),
      'east':{k:y for k,y in v['sectors'][3].items() if k!='estimate'},
      'east_result':{k:y for k,y in v['sectors'][3].get('estimate',{}).items() if k!='bootstrap_samples'}}
      for w,v in p.items()} for cond,p in z['spatial_strata'].items()}}
    (args.outdir/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
