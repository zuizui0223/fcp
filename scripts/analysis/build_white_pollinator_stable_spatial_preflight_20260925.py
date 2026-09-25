#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

EARTH_KM=6371.0088
POLL_TO_PLANT={'pollinates','visitsFlowersOf'}
PLANT_TO_POLL={'pollinatedBy','flowersVisitedBy'}
ITYPES=POLL_TO_PLANT|PLANT_TO_POLL
BEE_MARKERS=('apoidea','anthophila','apidae','megachilidae','halictidae','andrenidae','colletidae','melittidae','stenotritidae')

def md5(path: Path) -> str:
    h=hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def classify(path,name):
    s=(str(path)+' '+str(name)).lower()
    if 'sphingidae' in s:
        return 'lepidoptera','hawkmoth'
    if 'lepidoptera' in s:
        return 'lepidoptera','other_lepidoptera'
    if any(x in s for x in BEE_MARKERS):
        return 'bee','bee'
    if 'hymenoptera' in s:
        return 'other_hymenoptera','other_hymenoptera'
    if 'diptera' in s:
        return 'diptera','diptera'
    if 'coleoptera' in s:
        return 'coleoptera','coleoptera'
    if 'aves' in s or 'trochilidae' in s or 'nectariniidae' in s or 'meliphagidae' in s:
        return 'aves','aves'
    if 'chiroptera' in s or 'pteropodidae' in s or 'phyllostomidae' in s:
        return 'chiroptera','chiroptera'
    return 'other_unknown','other_unknown'

def fnum(x):
    try:
        v=float(x)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan

def xyz(lat,lon):
    lat=np.deg2rad(np.asarray(lat,float)); lon=np.deg2rad(np.asarray(lon,float))
    c=np.cos(lat)
    return np.column_stack([c*np.cos(lon),c*np.sin(lon),np.sin(lat)])

def nearest_km(photo_lat,photo_lon,point_lat,point_lon):
    if len(point_lat)==0:
        return np.full(len(photo_lat),np.nan)
    tree=cKDTree(xyz(point_lat,point_lon))
    chord,_=tree.query(xyz(photo_lat,photo_lon),k=1)
    chord=np.clip(chord,0,2)
    return EARTH_KM*2*np.arcsin(chord/2)

def scan_globi(path: Path, species_set: set[str]):
    csv.field_size_limit(sys.maxsize)
    rows=[]; total=0; poll_rows=0; matched=0
    with gzip.open(path,'rt',encoding='utf-8',errors='replace',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'sourceTaxonName','sourceTaxonPathNames','interactionTypeName','targetTaxonName','targetTaxonPathNames','decimalLatitude','decimalLongitude'}
        missing=required-set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f'missing stable GloBI columns: {sorted(missing)}')
        for r in reader:
            total+=1
            it=(r.get('interactionTypeName') or '').strip()
            if it not in ITYPES:
                continue
            poll_rows+=1
            if it in POLL_TO_PLANT:
                plant=(r.get('targetTaxonName') or '').strip()
                if plant not in species_set:
                    continue
                poll=(r.get('sourceTaxonName') or '').strip()
                ppath=(r.get('sourceTaxonPathNames') or r.get('sourceTaxonPath') or '').strip()
            else:
                plant=(r.get('sourceTaxonName') or '').strip()
                if plant not in species_set:
                    continue
                poll=(r.get('targetTaxonName') or '').strip()
                ppath=(r.get('targetTaxonPathNames') or r.get('targetTaxonPath') or '').strip()
            if not poll:
                continue
            guild,sub=classify(ppath,poll)
            rows.append({
                'plant_species':plant,
                'pollinator_name':poll,
                'pollinator_path':ppath,
                'guild':guild,
                'subguild':sub,
                'interaction_type':it,
                'decimalLatitude':fnum(r.get('decimalLatitude')),
                'decimalLongitude':fnum(r.get('decimalLongitude')),
                'referenceCitation':(r.get('referenceCitation') or '').strip(),
                'namespace':(r.get('namespace') or '').strip(),
            })
            matched+=1
            if total%1000000==0:
                print(json.dumps({'rows_scanned':total,'pollination_rows':poll_rows,'matched_rows':matched}),flush=True)
    d=pd.DataFrame(rows)
    if len(d):
        d=d.drop_duplicates(['plant_species','pollinator_name','interaction_type','decimalLatitude','decimalLongitude','referenceCitation','namespace']).reset_index(drop=True)
    return d,{'rows_scanned':total,'pollination_rows_scanned':poll_rows,'matched_rows_before_dedup':matched,'matched_rows_after_dedup':int(len(d))}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True)
    p.add_argument('--globi-tsv-gz',required=True)
    p.add_argument('--expected-md5',required=True)
    p.add_argument('--outdir',required=True)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    measured=pd.read_csv(a.measured,usecols=['photo_id','species','latitude','longitude'])
    species=sorted(measured.species.dropna().astype(str).unique())
    if len(species)!=499:
        raise SystemExit(f'expected 499 species, found {len(species)}')
    measured=measured.dropna(subset=['latitude','longitude']).copy()
    actual=md5(Path(a.globi_tsv_gz))
    if actual.lower()!=a.expected_md5.lower():
        raise SystemExit(f'GloBI MD5 mismatch: {actual}')
    inter,scan=scan_globi(Path(a.globi_tsv_gz),set(species))
    inter.to_csv(out/'stable_pollinator_interactions.csv.gz',index=False,compression='gzip')

    geo=inter.dropna(subset=['decimalLatitude','decimalLongitude']).copy()
    if len(geo):
        geo['cell1deg']=np.floor(geo.decimalLatitude).astype(int).astype(str)+':'+np.floor(geo.decimalLongitude).astype(int).astype(str)

    dist_parts=[]; species_rows=[]
    for sp,gp in measured.groupby('species'):
        gi=geo[geo.plant_species.eq(sp)].copy()
        all_lat=gi.decimalLatitude.to_numpy(float); all_lon=gi.decimalLongitude.to_numpy(float)
        bee=gi[gi.guild.eq('bee')]
        lep=gi[gi.guild.eq('lepidoptera')]
        hawk=gi[gi.subguild.eq('hawkmoth')]
        q=gp[['photo_id','species','latitude','longitude']].copy()
        q['nearest_any_km']=nearest_km(q.latitude,q.longitude,all_lat,all_lon)
        q['nearest_bee_km']=nearest_km(q.latitude,q.longitude,bee.decimalLatitude.to_numpy(float),bee.decimalLongitude.to_numpy(float))
        q['nearest_lepidoptera_km']=nearest_km(q.latitude,q.longitude,lep.decimalLatitude.to_numpy(float),lep.decimalLongitude.to_numpy(float))
        q['nearest_sphingidae_km']=nearest_km(q.latitude,q.longitude,hawk.decimalLatitude.to_numpy(float),hawk.decimalLongitude.to_numpy(float))
        finite=q.nearest_any_km.notna()&q.nearest_sphingidae_km.notna()
        q['relative_hawkmoth_distance']=np.where(
            finite,
            np.log1p(q.nearest_sphingidae_km)-np.log1p(q.nearest_any_km),
            np.nan
        )
        dist_parts.append(q)
        any_cells=int(gi.cell1deg.nunique()) if len(gi) else 0
        hawk_cells=int(hawk.cell1deg.nunique()) if len(hawk) else 0
        lep_cells=int(lep.cell1deg.nunique()) if len(lep) else 0
        eligible=bool(len(hawk)>=1 and len(gi)>=5 and any_cells>=2 and finite.all())
        species_rows.append({
            'species':sp,
            'n_third_cohort_photo_coords':int(len(q)),
            'n_georeferenced_any_interactions':int(len(gi)),
            'n_georeferenced_bee_interactions':int(len(bee)),
            'n_georeferenced_lepidoptera_interactions':int(len(lep)),
            'n_georeferenced_sphingidae_interactions':int(len(hawk)),
            'n_any_cells_1deg':any_cells,
            'n_lepidoptera_cells_1deg':lep_cells,
            'n_sphingidae_cells_1deg':hawk_cells,
            'hawk_spatial_eligible':eligible,
            'hawk_multisite_eligible':bool(eligible and hawk_cells>=2),
            'lepidoptera_spatial_eligible':bool(len(lep)>=1 and len(gi)>=5 and any_cells>=2 and q.nearest_lepidoptera_km.notna().all())
        })
    dist=pd.concat(dist_parts,ignore_index=True)
    scov=pd.DataFrame(species_rows)
    dist.to_csv(out/'photo_pollinator_distance_preflight.csv.gz',index=False,compression='gzip')
    scov.to_csv(out/'stable_spatial_coverage_by_plant.csv',index=False)

    thresholds=[50,100,250,500]
    proximity={}
    for col in ['nearest_any_km','nearest_bee_km','nearest_lepidoptera_km','nearest_sphingidae_km']:
        finite=dist[col].notna()
        proximity[col]={'finite_photo_rows':int(finite.sum()),'finite_species':int(dist.loc[finite,'species'].nunique())}
        for th in thresholds:
            proximity[col][f'fraction_within_{th}km_among_finite']=float((dist.loc[finite,col]<=th).mean()) if finite.any() else None

    hawk=int(scov.hawk_spatial_eligible.sum())
    multisite=int(scov.hawk_multisite_eligible.sum())
    lep=int(scov.lepidoptera_spatial_eligible.sum())
    summary={
        'schema':'fcp_white_pollinator_stable_spatial_preflight_v1',
        'status':'complete',
        'role':'outcome_blind_spatial_preflight',
        'stable_source':{
            'zenodo_record':'18064921','version':'v3','published':'2025-12-26',
            'file':'interactions.tsv.gz','expected_md5':a.expected_md5,'observed_md5':actual
        },
        'plant_frame_species':499,
        'scan':scan,
        'species_with_any_stable_interaction':int(inter.plant_species.nunique()) if len(inter) else 0,
        'species_with_any_georeferenced_stable_interaction':int(geo.plant_species.nunique()) if len(geo) else 0,
        'hawk_spatial_gate':{'pass':bool(hawk>=30),'eligible_species':hawk,'required_species':30},
        'hawk_multisite_sensitivity_gate':{'pass':bool(multisite>=20),'eligible_species':multisite,'required_species':20},
        'lepidoptera_spatial_gate':{'pass':bool(lep>=50),'eligible_species':lep,'required_species':50},
        'proximity_descriptives':proximity,
        'planned_primary_covariate':'log1p(nearest_sphingidae_km) - log1p(nearest_any_pollinator_km)',
        'hard_nonclaims':['preflight does not read colour outcomes','interaction distances describe documentation geography rather than abundance or true absence','preflight cannot establish hawkmoth causation']
    }
    (out/'stable_spatial_preflight_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':
    main()
