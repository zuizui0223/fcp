#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

EARTH_RADIUS_KM=6371.0

def xyz(lat,lon):
    latr=np.radians(np.asarray(lat,float)); lonr=np.radians(np.asarray(lon,float))
    return np.column_stack([np.cos(latr)*np.cos(lonr),np.cos(latr)*np.sin(lonr),np.sin(latr)])

def trim_to_radial_core(g,q):
    X=xyz(g.decimalLatitude,g.decimalLongitude)
    c=X.mean(0); c=c/np.linalg.norm(c)
    ang=np.arccos(np.clip(X@c,-1,1))
    cutoff=float(np.quantile(ang,q))
    return g.loc[ang<=cutoff].copy(), cutoff

def add_tc_grid(g,res):
    g=g.copy()
    g['tc_lat_idx']=np.floor((90-g.decimalLatitude.astype(float))/res).astype(int).clip(0,int(round(180/res))-1)
    g['tc_lon_idx']=np.floor((g.decimalLongitude.astype(float)+180)/res).astype(int).clip(0,int(round(360/res))-1)
    return g.sort_values(['tc_lat_idx','tc_lon_idx','decimalLatitude','decimalLongitude']).drop_duplicates(['tc_lat_idx','tc_lon_idx']).reset_index(drop=True)

def maximin(g,k):
    X=xyz(g.decimalLatitude,g.decimalLongitude)
    c=X.mean(0); c=c/np.linalg.norm(c)
    first=int(np.argmin(np.sum((X-c)**2,axis=1)))
    chosen=[first]; nearest=np.sum((X-X[first])**2,axis=1); nearest[first]=-1
    while len(chosen)<k:
        j=int(np.argmax(nearest)); chosen.append(j)
        nearest=np.minimum(nearest,np.sum((X-X[j])**2,axis=1)); nearest[chosen]=-1
    out=g.iloc[chosen].copy(); out['selection_rank']=np.arange(1,k+1)
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--occurrences',required=True); p.add_argument('--core-species',required=True); p.add_argument('--analysis',required=True)
    p.add_argument('--out',required=True); p.add_argument('--qc-out',required=True); p.add_argument('--species-qc-out',required=True)
    p.add_argument('--points-per-species',type=int,default=20); p.add_argument('--radial-quantile',type=float,default=.95)
    p.add_argument('--terraclimate-resolution',type=float,default=1/24)
    a=p.parse_args()
    occ=pd.read_csv(a.occurrences); core=pd.read_csv(a.core_species); analysis=pd.read_csv(a.analysis)
    names=sorted(set(core.canonical_name).intersection(set(analysis.canonical_name)))
    meta=core[['canonical_name','family','organization_state','C_local_coexistence_documented','S_spatial_segregation_documented']].drop_duplicates('canonical_name')
    pieces=[]; qrows=[]
    for name in names:
        g=occ.loc[occ.canonical_name.eq(name)].copy()
        if g.empty: raise SystemExit(f'No occurrences for {name}')
        trimmed,cutoff=trim_to_radial_core(g,a.radial_quantile)
        cells=add_tc_grid(trimmed,a.terraclimate_resolution)
        if len(cells)<a.points_per_species: raise SystemExit(f'{name}: only {len(cells)} unique TerraClimate cells after trim')
        pieces.append(maximin(cells,a.points_per_species))
        qrows.append({'canonical_name':name,'raw_occurrence_rows':len(g),'rows_after_radial_trim':len(trimmed),'unique_tc_cells_after_trim':len(cells),'radial_cutoff_km':cutoff*EARTH_RADIUS_KM})
    out=pd.concat(pieces,ignore_index=True)
    out['terraclimate_lat']=90-(out.tc_lat_idx.astype(float)+.5)*a.terraclimate_resolution
    out['terraclimate_lon']=-180+(out.tc_lon_idx.astype(float)+.5)*a.terraclimate_resolution
    out=out.drop(columns=['organization_state','C_local_coexistence_documented','S_spatial_segregation_documented'],errors='ignore').merge(meta,on=['canonical_name','family'],how='left',validate='many_to_one')
    cols=['canonical_name','family','organization_state','C_local_coexistence_documented','S_spatial_segregation_documented','selection_rank','decimalLatitude','decimalLongitude','terraclimate_lat','terraclimate_lon','tc_lat_idx','tc_lon_idx','gbif_key','year','basisOfRecord','datasetKey']
    out=out[cols].sort_values(['canonical_name','selection_rank'])
    op=Path(a.out); op.parent.mkdir(parents=True,exist_ok=True); out.to_csv(op,index=False)
    pd.DataFrame(qrows).to_csv(a.species_qc_out,index=False)
    qc={'protocol':'fixed-spatial-points-v1','core_species':len(names),'points_per_species':a.points_per_species,'rows':len(out),'radial_quantile':a.radial_quantile,'terraclimate_grid_resolution_degrees':a.terraclimate_resolution,'selection':'radial core -> TerraClimate-cell deduplication -> deterministic spherical maximin; start at cell nearest retained spherical centroid','source_sampling_guard':'inherits capped QC-filtered GBIF search sample; not a complete range census','min_unique_terraclimate_cells_after_trim':min(r['unique_tc_cells_after_trim'] for r in qrows),'states':out.drop_duplicates('canonical_name').organization_state.value_counts().to_dict()}
    Path(a.qc_out).write_text(json.dumps(qc,indent=2)+'\n'); print(json.dumps(qc,indent=2))
if __name__=='__main__': main()
