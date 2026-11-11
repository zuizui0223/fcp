"""Species-wise environment structure on the identical measured FCP photo IDs.

This is a descriptive, species-equal *within-species* layer. Supply complete
species partitions, not partial sampling shards. Climate/soil come from the
photo's own geocoordinate via fcp_direct_environment; no mois-photo proxy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .fcp_direct_environment import MOIS_36


def _bool(s:pd.Series)->pd.Series:
    return s.fillna('').astype(str).str.strip().str.lower().isin(('true','yes','1','y'))


def paired_profiles(photos:pd.DataFrame,environment:pd.DataFrame,
                    *, grid_degrees:float=0.1, min_colour:int=5,
                    min_sites:int=3) -> tuple[pd.DataFrame,pd.DataFrame]:
    if not (0 < grid_degrees <= 1):
        raise ValueError('grid must be a predefined valid geographic resolution')
    pp={'photo_id','inat_taxon_id','species','morph','global_classifiable','observer_id'}
    ee={'photo_id','inat_taxon_id','latitude','longitude','location_source'}
    if pp-set(photos) or ee-set(environment):
        raise ValueError('Missing measured flower or environmental photo-identity fields')
    p=photos.copy();e=environment.copy()
    if e.location_source.ne('direct_source_photo_coordinate').any():
        raise ValueError('Nearest mois photo is not a valid site-specific environment')
    p['photo_id']=p.photo_id.astype(str)
    e['photo_id']=e.photo_id.astype(str)
    if (p.photo_id.duplicated().any() or e.photo_id.duplicated().any() or
        p[['inat_taxon_id','photo_id']].duplicated().any()):
        raise ValueError('Photo identities must be unique before site join')
    if e.location_source.ne('direct_source_photo_coordinate').any():
        raise ValueError('Nearest mois photo is not a valid target-species environment')
    variables=sorted(set(MOIS_36)&set(e.columns))
    if not variables:
        raise ValueError('No matching mois-named source-raster features')
    joined=p.merge(e[['photo_id','inat_taxon_id','latitude','longitude',*variables]],
                   on=['photo_id','inat_taxon_id'],how='left',validate='one_to_one',indicator=True)
    if not joined['_merge'].eq('both').all():
        raise ValueError('Environmental table does not cover every supplied photo ID')
    joined['classifiable']=_bool(joined.global_classifiable)&joined.morph.isin(('white','yellow_orange','red_pink','blue_purple'))
    lat=pd.to_numeric(joined.latitude,errors='coerce')
    lon=pd.to_numeric(joined.longitude,errors='coerce')
    if not (lat.between(-90,90)&lon.between(-180,180)).all():
        raise ValueError('All photo coordinates must be complete for the submitted complete-species partition')
    joined['cell']=np.floor(lat.to_numpy(float)/grid_degrees).astype('int64').astype(str)+'|'+np.floor(lon.to_numpy(float)/grid_degrees).astype('int64').astype(str)
    profiles=[]
    sorting=[]
    for (taxon,species),group in joined.groupby(['inat_taxon_id','species'],sort=True):
        if group['photo_id'].nunique()!=len(group):
            raise ValueError('Duplicate photo ID within species')
        for variable in variables:
            take=group.copy()
            take[variable]=pd.to_numeric(take[variable],errors='coerce')
            take=take.loc[np.isfinite(take[variable])].copy()
            site_means=take.groupby('cell')[variable].mean()
            nsites=len(site_means)
            spread=float(site_means.std(ddof=0)) if nsites>=min_sites else np.nan
            profiles.append({'inat_taxon_id':int(taxon),'species':str(species),
                             'variable':variable, 'n_environment_photos':len(take),
                             'n_distinct_spatial_cells':nsites,
                             'environment_spread_site_equal':spread,
                             'location_source':'direct_source_photo_coordinate',
                             'grid_degrees':grid_degrees})
            eligible=take.loc[take.classifiable].copy()
            white=eligible[eligible.morph.eq('white')]
            coloured=eligible[eligible.morph.ne('white')]
            if min(len(white),len(coloured)) < min_colour:
                continue
            scale=float(eligible[variable].std(ddof=0))
            if not (np.isfinite(scale) and scale>0):
                continue
            # A signed within-species environmental contrast is SECONDARY to
            # W-vs-environment heterogeneity. It is never pigment fitness.
            effect=float((white[variable].mean()-coloured[variable].mean())/scale)
            observer_pairs=[]
            for _,obs in eligible.groupby('observer_id'):
                a=obs.loc[obs.morph.eq('white'),variable]
                b=obs.loc[obs.morph.ne('white'),variable]
                if len(a) and len(b): observer_pairs.append((a.mean()-b.mean())/scale)
            common_cells=[]
            for _,site in eligible.groupby('cell'):
                a=site.loc[site.morph.eq('white'),variable]
                b=site.loc[site.morph.ne('white'),variable]
                if len(a) and len(b): common_cells.append((a.mean()-b.mean())/scale)
            sorting.append({'inat_taxon_id':int(taxon),'species':str(species),
                            'variable':variable,'white_minus_nonwhite_sd':effect,
                            'n_white':len(white),'n_nonwhite':len(coloured),
                            'n_observer_pairs':len(observer_pairs),
                            'same_observer_effect_mean':float(np.mean(observer_pairs)) if len(observer_pairs)>=3 else np.nan,
                            'n_mixed_cells':len(common_cells),
                            'same_cell_effect_mean':float(np.mean(common_cells)) if len(common_cells)>=3 else np.nan,
                            'location_source':'direct_source_photo_coordinate'})
    return pd.DataFrame(profiles),pd.DataFrame(sorting)


def main(argv:list[str]|None=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--measured-complete-species-partition',type=Path,required=True)
    parser.add_argument('--environment',type=Path,required=True)
    parser.add_argument('--grid-degrees',type=float,default=.1)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args(argv)
    pp=pd.read_csv(args.measured_complete_species_partition,dtype={'photo_id':str},low_memory=False)
    ee=pd.read_csv(args.environment,dtype={'photo_id':str},low_memory=False)
    spread,sorting=paired_profiles(pp,ee,grid_degrees=args.grid_degrees)
    args.outdir.mkdir(parents=True,exist_ok=True)
    spread.to_csv(args.outdir/'species_environment_heterogeneity.csv',index=False)
    sorting.to_csv(args.outdir/'species_white_nonwhite_environment_sorting.csv',index=False)
    receipt={'schema':'fcp_species_environment_structure_v1',
             'n_input_photo_ids':len(pp),'n_species_with_environment':int(spread.inat_taxon_id.nunique()),
             'n_variables':int(spread.variable.nunique()),'within_species_environment_spread_is_site_equal':True,
             'within_species_white_environment_sorting_is_secondary':True,
             'identifies_common_environmental_causation':False,
             'requires_observer_exposure_and_geographic_negative_controls':True,
             'requires_constructed_H2_null_before_inference':True,
             'all_measured_photos_included_no_depth_cap':True}
    (args.outdir/'environment_structure_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
