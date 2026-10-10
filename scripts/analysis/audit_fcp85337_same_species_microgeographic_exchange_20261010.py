#!/usr/bin/env python3
"""Feasibility of within-species geographically matched flower-COLOUR nulls.

Original 85,337 source species×cell photo anchors, 39,075 existing classified
flower-photo states, 100,543 source geolocation/environment records. Before
building any source rain/solar/soil permutation claim, count same-species
photo sites genuinely close enough (50,100,250,500 km) for a geographically
matched colour-label exchange. Keep all source cohorts and missingness.
Eligibility is set only by fixed photo IDs and coordinates; whether labels
differ is reported separately and NEVER used to change group membership.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from compare_fcp_original85337_species_fixed_multiabiotic_20261010 import (
    match_original, populations,
)

DISTANCES_KM=(50.,100.,250.,500.)
EARTH_KM=6371.0088
N_CELLS_SOURCE=85337
CLASSIFIED_SOURCE=39075
N_ORIGINAL_SITES=100543
N_FOLDS=5
MIN_EXCHANGE_COLOUR_PHOTOS=100
MIN_EXCHANGE_COLOUR_FRACTION=.1
SCHEMA="fcp_original85337_within_species_nearby_photo_colour_exchangeability_v1"


def geometry(lat:np.ndarray,lon:np.ndarray)->np.ndarray:
    a=np.asarray(lat,float);b=np.asarray(lon,float)
    if len(a)!=len(b) or len(a)<1 or not (
        np.isfinite(a).all() and np.isfinite(b).all()
        and np.all(np.abs(a)<=90) and np.all(np.abs(b)<=180)
    ):raise ValueError("Original photo site coordinate missing or invalid")
    a=np.deg2rad(a);b=np.deg2rad(b)
    xyz=np.column_stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)))
    return EARTH_KM*np.arccos(np.clip(xyz@xyz.T,-1,1))


def source_heldout_mask(x:pd.DataFrame)->np.ndarray:
    from sklearn.model_selection import GroupKFold
    if not {"inat_taxon_id","cell_id"}.issubset(x):
        raise ValueError("Old photo species and 162-cell source fields missing")
    if x[["inat_taxon_id","cell_id"]].duplicated().any():
        raise ValueError("Repeated original species×region in same table")
    if x.cell_id.nunique()<N_FOLDS:
        raise ValueError("Original photo cell support cannot form 5 fold")
    eligibility=np.zeros(len(x),bool)
    for train,test in GroupKFold(n_splits=N_FOLDS).split(x,groups=x.cell_id):
        known=set(x.iloc[train].inat_taxon_id)
        eligibility[test]=x.iloc[test].inat_taxon_id.isin(known)
    return eligibility


def completed_source_microgroups(x:pd.DataFrame,distance_km:float)->list[np.ndarray]:
    if len(x)!=x.photo_id.nunique() or x[["inat_taxon_id","cell_id"]].duplicated().any():
        raise ValueError("One original photograph per species per original cell required")
    if distance_km<=0:raise ValueError("Photo source geographic distance cap positive")
    groups=[]
    for _,original in x.groupby("inat_taxon_id",sort=True):
        sorted_source=original.sort_values("photo_id",kind="stable")
        ids=sorted_source.index.to_numpy(int)
        d=geometry(sorted_source.latitude.to_numpy(float),sorted_source.longitude.to_numpy(float))
        clusters=[]
        for i in range(len(ids)):
            inserted=False
            for cluster in clusters:
                if all(d[i,k]<=distance_km+1e-8 for k in cluster):
                    cluster.append(i)
                    inserted=True
                    break
            if not inserted:clusters.append([i])
        for c in clusters:
            assert np.max(d[np.ix_(c,c)])<=distance_km+1e-7
            groups.append(ids[c])
    if len(groups)==0 or len(np.concatenate(groups))!=len(x):
        raise RuntimeError("Original source species-photo microgroups lost data")
    if len(set(np.concatenate(groups).tolist()))!=len(x):
        raise RuntimeError("Original source same photograph entered twice")
    return groups


def audit(x:pd.DataFrame)->dict:
    required={"inat_taxon_id","photo_id","cell_id","latitude","longitude",
              "morph","site_geo_status"}
    if not required.issubset(x):
        raise ValueError("Source labelled photo coordinate ledger incomplete")
    if x.index.tolist()!=list(range(len(x))):
        raise ValueError("Original source photo cohort must be zero-based index")
    if not x.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT").all():
        raise ValueError("Ungeolocated source photo entered geographically matched cohort")
    label=x.morph.to_numpy(str)
    allout={}
    for cap in DISTANCES_KM:
        clusters=completed_source_microgroups(x,cap)
        multiple=[ids for ids in clusters if len(ids)>=2]
        changes=[ids for ids in multiple if len(set(label[ids]))>=2]
        n_multi=int(sum(map(len,multiple)))
        n_changes=int(sum(map(len,changes)))
        taxa_multi={int(x.loc[ids[0],"inat_taxon_id"]) for ids in multiple}
        taxa_changes={int(x.loc[ids[0],"inat_taxon_id"]) for ids in changes}
        paircount=int(sum(len(ids)*(len(ids)-1)//2 for ids in multiple))
        heldout=x["eligible_heldout_original_species"].to_numpy(bool) if "eligible_heldout_original_species" in x else np.zeros(len(x),bool)
        test_nearby=int(sum(heldout[ids].sum() for ids in multiple))
        test_exchange=int(sum(heldout[ids].sum() for ids in changes))
        enough=(n_changes>=MIN_EXCHANGE_COLOUR_PHOTOS and
                n_changes/len(x)>=MIN_EXCHANGE_COLOUR_FRACTION)
        allout[str(int(cap))]={
            "max_site_pair_distance_km":cap,
            "n_source_species_photo_records":len(x),
            "n_unique_nominal_source_species":int(x.inat_taxon_id.nunique()),
            "n_same_species_photo_groups_including_singletons":len(clusters),
            "n_same_species_two_plus_source_photo_groups":len(multiple),
            "n_original_photos_in_two_plus_geographically_nearby_same_species_groups":n_multi,
            "n_same_species_taxa_with_two_plus_nearby_original_photo_sites":len(taxa_multi),
            "n_original_pairs_of_source_photos_within_complete_link_neighbourhoods":paircount,
            "n_groups_with_at_least_two_DIFFERENT_source_four_colours":len(changes),
            "n_source_photo_records_in_colour_exchangeable_nearby_groups":n_changes,
            "n_species_with_nearby_original_photo_colour_exchangeability":len(taxa_changes),
            "source_photo_fraction_in_colour_exchangeable_nearby_groups":float(n_changes/len(x)),
            "n_heldout_evaluable_original_photos_in_nearby_two_plus_groups":test_nearby,
            "n_heldout_evaluable_original_photos_in_exchangeable_groups":test_exchange,
            "group_membership_photo_colour_blind":True,
            "same_original_species_and_photo_cell_only":True,
            "complete_link_all_original_group_pair_distances_below_limit":True,
            "is_colour_permutation_null_identifiable":bool(enough),
            "readiness":"SOURCE_SAME_SPECIES_SPATIALLY_MATCHED_COLOUR_NULL_FEASIBLE" if enough
                        else "HOLD_INSUFFICIENT_NEARBY_INTRASPECIFIC_COLOUR_EXCHANGE",
        }
    return allout


def run(taxon_cell:pd.DataFrame,full_site:pd.DataFrame,*,strict:bool=True)->dict:
    old,ledger=match_original(taxon_cell,full_site,strict=strict)
    if strict and (len(old)!=N_CELLS_SOURCE or ledger["source_four_state_colour_classifiable"]!=CLASSIFIED_SOURCE):
        raise RuntimeError("Original full 85,337 photo-cell source or 39,075 label denominator drift")
    populations_,coverage=populations(old)
    cohorts={}
    for name,rows in populations_.items():
        x=rows.copy().reset_index(drop=True)
        x["eligible_heldout_original_species"]=source_heldout_mask(x)
        groups=audit(x)
        cohorts[name]={
            "n_original_classifiable_geolocated_environment_complete_photo_cells":len(x),
            "n_heldout_species_photo_records":int(x.eligible_heldout_original_species.sum()),
            "n_original_taxa_with_at_least_two_source_photo_cells":int(x.inat_taxon_id.value_counts().ge(2).sum()),
            "distance_audit":groups,
        }
    return {
        "schema":SCHEMA,"date_jst":"2026-10-10",
        "status":"SOURCE_ORIGINAL_SAME_SPECIES_SPATIALLY_MATCHED_PHOTO_COLOUR_EXCHANGEABILITY_PREFLIGHT_ONLY",
        "original_repeated_species_photo_ledger":ledger,
        "source_environment_coverage":coverage,
        "distance_caps_km_fixed_all_reported":list(DISTANCES_KM),
        "min_exchangeable_colour_photos_before_null":MIN_EXCHANGE_COLOUR_PHOTOS,
        "min_fraction_colour_exchangeable_before_null":MIN_EXCHANGE_COLOUR_FRACTION,
        "source_photo_colours_only_INSPECTED_FOR_POTENTIAL_EXCHANGE_NOT_RECLASSIFIED":True,
        "source_group_membership_no_outcome_selection":True,
        "no_additional_original_source_photo_pixels_or_environment_lookups":True,
        "same_species_source_photo_clusters_not_phylogenetic_species_turnover":True,
        "cohorts":cohorts,
        "hard_nonclaims":[
            "Multiple photos of one nominal species across geographic cells are not same individual, same population or proven colour genotypes",
            "An insufficient original nearby colour-exchange null is an identifiability bottleneck and not evidence of no phenotype–environment association",
            "A group contains at least two different original photographed coarse colours only after outcomes are examined for the feasibility count, never for inclusion",
            "Species and geography heldout predicted improvements in previous analyses are not geographically matched null p-values",
            "Source matched species × spatial cell is one photo per cell by design; nearby original photo availability may be severely limited",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-taxon-cell",required=True,type=Path)
    p.add_argument("--expanded-all-original-photo-sites",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    z=p.parse_args()
    report=run(pd.read_csv(z.original_taxon_cell,low_memory=False),
               pd.read_csv(z.expanded_all_original_photo_sites,low_memory=False))
    z.outdir.mkdir(parents=True,exist_ok=True)
    (z.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
