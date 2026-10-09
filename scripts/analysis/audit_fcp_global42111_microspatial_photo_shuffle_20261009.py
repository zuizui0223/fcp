#!/usr/bin/env python3
"""Microgeography-constrained colour-label permutation control for FCP 42111.

Restrict the already-frozen 250/500 km source congeneric photo cohorts to
null label exchanges among geographically NEAR original photographed sites,
rather than shuffle anywhere within a potentially 500-km genus×cell group.

No new photo pixels, taxonomic IDs, environmental rasters, morph labels,
genotypes, or prospective taxa. Retrospective conditional null only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort, outofsource_photo_moisture_gain,
    EXPECTED_ORIGINAL_N, EXPECTED_OBS_GAIN, SEED,
)
from audit_fcp_global42111_local_congeners_20261009 import (
    photo_populations, fold_plan, CLASSES,
)
from audit_fcp_global42111_local_congener_diameters_20261009 import EARTH_KM

GROUP_CAPS_KM=(250,500)
MICRO_DIAMETERS_KM=(50,100)
N_NULL=199
MIN_INFORMATIVE_SPECIES=100
MIN_INFORMATIVE_FRACTION=.1
SCHEMA="fcp_global42111_local_congener_microspatial_shuffle_null_v1"


def great_circle_matrix_km(lat: np.ndarray,lon: np.ndarray)->np.ndarray:
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    if lat.ndim!=1 or lon.ndim!=1 or len(lat)!=len(lon):
        raise ValueError("Invalid original photo coordinates")
    if not(np.isfinite(lat).all() and np.isfinite(lon).all() and
           np.all(np.abs(lat)<=90) and np.all(np.abs(lon)<=180)):
        raise ValueError("Original photograph coordinate missing or invalid")
    phi=np.deg2rad(lat)
    theta=np.deg2rad(lon)
    xyz=np.column_stack((np.cos(phi)*np.cos(theta),
                         np.cos(phi)*np.sin(theta),np.sin(phi)))
    cos=np.clip(xyz@xyz.T,-1,1)
    return np.arccos(cos)*EARTH_KM


def microgroups(source:pd.DataFrame, max_diameter_km:float)->list[np.ndarray]:
    """Greedy deterministic complete-link partition; all within-subgroup
    pairwise ORIGINAL photo distances <= fixed diameter.
    Order uses source taxon ID, never source flower colour or climatic value.
    """
    expected={"genus_cell_id","inat_taxon_id","latitude","longitude"}
    if not expected.issubset(source):
        raise ValueError("Missing original species genus-cell/source point")
    if source.inat_taxon_id.duplicated().any():
        raise ValueError("Repeated original photo species ID")
    if max_diameter_km<=0:
        raise ValueError("Microgeographic diameter must be positive")
    buckets=[]
    for _, group in source.groupby("genus_cell_id",sort=True):
        ordered=group.sort_values("inat_taxon_id",kind="stable")
        ids=ordered.index.to_numpy(dtype=int)
        dist=great_circle_matrix_km(
            ordered.latitude.to_numpy(float),ordered.longitude.to_numpy(float))
        clusters=[]
        for k in range(len(ids)):
            matched=False
            for existing in clusters:
                if all(dist[k,j]<=max_diameter_km+1e-8 for j in existing):
                    existing.append(k)
                    matched=True
                    break
            if not matched:
                clusters.append([k])
        for c in clusters:
            local=ids[c]
            if len(c)>1 and np.max(dist[np.ix_(c,c)])>max_diameter_km+1e-7:
                raise RuntimeError("Microgeographic connected-chain failed complete-link bound")
            buckets.append(local)
    flat=np.concatenate(buckets)
    if len(flat)!=len(source) or len(set(flat.tolist()))!=len(source):
        raise RuntimeError("Permutation neighborhood lost original photo species")
    return buckets


def assess(source:pd.DataFrame,diameter_km:int, nperm:int,seed:int)->dict:
    if source.index.tolist()!=list(range(len(source))):
        raise ValueError("Original source cohort must have stable zero-based indices")
    group_ids=microgroups(source,float(diameter_km))
    base=source.morph.to_numpy(str)
    if not np.isin(base,CLASSES).all():
        raise ValueError("Original source photo-coarse four-state labels missing")
    informative=[idx for idx in group_ids
                 if len(idx)>=2 and len(set(base[idx]))>=2]
    moved=np.unique(np.concatenate(informative)) if informative else np.array([],int)
    ninf=len(moved)
    gate=ninf>=MIN_INFORMATIVE_SPECIES and ninf/len(source)>=MIN_INFORMATIVE_FRACTION
    stats={
        "microgeographic_photo_diameter_km":diameter_km,
        "n_frozen_original_photo_species":len(source),
        "n_source_genus_cell_groups":int(source.genus_cell_id.nunique()),
        "n_microgeographic_complete_link_groups":len(group_ids),
        "n_groups_with_multiple_original_species":sum(len(q)>=2 for q in group_ids),
        "n_groups_with_multiple_distinct_photo_colours":len(informative),
        "n_source_species_in_informatively_shufflable_microgroups":ninf,
        "fraction_source_species_in_informatively_shufflable_microgroups":ninf/len(source),
        "minimum_effective_photo_exchange_count":MIN_INFORMATIVE_SPECIES,
        "minimum_effective_photo_exchange_fraction":MIN_INFORMATIVE_FRACTION,
        "all_source_photos_covered_by_exactly_one_complete_link_group":True,
        "group_membership_uses_only_original_species_ID_genus_and_public_photo_coordinates":True,
        "original_genus_cell_colour_counts_always_preserved":True,
        "microgroup_colour_counts_preserved_on_each_null":True,
        "source_test_fold_assignments_unchanged":True,
        "status":"READY_FOR_CONDITIONAL_SPATIAL_NULL" if gate else "HOLD_INSUFFICIENT_MICROSPATIAL_EXCHANGEABILITY",
    }
    if not gate:
        stats["p_value"]=None
        stats["n_refitted_null_models"]=0
        return stats
    folds=fold_plan(source)
    obs=outofsource_photo_moisture_gain(source,folds)
    rng=np.random.default_rng(seed)
    null=[]
    for _ in range(nperm):
        perm=base.copy()
        for ids in informative:
            perm[ids]=rng.permutation(base[ids])
        current=source.copy()
        current["morph"]=perm
        null.append(outofsource_photo_moisture_gain(current,folds))
    v=np.asarray(null,float)
    stats.update({
        "n_refitted_null_models":nperm,
        "observed_moisture_brier_gain":float(obs),
        "null_mean":float(v.mean()),
        "null_std":float(v.std(ddof=1)) if nperm>1 else None,
        "null_q025_q500_q975":[float(t) for t in np.quantile(v,[.025,.5,.975])],
        "null_count_ge_observed":int(np.sum(v>=obs-1e-12)),
        "p_value":float((1+np.sum(v>=obs-1e-12))/(nperm+1)),
    })
    return stats


def run(source:pd.DataFrame,*,strict:bool=True,nperm:int=N_NULL)->dict:
    pops,coverage=photo_populations(source,strict=strict)
    fixed=pops["CLIMATE_ALL"]
    cohorts={}
    for cap in GROUP_CAPS_KM:
        cohort=fixed_source_cohort(fixed,cap)
        if strict and len(cohort)!=EXPECTED_ORIGINAL_N[str(cap)]:
            raise RuntimeError("Photo source distance and taxon eligibility drift")
        fold=fold_plan(cohort)
        observed=outofsource_photo_moisture_gain(cohort,fold)
        if strict and abs(observed-EXPECTED_OBS_GAIN[str(cap)])>1e-8:
            raise RuntimeError("Original source observed moisture Brier check failed")
        diagnostics={}
        for diameter in MICRO_DIAMETERS_KM:
            diagnostics[str(diameter)]=assess(
                cohort,diameter,nperm,SEED+cap*10+diameter)
            if diagnostics[str(diameter)]["status"]=="READY_FOR_CONDITIONAL_SPATIAL_NULL":
                if abs(diagnostics[str(diameter)]["observed_moisture_brier_gain"]-observed)>1e-12:
                    raise RuntimeError("Observed photo source model changed")
        cohorts[str(cap)]={
            "n_original_source_photos":len(cohort),
            "n_original_genus_cell_groups":int(cohort.genus_cell_id.nunique()),
            "previously_frozen_observed_moisture_gain":float(observed),
            "microgeographic_diagnostics":diagnostics,
        }
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-09",
        "status":"SOURCE_ORIGINAL_PHOTOGRAPH_MICROGEOGRAPHY_EXCHANGEABILITY_AUDIT",
        "original_global_source_taxa":len(source),
        "original_colour_classifiable_taxa":coverage["source_classified_photo_taxa"],
        "original_climate_complete_photo_species":coverage["climate_eligible_classified_source_taxa"],
        "original_distance_caps_km":list(GROUP_CAPS_KM),
        "microgeographic_complete_link_diameter_caps_km":list(MICRO_DIAMETERS_KM),
        "null_permutations_each_estimable_comparison":nperm,
        "original_photos_and_climate_sources_unchanged":True,
        "all_groups_selected_without_photo_colour_outcomes":True,
        "cohorts":cohorts,
        "nonclaims":[
            "Matched-null label shuffling conditional on photographer-provided locations cannot remove taxon phylogenetic subgenus effects",
            "Limited exchangeable near-neighbor groups can make a microgeographic permutation test uninformative; HOLD is not a null phenotype result",
            "The original photo color may be atypical for a species or taxonomically inaccurate",
            "Even a rare photo-colour shuffle outcome is not causal evidence for climatic selection or fitness",
            "The original 250km and 500km source cohorts overlap and do not constitute independent replication",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    d=pd.read_csv(args.original_breadth_abiotic,low_memory=False)
    result=run(d)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
