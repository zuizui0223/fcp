#!/usr/bin/env python3
"""Conservative original FCP five-cohort H9 photo-ID reuse correction.

The first three 149900 actual pixel measurements are immutable. Both 256-worker
ROI jobs are already running on their frozen 49800 and 248100 photo-ID sets.

An independent real five-cohort source check found exactly three reused
discovery photo IDs and one reused observation ID in the final U100 H9 source.
To retain original per-species EXACTLY 100 distinct sampled photos AND equal
species weighting, never select another image, silently drop one photo to 99,
or change an old photo. Exclude the ENTIRE affected new taxon from statistical
inference, report the actual lost species/100-photo rows separately, and keep
all original source/ROI artifacts immutable. No ecological test is done here.
"""
from __future__ import annotations
import json
from typing import Any

import pandas as pd

PINNED_OLD_DISCOVERY_PHOTO_COLLISIONS = frozenset({
    45948761, 311814882, 517553374,
})
PINNED_OLD_DISCOVERY_OBSERVATION_COLLISIONS = frozenset({
    287821814,
})
EXPECTED_NEW_H9_COMPLETE_TAXA = 2481
EXPECTED_NEW_H9_SOURCE_PHOTOS = 248100


def inspect_reused_new_u100_source(
    new_source: pd.DataFrame,
    old_by_cohort: dict[str,pd.DataFrame],
    *,
    strict_original: bool = True,
) -> tuple[pd.DataFrame,dict[str,Any]]:
    """Return genuine new photo rows retained after dropping entire collided taxa.

    Every previously used old observation/photo ID is checked. The precise real
    discovery-only collision IDs are pinned so new unexpected source quality
    failures are NOT silently reconciled by automatic exclusion.
    """
    required={"photo_id","observation_id","inat_taxon_id","species"}
    if required - set(new_source):
        raise RuntimeError("original new H9 source identity keys missing")
    if set(old_by_cohort)!={"discovery","validation","prospective","wave4"}:
        raise RuntimeError("four original previous cohorts required for identity reuse")
    numrows=len(new_source)
    numeric_new=new_source.copy()
    for col in ("photo_id","observation_id","inat_taxon_id"):
        numeric_new[col]=pd.to_numeric(numeric_new[col],errors="raise").astype("int64")
        if numeric_new[col].le(0).any():
            raise RuntimeError("original photographic source identifier not positive")
    if (numeric_new.photo_id.nunique()!=numrows
        or numeric_new.observation_id.nunique()!=numrows):
        raise RuntimeError("duplicate photographic or observation IDs already in H9 new source")
    species_count=int(numeric_new.inat_taxon_id.nunique())
    if numeric_new.groupby("inat_taxon_id").size().ne(100).any():
        raise RuntimeError("original H9 new source must have exactly 100 photos per taxon")
    if numeric_new[["inat_taxon_id","species"]].drop_duplicates().shape[0]!=species_count:
        raise RuntimeError("original H9 source taxon labels not stable")
    if strict_original and (numrows!=EXPECTED_NEW_H9_SOURCE_PHOTOS
                            or species_count!=EXPECTED_NEW_H9_COMPLETE_TAXA):
        raise RuntimeError("not the SHA-authenticated original 248100 new H9 source")
    photo_reuse:set[int]=set()
    observation_reuse:set[int]=set()
    per_old={}
    for key in ("discovery","validation","prospective","wave4"):
        previous=old_by_cohort[key]
        if "photo_id" not in previous or "observation_id" not in previous:
            raise RuntimeError("cannot verify original source photograph + observation reuse")
        photo_old=set(pd.to_numeric(previous["photo_id"],errors="raise").astype("int64").tolist())
        observation_old=set(pd.to_numeric(previous["observation_id"],errors="raise").astype("int64").tolist())
        repeated_photos=photo_old.intersection(numeric_new.photo_id)
        repeated_observations=observation_old.intersection(numeric_new.observation_id)
        per_old[key]={
            "original_reused_photo_IDs":sorted(int(v) for v in repeated_photos),
            "original_reused_observation_IDs":sorted(int(v) for v in repeated_observations),
        }
        photo_reuse.update(repeated_photos)
        observation_reuse.update(repeated_observations)
    if strict_original:
        if (photo_reuse!=PINNED_OLD_DISCOVERY_PHOTO_COLLISIONS
            or observation_reuse!=PINNED_OLD_DISCOVERY_OBSERVATION_COLLISIONS
            or any(per_old[k]["original_reused_photo_IDs"] or per_old[k]["original_reused_observation_IDs"]
                   for k in ("validation","prospective","wave4"))):
            raise RuntimeError("the full real five-cohort reuse identity fingerprints drifted")
    mask=(numeric_new.photo_id.isin(photo_reuse)
          | numeric_new.observation_id.isin(observation_reuse))
    affected=numeric_new.loc[mask,["inat_taxon_id","species"]].drop_duplicates()
    excluded_taxon_ids=set(affected.inat_taxon_id.astype(int))
    # Entire original 100-photo species are conservatively excluded. No images
    # or outcomes from these taxa influence the equally-weighted statistics.
    retain=~numeric_new.inat_taxon_id.isin(excluded_taxon_ids)
    analysis_photo=numeric_new.loc[retain].copy()
    if (len(analysis_photo)!=numrows-100*len(excluded_taxon_ids)
        or analysis_photo.groupby("inat_taxon_id").size().ne(100).any()
        or analysis_photo.photo_id.isin(photo_reuse).any()
        or analysis_photo.observation_id.isin(observation_reuse).any()):
        raise RuntimeError("failed to conservatively preserve original 100-photo/taxon method")
    audit={
        "schema":"fcp_4730_original_H9_discovery_photo_reuse_complete_taxon_exclusion_v1",
        "original_new_H9_source_photo_rows_before_frozen_photo_ID_reuse":numrows,
        "original_new_H9_complete_taxa_before_frozen_photo_ID_reuse":species_count,
        "original_reused_photo_ID_count":len(photo_reuse),
        "original_reused_observation_ID_count":len(observation_reuse),
        "original_reused_photo_IDs":sorted(photo_reuse),
        "original_reused_observation_IDs":sorted(observation_reuse),
        "old_original_cohort_reuse_breakdown":per_old,
        "affected_new_H9_complete_species":len(excluded_taxon_ids),
        "affected_new_H9_taxon_IDs":sorted(excluded_taxon_ids),
        "affected_new_H9_species_names":sorted(affected.species.astype(str).unique().tolist()),
        "original_source_photo_rows_removed_whole_taxon":100*len(excluded_taxon_ids),
        "effective_new_U100_complete100_taxa":species_count-len(excluded_taxon_ids),
        "effective_new_U100_source_photo_rows":len(analysis_photo),
        "old_149900_unmodified":True,
        "original_wave4_49800_unmodified":True,
        "original_248100_pixel_measurement_jobs_unmodified":True,
        "never_use_99_photo_taxon":True,
        "no_missing_photo_replacement_or_colour_based_selection":True,
        "original_50km_H1_D_W_estimators_unchanged":True,
        "original_U100_4730_taxa_all_retain_auditable_status":True,
    }
    return analysis_photo,audit


def brief(report:dict)->str:
    return json.dumps({
        "excluded_taxa":report["affected_new_H9_taxon_IDs"],
        "excluded_n":report["affected_new_H9_complete_species"],
        "retained_taxa":report["effective_new_U100_complete100_taxa"],
        "retained_photo_rows":report["effective_new_U100_source_photo_rows"],
        "reused_photos":report["original_reused_photo_IDs"],
        "reused_observations":report["original_reused_observation_IDs"],
    },sort_keys=True)
