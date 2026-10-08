#!/usr/bin/env python3
"""Freeze one metadata-only original observation/photo identity for ALL FCP 42,111.

Uses *already archived* taxon/cell/observation/photo discovery indexes from
both independent iNaturalist global scans, not a new API call. A zero in the
later observer-capped capacity scan is NOT evidence no historical photo exists.
No image URL is guessed, no image bytes opened and no flower colour inferred.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SPECIES_SHA="5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
CAPACITY_SHA="12d94100d6343597aa2b87670555e801ce0681f4b69e3ad634fa6d80409008e8"
V1_SHA="b60a8b1b98fcf6745344d6acc394bb6ff905b9d5bded134031402892ff2efa2e"
V2_SHA="d7f08601ccc4d7dc9ff1943ea3e58ed903b9759277f7b536bc7423bc08189b3a"
N_WORLD_SPECIES=42111
SALT="FCP_ALL42111_SINGLE_HISTORIC_PHOTO_METADATA_20261008_V1"


def sha(path:Path)->str:
    dig=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):dig.update(chunk)
    return dig.hexdigest()


def read_inputs(species:Path, capacity:Path, v1:Path,v2:Path)->tuple[pd.DataFrame,pd.DataFrame]:
    for path,expected in ((species,SPECIES_SHA),(capacity,CAPACITY_SHA),(v1,V1_SHA),(v2,V2_SHA)):
        if sha(path)!=expected:
            raise RuntimeError(f"Historical source photo-ID metadata SHA256 drift: {path.name}")
    sp=pd.read_csv(species,usecols=["inat_taxon_id","species"],low_memory=False)
    cap=pd.read_csv(capacity,usecols=["inat_taxon_id","after_observer_cap"],low_memory=False)
    if len(sp)!=N_WORLD_SPECIES or len(cap)!=N_WORLD_SPECIES:
        raise RuntimeError("Not the complete historical 42111 source frame")
    if sp.inat_taxon_id.duplicated().any() or cap.inat_taxon_id.duplicated().any():
        raise RuntimeError("Original global photo discovery identities duplicated")
    joined=sp.merge(cap,on="inat_taxon_id",validate="one_to_one")
    if len(joined)!=N_WORLD_SPECIES:
        raise RuntimeError("Historical photo opportunity does not cover whole discovered frame")
    links=[]
    need=["inat_taxon_id","observation_id","photo_id","cell_id","species"]
    for tag,path in (("v1",v1),("v2",v2)):
        d=pd.read_csv(path,usecols=need,low_memory=False)
        d["source_scan"]=tag
        links.append(d)
    return joined,pd.concat(links,ignore_index=True)


def freeze_photo_identifiers(species:pd.DataFrame, evidence:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    must={"inat_taxon_id","observation_id","photo_id","cell_id","source_scan"}
    if not must.issubset(evidence):
        raise ValueError("Missing source-discovery photo-ID fields")
    if not {"inat_taxon_id","species","after_observer_cap"}.issubset(species):
        raise ValueError("Missing full species denominator or observer-capped capacity")
    all_taxa=set(species.inat_taxon_id.astype(int))
    e=evidence.copy().dropna(subset=["inat_taxon_id","observation_id","photo_id","cell_id"])
    for col in ("inat_taxon_id","observation_id","photo_id","cell_id"):
        e[col]=pd.to_numeric(e[col],errors="raise").astype("int64")
    if not set(e.inat_taxon_id).issubset(all_taxa):
        raise RuntimeError("Photo archive contains a taxon not in the fixed 42,111-species frame")
    if (e.photo_id<=0).any() or (e.observation_id<=0).any():
        raise RuntimeError("Photo or observation IDs are not valid")
    if (e.cell_id<0).any() or (e.cell_id>=162).any():
        raise RuntimeError("Photo discovery cell is outside fixed global equal-area range")

    # Never select image based on appearance, geographic suitability or
    # anything but immutable taxon/photo identities. Deduplicate rounds first.
    e=e.drop_duplicates(["inat_taxon_id","photo_id"],keep="first")
    e["photo_identity_hash"]=[
        hashlib.sha256(f"{SALT}|{int(t)}|{int(p)}".encode("utf-8")).hexdigest()
        for t,p in zip(e.inat_taxon_id,e.photo_id)
    ]
    e=e.sort_values(["inat_taxon_id","photo_identity_hash","photo_id"],
                    kind="stable")
    source_counts=e.groupby("inat_taxon_id",sort=False).size()
    if source_counts.size!=N_WORLD_SPECIES:
        missing=all_taxa-set(source_counts.index.astype(int))
        raise RuntimeError(
            f"Not every 42111 source-discovery species has an archived photo ID; n_missing={len(missing)}, first={sorted(missing)[:12]}"
        )
    # Rarest-source-first makes cross-taxon photo-ID uniqueness deterministic.
    ids_order=source_counts.sort_values(kind="stable").index.astype(int)
    grouped=e.groupby("inat_taxon_id",sort=False)
    used_obs:set[int]=set()
    used_photos:set[int]=set()
    chosen=[]
    unresolved=[]
    for taxon in ids_order:
        group=grouped.get_group(int(taxon))
        opt=group.loc[
            ~group.photo_id.isin(used_photos)&
            ~group.observation_id.isin(used_obs)
        ]
        if not len(opt):
            unresolved.append(int(taxon))
            continue
        selected=opt.iloc[0].to_dict()
        used_obs.add(int(selected["observation_id"]))
        used_photos.add(int(selected["photo_id"]))
        chosen.append(selected)
    if unresolved:
        raise RuntimeError(
            f"Cannot make all historical one-per-species photo identities globally disjoint by fixed rarity/hash ordering; n={len(unresolved)}, example={unresolved[:12]}. Preserve HOLD; do not select source photo by colour."
        )
    sampled=pd.DataFrame(chosen).drop(columns=["species"],errors="ignore")
    result=species.merge(sampled,on="inat_taxon_id",validate="one_to_one")
    if len(result)!=N_WORLD_SPECIES or result.photo_id.duplicated().any() or result.observation_id.duplicated().any():
        raise RuntimeError("Incomplete or duplicate one-species photo identity")
    # Compare source identity drift rather than guessing biological alias.
    original_name=evidence[["inat_taxon_id","species"]].drop_duplicates()
    names=original_name.groupby("inat_taxon_id").species.agg(lambda x:set(x.astype(str)))
    result["historical_source_name_exact_match"]=[
        bool(str(s) in names.get(int(t),set()))
        for t,s in zip(result.inat_taxon_id,result.species)
    ]
    result["new_photo_pixel_opened"]=False
    result["original_photo_URL_verified"]=False
    result["white_or_pigment_state_measured"]=False
    result["source_name_overwrites_current_species"]=False
    result=result.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)
    report={
        "n_all_world_species":len(result),
        "n_still_no_eligible_photo_by_later_capacity_scan":int(result.after_observer_cap.eq(0).sum()),
        "n_original_archived_photo_ids_even_when_later_capacity_zero":int(result.loc[result.after_observer_cap.eq(0),"photo_id"].notna().sum()),
        "n_distinct_archived_photo_IDs":int(result.photo_id.nunique()),
        "n_distinct_original_observation_IDs":int(result.observation_id.nunique()),
        "n_with_exact_discovery_species_name":int(result.historical_source_name_exact_match.sum()),
        "n_original_name_mismatch_flagged":int((~result.historical_source_name_exact_match).sum()),
        "n_v1_selected":int(result.source_scan.eq("v1").sum()),
        "n_v2_selected":int(result.source_scan.eq("v2").sum()),
        "n_historic_taxa_with_many_photo_candidates":int(source_counts.gt(1).sum()),
        "max_original_unique_photo_candidates_one_taxon":int(source_counts.max())
    }
    return result,report


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--species",required=True,type=Path)
    p.add_argument("--capacity",required=True,type=Path)
    p.add_argument("--links-v1",required=True,type=Path)
    p.add_argument("--links-v2",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    sp,frames=read_inputs(args.species,args.capacity,args.links_v1,args.links_v2)
    sample,stats=freeze_photo_identifiers(sp,frames)
    args.outdir.mkdir(parents=True,exist_ok=True)
    dest=args.outdir/"one_source_photo_id_per_each_42111_species.csv.gz"
    sample.to_csv(dest,index=False,lineterminator="\n",
                  compression={"method":"gzip","compresslevel":9,"mtime":0})
    record={
        "schema":"fcp_all42111_preexisting_archived_single_photo_id_manifest_v1",
        "status":"ORIGINAL_HISTORICAL_SOURCE_ID_42111_COMPLETE_NO_COLOUR_PIXELS",
        "date_jst":"2026-10-08",
        "selection_method":"one source photo ID per taxon; historical archived observation indexes, source-photo identity SHA, rarest-source-first deterministic collision handling",
        "source_sha256":{
            "species":SPECIES_SHA,"capacity":CAPACITY_SHA,
            "discovery_v1":V1_SHA,"discovery_v2":V2_SHA
        },
        "species_counts":stats,
        "archived_source_photo_id_manifest_sha256":sha(dest),
        "fresh_independent_new_photo_source":False,
        "new_image_pixels_opened":False,
        "photo_urls_retrieved":False,
        "colour_labels_or_pigment_processed":False,
        "historical_one-photo_not_fitness_or_ITV":"One candidate photo for a species does not demonstrate genetic colour polymorphism, its absence, or native site pigment frequency",
        "conclusion":"Complete original photo identity availability across all 42111 historically discovered species; previous 789 zeros denote later capped/exclusion stage, not total original photo absence",
        "confirmatory_decisions_changed":False
    }
    (args.outdir/"result.json").write_text(json.dumps(record,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(record,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
