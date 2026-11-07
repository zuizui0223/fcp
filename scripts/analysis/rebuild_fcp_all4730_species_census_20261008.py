#!/usr/bin/env python3
"""Rebuild the EXACT 4730-species historical U100 universe and allocated strata.

Metadata only. No photograph pixels, flower colour, previous model outcomes,
species replacement, environmental or geographic selection is involved.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

SPECIES_HASH="5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
CAPACITY_HASH="12d94100d6343597aa2b87670555e801ce0681f4b69e3ad634fa6d80409008e8"
P100_HASH="1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c"
P500_HASH="f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"
THIRD_HASH="16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59"
MAIN_HASH="e1301095533dfa755c1588cc550522bf31aabaef632083e265a138f0e3af6d2a"
REPLICATION_HASH="dbfaf954c1b82c569e9495006105639ad1d27b356a20b0ee3cfb5b4e240c00e8"
GROUPS={
    "historical_discovery_reserve":1000,
    "historical_p500":500,
    "historical_third":500,
    "new_primary":2000,
    "new_replication":730,
}
OUTCOME_FIREWALL={
    "all_photo_pixels_unopened_by_this_script":True,
    "new_colour_morph_measurement_performed":False,
    "previous_colour_outcome_used_in_selection":False,
    "environment_or_latitude_used_in_selection":False,
    "historical_p500_and_third_failures_reclassified_as_confirmation":False
}


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1<<20),b""):h.update(buf)
    return h.hexdigest()


def read_ids(path:Path, expected_hash:str, n:int,sep=",")->pd.DataFrame:
    if sha(path)!=expected_hash:raise RuntimeError(f"Source identity hash mismatch {path.name}")
    x=pd.read_csv(path,sep=sep,low_memory=False)
    if not {"inat_taxon_id","species"}.issubset(x.columns):raise ValueError(f"Identity fields absent {path.name}")
    if len(x)!=n:raise ValueError(f"Source allocation count mismatch {path.name} {len(x)}")
    x["inat_taxon_id"]=pd.to_numeric(x.inat_taxon_id,errors="raise").astype("int64")
    x["species"]=x.species.astype(str)
    if x.inat_taxon_id.duplicated().any() or x.species.duplicated().any():
        raise ValueError(f"Duplicated taxon or scientific binomial {path.name}")
    return x


def group_allocation(species:pd.DataFrame, capacity:pd.DataFrame,p100:pd.DataFrame,p500:pd.DataFrame,
                     third:pd.DataFrame,main:pd.DataFrame,replication:pd.DataFrame)->pd.DataFrame:
    if len(species)!=42111 or len(capacity)!=42111:
        raise RuntimeError("Frozen global species opportunity census should have 42111 entries")
    if not {"inat_taxon_id","species"}.issubset(species) or not {
        "inat_taxon_id","after_observer_cap","request_error"}.issubset(capacity):
        raise RuntimeError("Required original census fields absent")
    if species.inat_taxon_id.duplicated().any() or capacity.inat_taxon_id.duplicated().any():
        raise RuntimeError("Frozen source census repeated taxon ID")
    # Capacity at species rank is selected entirely without flower-colour or
    # population/latitudinal filtering, preserving original opportunity.
    both=species[["inat_taxon_id","species"]].merge(
        capacity[["inat_taxon_id","after_observer_cap","request_error"]],
        how="inner",on="inat_taxon_id",validate="one_to_one")
    if len(both)!=42111:raise RuntimeError("Capacity/species census ID universe disagreement")
    if both.request_error.fillna("").astype(str).str.len().gt(0).any():
        raise RuntimeError("Original capacity scan retained request errors")
    both["after_observer_cap"]=pd.to_numeric(both.after_observer_cap,errors="raise").astype(int)
    frame=both.loc[both.after_observer_cap>=100,["inat_taxon_id","species","after_observer_cap"]].copy()
    if len(frame)!=4730:raise RuntimeError(f"Historical U100 count mismatch: {len(frame)}")
    uni=dict(zip(frame.inat_taxon_id,frame.species))
    xgroup={}
    for g,ds in [
        ("historical_p500",p500),
        ("historical_third",third),
        ("new_primary",main),
        ("new_replication",replication)
    ]:
        for row in ds.itertuples(index=False):
            tid=int(row.inat_taxon_id)
            name=str(row.species)
            if tid not in uni or uni[tid]!=name:
                raise RuntimeError(f"Out-of-U100 or stale species identity in {g}: {tid}/{name}")
            if tid in xgroup:raise RuntimeError("Unexpected overlap among frozen selections")
            xgroup[tid]=g
    parent_ids=set(p100.inat_taxon_id.astype(int))
    if len(parent_ids)!=3730 or parent_ids != set(uni)-{
        int(k) for k in uni if k not in parent_ids
    }:
        raise RuntimeError("P100 source invalid")
    # P100 must comprise all four post-legacy pools; its complement is exactly
    # the original 1000 legacy resource. Ensures 4730 complete universe.
    old=set(uni)-parent_ids
    if len(old)!=1000:raise RuntimeError("Old legacy identity complement is not 1000")
    if set(xgroup)!=parent_ids:
        raise RuntimeError("New+P500+third do not partition the 3730 P100 species exactly")
    for tid in old:xgroup[tid]="historical_discovery_reserve"
    frame["allocation_group"]=frame.inat_taxon_id.map(xgroup)
    if frame.allocation_group.isna().any():raise RuntimeError("An original U100 taxon was not allocated")
    v=frame.allocation_group.value_counts().to_dict()
    if v!=GROUPS:raise RuntimeError(f"Allocation groups mismatch: {v}")
    return frame.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)


def main()->None:
    p=argparse.ArgumentParser()
    for n in ("species","capacity","p100","p500","third","primary","replication"):
        p.add_argument("--"+n,required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    frozen={
        "species":SPECIES_HASH,"capacity":CAPACITY_HASH,
        "p100":P100_HASH,"p500":P500_HASH,"third":THIRD_HASH,
        "primary":MAIN_HASH,"replication":REPLICATION_HASH
    }
    paths={k:getattr(args,k) for k in frozen}
    for k,path in paths.items():
        if sha(path)!=frozen[k]:raise RuntimeError(f"Historical {k} SHA256 drift")
    sp=pd.read_csv(args.species,low_memory=False)
    cap=pd.read_csv(args.capacity,low_memory=False)
    r=group_allocation(
        sp,cap,
        read_ids(args.p100,P100_HASH,3730),
        read_ids(args.p500,P500_HASH,500),
        read_ids(args.third,THIRD_HASH,500,"\t"),
        read_ids(args.primary,MAIN_HASH,2000,"\t"),
        read_ids(args.replication,REPLICATION_HASH,730,"\t")
    )
    args.outdir.mkdir(parents=True,exist_ok=True)
    fname=args.outdir/"entire_4730_species_allocation_ledger.csv"
    r.to_csv(fname,index=False,lineterminator="\n")
    result={
        "schema":"fcp_all4730_species_metadata_census_v1",
        "status":"HISTORICAL_U100_FULL_SPECIES_UNIVERSE_VERIFIED_NO_NEW_COLOUR",
        "original_metadata_discovery_species":42111,
        "u100_original_capacity_species":len(r),
        "group_counts":{k:int((r.allocation_group==k).sum()) for k in GROUPS},
        "historical_exposed_or_selected_species":2000,
        "all_new_unallocated_species_planned_for_acquisition":2730,
        "prospective_primary_species":2000,
        "prospective_separate_replication_species":730,
        "complete_taxon_id_and_name_uniqueness":bool(r.inat_taxon_id.is_unique and r.species.is_unique),
        "new_photo_pixels_opened":False,
        "new_biological_confirmation_result_available":False,
        "capacity_scan_not_new_image_measurement":True,
        "source_sha256":frozen,
        "result_ledger_sha256":sha(fname),
        "opportunity_definition":"historical iNaturalist metadata eligible photos after observer cap >=100",
        "outcome_firewall":OUTCOME_FIREWALL,
        "inference_boundaries":[
            "4730 is metadata capacity, not 4730 species actually photographed now",
            "historical P500 and third each included a species short of 100 fresh image records",
            "only the 2730 unexposed candidate species are prospective new photo sources",
            "reproductive fitness, gene variants, and biochemical white state not measured"
        ],
        "confirmatory_decisions_changed":False
    }
    (args.outdir/"result.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
