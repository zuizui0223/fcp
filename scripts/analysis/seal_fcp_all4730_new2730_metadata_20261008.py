#!/usr/bin/env python3
"""Seal the complete 2730-species metadata-only iNaturalist first draw.

Input = exactly twenty frozen deterministic source-species shards. Every
newly selected FCP species appears EXACTLY once, including API failures.
This script cannot open images and cannot recompute floristic colour.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

NUM_SHARDS=20
EXPECTED_SPECIES=2730
EXPECTED_PRIMARY=2000
EXPECTED_REPLICATION=730
PHOTO_TARGET=100
MAX_REQUEST_ERRORS_FRACTION=0.05
CAPACITY_PRIMARY_GATE_FULL100=1000
EXPECTED_HASH_MAIN="e1301095533dfa755c1588cc550522bf31aabaef632083e265a138f0e3af6d2a"
EXPECTED_HASH_REPLICATION="dbfaf954c1b82c569e9495006105639ad1d27b356a20b0ee3cfb5b4e240c00e8"


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()


def collected_receipts(root:Path)->tuple[pd.DataFrame,pd.DataFrame,list[dict]]:
    aud=[];photo=[];receipts=[]
    for i in range(NUM_SHARDS):
        tag=f"{i:02}"
        p=root/f"receipt_shard_{tag}.json"
        af=root/f"species_audit_shard_{tag}.csv.gz"
        pf=root/f"metadata_photos_shard_{tag}.csv.gz"
        if not p.exists() or not af.exists() or not pf.exists():
            raise RuntimeError(f"Shard {i} input incomplete; cannot seal all species")
        j=json.loads(p.read_text())
        if j.get("schema")!="fcp_all4730_new2730_metadata_shard_v1" or j.get("shard_index")!=i or j.get("shard_count")!=NUM_SHARDS:
            raise RuntimeError("Unexpected or duplicate shard identity")
        if j.get("source_selected_species_sha256")!=EXPECTED_HASH_MAIN or j.get("source_replication_species_sha256")!=EXPECTED_HASH_REPLICATION:
            raise RuntimeError("Prospective selected taxon identity changed during metadata draw")
        if j.get("observed_metadata_sha256")!=sha(pf) or j.get("species_audit_sha256")!=sha(af):
            raise RuntimeError("Source metadata or species audit file no longer matches source shard receipt")
        if any(j.get("outcome_firewall",{}).values()):
            raise RuntimeError("Source metadata shard has improperly opened photographic colour")
        a=pd.read_csv(af,low_memory=False)
        f=pd.read_csv(pf,low_memory=False)
        if j["species_query_attempts"]!=len(a) or j["total_metadata_photo_rows"]!=len(f):
            raise RuntimeError("Source rows or query attempts inconsistent with shard receipt")
        if len(f):
            if f.photo_id.duplicated().any() or f.observation_id.duplicated().any():
                raise RuntimeError("Source shard has repeated photo IDs")
        receipts.append(j);aud.append(a);photo.append(f)
    a=pd.concat(aud,ignore_index=True)
    f=pd.concat(photo,ignore_index=True)
    return a,f,receipts


def verify_unique_opportunities(a:pd.DataFrame,f:pd.DataFrame,receipts:list[dict])->dict:
    if len(a)!=EXPECTED_SPECIES or a.inat_taxon_id.nunique()!=EXPECTED_SPECIES or a.prospective_rank.nunique()!=EXPECTED_SPECIES:
        raise RuntimeError("Not every new species was attempted exactly once")
    ranks=sorted(pd.to_numeric(a.prospective_rank,errors="raise").astype(int))
    if ranks!=list(range(1,EXPECTED_SPECIES+1)):
        raise RuntimeError("Different/new source taxon ranks substituted")
    groups=a.allocation_group.value_counts().to_dict()
    if groups!={"new_primary":EXPECTED_PRIMARY,"new_replication":EXPECTED_REPLICATION}:
        raise RuntimeError(f"Missing primary or replication region: {groups}")
    if len(f):
        if f.observation_id.duplicated().any() or f.photo_id.duplicated().any():
            raise RuntimeError("Distinct taxon shards reused an iNaturalist observation or photo ID")
        # Every photo row must identify an attempted species with exact group.
        pairs=set(zip(a.inat_taxon_id.astype(str),a.allocation_group))
        if not set(zip(f.inat_taxon_id.astype(str),f.allocation_group)).issubset(pairs):
            raise RuntimeError("Photo metadata contains a species never attempted")
        counts=f.groupby("inat_taxon_id").size()
        ac=a.set_index("inat_taxon_id").retained.astype(int).sort_index()
        if not ac.reindex(counts.index).eq(counts).all():
            raise RuntimeError("Saved photo counts do not reproduce species-level opportunity audits")
        if int(counts.max())>PHOTO_TARGET:
            raise RuntimeError("More than frozen one hundred photo metadata rows per new species")
    errors=int(sum(int(x["request_errors"]) for x in receipts))
    assert errors==int(a.request_error.fillna("").astype(str).str.len().gt(0).sum())
    full=(a.retained.astype(int)==PHOTO_TARGET)
    nprimary=int(((a.allocation_group=="new_primary")&full).sum())
    nrep=int(((a.allocation_group=="new_replication")&full).sum())
    threshold_pass=(errors/EXPECTED_SPECIES<=MAX_REQUEST_ERRORS_FRACTION
                    and nprimary>=CAPACITY_PRIMARY_GATE_FULL100)
    summary={
        "n_species_attempted":EXPECTED_SPECIES,
        "n_primary_attempted":EXPECTED_PRIMARY,
        "n_replication_attempted":EXPECTED_REPLICATION,
        "n_primary_full100":nprimary,
        "n_replication_full100":nrep,
        "n_total_full100":nprimary+nrep,
        "n_any_photo_species":int((a.retained.astype(int)>0).sum()),
        "n_metadata_photo_records":int(len(f)),
        "n_request_errors":errors,
        "request_error_rate":errors/EXPECTED_SPECIES,
        "met_capacity_and_transport_gate":bool(threshold_pass),
        "n_at_least_40_photo_metadata":int((a.retained.astype(int)>=40).sum()),
        "n_at_least_40_primary":int(((a.allocation_group=="new_primary")&(a.retained.astype(int)>=40)).sum()),
        "n_at_least_40_replication":int(((a.allocation_group=="new_replication")&(a.retained.astype(int)>=40)).sum())
    }
    return summary


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--shard-root",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    audit,photos,receipts=collected_receipts(args.shard_root)
    m=verify_unique_opportunities(audit,photos,receipts)
    args.outdir.mkdir(parents=True,exist_ok=True)
    audit=audit.sort_values("prospective_rank",kind="stable").reset_index(drop=True)
    if len(photos):
        photos=photos.sort_values(["prospective_rank","h9_selection_order","observation_id","photo_id"],
                                  kind="stable").reset_index(drop=True)
    audit_out=args.outdir/"all2730_new_species_metadata_audit.csv.gz"
    p_out=args.outdir/"all2730_fresh_photo_metadata.csv.gz"
    audit.to_csv(audit_out,index=False,lineterminator="\n",
                 compression={"method":"gzip","compresslevel":9,"mtime":0})
    photos.to_csv(p_out,index=False,lineterminator="\n",
                  compression={"method":"gzip","compresslevel":9,"mtime":0})
    hashes={f"receipt_shard_{int(j['shard_index']):02}":hashlib.sha256(
        json.dumps(j,sort_keys=True).encode()).hexdigest() for j in receipts}
    out={
        "schema":"fcp_all4730_new2730_metadata_durable_denominator_v1",
        "date_jst":"2026-10-08",
        "status":"FROZEN_ONE_SHOT_NEW2730_OBSERVATION_METADATA_NO_PHOTO_PIXELS",
        "historical_u100_species":4730,
        "previously_allocated_species":2000,
        **m,
        "shards":NUM_SHARDS,
        "primary_species_manifest_SHA256":EXPECTED_HASH_MAIN,
        "replication_species_manifest_SHA256":EXPECTED_HASH_REPLICATION,
        "shard_receipt_canonical_json_sha256":hashes,
        "all_photo_metadata_sha256":sha(p_out),
        "all_species_audit_sha256":sha(audit_out),
        "biological_outcome_opened":False,
        "biological_colour_classification_authorized_by_this_run":False,
        "next_stage":"Independent freeze, synthetic qualification and one authorized biological photo measurement; no species substitution or refill",
        "capacity_failure_is_not_biological_monomorphism":True,
        "confirmatory_decisions_changed":False
    }
    (args.outdir/"result.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
