#!/usr/bin/env python3
"""One-pass prospective metadata draw for ALL 2730 previously unused U100 taxa.

Historical iNaturalist metadata sampling logic pinned from pre-outcome code.
Each rank is assigned to EXACTLY one fixed shard. No image pixels or colour
outcomes are opened. A technically failed shard reports failure, not replacement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

SELECTED_HASH="e1301095533dfa755c1588cc550522bf31aabaef632083e265a138f0e3af6d2a"
REPLICATION_HASH="dbfaf954c1b82c569e9495006105639ad1d27b356a20b0ee3cfb5b4e240c00e8"
OLD_EXCLUSIONS={
 "random_photo_first_h9_exclusion_ledger_v1.csv":"f9a6894740e9974399c055f92cba237be8ada41707d84e1807ba61b902c91b99",
 "random_photo_first_h9_fresh_metadata_v1.csv":"111d0f964618c0c3df749a6e4bd29f834214cda5d7a2d9bda7d43cdc9dbf4c6f",
 "rgfca_42111_species_breadth_measured.csv.gz":"38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605",
 "polymorphism_h2_p500_candidate_metadata_v1.csv.gz":"a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595",
 "polymorphism_h2_third_cohort_candidate_metadata_v1.csv.gz":"adec40e29e347b035872f2add95b67011906cb74e28510e6260ed1edbd075711",
}
PREVIOUS_EXCLUSION_UNION=228461  # 178462 legacy+P500 plus 49999 fresh third; disjoint required
EXACT_SELECTED=2000
EXACT_REPLICATION=730
SHARDS=20
PHOTO_TARGET=100
MIN_NEW_FULL100_PRIMARY=1000
OBSERVER_CAP=2
PER_PAGE=200
MAX_ACCURACY_M=5000
REQUEST_INTERVAL_SECONDS=1.05
REQUEST_TIMEOUT_SECONDS=45.0
REQUEST_RETRIES=0
ALLOWED_LICENSES=("cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa")


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for ch in iter(lambda:f.read(1<<20),b""):h.update(ch)
    return h.hexdigest()


def read_allocation(primary:Path, replication:Path)->pd.DataFrame:
    if sha(primary)!=SELECTED_HASH or sha(replication)!=REPLICATION_HASH:
        raise RuntimeError("Original untouched species allocation SHA256 drift")
    a=pd.read_csv(primary,sep="\t");b=pd.read_csv(replication,sep="\t")
    if len(a)!=EXACT_SELECTED or len(b)!=EXACT_REPLICATION:
        raise RuntimeError("Unmodified 2000+730 identity counts no longer valid")
    a["allocation_group"]="new_primary"
    b["allocation_group"]="new_replication"
    combined=pd.concat([a,b],ignore_index=True)
    req={"prospective_rank","inat_taxon_id","species","after_observer_cap","allocation_group"}
    if not req.issubset(combined):
        raise RuntimeError(f"Selected identity missing {sorted(req-set(combined))}")
    ranks=pd.to_numeric(combined.prospective_rank,errors="raise").astype(int)
    taxa=pd.to_numeric(combined.inat_taxon_id,errors="raise").astype(int)
    if sorted(ranks.tolist())!=list(range(1,2731)):
        raise RuntimeError("Not every one of the 2730 original selection ranks is present exactly once")
    if len(set(taxa))!=2730 or combined.species.nunique()!=2730:
        raise RuntimeError("Duplicated future species taxon or scientific name")
    if not (pd.to_numeric(combined.after_observer_cap,errors="raise")>=100).all():
        raise RuntimeError("Frozen U100 opportunity contains species with <100 eligible photos")
    return combined.sort_values("prospective_rank",kind="stable").reset_index(drop=True)


def shard_frame(pool:pd.DataFrame, index:int,count:int=SHARDS)->pd.DataFrame:
    if count!=SHARDS or not 0<=index<count:
        raise ValueError("Shard identity/budget cannot be retuned")
    mask=(pool.prospective_rank.astype(int)-1)%count==index
    frame=pool.loc[mask].copy().reset_index(drop=True)
    if not len(frame):
        raise RuntimeError("Invalid empty shard")
    # Preserve eligibility rank in audit after the prior engine sorts taxon ID.
    return frame


def exclusions(source_dir:Path)->tuple[set[int],set[int],dict]:
    obs=set(); photos=set(); audit={}
    for filename,expected_sha in OLD_EXCLUSIONS.items():
        file=source_dir/filename
        if not file.is_file() or sha(file)!=expected_sha:
            raise RuntimeError(f"Previous photo-ID exclusion source unavailable or mismatched: {filename}")
        frame=pd.read_csv(file,usecols=["observation_id","photo_id"],low_memory=False)
        a=set(pd.to_numeric(frame.observation_id,errors="raise").dropna().astype("int64"))
        b=set(pd.to_numeric(frame.photo_id,errors="raise").dropna().astype("int64"))
        obs|=a;photos|=b
        audit[filename]={"sha256":sha(file),"rows":len(frame),
                         "unique_obs":len(a),"unique_photo":len(b)}
    if len(obs)!=PREVIOUS_EXCLUSION_UNION or len(photos)!=PREVIOUS_EXCLUSION_UNION:
        raise RuntimeError(f"Historical union drift: obs={len(obs)} photos={len(photos)} expected={PREVIOUS_EXCLUSION_UNION}")
    return obs,photos,audit


def main()->None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--primary",required=True,type=Path)
    parser.add_argument("--replication",required=True,type=Path)
    parser.add_argument("--exclusion-dir",required=True,type=Path)
    parser.add_argument("--shard-index",required=True,type=int)
    parser.add_argument("--shard-count",required=True,type=int)
    parser.add_argument("--outdir",required=True,type=Path)
    args=parser.parse_args()

    if os.environ.get("GITHUB_ACTIONS")=="true" and os.environ.get("GITHUB_RUN_ATTEMPT","1")!="1":
        raise RuntimeError("Entire original one-pass metadata draw cannot be rerun")
    if os.environ.get("FCP_FRESH_METADATA_ALREADY_SEALED","")=="true":
        raise RuntimeError("Source programme already has a terminal recorded draw")

    pool=read_allocation(args.primary,args.replication)
    group=shard_frame(pool,args.shard_index,args.shard_count)
    obs,photo,source_audit=exclusions(args.exclusion_dir)

    # The source client is exactly the original source-frozen metadata-only
    # sampler used to qualify the previous independent third cohort. We open
    # no pixels even though photo URL strings are stored for future use.
    client=InaturalistObservationClient(
        user_agent=f"zuizui0223-fcp-all4730-meta-shard-{args.shard_index}/1.0",
        request_interval_seconds=REQUEST_INTERVAL_SECONDS,
        timeout_seconds=REQUEST_TIMEOUT_SECONDS,
        max_retries=REQUEST_RETRIES)
    frozen=freeze_h9_metadata(
        client=client,
        species_frame=group[["inat_taxon_id","species"]],
        exclusion_observation_ids=obs,
        exclusion_photo_ids=photo,
        per_page=PER_PAGE,
        observer_cap_n=OBSERVER_CAP,
        fixed_raw_photos=PHOTO_TARGET,
        maximum_positional_accuracy_m=MAX_ACCURACY_M,
        allowed_photo_licenses=ALLOWED_LICENSES)

    audit=frozen.species_audit.merge(
        group[["inat_taxon_id","prospective_rank","allocation_group"]],
        on="inat_taxon_id",how="left",validate="one_to_one")
    audit=audit.sort_values("prospective_rank",kind="stable").reset_index(drop=True)
    if len(audit)!=len(group) or audit.inat_taxon_id.nunique()!=len(group):
        raise RuntimeError("One or more species lost in metadata API opportunity audit")
    records=frozen.observations.copy()
    if len(records):
        records=records.merge(
            group[["inat_taxon_id","prospective_rank","allocation_group"]],
            on="inat_taxon_id",how="left",validate="many_to_one")
        if records.observation_id.duplicated().any() or records.photo_id.duplicated().any():
            raise RuntimeError("Within-shard reused photo or observation IDs")
        if set(records.observation_id.astype(int))&obs or set(records.photo_id.astype(int))&photo:
            raise RuntimeError("Previously used photo/observation IDs leaked")
        counts=records.groupby("inat_taxon_id").size()
        if int(counts.max())>PHOTO_TARGET:
            raise RuntimeError("New species exceeds predeclared max 100 photos")
        records=records.sort_values(
            ["prospective_rank","h9_selection_order","observation_id","photo_id"],
            kind="stable").reset_index(drop=True)
    args.outdir.mkdir(parents=True,exist_ok=True)
    audits=args.outdir/f"species_audit_shard_{args.shard_index:02}.csv.gz"
    photos=args.outdir/f"metadata_photos_shard_{args.shard_index:02}.csv.gz"
    audit.to_csv(audits,index=False,lineterminator="\n",
                 compression={"method":"gzip","compresslevel":9,"mtime":0})
    records.to_csv(photos,index=False,lineterminator="\n",
                   compression={"method":"gzip","compresslevel":9,"mtime":0})
    errors=int(audit.request_error.fillna("").astype(str).str.len().gt(0).sum())
    n_full=int((audit.retained.astype(int)==PHOTO_TARGET).sum())
    report={
        "schema":"fcp_all4730_new2730_metadata_shard_v1",
        "status":"ONE_SHOT_METADATA_DRAW_COMPLETE_BEFORE_IMAGE_PIXELS",
        "source_selected_species_sha256":SELECTED_HASH,
        "source_replication_species_sha256":REPLICATION_HASH,
        "shard_index":args.shard_index,"shard_count":SHARDS,
        "species_query_attempts":len(group),
        "species_eligible100":n_full,
        "source_species_with_any_photos":int((audit.retained.astype(int)>0).sum()),
        "total_metadata_photo_rows":len(records),
        "request_errors":errors,
        "query_contract":{
            "one_random_page_per_taxon":True,
            "per_page":PER_PAGE,"observer_cap":OBSERVER_CAP,
            "target_photo_records_per_species":PHOTO_TARGET,
            "request_interval_seconds":REQUEST_INTERVAL_SECONDS,
            "request_retries":REQUEST_RETRIES,
            "max_positional_accuracy_m":MAX_ACCURACY_M,
            "allowed_photo_licenses":list(ALLOWED_LICENSES)
        },
        "group_species":{z:int((audit.allocation_group==z).sum())
                         for z in ("new_primary","new_replication")},
        "group_full100":{z:int(((audit.allocation_group==z)&(audit.retained==PHOTO_TARGET)).sum())
                         for z in ("new_primary","new_replication")},
        "source_exclusion_id_union_size":len(obs),
        "exclusion_files":source_audit,
        "source_manifest_sha256":sha(args.primary),
        "observed_metadata_sha256":sha(photos),
        "species_audit_sha256":sha(audits),
        "outcome_firewall":{
            "new_photograph_pixels_opened":False,
            "flower_colour_morph_labels_computed":False,
            "candidate_colour_used_for_selection":False,
            "old_results_consulted_for_photo_replacement":False,
            "biological_test_run":False
        },
        "adjudication_note":"This is one immutable metadata opportunity draw. Failed or <100 species are NOT replaced; outcome code must never use photo pixel/colour labels from this stage.",
        "confirmatory_decisions_changed":False
    }
    target=args.outdir/f"receipt_shard_{args.shard_index:02}.json"
    target.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()
