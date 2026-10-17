#!/usr/bin/env python3
"""One-shot 4730 U100 expansion: exact original ROI-v4 and historical inference.

The extra 2733 species are outcome-blind original H9 metadata candidates. Only
those with an authenticated EXACT 100 original photo-ID gate are measured.
Every photo must have exactly one old classifier terminal including failures.
No source taxonomy, geography or original 149900-pixel classification is altered.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    ACQUISITION_FIELDS, WORKER_FIELDS, semantic_shard, compute_partition,
    reassemble_complete_measurement,
)
from fcp_pipeline.third_cohort_prospective_execution_gate import measurement_id, measurement_batch

from scripts.analysis.run_fcp_wave4_original_roi_and_statistics_20261017 import (
    analyze_original_methods, original_spatial_primary_50km_n,
    FRACTIONS, MORPHS,
)

SOURCE_RUN_ID = 38996743273
SOURCE_HEAD = "61abc4f71e74db3a6633a790e37d1f30e8e3bfcc"
N_BATCHES = 2
N_SHARDS = 32
N_PARTS = 4
TOTAL_PARTS = 256
TARGET = 100
SOURCE_SCHEMA = "fcp_u100_all4730_original_h9_only_v1"
FIREWALL_SCHEMA = "fcp_all4730_original_v4_blind_firewall_v1"
TERMINAL_SCHEMA = "fcp_all4730_new_original_method_pixel_terminal_v1"


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def authenticated_full_u100(source:Path)->tuple[dict,pd.DataFrame]:
    receipt_path=source/"all4730_metadata_result.json"
    photos_path=source/"remaining_U100_authorized_metadata.csv.gz"
    audit_path=source/"full_U100_4730_status_after_H9.csv"
    for p in (receipt_path,photos_path,audit_path):
        if not p.is_file(): raise RuntimeError("all4730 source metadata missing: "+str(p))
    receipt=json.loads(receipt_path.read_text())
    if (receipt.get("schema")!=SOURCE_SCHEMA
        or receipt.get("stage")!="all_4730_H9_metadata_terminal"
        or receipt.get("U100_identity_species")!=4730
        or receipt.get("old_measured_species")!=1499
        or receipt.get("wave4_original_photo_metadata_species")!=498
        or receipt.get("new_attempted_species")!=2733
        or receipt.get("verified_against_all_4_complete_shards") is not True
        or receipt.get("metadata_gate_pass") is not True
        or receipt.get("old_149900_method_and_results_changed") is not False
        or receipt.get("photo_pixels_opened") is not False):
        raise RuntimeError("unverified 4730 original H9 photo metadata source")
    for p in (photos_path,audit_path):
        if sha(p)!=receipt.get("outputs_sha256",{}).get(p.name):
            raise RuntimeError("original H9 source photo SHA drift: "+p.name)
    photo=pd.read_csv(photos_path,low_memory=False)
    n=int(receipt["new_complete_100_photo_species"])
    m=int(receipt["new_original_photo_IDs"])
    required={"inat_taxon_id","species","query_species","photo_id","observation_id",
              "photo_url_large","photo_license","latitude","longitude","observer_id"}
    if required-set(photo.columns):
        raise RuntimeError("original photo metadata missing required fields")
    if (n<1 or n>2733 or m!=n*TARGET or len(photo)!=m
        or photo.photo_id.nunique()!=m or photo.observation_id.nunique()!=m
        or photo.inat_taxon_id.nunique()!=n
        or not photo.groupby("inat_taxon_id").size().eq(TARGET).all()):
        raise RuntimeError("only full 100-photo per U100 taxon may be measured")
    if photo.photo_url_large.isna().any() or photo.photo_license.isna().any():
        raise RuntimeError("original photographic source missing URL/license")
    allx=pd.read_csv(audit_path)
    if (len(allx)!=4730 or allx.inat_taxon_id.nunique()!=4730
        or sum(allx.final_status.eq("original_149900_photo_measured"))!=1499
        or sum(allx.final_status.eq("wave4_49800_photo_id_authorized"))!=498
        or sum(allx.final_status.eq("new_metadata_100_of_100_photo_IDs"))!=n
        or sum(allx.final_status.eq("new_H9_metadata_under_100_or_request_error"))!=2733-n):
        raise RuntimeError("U100 attempted and successful identity statuses drift")
    return receipt,photo


def firewall(source:Path,dest:Path)->dict:
    receipt,photo=authenticated_full_u100(source)
    x=photo.reset_index(drop=True).copy()
    x.insert(0,"measurement_id",[measurement_id(v) for v in x.photo_id])
    x.insert(1,"measurement_batch",[measurement_batch(v) for v in x.measurement_id])
    if x.measurement_id.nunique()!=len(x) or not set(x.measurement_batch).issubset({0,1}):
        raise RuntimeError("original blind measurement IDs/batch changed")
    x["species"]=x["query_species"].astype(str)
    worker=pd.DataFrame({
        "measurement_id":x.measurement_id.astype(str),
        "image_filename":x.measurement_id.astype(str)+".jpg",
        "photo_license":x.photo_license.astype(str),
    })
    acquisition=pd.DataFrame({
        "measurement_id":worker.measurement_id,
        "image_filename":worker.image_filename,
        "photo_url_large":x.photo_url_large.astype(str),
        "photo_license":worker.photo_license,
    })
    if tuple(worker.columns)!=WORKER_FIELDS or tuple(acquisition.columns)!=ACQUISITION_FIELDS:
        raise RuntimeError("original location-blind worker interface changed")
    dest.mkdir(parents=True,exist_ok=True)
    (dest/"sealed_keys").mkdir(exist_ok=True)
    key=dest/"sealed_keys/metadata_join_key.csv"
    x.drop(columns=["photo_url_large"]).to_csv(key,index=False,lineterminator="\n")
    assignments=[]
    for b in range(N_BATCHES):
        folder=dest/f"batch_{b}"
        (folder/"worker_packet").mkdir(parents=True,exist_ok=True)
        (folder/"sealed_keys").mkdir(parents=True,exist_ok=True)
        mask=x.measurement_batch.eq(b)
        w=worker.loc[mask].reset_index(drop=True)
        a=acquisition.loc[mask].reset_index(drop=True)
        w.to_csv(folder/"worker_packet/measurement_manifest.csv",index=False,lineterminator="\n")
        a.to_csv(folder/"sealed_keys/acquisition_key.csv",index=False,lineterminator="\n")
        for mid in w.measurement_id:
            assignments.append((mid,b,semantic_shard(mid,N_SHARDS),compute_partition(mid,N_PARTS)))
    if len(assignments)!=len(x):
        raise RuntimeError("original 256 worker partition coverage incomplete")
    pd.DataFrame(assignments,columns=["measurement_id","batch","semantic_shard","compute_partition"]).to_csv(
        dest/"partition_assignments.csv",index=False,lineterminator="\n")
    out={
        "schema":FIREWALL_SCHEMA,"source_run_id":SOURCE_RUN_ID,
        "source_expected_commit":SOURCE_HEAD,"all_original_U100_taxa":4730,
        "old_frozen_photo_rows":149900,"wave4_original_source_photo_rows":49800,
        "new_100photo_species":int(receipt["new_complete_100_photo_species"]),
        "new_photo_rows":len(x),"new_under_100_species":int(receipt["new_unmet_100_photo_species"]),
        "metadata_photo_SHA256":receipt["outputs_sha256"]["remaining_U100_authorized_metadata.csv.gz"],
        "metadata_join_SHA256":sha(key),"partitions":TOTAL_PARTS,
        "blind_worker_columns":list(WORKER_FIELDS),"uses_unchanged_historical_worker":True,
        "coordinates_or_taxa_in_pixel_worker":False,"photo_pixels_opened":False,
        "existing_149900_results_modified":False,
    }
    (dest/"measurement_firewall_manifest.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    return out


def reassemble(firewall_root:Path,parts:Path,dest:Path)->dict:
    record=json.loads((firewall_root/"measurement_firewall_manifest.json").read_text())
    if (record.get("schema")!=FIREWALL_SCHEMA or record.get("partitions")!=TOTAL_PARTS
        or record.get("coordinates_or_taxa_in_pixel_worker") is not False):
        raise RuntimeError("untrusted 4730 original precolour firewall")
    worker=pd.concat([
        pd.read_csv(firewall_root/f"batch_{b}/worker_packet/measurement_manifest.csv",dtype=str).fillna("")
        for b in range(N_BATCHES)
    ],ignore_index=True)
    key_path=firewall_root/"sealed_keys/metadata_join_key.csv"
    join=pd.read_csv(key_path,dtype={"measurement_id":str},low_memory=False)
    if len(join)!=len(worker) or len(join)!=record["new_photo_rows"] or sha(key_path)!=record["metadata_join_SHA256"]:
        raise RuntimeError("full U100 original metadata/worker SHA changed")
    files={}
    for path in parts.glob("b*_partition_s*_p*"):
        match=re.fullmatch(r"b(\d+)_partition_s(\d{2})_p(\d{2})\.(csv|json)",path.name)
        if match is None:
            raise RuntimeError("unknown original worker receipt "+path.name)
        b,s,p,ext=int(match[1]),int(match[2]),int(match[3]),match[4]
        if not (0<=b<N_BATCHES and 0<=s<N_SHARDS and 0<=p<N_PARTS):
            raise RuntimeError("wrong original worker ID")
        if (b,s,p,ext) in files:raise RuntimeError("duplicate original worker result")
        files[b,s,p,ext]=path
    expected={(b,s,p) for b in range(N_BATCHES) for s in range(N_SHARDS) for p in range(N_PARTS)}
    for ext in ("csv","json"):
        if {v[:3] for v in files if v[3]==ext}!=expected:
            raise RuntimeError("all 256 original classifier partitions must terminate before data join")
    sealed=[]
    for b,s,p in sorted(expected):
        meta=json.loads(files[b,s,p,"json"].read_text())
        if (meta.get("status")!="complete_random_photo_first_terminal_partition"
            or meta.get("semantic_shard")!=s or meta.get("compute_partition")!=p
            or meta.get("source_urls_present") is not False
            or meta.get("species_present") is not False
            or meta.get("coordinates_present") is not False):
            raise RuntimeError("invalid old-style terminal ROI packet")
        sealed.append(pd.read_csv(files[b,s,p,"csv"],dtype={"measurement_id":str}).fillna(""))
    assembled=reassemble_complete_measurement(
        sealed,worker,join,expected_partition_receipts=TOTAL_PARTS)
    measured=assembled.joined_photos
    if len(measured)!=len(join) or measured.measurement_id.nunique()!=len(measured):
        raise RuntimeError("no missing, extra, or repeated measurement terminals allowed")
    measured["global_classifiable"]=measured.measurement_status.eq("classified_four_state_morph") & measured.morph.isin(MORPHS)
    if set(FRACTIONS)-set(measured.columns):
        raise RuntimeError("the original 9 fixed colour-fraction columns are missing")
    for col in FRACTIONS:
        measured[col]=pd.to_numeric(measured[col],errors="raise")
    denom=measured.loc[measured.global_classifiable,FRACTIONS].sum(axis=1)
    if len(denom) and not np.allclose(denom.to_numpy(float),1.0,atol=1e-8):
        raise RuntimeError("original nine-colour normalization changed")
    support=measured.groupby(["species","inat_taxon_id"]).agg(
        raw_photos=("measurement_id","size"),
        classifiable_photos=("global_classifiable","sum")).reset_index()
    if (len(support)!=record["new_100photo_species"] or
        not support.raw_photos.eq(TARGET).all()):
        raise RuntimeError("new original source species lost a full 100-photo denominator")
    support["measurement_evaluable"]=support.classifiable_photos.ge(40)
    dest.mkdir(parents=True,exist_ok=True)
    outphotos=dest/"all4730_new_original_measured_photos.csv.gz"
    supportpath=dest/"all4730_new_species_original_colour_support.csv"
    measured.to_csv(outphotos,index=False,lineterminator="\n",
                    compression={"method":"gzip","compresslevel":9,"mtime":0})
    support.to_csv(supportpath,index=False,lineterminator="\n")
    out={
        "schema":TERMINAL_SCHEMA,"measurement_photo_rows":len(measured),
        "new_complete_100_photo_species":len(support),
        "new_classifiable_photos":int(measured.global_classifiable.sum()),
        "new_nonclassifiable_photos":int((~measured.global_classifiable).sum()),
        "new_D_support_40photo_species":int(support.measurement_evaluable.sum()),
        "new_under_100_species_prior_H9":record["new_under_100_species"],
        "original_roi_v4_and_efficientsam_unchanged":True,
        "original_frozen_149900_results_changed":False,
        "measurement_terminal_partitions":TOTAL_PARTS,
        "measured_photo_file_sha256":sha(outphotos),
        "new_species_support_sha256":sha(supportpath),
        "source_metadata_sha256":record["metadata_photo_SHA256"],
        "technical_missing_photos_not_replaced":True,
    }
    (dest/"all4730_new_photo_terminal_result.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    return out


def original_inference(measured:Path, terminal:Path, dest:Path)->dict:
    p=json.loads(terminal.read_text())
    if (p.get("schema")!=TERMINAL_SCHEMA or p.get("measured_photo_file_sha256")!=sha(measured)
        or p.get("original_roi_v4_and_efficientsam_unchanged") is not True):
        raise RuntimeError("new original actual ROI colour terminal not source-authenticated")
    # Science uses the ORIGINAL existing worker and 50km/H1/W functions verbatim.
    old=analyze_original_methods(measured,dest)
    # Original 50km gate requires enough local photo pairs, not merely
    # n_classifiable>=40. Read the exact original spatial test result.
    n_50km=original_spatial_primary_50km_n(old["spatial_original_method"])
    if old["wave4_measured_table_sha256"]!=sha(measured):
        raise RuntimeError("original unchanged scientific estimator input SHA drift")
    out={
        "schema":"fcp_all4730_same_original_H1_D_50km_W_bulk_new_species_v1",
        "original_stats_called_without_rewriting":True,
        "historical_function":"run_fcp_wave4_original_roi_and_statistics_20261017.analyze_original_methods",
        "old_historical_result_schema":old["schema"],
        "source_measured_sha256":sha(measured),
        "n_new_species":p["new_complete_100_photo_species"],
        "n_new_D_evaluable":p["new_D_support_40photo_species"],
        "n_new_50km_evaluable":n_50km,
        "original_frozen_149900_result_modified":False,
        "original_H1_D":old["h1_original_method"],
        "original_50km_spatial":old["spatial_original_method"],
        "original_H2_W":old["h2_original_W_method"],
        "posthoc_50km_still_posthoc":True,
        "same_data_provider_not_independent_replication":True,
    }
    (dest/"all4730_new_original_method_inference.json").write_text(
        json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    return out


def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="action",required=True)
    a=sub.add_parser("firewall")
    a.add_argument("--source",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    b=sub.add_parser("reassemble")
    b.add_argument("--firewall",type=Path,required=True)
    b.add_argument("--parts",type=Path,required=True)
    b.add_argument("--out",type=Path,required=True)
    c=sub.add_parser("infer")
    c.add_argument("--measured",type=Path,required=True)
    c.add_argument("--terminal",type=Path,required=True)
    c.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=firewall(args.source,args.out) if args.action=="firewall" else (
        reassemble(args.firewall,args.parts,args.out) if args.action=="reassemble"
        else original_inference(args.measured,args.terminal,args.out))
    print("FCP_ALL4730_ORIGINAL_"+args.action.upper()+"_COMPLETE")
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
