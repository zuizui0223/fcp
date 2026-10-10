#!/usr/bin/env python3
"""Full 405 original-photo byte provenance from exact archived blind review queue.

Only iNaturalist openly licensed S3 public dataset; no raw images retained.
Neither source colour labels nor case group enters photo fetching/scoring.
"""
from __future__ import annotations
import argparse
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import time
import pandas as pd

PILOT=Path(__file__).resolve().parent/"probe_fcp_blinded_open_photo_bytes_20261010.py"
SPEC=importlib.util.spec_from_file_location("blind_openpilot",PILOT)
mod=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(mod)
EXPECTED_BLIND_SHA="8cc44a71dfb64601e23f68b68c25e8cbc534a332e0a17ff51a51bd8737a35c27"
EXPECTED_KEY_SHA="ef0275fedf7ab4e3e0aa5309fd80c9a9aee97de733b79a93bc2137a04f6f4b5c"
N_CASES=405

def digest_file(path):return sha256(Path(path).read_bytes()).hexdigest()

def manifest(blind,key):
    if digest_file(blind)!=EXPECTED_BLIND_SHA: raise RuntimeError("Original blinded 405 photo queue SHA mismatch")
    if digest_file(key)!=EXPECTED_KEY_SHA: raise RuntimeError("Original separately sealed review key SHA mismatch")
    d=pd.read_csv(blind,dtype=str)
    k=pd.read_csv(key,dtype=str)
    if len(d)!=N_CASES or len(k)!=N_CASES: raise RuntimeError("405 reviewer cases missing")
    if "original_algorithm_colour_label" in d or "kind" in d:
        raise RuntimeError("Blinding violated by review CSV")
    if not {"audit_case_id","photo_id"}.issubset(d):
        raise RuntimeError("Missing original blinded case identity")
    if not {"audit_case_id","photo_id","image_sha256"}.issubset(k):
        raise RuntimeError("Missing original separate image-byte hash")
    source=k[["audit_case_id","photo_id","image_sha256"]].copy()
    subset=d[["audit_case_id","photo_id"]].merge(source,on=["audit_case_id","photo_id"],
             how="inner",validate="one_to_one")
    if len(subset)!=N_CASES or subset.photo_id.duplicated().any() or subset.audit_case_id.duplicated().any():
        raise RuntimeError("Blinded-to-key exact original source photo IDs fail")
    if not subset.image_sha256.str.fullmatch("[0-9a-f]{64}").all():
        raise RuntimeError("Original per-photo byte SHA256 missing")
    return subset.sort_values("audit_case_id",kind="stable").reset_index(drop=True)

def run(blind,key,out,*,acquirer=mod.acquire_one,sleep=time.sleep):
    x=manifest(blind,key)
    p=Path(out);p.mkdir(parents=True,exist_ok=True)
    records=[]
    for idx,row in enumerate(x.itertuples(index=False),1):
        got,errors=acquirer(row.photo_id,row.image_sha256)
        rec={"audit_case_id":row.audit_case_id,"photo_id":row.photo_id,
             "historical_image_sha256":row.image_sha256,
             **got, "n_failed_variants":len(errors)}
        records.append(rec)
        if idx%25==0 or idx==len(x):
            z=Counter(t["status"] for t in records)
            print(f"recovered_and_checked={idx}/{N_CASES} states={dict(z)}",flush=True)
        sleep(.14)
    d=pd.DataFrame(records)
    d.to_csv(p/"photo_byte_provenance_and_whole_frame_qc.csv",index=False)
    z=Counter(d.status)
    r={"schema":"fcp_original_405_blinded_openphoto_provenance_v1",
       "n_original_blinded_photo_ids":N_CASES,
       "n_exact_historical_image_bytes":z["byte_match_verified"],
       "n_retrieved_original_photo_id_but_different_bytes":z["photo_id_retrieved_bytes_differ"],
       "n_unavailable_on_open_data_host":z["unavailable_on_open_data_host"],
       "statuses":dict(z),
       "original_blind_source_sha256":EXPECTED_BLIND_SHA,
       "original_sealed_key_sha256":EXPECTED_KEY_SHA,
       "source_colour_labels_seen_by_photo_fetcher":False,
       "image_pixels_retained_in_result":False,
       "image_licenses":"Only iNaturalist open-data S3; source image attribution still required for any redistribution",
       "label_reannotation_completed":False,
       "inferential_status":"EXACT_SOURCE_PIXEL_AVAILABILITY_ONLY",
       "not_claims":[
           "Open-data availability is not representative of all photo licenses or all flowering plants.",
           "Byte match validates source image identity, not biological flower organ, colour, or segmentation.",
           "No error-rate estimate or causal biological selection inference."
       ]}
    if sum(z.values())!=N_CASES:raise RuntimeError("One or more source photos unaccounted")
    (p/"result.json").write_text(json.dumps(r,indent=2)+"\n")
    return r

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--blind",required=True,type=Path)
    pa.add_argument("--sealed-key",required=True,type=Path)
    pa.add_argument("--outdir",required=True,type=Path)
    a=pa.parse_args()
    print(json.dumps(run(a.blind,a.sealed_key,a.outdir),indent=2))

if __name__=="__main__":main()
