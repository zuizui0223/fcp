#!/usr/bin/env python3
"""Bounded iNaturalist OPEN-DATA photo-byte recovery pilot; no colour annotation.

Only open-data S3 host is permitted. Never upload/release original image bytes.
Original image SHA256 is historical measurement evidence; a changed rendition
does not prove image substitution, but must never be called byte-verified.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from PIL import Image, ImageOps

SCHEMA="fcp_original_blind_openphoto_pilot_20261010_v1"
HOST="https://inaturalist-open-data.s3.amazonaws.com"
IMAGE_SIZE_VARIANTS=("large","original")
EXTENSIONS=("jpeg","jpg","png")
MAX_BYTES=8*1024*1024
PILOT_SHA256="06191b36eca9862f6241b5110ae13aee0ff7e7414a5caaaa1025daffca33eadd"
USER_AGENT="fcp-blind-photo-pilot/1.0 (github.com/zuizui0223/fcp; research provenance)"
SUCCESS_STATES=("byte_match_verified", "photo_id_retrieved_bytes_differ")


def candidate_urls(photo_id):
    if not str(photo_id).isdigit() or int(photo_id)<=0:
        raise ValueError("Invalid photo ID, refusing arbitrary URLs")
    return [f"{HOST}/photos/{int(photo_id)}/{sz}.{ext}"
            for sz in IMAGE_SIZE_VARIANTS for ext in EXTENSIONS]


def load_manifest(path):
    payload=Path(path).read_bytes()
    if sha256(payload).hexdigest()!=PILOT_SHA256:
        raise RuntimeError("24-photo manifest byte identity mismatch")
    d=pd.read_csv(BytesIO(payload),dtype=str)
    if list(d.columns)!=["audit_case_id","photo_id","frozen_original_image_sha256"]:
        raise RuntimeError("pilot schema/order drift")
    if len(d)!=24 or d.audit_case_id.duplicated().any() or d.photo_id.duplicated().any():
        raise RuntimeError("pilot size/unique case/photo drift")
    if not d.audit_case_id.str.fullmatch(r"FCQ-[0-9a-f]{16}").all():
        raise RuntimeError("unexpected blind audit IDs")
    if not d.photo_id.str.fullmatch(r"[0-9]+").all():
        raise RuntimeError("unexpected photo IDs")
    if not d.frozen_original_image_sha256.str.fullmatch(r"[0-9a-f]{64}").all():
        raise RuntimeError("missing original historical image hashes")
    return d


def pixel_quality(data):
    with Image.open(BytesIO(data)) as f:
        f.verify()
    with Image.open(BytesIO(data)) as raw:
        im=ImageOps.exif_transpose(raw).convert("RGB")
        im.thumbnail((384,384),Image.Resampling.BILINEAR)
        a=np.asarray(im,dtype=np.float32)/255.
        if a.ndim!=3 or a.shape[2]!=3: raise RuntimeError("RGB decode failed")
        lum=.2126*a[:,:,0]+.7152*a[:,:,1]+.0722*a[:,:,2]
        sat=a.max(axis=2)-a.min(axis=2)
        return {"width_qa_px":int(im.width),"height_qa_px":int(im.height),
                "bright_clip_fraction_whole_photo":round(float(np.mean(lum>=.98)),6),
                "dark_clip_fraction_whole_photo":round(float(np.mean(lum<=.02)),6),
                "whole_photo_rgb_channel_spread_mean":round(float(np.mean(sat)),6)}


def acquire_one(photo_id,historical_digest,opener=urlopen,*,sleep=time.sleep):
    tries=[]
    fallback=None
    for url in candidate_urls(photo_id):
        try:
            req=Request(url,headers={"User-Agent":USER_AGENT,"Accept":"image/jpeg,image/png"})
            with opener(req,timeout=12) as response:
                final_url=response.geturl()
                if not final_url.startswith(HOST+"/photos/"+str(int(photo_id))+"/"):
                    raise RuntimeError("unexpected redirect outside exact open-data photo ID")
                data=response.read(MAX_BYTES+1)
            if len(data)<1024 or len(data)>MAX_BYTES:
                raise RuntimeError("invalid photo size")
            q=pixel_quality(data)
            digest=sha256(data).hexdigest()
            photo={"status":"byte_match_verified" if digest==historical_digest
                           else "photo_id_retrieved_bytes_differ",
                   "retrieved_url":url,
                   "new_bytes_sha256":digest,
                   "bytes":len(data),**q}
            if digest==historical_digest:
                return photo,tries
            if fallback is None: fallback=photo
        except HTTPError as exc:
            tries.append({"candidate":url.rsplit("/",1)[-1],
                          "error_code":f"HTTP_{exc.code}"})
        except (URLError,OSError,RuntimeError,ValueError) as exc:
            tries.append({"candidate":url.rsplit("/",1)[-1],
                          "error_code":type(exc).__name__})
        sleep(.15)
    if fallback:
        return fallback,tries
    return {"status":"unavailable_on_open_data_host",
            "retrieved_url":"","new_bytes_sha256":"",
            "bytes":0},tries


def run(source,outdir,*,opener=urlopen,sleep=time.sleep):
    d=load_manifest(source)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for row in d.itertuples(index=False):
        got,tries=acquire_one(row.photo_id,row.frozen_original_image_sha256,
                              opener,sleep=sleep)
        record={"audit_case_id":row.audit_case_id,
                "photo_id":row.photo_id,
                "original_image_sha256":row.frozen_original_image_sha256,
                **got,
                "attempted_open_data_variants":len(candidate_urls(row.photo_id)),
                "failed_fetch_variants":len(tries)}
        rows.append(record)
        print(f"case={row.audit_case_id} status={record['status']} sha_match={record['status']=='byte_match_verified'}",
              flush=True)
        sleep(.20)
    r=pd.DataFrame(rows)
    r.to_csv(out/"byte_provenance_and_technical_qc.csv",index=False)
    status=Counter(r.status)
    result={
       "schema":SCHEMA,
       "pilot_manifest_sha256":PILOT_SHA256,
       "n_pilot_photo_ids":len(d),
       "n_byte_identical_to_historical_source":status["byte_match_verified"],
       "n_photo_id_retrieved_but_bytes_differ":status["photo_id_retrieved_bytes_differ"],
       "n_not_recovered_from_open_licensed_photo_host":status["unavailable_on_open_data_host"],
       "statuses":dict(status),
       "source_original_colour_labels_read":False,
       "image_pixels_persisted_or_uploaded":False,
       "photo_quality_metric_domain":"whole photograph, not flower/organ segmentation",
       "photo_colour_assignment_performed":False,
       "human_review_performed":False,
       "claim_ceiling":"Pilot image availability/byte-level provenance only; no biological flower colour classification, full 405 panel recovery, or empirical classifier error rate.",
       "per_case_file":"byte_provenance_and_technical_qc.csv",
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    print(json.dumps(run(args.manifest,args.outdir),indent=2))

if __name__=="__main__":
    main()
