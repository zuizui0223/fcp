#!/usr/bin/env python3
"""Produce an INTERNAL original-byte-verified visual review pilot, labels sealed.

Strict research-only temporary artifact: do not commit, publicly distribute, or
publish image pixels/contact sheets. Reviewers must follow original authors'
licences and attribution before any external presentation or distribution.
Neither species nor source machine colour/panel membership enters photo review.
"""
from __future__ import annotations
import argparse
import hashlib
from io import BytesIO
import json
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import time
import pandas as pd
from PIL import Image,ImageDraw,ImageFont,ImageOps

EXPECTED_CSV_SHA="06191b36eca9862f6241b5110ae13aee0ff7e7414a5caaaa1025daffca33eadd"
HOST="https://inaturalist-open-data.s3.amazonaws.com"
SIZES=("large","original")
EXTENSIONS=("jpeg","jpg","png")
MAX_BYTES=8*1024*1024
MAX_IMAGE_PIXELS=40_000_000
USER_AGENT="fcp-original-visual-pilot/1.0 (github.com/zuizui0223/fcp; internal original-photo review)"
Image.MAX_IMAGE_PIXELS=MAX_IMAGE_PIXELS

def read_blind(path):
    payload=Path(path).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=EXPECTED_CSV_SHA:
        raise RuntimeError("Original 24-image blinded pilot manifest SHA mismatch")
    d=pd.read_csv(BytesIO(payload),dtype=str)
    if list(d.columns)!=["audit_case_id","photo_id","frozen_original_image_sha256"]:
        raise RuntimeError("Only the original case/photo/hash blind columns may be loaded")
    if len(d)!=24 or d.photo_id.duplicated().any() or d.audit_case_id.duplicated().any():
        raise RuntimeError("Original 24 unique case IDs or photo IDs drift")
    if not d.photo_id.str.fullmatch(r"[1-9][0-9]*").all():
        raise RuntimeError("Unsafe source photo ID")
    if not d.frozen_original_image_sha256.str.fullmatch("[0-9a-f]{64}").all():
        raise RuntimeError("Source image SHA unavailable")
    return d.sort_values("audit_case_id",kind="stable").reset_index(drop=True)

def get_one(photo_id,original_sha,*,opener=urlopen,pause=time.sleep):
    for size in SIZES:
        for ext in EXTENSIONS:
            url=f"{HOST}/photos/{int(photo_id)}/{size}.{ext}"
            try:
                req=Request(url,headers={"User-Agent":USER_AGENT,"Accept":"image/jpeg,image/png"})
                with opener(req,timeout=20) as response:
                    final=response.geturl()
                    if not final.startswith(f"{HOST}/photos/{int(photo_id)}/"):
                        raise RuntimeError("Untrusted photo-host redirect")
                    payload=response.read(MAX_BYTES+1)
                if not (1024<=len(payload)<=MAX_BYTES):raise RuntimeError("Invalid photo bytes")
                if hashlib.sha256(payload).hexdigest()!=original_sha:
                    continue
                with Image.open(BytesIO(payload)) as original:
                    image=ImageOps.exif_transpose(original).convert("RGB")
                    image.load()
                    if image.width*image.height>MAX_IMAGE_PIXELS: raise RuntimeError("Oversized pixels")
                    return image, url
            except (HTTPError,URLError,RuntimeError,OSError,ValueError):
                pass
            pause(.08)
    return None,""

def render_sheet(entries,out):
    cols=4;rows=3;w,h=430,360
    sheet=Image.new("RGB",(cols*w,rows*h), "#ffffff")
    draw=ImageDraw.Draw(sheet)
    for ix,(case_id,im) in enumerate(entries):
        c,r=ix%cols,ix//cols
        x,y=c*w,r*h
        z=ImageOps.contain(im,(w-22,h-64),Image.Resampling.LANCZOS)
        sheet.paste(z,(x+(w-z.width)//2,y+(h-64-z.height)//2))
        draw.rectangle((x,y+h-44,x+w-1,y+h-1),fill="#f4f4f4",outline="#999999")
        draw.text((x+10,y+h-31),case_id,fill="#111111")
    out.parent.mkdir(parents=True,exist_ok=True)
    sheet.save(out,format="PNG",optimize=True)

def run(manifest,outdir,*,opener=urlopen,pause=time.sleep):
    d=read_blind(manifest)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    rows=[];picks=[]
    for row in d.itertuples(index=False):
        image,url=get_one(row.photo_id,row.frozen_original_image_sha256,opener=opener,pause=pause)
        status="source_byte_verified_visual_pilot" if image is not None else "unavailable_or_source_byte_mismatch"
        rows.append({"audit_case_id":row.audit_case_id,"status":status,
                     "original_photo_byte_verified":image is not None})
        if image is not None:
            picks.append((row.audit_case_id,image))
        print(f"blind_case={row.audit_case_id} status={status}",flush=True)
    manifest_out=pd.DataFrame(rows).sort_values("audit_case_id",kind="stable")
    manifest_out.to_csv(out/"blinded_pilot_visual_recovery_status.csv",index=False)
    for ix,start in enumerate(range(0,len(picks),12),1):
        render_sheet(picks[start:start+12],out/f"REVIEW_INTERNAL_blind_source_photos_{ix:02d}.png")
    result={
        "schema":"fcp_original24_internal_blind_visual_audit_v1",
        "n_expected_photo_cases":24,
        "n_original_source_byte_identical_visual_cases":len(picks),
        "n_failed_photo_cases":24-len(picks),
        "n_contact_sheets":(len(picks)+11)//12,
        "all_original_machine_colour_labels_and_species_sealed":True,
        "source_pilot_sha256":EXPECTED_CSV_SHA,
        "human_expert_reviews_completed":0,
        "automatic_colour_classifications_completed":0,
        "status":"INTERNAL_REVIEW_IMAGES_ONLY",
        "license_notice":"Images originate from iNaturalist openly licensed photo bucket. Do not publish/re-distribute pixels without original photograph creator attribution, confirmed exact licence and permitted use. Do not commit these pixels.",
        "hard_nonclaims":["Photo bytes matching prior source do not establish botanical organ or colour label correctness.",
                          "The original expert review still requires two independent blinded humans.",
                          "This is a purposefully small 24-image selected pilot, not a 405-photo measured confusion matrix."]
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.manifest,a.outdir),indent=2))

if __name__=="__main__": main()
