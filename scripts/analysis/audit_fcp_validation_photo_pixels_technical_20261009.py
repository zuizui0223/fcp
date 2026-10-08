#!/usr/bin/env python3
"""Open only frozen FCP Validation photo IDs for *technical* quality diagnostics.

This script never labels flower presence, species identity, pigments or colour
states. No classification or hypothesis test. 115 frozen photo IDs only;
no iNaturalist search/API queries; no pixels committed or uploaded.
Images are held temporarily in memory and discarded after metrics.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from PIL import Image, ImageFilter, ImageOps, UnidentifiedImageError

EXPECTED_CSV_GIT_BLOB = "467f1a1f51a8e14f243d25e7f62550521bee0b40"
N_PHOTOS = 115
N_SPECIES = 25
VALID_LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}
EXTENSIONS = ("jpg", "jpeg")  # The two frozen, common iNaturalist JPEG extensions.
PHOTO_HOST = "https://inaturalist-open-data.s3.amazonaws.com/photos"
USER_AGENT = "fcp-frozen-validation-photo-technical-quality/1.0 (115-fixed-ids)"
MAX_BYTES = 6_000_000
MIN_INTERVAL_SEC = 0.4
DOWNLOAD_TIMEOUT_SEC = 10


def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\\0" + raw).hexdigest()


def read_frozen_candidates(path: Path) -> list[dict]:
    raw = path.read_bytes()
    if git_blob_sha(raw) != EXPECTED_CSV_GIT_BLOB:
        raise ValueError("Frozen 115-photo Git blob SHA mismatch")
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    if len(rows) != N_PHOTOS or len({r["photo_id"] for r in rows}) != N_PHOTOS:
        raise ValueError("115 unique frozen photo IDs were not recovered")
    if len({r["observation_id"] for r in rows}) != N_PHOTOS:
        raise ValueError("Unexpected repeated observation IDs")
    if len({r["inat_taxon_id"] for r in rows}) != N_SPECIES:
        raise ValueError("Frozen 25 species changed")
    if Counter(int(r["gap_class"]) for r in rows) != {1:39,2:76}:
        raise ValueError("Frozen 39/76 photo split changed")
    expected_review_keys = (
        "flower_visible", "species_identity_confirmed", "target_flower_organ",
        "exposure_acceptable", "colour_measurable", "individual_identity_verified",
    )
    for r in rows:
        if not all(str(r[k]).strip() == "" for k in expected_review_keys):
            raise ValueError("Botanical/phenotypic outcome leaked into frozen input")
        if r["quality_review_status"] != "UNREVIEWED":
            raise ValueError("Previously reviewed photo cannot enter blinded technical audit")
        if r["photo_license"] not in VALID_LICENSES:
            raise ValueError("Non-open source photo license")
        if not all(r[key].isdecimal() for key in ("photo_id", "observation_id", "inat_taxon_id", "observer_id")):
            raise ValueError("Non-numeric photo/observation/taxon/observer identity")
    return rows


def make_urls(photo_id: str) -> tuple[str, str]:
    if not photo_id.isdecimal() or not photo_id:
        raise ValueError("Unsafe photo ID")
    return tuple(f"{PHOTO_HOST}/{photo_id}/medium.{ext}" for ext in EXTENSIONS)


def retrieve(photo_id: str, *, fetcher=None) -> tuple[bytes|None, str|None, str]:
    if fetcher is None:
        fetcher = download
    failure = []
    for url in make_urls(photo_id):
        try:
            raw = fetcher(url)
            if not (1024 <= len(raw) <= MAX_BYTES):
                raise ValueError("Downloaded image outside accepted byte range")
            return raw, url, "RETRIEVED"
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,
                ValueError,OSError) as exc:
            failure.append(type(exc).__name__)
    return None, None, "UNAVAILABLE_" + "_".join(failure)


def download(url: str) -> bytes:
    if not any(url.startswith(f"{PHOTO_HOST}/") for _ in (0,)):
        raise ValueError("Unexpected photo download host")
    req = urllib.request.Request(url, headers={
        "User-Agent": USER_AGENT, "Accept": "image/jpeg",
    })
    with urllib.request.urlopen(req, timeout=DOWNLOAD_TIMEOUT_SEC) as resp:
        typ = resp.headers.get("Content-Type", "")
        if not typ.startswith("image/"):
            raise ValueError("Not an image MIME type")
        size = resp.headers.get("Content-Length")
        if size is not None and int(size) > MAX_BYTES:
            raise ValueError("Oversize image response")
        raw = resp.read(MAX_BYTES+1)
    if len(raw)>MAX_BYTES:
        raise ValueError("Oversize image body")
    return raw


def technical_metrics(raw: bytes) -> dict:
    with Image.open(io.BytesIO(raw)) as img:
        if img.format not in ("JPEG","PNG","WEBP"):
            raise ValueError(f"Unexpected image container: {img.format}")
        img = ImageOps.exif_transpose(img)
        img.load()
        im = img.convert("RGB")
        w,h = im.size
        if min(w,h)<1 or max(w,h)>4096:
            raise ValueError("Image dimensions invalid")
        im.thumbnail((320,320), Image.Resampling.LANCZOS)
        x=np.asarray(im,dtype=np.uint8)
    lum=(0.2126*x[...,0].astype(float)+0.7152*x[...,1].astype(float)+0.0722*x[...,2].astype(float))
    clip=float(np.mean(np.max(x,axis=-1)>=250))
    near_white=float(np.mean((x>=245).all(axis=-1)))
    dark=float(np.mean(np.max(x,axis=-1)<=12))
    edge=np.asarray(Image.fromarray(x).convert("L").filter(ImageFilter.FIND_EDGES),dtype=float)
    return {
        "width":int(w),"height":int(h),
        "mean_image_luma_0_255":round(float(lum.mean()),4),
        "near_highlight_pixel_fraction":round(clip,6),
        "near_white_RGB_pixel_fraction":round(near_white,6),
        "very_dark_pixel_fraction":round(dark,6),
        "edge_filter_variance_diagnostic":round(float(edge.var()),4),
        "image_sha256":hashlib.sha256(raw).hexdigest(),
    }


def analyse(rows: list[dict], *, fetcher=None, sleeper=time.sleep, interval=MIN_INTERVAL_SEC):
    if len(rows)!=N_PHOTOS:
        raise ValueError("Incomplete frozen photo set")
    results=[]
    for index,row in enumerate(rows):
        pid=row["photo_id"]
        payload,url,status=retrieve(pid,fetcher=fetcher)
        item={
            "inat_taxon_id":row["inat_taxon_id"],
            "species":row["species"],
            "gap_class":int(row["gap_class"]),
            "target_year":int(row["target_year"]),
            "calendar_month":int(row["calendar_month"]),
            "observation_id":row["observation_id"],
            "photo_id":pid,
            "observer_id":row["observer_id"],
            "photo_license":row["photo_license"],
            "photo_page_url":f"https://www.inaturalist.org/photos/{pid}",
            "observation_page_url":f"https://www.inaturalist.org/observations/{row['observation_id']}",
            "source_image_url":url,
            "image_technical_status":status,
            "review_status":"UNREVIEWED",
            "flower_visible":"UNKNOWN",
            "species_identity_confirmed":"UNKNOWN",
            "target_flower_organ":"UNKNOWN",
            "exposure_acceptable":"UNKNOWN",
            "colour_measurable":"UNKNOWN",
            "individual_identity_verified":"UNKNOWN",
        }
        if payload is not None:
            try:
                item.update(technical_metrics(payload))
                item["image_technical_status"]="DECODED"
            except (UnidentifiedImageError,ValueError,OSError) as exc:
                item["image_technical_status"]="RETRIEVED_BUT_DECODE_FAILED"
                item["technical_error_class"]=type(exc).__name__
        results.append(item)
        if index+1<len(rows):
            sleeper(interval)
    if len(results)!=N_PHOTOS or len({x["photo_id"] for x in results})!=N_PHOTOS:
        raise RuntimeError("Any missing or duplicated original photo identity")
    sha=[r["image_sha256"] for r in results if r.get("image_sha256")]
    duplicate_hashes={h:n for h,n in Counter(sha).items() if n>1}
    for r in results:
        r["same_image_bytes_another_photo_id"]=r.get("image_sha256") in duplicate_hashes
    decoded=sum(r["image_technical_status"]=="DECODED" for r in results)
    report={
        "schema":"fcp_validation_115_frozen_photo_pixel_technical_v1",
        "status":"TECHNICAL_ONLY_IMAGE_BYTES_OPENED_NO_BOTANICAL_REVIEW",
        "source_csv_git_blob_sha":EXPECTED_CSV_GIT_BLOB,
        "n_photo_ids":N_PHOTOS,"n_species":N_SPECIES,
        "n_image_decode_success":decoded,
        "n_image_decode_or_transport_unavailable":N_PHOTOS-decoded,
        "n_duplicate_image_content_hashes":len(duplicate_hashes),
        "n_botanical_images_confirmed":None,
        "n_species_photo_quality_validated":None,
        "n_colour_labels_assigned":0,
        "original_10km_gate":"HOLD_UNCHANGED",
        "image_bytes_uploaded_or_committed":False,
        "technical_flags_are_never_automatic_botanical_fail_or_pass":True,
        "nonclaims":[
            "Brightness, edge energy and clipping are whole-image diagnostics not flower-petal ROI measurements",
            "Decodable photo can show leaves, soil, insects, garden signs, nonfocal flowers, or wrong species",
            "Image is never scored as botanically usable without explicit organ and identity review",
            "No published colour-frequency, climate, genotype, evolutionary or fitness result follows from these metrics",
            "Identities, photos and hypothesis gates are never replaced or changed based on technical values",
        ],
    }
    return report,results


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--queue",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    records=read_frozen_candidates(args.queue)
    report,items=analyse(records)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\\n",encoding="utf-8")
    (args.outdir/"photo_technical_diagnostics.json").write_text(json.dumps(items,indent=2,sort_keys=True)+"\\n",encoding="utf-8")
    with (args.outdir/"photo_technical_diagnostics.csv").open("w",newline="",encoding="utf-8") as f:
        names=sorted({k for row in items for k in row})
        w=csv.DictWriter(f,fieldnames=names)
        w.writeheader();w.writerows(items)
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
