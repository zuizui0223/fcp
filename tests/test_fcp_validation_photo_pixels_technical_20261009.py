"""Technical image preflight only; no new observations or biological outcomes."""
from __future__ import annotations
import hashlib
import io
import sys
import urllib.error
from pathlib import Path

import numpy as np
from PIL import Image
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from audit_fcp_validation_photo_pixels_technical_20261009 import (
    read_frozen_candidates,make_urls,retrieve,technical_metrics,analyse,git_blob_sha,
    N_PHOTOS, PHOTO_HOST,
)

ROOT=Path(__file__).resolve().parents[1]
QUEUE=ROOT/"results"/"fcp_validation_photo_quality_queue_20261009"/"candidate_photo_quality_review.csv"


def synthetic_image():
    array=np.random.default_rng(17).integers(0,255,size=(260,310,3),dtype=np.uint8)
    mem=io.BytesIO()
    Image.fromarray(array,"RGB").save(mem,format="JPEG",quality=80)
    return mem.getvalue()


def test_immutable_source_and_no_preexisting_qa():
    rows=read_frozen_candidates(QUEUE)
    assert len(rows)==115
    assert len({r["inat_taxon_id"] for r in rows})==25
    assert all(r["quality_review_status"]=="UNREVIEWED" for r in rows)
    assert len({r["photo_id"] for r in rows})==115
    assert git_blob_sha(QUEUE.read_bytes())=="467f1a1f51a8e14f243d25e7f62550521bee0b40"


def test_photo_url_is_exact_allowlisted_host_numeric_id():
    u=make_urls("123456789")
    assert u==(f"{PHOTO_HOST}/123456789/medium.jpg",f"{PHOTO_HOST}/123456789/medium.jpeg")
    with pytest.raises(ValueError,match="Unsafe"):
        make_urls("../etc/passwd")


def test_decoder_metrics_are_not_organ_classification():
    d=technical_metrics(synthetic_image())
    assert d["width"]==310 and d["height"]==260
    assert 0<=d["near_highlight_pixel_fraction"]<=1
    assert d["image_sha256"]==hashlib.sha256(synthetic_image()).hexdigest()
    assert "flower_visible" not in d
    assert "colour_measurable" not in d


def test_jpg_absent_jpeg_fallback_and_zero_image_replacement():
    valid=synthetic_image()
    seen=[]
    def fetcher(url):
        seen.append(url)
        if url.endswith(".jpg"):
            raise urllib.error.HTTPError(url,404,"not_found",{},None)
        return valid
    raw,url,status=retrieve("1234",fetcher=fetcher)
    assert raw==valid and url.endswith(".jpeg") and status=="RETRIEVED"
    assert len(seen)==2


def test_115_fixed_images_keep_all_botanical_statuses_unknown():
    rows=read_frozen_candidates(QUEUE)
    valid=synthetic_image()
    r,images=analyse(rows,fetcher=lambda url:valid,sleeper=lambda t:None,interval=0)
    assert r["n_image_decode_success"]==115
    assert r["n_botanical_images_confirmed"] is None
    assert r["original_10km_gate"]=="HOLD_UNCHANGED"
    assert len(images)==N_PHOTOS
    assert len({x["photo_id"] for x in images})==115
    assert all(x["colour_measurable"]=="UNKNOWN" for x in images)
    assert all(x["review_status"]=="UNREVIEWED" for x in images)
    assert r["n_duplicate_image_content_hashes"]==1


def test_unreadable_photo_is_unresolved_not_monochromatic_or_failed_morph():
    rows=read_frozen_candidates(QUEUE)
    fake=synthetic_image()
    calls=0
    def fetcher(url):
        nonlocal calls
        calls+=1
        if "/"+rows[0]["photo_id"]+"/" in url:
            raise urllib.error.URLError("simulated unavailable")
        return fake
    r, images=analyse(rows,fetcher=fetcher,sleeper=lambda t:None,interval=0)
    assert r["n_image_decode_success"]==114
    assert r["n_image_decode_or_transport_unavailable"]==1
    assert images[0]["image_technical_status"].startswith("UNAVAILABLE")
    assert images[0]["flower_visible"]=="UNKNOWN"
    assert images[0]["colour_measurable"]=="UNKNOWN"
    assert calls==116
