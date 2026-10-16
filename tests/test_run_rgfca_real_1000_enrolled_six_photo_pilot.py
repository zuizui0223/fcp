"""Mock provider only. Pilot source selection cannot read colour or become RGFCA result."""
from __future__ import annotations

import csv
import json

import pytest

from scripts.acquisition.run_rgfca_real_1000_enrolled_six_photo_pilot import (
    FINAL_FIELDS, SAMPLE_FIELDS, PILOT_SALT, TIERS,
    _trusted_s3_url, measure, pilot_taxa, precolour_metadata,
    sampling_hash, select_photo_rows,
)


def roster(tmp_path):
    path = tmp_path/"rgfca_selected_species.csv"
    with path.open("w",newline="") as h:
        writer=csv.writer(h)
        writer.writerow(["inat_taxon_id","species","source_photo_opportunities",
                         "opportunity_stratum","enrolment_hash"])
        for i,tier in enumerate(TIERS):
            for k in range(3):
                tid=1000+i*10+k
                writer.writerow([tid,f"Species {tid}",128,tier,"real_from_source_hash"])
    return path


def observation(oid: int, tid: int, user: int, pid: int):
    return {
        "id": oid,
        "taxon": {"id":tid,"ancestor_ids":[]},
        "user":{"id":user},
        "geojson":{"coordinates":[130.0+oid*.001,35.0]},
        "observed_on":"2026-08-01",
        "photos":[{"id":pid,"license_code":"cc-by",
                   "url":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/medium.JPG",
                   "attribution":"CC-BY biological photograph"}],
    }


def test_original_rgfca_taxon_enrolment_is_not_based_on_colour(tmp_path):
    r=roster(tmp_path)
    a=pilot_taxa(r,_synthetic=True)
    b=pilot_taxa(r,_synthetic=True)
    assert a==b
    assert len(a)==6
    assert all(sum(t["source_opportunity_tier"]==tier for t in a)==2 for tier in TIERS)
    assert sampling_hash("taxon",23)==sampling_hash("taxon",23)
    assert len({x["inat_taxon_id"] for x in a})==6
    assert all(x["original_20260927_source_photo_count"]==128 for x in a)


def test_provider_url_cannot_swap_untrusted_photo_ID_or_host():
    assert _trusted_s3_url("111","https://inaturalist-open-data.s3.amazonaws.com/photos/111/medium.JPG")==(
        "https://inaturalist-open-data.s3.amazonaws.com/photos/111/large.JPG")
    assert _trusted_s3_url("111","https://static.inaturalist.org/photos/111/square.jpg")==(
        "https://inaturalist-open-data.s3.amazonaws.com/photos/111/large.jpg")
    assert _trusted_s3_url("112","https://static.inaturalist.org/photos/111/square.jpg") is None
    assert _trusted_s3_url("111","https://evil.example/photos/111/medium.jpg") is None
    assert _trusted_s3_url("111","http://static.inaturalist.org/photos/111/medium.jpg") is None


def test_first_original_photo_source_is_a_bounded_live_API_not_archive(tmp_path):
    r=roster(tmp_path)
    count=[]
    def simulated_request(params):
        count.append(dict(params))
        ident=int(params["taxon_id"])
        return {"results":[observation(ident*100+j,ident,ident+j,ident*10000+j) for j in range(12)]}
    out=tmp_path/"technical-source"
    report=precolour_metadata(r,out,request=simulated_request,_synthetic=True)
    assert report["pilot_taxa"]==6
    assert len(count)==6  # An incomplete (<200) page halts each taxon.
    assert all(z["created_d2"]=="2026-09-27" for z in count)
    assert all(z["quality_grade"]=="research" for z in count)
    assert report["EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED"] is False
    assert report["NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS"] is True
    assert report["RGFCA_200_realization_photo_polymorphism_Ppoly_computed"] is False
    assert report["image_bytes_opened"] is False
    with (out/"rgfca_six_source_photo_candidates_before_colour.csv").open(newline="") as h:
        reader=csv.DictReader(h)
        assert tuple(reader.fieldnames)==SAMPLE_FIELDS
        photos=list(reader)
    assert len(photos)==72
    assert len({p["photo_id"] for p in photos})==72
    assert all(p["sample_hash"]==sampling_hash(
        "original_photo",p["inat_taxon_id"],p["photo_id"]) for p in photos)
    assert all("large.JPG" in p["photo_url_large"] for p in photos)


def test_photo_colour_or_prior_results_never_select_api_ids(tmp_path):
    taxon={"inat_taxon_id":1000,"species":"Species 1000",
           "source_opportunity_tier":TIERS[0]}
    entries=[{"inat_taxon_id":1000,"species":"Species 1000",
              "photo_id":str(1000+i),"observation_id":str(i),
              "observer_id":str(i),
              "latitude":35.0,"longitude":135.0,
              "photo_url_large":f"https://static.inaturalist.org/photos/{1000+i}/medium.jpg",
              "photo_license":"cc-by"} for i in range(100)]
    a=select_photo_rows(taxon,entries,per_fold=4)
    outcome=select_photo_rows(taxon,[{**x,"morph":"red_pink","soil_pH":9.9,
                                    "old_FCP_outcome":True} for x in entries],per_fold=4)
    assert a==outcome
    assert len(a)==8
    assert not ({x["photo_id"] for x in a}&{"-1000"})
    assert all(x["photo_url_large"].startswith(
        "https://inaturalist-open-data.s3.amazonaws.com/photos/") for x in a)


def test_provider_cursor_and_duplicate_observations_cannot_fake_independence(tmp_path):
    r=roster(tmp_path)
    def invalid_cursor(params):
        ident=int(params["taxon_id"])
        return {"results":[observation(51,ident,1,1234),
                           observation(50,ident,2,1235)]}
    with pytest.raises(ValueError,match="Nonmonotonic"):
        precolour_metadata(r,tmp_path/"bad-cursor",
                           request=invalid_cursor,_synthetic=True)


def test_real_pixel_worker_requires_precolour_receipt_before_method(tmp_path):
    selection=tmp_path/"source"
    selection.mkdir()
    meta=selection/"rgfca_six_source_photo_candidates_before_colour.csv"
    meta.write_text(",".join(SAMPLE_FIELDS)+"\n")
    (selection/"rgfca_six_source_precolour_receipt.json").write_text(json.dumps({
        "status":"Fake source selection",
        "selected_metadata_sha256":"WRONG",
        "EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED":True
    }))
    with pytest.raises(ValueError,match="Unverified"):
        measure(selection,tmp_path/"out",tmp_path/"missing-source",
                tmp_path/"missing-detector",tmp_path/"missing-sam")
