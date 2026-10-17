"""Second disjoint RGFCA six-taxon technical pilot: synthetic no-network gates."""
from __future__ import annotations
import csv
import json
import pytest
from scripts.acquisition.run_rgfca_real_1000_enrolled_six_photo_pilot import (
    pilot_taxa as first_taxa, TIERS, SAMPLE_FIELDS,
)
from scripts.acquisition.run_rgfca_disjoint_wave2_six_taxa_original_photo_pilot import (
    pilot_taxa as second_taxa, precolour_metadata, measure, sampling_hash,
)


def roster(tmp_path):
    path=tmp_path/"rgfca_first1000.csv"
    with path.open("w",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(["inat_taxon_id","species","source_photo_opportunities",
                     "opportunity_stratum","enrolment_hash"])
        for i,tier in enumerate(TIERS):
            for k in range(6):
                tax=1000+i*50+k
                wr.writerow([tax,f"Species {tax}",128+k,tier,"fake_source_no_colour"])
    return path


def fake_obs(oid,tid,user,pid):
    return {
        "id":oid,"taxon":{"id":tid,"ancestor_ids":[]},"user":{"id":user},
        "geojson":{"coordinates":[130.0,35.0+(oid%100)*.001]},
        "observed_on":"2026-08-01",
        "photos":[{"id":pid,"license_code":"cc-by",
                   "url":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/medium.JPG",
                   "attribution":"CC-BY original biological photograph"}],
    }


def test_second_six_are_deterministic_and_never_overlap_original_six(tmp_path):
    frame=roster(tmp_path)
    prior=first_taxa(frame,_synthetic=True)
    second=second_taxa(frame,_synthetic=True)
    assert second==second_taxa(frame,_synthetic=True)
    assert len(prior)==len(second)==6
    assert not ({x["inat_taxon_id"] for x in prior}&
                {x["inat_taxon_id"] for x in second})
    assert all(sum(x["source_opportunity_tier"]==tier for x in second)==2
               for tier in TIERS)
    assert all(x["original_20260927_source_photo_count"]>=128 for x in second)


def test_second_source_selection_keeps_precolour_fixed_observer_disjoint_folds(tmp_path):
    frame=roster(tmp_path)
    requested=[]
    def provider(params):
        requested.append(params.copy())
        tax=int(params["taxon_id"])
        return {"results":[fake_obs(tax*100+j,tax,tax+j,tax*10000+j)
                           for j in range(12)]}
    source=tmp_path/"new_disjoint_six"
    receipt=precolour_metadata(frame,source,request=provider,_synthetic=True)
    assert receipt["wave"]==2
    assert receipt["technical_only_disjoint_previous_six_taxa"] is True
    assert receipt["pilot_taxa"]==6
    assert len(requested)==6
    assert all(p["created_d2"]=="2026-09-27" and
               p["quality_grade"]=="research" for p in requested)
    assert receipt["NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS"] is True
    assert receipt["EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED"] is False
    assert receipt["image_bytes_opened"] is False
    assert receipt["RGFCA_200_realization_photo_polymorphism_Ppoly_computed"] is False
    photos=source/"rgfca_six_wave2_source_photo_candidates_before_colour.csv"
    with photos.open(newline="") as f:
        reader=csv.DictReader(f)
        assert tuple(reader.fieldnames)==SAMPLE_FIELDS
        saved=list(reader)
    assert len(saved)==72
    assert len({r["photo_id"] for r in saved})==72
    assert all(r["sample_hash"]==sampling_hash(
        "original_photo",r["inat_taxon_id"],r["photo_id"]) for r in saved)


def test_second_wave_unverified_source_fails_before_original_pixel_inference(tmp_path):
    root=tmp_path/"fake"
    root.mkdir()
    (root/"rgfca_six_wave2_source_photo_candidates_before_colour.csv").write_text(
        ",".join(SAMPLE_FIELDS)+"\n")
    (root/"rgfca_six_wave2_source_precolour_receipt.json").write_text(json.dumps({
        "status":"UNVERIFIED_FCP_PRECOLOUR_SOURCE",
        "selected_metadata_sha256":"not_a_digest",
        "EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED":True,
    }))
    with pytest.raises(ValueError,match="Unverified"):
        measure(root,tmp_path/"should_not_exist",tmp_path/"model",
                tmp_path/"detector",tmp_path/"sam")
