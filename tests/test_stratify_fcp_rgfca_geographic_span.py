"""Synthetic original dated metadata ONLY. No network, photo pixels, or colours."""
from __future__ import annotations

import csv
import json

import pytest

from fcp_pipeline.photo_polymorphism_screen import sha256
from scripts.acquisition.extract_rgfca_1000_original_monthly_photo_ids import CSV_FIELDS
from scripts.acquisition.stratify_fcp_rgfca_geographic_span import (
    ROSTER_FIELDS, SCHEMA, STATUS, _angle_km, _xyz, run,
)


def fixture(tmp_path):
    roster=tmp_path/"original_rgfca_selected_species.csv"
    raw=tmp_path/"rgfca_1000_complete_original_source_photo_IDs.csv"
    receipt=tmp_path/"rgfca_1000_targeted_monthly_photo_ID_receipt.json"
    taxa=[
        (101,"Widespread photographically","medium_256_1023"),
        (102,"Narrow photographically","deep_ge1024"),
        (103,"One faraway questionable point","shallow_128_255"),
        (104,"Across dateline but locally narrow","medium_256_1023"),
    ]
    photos={}
    for tid,name,tier in taxa:
        photos[tid]=[]
        if tid==101:
            coords=[(35.0+i*.1,139.0+i*.1) for i in range(4)]
            coords += [(51.0+i*.1,-0.5+i*.1) for i in range(4)]
            coords += [(37.0+i*.1,-122.0+i*.1) for i in range(4)]
        elif tid==102:
            coords=[(35+(i%4)*.04,139+(i%5)*.04) for i in range(20)]
        elif tid==103:
            coords=[(35+(i%4)*.02,139+(i%3)*.02) for i in range(11)]
            coords += [(-33.9,151.2)]
        else:
            coords=[(0+(i%4)*.03,179.92 if i<6 else -179.92) for i in range(12)]
        for j,(lat,lon) in enumerate(coords):
            obs=f"source_obs_{tid}_{j}"
            observer=f"observer_{tid}_{j%6 if tid!=101 else (j//4)*2+j%2}"
            photos[tid].append(dict(
                inat_taxon_id=tid,species=name,
                photo_id=str(tid*1000+j),observation_id=obs,
                observer_id=observer,
                latitude=lat,longitude=lon,observed_on="2026-08-20",
                positional_accuracy_m="250",
                photo_license="cc-by",photo_creator_id=observer,
                attribution="Verified CC-BY source creator",
                photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{tid*1000+j}/large.JPG",
                source_snapshot_date="2026-09-27"))
        if tid==101:
            repeated={**photos[tid][0],"photo_id":str(tid*1000+99)}
            repeated["photo_url_large"]=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{repeated['photo_id']}/large.JPG"
            photos[tid].append(repeated)
        if tid==102:
            # Inaccurate faraway observer position is never broad evidence.
            bad={**photos[tid][0],"photo_id":str(tid*1000+99),
                 "observation_id":"independent_with_low_accuracy",
                 "observer_id":"new_bad_observer",
                 "latitude":-20.0,"longitude":150.0,
                 "positional_accuracy_m":"100000"}
            bad["photo_url_large"]=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{bad['photo_id']}/large.JPG"
            photos[tid].append(bad)
    with roster.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=ROSTER_FIELDS,lineterminator="\n")
        w.writeheader()
        for tid,name,tier in taxa:
            w.writerow({"inat_taxon_id":tid,"species":name,"source_photo_opportunities":len(photos[tid]),
                        "opportunity_stratum":tier,"enrolment_hash":"source-only-placeholder"})
    with raw.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=CSV_FIELDS,lineterminator="\n")
        w.writeheader()
        for tid,_,_ in taxa:
            for x in photos[tid]:
                w.writerow(x)
    receipt.write_text(json.dumps({
        "schema":"fcp_rgfca_targeted_original_20260927_photo_ID_census_v1",
        "status":"COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS",
        "source_snapshot_date":"2026-09-27",
        "selected_all_photo_IDs_csv_sha256":sha256(raw),
        "frozen_1000_source_enrolment_sha256":sha256(roster),
        "all_enrolled_species_match_independent_full_source_photo_ID_counts":True,
        "original_image_pixels_opened":False,
        "historical_flower_colour_CSV_imported":False,
        "photos_NOT_classified_or_species_FCP_inferred":True,
        "n_selected_original_photo_IDs":sum(len(photos[tid]) for tid,_,_ in taxa),
    }),encoding="utf-8")
    return raw,receipt,roster


def results(out):
    with (out/"rgfca_1000_original_source_observed_geographic_span.csv").open(newline="") as f:
        return {int(row["inat_taxon_id"]):row for row in csv.DictReader(f)}


def test_true_observer_supported_broad_and_narrow_while_removing_single_point_outliers(tmp_path):
    raw,receipt,roster=fixture(tmp_path)
    out=tmp_path/"out"
    got=run(raw,receipt,roster,out,synthetic=True)
    assert got["schema"]==SCHEMA and got["status"]==STATUS
    assert got["source_species"]==4
    assert got["geographic_categories"]=={
        "broad_observed_1000km_supported":1,
        "narrow_observed_p95_le200km":2,
        "intermediate_or_unresolved_observed_span":1,
    }
    z=results(out)
    assert float(z[101]["max_supported_cluster_separation_km"])>1000
    assert int(z[101]["n_distinct_observations"])==12  # two photos same original observation
    assert int(z[101]["n_original_photo_metadata"])==13
    assert z[101]["geographic_sampling_category"]=="broad_observed_1000km_supported"
    assert z[102]["geographic_sampling_category"]=="narrow_observed_p95_le200km"
    assert int(z[102]["n_rejected_coordinate_accuracy"])==1
    assert z[103]["geographic_sampling_category"]=="intermediate_or_unresolved_observed_span"
    assert z[104]["geographic_sampling_category"]=="narrow_observed_p95_le200km"
    assert z[101]["is_native_range_confirmed"]=="False"
    assert all(x["original_image_pixels_opened"]=="False" for x in z.values())
    first=list(csv.DictReader((out/"rgfca_1000_precolour_broad_geography_first_priority.csv").open()))
    assert [int(first[i]["inat_taxon_id"]) for i in (0,1)]==[101,104]
    pairs=list(csv.DictReader((out/"rgfca_1000_original_source_broad_narrow_depthmatched_pairs.csv").open()))
    assert len(pairs)==1
    assert pairs[0]["broad_taxon_id"]=="101" and pairs[0]["narrow_taxon_id"]=="104"
    assert pairs[0]["opportunity_stratum"]=="medium_256_1023"
    assert float(pairs[0]["source_photo_count_ratio"])<=4
    assert pairs[0]["matched_before_any_original_photo_pixels"]=="True"
    assert got["n_broad_narrow_pairs_matched_on_original_photo_depth"]==1
    assert got["matching_does_not_replace_original_1000_species_cohort"] is True
    assert got["geographic_matching_pair_csv_sha256"]==sha256(
        out/"rgfca_1000_original_source_broad_narrow_depthmatched_pairs.csv")
    assert got["existing_200x250_repetition_schedule_unchanged"] is True
    assert got["n_original_new_pixels_measured"]==0


def test_dateline_metric_geodesic_does_not_invent_global_extent():
    assert _angle_km(_xyz(0,179.9),_xyz(0,-179.9))<30


def test_modified_source_photo_hash_fails_before_prioritization(tmp_path):
    raw,receipt,roster=fixture(tmp_path)
    with raw.open("a") as f:
        f.write("tampering\n")
    with pytest.raises(ValueError,match="Missing correct untouched"):
        run(raw,receipt,roster,tmp_path/"invalid",synthetic=True)


def test_missing_taxon_or_photos_fails_even_with_forged_hash_receipt(tmp_path):
    raw,receipt,roster=fixture(tmp_path)
    with raw.open(newline="") as f:
        rows=list(csv.DictReader(f))
    with raw.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows[:-1])
    d=json.loads(receipt.read_text())
    d["selected_all_photo_IDs_csv_sha256"]=sha256(raw)
    d["n_selected_original_photo_IDs"]-=1
    receipt.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="source photo ID count drifted"):
        run(raw,receipt,roster,tmp_path/"invalid",synthetic=True)


def test_no_destructive_rerun_or_image_data_inputs(tmp_path):
    raw,receipt,roster=fixture(tmp_path)
    out=tmp_path/"out"
    original=run(raw,receipt,roster,out,synthetic=True)
    with pytest.raises(ValueError,match="Refuse to overwrite"):
        run(raw,receipt,roster,out,synthetic=True)
    assert original["climate_soil_W_or_H2_estimated"] is False
    assert not set(["morph","pigment","flower_colour","measured_colour"]) & set(results(out)[101])



def test_geographic_narrow_controls_dont_overmatch_unrelated_source_count_tiers(tmp_path):
    raw,receipt,roster=fixture(tmp_path)
    with roster.open(newline="") as f:
        rows=list(csv.DictReader(f))
    # The only geographically narrow same-tier taxon is moved to the same
    # native-free test stratum as the second narrow species; no matching
    # across source-photo-depth strata is permitted.
    for row in rows:
        if row["inat_taxon_id"]=="104":
            row["opportunity_stratum"]="shallow_128_255"
    with roster.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=ROSTER_FIELDS)
        w.writeheader();w.writerows(rows)
    r=json.loads(receipt.read_text())
    r["frozen_1000_source_enrolment_sha256"]=sha256(roster)
    receipt.write_text(json.dumps(r))
    out=tmp_path/"unmatched"
    result=run(raw,receipt,roster,out,synthetic=True)
    assert result["n_broad_narrow_pairs_matched_on_original_photo_depth"]==0
    assert result["source_species"]==4
    assert list(csv.DictReader((out/"rgfca_1000_original_source_broad_narrow_depthmatched_pairs.csv").open()))==[]
