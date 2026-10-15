"""Source-only end-to-end no-loss guard between audited count and full ID freeze."""
import csv
import hashlib
import json

import pytest

from scripts.acquisition.verify_fcp_complete_photo_ids_vs_original_counts import (
    ARCHIVE_SHA, FRAME_SHA, SOURCE_DATE, validate,
)


def fixture(tmp_path):
    baseline=tmp_path/"original_count.csv"
    full_counts=tmp_path/"full_photo_ID_manifest_counts.csv"
    with baseline.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"])
        w.writerows([[11,"Species alba",2,SOURCE_DATE],[12,"Species rubra",0,SOURCE_DATE],
                     [13,"Species lutea",1,SOURCE_DATE],[14,"Species flava",5,SOURCE_DATE]])
    with full_counts.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_target_photo_ids","source_status","source_manifest_sha256"])
        w.writerows([[11,"Species alba",2,"complete","a"*64],
                     [12,"Species rubra",0,"complete","b"*64],
                     [13,"Species lutea",1,"complete","c"*64],
                     [14,"Species flava",5,"complete","d"*64]])
    audit={
        "schema":"fcp_42111_original_photo_opportunity_audit_v1",
        "status":"VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT",
        "source_frame_species":4,
        "source_frame_sha256":FRAME_SHA,
        "source_archive_sha256":ARCHIVE_SHA,
        "snapshot_date":SOURCE_DATE,
        "count_ledger_sha256":hashlib.sha256(baseline.read_bytes()).hexdigest(),
        "original_source_taxon_IDs_and_names_independently_matched":True,
        "all_source_photo_id_manifests_frozen":False,
        "all_original_image_pixels_measured":False,
        "W_D_H2_or_climate_soil_results":False,
        "all_source_photo_opportunity_ids":8,
    }
    source={
        "status":"COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE",
        "complete_metadata_census":True,
        "target_frame_species":4,
        "species_count_checkpoints_audited":4,
        "original_archive_sha256":ARCHIVE_SHA,
        "immutable_source_date":SOURCE_DATE,
        "full_pixel_measurement_complete":False,
        "biological_h2_and_environment_inference_ready":False,
        "old_measured_flower_colours_used":False,
        "no_per_species_photo_cap":True,
        "frozen_location_blind_worker_partitions":128,
        "n_all_available_licensed_original_photo_ids":8,
    }
    frozen={
        "status":"FRESH_BLIND_METADATA_FROZEN",
        "source_species":4,
        "worker_packets":128,
        "historical_measured_CSV_read":False,
        "all_available_licensed_photo_ids":8,
    }
    docs=[]
    for name,data in [("audit.json",audit),("source.json",source),("firewall.json",frozen)]:
        path=tmp_path/name
        path.write_text(json.dumps(data),encoding="utf-8")
        docs.append(path)
    return (baseline,*docs[:1],full_counts,*docs[1:])


def test_complete_42111_frame_count_comparison_synthetic_no_photo_pixels(tmp_path):
    data=fixture(tmp_path)
    result=validate(*data,_synthetic_test=True)
    assert result["status"]=="REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS"
    assert result["species"]==4
    assert result["distinct_original_photo_ids"]==8
    assert result["original_photo_IDs_covered_without_count_loss"] is True
    assert result["full_study_original_pixels_classified"] is False


def test_drops_one_source_photo_id_must_fail(tmp_path):
    data=fixture(tmp_path)
    path=data[2]
    path.write_text(path.read_text().replace("Species flava,5,complete","Species flava,4,complete"))
    with pytest.raises(ValueError, match="differs"):
        validate(*data,_synthetic_test=True)


def test_replaces_original_taxon_identity_must_fail(tmp_path):
    data=fixture(tmp_path)
    path=data[2]
    path.write_text(path.read_text().replace("Species flava","Fake replacement species"))
    with pytest.raises(ValueError, match="differs"):
        validate(*data,_synthetic_test=True)


def test_will_not_accept_metadata_only_source_receipt(tmp_path):
    data=fixture(tmp_path)
    path=data[3]
    source=json.loads(path.read_text())
    source["complete_metadata_census"]=False
    path.write_text(json.dumps(source))
    with pytest.raises(ValueError, match="NOT a genuinely complete"):
        validate(*data,_synthetic_test=True)


def test_will_not_accept_fake_pixel_outcomes_or_unsealed_manifest(tmp_path):
    data=fixture(tmp_path)
    path=data[4]
    frozen=json.loads(path.read_text())
    frozen["status"]="DRAFT_NOT_FROZEN"
    path.write_text(json.dumps(frozen))
    with pytest.raises(ValueError, match="NOT a genuinely complete"):
        validate(*data,_synthetic_test=True)
