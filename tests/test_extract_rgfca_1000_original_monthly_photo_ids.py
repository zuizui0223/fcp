"""Synthetic official four-gzip monthly snapshot, RGFCA-selected 1000 photo-ID audits.

The new smaller SQLite store MUST be byte-pinned, count-equivalent to the
original 42,111-source importer, globally collision-aware and outcome-blind.
No iNaturalist API, old measured colour CSV or original flower pixel is used.
"""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pytest

from fcp_pipeline.fcp_allphotos import sha256_file
from scripts.acquisition.extract_rgfca_1000_original_monthly_photo_ids import (
    CSV_FIELDS, process,
)


def gz_file(path, head, rows):
    with gzip.open(path,"wt",encoding="utf-8",newline="") as f:
        writer=csv.writer(f,delimiter="\t",lineterminator="\n")
        writer.writerow(head)
        writer.writerows(rows)
    return path


def setup(tmp_path, *, invalid_external_collision=False):
    frame=tmp_path/"original_42111_frame.csv"
    with frame.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","after_observer_cap"])
        w.writerows([[11,"Species alba",500],[12,"Species rubra",300],
                     [13,"Species lutea",250],[14,"Species cyanea",180]])
    # The FULL global photo count from the original full-41k importer:
    # 200 is reused cross-taxon and therefore globally quarantined. 202 is
    # duplicated IDENTICALLY, allowed once as original source ID.
    counts=tmp_path/"all_42111_species_photo_opportunity_counts.csv"
    with counts.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"])
        w.writerows([[11,"Species alba",3,"2026-09-27"],
                     [12,"Species rubra",2,"2026-09-27"],
                     [13,"Species lutea",1,"2026-09-27"],
                     [14,"Species cyanea",1,"2026-09-27"]])
    independent=tmp_path/"independent_audit.json"
    independent.write_text(json.dumps({
        "schema":"fcp_42111_original_photo_opportunity_audit_v1",
        "status":"VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT",
        "source_frame_species":4,
        "all_source_photo_opportunity_ids":7,
    }))
    roster=tmp_path/"rgfca_1000_species_only.csv"
    with roster.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","source_photo_opportunities",
                    "opportunity_stratum","enrolment_hash"])
        w.writerows([[11,"Species alba",3,"synthetic_tier","orig_A"],
                     [12,"Species rubra",2,"synthetic_tier","orig_B"]])
    receipt=tmp_path/"rgfca_receipt.json"
    receipt.write_text(json.dumps({
        "schema":"fcp_rgfca_photo_polymorphism_balanced_sample_v1",
        "status":"ORIGINAL_SOURCE_BALANCED_REPETITION_PLAN_NO_COLOUR_MEASURED",
        "selected_species":2,
        "selected_species_csv_sha256":sha256_file(roster),
        "new_flower_images_opened":False,
        "old_149900_flower_colours_used":False,
        "photographic_poly_probabilities_measured":False,
    }))
    taxa=gz_file(tmp_path/"taxa.csv.gz",["taxon_id","ancestry","rank"],[
        [11,"1/5","species"],[12,"1/5","species"],
        [13,"1/5","species"],[14,"1/5","species"],
        [111,"1/5/11","subspecies"],[99,"1/5","species"],
    ])
    obs_rows=[
        ["o-1",111,1,35.6,135.2,15,"Research Grade","2023-08-01"],
        ["o-2",12,2,36,136,25,"research","2022-08-05"],
        ["o-3",13,3,36.2,136.2,15,"Research Grade","2020-07-30"],
        ["o-4",14,4,37.1,137.1,15,"Research Grade","2021-08-03"],
        ["o-5",99,5,35.1,130,20,"Research Grade","2021-08-03"],
        ["o-6",12,6,36.3,136.3,20,"Needs ID","2021-08-03"],
    ]
    observations=gz_file(tmp_path/"observations.csv.gz",
         ["observation_uuid","taxon_id","observer_id","latitude","longitude",
          "positional_accuracy","quality_grade","observed_on"],obs_rows)
    photo_rows=[
        [100,"o-1",9,"JPG","CC-BY"],
        [101,"o-1",1,"jpeg","CC-BY-NC"],
        [103,"o-1",1,"png","CC0"],
        [200,"o-1",1,"jpg","CC-BY"],
        [102,"o-2",2,"jpg","CC-BY-SA"],
        [202,"o-2",2,"jpg","CC-BY"],
        [202,"o-2",2,"jpg","CC-BY"],
        [200,"o-3",3,"jpg","CC-BY"], # Important: same source ID, different OBSERVER/TAXON
        [300,"o-3",3,"jpg","CC-BY"],
        [400,"o-4",4,"jpg","CC-BY"],
        [900,"o-5",5,"jpg","CC-BY"], # Outside full 42111 species set
        [910,"o-6",6,"jpg","CC-BY"], # Not Research Grade
        [12345,"o-1",1,"jpg","CC-BY-ND"], # Invalid licence, excluded as original
    ]
    if invalid_external_collision:
        photo_rows.append([100,"o-6",6,"JPG","CC-BY"])
    photos=gz_file(tmp_path/"photos.csv.gz",
         ["photo_id","observation_uuid","observer_id","extension","license"],
         photo_rows)
    observers=gz_file(tmp_path/"observers.csv.gz",["observer_id","login","name"],[
        [1,"alice","Alice Flora"],[2,"bob","Bob Field"],
        [3,"c","C Specimen"],[4,"d","D Specimen"],
        [5,"e","E Specimen"],[6,"f","F Specimen"],
        [9,"photoauthor","Independent Photographer Nine"],
    ])
    tables={"taxa":taxa,"observations":observations,"photos":photos,"observers":observers}
    digest={k:sha256_file(v) for k,v in tables.items()}
    pin=tmp_path/"fcp_source_archive_pin.json"
    pin.write_text(json.dumps({
        "schema":"fcp_dated_inaturalist_metadata_bundle_pin_v1",
        "status":"FOUR_ORIGINAL_METADATA_TABLES_EXTRACTED_AND_PINNED",
        "snapshot_date":"2026-09-27",
        "s3_object_key":"metadata/inaturalist-open-data-20260927.tar.gz",
        "metadata_files_sha256":digest,
        "full_bundle_sha256":"a"*64,
        "source_images_opened":False,
        "historical_photo_colours_read":False,
    }))
    return frame,counts,independent,roster,receipt,pin,tables


def test_exact_date_selected_2_of_4_taxa_photo_IDs_match_globally_deduplicated_original(tmp_path):
    frame,counts,audit,roster,receipt,pin,tables=setup(tmp_path)
    destination=tmp_path/"selected_2"
    result=process(frame,counts,audit,roster,receipt,pin,tables,destination,synthetic=True)
    assert result["status"]=="COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"
    assert result["full_original_frame_taxa"]==4
    assert result["enrolled_outcome_blind_taxa"]==2
    assert result["n_selected_original_photo_IDs"]==5
    assert result["n_conservative_within_selected_photo_ID_conflicts_before_global_replay"]==0
    assert result["n_original_photo_rows_scanned_twice"]==13
    assert result["all_enrolled_species_match_independent_full_source_photo_ID_counts"] is True
    assert result["historical_flower_colour_CSV_imported"] is False
    assert result["original_image_pixels_opened"] is False
    with (destination/"rgfca_1000_complete_original_source_photo_IDs.csv").open(newline="") as f:
        rows=list(csv.DictReader(f))
    assert tuple(rows[0])==CSV_FIELDS
    assert {(r["inat_taxon_id"],r["photo_id"]) for r in rows}=={
        ("11","100"),("11","101"),("11","103"),("12","102"),("12","202")}
    assert len({r["photo_id"] for r in rows})==5
    assert "200" not in {r["photo_id"] for r in rows}
    creator=next(r for r in rows if r["photo_id"]=="100")
    assert "Independent Photographer Nine" in creator["attribution"]
    assert creator["observer_id"]=="1" and creator["photo_creator_id"]=="9"
    assert creator["photo_url_large"].endswith("/100/large.JPG")
    assert creator["photo_license"]=="cc-by"
    assert all(r["source_snapshot_date"]=="2026-09-27" for r in rows)


def test_all_other_original_photo_ID_conflicts_are_conservatively_quarantined_not_admitted(tmp_path):
    frame,counts,audit,roster,receipt,pin,tables=setup(
        tmp_path,invalid_external_collision=True)
    # Source photo 100 occurs in the selected source + outside legitimate
    # Research Grade observations; global 42111-source census would keep it,
    # but our conservative raw collision pass must REFUSE to overclaim.
    dest=tmp_path/"collision_not_fcp"
    with pytest.raises(ValueError,match="globally audited full universe"):
        process(frame,counts,audit,roster,receipt,pin,tables,dest,synthetic=True)
    doc=json.loads((dest/"FAILED_selected_source_taxon_count_mismatches.json").read_text())
    assert doc["status"]=="REJECTED_NOT_ORIGINAL_PHOTO_ID_FCP"
    assert doc["n_mismatched_selected_taxa"]>=1
    assert not (dest/"rgfca_1000_targeted_monthly_photo_ID_receipt.json").exists()


def test_rgfca_source_only_selection_metadata_cannot_be_replaced_by_colour_or_literature(tmp_path):
    frame,counts,audit,roster,receipt,pin,tables=setup(tmp_path)
    names=roster.read_text().replace("Species alba","Species fake")
    roster.write_text(names)
    with pytest.raises(ValueError,match="original source photo counts"):
        process(frame,counts,audit,roster,receipt,pin,tables,
                tmp_path/"bad",synthetic=True)


def test_immutable_original_dated_archive_and_all_four_table_sha_required(tmp_path):
    frame,counts,audit,roster,receipt,pin,tables=setup(tmp_path)
    doc=json.loads(pin.read_text())
    doc["s3_object_key"]="metadata/inaturalist-open-data-latest.tar.gz"
    pin.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="Dated immutable"):
        process(frame,counts,audit,roster,receipt,pin,tables,
                tmp_path/"not_monthly",synthetic=True)


def test_never_recreate_or_overwrite_completed_rgfca_photo_ID_ledger(tmp_path):
    args=setup(tmp_path)
    outdir=tmp_path/"done"
    process(*args,outdir,synthetic=True)
    with pytest.raises(ValueError,match="Never overwrite"):
        process(*args,outdir,synthetic=True)
