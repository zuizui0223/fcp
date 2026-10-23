"""End-to-end mini true-dated monthly metadata -> 1000-taxon source ID subset
-> RGFCA outcome-blind photo first screen -> 128 location-blind ROI worker files.

Tests have no network or old photo outcomes. The real frozen 2026-09-27
source tar is not available in the hosted test environment.
"""
from __future__ import annotations
import csv
import json

import pytest

from scripts.acquisition.extract_rgfca_1000_original_monthly_photo_ids import process
from scripts.acquisition.plan_rgfca_screen_from_targeted_original_monthly_photo_ids import create
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import freeze
from fcp_pipeline.photo_polymorphism_screen import sha256
from test_extract_rgfca_1000_original_monthly_photo_ids import setup


def fixture(tmp_path):
    frame,counts,audit,roster,receipt,pin,tables=setup(tmp_path)
    sourced=tmp_path/"exact_dated_2taxa_from_original"
    process(frame,counts,audit,roster,receipt,pin,tables,sourced,synthetic=True)
    return (
        sourced/"rgfca_1000_complete_original_source_photo_IDs.csv",
        sourced/"rgfca_1000_targeted_monthly_photo_ID_receipt.json",
        frame,counts,audit,roster,receipt
    )


def table(path):
    with path.open("r",encoding="utf-8",newline="") as src:
        return list(csv.DictReader(src))


def test_real_shaped_targeted_original_IDs_to_frozen_blind_ROI_worker(tmp_path):
    manifest,admission,frame,counts,audit,roster,roster_receipt=fixture(tmp_path)
    photo_screen=tmp_path/"photo_first_screen"
    result=create(manifest,admission,frame,counts,audit,roster,roster_receipt,
                  photo_screen,photos_per_fold=2,max_per_observer=2,
                  cell_candidate_limit=3,synthetic=True)
    assert result["source_species"]==4
    assert result["RGFCA_enrolment_species_count"]==2
    assert result["original_RGFCA_1000_full_taxon_photo_IDs_checked"]==5
    assert result["original_source_photo_ids_globally_considered"]==7
    assert result["outcome_blind_image_selection_no_floral_colour_opened"] is True
    assert result["selected_sample_photo_pixels_measured"] is False
    chosen=table(photo_screen/"photo_screen_original_metadata_selected.csv")
    assert len(chosen)==result["selected_photo_IDs"]
    assert chosen
    assert set(x["photo_id"] for x in chosen).issubset({"100","101","103","102","202"})
    assert not set(x["photo_id"] for x in chosen) & {"200","300","400"}
    assert all("morph" not in x and "colour" not in x for x in chosen)
    original=table(photo_screen/"photo_screen_all_species_opportunities.csv")
    assert len(original)==4
    for row in original:
        if row["inat_taxon_id"] in ("13","14"):
            assert row["total_selected_for_image_screen"]=="0"
            assert row["screen_unmeasured_is_not_monomorphism"]=="True"

    firewall=tmp_path/"blind_original_roi"
    fixed=freeze(photo_screen/"photo_screen_original_metadata_selected.csv",
                 photo_screen/"photo_screen_source_selection_receipt.json",
                 frame,counts,firewall,_synthetic=True)
    assert fixed["worker_partitions"]==128
    assert fixed["screen_selected_original_photo_IDs"]==len(chosen)
    assert fixed["image_worker_sees_only_three_fields"] is True
    assert fixed["original_image_pixels_opened"] is False
    assert fixed["no_old_photo_colour_rows_loaded"] is True
    assert len(list(firewall.glob("worker_partition_*.csv")))==128
    assert len(table(firewall/"metadata_join_key.csv"))==len(chosen)
    worker=next(p for p in firewall.glob("worker_partition_*.csv")
                if len(table(p)))
    assert list(table(worker)[0])==["measurement_id","image_filename","photo_license"]
    assert not set(table(worker)[0])&{"inat_taxon_id","species","latitude","longitude","morph"}


def test_target_monthly_source_manifest_digest_mutation_blocked(tmp_path):
    manifest,admission,frame,counts,audit,roster,receipt=fixture(tmp_path)
    manifest.write_text(manifest.read_text().replace("Species alba","Species faked"))
    with pytest.raises(ValueError,match="full original dated-source ID/collision provenance"):
        create(manifest,admission,frame,counts,audit,roster,receipt,tmp_path/"bad",
               photos_per_fold=2,max_per_observer=2,synthetic=True)


def test_forged_original_pixel_measurement_claim_blocks_source_photo_screen(tmp_path):
    manifest,admission,frame,counts,audit,roster,receipt=fixture(tmp_path)
    reported=json.loads(admission.read_text())
    reported["photos_NOT_classified_or_species_FCP_inferred"]=False
    admission.write_text(json.dumps(reported))
    with pytest.raises(ValueError,match="full original dated-source ID/collision provenance"):
        create(manifest,admission,frame,counts,audit,roster,receipt,tmp_path/"fake",
               photos_per_fold=2,max_per_observer=2,synthetic=True)


def test_not_source_selected_taxa_never_become_false_monomorphic(tmp_path):
    manifest,admission,frame,counts,audit,roster,receipt=fixture(tmp_path)
    out=tmp_path/"screen"
    create(manifest,admission,frame,counts,audit,roster,receipt,
           out,photos_per_fold=2,max_per_observer=2,synthetic=True)
    rows=table(out/"photo_screen_all_species_opportunities.csv")
    for row in rows:
        assert row["screen_unmeasured_is_not_monomorphism"]=="True"
    with pytest.raises(ValueError,match="Cannot replace"):
        create(manifest,admission,frame,counts,audit,roster,receipt,
               out,photos_per_fold=2,max_per_observer=2,synthetic=True)
