"""Synthetic contract tests: photo-first screen workers never see taxon or geography."""
from __future__ import annotations

import csv
import json

import pytest

from fcp_pipeline.photo_polymorphism_screen import sha256
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import (
    ACQUISITION, FREEZE_STATUS, METADATA, WORKER, freeze,
)
from scripts.acquisition.plan_fcp_photo_first_screen import plan
from test_plan_fcp_photo_first_screen import setup_source


def source_screen(tmp_path):
    stage, frame, counts, _ = setup_source(tmp_path)
    sample = tmp_path / "selection"
    plan(stage,frame,counts,sample,photos_per_fold=2,cell_candidate_limit=4,
         max_per_observer=2,_synthetic=True)
    return (
        sample / "photo_screen_original_metadata_selected.csv",
        sample / "photo_screen_source_selection_receipt.json",frame,counts
    )


def rows(path):
    with path.open("r",newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_128_original_roi_workers_are_exactly_location_blind(tmp_path):
    selected,receipt,frame,counts=source_screen(tmp_path)
    out=tmp_path/"blind"
    result=freeze(selected,receipt,frame,counts,out,_synthetic=True)
    assert result["status"]==FREEZE_STATUS
    assert result["worker_partitions"]==128
    assert result["screen_selected_original_photo_IDs"] > 0
    assert result["screen_selected_original_photo_IDs"] < 20
    assert result["source_screen_NOT_complete_photo_census"] is True
    assert result["original_image_pixels_opened"] is False
    assert result["colour_states_inferred"] is False
    assert result["no_literature_positive_labels_used"] is True
    assert result["no_old_photo_colour_rows_loaded"] is True
    assert len(list(out.glob("worker_partition_*.csv")))==128
    assert len(list(out.glob("acquisition_partition_*.csv")))==128
    sum_worker=0
    for n in range(128):
        worker=out/f"worker_partition_{n:03}.csv"
        acquisition=out/f"acquisition_partition_{n:03}.csv"
        assert worker.open().readline().strip().split(",")==list(WORKER)
        assert acquisition.open().readline().strip().split(",")==list(ACQUISITION)
        w=rows(worker)
        a=rows(acquisition)
        assert len(w)==len(a)
        assert {x["measurement_id"] for x in w}=={x["measurement_id"] for x in a}
        sum_worker+=len(w)
        assert result["frozen_file_sha256"][worker.name]==sha256(worker)
        assert result["frozen_file_sha256"][acquisition.name]==sha256(acquisition)
        assert not set(WORKER) & {
            "species","inat_taxon_id","photo_id","latitude","longitude","photo_url_large",
        }
    assert sum_worker==result["screen_selected_original_photo_IDs"]
    key=out/"metadata_join_key.csv"
    assert key.open().readline().strip().split(",")==list(METADATA)
    assert len(rows(key))==sum_worker
    assert result["photo_screen_fold_0"]+result["photo_screen_fold_1"]==sum_worker
    assert not (out/"fresh_study_freeze.json").exists()


def test_freeze_rejects_image_outcome_inserted_into_immutable_metadata(tmp_path):
    selected,receipt,frame,counts=source_screen(tmp_path)
    with selected.open(newline="") as f:
        data=list(csv.DictReader(f))
    keys=list(data[0])+["morph"]
    for rec in data:
        rec["morph"]="white"
    with selected.open("w",newline="") as f:
        w=csv.DictWriter(f,keys)
        w.writeheader()
        w.writerows(data)
    info=json.loads(receipt.read_text())
    info["photo_selection_csv_sha256"]=sha256(selected)
    receipt.write_text(json.dumps(info))
    with pytest.raises(ValueError,match="selection columns"):
        freeze(selected,receipt,frame,counts,tmp_path/"invalid",_synthetic=True)


def test_changed_original_photo_selection_digest_rejected(tmp_path):
    selected,receipt,frame,counts=source_screen(tmp_path)
    selected.write_text(selected.read_text().replace("Spec alba","Spec altered"))
    with pytest.raises(ValueError,match="no valid pre-colour"):
        freeze(selected,receipt,frame,counts,tmp_path/"invalid",_synthetic=True)


def test_tampered_source_colour_independence_cannot_freeze(tmp_path):
    selected,receipt,frame,counts=source_screen(tmp_path)
    x=json.loads(receipt.read_text())
    x["literature_polymorphism_labels_used"]=True
    receipt.write_text(json.dumps(x))
    with pytest.raises(ValueError,match="no valid pre-colour"):
        freeze(selected,receipt,frame,counts,tmp_path/"invalid",_synthetic=True)


def test_no_overwrite_of_existing_source_blind_packets(tmp_path):
    selected,receipt,frame,counts=source_screen(tmp_path)
    out=tmp_path/"blind"
    freeze(selected,receipt,frame,counts,out,_synthetic=True)
    with pytest.raises(ValueError,match="Cannot overwrite"):
        freeze(selected,receipt,frame,counts,out,_synthetic=True)
