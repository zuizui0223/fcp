"""Synthetic 128-worker original-image SCREEN reassembly, not biological data."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.acquisition.reassemble_fcp_photo_first_screen import join
from scripts.acquisition.run_fcp_photo_first_roi_work_unit import execute
from test_run_fcp_photo_first_roi_work_unit import screen_firewall, simulated_roi


def build_sealed_screen(tmp_path):
    firewall,_=screen_firewall(tmp_path)
    result_folder=tmp_path/"original_photo_screen_roi"
    result_folder.mkdir()
    completed=[]
    for part in range(128):
        worker=firewall/f"worker_partition_{part:03}.csv"
        with worker.open(newline="") as f:
            n=sum(1 for _ in csv.DictReader(f))
        if not n:
            continue
        output=result_folder/f"partition_{part:03}"
        got=execute(firewall,part,0,64,output,
                    source=Path("synthetic"),
                    detector=Path("synthetic_detector"),
                    sam=Path("synthetic_sam"),
                    _synthetic_loader=lambda url:b"dummy-image-payload",
                    _synthetic_evaluator=simulated_roi)
        assert got["source_photo_IDs"]==n
        completed.append((part,output))
    return firewall,result_folder,completed


def load(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def test_screen_photo_identity_full_selected_sample_reassembles_without_hiding_failures(tmp_path):
    firewall,folder,completed=build_sealed_screen(tmp_path)
    out=tmp_path/"new_original_screen"
    r=join(firewall,folder,out,_synthetic_test=True)
    assert r["status"]=="SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL"
    assert r["screen_not_entire_original_93million_census"] is True
    assert r["full_source_image_measurement_complete"] is False
    assert r["historical_measured_flower_colours_reused"] is False
    assert r["genetic_FCP_claim"] is False
    assert r["synthetic_test_only"] is True
    assert r["original_ROI_v4_method_attested"] is False
    assert r["n_terminal_measurements"]==r["n_classifiable_original_flower_photos"]
    original=load(out/"newly_remeasured_photo_screen.csv")
    assert len(original)==r["n_terminal_measurements"]
    assert all(x["global_classifiable"]=="true" for x in original)
    assert all(x["morph"]=="white" for x in original)
    assert all(x["screen_fold"] in {"0","1"} for x in original)
    assert len({x["photo_id"] for x in original})==len(original)
    assert (out/"sealed_photo_screen_terminal.json").is_file()


def test_screen_requires_every_original_selected_ROI_photo_not_partial(tmp_path):
    firewall,folder,completed=build_sealed_screen(tmp_path)
    part,output=completed[0]
    file=output/f"terminal_p{part:03}_b0000000.csv"
    file.unlink()
    with pytest.raises(ValueError,match="Extra/missing"):
        join(firewall,folder,tmp_path/"incomplete",_synthetic_test=True)
    assert not (tmp_path/"incomplete/newly_remeasured_photo_screen.csv").exists()
    assert not (tmp_path/"incomplete/sealed_photo_screen_terminal.json").exists()


def test_screen_rejects_faked_full_photo_census_terminal(tmp_path):
    firewall,folder,completed=build_sealed_screen(tmp_path)
    part,output=completed[0]
    receipt=output/f"terminal_p{part:03}_b0000000.json"
    doc=json.loads(receipt.read_text(encoding="utf-8"))
    doc["full_study_complete"]=True
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="work unit is incomplete or changed"):
        join(firewall,folder,tmp_path/"not_global",_synthetic_test=True)


def test_screen_rejects_sample_source_identity_tampering(tmp_path):
    firewall,folder,completed=build_sealed_screen(tmp_path)
    data=firewall/"metadata_join_key.csv"
    data.write_text(data.read_text().replace("Spec alba","Another species"))
    with pytest.raises(ValueError,match="Frozen screen worker"):
        join(firewall,folder,tmp_path/"changed",_synthetic_test=True)


def test_screen_refuses_reassembly_overwrite(tmp_path):
    firewall,folder,completed=build_sealed_screen(tmp_path)
    out=tmp_path/"existing"
    join(firewall,folder,out,_synthetic_test=True)
    with pytest.raises(ValueError,match="Do not overwrite"):
        join(firewall,folder,out,_synthetic_test=True)
