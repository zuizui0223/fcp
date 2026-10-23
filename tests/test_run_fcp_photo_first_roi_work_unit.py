"""Photo-first ORIGINAL ROI-v4 work-unit tests with synthetic image evaluator.

No actual provider requests or biological polymorphism measurements.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from fcp_pipeline.fcp_fresh_study import BIO
from scripts.acquisition.run_fcp_photo_first_roi_work_unit import (
    ScreenSourceIndex, execute,
)
from test_freeze_fcp_photo_first_blind_workers import source_screen
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import freeze


def screen_firewall(tmp_path):
    selected,receipt,frame,counts=source_screen(tmp_path)
    output=tmp_path/"blind_original_screen"
    frozen=freeze(selected,receipt,frame,counts,output,_synthetic=True)
    assert frozen["screen_selected_original_photo_IDs"] > 0
    part=next(n for n in range(128) if sum(1 for _ in
              (output/f"worker_partition_{n:03}.csv").open())>1)
    return output,part


def simulated_roi(worker, images, source, detector, sam, out):
    data=pd.read_csv(worker,dtype=str)
    result=pd.DataFrame({
        "measurement_id":data.measurement_id,
        "morph":["white"]*len(data),
        "measurement_status":["classified_four_state_morph"]*len(data),
    })
    for c in BIO:
        result[f"flower_fraction_{c}"]=1.0 if c=="white" else 0.0
    return result


def test_screen_worker_downloads_selected_originals_only_not_global_full_census(tmp_path):
    fw,part=screen_firewall(tmp_path)
    index=ScreenSourceIndex(fw,part,synthetic=True)
    worker,acq,n_batch=index.get_batch(part,0,64)
    assert len(worker)==len(acq)>0 and n_batch==1
    assert list(worker.columns)==["measurement_id","image_filename","photo_license"]
    assert "photo_url_large" not in worker
    assert acq.photo_url_large.str.startswith(
        "https://inaturalist-open-data.s3.amazonaws.com/photos/").all()
    output=tmp_path/"test_ROI_batch"
    result=execute(fw,part,0,64,output,
                   source=Path("no_original_source"),
                   detector=Path("no_real_detector"),
                   sam=Path("no_real_sam"),
                   _synthetic_loader=lambda url:b"synthetic-image-bytes",
                   _synthetic_evaluator=simulated_roi)
    assert result["status"]=="ORIGINAL_PHOTO_SCREEN_UNIT_MEASURED_NOT_ALL_FCP"
    assert result["source_photo_IDs"]==len(worker)
    assert result["n_original_photo_classified"]==len(worker)
    assert result["synthetic_test_only"] is True
    assert result["original_photo_screen_complete"] is False
    assert result["full_93million_original_photo_census_complete"] is False
    assert result["W_H2_or_environment_results_computed"] is False
    assert (output/f"screen_p{part:03}_b0000000.json").is_file()
    assert (output/f"terminal_p{part:03}_b0000000.csv").is_file()


def test_bounded_screen_refuses_forged_full_study_or_wrong_partition(tmp_path):
    fw,part=screen_firewall(tmp_path)
    index=ScreenSourceIndex(fw,part,synthetic=True)
    with pytest.raises(ValueError,match="partition/batch"):
        index.get_batch((part+1)%128,0,64)
    with pytest.raises(ValueError,match="no work unit"):
        index.get_batch(part,1000,64)
    doc=json.loads((fw/"photo_screen_freeze.json").read_text())
    doc["status"]="FRESH_BLIND_METADATA_FROZEN"
    (fw/"photo_screen_freeze.json").write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="photo-screen-only"):
        ScreenSourceIndex(fw,part,synthetic=True)


def test_screen_requires_immutable_source_provenance_and_bounds(tmp_path):
    fw,part=screen_firewall(tmp_path)
    path=fw/f"worker_partition_{part:03}.csv"
    path.write_text(path.read_text()+"TAMPERED_ROW\n")
    with pytest.raises(ValueError,match="worker mutated"):
        ScreenSourceIndex(fw,part,synthetic=True)


def test_roi_cannot_be_relabelled_complete_after_one_64_image_batch(tmp_path):
    fw,part=screen_firewall(tmp_path)
    out=tmp_path/"ROI"
    kwargs={
        "source":Path("no_original_source"),
        "detector":Path("no_real_detector"),
        "sam":Path("no_real_sam"),
        "_synthetic_loader":lambda url:b"synthetic-image-bytes",
        "_synthetic_evaluator":simulated_roi,
    }
    execute(fw,part,0,64,out,**kwargs)
    with pytest.raises(ValueError,match="already exists"):
        execute(fw,part,0,64,out,**kwargs)
    with pytest.raises(ValueError,match="fixed 64-photo"):
        execute(fw,part,0,128,tmp_path/"bad",**kwargs)
