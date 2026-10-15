#!/usr/bin/env python3
"""Run one pre-colour 42,111-candidate photo-screen ROI work unit.

Exact original ROI-v4/EfficientSAM worker is reused UNCHANGED through the
existing image-partition runtime; only its metadata preflight targets a
separately frozen bounded discovery/holdout screen, not all 93M photographs.
A completed screen batch is never an all-species biological inference.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_polymorphism_screen import (
    FRAME_SHA, SOURCE_COUNT_SHA, SOURCE_PHOTOS, SNAPSHOT, sha256,
)
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import (
    ACQUISITION, FREEZE_STATUS, WORKER,
)
from scripts.acquisition.run_fcp_fresh_image_partition import run_batch


class ScreenSourceIndex:
    def __init__(self, firewall: Path, part: int, *, synthetic: bool = False):
        self.firewall = firewall
        self.part = part
        self.freeze_path = firewall / "photo_screen_freeze.json"
        self.receipt = json.loads(self.freeze_path.read_text(encoding="utf-8"))
        doc = self.receipt
        if (doc.get("schema") != "fcp_original_photo_screen_blind_firewall_v1"
            or doc.get("status") != FREEZE_STATUS
            or doc.get("source_snapshot_date") != SNAPSHOT
            or doc.get("worker_partitions") != 128
            or not 0 <= part < 128
            or doc.get("image_worker_sees_only_three_fields") is not True
            or doc.get("source_screen_NOT_complete_photo_census") is not True
            or doc.get("original_image_pixels_opened") is not False
            or doc.get("colour_states_inferred") is not False
            or doc.get("no_old_photo_colour_rows_loaded") is not True
            or doc.get("no_literature_positive_labels_used") is not True
            or (not synthetic and
                (doc.get("source_species") != 42111
                 or doc.get("source_frame_sha256") != FRAME_SHA
                 or doc.get("source_photo_counts_sha256") != SOURCE_COUNT_SHA
                 or doc.get("source_all_provider_photo_IDs") != SOURCE_PHOTOS))):
            raise ValueError("No valid original photo-screen-only frozen ROI firewall")
        worker_path = firewall / f"worker_partition_{part:03}.csv"
        acquisition_path = firewall / f"acquisition_partition_{part:03}.csv"
        for path in (worker_path, acquisition_path):
            if (doc.get("frozen_file_sha256",{}).get(path.name) != sha256(path)):
                raise ValueError("Source/photo URL worker mutated since original pre-colour screen freeze")
        self.workers = pd.read_csv(worker_path, dtype=str).fillna("")
        self.keys = pd.read_csv(acquisition_path, dtype=str).fillna("")
        if (tuple(self.workers.columns) != WORKER
            or tuple(self.keys.columns) != ACQUISITION
            or self.workers.measurement_id.duplicated().any()
            or self.keys.measurement_id.duplicated().any()
            or len(self.workers) != len(self.keys)
            or set(self.workers.measurement_id) != set(self.keys.measurement_id)):
            raise ValueError("Original selected blind-worker packet/source acquisition key invalid")
        self.freeze_sha256 = sha256(self.freeze_path)

    def get_batch(self, part: int, batch: int, size: int):
        if part != self.part or not 0 <= batch or not 1 <= size <= 1024:
            raise ValueError("Screen ROI partition/batch identity invalid")
        start = batch * size
        if start >= len(self.workers):
            raise ValueError("Screen partition has no work unit at requested index")
        selected = self.workers.iloc[start:start+size].copy()
        original = selected.merge(
            self.keys, on=["measurement_id","image_filename","photo_license"],
            how="left", validate="one_to_one")
        if len(original) != len(selected) or original.photo_url_large.isna().any():
            raise ValueError("Frozen original S3 photo URL absent")
        return selected, original, math.ceil(len(self.workers) / size)


def execute(firewall: Path, part: int, batch: int, batch_size: int,
            outdir: Path, *, source: Path, detector: Path, sam: Path,
            _synthetic_loader=None, _synthetic_evaluator=None) -> dict:
    if batch_size != 64:
        raise ValueError("Image screen runs fixed 64-photo ROI compute work units, NOT taxon caps")
    test_only = _synthetic_loader is not None or _synthetic_evaluator is not None
    if test_only and (_synthetic_loader is None or _synthetic_evaluator is None):
        raise ValueError("Synthetic test requires paired loader and evaluator")
    index = ScreenSourceIndex(firewall, part, synthetic=test_only)
    record = run_batch(
        firewall,part,batch,batch_size,outdir,
        source=source,detector=detector,sam=sam,
        indexed_source=index,
        loader=_synthetic_loader,evaluator=_synthetic_evaluator,
        _test_only_skip_method_check=test_only)
    if (record.get("status") != "COMPLETE_WORK_UNIT"
        or record.get("full_study_complete") is not False
        or record.get("reused_historical_measurement_csv") is not False
        or record.get("image_bytes_persisted") is not False):
        raise ValueError("Original screen terminal incorrectly claims all-photo completion")
    receipt = {
        "schema": "fcp_photo_screen_bounded_original_roi_work_unit_v1",
        "status": "ORIGINAL_PHOTO_SCREEN_UNIT_MEASURED_NOT_ALL_FCP",
        "screen_firewall_sha256": index.freeze_sha256,
        "partition": part, "batch": batch, "batch_size": batch_size,
        "source_photo_IDs": record["n_source_photo_ids"],
        "n_original_photo_classified": record["n_classifiable"],
        "terminal_sha256": record["terminal_csv_sha256"],
        "frozen_original_roi_method": record["method"],
        "historical_photo_colours_reused": False,
        "original_photo_screen_complete": False,
        "full_93million_original_photo_census_complete": False,
        "W_H2_or_environment_results_computed": False,
        "synthetic_test_only": test_only,
    }
    name = f"screen_p{part:03}_b{batch:07}.json"
    target = outdir / name
    if target.exists():
        raise ValueError("Cannot overwrite original photo-screen work-unit receipt")
    target.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return receipt


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--firewall",required=True,type=Path)
    p.add_argument("--partition",type=int,required=True)
    p.add_argument("--batch",type=int,required=True)
    p.add_argument("--batch-size",type=int,default=64)
    p.add_argument("--original-measure-script",required=True,type=Path)
    p.add_argument("--original-detector",required=True,type=Path)
    p.add_argument("--original-efficient-sam",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args(argv)
    print(json.dumps(execute(args.firewall,args.partition,args.batch,
        args.batch_size,args.output,source=args.original_measure_script,
        detector=args.original_detector,sam=args.original_efficient_sam),
        sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
