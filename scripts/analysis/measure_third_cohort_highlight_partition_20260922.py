#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageOps

def _srgb_to_linear(x: np.ndarray) -> np.ndarray:
    return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)

def highlight_metrics(rgb_pixels: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(rgb_pixels, dtype=np.uint8)
    if x.ndim != 2 or x.shape[1] != 3 or len(x) == 0:
        raise ValueError("flower-mask RGB pixels are unavailable")
    m = x.max(axis=1)
    clip = float(np.mean(m == 255))
    near = float(np.mean(m >= 250))
    srgb = x.astype(np.float64) / 255.0
    lin = _srgb_to_linear(srgb)
    lum = 0.2126729 * lin[:,0] + 0.7151522 * lin[:,1] + 0.0721750 * lin[:,2]
    return clip, near, float(np.quantile(lum, 0.99))

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--worker-manifest", type=Path, required=True)
    p.add_argument("--images-dir", type=Path, required=True)
    p.add_argument("--frozen-source-root", type=Path, required=True)
    p.add_argument("--efficient-sam-dir", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    args = p.parse_args()

    sys.path.insert(0, str(args.frozen_source_root))
    from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
    from fcp_pipeline.flower_roi_v4_runtime import (
        FrozenFlowerColourEstimator,
        file_sha256,
        validate_scaleout_authorization,
    )

    roi_path = args.frozen_source_root / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
    locked_path = args.frozen_source_root / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"
    detector = args.frozen_source_root / "data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    roi = json.loads(roi_path.read_text(encoding="utf-8"))
    locked = json.loads(locked_path.read_text(encoding="utf-8"))
    validate_roi_v4_contract(roi)
    detector_sha = file_sha256(detector)
    validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    if list(worker.columns) != ["measurement_id","image_filename","photo_license"]:
        raise RuntimeError("blind worker schema drift")
    estimator = FrozenFlowerColourEstimator(
        detector, args.efficient_sam_dir, roi, torch_threads=1
    )

    rows = []
    for i, row in enumerate(worker.itertuples(index=False), start=1):
        path = args.images_dir / str(row.image_filename)
        image_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        rec = {
            "measurement_id": str(row.measurement_id),
            "reacquired_image_sha256": image_hash,
            "technical_status": "",
            "roi_status": "",
            "failure_reason": "",
            "mask_pixels": 0,
            "clip_fraction": np.nan,
            "near_clip_fraction": np.nan,
            "luminance_q99": np.nan,
        }
        try:
            with Image.open(path) as source:
                measured = estimator.measure(source)
                rgb = np.asarray(ImageOps.exif_transpose(source).convert("RGB"), dtype=np.uint8)
            mask = np.asarray(measured["flower_mask"], dtype=bool)
            rec["roi_status"] = str(measured.get("automated_colour_state_status") or "")
            rec["mask_pixels"] = int(mask.sum())
            if mask.shape != rgb.shape[:2] or not mask.any():
                raise RuntimeError("flower_mask_unavailable")
            clip, near, lum = highlight_metrics(rgb[mask])
            rec.update({
                "technical_status": "highlight_metrics_available",
                "clip_fraction": clip,
                "near_clip_fraction": near,
                "luminance_q99": lum,
            })
        except Exception as exc:
            rec["technical_status"] = "highlight_metrics_unavailable"
            rec["failure_reason"] = f"{type(exc).__name__}:{str(exc)[:300]}"
        rows.append(rec)
        if i % 10 == 0 or i == len(worker):
            print(f"highlight_measured={i}/{len(worker)}", flush=True)

    out = pd.DataFrame(rows)
    if out["measurement_id"].nunique() != len(out):
        raise RuntimeError("duplicate technical IDs")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False, lineterminator="\n")

if __name__ == "__main__":
    main()
