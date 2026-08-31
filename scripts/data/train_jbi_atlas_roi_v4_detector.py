#!/usr/bin/env python3
"""Train the frozen flower-specific detector without validation or test access."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-dir", type=Path, required=True)
    parser.add_argument("--upstream-weight", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_roi_v4_contract(contract)
    detector = contract["detector"]
    training = detector["training"]
    materialization = json.loads(
        (args.training_dir / "training_manifest.json").read_text(encoding="utf-8")
    )
    if (
        materialization.get("status") != "pass_roi_v4_training_materialization"
        or materialization.get("images") != 400
        or materialization.get("evaluable_training_boxes") != 6991
        or materialization.get("jrc_test_directory_read") is not False
    ):
        raise RuntimeError("ROI v4 training materialization changed")
    if sha256(args.upstream_weight) != detector["upstream_weight_sha256"]:
        raise RuntimeError("YOLO11n upstream weight hash changed")

    import torch
    import ultralytics
    from ultralytics import YOLO

    if ultralytics.__version__ != detector["ultralytics_version"]:
        raise RuntimeError("Ultralytics version changed")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    started = time.time()
    model = YOLO(str(args.upstream_weight))
    model.train(
        data=str(args.training_dir / "data.yaml"),
        epochs=int(training["epochs"]),
        imgsz=int(training["image_size"]),
        batch=int(training["batch"]),
        device="cpu",
        workers=0,
        optimizer=training["optimizer"],
        lr0=float(training["initial_learning_rate"]),
        lrf=float(training["final_learning_rate_fraction"]),
        momentum=float(training["momentum"]),
        weight_decay=float(training["weight_decay"]),
        warmup_epochs=float(training["warmup_epochs"]),
        warmup_momentum=float(training["warmup_momentum"]),
        warmup_bias_lr=float(training["warmup_bias_learning_rate"]),
        box=float(training["box_loss_gain"]),
        cls=float(training["class_loss_gain"]),
        dfl=float(training["distribution_focal_loss_gain"]),
        seed=int(training["seed"]),
        deterministic=True,
        patience=0,
        val=False,
        save=True,
        save_period=10,
        plots=False,
        amp=False,
        rect=False,
        multi_scale=0.0,
        single_cls=True,
        hsv_h=float(training["hsv_h"]),
        hsv_s=float(training["hsv_s"]),
        hsv_v=float(training["hsv_v"]),
        degrees=float(training["degrees"]),
        translate=float(training["translate"]),
        scale=float(training["scale"]),
        shear=float(training["shear"]),
        perspective=float(training["perspective"]),
        flipud=float(training["vertical_flip_probability"]),
        fliplr=float(training["horizontal_flip_probability"]),
        mosaic=float(training["mosaic_probability"]),
        close_mosaic=int(training["close_mosaic_epochs"]),
        mixup=float(training["mixup_probability"]),
        cutmix=float(training["cutmix_probability"]),
        copy_paste=float(training["copy_paste_probability"]),
        project=str(args.output_dir),
        name="frozen_train",
        exist_ok=False,
        verbose=True,
    )
    source_last = args.output_dir / "frozen_train/weights/last.pt"
    if not source_last.is_file():
        raise RuntimeError("frozen last-epoch detector was not produced")
    frozen_weight = args.output_dir / "jrc_yolo11n_last_v4.pt"
    shutil.copy2(source_last, frozen_weight)
    results = args.output_dir / "frozen_train/results.csv"
    manifest = {
        "protocol": contract["protocol"],
        "status": "complete_roi_v4_detector_training_not_yet_qualified",
        "epochs": int(training["epochs"]),
        "weight_selection": training["weight_selection"],
        "trained_weight_path": frozen_weight.name,
        "trained_weight_sha256": sha256(frozen_weight),
        "training_results_sha256": sha256(results),
        "elapsed_seconds": time.time() - started,
        "jrc_test_images_decoded_or_scored": False,
        "scaleout_candidate_pixels_opened": False,
        "environment": {
            "python": platform.python_version(),
            "torch": torch.__version__,
            "ultralytics": ultralytics.__version__,
            "device": "cpu",
        },
    }
    (args.output_dir / "training_result_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
