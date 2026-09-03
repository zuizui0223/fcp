#!/usr/bin/env python3
"""Measure one location-blind fresh photo partition with frozen ROI-v4 + palette.

The process receives only blinded image filenames and licences. It never receives
source URLs, species, coordinates, dates, observers, climate, or pollinator data.
Every input row reaches one terminal morph/status; ROI failures are retained as
``mixed_uncertain`` and never replaced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageOps

from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
from fcp_pipeline.flower_roi_v4_runtime import (
    FrozenFlowerColourEstimator,
    file_sha256,
    validate_scaleout_authorization,
)
from fcp_pipeline.photo_first_measurement import (
    BIOLOGICAL_PALETTE,
    NUISANCE,
    REFERENCE_RGB,
    classify_masked_rgb,
)
from fcp_pipeline.photo_first_measurement_execution import (
    WORKER_FIELDS,
    validate_execution_contract,
    validate_terminal_partition_results,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MEASUREMENT = ROOT / "docs/supporting/random_photo_first_measurement_contract_v1.json"
DEFAULT_EXECUTION = ROOT / "docs/supporting/random_photo_first_measurement_execution_v1.json"
DEFAULT_ROI = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
DEFAULT_LOCKED = ROOT / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-manifest", type=Path, required=True)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--detector-weight", type=Path, required=True)
    parser.add_argument("--efficient-sam-dir", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT)
    parser.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    parser.add_argument("--roi-contract", type=Path, default=DEFAULT_ROI)
    parser.add_argument("--locked-result", type=Path, default=DEFAULT_LOCKED)
    return parser.parse_args()


def _empty_palette() -> tuple[dict[str, int], dict[str, float]]:
    return (
        {name: 0 for name in REFERENCE_RGB},
        {name: 0.0 for name in BIOLOGICAL_PALETTE},
    )


def main() -> int:
    args = parse_args()
    measurement_contract = json.loads(args.measurement_contract.read_text(encoding="utf-8"))
    execution_contract = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    roi_contract = json.loads(args.roi_contract.read_text(encoding="utf-8"))
    locked = json.loads(args.locked_result.read_text(encoding="utf-8"))
    validate_execution_contract(measurement_contract, execution_contract)
    validate_roi_v4_contract(roi_contract)
    detector_sha = file_sha256(args.detector_weight)
    expected_detector = execution_contract["runtime"]["detector_sha256"]
    if detector_sha != expected_detector:
        raise RuntimeError("fresh measurement detector hash changed")
    validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    if tuple(worker.columns) != WORKER_FIELDS:
        raise ValueError("blind worker manifest fields changed or leaked")
    if worker["measurement_id"].nunique() != len(worker):
        raise ValueError("blind worker manifest contains duplicate IDs")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result_columns = [
        "measurement_id",
        "morph",
        "measurement_status",
        "roi_status",
        "failure_reasons",
        "image_sha256",
        "flower_effective_pixels",
        "background_effective_pixels",
        "horizontal_flip_mask_iou",
        "horizontal_flip_colour_delta_e",
        "mask_pixels",
        "nuisance_pixel_fraction",
        *[f"palette_count_{name}" for name in REFERENCE_RGB],
        *[f"flower_fraction_{name}" for name in BIOLOGICAL_PALETTE],
    ]
    if len(worker) == 0:
        pd.DataFrame(columns=result_columns).to_csv(
            args.output_csv, index=False, lineterminator="\n"
        )
        print(json.dumps({"status": "complete_empty_blind_measurement", "rows": 0}, indent=2))
        return 0

    estimator = FrozenFlowerColourEstimator(
        args.detector_weight,
        args.efficient_sam_dir,
        roi_contract,
        torch_threads=int(execution_contract["runtime"]["torch_threads_per_worker"]),
    )
    threshold = measurement_contract["coarse_colour_state"]
    minimum_pixels = int(execution_contract["measurement"]["minimum_flower_mask_pixels"])
    results = []

    for position, row in enumerate(worker.itertuples(index=False), start=1):
        image_path = args.images_dir / str(row.image_filename)
        if not image_path.is_file():
            raise RuntimeError(f"blind acquisition packet is missing {row.image_filename}")
        image_hash = hashlib.sha256(image_path.read_bytes()).hexdigest()
        try:
            with Image.open(image_path) as source:
                measured = estimator.measure(source)
                oriented = ImageOps.exif_transpose(source).convert("RGB")
                rgb = np.asarray(oriented, dtype=np.uint8)
        except Exception as exc:
            counts, fractions = _empty_palette()
            result = {
                "measurement_id": str(row.measurement_id),
                "morph": "mixed_uncertain",
                "measurement_status": "not_evaluable_roi_or_flip_gate",
                "roi_status": "roi_runtime_failure",
                "failure_reasons": f"{type(exc).__name__}:{str(exc)[:350]}",
                "image_sha256": image_hash,
                "flower_effective_pixels": 0,
                "background_effective_pixels": 0,
                "horizontal_flip_mask_iou": np.nan,
                "horizontal_flip_colour_delta_e": np.nan,
                "mask_pixels": 0,
                "nuisance_pixel_fraction": np.nan,
            }
        else:
            roi_status = str(measured.get("automated_colour_state_status") or "")
            flower_pixels = int(measured.get("flower_effective_pixels") or 0)
            background_pixels = int(measured.get("background_effective_pixels") or 0)
            flip_iou = measured.get("horizontal_flip_mask_iou")
            flip_delta = measured.get("horizontal_flip_colour_delta_e")
            roi_failures = str(measured.get("failure_reasons") or "")
            if roi_status != "automated_colour_state_admitted":
                counts, fractions = _empty_palette()
                result = {
                    "measurement_id": str(row.measurement_id),
                    "morph": "mixed_uncertain",
                    "measurement_status": "not_evaluable_roi_or_flip_gate",
                    "roi_status": roi_status or "automated_colour_state_not_evaluable",
                    "failure_reasons": roi_failures or "roi_or_flip_gate_failed",
                    "image_sha256": image_hash,
                    "flower_effective_pixels": flower_pixels,
                    "background_effective_pixels": background_pixels,
                    "horizontal_flip_mask_iou": flip_iou,
                    "horizontal_flip_colour_delta_e": flip_delta,
                    "mask_pixels": flower_pixels,
                    "nuisance_pixel_fraction": np.nan,
                }
            else:
                mask = np.asarray(measured["flower_mask"], dtype=bool)
                if mask.shape != rgb.shape[:2]:
                    raise RuntimeError("ROI flower mask and oriented RGB shape drifted")
                classified = classify_masked_rgb(
                    rgb[mask],
                    minimum_mask_pixels=minimum_pixels,
                    minimum_dominant_fraction=float(threshold["minimum_dominant_fraction"]),
                    minimum_margin=float(threshold["minimum_margin_over_second_group"]),
                )
                counts = dict(classified["palette_counts"])
                fractions = dict(classified["flower_only_fractions"])
                nuisance = sum(int(counts[name]) for name in NUISANCE)
                total = int(sum(int(value) for value in counts.values()))
                result = {
                    "measurement_id": str(row.measurement_id),
                    "morph": str(classified["morph"]),
                    "measurement_status": str(classified["measurement_status"]),
                    "roi_status": roi_status,
                    "failure_reasons": (
                        "" if classified["measurement_status"] == "classified_four_state_morph"
                        else str(classified["measurement_status"])
                    ),
                    "image_sha256": image_hash,
                    "flower_effective_pixels": flower_pixels,
                    "background_effective_pixels": background_pixels,
                    "horizontal_flip_mask_iou": flip_iou,
                    "horizontal_flip_colour_delta_e": flip_delta,
                    "mask_pixels": int(classified["mask_pixels"]),
                    "nuisance_pixel_fraction": nuisance / total if total else np.nan,
                }
        result.update({f"palette_count_{name}": int(counts[name]) for name in REFERENCE_RGB})
        result.update(
            {
                f"flower_fraction_{name}": float(fractions[name])
                for name in BIOLOGICAL_PALETTE
            }
        )
        results.append(result)
        if position % 10 == 0 or position == len(worker):
            print(f"measured={position}/{len(worker)}", flush=True)

    frame = pd.DataFrame(results).reindex(columns=result_columns)
    validate_terminal_partition_results(
        frame, worker["measurement_id"].astype(str).tolist()
    )
    frame.to_csv(args.output_csv, index=False, lineterminator="\n")
    print(
        json.dumps(
            {
                "status": "complete_location_blind_random_photo_first_measurement",
                "rows": int(len(frame)),
                "classified": int(frame["morph"].isin(["white", "yellow_orange", "red_pink", "blue_purple"]).sum()),
                "mixed_uncertain": int(frame["morph"].eq("mixed_uncertain").sum()),
                "coordinates_opened": False,
                "species_opened": False,
                "source_urls_opened": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
