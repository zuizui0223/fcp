#!/usr/bin/env python3
"""Recover frozen ROI-v4 background palette counts for one blinded matched partition."""
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
from fcp_pipeline.photo_first_measurement import REFERENCE_RGB, nearest_palette_counts

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXECUTION = ROOT / "docs/supporting/random_photo_first_measurement_execution_v1.json"
DEFAULT_ROI = ROOT / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
DEFAULT_LOCKED = ROOT / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected-blind-manifest", type=Path, required=True)
    ap.add_argument("--acquisition-receipt", type=Path, required=True)
    ap.add_argument("--images-dir", type=Path, required=True)
    ap.add_argument("--expected-key", type=Path, required=True)
    ap.add_argument("--detector-weight", type=Path, required=True)
    ap.add_argument("--efficient-sam-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-manifest", type=Path, required=True)
    ap.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    ap.add_argument("--roi-contract", type=Path, default=DEFAULT_ROI)
    ap.add_argument("--locked-result", type=Path, default=DEFAULT_LOCKED)
    return ap.parse_args()


def _expected_integer(value: object, *, label: str) -> int:
    """Parse integer-like CSV values without changing the frozen numeric value.

    Original measurement columns can be materialized by pandas as strings such as
    ``56230.0`` because terminal tables also contain missing values.  Require a
    finite, exactly integral number rather than silently rounding it.
    """
    numeric = float(value)
    if not np.isfinite(numeric) or numeric < 0 or not numeric.is_integer():
        raise RuntimeError(f"frozen expected {label} is not a nonnegative integer: {value!r}")
    return int(numeric)


def main() -> int:
    args = parse_args()
    execution = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    roi_contract = json.loads(args.roi_contract.read_text(encoding="utf-8"))
    locked = json.loads(args.locked_result.read_text(encoding="utf-8"))
    validate_roi_v4_contract(roi_contract)
    detector_sha = file_sha256(args.detector_weight)
    if detector_sha != execution["runtime"]["detector_sha256"]:
        raise RuntimeError("background recovery detector SHA differs from frozen runtime")
    validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)

    selected = pd.read_csv(args.selected_blind_manifest, dtype=str).fillna("")
    receipt = pd.read_csv(args.acquisition_receipt, dtype=str).fillna("")
    expected = pd.read_csv(args.expected_key, dtype=str).fillna("")
    if selected["measurement_id"].nunique() != len(selected):
        raise RuntimeError("selected blind partition contains duplicate IDs")
    if receipt["measurement_id"].nunique() != len(receipt) or set(receipt["measurement_id"]) != set(selected["measurement_id"]):
        raise RuntimeError("acquisition receipt does not exactly cover selected blind partition")
    if expected["measurement_id"].nunique() != len(expected):
        raise RuntimeError("safe expected key contains duplicate IDs")
    subset = selected[["measurement_id", "image_filename"]].merge(
        receipt[["measurement_id", "acquisition_status", "image_sha256"]],
        on="measurement_id", how="left", validate="one_to_one"
    ).merge(expected, on="measurement_id", how="left", validate="one_to_one")
    if len(subset) != len(selected) or subset["expected_image_sha256"].astype(str).str.len().eq(0).any():
        raise RuntimeError("safe expected reproduction key does not cover selected partition")

    estimator = FrozenFlowerColourEstimator(
        args.detector_weight,
        args.efficient_sam_dir,
        roi_contract,
        torch_threads=int(execution["runtime"]["torch_threads_per_worker"]),
    )
    names = tuple(REFERENCE_RGB)
    output_rows = []
    exact_pass = 0
    acquisition_fail = 0
    sha_fail = 0
    roi_fail = 0
    reproduction_fail = 0

    for position, row in enumerate(subset.itertuples(index=False), start=1):
        measurement = str(row.measurement_id)
        result = {
            "measurement_id": measurement,
            "recovery_status": "",
            "exact_image_sha_match": False,
            "roi_admission_reproduced": False,
            "flower_mask_pixels_exact": False,
            "flower_palette_exact": False,
            "background_pixels_exact": False,
            "background_effective_pixels": 0,
            **{f"background_palette_count_{name}": 0 for name in names},
        }
        if str(row.acquisition_status) != "acquired_and_decode_verified":
            result["recovery_status"] = "acquisition_failed_no_replacement"
            acquisition_fail += 1
            output_rows.append(result)
            continue
        image_path = args.images_dir / str(row.image_filename)
        if not image_path.is_file():
            raise RuntimeError(f"acquired image file is missing for {measurement}")
        actual_sha = hashlib.sha256(image_path.read_bytes()).hexdigest()
        expected_sha = str(row.expected_image_sha256)
        receipt_sha = str(row.image_sha256)
        if actual_sha != expected_sha or receipt_sha != expected_sha:
            result["recovery_status"] = "image_sha_mismatch_no_replacement"
            sha_fail += 1
            output_rows.append(result)
            continue
        result["exact_image_sha_match"] = True
        try:
            with Image.open(image_path) as source:
                measured = estimator.measure(source)
                oriented = ImageOps.exif_transpose(source).convert("RGB")
                rgb = np.asarray(oriented, dtype=np.uint8)
        except Exception:
            result["recovery_status"] = "roi_runtime_reproduction_failure"
            roi_fail += 1
            output_rows.append(result)
            continue
        admitted = str(measured.get("automated_colour_state_status") or "") == "automated_colour_state_admitted"
        result["roi_admission_reproduced"] = admitted
        if not admitted:
            result["recovery_status"] = "roi_admission_not_reproduced"
            roi_fail += 1
            output_rows.append(result)
            continue
        flower_mask = np.asarray(measured["flower_mask"], dtype=bool)
        background_mask = np.asarray(measured["background_mask"], dtype=bool)
        if flower_mask.shape != rgb.shape[:2] or background_mask.shape != rgb.shape[:2]:
            raise RuntimeError("recovered ROI masks and oriented RGB have incompatible shapes")
        flower_counts = nearest_palette_counts(rgb[flower_mask])
        background_counts = nearest_palette_counts(rgb[background_mask])
        flower_pixels = int(np.count_nonzero(flower_mask))
        background_pixels = int(np.count_nonzero(background_mask))
        expected_flower_pixels = _expected_integer(
            row.expected_flower_mask_pixels, label="flower mask pixels"
        )
        expected_background_pixels = _expected_integer(
            row.expected_background_effective_pixels, label="background pixels"
        )
        result["flower_mask_pixels_exact"] = flower_pixels == expected_flower_pixels
        result["background_pixels_exact"] = background_pixels == expected_background_pixels
        palette_exact = all(
            int(flower_counts[name])
            == _expected_integer(
                getattr(row, f"expected_flower_palette_count_{name}"),
                label=f"flower palette count {name}",
            )
            for name in names
        )
        result["flower_palette_exact"] = palette_exact
        result["background_effective_pixels"] = background_pixels
        if not (
            result["flower_mask_pixels_exact"]
            and result["background_pixels_exact"]
            and palette_exact
        ):
            result["recovery_status"] = "original_roi_or_flower_palette_not_exactly_reproduced"
            reproduction_fail += 1
            output_rows.append(result)
            continue
        if sum(background_counts.values()) != background_pixels or background_pixels <= 0:
            raise RuntimeError("recovered background palette census is internally inconsistent")
        result.update({f"background_palette_count_{name}": int(background_counts[name]) for name in names})
        result["recovery_status"] = "exact_matched_background_recovered"
        exact_pass += 1
        output_rows.append(result)
        if position % 10 == 0 or position == len(subset):
            print(f"background_recovered_or_failed={position}/{len(subset)}", flush=True)

    frame = pd.DataFrame(output_rows)
    if len(frame) != len(selected) or frame["measurement_id"].nunique() != len(selected):
        raise RuntimeError("background recovery output does not exactly cover selected partition")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_csv, index=False, lineterminator="\n")
    manifest = {
        "status": "complete_background_control_partition_recovery",
        "selected_rows": int(len(frame)),
        "exact_matched_background_recovered": int(exact_pass),
        "acquisition_failed_no_replacement": int(acquisition_fail),
        "image_sha_mismatch_no_replacement": int(sha_fail),
        "roi_runtime_or_admission_failure": int(roi_fail),
        "original_roi_or_flower_palette_not_exactly_reproduced": int(reproduction_fail),
        "background_colour_pixels_opened": bool(exact_pass or reproduction_fail),
        "species_opened_to_worker": False,
        "coordinates_opened_to_worker": False,
        "replacement_photos_used": False,
    }
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
