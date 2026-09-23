#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageOps

# Load current v2 perturbation helpers without importing the current
# fcp_pipeline package. The biological classifier/ROI package is loaded later
# from the exact historical source root.
_HELPER_PATH = (
    Path(__file__).resolve().parents[2]
    / "fcp_pipeline"
    / "measurement_validity_v2.py"
)
_HELPER_SPEC = importlib.util.spec_from_file_location(
    "fcp_v2_measurement_helper_bio", _HELPER_PATH
)
if _HELPER_SPEC is None or _HELPER_SPEC.loader is None:
    raise RuntimeError("cannot load v2 measurement helper")
_HELPER = importlib.util.module_from_spec(_HELPER_SPEC)
sys.modules[_HELPER_SPEC.name] = _HELPER
_HELPER_SPEC.loader.exec_module(_HELPER)

EV_LEVELS = tuple(float(v) for v in _HELPER.EV_LEVELS)
apply_exposure_ev = _HELPER.apply_exposure_ev
frozen_prompt_jitter_set = _HELPER.frozen_prompt_jitter_set
neutralize_background = _HELPER.neutralize_background

HEAVY_EV_LEVELS = (-1.0, -0.5, 0.5, 1.0)
JITTER_LABELS = (
    "base",
    "x_m5",
    "x_p5",
    "y_m5",
    "y_p5",
    "scale_m10",
    "scale_p10",
)
BIOLOGICAL_PALETTE = (
    "white",
    "yellow",
    "orange",
    "red",
    "pink",
    "magenta",
    "purple",
    "blue",
    "bronze",
)


def _as_bool(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes"}


def _tag_ev(prefix: str, ev: float) -> str:
    value = f"{float(ev):+.1f}".replace("+", "p").replace("-", "m").replace(".", "_")
    return f"{prefix}_{value}"


def _empty_fractions() -> dict[str, float]:
    return {name: 0.0 for name in BIOLOGICAL_PALETTE}


def _row(
    *,
    measurement_id: str,
    condition_id: str,
    condition_family: str,
    condition_value: str,
    heavy_counterfactual: bool,
    morph: str,
    measurement_status: str,
    counterfactual_status: str,
    roi_status: str,
    failure_reasons: str,
    image_sha256: str,
    mask_pixels: int,
    flower_effective_pixels: int,
    background_effective_pixels: int,
    horizontal_flip_mask_iou: float,
    horizontal_flip_colour_delta_e: float,
    fractions: dict[str, float],
) -> dict[str, object]:
    out: dict[str, object] = {
        "measurement_id": measurement_id,
        "condition_id": condition_id,
        "condition_family": condition_family,
        "condition_value": condition_value,
        "heavy_counterfactual": bool(heavy_counterfactual),
        "morph": morph,
        "measurement_status": measurement_status,
        "counterfactual_status": counterfactual_status,
        "roi_status": roi_status,
        "failure_reasons": failure_reasons,
        "image_sha256": image_sha256,
        "mask_pixels": int(mask_pixels),
        "flower_effective_pixels": int(flower_effective_pixels),
        "background_effective_pixels": int(background_effective_pixels),
        "horizontal_flip_mask_iou": horizontal_flip_mask_iou,
        "horizontal_flip_colour_delta_e": horizontal_flip_colour_delta_e,
    }
    for name in BIOLOGICAL_PALETTE:
        out[f"flower_fraction_{name}"] = float(fractions.get(name, 0.0))
    return out


def _classify_mask(
    *,
    rgb: np.ndarray,
    mask: np.ndarray,
    classify_masked_rgb: object,
    minimum_pixels: int,
    minimum_dominant_fraction: float,
    minimum_margin: float,
) -> dict[str, object]:
    m = np.asarray(mask, dtype=bool)
    if m.shape != rgb.shape[:2]:
        raise RuntimeError("flower mask and RGB shape drift")
    classified = classify_masked_rgb(
        rgb[m],
        minimum_mask_pixels=minimum_pixels,
        minimum_dominant_fraction=minimum_dominant_fraction,
        minimum_margin=minimum_margin,
    )
    return {
        "morph": str(classified["morph"]),
        "measurement_status": str(classified["measurement_status"]),
        "mask_pixels": int(classified["mask_pixels"]),
        "fractions": {
            name: float(classified["flower_only_fractions"][name])
            for name in BIOLOGICAL_PALETTE
        },
    }


def _classify_full_pipeline(
    *,
    image: Image.Image,
    rgb: np.ndarray,
    estimator: object,
    classify_masked_rgb: object,
    minimum_pixels: int,
    minimum_dominant_fraction: float,
    minimum_margin: float,
) -> tuple[dict[str, object], dict[str, object]]:
    measured = estimator.measure(image)
    roi_status = str(measured.get("automated_colour_state_status") or "")
    flower_pixels = int(measured.get("flower_effective_pixels") or 0)
    background_pixels = int(measured.get("background_effective_pixels") or 0)
    flip_iou = (
        float(measured["horizontal_flip_mask_iou"])
        if measured.get("horizontal_flip_mask_iou") is not None
        else float("nan")
    )
    flip_delta = (
        float(measured["horizontal_flip_colour_delta_e"])
        if measured.get("horizontal_flip_colour_delta_e") is not None
        else float("nan")
    )
    roi_failures = str(measured.get("failure_reasons") or "")
    if roi_status != "automated_colour_state_admitted":
        classified = {
            "morph": "mixed_uncertain",
            "measurement_status": "not_evaluable_roi_or_flip_gate",
            "mask_pixels": flower_pixels,
            "fractions": _empty_fractions(),
        }
    else:
        mask = np.asarray(measured["flower_mask"], dtype=bool)
        classified = _classify_mask(
            rgb=rgb,
            mask=mask,
            classify_masked_rgb=classify_masked_rgb,
            minimum_pixels=minimum_pixels,
            minimum_dominant_fraction=minimum_dominant_fraction,
            minimum_margin=minimum_margin,
        )
    meta = {
        "roi_status": roi_status or "automated_colour_state_not_evaluable",
        "failure_reasons": (
            ""
            if classified["measurement_status"] == "classified_four_state_morph"
            else (
                roi_failures
                if classified["measurement_status"]
                == "not_evaluable_roi_or_flip_gate"
                else str(classified["measurement_status"])
            )
        ),
        "flower_effective_pixels": flower_pixels,
        "background_effective_pixels": background_pixels,
        "horizontal_flip_mask_iou": flip_iou,
        "horizontal_flip_colour_delta_e": flip_delta,
        "measurement": measured,
    }
    return classified, meta


def _jitter_masks(
    *,
    estimator: object,
    oriented: Image.Image,
    boxes: list[list[float]],
    box_to_canvas: object,
    select_prompt_mask: object,
    letterboxed_rgb: object,
    canvas_mask_to_original: object,
    canvas_size: int,
) -> list[np.ndarray | None]:
    if not boxes:
        return [None] * len(JITTER_LABELS)
    width, height = oriented.size
    canvas, geometry = letterboxed_rgb(oriented)
    embeddings = estimator.encoder.run(
        None,
        {
            "batched_images": (
                canvas.transpose(2, 0, 1)[None].astype(np.float32) / 255.0
            )
        },
    )[0]
    box_variants = [frozen_prompt_jitter_set(box) for box in boxes]
    outputs: list[np.ndarray | None] = []
    for variant_index in range(len(JITTER_LABELS)):
        union_canvas = np.zeros((canvas_size, canvas_size), dtype=bool)
        valid_prompt = False
        for variants in box_variants:
            jittered = variants[variant_index]
            try:
                canvas_box = box_to_canvas(
                    jittered,
                    width=width,
                    height=height,
                )
            except ValueError:
                continue
            points = np.array(
                [
                    [
                        [
                            [canvas_box[0], canvas_box[1]],
                            [canvas_box[2], canvas_box[3]],
                        ]
                    ]
                ],
                dtype=np.float32,
            )
            labels = np.array([[[2.0, 3.0]]], dtype=np.float32)
            output_masks, predicted_iou, _ = estimator.decoder.run(
                None,
                {
                    "image_embeddings": embeddings,
                    "batched_point_coords": points,
                    "batched_point_labels": labels,
                    "orig_im_size": np.array(
                        [canvas_size, canvas_size], dtype=np.int64
                    ),
                },
            )
            selected = select_prompt_mask(
                output_masks[0, 0],
                predicted_iou[0, 0],
                canvas_box,
            )
            if selected.any():
                valid_prompt = True
                union_canvas |= selected
        if not valid_prompt or not union_canvas.any():
            outputs.append(None)
            continue
        outputs.append(
            canvas_mask_to_original(
                union_canvas,
                geometry,
                width=width,
                height=height,
            )
        )
    return outputs


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--worker-manifest", type=Path, required=True)
    p.add_argument("--images-dir", type=Path, required=True)
    p.add_argument("--frozen-source-root", type=Path, required=True)
    p.add_argument("--efficient-sam-dir", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    args = p.parse_args()

    # Ensure the historical package is the one imported for ROI and biological
    # palette classification.
    for key in list(sys.modules):
        if key == "fcp_pipeline" or key.startswith("fcp_pipeline."):
            del sys.modules[key]
    sys.path.insert(0, str(args.frozen_source_root))

    from fcp_pipeline.flower_roi_v4 import (
        CANVAS_SIZE,
        box_to_canvas,
        select_prompt_mask,
        validate_roi_v4_contract,
    )
    from fcp_pipeline.flower_roi_v4_runtime import (
        FrozenFlowerColourEstimator,
        _canvas_mask_to_original,
        _letterboxed_rgb,
        file_sha256,
        validate_scaleout_authorization,
    )
    from fcp_pipeline.photo_first_measurement import classify_masked_rgb

    roi_path = (
        args.frozen_source_root
        / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json"
    )
    locked_path = (
        args.frozen_source_root
        / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json"
    )
    detector = (
        args.frozen_source_root
        / "data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    )
    measurement_contract_path = (
        args.frozen_source_root
        / "docs/supporting/random_photo_first_measurement_contract_v1.json"
    )

    roi_contract = json.loads(roi_path.read_text(encoding="utf-8"))
    measurement_contract = json.loads(
        measurement_contract_path.read_text(encoding="utf-8")
    )
    locked = json.loads(locked_path.read_text(encoding="utf-8"))
    validate_roi_v4_contract(roi_contract)
    detector_sha = file_sha256(detector)
    validate_scaleout_authorization(
        locked,
        trained_weight_sha256=detector_sha,
    )

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    expected_columns = [
        "measurement_id",
        "image_filename",
        "photo_license",
        "heavy_counterfactual",
        "expected_source_sha256",
    ]
    if list(worker.columns) != expected_columns:
        raise RuntimeError(
            f"Pass-B worker schema drift: {list(worker.columns)}"
        )
    if worker["measurement_id"].nunique() != len(worker):
        raise RuntimeError("Pass-B worker IDs are not unique")

    estimator = FrozenFlowerColourEstimator(
        detector,
        args.efficient_sam_dir,
        roi_contract,
        torch_threads=1,
    )
    threshold = measurement_contract["coarse_colour_state"]
    minimum_pixels = int(
        measurement_contract["flower_roi"]["minimum_union_flower_pixels"]
    )
    minimum_dominant_fraction = float(
        threshold["minimum_dominant_fraction"]
    )
    minimum_margin = float(
        threshold["minimum_margin_over_second_group"]
    )

    rows: list[dict[str, object]] = []
    for position, item in enumerate(worker.itertuples(index=False), start=1):
        measurement_id = str(item.measurement_id)
        heavy = _as_bool(item.heavy_counterfactual)
        expected_sha = str(item.expected_source_sha256).strip().lower()
        image_path = args.images_dir / str(item.image_filename)
        if not image_path.is_file():
            raise RuntimeError(
                f"Pass-B acquisition packet missing {item.image_filename}"
            )
        image_sha = hashlib.sha256(image_path.read_bytes()).hexdigest()
        if image_sha != expected_sha:
            raise RuntimeError(
                "source SHA drift reached biological measurement worker"
            )

        with Image.open(image_path) as source:
            oriented = ImageOps.exif_transpose(source).convert("RGB")
        rgb = np.asarray(oriented, dtype=np.uint8)

        baseline_error = ""
        try:
            baseline_class, baseline_meta = _classify_full_pipeline(
                image=oriented,
                rgb=rgb,
                estimator=estimator,
                classify_masked_rgb=classify_masked_rgb,
                minimum_pixels=minimum_pixels,
                minimum_dominant_fraction=minimum_dominant_fraction,
                minimum_margin=minimum_margin,
            )
        except Exception as exc:
            baseline_error = f"{type(exc).__name__}:{str(exc)[:350]}"
            baseline_class = {
                "morph": "mixed_uncertain",
                "measurement_status": "not_evaluable_roi_or_flip_gate",
                "mask_pixels": 0,
                "fractions": _empty_fractions(),
            }
            baseline_meta = {
                "roi_status": "roi_runtime_failure",
                "failure_reasons": baseline_error,
                "flower_effective_pixels": 0,
                "background_effective_pixels": 0,
                "horizontal_flip_mask_iou": float("nan"),
                "horizontal_flip_colour_delta_e": float("nan"),
                "measurement": None,
            }

        rows.append(
            _row(
                measurement_id=measurement_id,
                condition_id="baseline",
                condition_family="baseline",
                condition_value="baseline",
                heavy_counterfactual=heavy,
                morph=str(baseline_class["morph"]),
                measurement_status=str(
                    baseline_class["measurement_status"]
                ),
                counterfactual_status="available",
                roi_status=str(baseline_meta["roi_status"]),
                failure_reasons=str(
                    baseline_meta["failure_reasons"]
                ),
                image_sha256=image_sha,
                mask_pixels=int(baseline_class["mask_pixels"]),
                flower_effective_pixels=int(
                    baseline_meta["flower_effective_pixels"]
                ),
                background_effective_pixels=int(
                    baseline_meta["background_effective_pixels"]
                ),
                horizontal_flip_mask_iou=float(
                    baseline_meta["horizontal_flip_mask_iou"]
                ),
                horizontal_flip_colour_delta_e=float(
                    baseline_meta["horizontal_flip_colour_delta_e"]
                ),
                fractions=dict(baseline_class["fractions"]),
            )
        )

        measured = baseline_meta.get("measurement")
        baseline_mask = (
            np.asarray(measured["flower_mask"], dtype=bool)
            if isinstance(measured, dict)
            and measured.get("flower_mask") is not None
            else None
        )
        baseline_roi_admitted = (
            str(baseline_meta["roi_status"])
            == "automated_colour_state_admitted"
            and baseline_mask is not None
            and baseline_mask.any()
        )

        # All-row fixed-mask exposure sensitivity.
        for ev in EV_LEVELS:
            condition_id = _tag_ev("fixed_ev", ev)
            if baseline_roi_admitted:
                ev_rgb = apply_exposure_ev(rgb, ev)
                classified = _classify_mask(
                    rgb=ev_rgb,
                    mask=baseline_mask,
                    classify_masked_rgb=classify_masked_rgb,
                    minimum_pixels=minimum_pixels,
                    minimum_dominant_fraction=minimum_dominant_fraction,
                    minimum_margin=minimum_margin,
                )
                cf_status = "available"
                failure = (
                    ""
                    if classified["measurement_status"]
                    == "classified_four_state_morph"
                    else str(classified["measurement_status"])
                )
            else:
                classified = {
                    "morph": "mixed_uncertain",
                    "measurement_status": "not_evaluable_roi_or_flip_gate",
                    "mask_pixels": int(
                        baseline_mask.sum()
                        if baseline_mask is not None
                        else 0
                    ),
                    "fractions": _empty_fractions(),
                }
                cf_status = "unavailable_baseline_roi"
                failure = (
                    baseline_error
                    or str(baseline_meta["failure_reasons"])
                    or "baseline_roi_not_admitted"
                )
            rows.append(
                _row(
                    measurement_id=measurement_id,
                    condition_id=condition_id,
                    condition_family="fixed_mask_exposure",
                    condition_value=f"{ev:+.1f}",
                    heavy_counterfactual=heavy,
                    morph=str(classified["morph"]),
                    measurement_status=str(
                        classified["measurement_status"]
                    ),
                    counterfactual_status=cf_status,
                    roi_status="baseline_mask_fixed",
                    failure_reasons=failure,
                    image_sha256=image_sha,
                    mask_pixels=int(classified["mask_pixels"]),
                    flower_effective_pixels=int(
                        baseline_meta["flower_effective_pixels"]
                    ),
                    background_effective_pixels=int(
                        baseline_meta["background_effective_pixels"]
                    ),
                    horizontal_flip_mask_iou=float(
                        baseline_meta["horizontal_flip_mask_iou"]
                    ),
                    horizontal_flip_colour_delta_e=float(
                        baseline_meta["horizontal_flip_colour_delta_e"]
                    ),
                    fractions=dict(classified["fractions"]),
                )
            )

        if heavy:
            # Full-pipeline exposure perturbations.
            for ev in HEAVY_EV_LEVELS:
                condition_id = _tag_ev("full_ev", ev)
                try:
                    ev_rgb = apply_exposure_ev(rgb, ev)
                    classified, meta = _classify_full_pipeline(
                        image=Image.fromarray(ev_rgb, mode="RGB"),
                        rgb=ev_rgb,
                        estimator=estimator,
                        classify_masked_rgb=classify_masked_rgb,
                        minimum_pixels=minimum_pixels,
                        minimum_dominant_fraction=minimum_dominant_fraction,
                        minimum_margin=minimum_margin,
                    )
                    cf_status = "available"
                    failure = str(meta["failure_reasons"])
                except Exception as exc:
                    classified = {
                        "morph": "mixed_uncertain",
                        "measurement_status": "not_evaluable_roi_or_flip_gate",
                        "mask_pixels": 0,
                        "fractions": _empty_fractions(),
                    }
                    meta = {
                        "roi_status": "roi_runtime_failure",
                        "flower_effective_pixels": 0,
                        "background_effective_pixels": 0,
                        "horizontal_flip_mask_iou": float("nan"),
                        "horizontal_flip_colour_delta_e": float("nan"),
                    }
                    cf_status = "runtime_failure"
                    failure = f"{type(exc).__name__}:{str(exc)[:350]}"
                rows.append(
                    _row(
                        measurement_id=measurement_id,
                        condition_id=condition_id,
                        condition_family="full_pipeline_exposure",
                        condition_value=f"{ev:+.1f}",
                        heavy_counterfactual=True,
                        morph=str(classified["morph"]),
                        measurement_status=str(
                            classified["measurement_status"]
                        ),
                        counterfactual_status=cf_status,
                        roi_status=str(meta["roi_status"]),
                        failure_reasons=failure,
                        image_sha256=image_sha,
                        mask_pixels=int(classified["mask_pixels"]),
                        flower_effective_pixels=int(
                            meta["flower_effective_pixels"]
                        ),
                        background_effective_pixels=int(
                            meta["background_effective_pixels"]
                        ),
                        horizontal_flip_mask_iou=float(
                            meta["horizontal_flip_mask_iou"]
                        ),
                        horizontal_flip_colour_delta_e=float(
                            meta["horizontal_flip_colour_delta_e"]
                        ),
                        fractions=dict(classified["fractions"]),
                    )
                )

            # Background neutralization.
            if baseline_mask is not None and baseline_mask.any():
                try:
                    neutral_rgb = neutralize_background(
                        rgb,
                        baseline_mask,
                    )
                    classified, meta = _classify_full_pipeline(
                        image=Image.fromarray(neutral_rgb, mode="RGB"),
                        rgb=neutral_rgb,
                        estimator=estimator,
                        classify_masked_rgb=classify_masked_rgb,
                        minimum_pixels=minimum_pixels,
                        minimum_dominant_fraction=minimum_dominant_fraction,
                        minimum_margin=minimum_margin,
                    )
                    cf_status = "available"
                    failure = str(meta["failure_reasons"])
                except Exception as exc:
                    classified = {
                        "morph": "mixed_uncertain",
                        "measurement_status": "not_evaluable_roi_or_flip_gate",
                        "mask_pixels": 0,
                        "fractions": _empty_fractions(),
                    }
                    meta = {
                        "roi_status": "roi_runtime_failure",
                        "flower_effective_pixels": 0,
                        "background_effective_pixels": 0,
                        "horizontal_flip_mask_iou": float("nan"),
                        "horizontal_flip_colour_delta_e": float("nan"),
                    }
                    cf_status = "runtime_failure"
                    failure = f"{type(exc).__name__}:{str(exc)[:350]}"
            else:
                classified = {
                    "morph": "mixed_uncertain",
                    "measurement_status": "not_evaluable_roi_or_flip_gate",
                    "mask_pixels": 0,
                    "fractions": _empty_fractions(),
                }
                meta = {
                    "roi_status": "baseline_mask_unavailable",
                    "flower_effective_pixels": 0,
                    "background_effective_pixels": 0,
                    "horizontal_flip_mask_iou": float("nan"),
                    "horizontal_flip_colour_delta_e": float("nan"),
                }
                cf_status = "unavailable_baseline_mask"
                failure = "baseline_flower_mask_unavailable"

            rows.append(
                _row(
                    measurement_id=measurement_id,
                    condition_id="neutral_bg",
                    condition_family="background_neutralization",
                    condition_value="srgb_128_128_128",
                    heavy_counterfactual=True,
                    morph=str(classified["morph"]),
                    measurement_status=str(
                        classified["measurement_status"]
                    ),
                    counterfactual_status=cf_status,
                    roi_status=str(meta["roi_status"]),
                    failure_reasons=failure,
                    image_sha256=image_sha,
                    mask_pixels=int(classified["mask_pixels"]),
                    flower_effective_pixels=int(
                        meta["flower_effective_pixels"]
                    ),
                    background_effective_pixels=int(
                        meta["background_effective_pixels"]
                    ),
                    horizontal_flip_mask_iou=float(
                        meta["horizontal_flip_mask_iou"]
                    ),
                    horizontal_flip_colour_delta_e=float(
                        meta["horizontal_flip_colour_delta_e"]
                    ),
                    fractions=dict(classified["fractions"]),
                )
            )

            # Fixed-detector-box prompt-jitter masks.
            boxes = (
                measured.get("boxes") or []
                if isinstance(measured, dict)
                else []
            )
            try:
                jitter_masks = _jitter_masks(
                    estimator=estimator,
                    oriented=oriented,
                    boxes=boxes,
                    box_to_canvas=box_to_canvas,
                    select_prompt_mask=select_prompt_mask,
                    letterboxed_rgb=_letterboxed_rgb,
                    canvas_mask_to_original=_canvas_mask_to_original,
                    canvas_size=CANVAS_SIZE,
                )
            except Exception:
                jitter_masks = [None] * len(JITTER_LABELS)

            for label, mask in zip(
                JITTER_LABELS,
                jitter_masks,
                strict=True,
            ):
                if mask is None or not np.asarray(mask, dtype=bool).any():
                    classified = {
                        "morph": "mixed_uncertain",
                        "measurement_status": "not_evaluable_roi_or_flip_gate",
                        "mask_pixels": 0,
                        "fractions": _empty_fractions(),
                    }
                    cf_status = "unavailable_prompt_mask"
                    failure = "prompt_mask_unavailable"
                else:
                    classified = _classify_mask(
                        rgb=rgb,
                        mask=np.asarray(mask, dtype=bool),
                        classify_masked_rgb=classify_masked_rgb,
                        minimum_pixels=minimum_pixels,
                        minimum_dominant_fraction=minimum_dominant_fraction,
                        minimum_margin=minimum_margin,
                    )
                    cf_status = "available"
                    failure = (
                        ""
                        if classified["measurement_status"]
                        == "classified_four_state_morph"
                        else str(classified["measurement_status"])
                    )
                rows.append(
                    _row(
                        measurement_id=measurement_id,
                        condition_id=f"jitter_{label}",
                        condition_family="roi_prompt_jitter",
                        condition_value=label,
                        heavy_counterfactual=True,
                        morph=str(classified["morph"]),
                        measurement_status=str(
                            classified["measurement_status"]
                        ),
                        counterfactual_status=cf_status,
                        roi_status="fixed_detector_box_prompt_jitter",
                        failure_reasons=failure,
                        image_sha256=image_sha,
                        mask_pixels=int(classified["mask_pixels"]),
                        flower_effective_pixels=int(
                            classified["mask_pixels"]
                        ),
                        background_effective_pixels=0,
                        horizontal_flip_mask_iou=float("nan"),
                        horizontal_flip_colour_delta_e=float("nan"),
                        fractions=dict(classified["fractions"]),
                    )
                )

        if position % 5 == 0 or position == len(worker):
            print(
                f"pass_b_biological_measured={position}/{len(worker)}",
                flush=True,
            )

    frame = pd.DataFrame(rows)
    if frame.empty and len(worker):
        raise RuntimeError("Pass-B biological worker produced no rows")
    if len(frame):
        if frame.duplicated(["measurement_id", "condition_id"]).any():
            raise RuntimeError("duplicate biological condition rows")
        baseline = frame.loc[frame["condition_id"].eq("baseline")]
        if baseline["measurement_id"].nunique() != len(worker):
            raise RuntimeError("Pass-B baseline coverage is incomplete")
        fixed = frame.loc[
            frame["condition_family"].eq("fixed_mask_exposure")
        ]
        if len(fixed) != len(worker) * len(EV_LEVELS):
            raise RuntimeError("fixed-mask exposure coverage is incomplete")
        heavy_n = int(
            worker["heavy_counterfactual"]
            .map(_as_bool)
            .sum()
        )
        heavy_expected = heavy_n * (
            len(HEAVY_EV_LEVELS) + 1 + len(JITTER_LABELS)
        )
        heavy_rows = frame.loc[
            frame["condition_family"].isin(
                {
                    "full_pipeline_exposure",
                    "background_neutralization",
                    "roi_prompt_jitter",
                }
            )
        ]
        if len(heavy_rows) != heavy_expected:
            raise RuntimeError(
                f"heavy biological coverage drift: {len(heavy_rows)} != {heavy_expected}"
            )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_csv, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
