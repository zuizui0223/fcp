#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageOps
from skimage.color import rgb2lab

from fcp_pipeline.measurement_validity_v2 import (
    EV_LEVELS,
    apply_exposure_ev,
    exposure_metrics,
    frozen_prompt_jitter_set,
    mask_boundary_fraction,
    mask_iou,
    neutralize_background,
    relative_luminance,
)

HEAVY_EV_LEVELS = (-1.0, -0.5, 0.5, 1.0)
HIGH_BACKGROUND_LUMINANCE = 0.90


def _as_bool(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes"}


def _lab_mean(rgb: np.ndarray, mask: np.ndarray) -> np.ndarray | None:
    x = np.asarray(rgb, dtype=np.uint8)
    m = np.asarray(mask, dtype=bool)
    if x.ndim != 3 or x.shape[2] != 3 or m.shape != x.shape[:2] or not m.any():
        return None
    lab = rgb2lab(x.astype(np.float32) / 255.0)
    out = np.mean(lab[m].astype(np.float64), axis=0)
    return out if np.all(np.isfinite(out)) else None


def _lab_delta(a: np.ndarray | None, b: np.ndarray | None) -> float:
    if a is None or b is None:
        return float("nan")
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def _prefix_metrics(prefix: str, values: dict[str, float]) -> dict[str, float]:
    return {f"{prefix}_{key}": float(value) for key, value in values.items()}


def _mask_exposure_metrics(rgb: np.ndarray, mask: np.ndarray, prefix: str) -> dict[str, float]:
    m = np.asarray(mask, dtype=bool)
    if m.shape != rgb.shape[:2] or not m.any():
        return {
            f"{prefix}_{name}": float("nan")
            for name in (
                "clip_fraction",
                "near_clip_fraction",
                "clip_fraction_r",
                "clip_fraction_g",
                "clip_fraction_b",
                "luminance_q01",
                "luminance_q10",
                "luminance_q50",
                "luminance_q90",
                "luminance_q99",
                "dynamic_range_q99_minus_q01",
                "black_fraction",
            )
        }
    return _prefix_metrics(prefix, exposure_metrics(rgb[m]))


def _flower_background_lab_distance(rgb: np.ndarray, flower: np.ndarray, background: np.ndarray) -> float:
    return _lab_delta(_lab_mean(rgb, flower), _lab_mean(rgb, background))


def _boundary_adjacent_high_background_fraction(
    rgb: np.ndarray,
    flower_mask: np.ndarray,
    background_mask: np.ndarray,
) -> float:
    flower = np.asarray(flower_mask, dtype=bool)
    background = np.asarray(background_mask, dtype=bool)
    if flower.shape != background.shape or flower.shape != rgb.shape[:2]:
        return float("nan")
    if not flower.any() or not background.any():
        return float("nan")

    # 8-neighbour dilation without scipy; only background pixels adjacent to the
    # flower boundary are eligible.
    p = np.pad(flower, 1, constant_values=False)
    adjacent = np.zeros_like(flower, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0:
                continue
            adjacent |= p[1 + dy : 1 + dy + flower.shape[0], 1 + dx : 1 + dx + flower.shape[1]]
    adjacent &= ~flower
    adjacent &= background
    n = int(adjacent.sum())
    if n == 0:
        return float("nan")
    lum = relative_luminance(rgb)
    return float(np.mean(lum[adjacent] >= HIGH_BACKGROUND_LUMINANCE))


def _technical_summary_from_measurement(
    rgb: np.ndarray,
    measured: dict[str, object],
    prefix: str,
    baseline_mask: np.ndarray | None = None,
    baseline_lab: np.ndarray | None = None,
) -> dict[str, object]:
    flower = np.asarray(measured["flower_mask"], dtype=bool)
    background = np.asarray(measured["background_mask"], dtype=bool)
    out: dict[str, object] = {
        f"{prefix}_roi_status": str(measured.get("automated_colour_state_status") or ""),
        f"{prefix}_failure_reasons": str(measured.get("failure_reasons") or ""),
        f"{prefix}_detector_box_count": int(len(measured.get("boxes") or [])),
        f"{prefix}_retained_instances": int(measured.get("retained_instances") or 0),
        f"{prefix}_flower_effective_pixels": int(flower.sum()),
        f"{prefix}_background_effective_pixels": int(background.sum()),
        f"{prefix}_mask_area_fraction": float(flower.mean()),
        f"{prefix}_mask_boundary_fraction": float(mask_boundary_fraction(flower)) if flower.any() else float("nan"),
        f"{prefix}_horizontal_flip_mask_iou": float(measured.get("horizontal_flip_mask_iou") or 0.0),
        f"{prefix}_horizontal_flip_colour_delta_e": (
            float(measured["horizontal_flip_colour_delta_e"])
            if measured.get("horizontal_flip_colour_delta_e") is not None
            else float("nan")
        ),
        f"{prefix}_flower_background_lab_distance": _flower_background_lab_distance(
            rgb, flower, background
        ),
        f"{prefix}_boundary_adjacent_high_background_fraction": _boundary_adjacent_high_background_fraction(
            rgb, flower, background
        ),
    }
    out.update(_mask_exposure_metrics(rgb, flower, f"{prefix}_flower"))
    out.update(_mask_exposure_metrics(rgb, background, f"{prefix}_background"))

    for channel in ("L", "a", "b"):
        for metric in ("mean", "sd", "q10", "q50", "q90"):
            key = f"flower_{channel}_{metric}"
            out[f"{prefix}_{key}"] = (
                float(measured[key]) if measured.get(key) is not None else float("nan")
            )
            bkey = f"background_{channel}_{metric}"
            out[f"{prefix}_{bkey}"] = (
                float(measured[bkey]) if measured.get(bkey) is not None else float("nan")
            )

    if baseline_mask is not None:
        out[f"{prefix}_mask_iou_to_baseline"] = mask_iou(flower, baseline_mask)
    if baseline_lab is not None:
        out[f"{prefix}_flower_lab_delta_e_to_baseline"] = _lab_delta(
            _lab_mean(rgb, flower), baseline_lab
        )
    return out


def _prompt_jitter_summary(
    estimator: object,
    oriented: Image.Image,
    rgb: np.ndarray,
    boxes: list[list[float]],
    baseline_mask: np.ndarray,
    baseline_lab: np.ndarray | None,
    *,
    box_to_canvas: object,
    select_prompt_mask: object,
    letterboxed_rgb: object,
    canvas_mask_to_original: object,
    canvas_size: int,
) -> dict[str, object]:
    if not boxes:
        return {
            "jitter_valid_variants": 0,
            "jitter_mask_iou_min": float("nan"),
            "jitter_mask_iou_median": float("nan"),
            "jitter_colour_delta_e_max": float("nan"),
            "jitter_colour_delta_e_median": float("nan"),
        }

    width, height = oriented.size
    canvas, geometry = letterboxed_rgb(oriented)
    embeddings = estimator.encoder.run(
        None,
        {"batched_images": canvas.transpose(2, 0, 1)[None].astype(np.float32) / 255.0},
    )[0]

    variant_boxes = [frozen_prompt_jitter_set(box) for box in boxes]
    ious: list[float] = []
    deltas: list[float] = []
    for variant_index in range(7):
        union_canvas = np.zeros((canvas_size, canvas_size), dtype=bool)
        for box_set in variant_boxes:
            jittered = box_set[variant_index]
            try:
                canvas_box = box_to_canvas(jittered, width=width, height=height)
            except ValueError:
                continue
            points = np.array(
                [[[[canvas_box[0], canvas_box[1]], [canvas_box[2], canvas_box[3]]]]],
                dtype=np.float32,
            )
            labels = np.array([[[2.0, 3.0]]], dtype=np.float32)
            output_masks, predicted_iou, _ = estimator.decoder.run(
                None,
                {
                    "image_embeddings": embeddings,
                    "batched_point_coords": points,
                    "batched_point_labels": labels,
                    "orig_im_size": np.array([canvas_size, canvas_size], dtype=np.int64),
                },
            )
            selected = select_prompt_mask(
                output_masks[0, 0], predicted_iou[0, 0], canvas_box
            )
            union_canvas |= selected

        restored = canvas_mask_to_original(
            union_canvas, geometry, width=width, height=height
        )
        if not restored.any():
            continue
        ious.append(mask_iou(restored, baseline_mask))
        deltas.append(_lab_delta(_lab_mean(rgb, restored), baseline_lab))

    finite_deltas = [v for v in deltas if math.isfinite(v)]
    return {
        "jitter_valid_variants": int(len(ious)),
        "jitter_mask_iou_min": float(np.min(ious)) if ious else float("nan"),
        "jitter_mask_iou_median": float(np.median(ious)) if ious else float("nan"),
        "jitter_colour_delta_e_max": (
            float(np.max(finite_deltas)) if finite_deltas else float("nan")
        ),
        "jitter_colour_delta_e_median": (
            float(np.median(finite_deltas)) if finite_deltas else float("nan")
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--worker-manifest", type=Path, required=True)
    p.add_argument("--images-dir", type=Path, required=True)
    p.add_argument("--frozen-source-root", type=Path, required=True)
    p.add_argument("--efficient-sam-dir", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    args = p.parse_args()

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
    roi = json.loads(roi_path.read_text(encoding="utf-8"))
    locked = json.loads(locked_path.read_text(encoding="utf-8"))
    validate_roi_v4_contract(roi)
    detector_sha = file_sha256(detector)
    validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    expected = [
        "measurement_id",
        "image_filename",
        "photo_license",
        "heavy_counterfactual",
    ]
    if list(worker.columns) != expected:
        raise RuntimeError(f"blind worker schema drift: {list(worker.columns)}")

    estimator = FrozenFlowerColourEstimator(
        detector,
        args.efficient_sam_dir,
        roi,
        torch_threads=1,
    )

    rows: list[dict[str, object]] = []
    for i, row in enumerate(worker.itertuples(index=False), start=1):
        path = args.images_dir / str(row.image_filename)
        heavy = _as_bool(row.heavy_counterfactual)
        rec: dict[str, object] = {
            "measurement_id": str(row.measurement_id),
            "source_image_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "heavy_counterfactual": heavy,
            "technical_status": "",
            "technical_failure_reason": "",
        }
        try:
            with Image.open(path) as source:
                oriented = ImageOps.exif_transpose(source).convert("RGB")
            rgb = np.asarray(oriented, dtype=np.uint8)
            measured = estimator.measure(oriented)
            baseline_mask = np.asarray(measured["flower_mask"], dtype=bool)
            baseline_lab = _lab_mean(rgb, baseline_mask)
            if not baseline_mask.any():
                raise RuntimeError("baseline_flower_mask_unavailable")

            rec.update(
                _technical_summary_from_measurement(
                    rgb,
                    measured,
                    "base",
                )
            )

            # All-row fixed-mask photometric diagnostics. No biological palette
            # or morph classification is computed in this pass.
            for ev in EV_LEVELS:
                ev_rgb = apply_exposure_ev(rgb, ev)
                tag = f"fixed_ev_{ev:+.1f}".replace("+", "p").replace("-", "m").replace(".", "_")
                rec.update(
                    _mask_exposure_metrics(
                        ev_rgb,
                        baseline_mask,
                        f"{tag}_flower",
                    )
                )
                rec[f"{tag}_flower_lab_delta_e_to_base"] = _lab_delta(
                    _lab_mean(ev_rgb, baseline_mask),
                    baseline_lab,
                )

            if heavy:
                # Full-pipeline exposure counterfactuals.
                for ev in HEAVY_EV_LEVELS:
                    ev_rgb = apply_exposure_ev(rgb, ev)
                    ev_image = Image.fromarray(ev_rgb, mode="RGB")
                    perturbed = estimator.measure(ev_image)
                    tag = f"full_ev_{ev:+.1f}".replace("+", "p").replace("-", "m").replace(".", "_")
                    rec.update(
                        _technical_summary_from_measurement(
                            ev_rgb,
                            perturbed,
                            tag,
                            baseline_mask=baseline_mask,
                            baseline_lab=baseline_lab,
                        )
                    )

                # Full-pipeline background-neutralization counterfactual.
                neutral_rgb = neutralize_background(rgb, baseline_mask)
                neutral = estimator.measure(Image.fromarray(neutral_rgb, mode="RGB"))
                rec.update(
                    _technical_summary_from_measurement(
                        neutral_rgb,
                        neutral,
                        "neutral_bg",
                        baseline_mask=baseline_mask,
                        baseline_lab=baseline_lab,
                    )
                )

                # Segmentation prompt sensitivity with the detector boxes fixed.
                rec.update(
                    _prompt_jitter_summary(
                        estimator,
                        oriented,
                        rgb,
                        measured.get("boxes") or [],
                        baseline_mask,
                        baseline_lab,
                        box_to_canvas=box_to_canvas,
                        select_prompt_mask=select_prompt_mask,
                        letterboxed_rgb=_letterboxed_rgb,
                        canvas_mask_to_original=_canvas_mask_to_original,
                        canvas_size=CANVAS_SIZE,
                    )
                )

            rec["technical_status"] = "technical_metrics_available"
        except Exception as exc:
            rec["technical_status"] = "technical_metrics_unavailable"
            rec["technical_failure_reason"] = f"{type(exc).__name__}:{str(exc)[:400]}"
        rows.append(rec)
        if i % 5 == 0 or i == len(worker):
            print(f"technical_measured={i}/{len(worker)}", flush=True)

    out = pd.DataFrame(rows)
    if out["measurement_id"].nunique() != len(out):
        raise RuntimeError("duplicate technical IDs")
    forbidden = (
        "species",
        "taxon",
        "photo_id",
        "observation_id",
        "latitude",
        "longitude",
        "observer",
        "morph",
        "q_white",
        "spatial_outcome",
    )
    leaked = [
        c for c in out.columns
        if any(token in c.casefold() for token in forbidden)
    ]
    if leaked:
        raise RuntimeError(f"technical output leaked protected columns: {leaked}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
