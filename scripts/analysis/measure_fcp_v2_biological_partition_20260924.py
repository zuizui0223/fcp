#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageOps

ROOT = Path(__file__).resolve().parents[2]
_HELPER_PATH = ROOT / "fcp_pipeline" / "measurement_validity_v2.py"
_SPEC = importlib.util.spec_from_file_location("fcp_v2_measurement_helper", _HELPER_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError("cannot load v2 technical helper")
_HELPER = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _HELPER
_SPEC.loader.exec_module(_HELPER)

apply_exposure_ev = _HELPER.apply_exposure_ev
frozen_prompt_jitter_set = _HELPER.frozen_prompt_jitter_set
neutralize_background = _HELPER.neutralize_background

BIO = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
ZERO_FRACTIONS = {name: 0.0 for name in BIO}
HEAVY_EV_LEVELS = (-1.0, -0.5, 0.5, 1.0)


def _as_bool(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes"}


def _tag_ev(prefix: str, ev: float) -> str:
    return f"{prefix}_{ev:+.1f}".replace("+", "p").replace("-", "m").replace(".", "_")


def _empty(prefix: str, status: str) -> dict[str, object]:
    out: dict[str, object] = {
        f"{prefix}_morph": "mixed_uncertain",
        f"{prefix}_status": status,
        f"{prefix}_mask_pixels": 0,
    }
    out.update({f"{prefix}_fraction_{name}": 0.0 for name in BIO})
    return out


def _classify_mask(
    rgb: np.ndarray,
    mask: np.ndarray,
    *,
    prefix: str,
    classify_masked_rgb: object,
    minimum_pixels: int,
    minimum_dominant_fraction: float,
    minimum_margin: float,
) -> dict[str, object]:
    m = np.asarray(mask, dtype=bool)
    if m.shape != rgb.shape[:2] or not m.any():
        return _empty(prefix, "not_evaluable_roi_or_flip_gate")
    classified = classify_masked_rgb(
        rgb[m],
        minimum_mask_pixels=minimum_pixels,
        minimum_dominant_fraction=minimum_dominant_fraction,
        minimum_margin=minimum_margin,
    )
    out: dict[str, object] = {
        f"{prefix}_morph": str(classified["morph"]),
        f"{prefix}_status": str(classified["measurement_status"]),
        f"{prefix}_mask_pixels": int(classified["mask_pixels"]),
    }
    fractions = classified["flower_only_fractions"]
    out.update({f"{prefix}_fraction_{name}": float(fractions[name]) for name in BIO})
    return out


def _classify_pipeline(
    rgb: np.ndarray,
    measured: dict[str, object],
    *,
    prefix: str,
    classify_masked_rgb: object,
    minimum_pixels: int,
    minimum_dominant_fraction: float,
    minimum_margin: float,
) -> tuple[dict[str, object], np.ndarray | None]:
    roi_status = str(measured.get("automated_colour_state_status") or "")
    mask = np.asarray(measured.get("flower_mask"), dtype=bool) if measured.get("flower_mask") is not None else None
    if roi_status != "automated_colour_state_admitted" or mask is None:
        return _empty(prefix, "not_evaluable_roi_or_flip_gate"), None
    return (
        _classify_mask(
            rgb,
            mask,
            prefix=prefix,
            classify_masked_rgb=classify_masked_rgb,
            minimum_pixels=minimum_pixels,
            minimum_dominant_fraction=minimum_dominant_fraction,
            minimum_margin=minimum_margin,
        ),
        mask,
    )


def _prompt_variants(
    estimator: object,
    oriented: Image.Image,
    rgb: np.ndarray,
    boxes: list[list[float]],
    *,
    classify_masked_rgb: object,
    minimum_pixels: int,
    minimum_dominant_fraction: float,
    minimum_margin: float,
    box_to_canvas: object,
    select_prompt_mask: object,
    letterboxed_rgb: object,
    canvas_mask_to_original: object,
    canvas_size: int,
) -> dict[str, object]:
    out: dict[str, object] = {}
    if not boxes:
        for i in range(7):
            out.update(_empty(f"jitter_{i}", "not_evaluable_no_detector_box"))
        return out

    width, height = oriented.size
    canvas, geometry = letterboxed_rgb(oriented)
    embeddings = estimator.encoder.run(
        None,
        {"batched_images": canvas.transpose(2, 0, 1)[None].astype(np.float32) / 255.0},
    )[0]
    variant_boxes = [frozen_prompt_jitter_set(box) for box in boxes]

    for variant_index in range(7):
        union_canvas = np.zeros((canvas_size, canvas_size), dtype=bool)
        valid_prompt = False
        for box_set in variant_boxes:
            try:
                canvas_box = box_to_canvas(
                    box_set[variant_index],
                    width=width,
                    height=height,
                )
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
                output_masks[0, 0],
                predicted_iou[0, 0],
                canvas_box,
            )
            union_canvas |= selected
            valid_prompt = True
        prefix = f"jitter_{variant_index}"
        if not valid_prompt or not union_canvas.any():
            out.update(_empty(prefix, "not_evaluable_prompt_mask"))
            continue
        restored = canvas_mask_to_original(
            union_canvas,
            geometry,
            width=width,
            height=height,
        )
        out.update(
            _classify_mask(
                rgb,
                restored,
                prefix=prefix,
                classify_masked_rgb=classify_masked_rgb,
                minimum_pixels=minimum_pixels,
                minimum_dominant_fraction=minimum_dominant_fraction,
                minimum_margin=minimum_margin,
            )
        )
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--worker-manifest", type=Path, required=True)
    p.add_argument("--acquisition-receipt", type=Path, required=True)
    p.add_argument("--expected-sha", type=Path, required=True)
    p.add_argument("--images-dir", type=Path, required=True)
    p.add_argument("--frozen-source-root", type=Path, required=True)
    p.add_argument("--efficient-sam-dir", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    args = p.parse_args()

    sys.path.insert(0, str(args.frozen_source_root))
    import fcp_pipeline
    from fcp_pipeline.flower_roi_v4 import CANVAS_SIZE, box_to_canvas, select_prompt_mask
    from fcp_pipeline.flower_roi_v4_runtime import (
        FrozenFlowerColourEstimator,
        _canvas_mask_to_original,
        _letterboxed_rgb,
        file_sha256,
        validate_scaleout_authorization,
    )
    from fcp_pipeline.photo_first_measurement import classify_masked_rgb

    frozen_root = args.frozen_source_root.resolve()
    module_path = Path(fcp_pipeline.__file__).resolve()
    if not module_path.is_relative_to(frozen_root):
        raise RuntimeError(f"historical biological package shadowed: {module_path}")

    import json
    roi = json.loads(
        (frozen_root / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json").read_text()
    )
    locked = json.loads(
        (frozen_root / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json").read_text()
    )
    measurement_contract = json.loads(
        (frozen_root / "docs/supporting/random_photo_first_measurement_contract_v1.json").read_text()
    )
    execution_contract = json.loads(
        (frozen_root / "docs/supporting/random_photo_first_measurement_execution_v1.json").read_text()
    )
    detector = frozen_root / "data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    detector_sha = file_sha256(detector)
    if detector_sha != execution_contract["runtime"]["detector_sha256"]:
        raise RuntimeError("historical detector hash drift")
    validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)

    threshold = measurement_contract["coarse_colour_state"]
    minimum_pixels = int(execution_contract["measurement"]["minimum_flower_mask_pixels"])
    minimum_dominant_fraction = float(threshold["minimum_dominant_fraction"])
    minimum_margin = float(threshold["minimum_margin_over_second_group"])

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    receipt = pd.read_csv(args.acquisition_receipt, dtype=str).fillna("")
    expected = pd.read_csv(args.expected_sha, dtype=str).fillna("")
    required_worker = ["measurement_id", "image_filename", "photo_license", "heavy_counterfactual"]
    if list(worker.columns) != required_worker:
        raise RuntimeError(f"biological worker schema drift: {list(worker.columns)}")
    if worker["measurement_id"].nunique() != len(worker):
        raise RuntimeError("duplicate worker measurement IDs")
    joined = worker.merge(
        receipt[["measurement_id", "acquisition_status", "source_image_sha256", "failure_reason"]],
        on="measurement_id",
        how="left",
        validate="one_to_one",
    ).merge(
        expected,
        on="measurement_id",
        how="left",
        validate="one_to_one",
        suffixes=("_pass_b", "_pass_t"),
    )
    if len(joined) != len(worker):
        raise RuntimeError("biological acquisition/SHA join lost rows")

    estimator = FrozenFlowerColourEstimator(
        detector,
        args.efficient_sam_dir,
        roi,
        torch_threads=1,
    )

    rows: list[dict[str, object]] = []
    for pos, row in enumerate(joined.itertuples(index=False), start=1):
        mid = str(row.measurement_id)
        heavy = _as_bool(row.heavy_counterfactual)
        rec: dict[str, object] = {
            "measurement_id": mid,
            "heavy_counterfactual": heavy,
            "pass_b_source_image_sha256": str(row.source_image_sha256_pass_b),
            "pass_t_source_image_sha256": str(row.source_image_sha256_pass_t),
            "source_identity_status": "",
        }

        if str(row.acquisition_status) != "acquired_and_decode_verified":
            rec["source_identity_status"] = "pass_b_acquisition_failed"
            rec.update(_empty("base", "image_acquisition_failed"))
            rows.append(rec)
            continue
        if str(row.source_image_sha256_pass_b) != str(row.source_image_sha256_pass_t):
            rec["source_identity_status"] = "source_byte_drift"
            rec.update(_empty("base", "source_byte_drift"))
            rows.append(rec)
            continue

        rec["source_identity_status"] = "exact_source_sha_match"
        image_path = args.images_dir / str(row.image_filename)
        try:
            with Image.open(image_path) as source:
                oriented = ImageOps.exif_transpose(source).convert("RGB")
            rgb = np.asarray(oriented, dtype=np.uint8)
            measured = estimator.measure(oriented)
            base_out, base_mask = _classify_pipeline(
                rgb,
                measured,
                prefix="base",
                classify_masked_rgb=classify_masked_rgb,
                minimum_pixels=minimum_pixels,
                minimum_dominant_fraction=minimum_dominant_fraction,
                minimum_margin=minimum_margin,
            )
            rec.update(base_out)

            if base_mask is not None:
                for ev in _HELPER.EV_LEVELS:
                    ev_rgb = apply_exposure_ev(rgb, ev)
                    rec.update(
                        _classify_mask(
                            ev_rgb,
                            base_mask,
                            prefix=_tag_ev("fixed_ev", ev),
                            classify_masked_rgb=classify_masked_rgb,
                            minimum_pixels=minimum_pixels,
                            minimum_dominant_fraction=minimum_dominant_fraction,
                            minimum_margin=minimum_margin,
                        )
                    )

            if heavy:
                for ev in HEAVY_EV_LEVELS:
                    ev_rgb = apply_exposure_ev(rgb, ev)
                    ev_measured = estimator.measure(Image.fromarray(ev_rgb, mode="RGB"))
                    out, _ = _classify_pipeline(
                        ev_rgb,
                        ev_measured,
                        prefix=_tag_ev("full_ev", ev),
                        classify_masked_rgb=classify_masked_rgb,
                        minimum_pixels=minimum_pixels,
                        minimum_dominant_fraction=minimum_dominant_fraction,
                        minimum_margin=minimum_margin,
                    )
                    rec.update(out)

                if base_mask is not None:
                    neutral_rgb = neutralize_background(rgb, base_mask)
                    neutral_measured = estimator.measure(Image.fromarray(neutral_rgb, mode="RGB"))
                    out, _ = _classify_pipeline(
                        neutral_rgb,
                        neutral_measured,
                        prefix="neutral_bg",
                        classify_masked_rgb=classify_masked_rgb,
                        minimum_pixels=minimum_pixels,
                        minimum_dominant_fraction=minimum_dominant_fraction,
                        minimum_margin=minimum_margin,
                    )
                    rec.update(out)
                    rec.update(
                        _prompt_variants(
                            estimator,
                            oriented,
                            rgb,
                            measured.get("boxes") or [],
                            classify_masked_rgb=classify_masked_rgb,
                            minimum_pixels=minimum_pixels,
                            minimum_dominant_fraction=minimum_dominant_fraction,
                            minimum_margin=minimum_margin,
                            box_to_canvas=box_to_canvas,
                            select_prompt_mask=select_prompt_mask,
                            letterboxed_rgb=_letterboxed_rgb,
                            canvas_mask_to_original=_canvas_mask_to_original,
                            canvas_size=CANVAS_SIZE,
                        )
                    )
        except Exception as exc:
            rec["source_identity_status"] = "exact_source_sha_match_measurement_failed"
            rec.update(_empty("base", f"measurement_runtime_failure:{type(exc).__name__}:{str(exc)[:240]}"))

        rows.append(rec)
        if pos % 5 == 0 or pos == len(joined):
            print(f"biological_measured_or_terminal={pos}/{len(joined)}", flush=True)

    out = pd.DataFrame(rows)
    if len(out) != len(worker) or out["measurement_id"].nunique() != len(worker):
        raise RuntimeError("biological terminal partition coverage drift")
    forbidden = ("species", "taxon", "photo_id", "observation_id", "latitude", "longitude", "observer")
    leaked = [c for c in out.columns if any(token in c.casefold() for token in forbidden)]
    if leaked:
        raise RuntimeError(f"biological worker output leaked metadata: {leaked}")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
