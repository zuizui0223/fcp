#!/usr/bin/env python3
"""Run a blinded SAM2 rescue pilot on reviewer-1 calibration ROIs.

The pilot uses exactly one rescue_segment record and one usable-fresh control for each
species that has rescue records. Selection is based only on the already-frozen blinded
review order. Colour scores, geography, observer metadata and dates are never used.

SAM2 receives only the previously selected Florence bounding box. Its output remains a
technical ROI-recovery diagnostic; this script never emits a biological colour state or
opens any evaluation image.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageOps

PROTOCOL = "jbi-ch1-sam2-roi-rescue-pilot-v1"
MODEL_ID = "facebook/sam2.1-hiera-tiny"
TRANSFORMERS_VERSION = "5.16.1"
USER_AGENT = "zuizui0223-fcp-jbi-ch1-sam2-rescue/1.0 (research reproducibility)"
EXPECTED_SPECIES = [
    "Antirrhinum majus",
    "Dactylorhiza sambucina",
    "Gentiana lutea",
    "Ipomoea purpurea",
    "Raphanus sativus",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_features(path: Path) -> list[dict]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 480:
        raise RuntimeError(f"expected 480 calibration features, found {len(rows)}")
    if len({str(row["blind_id"]) for row in rows}) != 480:
        raise RuntimeError("duplicate blind IDs")
    if len({str(row["photo_id"]) for row in rows}) != 480:
        raise RuntimeError("duplicate photo IDs")
    if any(row.get("evaluation_row") is not False for row in rows):
        raise RuntimeError("evaluation row present in feature input")
    if any(row.get("calibration_only") is not True or row.get("final_label") is not False for row in rows):
        raise RuntimeError("feature calibration/final-label firewall violation")
    return rows


def decision_at(spec: dict, ordinal: int) -> tuple[str, str]:
    roi = "usable"
    if ordinal in set(spec.get("rescue_segment", [])):
        roi = "rescue_segment"
    if ordinal in set(spec.get("invalid", [])):
        roi = "invalid"
    if ordinal in set(spec.get("ambiguous", [])):
        roi = "ambiguous"

    condition = "fresh"
    if ordinal in set(spec.get("senescent", [])):
        condition = "senescent"
    if ordinal in set(spec.get("damaged", [])):
        condition = "damaged"
    if ordinal in set(spec.get("mixed_or_ambiguous", [])):
        condition = "mixed_or_ambiguous"
    if ordinal in set(spec.get("not_evaluable", [])):
        condition = "not_evaluable"
    return roi, condition


def select_pilot_rows(features: list[dict], review: dict) -> list[dict]:
    if review.get("calibration_only") is not True or review.get("evaluation_rows_opened") is not False:
        raise RuntimeError("review firewall violation")
    if review.get("final_label") is not False:
        raise RuntimeError("review unexpectedly contains final labels")

    selected: list[dict] = []
    by_species: dict[str, list[dict]] = {}
    for species in sorted({str(row["species"]) for row in features}):
        group = sorted([row for row in features if str(row["species"]) == species], key=lambda row: str(row["blind_id"]))
        if len(group) != 80:
            raise RuntimeError(f"{species}: expected 80 feature rows, found {len(group)}")
        by_species[species] = group

    review_species = review.get("species", {})
    species_with_rescue = sorted(
        species for species, spec in review_species.items() if spec.get("rescue_segment")
    )
    if species_with_rescue != sorted(EXPECTED_SPECIES):
        raise RuntimeError(f"unexpected rescue species: {species_with_rescue}")

    for species in EXPECTED_SPECIES:
        spec = review_species[species]
        group = by_species[species]
        rescue_ordinal = min(int(x) for x in spec["rescue_segment"])
        rescue_roi, rescue_condition = decision_at(spec, rescue_ordinal)
        if rescue_roi != "rescue_segment":
            raise RuntimeError(f"{species}: selected rescue ordinal is not rescue_segment")
        rescue = dict(group[rescue_ordinal - 1])
        rescue.update({
            "pilot_arm": "rescue",
            "review_order_within_species": rescue_ordinal,
            "reviewer1_roi_validity": rescue_roi,
            "reviewer1_condition": rescue_condition,
        })
        selected.append(rescue)

        control = None
        for ordinal in range(1, 81):
            roi, condition = decision_at(spec, ordinal)
            if roi == "usable" and condition == "fresh":
                control = dict(group[ordinal - 1])
                control.update({
                    "pilot_arm": "control",
                    "review_order_within_species": ordinal,
                    "reviewer1_roi_validity": roi,
                    "reviewer1_condition": condition,
                })
                break
        if control is None:
            raise RuntimeError(f"{species}: no usable-fresh control")
        selected.append(control)

    if len(selected) != 10:
        raise RuntimeError(f"expected 10 pilot rows, found {len(selected)}")
    if len({str(row["blind_id"]) for row in selected}) != 10:
        raise RuntimeError("pilot selection contains duplicate blind IDs")

    # Sheet order is deterministic but hides species/arm. The hash is based only on blind ID.
    selected.sort(key=lambda row: hashlib.sha256(f"sam2-pilot-v1|{row['blind_id']}".encode()).hexdigest())
    for slot, row in enumerate(selected, start=1):
        row["pilot_slot"] = slot
    return selected


def download_image(url: str) -> Image.Image:
    errors: list[str] = []
    for attempt in range(1, 4):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
            with urlopen(req, timeout=60) as response:
                payload = response.read()
            if len(payload) < 1024:
                raise RuntimeError(f"response too small: {len(payload)} bytes")
            with Image.open(io.BytesIO(payload)) as image:
                return ImageOps.exif_transpose(image).convert("RGB").copy()
        except Exception as exc:
            errors.append(f"attempt={attempt}: {type(exc).__name__}: {exc}")
            time.sleep(0.5 * attempt)
    raise RuntimeError("; ".join(errors[-3:]))


def clip_box(box: list[float], size: tuple[int, int]) -> list[float]:
    if not isinstance(box, list) or len(box) != 4:
        raise ValueError("selected_bbox must contain four coordinates")
    width, height = size
    x0, y0, x1, y1 = map(float, box)
    x0 = max(0.0, min(float(width), x0))
    x1 = max(0.0, min(float(width), x1))
    y0 = max(0.0, min(float(height), y0))
    y1 = max(0.0, min(float(height), y1))
    if x1 <= x0 or y1 <= y0:
        raise ValueError("invalid clipped bounding box")
    return [x0, y0, x1, y1]


def mask_metrics(mask: np.ndarray, box: list[float]) -> dict[str, float | list[int]]:
    if mask.ndim != 2:
        raise ValueError("mask must be 2D")
    binary = mask.astype(bool)
    height, width = binary.shape
    pixels = int(binary.sum())
    if pixels == 0:
        raise ValueError("empty SAM2 mask")
    x0, y0, x1, y1 = box
    ix0 = max(0, min(width, int(np.floor(x0))))
    iy0 = max(0, min(height, int(np.floor(y0))))
    ix1 = max(ix0 + 1, min(width, int(np.ceil(x1))))
    iy1 = max(iy0 + 1, min(height, int(np.ceil(y1))))
    prompt_area = max(1, (ix1 - ix0) * (iy1 - iy0))
    inside = int(binary[iy0:iy1, ix0:ix1].sum())
    ys, xs = np.where(binary)
    mask_box = [int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1]
    return {
        "mask_area_fraction_image": round(pixels / float(width * height), 6),
        "mask_to_prompt_area_ratio": round(pixels / float(prompt_area), 6),
        "mask_pixels_inside_prompt_fraction": round(inside / float(pixels), 6),
        "mask_bbox": mask_box,
    }


def crop_with_context(image: Image.Image, box: list[float], pad_fraction: float = 0.12) -> Image.Image:
    width, height = image.size
    x0, y0, x1, y1 = map(float, box)
    bw, bh = x1 - x0, y1 - y0
    pad = max(bw, bh) * pad_fraction
    left = max(0, int(np.floor(x0 - pad)))
    top = max(0, int(np.floor(y0 - pad)))
    right = min(width, int(np.ceil(x1 + pad)))
    bottom = min(height, int(np.ceil(y1 + pad)))
    return image.crop((left, top, right, bottom)).convert("RGB")


def masked_image(image: Image.Image, mask: np.ndarray) -> Image.Image:
    arr = np.asarray(image, dtype=np.uint8)
    binary = mask.astype(bool)
    out = np.full_like(arr, 255)
    out[binary] = arr[binary]
    return Image.fromarray(out, mode="RGB")


def sam2_segment(model, processor, image: Image.Image, box: list[float], torch_module) -> tuple[np.ndarray, float]:
    inputs = processor(images=image, input_boxes=[[box]], return_tensors="pt")
    with torch_module.no_grad():
        outputs = model(**inputs, multimask_output=True)
    scores = outputs.iou_scores.detach().cpu().numpy()[0, 0]
    best = int(np.argmax(scores))
    masks = processor.post_process_masks(outputs.pred_masks.detach().cpu(), inputs["original_sizes"])[0]
    mask_tensor = masks
    while getattr(mask_tensor, "ndim", 0) > 3 and mask_tensor.shape[0] == 1:
        mask_tensor = mask_tensor[0]
    if mask_tensor.ndim == 3:
        mask_tensor = mask_tensor[best]
    elif mask_tensor.ndim != 2:
        raise RuntimeError(f"unexpected SAM2 mask shape: {tuple(mask_tensor.shape)}")
    mask = mask_tensor.detach().cpu().numpy().astype(bool)
    return mask, float(scores[best])


def pair_tile(before: Image.Image, masked: Image.Image, slot: int, blind_id: str) -> Image.Image:
    panel_w, panel_h, header_h = 360, 300, 50
    tile = Image.new("RGB", (panel_w * 2, panel_h + header_h), "white")
    font = ImageFont.load_default()
    draw = ImageDraw.Draw(tile)
    draw.text((8, 7), f"slot {slot:02d} | {blind_id}", fill="black", font=font)
    draw.text((8, 27), "before", fill="black", font=font)
    draw.text((panel_w + 8, 27), "SAM2 masked", fill="black", font=font)
    for col, image in enumerate((before, masked)):
        contained = ImageOps.contain(image, (panel_w, panel_h), Image.Resampling.LANCZOS)
        x = col * panel_w + (panel_w - contained.width) // 2
        y = header_h + (panel_h - contained.height) // 2
        tile.paste(contained, (x, y))
    return tile


def save_sheets(tiles: list[Image.Image], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    per_page = 5
    for start in range(0, len(tiles), per_page):
        subset = tiles[start:start + per_page]
        page = Image.new("RGB", (720, len(subset) * 350), "white")
        for index, tile in enumerate(subset):
            page.paste(tile, (0, index * 350))
        path = output_dir / f"sam2_rescue_pilot_page_{start // per_page + 1:02d}.jpg"
        page.save(path, quality=92)
        paths.append(str(path))
    return paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"))
    parser.add_argument("--review", type=Path, default=Path("docs/supporting/jbi_ch1_blind_roi_condition_review_r1_v1.json"))
    parser.add_argument("--protocol", type=Path, default=Path("docs/supporting/jbi_ch1_sam2_rescue_pilot_protocol_v1.json"))
    parser.add_argument("--output-jsonl", type=Path, default=Path("data/calibration/jbi_ch1_sam2_rescue_pilot_v1.jsonl"))
    parser.add_argument("--manifest", type=Path, default=Path("docs/supporting/jbi_ch1_sam2_rescue_pilot_manifest_v1.json"))
    parser.add_argument("--sheets-dir", type=Path, default=Path("artifacts/jbi_ch1_sam2_rescue_pilot_v1/contact_sheets"))
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if protocol.get("protocol") != PROTOCOL or protocol.get("status") != "frozen_before_pilot_execution":
        raise RuntimeError("SAM2 pilot protocol is not frozen")
    if protocol.get("evaluation_rows_opened") is not False or protocol.get("final_label") is not False:
        raise RuntimeError("SAM2 pilot protocol firewall violation")
    if protocol.get("model") != MODEL_ID or protocol.get("transformers_version") != TRANSFORMERS_VERSION:
        raise RuntimeError("SAM2 model/version differs from frozen protocol")

    features = load_features(args.features)
    review = json.loads(args.review.read_text(encoding="utf-8"))
    rows = select_pilot_rows(features, review)

    import torch
    import transformers
    from transformers import Sam2Model, Sam2Processor

    if transformers.__version__ != TRANSFORMERS_VERSION:
        raise RuntimeError(f"expected transformers {TRANSFORMERS_VERSION}, found {transformers.__version__}")
    processor = Sam2Processor.from_pretrained(MODEL_ID)
    model = Sam2Model.from_pretrained(MODEL_ID)
    model.eval()

    output_rows: list[dict] = []
    tiles: list[Image.Image] = []
    failures: list[dict] = []
    for row in rows:
        blind_id = str(row["blind_id"])
        slot = int(row["pilot_slot"])
        try:
            image_url = str(row.get("downloaded_from", "")).strip()
            if not image_url:
                raise RuntimeError("feature row has no frozen downloaded_from URL")
            image = download_image(image_url)
            box = clip_box(list(row.get("selected_bbox", [])), image.size)
            mask, score = sam2_segment(model, processor, image, box, torch)
            metrics = mask_metrics(mask, box)
            before = crop_with_context(image, box)
            after = crop_with_context(masked_image(image, mask), box)
            tiles.append(pair_tile(before, after, slot, blind_id))
            record = {
                "protocol": PROTOCOL,
                "pilot_slot": slot,
                "blind_id": blind_id,
                "photo_id": str(row["photo_id"]),
                "species": str(row["species"]),
                "pilot_arm": str(row["pilot_arm"]),
                "review_order_within_species": int(row["review_order_within_species"]),
                "reviewer1_roi_validity": str(row["reviewer1_roi_validity"]),
                "reviewer1_condition": str(row["reviewer1_condition"]),
                "source_florence_bbox": [round(float(x), 2) for x in box],
                "sam2_predicted_iou": round(score, 6),
                **metrics,
                "sam2_visual_decision": "",
                "calibration_only": True,
                "evaluation_row": False,
                "final_label": False,
                "colour_state_label": None,
                "pilot_only": True,
            }
        except Exception as exc:
            failures.append({"pilot_slot": slot, "blind_id": blind_id, "error": f"{type(exc).__name__}: {exc}"})
            blank = Image.new("RGB", (360, 260), "white")
            ImageDraw.Draw(blank).text((10, 10), "SAM2 pilot failure", fill="black")
            tiles.append(pair_tile(blank, blank, slot, blind_id))
            record = {
                "protocol": PROTOCOL,
                "pilot_slot": slot,
                "blind_id": blind_id,
                "photo_id": str(row["photo_id"]),
                "species": str(row["species"]),
                "pilot_arm": str(row["pilot_arm"]),
                "review_order_within_species": int(row["review_order_within_species"]),
                "reviewer1_roi_validity": str(row["reviewer1_roi_validity"]),
                "reviewer1_condition": str(row["reviewer1_condition"]),
                "sam2_visual_decision": "",
                "calibration_only": True,
                "evaluation_row": False,
                "final_label": False,
                "colour_state_label": None,
                "pilot_only": True,
                "error": failures[-1]["error"],
            }
        output_rows.append(record)
        print(f"slot {slot:02d} {blind_id}: {'ok' if 'error' not in record else 'failed'}", flush=True)

    output_rows.sort(key=lambda row: int(row["pilot_slot"]))
    if len(output_rows) != 10 or len({row["blind_id"] for row in output_rows}) != 10:
        raise RuntimeError("SAM2 pilot output row contract violated")
    if any(row["evaluation_row"] or row["final_label"] for row in output_rows):
        raise RuntimeError("SAM2 pilot output firewall violated")
    if any(row.get("colour_state_label") is not None for row in output_rows):
        raise RuntimeError("SAM2 pilot unexpectedly created a colour state")

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.output_jsonl.write_text(
        "".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in output_rows),
        encoding="utf-8",
    )
    sheets = save_sheets(sorted(tiles, key=lambda tile: 0), args.sheets_dir)
    # tiles were appended in deterministic pilot-slot order already.

    manifest = {
        "protocol": PROTOCOL,
        "status": "pilot_images_generated_pending_blind_visual_decision",
        "n_rows": 10,
        "n_rescue": sum(row["pilot_arm"] == "rescue" for row in output_rows),
        "n_control": sum(row["pilot_arm"] == "control" for row in output_rows),
        "n_failures": len(failures),
        "failures": failures,
        "model": MODEL_ID,
        "transformers_version": TRANSFORMERS_VERSION,
        "source_feature_sha256": sha256(args.features),
        "source_review_sha256": sha256(args.review),
        "frozen_protocol_sha256": sha256(args.protocol),
        "output_jsonl": str(args.output_jsonl),
        "contact_sheets": sheets,
        "sheet_exposes_species": False,
        "sheet_exposes_pilot_arm": False,
        "sheet_exposes_colour_scores": False,
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "visual_decisions_completed": False,
        "scaleup_rule": protocol["scaleup_rule"],
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
