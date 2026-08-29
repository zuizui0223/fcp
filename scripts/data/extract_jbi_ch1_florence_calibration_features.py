#!/usr/bin/env python3
"""Extract quota-independent flower ROI colour features for one calibration species.

This worker processes only frozen calibration rows. It reuses the exact Florence model,
fixed palette, species score mapping, and geometry-only ROI selector from the validated
six-image v2 pilot. It never opens evaluation rows and never emits a final label.
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import time
from urllib.request import Request, urlopen

import pandas as pd
from PIL import Image, ImageOps
import torch
from transformers import AutoModelForMultimodalLM, AutoProcessor

from fcp_pipeline.florence_roi import box_area_fraction, choose_flower_box

PROTOCOL = "jbi-ch1-florence-calibration-features-v1"
USER_AGENT = "zuizui0223-fcp-jbi-ch1-florence-features/1.0 (research reproducibility)"


def load_pilot_module():
    path = Path(__file__).with_name("run_jbi_ch1_florence_colour_pilot.py")
    spec = importlib.util.spec_from_file_location("jbi_florence_pilot_v2", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load frozen Florence pilot implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.PROTOCOL != "jbi-ch1-florence-colour-pilot-v2":
        raise RuntimeError(f"expected validated pilot v2 implementation, found {module.PROTOCOL}")
    return module


def candidate_urls(row: pd.Series) -> list[str]:
    urls: list[str] = []
    for col in ("photo_url", "photo_url_api"):
        if col in row.index and not pd.isna(row[col]):
            value = str(row[col]).strip()
            if value and value not in urls:
                urls.append(value)
    return urls


def download_image(row: pd.Series, output: Path) -> str:
    errors: list[str] = []
    for url in candidate_urls(row):
        for attempt in range(1, 4):
            try:
                req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
                with urlopen(req, timeout=60) as response:
                    payload = response.read()
                if len(payload) < 1024:
                    raise RuntimeError(f"response too small: {len(payload)} bytes")
                with Image.open(io.BytesIO(payload)) as image:
                    image = ImageOps.exif_transpose(image)
                    image.load()
                    image.convert("RGB").save(output, format="JPEG", quality=95)
                return url
            except Exception as exc:
                errors.append(f"{url} attempt={attempt}: {type(exc).__name__}: {exc}")
                time.sleep(0.5 * attempt)
    raise RuntimeError("; ".join(errors[-8:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument(
        "--split",
        type=Path,
        default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    pilot = load_pilot_module()
    split = pd.read_csv(args.split)
    required = {"species", "photo_id", "split", "photo_url"}
    missing = required - set(split.columns)
    if missing:
        raise RuntimeError(f"split missing required columns: {sorted(missing)}")

    evaluation_ids = set(
        split.loc[split["split"].astype(str).eq("evaluation"), "photo_id"].astype(str)
    )
    rows = split.loc[
        split["split"].astype(str).eq("calibration")
        & split["species"].astype(str).eq(args.species)
    ].copy()
    if len(rows) != 80:
        raise RuntimeError(f"{args.species}: expected 80 calibration rows, found {len(rows)}")
    if set(rows["photo_id"].astype(str)) & evaluation_ids:
        raise RuntimeError("evaluation photo leaked into calibration worker")

    processor = AutoProcessor.from_pretrained(pilot.MODEL_ID)
    model = AutoModelForMultimodalLM.from_pretrained(pilot.MODEL_ID)
    model.eval()

    records = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        rows = rows.sort_values("split_rank_hash", kind="mergesort")
        for index, (_, row) in enumerate(rows.iterrows(), start=1):
            photo_id = str(row["photo_id"])
            blind = __import__("hashlib").sha256(
                f"jbi-ch1-calibration-v1\x1f{args.species}\x1f{photo_id}".encode()
            ).hexdigest()[:16]
            image_path = tmpdir / f"{blind}.jpg"
            used_url = download_image(row, image_path)
            image = Image.open(image_path).convert("RGB")

            detection = pilot.run_task(
                model,
                processor,
                image,
                "<OPEN_VOCABULARY_DETECTION>",
                "flower",
            )
            boxes = pilot.extract_boxes(detection, "<OPEN_VOCABULARY_DETECTION>")
            prompt_used = "flower"
            if not boxes:
                detection = pilot.run_task(
                    model,
                    processor,
                    image,
                    "<OPEN_VOCABULARY_DETECTION>",
                    "flower petals",
                )
                boxes = pilot.extract_boxes(detection, "<OPEN_VOCABULARY_DETECTION>")
                prompt_used = "flower petals"

            box = choose_flower_box(boxes, image.size)
            if box is None:
                record = {
                    "protocol": PROTOCOL,
                    "model": pilot.MODEL_ID,
                    "species": args.species,
                    "photo_id": photo_id,
                    "blind_id": blind,
                    "evaluation_row": False,
                    "calibration_only": True,
                    "final_label": False,
                    "localization_status": "no_valid_flower_box",
                    "n_detected_boxes": len(boxes),
                    "detection_prompt": prompt_used,
                    "candidate_state_numeric": "unresolved",
                    "feature_status": "localization_failed",
                    "downloaded_from": used_url,
                }
            else:
                counts = pilot.nearest_palette_counts(image, box)
                fractions = pilot.flower_only_fractions(counts)
                candidate, scores = pilot.deterministic_candidate(args.species, fractions)
                area_fractions = [box_area_fraction(b, image.size) for b in boxes]
                record = {
                    "protocol": PROTOCOL,
                    "model": pilot.MODEL_ID,
                    "species": args.species,
                    "photo_id": photo_id,
                    "blind_id": blind,
                    "evaluation_row": False,
                    "calibration_only": True,
                    "final_label": False,
                    "localization_status": "florence_open_vocab_box",
                    "n_detected_boxes": len(boxes),
                    "detection_prompt": prompt_used,
                    "selected_bbox": [round(float(x), 2) for x in box],
                    "selected_bbox_area_fraction": round(
                        float(box_area_fraction(box, image.size)), 6
                    ),
                    "detected_box_area_fraction_min": round(min(area_fractions), 6),
                    "detected_box_area_fraction_max": round(max(area_fractions), 6),
                    "palette_counts": counts,
                    "flower_only_fractions": {
                        k: round(float(v), 6) for k, v in fractions.items()
                    },
                    "candidate_scores": {
                        k: round(float(v), 6) for k, v in scores.items()
                    },
                    "candidate_state_numeric": candidate,
                    "feature_status": "ok",
                    "downloaded_from": used_url,
                }
            records.append(record)
            if index % 10 == 0:
                print(
                    f"{args.species}: {index}/80 "
                    f"{record['localization_status']} {record['candidate_state_numeric']}",
                    flush=True,
                )

    if len(records) != 80 or any(r["evaluation_row"] for r in records):
        raise RuntimeError("calibration feature output violated frozen row contract")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
