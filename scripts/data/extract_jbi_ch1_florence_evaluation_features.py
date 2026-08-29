#!/usr/bin/env python3
"""Extract frozen Florence colour features from the 720 evaluation photographs.

The worker is evaluation-only. It imports and reuses the validated calibration worker's
image download, sharding, biological-axis helper, Florence pilot implementation and
geometry-only ROI selector. It emits numeric features but no final biological labels.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile

import pandas as pd
from PIL import Image
from transformers import AutoModelForMultimodalLM, AutoProcessor

PROTOCOL = "jbi-ch1-florence-evaluation-features-v1"


def load_calibration_worker():
    path = Path(__file__).with_name("extract_jbi_ch1_florence_calibration_features.py")
    spec = importlib.util.spec_from_file_location("jbi_calibration_worker", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load frozen calibration worker: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.PROTOCOL != "jbi-ch1-florence-calibration-features-v1":
        raise RuntimeError(f"unexpected calibration worker protocol: {module.PROTOCOL}")
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument("--representation", type=Path, default=Path("docs/supporting/jbi_ch1_continuous_colour_representation_v1.json"))
    parser.add_argument("--opening-contract", type=Path, default=Path("docs/supporting/jbi_ch1_evaluation_opening_contract_v1.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    args = parser.parse_args()

    calibration = load_calibration_worker()
    pilot = calibration.load_pilot_module()
    representation = json.loads(args.representation.read_text(encoding="utf-8"))
    opening = json.loads(args.opening_contract.read_text(encoding="utf-8"))
    if representation.get("status") != "frozen_before_evaluation_values_inspected":
        raise RuntimeError("continuous representation was not frozen before evaluation")
    if opening.get("user_authorized_evaluation_opening") is not True:
        raise RuntimeError("evaluation opening lacks explicit authorization")
    if args.species not in representation.get("per_species", {}):
        raise RuntimeError(f"species absent from frozen representation: {args.species}")

    split = pd.read_csv(args.split)
    required = {"species", "photo_id", "split", "photo_url", "split_rank_hash"}
    missing = required - set(split.columns)
    if missing:
        raise RuntimeError(f"split missing required columns: {sorted(missing)}")
    calibration_ids = set(split.loc[split["split"].astype(str).eq("calibration"), "photo_id"].astype(str))
    all_rows = split.loc[
        split["split"].astype(str).eq("evaluation")
        & split["species"].astype(str).eq(args.species)
    ].copy()
    if len(all_rows) != 120:
        raise RuntimeError(f"{args.species}: expected 120 evaluation rows, found {len(all_rows)}")
    if set(all_rows["photo_id"].astype(str)) & calibration_ids:
        raise RuntimeError("calibration photo leaked into evaluation worker")

    rows = calibration.select_shard(all_rows, args.shard_index, args.shard_count)
    expected = len(range(args.shard_index, 120, args.shard_count))
    if len(rows) != expected:
        raise RuntimeError(f"shard row count mismatch: expected {expected}, found {len(rows)}")

    processor = AutoProcessor.from_pretrained(pilot.MODEL_ID)
    model = AutoModelForMultimodalLM.from_pretrained(pilot.MODEL_ID)
    model.eval()

    records = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        for index, (_, row) in enumerate(rows.iterrows(), start=1):
            photo_id = str(row["photo_id"])
            blind = hashlib.sha256(
                f"jbi-ch1-evaluation-v1\x1f{args.species}\x1f{photo_id}".encode("utf-8")
            ).hexdigest()[:16]
            image_path = tmpdir / f"{blind}.jpg"
            used_url = calibration.download_image(row, image_path)
            image = Image.open(image_path).convert("RGB")

            detection = pilot.run_task(model, processor, image, "<OPEN_VOCABULARY_DETECTION>", "flower")
            boxes = pilot.extract_boxes(detection, "<OPEN_VOCABULARY_DETECTION>")
            prompt_used = "flower"
            if not boxes:
                detection = pilot.run_task(model, processor, image, "<OPEN_VOCABULARY_DETECTION>", "flower petals")
                boxes = pilot.extract_boxes(detection, "<OPEN_VOCABULARY_DETECTION>")
                prompt_used = "flower petals"

            box = calibration.choose_flower_box(boxes, image.size)
            base = {
                "protocol": PROTOCOL,
                "model": pilot.MODEL_ID,
                "species": args.species,
                "photo_id": photo_id,
                "blind_id": blind,
                "evaluation_row": True,
                "calibration_only": False,
                "final_label": False,
                "evaluation_feature_measurement": True,
                "frozen_representation_protocol": representation["protocol"],
                "detection_prompt": prompt_used,
                "n_detected_boxes": len(boxes),
                "downloaded_from": used_url,
                "compute_shard_index": args.shard_index,
                "compute_shard_count": args.shard_count,
                "direct_palette_candidate_is_final_state": False,
            }
            if box is None:
                record = {
                    **base,
                    "localization_status": "no_valid_flower_box",
                    "candidate_state_numeric": "unresolved",
                    "feature_status": "localization_failed",
                }
            else:
                counts = pilot.nearest_palette_counts(image, box)
                fractions = pilot.flower_only_fractions(counts)
                candidate, scores = pilot.deterministic_candidate(args.species, fractions)
                area_fractions = [calibration.box_area_fraction(b, image.size) for b in boxes]
                record = {
                    **base,
                    "localization_status": "florence_open_vocab_box",
                    "selected_bbox": [round(float(x), 2) for x in box],
                    "selected_bbox_area_fraction": round(float(calibration.box_area_fraction(box, image.size)), 6),
                    "detected_box_area_fraction_min": round(min(area_fractions), 6),
                    "detected_box_area_fraction_max": round(max(area_fractions), 6),
                    "palette_counts": counts,
                    "flower_only_fractions": {k: round(float(v), 6) for k, v in fractions.items()},
                    "candidate_scores": {k: round(float(v), 6) for k, v in scores.items()},
                    "candidate_state_numeric": candidate,
                    "feature_status": "ok",
                    **calibration.biological_axis_features(args.species, fractions),
                }
            records.append(record)
            if index % 5 == 0 or index == len(rows):
                print(
                    f"{args.species} evaluation shard {args.shard_index}/{args.shard_count}: "
                    f"{index}/{len(rows)} {record['localization_status']}",
                    flush=True,
                )

    if len(records) != expected:
        raise RuntimeError("evaluation shard emitted wrong number of records")
    if any(r["evaluation_row"] is not True or r["calibration_only"] is not False or r["final_label"] is not False for r in records):
        raise RuntimeError("evaluation feature output violated frozen opening contract")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
