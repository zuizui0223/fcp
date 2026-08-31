#!/usr/bin/env python3
"""Prepare only the public JRC train split for the frozen ROI v4 detector."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract


SOURCE_INVENTORY = ROOT / "data/atlas/qualification/roi_v3_sources/jrc_flower_detection_source_inventory_v1.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jrc-root", type=Path, required=True)
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
    annotation_path = args.jrc_root / "annotations/train.json"
    if sha256(annotation_path) != contract["jrc_source"]["train_annotation_sha256"]:
        raise RuntimeError("JRC train annotation hash changed")
    annotation = json.loads(annotation_path.read_text(encoding="utf-8"))
    images = {int(row["id"]): row for row in annotation["images"]}
    boxes: dict[int, list[dict[str, object]]] = {image_id: [] for image_id in images}
    source_boxes = 0
    not_evaluable = 0
    for row in annotation["annotations"]:
        source_boxes += 1
        image = images[int(row["image_id"])]
        x, y, width, height = (float(value) for value in row["bbox"])
        x0 = max(0.0, min(float(image["width"]), x))
        y0 = max(0.0, min(float(image["height"]), y))
        x1 = max(0.0, min(float(image["width"]), x + width))
        y1 = max(0.0, min(float(image["height"]), y + height))
        if x1 <= x0 or y1 <= y0:
            not_evaluable += 1
            continue
        boxes[int(row["image_id"])].append(
            {"annotation_id": int(row["id"]), "xyxy": [x0, y0, x1, y1]}
        )
    if (len(images), source_boxes, not_evaluable) != (400, 6992, 1):
        raise RuntimeError("JRC training denominator changed")

    with SOURCE_INVENTORY.open(encoding="utf-8", newline="") as handle:
        inventory = {
            row["file_name"]: row
            for row in csv.DictReader(handle)
            if row["split"] == "train"
        }
    if len(inventory) != 400:
        raise RuntimeError("committed JRC train inventory changed")

    image_dir = args.output_dir / "images/train"
    label_dir = args.output_dir / "labels/train"
    image_dir.mkdir(parents=True, exist_ok=False)
    label_dir.mkdir(parents=True, exist_ok=False)
    image_rows = []
    for image_id in sorted(images):
        image = images[image_id]
        name = str(image["file_name"])
        source = args.jrc_root / "images/train" / name
        expected = inventory.get(name)
        if expected is None or sha256(source) != expected["image_sha256"]:
            raise RuntimeError(f"JRC train image identity changed: {name}")
        destination = image_dir / name
        try:
            destination.hardlink_to(source)
        except OSError:
            shutil.copy2(source, destination)
        label_lines = []
        for target in boxes[image_id]:
            x0, y0, x1, y1 = target["xyxy"]
            width = float(image["width"])
            height = float(image["height"])
            label_lines.append(
                "0 {:.12f} {:.12f} {:.12f} {:.12f}".format(
                    (x0 + x1) / (2 * width),
                    (y0 + y1) / (2 * height),
                    (x1 - x0) / width,
                    (y1 - y0) / height,
                )
            )
        (label_dir / f"{Path(name).stem}.txt").write_text(
            "\n".join(label_lines) + "\n", encoding="utf-8"
        )
        image_rows.append(
            {
                "image_id": image_id,
                "file_name": name,
                "image_sha256": expected["image_sha256"],
                "evaluable_boxes": len(boxes[image_id]),
            }
        )
    (args.output_dir / "data.yaml").write_text(
        f"path: {args.output_dir.as_posix()}\ntrain: images/train\nnames:\n  0: flower\n",
        encoding="utf-8",
    )
    rows_path = args.output_dir / "training_images.json"
    rows_path.write_text(json.dumps(image_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "protocol": contract["protocol"],
        "status": "pass_roi_v4_training_materialization",
        "images": 400,
        "source_annotation_boxes": 6992,
        "evaluable_training_boxes": 6991,
        "source_not_evaluable_boxes": 1,
        "training_images_sha256": sha256(rows_path),
        "jrc_test_directory_read": False,
        "jrc_test_images_decoded_or_scored": False,
        "scaleout_candidate_pixels_opened": False,
    }
    (args.output_dir / "training_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
