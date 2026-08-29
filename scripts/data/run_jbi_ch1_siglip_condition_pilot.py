#!/usr/bin/env python3
"""Six-image quota-independent flower-condition diagnostic with SigLIP.

The flower ROI is taken from the already-computed Florence v2 result. SigLIP receives
only that crop and a fixed set of condition descriptions. Highest similarity is recorded
without fitting a threshold. This is calibration-only diagnostic evidence, not a final
condition label for the 480 or 720 images.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from PIL import Image
from transformers import pipeline

MODEL_ID = "google/siglip-base-patch16-224"
PROTOCOL = "jbi-ch1-siglip-condition-pilot-v1"

LABELS = {
    "fresh": "a fresh healthy flower with intact non-wilted petals",
    "senescent": "a wilted dried senescent old flower with drying petals",
    "damaged": "a damaged flower with visibly injured torn or severely degraded petals",
}


def crop_box(image: Image.Image, box: list[float]) -> Image.Image:
    if len(box) != 4:
        raise ValueError("bbox must contain four coordinates")
    w, h = image.size
    x0, y0, x1, y1 = (float(v) for v in box)
    x0 = max(0, min(w - 1, int(x0)))
    y0 = max(0, min(h - 1, int(y0)))
    x1 = max(x0 + 1, min(w, int(x1 + 0.999)))
    y1 = max(y0 + 1, min(h, int(y1 + 0.999)))
    return image.crop((x0, y0, x1, y1)).convert("RGB")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-manifest", type=Path, required=True)
    parser.add_argument("--florence-v2", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.pilot_manifest)
    if len(manifest) != 6 or manifest["species"].nunique() != 6:
        raise RuntimeError("condition pilot requires the same six calibration images")
    if manifest.get("evaluation_row", pd.Series([False] * len(manifest))).astype(bool).any():
        raise RuntimeError("evaluation leakage")

    florence = json.loads(args.florence_v2.read_text(encoding="utf-8"))
    if florence.get("protocol") != "jbi-ch1-florence-colour-pilot-v2":
        raise RuntimeError("condition pilot requires Florence colour pilot v2")
    if florence.get("evaluation_rows_opened") is not False:
        raise RuntimeError("Florence pilot violated evaluation firewall")
    by_id = {r["blind_id"]: r for r in florence["results"]}
    expected = json.loads(args.expected.read_text(encoding="utf-8"))["expectations"]
    if set(by_id) != set(expected):
        raise RuntimeError("Florence and condition expectation blind IDs differ")

    classifier = pipeline(
        task="zero-shot-image-classification",
        model=MODEL_ID,
        device=-1,
    )
    label_texts = list(LABELS.values())
    reverse = {value: key for key, value in LABELS.items()}

    results = []
    for _, row in manifest.sort_values("species", kind="mergesort").iterrows():
        blind_id = str(row["blind_id"])
        image = Image.open(str(row["image_path"])).convert("RGB")
        bbox = by_id[blind_id].get("selected_bbox")
        if bbox is None or by_id[blind_id].get("localization_status") != "florence_open_vocab_box":
            observed = "unresolved"
            scores = {}
        else:
            crop = crop_box(image, bbox)
            raw = classifier(crop, candidate_labels=label_texts)
            scores = {reverse[item["label"]]: float(item["score"]) for item in raw}
            observed = max(scores, key=scores.get)
        exp = expected[blind_id]["expected_condition"]
        results.append({
            "protocol": PROTOCOL,
            "model": MODEL_ID,
            "species": str(row["species"]),
            "blind_id": blind_id,
            "evaluation_row": False,
            "pilot_only": True,
            "final_label": False,
            "expected_condition": exp,
            "observed_condition": observed,
            "match": observed == exp,
            "scores": {k: round(v, 6) for k, v in scores.items()},
            "numeric_threshold_used": False,
        })
        print(json.dumps(results[-1], ensure_ascii=False), flush=True)

    payload = {
        "protocol": PROTOCOL,
        "status": "pilot_complete_not_final",
        "model": MODEL_ID,
        "candidate_condition_prompts": LABELS,
        "n_images": 6,
        "n_predeclared_checks": 6,
        "n_matches": sum(r["match"] for r in results),
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "numeric_threshold_used": False,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
