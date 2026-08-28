#!/usr/bin/env python3
"""Run a six-image Copilot CLI vision pilot and validate structured outputs.

Diagnostic only: these six calibration images do not become final labels and no
held-out evaluation image is opened.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess

import pandas as pd

PROTOCOL = "jbi-ch1-vision-pilot-v1"
# GitHub Copilot CLI currently documents this exact model string and GitHub's image
# input documentation uses gpt-5.4 as a vision-session example.
MODEL = "gpt-5.4"

VISIBILITY = {"evaluable", "not_evaluable"}
FAILURE_CODES = {"", "occlusion", "distance", "blur", "overexposure", "underexposure", "non_target_organ", "insufficient_flower_area", "other"}
REGION = {"single_target_clear", "multiple_flowers_clear", "ambiguous", "not_applicable"}
SEGMENTATION = {"feasible", "uncertain", "not_feasible", "not_applicable"}
PATTERN = {"approximately_uniform", "multi_colour_pattern", "unresolved", "not_applicable"}


def make_prompt(species: str, blind_id: str) -> str:
    return f'''You are performing a blinded calibration pilot for a flower-colour spatial ecology dataset.
Examine ONLY the attached photograph. Species identity is provided only because later colour coding is species-specific: {species}.
Do not infer geography, environment, population membership, adaptation, pollination, or evolutionary mechanism.
Return ONLY one valid JSON object, with no markdown fences or surrounding prose.

Required schema:
{{
  "blind_id": "{blind_id}",
  "flower_visibility": "evaluable|not_evaluable",
  "visibility_failure_code": "|occlusion|distance|blur|overexposure|underexposure|non_target_organ|insufficient_flower_area|other",
  "flower_region": "single_target_clear|multiple_flowers_clear|ambiguous|not_applicable",
  "segmentation_feasibility": "feasible|uncertain|not_feasible|not_applicable",
  "target_flower_bbox_pct": [x0,y0,x1,y1] or null,
  "primary_petals_visible": true or false,
  "apparent_petals_colour_terms": [short visual colour terms],
  "colour_pattern": "approximately_uniform|multi_colour_pattern|unresolved|not_applicable",
  "confidence": number from 0 to 1,
  "notes": "maximum 20 words"
}}
Rules:
- evaluable means petal/flower colour can be judged with enough visible target-flower area.
- Ignore leaves, stems, soil, sky, labels, and background objects when describing petal colour.
- Do not force a colour or segmentation decision; use unresolved/uncertain when needed.
- If not_evaluable, set colour_pattern=not_applicable and target_flower_bbox_pct=null.
- Bounding-box values are rough percentages from 0 to 100.
- This is calibration, not biological inference.'''


def parse_json(text: str) -> dict:
    text = text.strip()
    if text.startswith("```"):
        lines = text.splitlines()[1:]
        if lines and lines[-1].strip() == "```":
            lines.pop()
        text = "\n".join(lines).strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        start, end = text.find("{"), text.rfind("}")
        if start >= 0 and end > start:
            return json.loads(text[start:end + 1])
        raise


def validate(result: dict, blind_id: str) -> None:
    required = {
        "blind_id", "flower_visibility", "visibility_failure_code", "flower_region",
        "segmentation_feasibility", "target_flower_bbox_pct", "primary_petals_visible",
        "apparent_petals_colour_terms", "colour_pattern", "confidence", "notes"
    }
    missing = required - set(result)
    if missing:
        raise ValueError(f"missing keys: {sorted(missing)}")
    if str(result["blind_id"]) != blind_id:
        raise ValueError("blind_id mismatch")
    if result["flower_visibility"] not in VISIBILITY:
        raise ValueError("invalid flower_visibility")
    if result["visibility_failure_code"] not in FAILURE_CODES:
        raise ValueError("invalid visibility_failure_code")
    if result["flower_region"] not in REGION:
        raise ValueError("invalid flower_region")
    if result["segmentation_feasibility"] not in SEGMENTATION:
        raise ValueError("invalid segmentation_feasibility")
    if result["colour_pattern"] not in PATTERN:
        raise ValueError("invalid colour_pattern")
    if not isinstance(result["primary_petals_visible"], bool):
        raise ValueError("primary_petals_visible must be boolean")
    if not isinstance(result["apparent_petals_colour_terms"], list):
        raise ValueError("colour terms must be a list")
    confidence = float(result["confidence"])
    if not 0 <= confidence <= 1:
        raise ValueError("confidence outside [0,1]")
    bbox = result["target_flower_bbox_pct"]
    if bbox is not None:
        if not isinstance(bbox, list) or len(bbox) != 4:
            raise ValueError("bbox must be null or four values")
        values = [float(v) for v in bbox]
        if any(v < 0 or v > 100 for v in values) or values[2] <= values[0] or values[3] <= values[1]:
            raise ValueError("invalid bbox")
    if result["flower_visibility"] == "not_evaluable":
        if bbox is not None or result["colour_pattern"] != "not_applicable":
            raise ValueError("not_evaluable output is internally inconsistent")


def run_one(image: Path, species: str, blind_id: str) -> dict:
    command = [
        "copilot", "--experimental", "--no-custom-instructions", "--no-ask-user",
        "--no-remote", "--no-remote-export", "--model", MODEL,
        "--attachment", str(image.resolve()), "-s", "-p", make_prompt(species, blind_id),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=240, env=os.environ.copy())
    if completed.returncode != 0:
        raise RuntimeError(f"Copilot CLI failed rc={completed.returncode}: {completed.stderr[-4000:]}")
    result = parse_json(completed.stdout)
    validate(result, blind_id)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-manifest", type=Path, default=Path("artifacts/jbi_ch1_vision_pilot_v1/pilot_manifest.csv"))
    parser.add_argument("--output-jsonl", type=Path, default=Path("data/calibration/jbi_ch1_vision_pilot_v1.jsonl"))
    parser.add_argument("--manifest-json", type=Path, default=Path("docs/supporting/jbi_ch1_vision_pilot_manifest_v1.json"))
    args = parser.parse_args()

    frame = pd.read_csv(args.pilot_manifest)
    if len(frame) != 6 or frame["species"].nunique() != 6:
        raise RuntimeError("pilot must contain exactly six species/images")
    if frame.get("evaluation_row", pd.Series([False] * len(frame))).astype(bool).any():
        raise RuntimeError("evaluation leakage detected")

    records = []
    for _, row in frame.sort_values("species", kind="mergesort").iterrows():
        species, bid = str(row["species"]), str(row["blind_id"])
        result = run_one(Path(str(row["image_path"])), species, bid)
        result.update({
            "species": species,
            "model": MODEL,
            "protocol": PROTOCOL,
            "pilot_only": True,
            "used_for_final_calibration_rule": False,
        })
        records.append(result)
        print(json.dumps(result, ensure_ascii=False), flush=True)

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.output_jsonl.write_text("".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records), encoding="utf-8")
    summary = {
        "protocol": PROTOCOL,
        "status": "pass",
        "model": MODEL,
        "n_images": 6,
        "species_count": 6,
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "pilot_only": True,
        "used_for_final_calibration_rule": False,
        "output_jsonl": str(args.output_jsonl),
        "results": records,
    }
    args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
