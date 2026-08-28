#!/usr/bin/env python3
"""Run a six-image Copilot CLI vision pilot and validate structured outputs.

The pilot is diagnostic only. It does not write final calibration labels and must not
be scaled to all 480 images until the output schema and agreement are reviewed.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess

import pandas as pd


PROTOCOL = "jbi-ch1-vision-pilot-v1"
MODEL = "gpt-5-mini"

VISIBILITY = {"evaluable", "not_evaluable"}
FAILURE_CODES = {
    "",
    "occlusion",
    "distance",
    "blur",
    "overexposure",
    "underexposure",
    "non_target_organ",
    "insufficient_flower_area",
    "other",
}
REGION = {"single_target_clear", "multiple_flowers_clear", "ambiguous", "not_applicable"}
SEGMENTATION = {"feasible", "uncertain", "not_feasible", "not_applicable"}
PATTERN = {"approximately_uniform", "multi_colour_pattern", "unresolved", "not_applicable"}


def prompt(species: str, blind_id: str) -> str:
    return f"""You are performing a blinded calibration pilot for a flower-colour spatial ecology dataset.
Examine ONLY the attached photograph. Species identity is provided only because later colour coding is species-specific: {species}.
Do not infer geography, environment, population membership, adaptation, pollination, or evolutionary mechanism.
Return ONLY one valid JSON object, with no markdown fences and no prose before or after it.

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
- evaluable means petal/flower colour can be judged from the photograph with enough visible target-flower area.
- Ignore leaves, stems, soil, sky, labels, and other background objects when describing petal colour.
- Do not force a colour or segmentation decision if ambiguous; use unresolved/uncertain.
- If flower_visibility is not_evaluable, set colour_pattern to not_applicable and target_flower_bbox_pct to null.
- Bounding box values, when used, are rough percentages of image width/height from 0 to 100.
- This is a calibration pilot, not a biological inference.
"""


def parse_json_response(text: str) -> dict:
    stripped = text.strip()
    if stripped.startswith("```"):
        lines = stripped.splitlines()
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip() == "```":
            lines = lines[:-1]
        stripped = "\n".join(lines).strip()
    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        start = stripped.find("{")
        end = stripped.rfind("}")
        if start >= 0 and end > start:
            return json.loads(stripped[start : end + 1])
        raise


def validate_result(result: dict, expected_blind_id: str) -> None:
    required = {
        "blind_id",
        "flower_visibility",
        "visibility_failure_code",
        "flower_region",
        "segmentation_feasibility",
        "target_flower_bbox_pct",
        "primary_petals_visible",
        "apparent_petals_colour_terms",
        "colour_pattern",
        "confidence",
        "notes",
    }
    missing = required - set(result)
    if missing:
        raise ValueError(f"vision result missing keys: {sorted(missing)}")
    if str(result["blind_id"]) != expected_blind_id:
        raise ValueError("vision result blind_id mismatch")
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
        raise ValueError("apparent_petals_colour_terms must be a list")
    confidence = float(result["confidence"])
    if not 0 <= confidence <= 1:
        raise ValueError("confidence must be in [0,1]")
    bbox = result["target_flower_bbox_pct"]
    if bbox is not None:
        if not isinstance(bbox, list) or len(bbox) != 4:
            raise ValueError("target_flower_bbox_pct must be null or four numbers")
        values = [float(v) for v in bbox]
        if any(v < 0 or v > 100 for v in values):
            raise ValueError("bbox values must be in [0,100]")
        if values[2] <= values[0] or values[3] <= values[1]:
            raise ValueError("bbox must have positive width and height")
    if result["flower_visibility"] == "not_evaluable":
        if result["target_flower_bbox_pct"] is not None:
            raise ValueError("not_evaluable result must have null bbox")
        if result["colour_pattern"] != "not_applicable":
            raise ValueError("not_evaluable result must have colour_pattern=not_applicable")


def run_one(image_path: Path, species: str, blind_id: str, *, timeout: int = 240) -> tuple[dict, str]:
    command = [
        "copilot",
        "--experimental",
        "--no-custom-instructions",
        "--no-ask-user",
        "--no-remote",
        "--no-remote-export",
        "--model",
        MODEL,
        "--attachment",
        str(image_path.resolve()),
        "-s",
        "-p",
        prompt(species, blind_id),
    ]
    env = os.environ.copy()
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Copilot CLI failed rc={completed.returncode}: stderr={completed.stderr[-4000:]}"
        )
    result = parse_json_response(completed.stdout)
    validate_result(result, blind_id)
    return result, completed.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pilot-manifest",
        type=Path,
        default=Path("artifacts/jbi_ch1_vision_pilot_v1/pilot_manifest.csv"),
    )
    parser.add_argument(
        "--output-jsonl",
        type=Path,
        default=Path("data/calibration/jbi_ch1_vision_pilot_v1.jsonl"),
    )
    parser.add_argument(
        "--manifest-json",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_vision_pilot_manifest_v1.json"),
    )
    args = parser.parse_args()

    frame = pd.read_csv(args.pilot_manifest)
    if len(frame) != 6 or frame["species"].nunique() != 6:
        raise RuntimeError("vision pilot must contain exactly six species/images")
    if frame.get("evaluation_row", pd.Series([False] * len(frame))).astype(bool).any():
        raise RuntimeError("evaluation row leaked into vision pilot")

    records = []
    for _, row in frame.sort_values("species", kind="mergesort").iterrows():
        species = str(row["species"])
        blind_id = str(row["blind_id"])
        image_path = Path(str(row["image_path"]))
        result, raw = run_one(image_path, species, blind_id)
        result["species"] = species
        result["model"] = MODEL
        result["protocol"] = PROTOCOL
        result["pilot_only"] = True
        result["used_for_final_calibration_rule"] = False
        records.append(result)
        print(json.dumps(result, ensure_ascii=False), flush=True)

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.output_jsonl.write_text(
        "".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records),
        encoding="utf-8",
    )
    summary = {
        "protocol": PROTOCOL,
        "status": "pass",
        "model": MODEL,
        "n_images": len(records),
        "species_count": len({r["species"] for r in records}),
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "pilot_only": True,
        "used_for_final_calibration_rule": False,
        "output_jsonl": str(args.output_jsonl),
        "results": records,
    }
    args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_json.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
