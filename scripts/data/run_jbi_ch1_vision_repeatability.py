#!/usr/bin/env python3
"""Measure independent repeatability of the six-image calibration vision pilot.

Three fresh Copilot CLI sessions are run for each of six calibration-only images.
The purpose is to measure operational repeatability, not to assign final colour states.
No evaluation image is accessed and no previous response is supplied to a later pass.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
from typing import Any

import pandas as pd

PROTOCOL = "jbi-ch1-vision-repeatability-v2"
MODEL_REQUESTED = "auto"
N_PASSES = 3

ENUMS = {
    "flower_visibility": {"evaluable", "not_evaluable"},
    "visibility_failure_code": {"", "occlusion", "distance", "blur", "overexposure", "underexposure", "non_target_organ", "insufficient_flower_area", "other"},
    "flower_condition": {"fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_applicable"},
    "flower_region": {"single_target_clear", "multiple_flowers_clear", "ambiguous", "not_applicable"},
    "within_photo_flower_consistency": {"consistent", "variable_between_flowers", "single_flower", "unresolved", "not_applicable"},
    "segmentation_feasibility": {"feasible", "uncertain", "not_feasible", "not_applicable"},
    "colour_pattern": {"approximately_uniform", "multi_colour_pattern", "unresolved", "not_applicable"},
    "diagnostic_colour_scope": {"petal_background", "pattern_markings", "both", "unresolved", "not_applicable"},
}
AGREEMENT_FIELDS = (
    "flower_visibility",
    "flower_condition",
    "flower_region",
    "within_photo_flower_consistency",
    "segmentation_feasibility",
    "colour_pattern",
    "diagnostic_colour_scope",
)


def make_prompt(species: str, blind_id: str, pass_index: int) -> str:
    return f'''You are an independent blinded reviewer in pass {pass_index} of a calibration repeatability test for flower-colour spatial ecology.
Examine ONLY the attached photograph. Species identity is provided because the later diagnostic colour rule is species-specific: {species}.
Do not infer geography, environment, adaptation, population membership, pollination, or evolutionary mechanism.
Do not assume this photograph resembles any prior photograph or prior answer.
Return ONLY one valid JSON object, no markdown and no surrounding prose.

Required schema:
{{
  "blind_id": "{blind_id}",
  "flower_visibility": "evaluable|not_evaluable",
  "visibility_failure_code": "|occlusion|distance|blur|overexposure|underexposure|non_target_organ|insufficient_flower_area|other",
  "flower_condition": "fresh|senescent|damaged|mixed_or_ambiguous|not_applicable",
  "flower_region": "single_target_clear|multiple_flowers_clear|ambiguous|not_applicable",
  "within_photo_flower_consistency": "consistent|variable_between_flowers|single_flower|unresolved|not_applicable",
  "segmentation_feasibility": "feasible|uncertain|not_feasible|not_applicable",
  "target_flower_bbox_pct": [x0,y0,x1,y1] or null,
  "primary_petals_visible": true or false,
  "apparent_petals_colour_terms": [short visual colour terms],
  "colour_pattern": "approximately_uniform|multi_colour_pattern|unresolved|not_applicable",
  "diagnostic_colour_scope": "petal_background|pattern_markings|both|unresolved|not_applicable",
  "confidence": number from 0 to 1,
  "notes": "maximum 20 words"
}}

Rules:
- flower_visibility asks only whether petal colour can visually be assessed.
- flower_condition=fresh only when the visible target petals are not clearly senescent/dried, badly damaged, or mixed in condition.
- Senescent/dried flowers can still be visible but are not suitable for defining the normal fresh-flower colour codebook.
- If multiple flowers are visible, assess whether their diagnostic colour appearance is mutually consistent. Never majority-vote discordant flowers.
- A multi-colour pattern can be normal pigmentation within one flower. Do not interpret within-flower spots, rays, throats, centers, or markings as different morphs.
- diagnostic_colour_scope asks which region would carry a species colour-state contrast if visible; choose unresolved if the image alone cannot establish that.
- For multiple consistent flowers, bbox one representative well-exposed target flower rather than the entire cluster.
- If flower_visibility=not_evaluable, set flower_condition, consistency, segmentation, colour_pattern, diagnostic_colour_scope to not_applicable and bbox=null.
- This is calibration repeatability only, not a biological inference or final label.'''


def parse_json(text: str) -> dict[str, Any]:
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


def validate(result: dict[str, Any], blind_id: str) -> None:
    required = {
        "blind_id", *ENUMS.keys(), "target_flower_bbox_pct", "primary_petals_visible",
        "apparent_petals_colour_terms", "confidence", "notes"
    }
    missing = required - set(result)
    if missing:
        raise ValueError(f"missing keys: {sorted(missing)}")
    if str(result["blind_id"]) != blind_id:
        raise ValueError("blind_id mismatch")
    for field, allowed in ENUMS.items():
        if result[field] not in allowed:
            raise ValueError(f"invalid {field}: {result[field]!r}")
    if not isinstance(result["primary_petals_visible"], bool):
        raise ValueError("primary_petals_visible must be boolean")
    if not isinstance(result["apparent_petals_colour_terms"], list):
        raise ValueError("apparent_petals_colour_terms must be a list")
    confidence = float(result["confidence"])
    if not 0 <= confidence <= 1:
        raise ValueError("confidence outside [0,1]")
    bbox = result["target_flower_bbox_pct"]
    if bbox is not None:
        if not isinstance(bbox, list) or len(bbox) != 4:
            raise ValueError("bbox must be null or four numbers")
        vals = [float(v) for v in bbox]
        if any(v < 0 or v > 100 for v in vals) or vals[2] <= vals[0] or vals[3] <= vals[1]:
            raise ValueError("invalid bbox")
    if result["flower_visibility"] == "not_evaluable":
        for field in ("flower_condition", "within_photo_flower_consistency", "segmentation_feasibility", "colour_pattern", "diagnostic_colour_scope"):
            if result[field] != "not_applicable":
                raise ValueError(f"not_evaluable requires {field}=not_applicable")
        if bbox is not None:
            raise ValueError("not_evaluable requires bbox=null")


def run_one(image: Path, species: str, blind_id: str, pass_index: int) -> dict[str, Any]:
    command = [
        "copilot", "--experimental", "--no-custom-instructions", "--no-ask-user",
        "--no-remote", "--no-remote-export", "--model", MODEL_REQUESTED,
        "--attachment", str(image.resolve()), "-s", "-p", make_prompt(species, blind_id, pass_index),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=240, env=os.environ.copy())
    if completed.returncode != 0:
        raise RuntimeError(f"Copilot CLI failed rc={completed.returncode}: {completed.stderr[-4000:]}")
    result = parse_json(completed.stdout)
    validate(result, blind_id)
    result.update({
        "species": species,
        "pass_index": pass_index,
        "model_requested": MODEL_REQUESTED,
        "protocol": PROTOCOL,
        "pilot_only": True,
        "used_for_final_calibration_rule": False,
    })
    return result


def normalize_terms(values: list[Any]) -> set[str]:
    return {str(x).strip().casefold() for x in values if str(x).strip()}


def pairwise_jaccard(sets: list[set[str]]) -> float | None:
    values = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            union = sets[i] | sets[j]
            values.append(1.0 if not union else len(sets[i] & sets[j]) / len(union))
    return sum(values) / len(values) if values else None


def summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    by_image: dict[str, list[dict[str, Any]]] = {}
    for row in records:
        by_image.setdefault(str(row["blind_id"]), []).append(row)

    image_rows = []
    for blind_id, rows in sorted(by_image.items()):
        if len(rows) != N_PASSES:
            raise ValueError(f"{blind_id}: expected {N_PASSES} passes, got {len(rows)}")
        item: dict[str, Any] = {
            "blind_id": blind_id,
            "species": rows[0]["species"],
            "n_passes": len(rows),
        }
        for field in AGREEMENT_FIELDS:
            values = [str(r[field]) for r in rows]
            counts = Counter(values)
            mode, mode_n = counts.most_common(1)[0]
            item[f"{field}_unanimous"] = len(counts) == 1
            item[f"{field}_mode"] = mode
            item[f"{field}_mode_fraction"] = mode_n / len(values)
            item[f"{field}_values"] = values
        term_sets = [normalize_terms(r["apparent_petals_colour_terms"]) for r in rows]
        item["colour_terms_pairwise_jaccard"] = pairwise_jaccard(term_sets)
        item["confidence_mean"] = sum(float(r["confidence"]) for r in rows) / len(rows)
        image_rows.append(item)

    field_summary = {}
    for field in AGREEMENT_FIELDS:
        unanimous = [bool(row[f"{field}_unanimous"]) for row in image_rows]
        field_summary[field] = {
            "n_images": len(unanimous),
            "n_unanimous": sum(unanimous),
            "unanimous_fraction": sum(unanimous) / len(unanimous),
        }
    return {
        "protocol": PROTOCOL,
        "status": "measured_not_yet_accepted_for_scaleup",
        "n_images": len(image_rows),
        "passes_per_image": N_PASSES,
        "n_valid_responses": len(records),
        "model_requested": MODEL_REQUESTED,
        "independent_sessions": True,
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "used_for_final_calibration_rule": False,
        "field_repeatability": field_summary,
        "per_image": image_rows,
        "scaleup_decision": "requires_review_of_measured_repeatability; no numeric acceptance threshold was invented post hoc"
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-manifest", type=Path, default=Path("artifacts/jbi_ch1_vision_repeatability_v2/pilot_manifest.csv"))
    parser.add_argument("--output-jsonl", type=Path, default=Path("data/calibration/jbi_ch1_vision_repeatability_v2.jsonl"))
    parser.add_argument("--summary-json", type=Path, default=Path("docs/supporting/jbi_ch1_vision_repeatability_v2.json"))
    args = parser.parse_args()

    frame = pd.read_csv(args.pilot_manifest)
    if len(frame) != 6 or frame["species"].nunique() != 6 or frame["blind_id"].nunique() != 6:
        raise RuntimeError("repeatability pilot requires exactly one calibration image per six species")
    if frame.get("evaluation_row", pd.Series([False] * len(frame))).astype(bool).any():
        raise RuntimeError("evaluation leakage detected")

    records = []
    for pass_index in range(1, N_PASSES + 1):
        for _, row in frame.sort_values("species", kind="mergesort").iterrows():
            result = run_one(Path(str(row["image_path"])), str(row["species"]), str(row["blind_id"]), pass_index)
            records.append(result)
            print(json.dumps(result, ensure_ascii=False), flush=True)

    summary = summarize(records)
    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.output_jsonl.write_text("".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records), encoding="utf-8")
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
