#!/usr/bin/env python3
"""Independently re-score state-bearing Chapter 1 calibration records.

The second reviewer never receives the first-pass candidate state or confidence. It sees
only the frozen photograph, species identity, and the same predeclared literature-backed
candidate codebook. Final calibration consensus is created only after both passes are
compared. Evaluation rows are forbidden.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
from typing import Any
from urllib.request import Request, urlopen

import pandas as pd
from PIL import Image, ImageOps

PROTOCOL = "jbi-ch1-calibration-independent-rescore-v1"
MODEL_REQUESTED = "auto"
USER_AGENT = "zuizui0223-fcp-jbi-ch1-rescore/1.0 (research reproducibility)"

VISIBILITY = {"evaluable", "not_evaluable"}
FAILURE = {"", "occlusion", "distance", "blur", "overexposure", "underexposure", "non_target_organ", "insufficient_flower_area", "other"}
CONDITION = {"fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_applicable"}
REGION = {"single_target_clear", "multiple_flowers_clear", "ambiguous", "not_applicable"}
CONSISTENCY = {"consistent", "variable_between_flowers", "single_flower", "unresolved", "not_applicable"}
SEGMENTATION = {"feasible", "uncertain", "not_feasible", "not_applicable"}


def blind_id(species: str, photo_id: str) -> str:
    return hashlib.sha256(f"jbi-ch1-calibration-v1\x1f{species}\x1f{photo_id}".encode()).hexdigest()[:16]


def candidate_urls(row: pd.Series) -> list[str]:
    urls: list[str] = []
    for col in ("photo_url", "photo_url_api"):
        if col in row.index and not pd.isna(row[col]):
            value = str(row[col]).strip()
            if value and value not in urls:
                urls.append(value)
    return urls


def download_image(row: pd.Series, output: Path) -> None:
    errors = []
    for url in candidate_urls(row):
        for attempt in range(1, 4):
            try:
                req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
                with urlopen(req, timeout=60) as response:
                    payload = response.read()
                if len(payload) < 1024:
                    raise RuntimeError(f"response too small: {len(payload)}")
                with Image.open(io.BytesIO(payload)) as image:
                    image = ImageOps.exif_transpose(image)
                    image.load()
                    image.convert("RGB").save(output, format="JPEG", quality=95)
                return
            except Exception as exc:
                errors.append(f"{url} attempt={attempt}: {type(exc).__name__}: {exc}")
                time.sleep(attempt * 0.5)
    raise RuntimeError("; ".join(errors[-8:]))


def prompt(species: str, bid: str, cfg: dict[str, Any]) -> str:
    definitions = cfg.get("state_definitions", {})
    lines = "\n".join(f"- {state}: {definitions.get(state, '')}" for state in cfg["candidate_states"])
    return f'''You are an INDEPENDENT second blinded reviewer of one calibration photograph for a flower-colour spatial ecology study.
You have NOT been given the first reviewer's label, confidence, or notes. Judge only this image using the predeclared species codebook.
Species: {species}.
Do not infer geography, environment, population, adaptation, pollination, or evolution.
Do not invent a state. Return ONLY one valid JSON object.

Candidate states:
{lines}
Diagnostic scope: {cfg.get('diagnostic_scope', '')}
Species note: {cfg.get('structure_note', '')}

Required schema:
{{
  "blind_id": "{bid}",
  "flower_visibility": "evaluable|not_evaluable",
  "visibility_failure_code": "|occlusion|distance|blur|overexposure|underexposure|non_target_organ|insufficient_flower_area|other",
  "flower_condition": "fresh|senescent|damaged|mixed_or_ambiguous|not_applicable",
  "flower_region": "single_target_clear|multiple_flowers_clear|ambiguous|not_applicable",
  "within_photo_flower_consistency": "consistent|variable_between_flowers|single_flower|unresolved|not_applicable",
  "segmentation_feasibility": "feasible|uncertain|not_feasible|not_applicable",
  "candidate_state": "one exact candidate state",
  "candidate_state_confidence": 0.0,
  "notes": "maximum 20 words"
}}

Rules:
- Only fresh, visually evaluable, feasibly segmentable flowers can receive a non-unresolved state.
- Senescent/damaged/mixed-condition flowers => unresolved.
- Multiple flowers with discordant diagnostic colour => unresolved; never majority-vote.
- Within-flower markings are not separate morphs unless the species codebook defines whole-pattern distribution as diagnostic.
- For Antirrhinum, classify whole-corolla pigment distribution, not hue alone.
- If any doubt prevents a clean codebook match, choose unresolved.
- This is calibration review, not a biological inference.'''


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


def validate(result: dict[str, Any], bid: str, allowed_states: set[str]) -> None:
    required = {
        "blind_id", "flower_visibility", "visibility_failure_code", "flower_condition",
        "flower_region", "within_photo_flower_consistency", "segmentation_feasibility",
        "candidate_state", "candidate_state_confidence", "notes"
    }
    missing = required - set(result)
    if missing:
        raise ValueError(f"missing keys: {sorted(missing)}")
    if str(result["blind_id"]) != bid:
        raise ValueError("blind_id mismatch")
    enums = {
        "flower_visibility": VISIBILITY,
        "visibility_failure_code": FAILURE,
        "flower_condition": CONDITION,
        "flower_region": REGION,
        "within_photo_flower_consistency": CONSISTENCY,
        "segmentation_feasibility": SEGMENTATION,
    }
    for field, allowed in enums.items():
        if result[field] not in allowed:
            raise ValueError(f"invalid {field}: {result[field]!r}")
    if result["candidate_state"] not in allowed_states:
        raise ValueError("candidate_state outside frozen codebook")
    confidence = float(result["candidate_state_confidence"])
    if not 0 <= confidence <= 1:
        raise ValueError("confidence outside [0,1]")
    gate_fail = (
        result["flower_visibility"] != "evaluable"
        or result["flower_condition"] != "fresh"
        or result["segmentation_feasibility"] != "feasible"
        or result["within_photo_flower_consistency"] in {"variable_between_flowers", "unresolved", "not_applicable"}
    )
    if gate_fail and result["candidate_state"] != "unresolved":
        raise ValueError("measurement gate requires unresolved")


def run_copilot(image: Path, species: str, bid: str, cfg: dict[str, Any]) -> dict[str, Any]:
    command = [
        "copilot", "--experimental", "--no-custom-instructions", "--no-ask-user",
        "--no-remote", "--no-remote-export", "--model", MODEL_REQUESTED,
        "--attachment", str(image.resolve()), "-s", "-p", prompt(species, bid, cfg),
    ]
    errors = []
    for attempt in range(1, 4):
        completed = subprocess.run(command, capture_output=True, text=True, timeout=240, env=os.environ.copy())
        if completed.returncode == 0:
            try:
                result = parse_json(completed.stdout)
                validate(result, bid, set(cfg["candidate_states"]))
                return result
            except Exception as exc:
                errors.append(f"attempt={attempt} invalid: {type(exc).__name__}: {exc}")
        else:
            errors.append(f"attempt={attempt} rc={completed.returncode}: {completed.stderr[-800:]}")
        time.sleep(attempt * 2)
    raise RuntimeError(" | ".join(errors[-3:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=2)
    parser.add_argument("--screening", type=Path, default=Path("data/calibration/jbi_ch1_semantic_screening_v1.jsonl"))
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument("--codebook", type=Path, default=Path("docs/supporting/jbi_ch1_species_colour_candidate_codebook_v1.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    first = [json.loads(x) for x in args.screening.read_text(encoding="utf-8").splitlines() if x.strip()]
    candidates = [
        r for r in first
        if r["species"] == args.species
        and bool(r.get("state_usable_for_calibration"))
        and r.get("candidate_state") != "unresolved"
        and r.get("evaluation_row") is False
    ]
    candidates = sorted(candidates, key=lambda r: r["blind_id"])
    candidates = [r for i, r in enumerate(candidates) if i % args.shard_count == args.shard_index]

    split = pd.read_csv(args.split)
    eval_ids = set(split.loc[split["split"].astype(str).eq("evaluation"), "photo_id"].astype(str))
    calibration = split.loc[split["split"].astype(str).eq("calibration")].copy()
    by_photo = {str(row["photo_id"]): row for _, row in calibration.iterrows()}
    codebook = json.loads(args.codebook.read_text(encoding="utf-8"))
    cfg = codebook["species"][args.species]

    records = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for n, first_row in enumerate(candidates, start=1):
            photo_id = str(first_row["photo_id"])
            if photo_id in eval_ids or photo_id not in by_photo:
                raise RuntimeError("evaluation leakage or missing calibration photo")
            bid = blind_id(args.species, photo_id)
            if bid != str(first_row["blind_id"]):
                raise RuntimeError("blind ID mismatch with frozen split")
            image = tmp / f"{bid}.jpg"
            download_image(by_photo[photo_id], image)
            second = run_copilot(image, args.species, bid, cfg)
            second.update({
                "species": args.species,
                "photo_id": photo_id,
                "protocol": PROTOCOL,
                "model_requested": MODEL_REQUESTED,
                "independent_of_first_candidate_state": True,
                "first_candidate_state_not_supplied_to_reviewer": True,
                "evaluation_row": False,
                "final_label": False,
                "shard_index": args.shard_index,
                "shard_count": args.shard_count,
            })
            records.append(second)
            print(f"{args.species} independent rescore {n}/{len(candidates)} {bid} {second['candidate_state']}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
