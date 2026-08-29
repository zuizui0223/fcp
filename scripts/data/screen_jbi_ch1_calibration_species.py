#!/usr/bin/env python3
"""One-pass semantic screening of a frozen Chapter 1 calibration shard.

This is calibration-only and screening-only. It never opens evaluation rows and never
creates new biological colour states. Candidate states are constrained by the frozen
literature-backed species codebook.
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

PROTOCOL = "jbi-ch1-semantic-screening-v1"
MODEL_REQUESTED = "auto"
USER_AGENT = "zuizui0223-fcp-jbi-ch1-screen/1.0 (research reproducibility)"

VISIBILITY = {"evaluable", "not_evaluable"}
FAILURE = {"", "occlusion", "distance", "blur", "overexposure", "underexposure", "non_target_organ", "insufficient_flower_area", "other"}
CONDITION = {"fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_applicable"}
REGION = {"single_target_clear", "multiple_flowers_clear", "ambiguous", "not_applicable"}
CONSISTENCY = {"consistent", "variable_between_flowers", "single_flower", "unresolved", "not_applicable"}
SEGMENTATION = {"feasible", "uncertain", "not_feasible", "not_applicable"}


def blind_id(species: str, photo_id: str) -> str:
    return hashlib.sha256(f"jbi-ch1-calibration-v1\x1f{species}\x1f{photo_id}".encode()).hexdigest()[:16]


def candidate_urls(row: pd.Series) -> list[str]:
    out: list[str] = []
    for col in ("photo_url", "photo_url_api"):
        if col in row.index and not pd.isna(row[col]):
            value = str(row[col]).strip()
            if value and value not in out:
                out.append(value)
    return out


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
                    rgb = image.convert("RGB")
                    # Normalize container encoding only; do not change visible colour.
                    rgb.save(output, format="JPEG", quality=95)
                return
            except Exception as exc:
                errors.append(f"{url} attempt={attempt}: {type(exc).__name__}: {exc}")
                time.sleep(attempt * 0.5)
    raise RuntimeError("; ".join(errors[-8:]))


def prompt(species: str, bid: str, cfg: dict[str, Any]) -> str:
    states = cfg["candidate_states"]
    definitions = cfg.get("state_definitions", {})
    state_lines = "\n".join(
        f"- {state}: {definitions.get(state, '')}" for state in states
    )
    return f'''You are screening one BLINDED CALIBRATION photograph for a flower-colour spatial ecology study.
Use ONLY the attached image plus the predeclared species codebook below. Species: {species}.
Do not infer geography, environment, population, adaptation, pollination, or evolution.
Do not invent a new colour state. Return ONLY one valid JSON object and no markdown.

Predeclared candidate states:
{state_lines}
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
  "candidate_state": "one exact predeclared candidate state",
  "candidate_state_confidence": 0.0,
  "apparent_petals_colour_terms": ["short visual terms"],
  "notes": "maximum 20 words"
}}

Rules:
- flower_visibility concerns whether flower/petal colour can be judged visually.
- Only fresh flowers may receive a non-unresolved candidate colour state.
- Senescent, damaged, or mixed-condition flowers must have candidate_state=unresolved.
- If multiple visible flowers disagree in diagnostic colour appearance, set within_photo_flower_consistency=variable_between_flowers and candidate_state=unresolved. Never majority-vote.
- Within-flower spots, rays, throats, veins, centers, or markings are not separate morphs unless the species codebook explicitly defines the whole pattern as diagnostic.
- If visibility is not_evaluable or segmentation is not feasible/uncertain, candidate_state=unresolved.
- For Antirrhinum, match the whole-corolla pigment distribution pattern; do not classify background hue alone.
- `unresolved` is a valid and preferred state whenever the photo does not support a clean codebook assignment.
- This output is screening-only and is not a final biological label.'''


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
        "candidate_state", "candidate_state_confidence", "apparent_petals_colour_terms", "notes"
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
        raise ValueError(f"candidate_state not in frozen codebook: {result['candidate_state']!r}")
    confidence = float(result["candidate_state_confidence"])
    if not 0 <= confidence <= 1:
        raise ValueError("candidate_state_confidence outside [0,1]")
    if not isinstance(result["apparent_petals_colour_terms"], list):
        raise ValueError("apparent_petals_colour_terms must be list")

    must_unresolve = (
        result["flower_visibility"] != "evaluable"
        or result["flower_condition"] != "fresh"
        or result["segmentation_feasibility"] != "feasible"
        or result["within_photo_flower_consistency"] in {"variable_between_flowers", "unresolved", "not_applicable"}
    )
    if must_unresolve and result["candidate_state"] != "unresolved":
        raise ValueError("measurement gate requires candidate_state=unresolved")


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
                errors.append(f"attempt={attempt} invalid_response: {type(exc).__name__}: {exc}; stdout={completed.stdout[-1200:]}")
        else:
            errors.append(f"attempt={attempt} copilot_rc={completed.returncode}: {completed.stderr[-1200:]}")
        time.sleep(attempt * 2)
    raise RuntimeError(" | ".join(errors[-3:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=2)
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument("--codebook", type=Path, default=Path("docs/supporting/jbi_ch1_species_colour_candidate_codebook_v1.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.shard_index < 0 or args.shard_index >= args.shard_count:
        raise ValueError("invalid shard index")
    split = pd.read_csv(args.split)
    codebook = json.loads(args.codebook.read_text(encoding="utf-8"))
    if args.species not in codebook["species"]:
        raise ValueError(f"species absent from codebook: {args.species}")
    cfg = codebook["species"][args.species]

    rows = split.loc[(split["split"].astype(str) == "calibration") & (split["species"].astype(str) == args.species)].copy()
    if len(rows) != 80:
        raise RuntimeError(f"expected 80 calibration rows for {args.species}, found {len(rows)}")
    rows["_blind"] = [blind_id(args.species, str(p)) for p in rows["photo_id"].astype(str)]
    rows = rows.sort_values("_blind", kind="mergesort").reset_index(drop=True)
    rows = rows.loc[rows.index % args.shard_count == args.shard_index].copy()

    records = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        for n, (_, row) in enumerate(rows.iterrows(), start=1):
            photo_id = str(row["photo_id"])
            bid = blind_id(args.species, photo_id)
            image_path = tmpdir / f"{bid}.jpg"
            download_image(row, image_path)
            result = run_copilot(image_path, args.species, bid, cfg)
            usable = (
                result["flower_visibility"] == "evaluable"
                and result["flower_condition"] == "fresh"
                and result["segmentation_feasibility"] == "feasible"
                and result["within_photo_flower_consistency"] in {"consistent", "single_flower"}
                and result["candidate_state"] != "unresolved"
            )
            result.update({
                "species": args.species,
                "photo_id": photo_id,
                "protocol": PROTOCOL,
                "model_requested": MODEL_REQUESTED,
                "shard_index": args.shard_index,
                "shard_count": args.shard_count,
                "screening_only": True,
                "final_label": False,
                "state_usable_for_calibration": usable,
                "evaluation_row": False,
            })
            records.append(result)
            print(f"{args.species} shard {args.shard_index}: {n}/{len(rows)} {bid} {result['candidate_state']}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in records), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
