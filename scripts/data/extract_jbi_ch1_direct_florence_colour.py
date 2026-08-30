#!/usr/bin/env python3
"""Direct, symmetric Florence-2 + continuous CIELAB measurement for Chapter 1.

The identical frozen procedure is applied independently to every calibration and
evaluation photograph:

1. Florence-2 open-vocabulary detection with the literal query ``flower``;
2. choose the largest valid returned box, with a coordinate tuple tie-break;
3. fixed inner-margin ellipse inside that box;
4. component-wise 10% trimmed mean in CIE L*a*b* (D65).

The procedure never reads latitude, longitude, graph edges, colour outcomes, species
contrasts, or Stage-A results when locating/measuring a flower.  Species is used only to
select the pre-frozen shard and to preserve the split's blocking identity.  Empty or
invalid detections fail closed; there is no whole-image fallback and no row deletion.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import torch
from PIL import Image
from transformers import AutoModelForCausalLM, AutoProcessor

from scripts.data.extract_jbi_ch1_roi_colour_features import (
    download_image,
    measure,
    normalize_box,
)

DEFAULT_MODEL_ID = "florence-community/Florence-2-base-ft"
TASK_PROMPT = "<OPEN_VOCABULARY_DETECTION>"
TEXT_QUERY = "flower"
EXPECTED_COMPONENTS = ("L_star", "a_star", "b_star")
EXPECTED_PER_SPLIT = {"calibration": 80, "evaluation": 120}
EXPECTED_SHARD_COUNTS = {"calibration": 4, "evaluation": 6}
EXPECTED_ROWS_PER_SHARD = 20


def scalar_id(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            pass
    return text


def stable_blind_id(species: str, split: str, photo_id: str) -> str:
    payload = f"jbi_ch1_direct_florence_colour_v1\t{species}\t{split}\t{photo_id}".encode()
    return hashlib.sha256(payload).hexdigest()[:24]


def read_frozen_shard(
    path: Path,
    *,
    species: str,
    split: str,
    shard_index: int,
    shard_count: int,
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("frozen split has no CSV header")
        required = {"species", "split", "photo_id", "photo_url"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"frozen split is missing {sorted(missing)}")
        rows = [
            dict(row)
            for row in reader
            if str(row.get("species", "")).strip() == species
            and str(row.get("split", "")).strip().lower() == split
        ]
    expected = EXPECTED_PER_SPLIT[split]
    if len(rows) != expected:
        raise ValueError(f"{species} {split}: expected {expected} frozen rows, found {len(rows)}")
    expected_shards = EXPECTED_SHARD_COUNTS[split]
    if shard_count != expected_shards:
        raise ValueError(
            f"{split} contract requires shard_count={expected_shards}, received {shard_count}"
        )
    rows.sort(
        key=lambda row: (
            str(row.get("split_rank_hash", "")),
            scalar_id(row.get("photo_id")),
        )
    )
    if len(rows) % shard_count:
        raise ValueError(f"{species} {split}: rows are not divisible by shard_count")
    size = len(rows) // shard_count
    if size != EXPECTED_ROWS_PER_SHARD:
        raise ValueError(f"expected {EXPECTED_ROWS_PER_SHARD} rows/shard, found {size}")
    start = shard_index * size
    return rows[start : start + size]


def _box_array(value: Any) -> np.ndarray:
    if isinstance(value, dict):
        for keys in (
            ("x1", "y1", "x2", "y2"),
            ("xmin", "ymin", "xmax", "ymax"),
            ("left", "top", "right", "bottom"),
        ):
            if all(key in value for key in keys):
                value = [value[key] for key in keys]
                break
        else:
            raise ValueError("unrecognized box object")
    array = np.asarray(value, dtype=float)
    if array.shape != (4,) or not np.isfinite(array).all():
        raise ValueError("box must be a finite four-vector")
    return array


def iter_box_candidates(value: Any, path: str = "") -> Iterable[tuple[str, np.ndarray]]:
    if isinstance(value, dict):
        # Conventional list fields get first-class treatment.
        for key in ("bboxes", "boxes"):
            if key in value and isinstance(value[key], (list, tuple)):
                for index, item in enumerate(value[key]):
                    try:
                        yield f"{path}.{key}[{index}]" if path else f"{key}[{index}]", _box_array(item)
                    except (TypeError, ValueError):
                        continue
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else str(key)
            if key in {"bboxes", "boxes"}:
                continue
            yield from iter_box_candidates(child, child_path)
    elif isinstance(value, (list, tuple)):
        # A bare four-vector can occur in processor-version-specific output.
        try:
            yield path, _box_array(value)
            return
        except (TypeError, ValueError):
            pass
        for index, child in enumerate(value):
            yield from iter_box_candidates(child, f"{path}[{index}]")


def select_largest_valid_box(
    parsed: Any,
    *,
    image_width: int,
    image_height: int,
) -> tuple[tuple[int, int, int, int], str, list[list[float]]]:
    """Select by geometry only; colour/content scores never enter the decision."""

    candidates: list[tuple[float, tuple[int, int, int, int], str, list[float]]] = []
    raw_boxes: list[list[float]] = []
    for path, raw in iter_box_candidates(parsed):
        raw_list = [float(value) for value in raw]
        raw_boxes.append(raw_list)
        try:
            normalized = normalize_box(raw, image_width, image_height, path)
        except ValueError:
            continue
        x1, y1, x2, y2 = normalized
        area = float((x2 - x1) * (y2 - y1))
        if area <= 0:
            continue
        candidates.append((area, normalized, path, raw_list))
    if not candidates:
        raise ValueError("Florence returned no valid flower box")
    # Max area; coordinate tuple and then source path make ties deterministic.
    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    _, box, path, _ = candidates[0]
    return box, path, raw_boxes


def load_model(model_id: str) -> tuple[Any, Any, dict[str, Any]]:
    processor = AutoProcessor.from_pretrained(model_id, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        model_id,
        trust_remote_code=True,
        torch_dtype=torch.float32,
    )
    model.eval()
    config = getattr(model, "config", None)
    metadata = {
        "model_id": model_id,
        "model_commit_hash": getattr(config, "_commit_hash", None),
        "model_architectures": getattr(config, "architectures", None),
        "processor_class": type(processor).__name__,
        "model_class": type(model).__name__,
        "torch_version": torch.__version__,
        "transformers_version": __import__("transformers").__version__,
        "python_version": platform.python_version(),
    }
    return processor, model, metadata


def florence_detect(processor: Any, model: Any, image: Image.Image) -> tuple[Any, str]:
    prompt = TASK_PROMPT + TEXT_QUERY
    inputs = processor(text=prompt, images=image, return_tensors="pt")
    # All tensors stay on CPU.  Generation is deterministic (no sampling).
    with torch.inference_mode():
        generated_ids = model.generate(
            **inputs,
            max_new_tokens=1024,
            num_beams=3,
            do_sample=False,
        )
    generated_text = processor.batch_decode(
        generated_ids,
        skip_special_tokens=False,
        clean_up_tokenization_spaces=False,
    )[0]
    parsed = processor.post_process_generation(
        generated_text,
        task=TASK_PROMPT,
        image_size=(image.width, image.height),
    )
    return parsed, generated_text


def process_row(
    source: dict[str, str],
    *,
    species: str,
    split: str,
    processor: Any,
    model: Any,
    model_metadata: dict[str, Any],
) -> dict[str, Any]:
    photo_id = scalar_id(source.get("photo_id"))
    output: dict[str, Any] = {
        **source,
        "photo_id": photo_id,
        "species": species,
        "split": split,
        "blind_id": stable_blind_id(species, split, photo_id),
        "protocol": "jbi_ch1_direct_florence_colour_v1",
        "feature_status": "pending",
        "colour_feature_status": "pending",
        "evaluation_row": split == "evaluation",
        "calibration_only": split == "calibration",
        "final_label": False,
        "continuous_colour_representation": "fixed_inner_ellipse_trimmed_mean_cielab_d65",
        "continuous_colour_dimension": 3,
        "continuous_colour_component_names": list(EXPECTED_COMPONENTS),
        "florence_task_prompt": TASK_PROMPT,
        "florence_text_query": TEXT_QUERY,
        "box_selection_rule": "largest_valid_area_then_coordinate_tuple_then_path",
        "measurement_uses_species_outcome": False,
        "measurement_uses_coordinates": False,
        "measurement_uses_graph": False,
        "measurement_rule_tuned_on_evaluation": False,
        "model_metadata": model_metadata,
    }
    try:
        if not photo_id:
            raise ValueError("photo_id is empty")
        image, image_sha = download_image(str(source.get("photo_url", "")))
        parsed, generated_text = florence_detect(processor, model, image)
        selected_box, selected_path, raw_boxes = select_largest_valid_box(
            parsed,
            image_width=image.width,
            image_height=image.height,
        )
        vector, diagnostics = measure(image, selected_box)
        output.update(
            {
                "feature_status": "ok",
                "feature_method": "florence_open_vocab_box_plus_fixed_roi_cielab",
                "colour_feature_status": "ok",
                "continuous_colour_vector": vector.tolist(),
                "selected_box_xyxy": list(selected_box),
                "selected_box_source_path": selected_path,
                "florence_candidate_boxes": raw_boxes,
                "florence_candidate_box_count": len(raw_boxes),
                "florence_generated_text_sha256": hashlib.sha256(
                    generated_text.encode("utf-8", "replace")
                ).hexdigest(),
                "source_image_width_px": image.width,
                "source_image_height_px": image.height,
                "source_image_sha256": image_sha,
                **diagnostics,
            }
        )
    except Exception as exc:
        output.update(
            {
                "feature_status": "failed",
                "colour_feature_status": "failed",
                "feature_error_type": type(exc).__name__,
                "feature_error": str(exc),
                "colour_feature_error_type": type(exc).__name__,
                "colour_feature_error": str(exc),
            }
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen-split", type=Path, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--split", choices=("calibration", "evaluation"), required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output-jsonl", type=Path, required=True)
    parser.add_argument("--model-id", default=os.environ.get("FLORENCE_MODEL_ID", DEFAULT_MODEL_ID))
    args = parser.parse_args()

    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index/count")
    rows = read_frozen_shard(
        args.frozen_split,
        species=args.species,
        split=args.split,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )
    processor, model, model_metadata = load_model(args.model_id)
    output: list[dict[str, Any]] = []
    failures = 0
    for index, source in enumerate(rows, start=1):
        result = process_row(
            source,
            species=args.species,
            split=args.split,
            processor=processor,
            model=model,
            model_metadata=model_metadata,
        )
        output.append(result)
        failures += result["colour_feature_status"] != "ok"
        if index % 5 == 0 or index == len(rows):
            print(
                f"{args.species} {args.split} shard {args.shard_index}/{args.shard_count}: "
                f"{index}/{len(rows)}; failures={failures}",
                flush=True,
            )

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.output_jsonl.open("w", encoding="utf-8") as handle:
        for row in output:
            handle.write(json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n")
    if len(output) != EXPECTED_ROWS_PER_SHARD:
        raise RuntimeError(f"expected {EXPECTED_ROWS_PER_SHARD} output rows, found {len(output)}")
    if failures:
        raise SystemExit(f"{failures}/{len(output)} direct Florence colour measurements failed")
    print(f"validated {len(output)} direct Florence continuous-colour rows")


if __name__ == "__main__":
    main()
