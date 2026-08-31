#!/usr/bin/env python3
"""Extract label-free, model-consensus flower-candidate colour states.

This module intentionally reads only the location-free development review
packet.  It never reads coordinates, dates, observers, environments, prior
C/G labels, or the locked partition.  Output from a limited run is feasibility
evidence only and must not be joined to geography.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from PIL import Image
from skimage.color import rgb2lab


PROTOCOL_VERSION = "fcp-inaturalist-automated-colour-state-v2"
DEFAULT_SEED = "fcp-inaturalist-automated-colour-state-v2"
MODEL_ID = "CIDAS/clipseg-rd64-refined"
MODEL_REVISION = "999e0328d9e10b484360c477313983f9afdd7050"
MODEL_SIZE = 352
POSITIVE_PROMPTS = ("flower", "petals", "blossom")
NEGATIVE_PROMPTS = ("leaves", "background")
PROMPTS = POSITIVE_PROMPTS + NEGATIVE_PROMPTS
MIN_EFFECTIVE_PIXELS = 100.0
MIN_VALID_POSITIVE_PROMPTS = 2
MAX_PROMPT_DELTA_E = 10.0
MAX_FLIP_DELTA_E = 5.0
MIN_FLIP_SOFT_IOU = 0.50
BACKGROUND_REQUIRED_FOR_PHOTO_ADMISSION = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_rank(seed: str, *values: object) -> str:
    token = "\x1f".join([seed, *(str(value) for value in values)])
    return hashlib.sha256(token.encode("utf-8")).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: Any) -> Any:
    """Replace non-finite diagnostics with explicit JSON null values.

    A negative-control mask may legitimately have zero effective weight when
    the flower prompts dominate.  Such a missing diagnostic must remain
    distinguishable from a measured zero, while never producing non-standard
    NaN/Infinity JSON.
    """

    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    return value


def contract_sha256(model_sha256: str) -> str:
    contract = {
        "protocol": PROTOCOL_VERSION,
        "model_id": MODEL_ID,
        "model_revision": MODEL_REVISION,
        "model_sha256": model_sha256,
        "model_size": MODEL_SIZE,
        "positive_prompts": POSITIVE_PROMPTS,
        "negative_prompts": NEGATIVE_PROMPTS,
        "min_effective_pixels": MIN_EFFECTIVE_PIXELS,
        "min_valid_positive_prompts": MIN_VALID_POSITIVE_PROMPTS,
        "max_prompt_delta_e": MAX_PROMPT_DELTA_E,
        "max_flip_delta_e": MAX_FLIP_DELTA_E,
        "min_flip_soft_iou": MIN_FLIP_SOFT_IOU,
        "background_required_for_photo_admission": BACKGROUND_REQUIRED_FOR_PHOTO_ADMISSION,
    }
    return hashlib.sha256(
        json.dumps(contract, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def select_encounters(
    rows: Iterable[dict[str, str]], limit_per_species: int, seed: str
) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["canonical_name"]].append(row)
    selected: list[dict[str, str]] = []
    for species in sorted(grouped):
        ranked = sorted(
            grouped[species],
            key=lambda row: stable_rank(seed, species, row["encounter_blind_id"]),
        )
        selected.extend(ranked if limit_per_species <= 0 else ranked[:limit_per_species])
    return selected


def sigmoid(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, -40.0, 40.0)
    return 1.0 / (1.0 + np.exp(-clipped))


def candidate_weights(logits: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-positive, positive ensemble, and negative weights."""

    if logits.shape[0] != len(PROMPTS):
        raise ValueError(f"expected {len(PROMPTS)} prompt maps, got {logits.shape}")
    positive = logits[: len(POSITIVE_PROMPTS)]
    negative = logits[len(POSITIVE_PROMPTS) :]
    max_negative = np.max(negative, axis=0)
    positive_weights = np.clip(2.0 * (sigmoid(positive - max_negative) - 0.5), 0.0, 1.0)
    ensemble = np.mean(positive_weights, axis=0)
    max_positive = np.max(positive, axis=0)
    negative_weight = np.clip(
        2.0 * (sigmoid(max_negative - max_positive) - 0.5), 0.0, 1.0
    )
    return positive_weights, ensemble, negative_weight


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(keep):
        return math.nan
    values = values[keep]
    weights = weights[keep]
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    target = probability * cumulative[-1]
    return float(values[min(np.searchsorted(cumulative, target, side="left"), len(values) - 1)])


def lab_summary(rgb: np.ndarray, weights: np.ndarray, prefix: str) -> dict[str, float]:
    flat_weights = weights.reshape(-1).astype(float)
    total = float(np.sum(flat_weights))
    result: dict[str, float] = {f"{prefix}_effective_pixels": total}
    if total <= 0 or not math.isfinite(total):
        for channel in ("L", "a", "b"):
            for suffix in ("mean", "sd", "q10", "q50", "q90"):
                result[f"{prefix}_{channel}_{suffix}"] = math.nan
        return result
    lab = rgb2lab(rgb.astype(np.float32) / 255.0).reshape(-1, 3)
    for index, channel in enumerate(("L", "a", "b")):
        values = lab[:, index].astype(float)
        mean = float(np.sum(values * flat_weights) / total)
        variance = float(np.sum(((values - mean) ** 2) * flat_weights) / total)
        result[f"{prefix}_{channel}_mean"] = mean
        result[f"{prefix}_{channel}_sd"] = math.sqrt(max(variance, 0.0))
        for suffix, probability in (("q10", 0.1), ("q50", 0.5), ("q90", 0.9)):
            result[f"{prefix}_{channel}_{suffix}"] = weighted_quantile(
                values, flat_weights, probability
            )
    return result


def mean_lab(summary: dict[str, float], prefix: str) -> np.ndarray:
    return np.array(
        [summary[f"{prefix}_L_mean"], summary[f"{prefix}_a_mean"], summary[f"{prefix}_b_mean"]],
        dtype=float,
    )


def delta_e(first: np.ndarray, second: np.ndarray) -> float:
    if not np.all(np.isfinite(first)) or not np.all(np.isfinite(second)):
        return math.inf
    return float(np.linalg.norm(first - second))


def summarize_photo(
    rgb: np.ndarray,
    original_logits: np.ndarray,
    flipped_rgb: np.ndarray,
    flipped_logits: np.ndarray,
) -> dict[str, Any]:
    prompt_weights, ensemble, negative = candidate_weights(original_logits)
    _flip_prompts, flip_ensemble, _flip_negative = candidate_weights(flipped_logits)
    ensemble_summary = lab_summary(rgb, ensemble, "flower")
    flip_summary = lab_summary(flipped_rgb, flip_ensemble, "flip")
    background_summary = lab_summary(rgb, negative, "background")
    prompt_summaries = [lab_summary(rgb, weights, f"prompt_{index}") for index, weights in enumerate(prompt_weights)]
    valid_prompt_indices = [
        index
        for index, summary in enumerate(prompt_summaries)
        if summary[f"prompt_{index}_effective_pixels"] >= MIN_EFFECTIVE_PIXELS
        and np.all(np.isfinite(mean_lab(summary, f"prompt_{index}")))
    ]
    prompt_deltas = [
        delta_e(
            mean_lab(prompt_summaries[first], f"prompt_{first}"),
            mean_lab(prompt_summaries[second], f"prompt_{second}"),
        )
        for first, second in itertools.combinations(valid_prompt_indices, 2)
    ]
    prompt_max_delta = max(prompt_deltas, default=math.inf)
    flip_delta = delta_e(mean_lab(ensemble_summary, "flower"), mean_lab(flip_summary, "flip"))
    unflipped = np.fliplr(flip_ensemble)
    denominator = float(np.sum(np.maximum(ensemble, unflipped)))
    soft_iou = (
        float(np.sum(np.minimum(ensemble, unflipped))) / denominator
        if denominator > 0
        else 0.0
    )
    flower_features_finite = all(
        math.isfinite(float(value))
        for mapping in (ensemble_summary, flip_summary)
        for value in mapping.values()
    )
    background_features_available = all(
        math.isfinite(float(value)) for value in background_summary.values()
    )
    reasons: list[str] = []
    if ensemble_summary["flower_effective_pixels"] < MIN_EFFECTIVE_PIXELS:
        reasons.append("insufficient_effective_weight")
    if len(valid_prompt_indices) < MIN_VALID_POSITIVE_PROMPTS:
        reasons.append("insufficient_valid_positive_prompts")
    if prompt_max_delta > MAX_PROMPT_DELTA_E:
        reasons.append("prompt_colour_instability")
    if flip_delta > MAX_FLIP_DELTA_E:
        reasons.append("reflection_colour_instability")
    if soft_iou < MIN_FLIP_SOFT_IOU:
        reasons.append("reflection_mask_instability")
    if not flower_features_finite:
        reasons.append("nonfinite_flower_feature")
    return {
        **ensemble_summary,
        **background_summary,
        "valid_positive_prompts": len(valid_prompt_indices),
        "prompt_max_delta_e": prompt_max_delta,
        "flip_delta_e": flip_delta,
        "flip_soft_iou": soft_iou,
        "background_features_available": background_features_available,
        "automated_colour_state_status": (
            "automated_colour_state_admitted"
            if not reasons
            else "automated_colour_state_not_evaluable"
        ),
        "failure_reasons": ";".join(reasons),
    }


def image_array(image: Image.Image) -> np.ndarray:
    resized = image.convert("RGB").resize((MODEL_SIZE, MODEL_SIZE), Image.Resampling.BILINEAR)
    return np.asarray(resized, dtype=np.uint8)


def inference_batch(model: Any, processor: Any, images: list[Image.Image]) -> np.ndarray:
    import torch

    texts: list[str] = []
    repeated_images: list[Image.Image] = []
    for image in images:
        for prompt in PROMPTS:
            texts.append(prompt)
            repeated_images.append(image)
    inputs = processor(text=texts, images=repeated_images, return_tensors="pt")
    with torch.inference_mode():
        logits = model(**inputs).logits.detach().cpu().numpy()
    return logits.reshape(len(images), len(PROMPTS), logits.shape[-2], logits.shape[-1])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review-artifact", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--public-manifest", type=Path, required=True)
    parser.add_argument("--limit-encounters-per-species", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--seed", default=DEFAULT_SEED)
    parser.add_argument("--torch-threads", type=int, default=8)
    parser.add_argument(
        "--species",
        action="append",
        default=[],
        help="Optional exact species subset for disjoint cache-only workers.",
    )
    parser.add_argument(
        "--cache-only",
        action="store_true",
        help="Populate validated per-photo cache records without aggregating outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    review_artifact = args.review_artifact.resolve()
    model_dir = args.model_dir.resolve()
    output_dir = args.output_dir.resolve()
    cache_dir = output_dir / "photo_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    model_path = model_dir / "model.safetensors"
    if sha256(model_path) != "d00ca85d6b859f9d07b7cfb8ef26fe9771cb275b34c9368f2ecf603139307f55":
        raise RuntimeError("pinned CLIPSeg safetensors hash mismatch")
    contract_hash = contract_sha256(sha256(model_path))

    reviewer_rows = read_csv(review_artifact / "reviewer_A_annotation_sheet.csv")
    selected = select_encounters(reviewer_rows, args.limit_encounters_per_species, args.seed)
    if args.species:
        unknown = sorted(set(args.species) - {row["canonical_name"] for row in reviewer_rows})
        if unknown:
            raise RuntimeError(f"unknown requested species: {unknown}")
        selected = [row for row in selected if row["canonical_name"] in set(args.species)]
    photo_inputs: list[dict[str, str]] = []
    for encounter in selected:
        for image_file in encounter["image_files"].split("|"):
            photo_inputs.append(
                {
                    "canonical_name": encounter["canonical_name"],
                    "encounter_blind_id": encounter["encounter_blind_id"],
                    "photo_blind_id": Path(image_file).stem,
                    "image_file": image_file,
                }
            )
    photo_inputs.sort(key=lambda row: (row["canonical_name"], row["encounter_blind_id"], row["photo_blind_id"]))

    pending: list[dict[str, str]] = []
    for row in photo_inputs:
        image_hash = sha256(review_artifact / row["image_file"])
        cache_path = cache_dir / f"{row['photo_blind_id']}.json"
        if cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("image_sha256") == image_hash and cached.get("contract_sha256") == contract_hash:
                continue
        pending.append(row | {"image_sha256": image_hash})

    os.environ.setdefault("USE_TF", "0")
    os.environ.setdefault("TRANSFORMERS_NO_TF", "1")
    import torch
    import transformers
    from transformers import CLIPSegForImageSegmentation, CLIPSegProcessor

    torch.set_num_threads(args.torch_threads)
    processor = CLIPSegProcessor.from_pretrained(model_dir, local_files_only=True, use_fast=False)
    model = CLIPSegForImageSegmentation.from_pretrained(model_dir, local_files_only=True).eval()
    started = time.time()
    for offset in range(0, len(pending), args.batch_size):
        batch = pending[offset : offset + args.batch_size]
        images = [Image.open(review_artifact / row["image_file"]).convert("RGB") for row in batch]
        flipped = [image.transpose(Image.Transpose.FLIP_LEFT_RIGHT) for image in images]
        batch_started = time.time()
        original_logits = inference_batch(model, processor, images)
        flipped_logits = inference_batch(model, processor, flipped)
        elapsed = time.time() - batch_started
        for index, row in enumerate(batch):
            summary = summarize_photo(
                image_array(images[index]),
                original_logits[index],
                image_array(flipped[index]),
                flipped_logits[index],
            )
            record = {
                **row,
                **summary,
                "image_sha256": row["image_sha256"],
                "contract_sha256": contract_hash,
                "model_id": MODEL_ID,
                "model_revision": MODEL_REVISION,
                "batch_elapsed_seconds": elapsed,
            }
            cache_path = cache_dir / f"{row['photo_blind_id']}.json"
            cache_path.write_text(
                json.dumps(json_safe(record), indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
        print(
            f"processed={min(offset + len(batch), len(pending))}/{len(pending)} "
            f"batch_seconds={elapsed:.3f}",
            flush=True,
        )

    if args.cache_only:
        print(
            json.dumps(
                {
                    "status": "complete_cache_only_worker",
                    "protocol": PROTOCOL_VERSION,
                    "species": sorted(args.species),
                    "selected_encounters": len(selected),
                    "selected_photos": len(photo_inputs),
                    "contract_sha256": contract_hash,
                    "cache_records_present": sum(
                        (cache_dir / f"{row['photo_blind_id']}.json").exists()
                        for row in photo_inputs
                    ),
                },
                indent=2,
                sort_keys=True,
            ),
            flush=True,
        )
        return

    photo_rows = [
        json.loads((cache_dir / f"{row['photo_blind_id']}.json").read_text(encoding="utf-8"))
        for row in photo_inputs
    ]
    encounter_rows: list[dict[str, Any]] = []
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in photo_rows:
        grouped[(row["canonical_name"], row["encounter_blind_id"])].append(row)
    feature_fields = [
        f"flower_{channel}_{metric}"
        for channel in ("L", "a", "b")
        for metric in ("mean", "sd", "q10", "q50", "q90")
    ]
    background_feature_fields = [
        f"background_{channel}_{metric}"
        for channel in ("L", "a", "b")
        for metric in ("mean", "sd", "q10", "q50", "q90")
    ]
    for (species, encounter_id), rows_for_encounter in sorted(grouped.items()):
        admitted = [
            row
            for row in rows_for_encounter
            if row["automated_colour_state_status"] == "automated_colour_state_admitted"
        ]
        aggregated: dict[str, Any] = {
            "canonical_name": species,
            "encounter_blind_id": encounter_id,
            "n_photos": len(rows_for_encounter),
            "n_admitted_photos": len(admitted),
            "encounter_status": (
                "automated_colour_state_admitted"
                if admitted
                else "automated_colour_state_not_evaluable"
            ),
        }
        for field in feature_fields:
            aggregated[field] = statistics.median(float(row[field]) for row in admitted) if admitted else ""
        background_available = [
            row for row in admitted if row["background_features_available"] is True
        ]
        aggregated["n_background_control_photos"] = len(background_available)
        aggregated["background_control_status"] = (
            "background_control_available"
            if background_available
            else "background_control_not_evaluable"
        )
        for field in background_feature_fields:
            aggregated[field] = (
                statistics.median(float(row[field]) for row in background_available)
                if background_available
                else ""
            )
        pairwise = [
            delta_e(mean_lab(first, "flower"), mean_lab(second, "flower"))
            for first, second in itertools.combinations(admitted, 2)
        ]
        aggregated["within_encounter_photo_delta_e_median"] = (
            statistics.median(pairwise) if pairwise else ""
        )
        encounter_rows.append(aggregated)

    output_dir.mkdir(parents=True, exist_ok=True)
    photo_path = output_dir / "photo_features.csv"
    encounter_path = output_dir / "encounter_features.csv"
    write_csv(photo_path, photo_rows)
    write_csv(encounter_path, encounter_rows)
    by_species: list[dict[str, Any]] = []
    for species in sorted({row["canonical_name"] for row in encounter_rows}):
        species_encounters = [row for row in encounter_rows if row["canonical_name"] == species]
        species_photos = [row for row in photo_rows if row["canonical_name"] == species]
        admitted_encounters = sum(
            row["encounter_status"] == "automated_colour_state_admitted"
            for row in species_encounters
        )
        background_encounters = sum(
            row["background_control_status"] == "background_control_available"
            for row in species_encounters
        )
        by_species.append(
            {
                "canonical_name": species,
                "selected_encounters": len(species_encounters),
                "selected_photos": len(species_photos),
                "admitted_photos": sum(
                    row["automated_colour_state_status"] == "automated_colour_state_admitted"
                    for row in species_photos
                ),
                "admitted_encounters": admitted_encounters,
                "admitted_encounter_share": admitted_encounters / len(species_encounters),
                "background_control_encounters": background_encounters,
                "background_control_encounter_share": (
                    background_encounters / admitted_encounters if admitted_encounters else 0.0
                ),
            }
        )

    private_manifest = {
        "status": "complete_automated_colour_state_development_feasibility_not_spatial",
        "protocol": PROTOCOL_VERSION,
        "source_artifact": str(review_artifact),
        "selection": (
            "all development encounters"
            if args.limit_encounters_per_species <= 0
            else f"stable outcome-blind {args.limit_encounters_per_species} encounters per species"
        ),
        "limit_encounters_per_species": args.limit_encounters_per_species,
        "seed": args.seed,
        "model_id": MODEL_ID,
        "model_revision": MODEL_REVISION,
        "model_file_sha256": {path.name: sha256(path) for path in sorted(model_dir.glob("*")) if path.is_file()},
        "contract_sha256": contract_hash,
        "extractor_source_sha256": sha256(Path(__file__).resolve()),
        "prompts": {"positive": POSITIVE_PROMPTS, "negative": NEGATIVE_PROMPTS},
        "gates": {
            "min_effective_pixels": MIN_EFFECTIVE_PIXELS,
            "min_valid_positive_prompts": MIN_VALID_POSITIVE_PROMPTS,
            "max_prompt_delta_e": MAX_PROMPT_DELTA_E,
            "max_flip_delta_e": MAX_FLIP_DELTA_E,
            "min_flip_soft_iou": MIN_FLIP_SOFT_IOU,
            "background_required_for_photo_admission": BACKGROUND_REQUIRED_FOR_PHOTO_ADMISSION,
        },
        "selected_encounters": len(encounter_rows),
        "selected_photos": len(photo_rows),
        "species_summary": by_species,
        "elapsed_seconds_new_inference": time.time() - started,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "torch": torch.__version__,
            "transformers": transformers.__version__,
            "torch_threads": args.torch_threads,
            "batch_size": args.batch_size,
            "device": "cpu",
        },
        "private_output_sha256": {
            photo_path.name: sha256(photo_path),
            encounter_path.name: sha256(encounter_path),
        },
        "claim_ceiling": (
            "Image-only automated measurement feasibility. No coordinates were read; this run cannot "
            "estimate spatial randomness, morph frequencies, botanical morph identity, boundaries or universality."
        ),
    }
    private_manifest_path = output_dir / "run_manifest.json"
    private_manifest_path.write_text(
        json.dumps(private_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    public = {
        key: value
        for key, value in private_manifest.items()
        if key not in {"source_artifact", "private_output_sha256"}
    }
    public["source_artifact_reference"] = "private location-free development review packet"
    public["run_manifest_sha256"] = sha256(private_manifest_path)
    args.public_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.public_manifest.write_text(
        json.dumps(public, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(public, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
