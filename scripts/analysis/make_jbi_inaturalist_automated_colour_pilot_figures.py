#!/usr/bin/env python3
"""Build the species-free atlas display and locked pilot result figure."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from PIL import Image, ImageOps
from skimage.color import lab2rgb

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.data.extract_inaturalist_automated_colour_states import (
    candidate_weights,
    inference_batch,
)


ADMITTED = "automated_colour_state_admitted"
PHOTO_BAR_COUNT = 24


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def select_photo_bar_encounters(
    rows: list[dict[str, str]], count: int = PHOTO_BAR_COUNT
) -> list[dict[str, str]]:
    admitted = sorted(
        (row for row in rows if row["encounter_status"] == ADMITTED),
        key=lambda row: (float(row["longitude"]), row["encounter_blind_id"]),
    )
    if len(admitted) < count:
        raise ValueError("not enough admitted encounters for photo bar")
    positions = np.rint(np.linspace(0, len(admitted) - 1, count)).astype(int)
    if len(set(positions.tolist())) != count:
        raise RuntimeError("photo-bar quantile selection produced duplicate positions")
    return [admitted[index] for index in positions]


def mass_box(weights: np.ndarray, quantile: float = 0.05, pad: float = 0.12) -> tuple[int, int, int, int]:
    """Return a padded box containing the central soft-mask mass."""

    values = np.asarray(weights, dtype=float)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or values.sum() <= 0:
        raise ValueError("invalid soft mask")
    y_mass = values.sum(axis=1)
    x_mass = values.sum(axis=0)
    y_cdf = np.cumsum(y_mass) / y_mass.sum()
    x_cdf = np.cumsum(x_mass) / x_mass.sum()
    y0 = int(np.searchsorted(y_cdf, quantile, side="left"))
    y1 = int(np.searchsorted(y_cdf, 1.0 - quantile, side="left")) + 1
    x0 = int(np.searchsorted(x_cdf, quantile, side="left"))
    x1 = int(np.searchsorted(x_cdf, 1.0 - quantile, side="left")) + 1
    height, width = values.shape
    extra_y = int(round((y1 - y0) * pad))
    extra_x = int(round((x1 - x0) * pad))
    return (
        max(0, x0 - extra_x),
        max(0, y0 - extra_y),
        min(width, x1 + extra_x),
        min(height, y1 + extra_y),
    )


def display_rgb(rows: list[dict[str, str]]) -> np.ndarray:
    lab = np.asarray(
        [
            [float(row["flower_L_mean"]), float(row["flower_a_mean"]), float(row["flower_b_mean"])]
            for row in rows
        ],
        dtype=float,
    )
    return np.clip(lab2rgb(lab.reshape(-1, 1, 3)).reshape(-1, 3), 0.0, 1.0)


def crop_from_weight(image: Image.Image, weights: np.ndarray) -> Image.Image:
    x0, y0, x1, y1 = mass_box(weights)
    width, height = image.size
    mask_height, mask_width = weights.shape
    box = (
        int(np.floor(x0 * width / mask_width)),
        int(np.floor(y0 * height / mask_height)),
        int(np.ceil(x1 * width / mask_width)),
        int(np.ceil(y1 * height / mask_height)),
    )
    crop = image.convert("RGB").crop(box)
    return ImageOps.fit(crop, (320, 320), method=Image.Resampling.LANCZOS)


def load_photo_bar(
    selected: list[dict[str, str]],
    photo_rows: list[dict[str, str]],
    image_root: Path,
    model_dir: Path,
    provenance_rows: list[dict[str, str]],
) -> tuple[list[Image.Image], list[dict[str, Any]]]:
    by_encounter: dict[str, list[dict[str, str]]] = {}
    for row in photo_rows:
        if row["automated_colour_state_status"] == ADMITTED:
            by_encounter.setdefault(row["encounter_blind_id"], []).append(row)
    chosen: list[dict[str, str]] = []
    images: list[Image.Image] = []
    for encounter in selected:
        candidates = sorted(
            by_encounter.get(encounter["encounter_blind_id"], []),
            key=lambda row: row["photo_blind_id"],
        )
        if not candidates:
            raise RuntimeError("selected admitted encounter has no admitted photograph")
        row = candidates[0]
        path = image_root / row["image_file"]
        if sha256(path) != row["image_sha256"]:
            raise RuntimeError(f"image hash mismatch: {path}")
        chosen.append(row)
        images.append(Image.open(path).convert("RGB"))

    from transformers import CLIPSegForImageSegmentation, CLIPSegProcessor

    processor = CLIPSegProcessor.from_pretrained(
        model_dir, local_files_only=True, use_fast=False
    )
    model = CLIPSegForImageSegmentation.from_pretrained(model_dir, local_files_only=True)
    model.eval()
    crops: list[Image.Image] = []
    batch_size = 4
    for start in range(0, len(images), batch_size):
        batch_images = images[start : start + batch_size]
        logits = inference_batch(model, processor, batch_images)
        for image, maps in zip(batch_images, logits):
            _prompts, ensemble, _negative = candidate_weights(maps)
            crops.append(crop_from_weight(image, ensemble))

    provenance_by_photo = {row["photo_blind_id"]: row for row in provenance_rows}
    audit: list[dict[str, Any]] = []
    for encounter, photo in zip(selected, chosen):
        provenance = provenance_by_photo.get(photo["photo_blind_id"])
        if provenance is None:
            raise RuntimeError("photo attribution is missing")
        audit.append(
            {
                "encounter_blind_id": encounter["encounter_blind_id"],
                "photo_blind_id": photo["photo_blind_id"],
                "image_sha256": photo["image_sha256"],
                "longitude": float(encounter["longitude"]),
                "photo_license": provenance["photo_license"],
                "photo_attribution": provenance["photo_attribution"],
            }
        )
    return crops, audit


def make_atlas_figure(
    rows: list[dict[str, str]], crops: list[Image.Image], output: Path
) -> None:
    admitted = [row for row in rows if row["encounter_status"] == ADMITTED]
    colours = display_rgb(admitted)
    longitude = np.radians([float(row["longitude"]) for row in admitted])
    latitude = np.radians([float(row["latitude"]) for row in admitted])
    figure = plt.figure(figsize=(15, 8.8), constrained_layout=True)
    grid = GridSpec(3, 12, figure=figure, height_ratios=(3.6, 1.0, 1.0))
    axis = figure.add_subplot(grid[0, :], projection="mollweide")
    axis.scatter(longitude, latitude, c=colours, s=25, alpha=0.88, linewidths=0)
    axis.grid(True, color="#d7d7d7", linewidth=0.6)
    axis.set_title(
        "Species-free display of 306 locked model-consensus flower-candidate colours",
        loc="left",
        fontsize=14,
        weight="bold",
    )
    axis.text(
        0.0,
        -0.12,
        "Display only: species labels are hidden here but retained in every inferential step.",
        transform=axis.transAxes,
        fontsize=9,
    )
    for index, crop in enumerate(crops):
        photo_axis = figure.add_subplot(grid[1 + index // 12, index % 12])
        photo_axis.imshow(crop)
        photo_axis.set_xticks([])
        photo_axis.set_yticks([])
        for spine in photo_axis.spines.values():
            spine.set_linewidth(0.45)
            spine.set_color("#777777")
    figure.savefig(output, dpi=300, facecolor="white")
    plt.close(figure)


def make_result_figure(
    development: dict[str, Any], spatial_rows: list[dict[str, str]], output: Path
) -> None:
    development_rows = development["species_results"]
    names = [row["canonical_name"] for row in development_rows]
    shares = np.asarray([float(row["admitted_encounter_share"]) for row in development_rows])
    passed = [row["development_gate_status"] == "pass" for row in development_rows]
    result_names = [row["canonical_name"] for row in spatial_rows]
    rho = np.asarray([float(row["primary_rho"]) for row in spatial_rows])
    q_value = np.asarray([float(row["primary_bh_q"]) for row in spatial_rows])
    contrast = np.asarray([float(row["flower_minus_background_rho"]) for row in spatial_rows])
    contrast_q = np.asarray([float(row["contrast_bh_q"]) for row in spatial_rows])

    figure, axes = plt.subplots(1, 3, figsize=(15, 5.6), constrained_layout=True)
    y = np.arange(len(names))
    axes[0].barh(y, shares, color=["#2878b5" if value else "#b8b8b8" for value in passed])
    axes[0].axvline(0.70, color="#b22222", linestyle="--", linewidth=1.2)
    axes[0].set_yticks(y, labels=names, fontstyle="italic")
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Admitted development encounters")
    axes[0].set_title("a  Location-free image gate", loc="left", weight="bold")
    for position, value in zip(y, shares):
        axes[0].text(value + 0.015, position, f"{100 * value:.1f}%", va="center", fontsize=8)

    result_y = np.arange(len(result_names))
    axes[1].axvline(0, color="#555555", linewidth=0.8)
    axes[1].scatter(rho, result_y, color="#2878b5", s=65, zorder=3)
    axes[1].set_yticks(result_y, labels=result_names, fontstyle="italic")
    axes[1].invert_yaxis()
    axes[1].set_xlim(-0.125, 0.125)
    axes[1].set_xlabel("Primary spatial rho")
    axes[1].set_title("b  Locked random-mark test", loc="left", weight="bold")
    for position, value, q in zip(result_y, rho, q_value):
        axes[1].text(0.12, position, f"q={q:.3f}", ha="right", va="center", fontsize=8)

    axes[2].axvline(0, color="#555555", linewidth=0.8)
    axes[2].scatter(contrast, result_y, color="#d95f02", s=65, zorder=3)
    axes[2].set_yticks(result_y, labels=result_names, fontstyle="italic")
    axes[2].invert_yaxis()
    axes[2].set_xlim(-0.15, 0.15)
    axes[2].set_xlabel("Flower minus background rho")
    axes[2].set_title("c  Negative control", loc="left", weight="bold")
    for position, value, q in zip(result_y, contrast, contrast_q):
        axes[2].text(0.145, position, f"q={q:.3f}", ha="right", va="center", fontsize=8)

    figure.suptitle(
        "Image measurement passed for three species; spatial organization was not detected",
        fontsize=15,
        weight="bold",
    )
    figure.savefig(output, dpi=300, facecolor="white")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-table", type=Path, required=True)
    parser.add_argument("--development-gate", type=Path, required=True)
    parser.add_argument("--spatial-results", type=Path, required=True)
    parser.add_argument("--photo-features", type=Path, required=True)
    parser.add_argument("--image-root", type=Path, required=True)
    parser.add_argument("--photo-provenance", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("docs/figures"))
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_csv(args.analysis_table)
    selected = select_photo_bar_encounters(rows)
    crops, crop_audit = load_photo_bar(
        selected,
        read_csv(args.photo_features),
        args.image_root,
        args.model_dir,
        read_csv(args.photo_provenance),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    atlas_path = args.output_dir / "jbi_inaturalist_automated_colour_pilot_figure1_atlas.png"
    result_path = args.output_dir / "jbi_inaturalist_automated_colour_pilot_figure2_results.png"
    make_atlas_figure(rows, crops, atlas_path)
    development = json.loads(args.development_gate.read_text(encoding="utf-8"))
    make_result_figure(development, read_csv(args.spatial_results), result_path)
    manifest = {
        "status": "complete_automated_colour_pilot_figures",
        "display_species_labels_removed": True,
        "inference_species_conditioned": True,
        "photo_bar_selection": "24 equally spaced positions after sorting admitted encounters by longitude; first admitted photo blind ID per encounter",
        "photo_bar_crop": "central 90% CLIPSeg flower-candidate soft-mask mass plus 12% padding; display only",
        "source_sha256": {
            "analysis_table": sha256(args.analysis_table),
            "development_gate": sha256(args.development_gate),
            "spatial_results": sha256(args.spatial_results),
            "photo_features": sha256(args.photo_features),
            "figure_script": sha256(Path(__file__).resolve()),
        },
        "output_sha256": {
            atlas_path.name: sha256(atlas_path),
            result_path.name: sha256(result_path),
        },
        "displayed_photo_audit": crop_audit,
        "claim_ceiling": "Species-free display and frozen three-species result summary; the map is not inferential evidence and crops are model-consensus flower candidates, not verified flower tissue.",
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
