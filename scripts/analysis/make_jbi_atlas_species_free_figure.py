#!/usr/bin/env python3
"""Build the frozen species-free 200-species atlas map and photo bar."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np
from PIL import Image, ImageOps


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
from fcp_pipeline.flower_roi_v4_runtime import (
    FrozenFlowerColourEstimator,
    file_sha256,
    validate_scaleout_authorization,
)
from fcp_pipeline.atlas_measurement import validate_measurement_result_rows
from scripts.data.validate_jbi_atlas_roi_v4_gate_evidence import (
    load_committed_locked_scaleout_result,
)
from scripts.data.validate_jbi_atlas_species_free_visualization_contract import (
    CONTRACT as VISUALIZATION_CONTRACT,
    validate_visualization_contract,
)


ADMITTED = "automated_colour_state_admitted"
LAB_FIELDS = ("flower_L_mean", "flower_a_mean", "flower_b_mean")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def select_photo_bar_rows(
    rows: Sequence[Mapping[str, str]], *, count: int = 48
) -> list[dict[str, str]]:
    """Select frozen equally spaced longitude ranks without colour or taxonomy."""

    if count < 2:
        raise ValueError("photo bar requires at least two images")
    ordered = sorted(
        (dict(row) for row in rows),
        key=lambda row: (
            float(row["longitude"]),
            float(row["latitude"]),
            str(row["measurement_id"]),
        ),
    )
    if len(ordered) < count:
        raise ValueError("not enough admitted atlas rows for the frozen photo bar")
    positions = [
        math.floor(index * (len(ordered) - 1) / (count - 1) + 0.5)
        for index in range(count)
    ]
    if len(set(positions)) != count:
        raise RuntimeError("frozen photo-bar ranks are not unique")
    return [ordered[position] for position in positions]


def central_mask_box(
    mask: np.ndarray, *, quantile: float = 0.05, padding: float = 0.12
) -> tuple[int, int, int, int]:
    values = np.asarray(mask, dtype=bool)
    if values.ndim != 2 or not values.any():
        raise ValueError("photo-bar flower mask is empty")
    y_mass = values.sum(axis=1).astype(float)
    x_mass = values.sum(axis=0).astype(float)
    y_cdf = np.cumsum(y_mass) / y_mass.sum()
    x_cdf = np.cumsum(x_mass) / x_mass.sum()
    y0 = int(np.searchsorted(y_cdf, quantile, side="left"))
    y1 = int(np.searchsorted(y_cdf, 1.0 - quantile, side="left")) + 1
    x0 = int(np.searchsorted(x_cdf, quantile, side="left"))
    x1 = int(np.searchsorted(x_cdf, 1.0 - quantile, side="left")) + 1
    height, width = values.shape
    pad_y = int(round((y1 - y0) * padding))
    pad_x = int(round((x1 - x0) * padding))
    return (
        max(0, x0 - pad_x),
        max(0, y0 - pad_y),
        min(width, x1 + pad_x),
        min(height, y1 + pad_y),
    )


def crop_from_mask(image: Image.Image, mask: np.ndarray) -> tuple[Image.Image, tuple[int, int, int, int]]:
    oriented = ImageOps.exif_transpose(image).convert("RGB")
    if mask.shape != (oriented.height, oriented.width):
        raise ValueError("photo-bar mask and oriented image dimensions differ")
    box = central_mask_box(mask)
    crop = ImageOps.fit(
        oriented.crop(box),
        (320, 320),
        method=Image.Resampling.LANCZOS,
    )
    return crop, box


def display_rgb(rows: Sequence[Mapping[str, str]]) -> np.ndarray:
    from skimage.color import lab2rgb

    lab = np.asarray(
        [[float(row[field]) for field in LAB_FIELDS] for row in rows], dtype=float
    )
    if lab.ndim != 2 or lab.shape[1] != 3 or not np.isfinite(lab).all():
        raise ValueError("atlas display Lab values are invalid")
    return np.clip(lab2rgb(lab.reshape(-1, 1, 3)).reshape(-1, 3), 0.0, 1.0)


def make_figure(
    points: Sequence[Mapping[str, str]],
    crops: Sequence[Image.Image],
    *,
    png_path: Path,
    pdf_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    ordered = sorted(points, key=lambda row: str(row["measurement_id"]))
    colours = display_rgb(ordered)
    longitude = np.radians([float(row["longitude"]) for row in ordered])
    latitude = np.radians([float(row["latitude"]) for row in ordered])
    figure = plt.figure(figsize=(16, 11.2), constrained_layout=True)
    grid = GridSpec(5, 12, figure=figure, height_ratios=(4.8, 1, 1, 1, 1))
    axis = figure.add_subplot(grid[0, :], projection="mollweide")
    axis.scatter(longitude, latitude, c=colours, s=5, alpha=0.55, linewidths=0)
    axis.grid(True, color="#d8d8d8", linewidth=0.45)
    axis.set_title(
        f"Species-free global display of {len(ordered):,} admitted flower-candidate colours",
        loc="left",
        fontsize=14,
        weight="bold",
    )
    axis.text(
        0.0,
        -0.10,
        "Display only: species are hidden here but retained in every transition surface and null.",
        transform=axis.transAxes,
        fontsize=9,
    )
    for index, crop in enumerate(crops):
        photo_axis = figure.add_subplot(grid[1 + index // 12, index % 12])
        photo_axis.imshow(crop)
        photo_axis.set_xticks([])
        photo_axis.set_yticks([])
        for spine in photo_axis.spines.values():
            spine.set_linewidth(0.4)
            spine.set_color("#777777")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        png_path,
        dpi=300,
        facecolor="white",
        metadata={"Software": "FCP frozen atlas visualization v1"},
    )
    figure.savefig(
        pdf_path,
        facecolor="white",
        metadata={
            "Creator": "FCP frozen atlas visualization v1",
            "CreationDate": datetime(2026, 9, 1, tzinfo=timezone.utc),
            "ModDate": datetime(2026, 9, 1, tzinfo=timezone.utc),
        },
    )
    plt.close(figure)


def load_measurements(directory: Path) -> dict[str, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(directory.glob("measurement_shard_*.csv")):
        for raw in read_csv(path):
            row: dict[str, Any] = dict(raw)
            folded = str(row["background_features_available"]).casefold()
            if folded not in {"true", "false"}:
                raise ValueError("invalid measurement background boolean")
            row["background_features_available"] = folded == "true"
            rows.append(row)
    validate_measurement_result_rows(rows)
    by_id = {str(row["measurement_id"]): row for row in rows}
    if len(by_id) != len(rows):
        raise RuntimeError("measurement rows are duplicated")
    return by_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas-points", type=Path, required=True)
    parser.add_argument("--measurement-results-dir", type=Path, required=True)
    parser.add_argument("--sealed-coordinate-key", type=Path, required=True)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--roi-evidence-dir", type=Path, required=True)
    parser.add_argument("--trained-weight", type=Path, required=True)
    parser.add_argument("--efficient-sam-weights-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--torch-threads", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    visualization = json.loads(VISUALIZATION_CONTRACT.read_text(encoding="utf-8"))
    validate_visualization_contract(visualization)
    roi_contract_path = ROOT / visualization["parents"]["roi_v4"]["path"]
    roi_contract = json.loads(roi_contract_path.read_text(encoding="utf-8"))
    validate_roi_v4_contract(roi_contract)
    trained_weight_sha = file_sha256(args.trained_weight)
    locked = load_committed_locked_scaleout_result(args.roi_evidence_dir)
    validate_scaleout_authorization(locked, trained_weight_sha256=trained_weight_sha)

    points = read_csv(args.atlas_points)
    if not points or any(set(row).intersection({"species", "taxon", "cohort_id"}) for row in points):
        raise RuntimeError("atlas display points are empty or expose inferential labels")
    selected = select_photo_bar_rows(
        points, count=int(visualization["photo_bar"]["count"])
    )
    coordinates = {
        row["measurement_id"]: row for row in read_csv(args.sealed_coordinate_key)
    }
    measurements = load_measurements(args.measurement_results_dir)
    estimator = FrozenFlowerColourEstimator(
        args.trained_weight,
        args.efficient_sam_weights_dir,
        roi_contract,
        torch_threads=args.torch_threads,
    )
    crops: list[Image.Image] = []
    audits: list[dict[str, Any]] = []
    crop_dir = args.output_dir / "photo_bar_crops"
    crop_dir.mkdir(parents=True, exist_ok=True)
    for order, point in enumerate(selected, start=1):
        measurement_id = point["measurement_id"]
        coordinate = coordinates.get(measurement_id)
        measured = measurements.get(measurement_id)
        if coordinate is None or measured is None or measured["automated_colour_state_status"] != ADMITTED:
            raise RuntimeError("selected photo-bar row lacks an admitted sealed parent")
        image_path = args.images_dir / f"{measurement_id}.jpg"
        if not image_path.is_file() or sha256(image_path) != measured["image_sha256"]:
            raise RuntimeError("selected photo-bar image identity changed")
        with Image.open(image_path) as image:
            rerun = estimator.measure(image)
            observed_lab = np.asarray([float(measured[field]) for field in LAB_FIELDS])
            rerun_lab = np.asarray([float(rerun[field]) for field in LAB_FIELDS])
            if (
                rerun["automated_colour_state_status"] != ADMITTED
                or not np.allclose(observed_lab, rerun_lab, rtol=0.0, atol=1e-9)
            ):
                raise RuntimeError("selected photo-bar ROI rerun did not reproduce")
            crop, box = crop_from_mask(image, rerun["flower_mask"])
        crop_path = crop_dir / f"photo_bar_{order:02d}.png"
        crop.save(crop_path, format="PNG", optimize=False)
        crops.append(crop)
        audits.append(
            {
                "photo_bar_order": order,
                "measurement_id": measurement_id,
                "image_sha256": sha256(image_path),
                "longitude": float(point["longitude"]),
                "latitude": float(point["latitude"]),
                "photo_license": coordinate["photo_license"],
                "attribution": coordinate["attribution"],
                "crop_box_xyxy": list(box),
                "crop_sha256": sha256(crop_path),
            }
        )

    png = args.output_dir / "jbi_atlas_species_free_map_photo_bar.png"
    pdf = args.output_dir / "jbi_atlas_species_free_map_photo_bar.pdf"
    make_figure(points, crops, png_path=png, pdf_path=pdf)
    manifest = {
        "protocol": visualization["protocol"],
        "status": "complete_species_free_atlas_visualization",
        "display_species_labels_removed": True,
        "inference_remains_species_conditioned": True,
        "map_points": len(points),
        "photo_bar_images": len(crops),
        "source_sha256": {
            "atlas_points": sha256(args.atlas_points),
            "sealed_coordinate_key": sha256(args.sealed_coordinate_key),
            "visualization_contract": sha256(VISUALIZATION_CONTRACT),
            "roi_evidence_manifest": sha256(
                args.roi_evidence_dir / "gate_evidence_manifest.json"
            ),
            "trained_weight": trained_weight_sha,
            "executable": sha256(Path(__file__).resolve()),
        },
        "output_sha256": {png.name: sha256(png), pdf.name: sha256(pdf)},
        "displayed_photo_audit": audits,
        "claim_ceiling": visualization["claim_ceiling"],
    }
    manifest_path = args.output_dir / "species_free_visualization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
