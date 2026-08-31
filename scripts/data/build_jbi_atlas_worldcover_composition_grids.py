#!/usr/bin/env python3
"""Systematically sample all frozen ESA WorldCover COGs into atlas grids."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
from typing import Any, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import (
    LAND_COVER_CODES,
    validate_environment_contract,
)
from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    equal_area_cell_centers,
    equal_area_cell_ids,
)


PRIMARY_SAMPLES = 60
SENSITIVITY_SAMPLES = 30


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty composition: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def valid_cache(path: Path, row: Mapping[str, str]) -> bool:
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as cached:
            return (
                str(cached["etag"].item()) == row["etag"]
                and int(cached["size_bytes"].item()) == int(row["size_bytes"])
                and cached["primary"].shape == (PRIMARY_SAMPLES, PRIMARY_SAMPLES)
                and cached["sensitivity"].shape
                == (SENSITIVITY_SAMPLES, SENSITIVITY_SAMPLES)
            )
    except (OSError, ValueError, KeyError):
        return False


def sample_tile(row: Mapping[str, str], cache_dir: Path, retries: int) -> Path:
    import rasterio
    from rasterio.enums import Resampling

    name = Path(row["key"]).stem + ".npz"
    path = cache_dir / name
    if valid_cache(path, row):
        return path
    cache_dir.mkdir(parents=True, exist_ok=True)
    error: Exception | None = None
    for attempt in range(retries):
        try:
            with rasterio.Env(
                GDAL_HTTP_MAX_RETRY="4",
                GDAL_HTTP_RETRY_DELAY="1",
                GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
            ):
                with rasterio.open(row["https_url"]) as dataset:
                    if (
                        dataset.count != 1
                        or dataset.nodata != 0
                        or dataset.crs.to_epsg() != 4326
                        or not math.isclose(dataset.bounds.right - dataset.bounds.left, 3.0)
                        or not math.isclose(dataset.bounds.top - dataset.bounds.bottom, 3.0)
                    ):
                        raise RuntimeError("WorldCover COG structure changed")
                    primary = dataset.read(
                        1,
                        out_shape=(PRIMARY_SAMPLES, PRIMARY_SAMPLES),
                        resampling=Resampling.nearest,
                    )
                    sensitivity = dataset.read(
                        1,
                        out_shape=(SENSITIVITY_SAMPLES, SENSITIVITY_SAMPLES),
                        resampling=Resampling.nearest,
                    )
                    bounds = np.asarray(dataset.bounds, dtype=float)
            observed = set(np.unique(primary)) | set(np.unique(sensitivity))
            if not observed.issubset({0, *LAND_COVER_CODES}):
                raise RuntimeError(f"unexpected WorldCover classes: {sorted(observed)}")
            partial = path.with_suffix(".partial")
            with partial.open("wb") as handle:
                np.savez_compressed(
                    handle,
                    primary=primary,
                    sensitivity=sensitivity,
                    bounds=bounds,
                    etag=np.asarray(row["etag"]),
                    size_bytes=np.asarray(int(row["size_bytes"])),
                )
            partial.replace(path)
            return path
        except Exception as exc:
            error = exc
            if attempt + 1 < retries:
                time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"{row['key']}: {error}")


def sample_coordinates(bounds: np.ndarray, size: int) -> tuple[np.ndarray, np.ndarray]:
    left, bottom, right, top = (float(value) for value in bounds)
    longitude = left + (np.arange(size, dtype=float) + 0.5) * (right - left) / size
    latitude = top - (np.arange(size, dtype=float) + 0.5) * (top - bottom) / size
    lon_grid, lat_grid = np.meshgrid(longitude, latitude)
    return lat_grid.reshape(-1), lon_grid.reshape(-1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--inventory-manifest", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--environment-contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json"),
    )
    parser.add_argument(
        "--atlas-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--maximum-tiles", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    environment = json.loads(args.environment_contract.read_text(encoding="utf-8"))
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    inventory_manifest = json.loads(args.inventory_manifest.read_text(encoding="utf-8"))
    validate_environment_contract(environment)
    if (
        inventory_manifest.get("status") != "pass_worldcover_inventory_freeze"
        or sha256(args.inventory) != inventory_manifest.get("inventory_sha256")
    ):
        raise RuntimeError("WorldCover inventory identity mismatch")
    rows = read_csv(args.inventory)
    formal = args.maximum_tiles <= 0 or args.maximum_tiles >= len(rows)
    selected = rows if formal else rows[: args.maximum_tiles]
    failures: list[str] = []
    completed = 0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(sample_tile, row, args.cache_dir, args.retries): row
            for row in selected
        }
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                failures.append(str(exc))
            completed += 1
            if completed % 25 == 0 or completed == len(selected):
                print(
                    f"worldcover_tiles={completed}/{len(selected)} failures={len(failures)}",
                    flush=True,
                )
    if failures:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "worldcover_tile_failures.json").write_text(
            json.dumps({"failures": failures}, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        raise SystemExit(2)

    candidates = {
        int(row["scale_km"]): (int(row["n_lon"]), int(row["n_sinlat"]))
        for row in atlas["geometry_only_scale_selection"]["candidates"]
    }
    accumulators: dict[tuple[int, int], dict[str, np.ndarray]] = {}
    for scale, (n_lon, n_sinlat) in candidates.items():
        grid = EqualAreaGrid(n_lon=n_lon, n_sinlat=n_sinlat)
        for sample_size in (PRIMARY_SAMPLES, SENSITIVITY_SAMPLES):
            accumulators[(scale, sample_size)] = {
                "weighted": np.zeros((grid.n_cells, len(LAND_COVER_CODES))),
                "counts": np.zeros(grid.n_cells, dtype=int),
            }

    for row in selected:
        path = args.cache_dir / (Path(row["key"]).stem + ".npz")
        with np.load(path, allow_pickle=False) as cached:
            bounds = cached["bounds"]
            for sample_size, field in (
                (PRIMARY_SAMPLES, "primary"),
                (SENSITIVITY_SAMPLES, "sensitivity"),
            ):
                classes = cached[field].reshape(-1).astype(int)
                latitude, longitude = sample_coordinates(bounds, sample_size)
                valid = classes != 0
                classes = classes[valid]
                latitude = latitude[valid]
                longitude = longitude[valid]
                weights = np.cos(np.deg2rad(latitude))
                for scale, (n_lon, n_sinlat) in candidates.items():
                    grid = EqualAreaGrid(n_lon=n_lon, n_sinlat=n_sinlat)
                    cell_ids = equal_area_cell_ids(latitude, longitude, grid)
                    target = accumulators[(scale, sample_size)]
                    np.add.at(target["counts"], cell_ids, 1)
                    for column, code in enumerate(LAND_COVER_CODES):
                        keep = classes == code
                        np.add.at(
                            target["weighted"][:, column],
                            cell_ids[keep],
                            weights[keep],
                        )

    outputs: dict[str, str] = {}
    scale_results: list[dict[str, Any]] = []
    minimum = {100: 25, 250: 100, 500: 400}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for scale, (n_lon, n_sinlat) in candidates.items():
        grid = EqualAreaGrid(n_lon=n_lon, n_sinlat=n_sinlat)
        cell_ids, latitude, longitude = equal_area_cell_centers(grid)
        for sample_size, label in (
            (PRIMARY_SAMPLES, "primary"),
            (SENSITIVITY_SAMPLES, "sensitivity"),
        ):
            target = accumulators[(scale, sample_size)]
            denominator = target["weighted"].sum(axis=1)
            composition = np.full_like(target["weighted"], np.nan)
            present = denominator > 0
            composition[present] = target["weighted"][present] / denominator[present, None]
            keep = target["counts"] >= minimum[scale]
            output_rows = []
            for cell in cell_ids[keep]:
                cell = int(cell)
                output_rows.append(
                    {
                        "scale_km": scale,
                        "cell_id": cell,
                        "latitude": float(latitude[cell]),
                        "longitude": float(longitude[cell]),
                        "valid_samples": int(target["counts"][cell]),
                        **{
                            f"worldcover_{code}": float(composition[cell, index])
                            for index, code in enumerate(LAND_COVER_CODES)
                        },
                    }
                )
            path = args.output_dir / f"worldcover_composition_{label}_{scale}km.csv"
            write_csv(path, output_rows)
            outputs[path.name] = sha256(path)
            scale_results.append(
                {
                    "scale_km": scale,
                    "sampling": label,
                    "sample_points_per_source_tile_axis": sample_size,
                    "cells_passing_minimum_samples": len(output_rows),
                    "minimum_valid_samples": minimum[scale],
                }
            )

    import rasterio

    manifest = {
        "status": (
            "pass_worldcover_composition_freeze"
            if formal
            else "smoke_only_worldcover_composition"
        ),
        "formal_complete_inventory": formal,
        "inventory_sha256": sha256(args.inventory),
        "inventory_objects": len(rows),
        "tiles_processed": len(selected),
        "tile_failures": 0,
        "primary_samples_per_tile_axis": PRIMARY_SAMPLES,
        "sensitivity_samples_per_tile_axis": SENSITIVITY_SAMPLES,
        "output_sha256": outputs,
        "scale_results": scale_results,
        "scaleout_colour_opened": False,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "rasterio": rasterio.__version__,
        },
        "claim_ceiling": "Pre-colour systematic WorldCover composition only; not exhaustive 10 m aggregation and no flower concordance.",
    }
    (args.output_dir / "worldcover_composition_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
