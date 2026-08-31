#!/usr/bin/env python3
"""Area-weight WorldClim and RESOLVE onto every frozen equal-area grid."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import (
    validate_environment_contract,
    weighted_cell_means,
    weighted_dominant_labels,
)
from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    equal_area_cell_centers,
    equal_area_cell_ids,
)


BIO_NUMBERS = (1, 4, 12, 15)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty grid: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worldclim-archive", type=Path, required=True)
    parser.add_argument("--worldclim-dir", type=Path, required=True)
    parser.add_argument("--ecoregion-archive", type=Path, required=True)
    parser.add_argument("--ecoregion-shapefile", type=Path, required=True)
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
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    environment = json.loads(args.environment_contract.read_text(encoding="utf-8"))
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    validate_environment_contract(environment)
    expected_climate = environment["sources"]["macroclimate"]["archive_sha256"]
    expected_ecoregion = environment["sources"]["ecoregion"]["archive_sha256"]
    if sha256(args.worldclim_archive) != expected_climate:
        raise RuntimeError("WorldClim archive SHA-256 mismatch")
    if sha256(args.ecoregion_archive) != expected_ecoregion:
        raise RuntimeError("RESOLVE archive SHA-256 mismatch")

    import rasterio
    from rasterio.features import rasterize
    import shapefile

    raster_paths = [
        args.worldclim_dir / f"wc2.1_10m_bio_{number}.tif"
        for number in BIO_NUMBERS
    ]
    datasets = [rasterio.open(path) for path in raster_paths]
    reference = datasets[0]
    if any(
        dataset.shape != reference.shape
        or dataset.transform != reference.transform
        or dataset.crs != reference.crs
        for dataset in datasets[1:]
    ):
        raise RuntimeError("WorldClim rasters do not share one grid")
    arrays = [dataset.read(1, masked=True) for dataset in datasets]
    climate_valid = ~np.logical_or.reduce([np.ma.getmaskarray(array) for array in arrays])

    reader = shapefile.Reader(str(args.ecoregion_shapefile), encoding="cp1252")
    attributes: dict[int, dict[str, str]] = {}
    shapes: list[tuple[dict[str, Any], int]] = []
    for item in reader.iterShapeRecords():
        record = item.record.as_dict()
        eco_id = int(record["ECO_ID"])
        current = {
            "ecoregion": str(record["ECO_NAME"]),
            "biome": str(record["BIOME_NAME"]),
            "realm": str(record["REALM"]),
        }
        previous = attributes.setdefault(eco_id, current)
        if previous != current:
            raise RuntimeError(f"ECO_ID {eco_id} has inconsistent attributes")
        shapes.append((item.shape.__geo_interface__, eco_id))
    eco_raster = rasterize(
        shapes,
        out_shape=reference.shape,
        transform=reference.transform,
        fill=0,
        all_touched=False,
        dtype="int32",
    )
    valid = climate_valid & (eco_raster > 0)
    source_row, source_column = np.nonzero(valid)
    longitude = reference.transform.c + (source_column + 0.5) * reference.transform.a
    latitude = reference.transform.f + (source_row + 0.5) * reference.transform.e
    area_weight = np.cos(np.deg2rad(latitude))
    climate_values = np.column_stack(
        [np.asarray(array.data[source_row, source_column], dtype=float) for array in arrays]
    )
    eco_values = eco_raster[source_row, source_column]

    output_files: dict[str, str] = {}
    scale_results: list[dict[str, Any]] = []
    candidates = {
        int(row["scale_km"]): (int(row["n_lon"]), int(row["n_sinlat"]))
        for row in atlas["geometry_only_scale_selection"]["candidates"]
    }
    for scale in environment["grid"]["scales_km"]:
        n_lon, n_sinlat = candidates[scale]
        grid = EqualAreaGrid(n_lon=n_lon, n_sinlat=n_sinlat)
        source_cell = equal_area_cell_ids(latitude, longitude, grid)
        climate_mean = weighted_cell_means(
            source_cell,
            area_weight,
            climate_values,
            n_cells=grid.n_cells,
        )
        dominant = weighted_dominant_labels(source_cell, area_weight, eco_values)
        cell_ids, center_latitude, center_longitude = equal_area_cell_centers(grid)
        selected = [
            int(cell)
            for cell in cell_ids
            if int(cell) in dominant and np.isfinite(climate_mean[int(cell)]).all()
        ]
        rows: list[dict[str, Any]] = []
        for cell in selected:
            label = attributes[dominant[cell]]
            rows.append(
                {
                    "scale_km": scale,
                    "cell_id": cell,
                    "latitude": float(center_latitude[cell]),
                    "longitude": float(center_longitude[cell]),
                    "bio1": float(climate_mean[cell, 0]),
                    "bio4": float(climate_mean[cell, 1]),
                    "bio12": float(climate_mean[cell, 2]),
                    "bio15": float(climate_mean[cell, 3]),
                    "realm": label["realm"],
                    "biome": label["biome"],
                    "ecoregion": label["ecoregion"],
                    "dominant_eco_id": dominant[cell],
                }
            )
        path = args.output_dir / f"climate_ecoregion_grid_{scale}km.csv"
        write_csv(path, rows)
        output_files[path.name] = sha256(path)
        scale_results.append(
            {
                "scale_km": scale,
                "global_cells": grid.n_cells,
                "terrestrial_cells": len(rows),
                "source_pixels": len(source_cell),
            }
        )

    source_files = [
        *raster_paths,
        args.ecoregion_shapefile,
        args.ecoregion_shapefile.with_suffix(".dbf"),
        args.ecoregion_shapefile.with_suffix(".shx"),
        args.ecoregion_shapefile.with_suffix(".prj"),
    ]
    manifest = {
        "status": "pass_climate_ecoregion_source_freeze",
        "protocol": environment["protocol"],
        "scaleout_colour_opened": False,
        "available_primary_families": ["macroclimate", "ecoregion"],
        "not_yet_frozen_primary_families": ["terrain", "land_cover"],
        "archive_sha256": {
            args.worldclim_archive.name: sha256(args.worldclim_archive),
            args.ecoregion_archive.name: sha256(args.ecoregion_archive),
        },
        "extracted_source_sha256": {path.name: sha256(path) for path in source_files},
        "output_sha256": output_files,
        "scale_results": scale_results,
        "source_pixels_with_climate_and_ecoregion": int(valid.sum()),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "rasterio": rasterio.__version__,
            "pyshp": shapefile.__version__,
        },
        "claim_ceiling": "Independent pre-colour climate and ecoregion grids only; terrain, land cover and flower concordance remain unopened.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "climate_ecoregion_source_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    for dataset in datasets:
        dataset.close()


if __name__ == "__main__":
    main()
