"""Small synthetic-grid tests for the outcome-blind SoilGrids proxy search."""
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))

from run_polymorphism_soil_neighbourhood_coverage_20261008 import (
    nearest_common_valid, valid_grid,
)
from run_polymorphism_full_gradient_partition_20261007 import SOIL_RAW, SOIL_DEPTHS


def test_nearest_valid_within_ten_km():
    grid = np.zeros((7, 7), dtype=bool)
    grid[0, 0] = True
    ref = (from_origin(0.0, 0.2, 0.05, 0.05), CRS.from_epsg(4326), 7, 7)
    lon = np.array([0.025, 0.075, 0.275, np.nan])
    lat = np.array([0.175, 0.175, 0.025, 0.175])
    out = nearest_common_valid(lon, lat, grid, ref)
    assert out["status"].tolist() == [
        "original_valid",
        "proxy_within_10km",
        "no_valid_pixel_within_10km",
        "invalid_coordinate",
    ]
    assert out["substitutions"] == 1
    assert out["distance_km"][0] == 0
    assert 5.0 < out["distance_km"][1] < 6.0
    assert np.isnan(out["distance_km"][2])
    assert out["row"][1] == 0 and out["col"][1] == 0


def test_common_grid_rejects_missing_and_unphysical_awc(tmp_path):
    affine = from_origin(0.0, 0.2, 0.05, 0.05)
    for prop in SOIL_RAW:
        for depth, _ in SOIL_DEPTHS:
            file = tmp_path / prop / f"{prop}_{depth}_mean_5000.tif"
            file.parent.mkdir(parents=True, exist_ok=True)
            a = np.full((4, 4), 50, dtype=np.int16)
            if prop == "cec" and depth == "5-15cm":
                a[0, 2] = -32768
            if prop == "wv1500":
                a[2, 2] = 90
            with rasterio.open(
                file, "w", driver="GTiff", height=4, width=4,
                count=1, dtype="int16", nodata=-32768,
                crs="EPSG:4326", transform=affine,
            ) as ds:
                ds.write(a, 1)
    common, ref = valid_grid(tmp_path)
    assert ref[2:] == (4, 4)
    assert common.shape == (4, 4)
    assert bool(common[1, 1])
    assert not bool(common[0, 2])
    assert not bool(common[2, 2])
    assert int(common.sum()) == 14
