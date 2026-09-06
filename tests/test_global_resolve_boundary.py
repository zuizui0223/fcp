from __future__ import annotations

import numpy as np
import pytest

from fcp_pipeline.global_resolve_boundary import (
    great_circle_samples,
    load_resolve_boundary_index,
    score_resolve_boundary_pairs,
)


def test_great_circle_samples_respect_interval_and_endpoints():
    lat, lon, distance = great_circle_samples(0.0, 0.0, 0.0, 1.0, maximum_interval_km=10.0)
    assert distance > 100.0
    assert len(lat) == int(np.ceil(distance / 10.0)) + 1
    assert lat[0] == pytest.approx(0.0, abs=1e-12)
    assert lon[0] == pytest.approx(0.0, abs=1e-12)
    assert lat[-1] == pytest.approx(0.0, abs=1e-12)
    assert lon[-1] == pytest.approx(1.0, abs=1e-12)
    assert distance / (len(lat) - 1) <= 10.0 + 1e-12


def _write_two_polygon_resolve_fixture(path):
    shapefile = pytest.importorskip("shapefile")
    writer = shapefile.Writer(str(path))
    writer.field("REALM", "C", size=32)
    writer.field("BIOME_NUM", "N", size=8, decimal=0)
    writer.field("ECO_ID", "N", size=8, decimal=0)
    writer.poly([[[-2.0, -1.0], [0.0, -1.0], [0.0, 1.0], [-2.0, 1.0], [-2.0, -1.0]]])
    writer.record("RealmA", 1, 1)
    writer.poly([[[0.0, -1.0], [2.0, -1.0], [2.0, 1.0], [0.0, 1.0], [0.0, -1.0]]])
    writer.record("RealmB", 2, 2)
    writer.close()


def test_resolve_boundary_score_and_external_missingness(tmp_path):
    pytest.importorskip("shapely")
    shp = tmp_path / "fixture.shp"
    _write_two_polygon_resolve_fixture(shp)
    index = load_resolve_boundary_index(shp)

    # Row 0 -> row 1 crosses exactly one realm+biome+ecoregion boundary = 3+2+1.
    # Row 0 -> row 2 stays in the left polygon = 0.
    # Row 1 -> row 3 leaves the supplied terrestrial polygons and must be NaN, not 0.
    latitude = np.asarray([0.0, 0.0, 0.25, 0.0])
    longitude = np.asarray([-1.0, 1.0, -1.25, 3.0])
    pairs = np.asarray([[0, 1], [0, 2], [1, 3]], dtype=np.int64)
    scores = score_resolve_boundary_pairs(
        pairs,
        latitude,
        longitude,
        index,
        maximum_interval_km=20.0,
        maximum_points_per_batch=1000,
    )
    assert scores[0] == pytest.approx(6.0)
    assert scores[1] == pytest.approx(0.0)
    assert np.isnan(scores[2])
