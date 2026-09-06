from __future__ import annotations

import numpy as np

from fcp_pipeline.global_colour_zones import (
    connected_components,
    extract_persistent_colour_zones,
    persistent_hotspots,
)


def test_persistent_hotspot_requires_repeated_membership_not_one_large_peak():
    rng = np.random.default_rng(5)
    fields = rng.normal(0.0, 0.05, size=(200, 100))
    fields[:, 40:44] += 5.0
    # Cell 80 is extreme only in 50/200 resamples and must not enter the 60% seed.
    fields[:50, 80] += 10.0
    persistence, _, seed, opportunities = persistent_hotspots(
        fields,
        hotspot_quantile=0.90,
        minimum_persistence=0.60,
        minimum_evaluable_resamples=100,
    )
    assert np.all(seed[40:44])
    assert not seed[80]
    assert persistence[80] < 0.60
    assert np.all(opportunities == 200)


def test_missing_support_is_not_counted_as_non_hotspot():
    fields = np.full((200, 20), np.nan)
    fields[:120, :10] = 0.0
    fields[:120, 3] = 10.0
    persistence, _, seed, opportunities = persistent_hotspots(
        fields,
        hotspot_quantile=0.90,
        minimum_persistence=0.60,
        minimum_evaluable_resamples=100,
    )
    assert opportunities[3] == 120
    assert persistence[3] == 1.0
    assert seed[3]
    assert opportunities[15] == 0
    assert np.isnan(persistence[15])


def test_dateline_cells_are_connected_under_longitude_wrap():
    # 4 longitudes x 2 sin-lat rows. Cell 0 and cell 3 sit at opposite array
    # edges but are geographic neighbours across the dateline.
    mask = np.zeros(8, dtype=bool)
    mask[[0, 3]] = True
    components = connected_components(mask, n_lon=4, n_sinlat=2)
    assert components == ((0, 3),)


def test_extract_zones_assigns_neutral_ids_by_integrated_intensity():
    rng = np.random.default_rng(42)
    n_lon, n_sinlat = 10, 5
    n_cells = n_lon * n_sinlat
    fields = rng.normal(0.0, 0.02, size=(200, n_cells))
    # Two persistent three-cell zones far apart; the first is stronger and must
    # become Z01 without reference to any named geography.
    zone_a = [11, 12, 13]
    zone_b = [36, 37, 38]
    fields[:, zone_a] += 8.0
    fields[:, zone_b] += 5.0
    result = extract_persistent_colour_zones(
        fields,
        n_lon=n_lon,
        n_sinlat=n_sinlat,
        opportunity=np.ones(n_cells),
        hotspot_quantile=0.90,
        minimum_persistence=0.60,
        minimum_evaluable_resamples=100,
        minimum_zone_cells=3,
    )
    assert len(result.zones) == 2
    assert result.zones[0].zone_id == "Z01"
    assert set(result.zones[0].cell_indices) == set(zone_a)
    assert result.zones[1].zone_id == "Z02"
    assert set(result.zones[1].cell_indices) == set(zone_b)


def test_components_below_minimum_zone_size_are_removed():
    rng = np.random.default_rng(4)
    n_lon, n_sinlat = 10, 5
    n_cells = n_lon * n_sinlat
    fields = rng.normal(0.0, 0.01, size=(200, n_cells))
    fields[:, [10, 11]] += 10.0
    result = extract_persistent_colour_zones(
        fields,
        n_lon=n_lon,
        n_sinlat=n_sinlat,
        opportunity=np.ones(n_cells),
        hotspot_quantile=0.90,
        minimum_persistence=0.60,
        minimum_evaluable_resamples=100,
        minimum_zone_cells=3,
    )
    assert result.zones == ()
