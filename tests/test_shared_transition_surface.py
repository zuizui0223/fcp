import numpy as np
import pytest

from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    build_edge_cell_geometry,
    cell_mean_intensity,
    equal_area_cell_centers,
    equal_area_cell_ids,
    geometry_opportunity_summary,
    spherical_edge_midpoints,
)


def test_equal_area_grid_uses_uniform_longitude_and_sin_latitude_bins():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    ids = equal_area_cell_ids(
        np.array([-60.0, -10.0, 10.0, 60.0]),
        np.array([-170.0, -80.0, 10.0, 100.0]),
        grid,
    )
    np.testing.assert_array_equal(ids, np.array([0, 1, 6, 7]))
    assert grid.n_cells == 8
    assert grid.cell_area_km2 > 0

    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    np.testing.assert_array_equal(cell_id, np.arange(8))
    assert len(latitude) == len(longitude) == 8
    # Equal-area latitude bands are centred at sin(latitude) = -0.5 and +0.5.
    assert latitude[0] == pytest.approx(-30.0)
    assert latitude[4] == pytest.approx(30.0)


def test_spherical_midpoint_handles_dateline_without_wrapping_to_greenwich():
    lat, lon = spherical_edge_midpoints(
        np.array([10.0, 10.0]),
        np.array([179.0, -179.0]),
        np.array([[0, 1]]),
    )
    assert lat[0] == pytest.approx(10.0015, abs=0.01)
    assert abs(abs(lon[0]) - 180.0) < 1e-6


def test_geometry_filters_edges_and_freezes_detectability_without_colour():
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    latitude = np.array([0.0, 0.0, 0.1, 0.1])
    longitude = np.array([0.0, 0.1, 0.0, 0.1])
    edges = np.array([[0, 1], [0, 2], [2, 3], [1, 3]])
    distance = np.array([10.0, 20.0, 30.0, 900.0])

    geometry = build_edge_cell_geometry(
        latitude,
        longitude,
        edges,
        distance,
        grid=grid,
        max_edge_km=100.0,
        min_edges_per_cell=2,
    )

    np.testing.assert_array_equal(geometry.retained_edge_indices, np.array([0, 1, 2]))
    assert len(geometry.retained_edges) == 3
    assert geometry.cell_edge_count.sum() == 3
    assert np.count_nonzero(geometry.detectable) == 1

    # Colour values are supplied only after the geometry is frozen.
    cell_a = cell_mean_intensity(np.array([0.0, 0.5, 1.0]), geometry)
    cell_b = cell_mean_intensity(np.array([1.0, 0.5, 0.0]), geometry)
    np.testing.assert_array_equal(np.isfinite(cell_a), geometry.detectable)
    np.testing.assert_array_equal(np.isfinite(cell_b), geometry.detectable)
    assert cell_a[geometry.detectable][0] == pytest.approx(0.5)
    assert cell_b[geometry.detectable][0] == pytest.approx(0.5)


def test_opportunity_summary_counts_shared_cells_and_species_support():
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    latitude = np.array([0.0, 0.0, 0.1, 0.1])
    longitude = np.array([0.0, 0.1, 0.0, 0.1])
    edges = np.array([[0, 1], [0, 2], [2, 3]])
    distance = np.array([10.0, 20.0, 30.0])

    geometry_a = build_edge_cell_geometry(
        latitude,
        longitude,
        edges,
        distance,
        grid=grid,
        max_edge_km=100.0,
        min_edges_per_cell=2,
    )
    geometry_b = build_edge_cell_geometry(
        latitude,
        longitude,
        edges,
        distance,
        grid=grid,
        max_edge_km=100.0,
        min_edges_per_cell=2,
    )
    summary = geometry_opportunity_summary(
        [geometry_a, geometry_b],
        min_detectable_species=2,
    )
    assert summary["n_cells"] == grid.n_cells
    assert summary["n_cells_A_ge_2"] == 1
    assert summary["max_A"] == 2
    assert summary["species_with_any_shared_opportunity"] == 2
    assert summary["retained_edges_per_species"] == [3, 3]


def test_invalid_geometry_inputs_fail_explicitly():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    with pytest.raises(ValueError, match="removed every edge"):
        build_edge_cell_geometry(
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([[0, 1]]),
            np.array([1000.0]),
            grid=grid,
            max_edge_km=10.0,
            min_edges_per_cell=1,
        )

    with pytest.raises(ValueError, match="match the retained geometry edges"):
        geometry = build_edge_cell_geometry(
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([[0, 1]]),
            np.array([10.0]),
            grid=grid,
            max_edge_km=100.0,
            min_edges_per_cell=1,
        )
        cell_mean_intensity(np.array([0.2, 0.3]), geometry)
