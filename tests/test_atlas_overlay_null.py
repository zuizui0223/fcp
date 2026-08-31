from __future__ import annotations

import numpy as np

from fcp_pipeline.atlas_overlay_null import (
    build_moran_sign_basis,
    equal_area_rook_adjacency,
    geographic_design,
    moran_quadratic,
    residual_coefficients,
    spectral_family_test,
)


def basis():
    n_lon, n_sinlat = 8, 4
    cell_ids = np.arange(n_lon * n_sinlat)
    row, column = np.divmod(cell_ids, n_lon)
    latitude = np.rad2deg(np.arcsin(-1.0 + (row + 0.5) * 2.0 / n_sinlat))
    longitude = -180.0 + (column + 0.5) * 360.0 / n_lon
    adjacency = equal_area_rook_adjacency(
        cell_ids, n_lon=n_lon, n_sinlat=n_sinlat
    )
    return build_moran_sign_basis(
        adjacency,
        np.linspace(1.0, 3.0, cell_ids.size),
        geographic_design(latitude, longitude),
    )


def test_sign_randomization_preserves_variance_trends_and_moran_energy() -> None:
    item = basis()
    values = np.sin(np.linspace(0, 4 * np.pi, item.sqrt_weights.size))
    coefficients = residual_coefficients(values, item)
    signs = np.where(np.arange(coefficients.size) % 2, -1.0, 1.0)
    randomized = coefficients * signs
    assert np.linalg.norm(randomized) == np.linalg.norm(coefficients)
    assert moran_quadratic(randomized, item) == moran_quadratic(coefficients, item)


def test_spectral_family_test_detects_alignment_and_uses_family_max() -> None:
    item = basis()
    rng = np.random.default_rng(8)
    aligned = rng.normal(size=item.sqrt_weights.size)
    unrelated = rng.normal(size=item.sqrt_weights.size)
    result = spectral_family_test(
        aligned,
        {"aligned": aligned, "unrelated": unrelated},
        item,
        randomizations=199,
        rng=np.random.default_rng(10),
    )
    assert result["family_statistic"] == result["observed_by_overlay"]["aligned"]
    assert result["p_value"] <= 0.01
    assert len(result["null_statistics"]) == 199


def test_equal_area_adjacency_wraps_longitude() -> None:
    adjacency = equal_area_rook_adjacency(
        np.arange(12), n_lon=4, n_sinlat=3
    )
    assert adjacency[0, 3] == 1.0
    assert adjacency[0, 4] == 1.0


def test_equal_area_adjacency_repairs_isolated_subset_cells() -> None:
    adjacency = equal_area_rook_adjacency(
        np.array([0, 6, 15]), n_lon=8, n_sinlat=4
    )
    assert np.all(adjacency.sum(axis=1) >= 1)
    assert np.allclose(adjacency, adjacency.T)
