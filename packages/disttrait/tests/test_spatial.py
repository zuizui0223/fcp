import numpy as np
import pytest

from disttrait import (
    distribution_spatial_association,
    great_circle_pairwise_km,
    jensen_shannon_pairwise,
    matched_difference_spatial_rho,
    partial_distribution_spatial_association,
    spatial_permutation_null,
    spatial_rho,
)


def test_pairwise_distance_and_jsd_shapes():
    lat = [0.0, 0.0, 0.0, 0.0]
    lon = [0.0, 1.0, 2.0, 3.0]
    traits = np.array([
        [1.0, 0.0],
        [0.8, 0.2],
        [0.2, 0.8],
        [0.0, 1.0],
    ])
    assert great_circle_pairwise_km(lat, lon).shape == (6,)
    assert jensen_shannon_pairwise(traits).shape == (6,)
    assert spatial_rho(lat, lon, traits) > 0.8


def test_spatial_null_is_deterministic_and_geometry_preserving():
    lat = [0.0, 0.0, 0.0, 0.0]
    lon = [0.0, 1.0, 2.0, 3.0]
    traits = np.array([
        [1.0, 0.0],
        [0.8, 0.2],
        [0.2, 0.8],
        [0.0, 1.0],
    ])
    obs1, null1 = spatial_permutation_null(lat, lon, traits, n_permutations=19, seed=7, key="sp")
    obs2, null2 = spatial_permutation_null(lat, lon, traits, n_permutations=19, seed=7, key="sp")
    assert obs1 == obs2
    assert np.array_equal(null1, null2)
    assert len(null1) == 19


def test_matched_background_uses_pairwise_difference_not_rho_difference():
    lat = [0.0, 0.0, 0.0, 0.0]
    lon = [0.0, 1.0, 2.0, 3.0]
    focal = np.array([[1, 0], [0.8, 0.2], [0.2, 0.8], [0, 1]], dtype=float)
    background = np.tile([0.5, 0.5], (4, 1))
    assert matched_difference_spatial_rho(lat, lon, focal, background) == spatial_rho(lat, lon, focal)


def test_distribution_spatial_association_and_partial_null():
    d = np.array([0.0, 0.2, 0.4, 0.6, 0.8])
    observed = np.array([-0.2, -0.1, 0.0, 0.2, 0.4])
    null = np.column_stack([
        observed[::-1],
        np.roll(observed, 1),
        np.roll(observed, 2),
    ])
    res = distribution_spatial_association(d, observed, spatial_null=null)
    assert res.observed > 0.9
    assert res.p_upper is not None
    partial = partial_distribution_spatial_association(
        d,
        observed,
        controls=[np.array([1, 2, 1, 2, 1], dtype=float)],
        spatial_null=null,
    )
    assert np.isfinite(partial.observed)
    assert partial.null is not None


def test_spatial_edge_cases_match_frozen_fcp_semantics():
    traits_constant = np.tile([0.5, 0.5], (4, 1))
    lat = [0.0, 0.0, 0.0, 0.0]
    lon = [0.0, 1.0, 2.0, 3.0]
    assert spatial_rho(lat, lon, traits_constant) == 0.0

    with pytest.raises(ValueError, match="not_evaluable_pair_geometry"):
        spatial_permutation_null(
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
            np.array([[1.0, 0.0], [0.8, 0.2], [0.2, 0.8], [0.0, 1.0]]),
            n_permutations=3,
        )
