from __future__ import annotations

import numpy as np

from fcp_pipeline.global_g3 import species_distance_colour_rho


def test_constant_colour_divergence_is_structural_zero_not_missing():
    latitude = np.zeros(6)
    longitude = np.array([-50.0, -30.0, -10.0, 10.0, 30.0, 50.0])
    colours = np.tile(np.array([[1.0, 0.0, 0.0, 0.0]]), (6, 1))
    assert species_distance_colour_rho(latitude, longitude, colours) == 0.0


def test_constant_geographic_distance_information_is_not_evaluable():
    latitude = np.zeros(5)
    longitude = np.zeros(5)
    colours = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.8, 0.2, 0.0, 0.0],
            [0.5, 0.5, 0.0, 0.0],
            [0.2, 0.8, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
        ]
    )
    assert np.isnan(species_distance_colour_rho(latitude, longitude, colours))


def test_broad_geographic_colour_split_gives_positive_rho():
    latitude = np.zeros(8)
    longitude = np.array([-70.0, -55.0, -40.0, -25.0, 25.0, 40.0, 55.0, 70.0])
    colours = np.array(
        [[1.0, 0.0, 0.0, 0.0]] * 4 + [[0.0, 0.0, 1.0, 0.0]] * 4,
        dtype=float,
    )
    rho = species_distance_colour_rho(latitude, longitude, colours)
    assert np.isfinite(rho)
    assert rho > 0.5
