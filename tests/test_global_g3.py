from __future__ import annotations

import numpy as np
import pandas as pd

from fcp_pipeline.global_g3 import run_g3_prevalence, species_distance_colour_rho


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


def _small_pool(n_species: int = 6, photos_per_species: int = 8) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    longitudes = np.linspace(-70.0, 70.0, photos_per_species)
    for s in range(n_species):
        species = f"species_{s:02d}"
        for p, lon in enumerate(longitudes):
            x = p / (photos_per_species - 1)
            rows.append(
                {
                    "photo_id": 50_000 + s * 100 + p,
                    "species": species,
                    "latitude": float(s - n_species / 2),
                    "longitude": float(lon),
                    "colour_white": float(1.0 - x),
                    "colour_yellow_orange": 0.0,
                    "colour_red_pink": float(x),
                    "colour_blue_purple": 0.0,
                }
            )
    return pd.DataFrame(rows)


def test_repeated_g3_schedule_is_deterministic_and_heterogeneity_is_finite():
    pool = _small_pool()
    kwargs = dict(
        n_outer=6,
        species_per_outer=4,
        photos_per_species=5,
        minimum_pool_photos_per_species=7,
        species_seed=123,
        photo_master_seed=456,
        variance_floor=1e-6,
    )
    a = run_g3_prevalence(pool, **kwargs)
    b = run_g3_prevalence(pool, **kwargs)
    pd.testing.assert_frame_equal(a.outer, b.outer)
    pd.testing.assert_frame_equal(a.species, b.species)
    assert len(a.outer) == 6
    assert len(a.species) == 6
    assert (a.outer["scheduled_species"] == 4).all()
    assert (a.outer["evaluable_species"] == 4).all()
    assert np.isfinite(a.outer[["mean_rho", "median_rho", "positive_fraction"]].to_numpy()).all()
    assert np.isfinite(a.tau2_fisher_z)
    assert a.tau2_species_used >= 2
    assert np.isclose(a.tau2_fisher_z, b.tau2_fisher_z)
