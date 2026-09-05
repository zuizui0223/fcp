from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.global_barrier_field import (
    edge_jensen_shannon_divergence,
    equal_area_grid_centers,
    node_colour_concentration_permutation_test,
    permute_colour_vectors_within_species,
    prepare_barrier_geometry,
    spherical_midpoint,
)


def _chain_fixture(n_species: int = 40, *, planted: bool) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Repeated eight-photo species chains with fixed geometry.

    In the planted fixture every species changes from colour state A to B across
    the same central geographic edge. In the null fixture every photo has the
    same colour vector, so no geographic colour discontinuity exists.
    """
    longitude = np.array([-140.0, -100.0, -60.0, -20.0, 20.0, 60.0, 100.0, 140.0])
    latitude = np.zeros_like(longitude)
    photo_species: list[str] = []
    colour_rows: list[list[float]] = []
    edge_nodes: list[tuple[int, int]] = []
    edge_species: list[str] = []
    midpoint_lat: list[float] = []
    midpoint_lon: list[float] = []

    for species_id in range(n_species):
        label = f"species_{species_id:03d}"
        offset = len(photo_species)
        photo_species.extend([label] * len(longitude))
        if planted:
            colour_rows.extend([[1.0, 0.0, 0.0, 0.0]] * 4)
            colour_rows.extend([[0.0, 0.0, 1.0, 0.0]] * 4)
        else:
            colour_rows.extend([[0.25, 0.25, 0.25, 0.25]] * len(longitude))

        for local in range(len(longitude) - 1):
            a = offset + local
            b = offset + local + 1
            edge_nodes.append((a, b))
            lat_mid, lon_mid = spherical_midpoint(
                [latitude[local]],
                [longitude[local]],
                [latitude[local + 1]],
                [longitude[local + 1]],
            )
            midpoint_lat.append(float(lat_mid[0]))
            midpoint_lon.append(float(lon_mid[0]))
            edge_species.append(label)

    edge_frame = pd.DataFrame(
        {
            "species": edge_species,
            "midpoint_latitude": midpoint_lat,
            "midpoint_longitude": midpoint_lon,
        }
    )
    return (
        edge_frame,
        np.asarray(colour_rows, dtype=float),
        np.asarray(photo_species),
        np.asarray(edge_nodes, dtype=np.int64),
    )


def test_jsd_respects_complete_soft_colour_vectors():
    colours = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )
    score = edge_jensen_shannon_divergence(colours, [[0, 1], [0, 2]])
    assert score[0] == pytest.approx(0.0)
    assert score[1] == pytest.approx(1.0)


def test_photo_vector_permutation_preserves_species_row_multisets():
    colours = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    species = np.array(["a", "a", "b", "b"])
    shuffled = permute_colour_vectors_within_species(colours, species, rng=np.random.default_rng(7))
    for label in ["a", "b"]:
        idx = np.flatnonzero(species == label)
        before = sorted(map(tuple, colours[idx].tolist()))
        after = sorted(map(tuple, shuffled[idx].tolist()))
        assert before == after


def test_planted_shared_node_colour_boundary_is_detected_by_primary_null():
    edge_frame, colours, photo_species, edge_nodes = _chain_fixture(planted=True)
    geometry = prepare_barrier_geometry(
        edge_frame,
        grid=equal_area_grid_centers(36, 18),
        kernel_km=500.0,
        cutoff_multiplier=3.0,
    )
    result = node_colour_concentration_permutation_test(
        geometry,
        colours,
        photo_species,
        edge_nodes,
        minimum_distinct_species=5,
        permutations=99,
        seed=20260904,
    )
    assert result["primary_rgfca_allowed"] is True
    assert result["null_unit"] == "complete_photo_colour_vector_within_species"
    assert result["observed"].concentration > result["null_mean"]
    assert result["p_upper"] <= 0.02


def test_uniform_colour_fixture_is_not_false_positive_under_node_null():
    edge_frame, colours, photo_species, edge_nodes = _chain_fixture(planted=False)
    geometry = prepare_barrier_geometry(
        edge_frame,
        grid=equal_area_grid_centers(36, 18),
        kernel_km=500.0,
        cutoff_multiplier=3.0,
    )
    result = node_colour_concentration_permutation_test(
        geometry,
        colours,
        photo_species,
        edge_nodes,
        minimum_distinct_species=5,
        permutations=39,
        seed=20260904,
    )
    assert result["observed"].concentration == pytest.approx(0.0, abs=1e-15)
    assert result["null_mean"] == pytest.approx(0.0, abs=1e-15)
    assert result["p_upper"] == pytest.approx(1.0)
