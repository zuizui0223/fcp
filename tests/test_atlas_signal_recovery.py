from __future__ import annotations

import numpy as np

from fcp_pipeline.atlas_signal_recovery import (
    AtlasSpeciesGeometry,
    atlas_shared_concentration,
    batched_permutation_p_value,
    permutation_p_value,
    synthetic_colour_vectors,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, build_edge_cell_geometry
from fcp_pipeline.spatial_graph import spherical_knn_edges


def geometries() -> list[AtlasSpeciesGeometry]:
    grid = EqualAreaGrid(n_lon=12, n_sinlat=6)
    base_lat = np.array([-2.0, -1.0, -0.5, 0.5, 1.0, 2.0, -1.5, 1.5])
    base_lon = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0, -0.2, 0.2])
    out = []
    for index in range(4):
        lat = base_lat + 0.05 * index
        lon = base_lon + 0.05 * index
        edges, distance = spherical_knn_edges(lat, lon, k=3)
        geometry = build_edge_cell_geometry(
            lat,
            lon,
            edges,
            distance,
            grid=grid,
            max_edge_km=1000,
            min_edges_per_cell=1,
        )
        out.append(AtlasSpeciesGeometry(f"species-{index}", lat, lon, geometry))
    return out


def test_shared_signal_raises_exact_surface_concentration() -> None:
    items = geometries()
    null = synthetic_colour_vectors(
        items,
        effect_size=0.0,
        scenario="null_stationary",
        shared_anchor_latitude=0.0,
        shared_anchor_longitude=0.0,
        rng=np.random.default_rng(4),
    )
    signal = synthetic_colour_vectors(
        items,
        effect_size=8.0,
        scenario="shared_geographic_boundary",
        shared_anchor_latitude=0.0,
        shared_anchor_longitude=0.0,
        rng=np.random.default_rng(4),
    )
    null_stat = atlas_shared_concentration(null, items, min_detectable_species=3)
    signal_stat = atlas_shared_concentration(signal, items, min_detectable_species=3)
    assert signal_stat > null_stat


def test_signal_recovery_permutation_p_value_is_nonzero() -> None:
    items = geometries()
    signal = synthetic_colour_vectors(
        items,
        effect_size=8.0,
        scenario="shared_geographic_boundary",
        shared_anchor_latitude=0.0,
        shared_anchor_longitude=0.0,
        rng=np.random.default_rng(9),
    )
    observed, p_value = permutation_p_value(
        signal,
        items,
        min_detectable_species=3,
        permutations=19,
        rng=np.random.default_rng(11),
    )
    assert observed > 0
    assert 0.05 <= p_value <= 1.0

    batched_observed, batched_p = batched_permutation_p_value(
        signal,
        items,
        min_detectable_species=3,
        permutations=19,
        rng=np.random.default_rng(11),
        batch_size=7,
    )
    assert batched_observed == observed
    assert 0.05 <= batched_p <= 1.0
