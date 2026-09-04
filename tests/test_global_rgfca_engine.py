from __future__ import annotations

import numpy as np
import pandas as pd

from fcp_pipeline.global_barrier_field import edge_jensen_shannon_divergence
from fcp_pipeline.global_rgfca_engine import (
    _raw_jsd_from_source_rows,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
    prepare_sparse_outer_geometry,
    run_g1_shard,
)
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule


def _synthetic_pool(n_species: int = 12, photos_per_species: int = 16, *, planted: bool = True) -> pd.DataFrame:
    rows = []
    longitudes = np.linspace(-150.0, 150.0, photos_per_species)
    for s in range(n_species):
        species = f"species_{s:03d}"
        for p, lon in enumerate(longitudes):
            if planted:
                colour = [1.0, 0.0, 0.0, 0.0] if lon < 0 else [0.0, 0.0, 1.0, 0.0]
            else:
                colour = [0.25, 0.25, 0.25, 0.25]
            rows.append(
                {
                    "photo_id": 1_000_000 + s * 1000 + p,
                    "species": species,
                    "latitude": float((s % 3 - 1) * 0.5),
                    "longitude": float(lon),
                    "colour_white": colour[0],
                    "colour_yellow_orange": colour[1],
                    "colour_red_pink": colour[2],
                    "colour_blue_purple": colour[3],
                }
            )
    return pd.DataFrame(rows)


def test_pairwise_cache_matches_direct_jsd_on_fixed_edges():
    pool = canonical_colour_pool(_synthetic_pool(n_species=3, photos_per_species=10))
    cache = build_pairwise_jsd_cache(pool)
    rows = np.arange(10, dtype=np.int64)
    edges = np.column_stack([rows[:-1], rows[1:]])
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    cached = _raw_jsd_from_source_rows(cache, edges, identity)
    direct = edge_jensen_shannon_divergence(pool.colours, edges)
    assert np.allclose(cached, direct)


def test_null_source_rows_stay_inside_species_blocks_and_are_reproducible():
    pool = canonical_colour_pool(_synthetic_pool(n_species=5, photos_per_species=12))
    a = null_source_row_matrix(pool, [0, 7, 11], master_seed=222)
    b = null_source_row_matrix(pool, [0, 7, 11], master_seed=222)
    assert np.array_equal(a, b)
    for row in a:
        assert np.array_equal(pool.species[row], pool.species)
        for start, stop in pool.species_slices:
            assert sorted(row[start:stop].tolist()) == list(range(start, stop))


def test_sparse_outer_geometry_keeps_species_equal_total_kernel_weight():
    frame = _synthetic_pool(n_species=10, photos_per_species=14)
    pool = canonical_colour_pool(frame)
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=4,
        species_per_outer=8,
        photos_per_species=10,
        minimum_pool_photos_per_species=12,
        species_seed=5,
        photo_master_seed=6,
    )
    geometry = prepare_sparse_outer_geometry(
        pool,
        schedule.outer_photo_ids[0],
        schedule.outer_species[0],
        k=3,
        kernel_km=800.0,
        minimum_distinct_species=3,
    )
    assert geometry.weighted_kernel.shape[1] == 36 * 18
    for species_index in range(8):
        total = float(geometry.weighted_kernel[geometry.edge_species_index == species_index].sum())
        assert np.isclose(total, 1.0)
    assert np.isclose(float(geometry.opportunity.sum()), 8.0)
    assert geometry.evaluable.any()


def test_end_to_end_uniform_fixture_is_null_with_all_consensus_fields_zero_concentration():
    result = run_g1_shard(
        _synthetic_pool(n_species=10, photos_per_species=14, planted=False),
        null_indices=[0, 1, 2, 3],
        n_outer=6,
        species_per_outer=8,
        photos_per_species=10,
        minimum_pool_photos_per_species=12,
        k=3,
        n_lon=18,
        n_sinlat=9,
        kernel_km=1000.0,
        minimum_distinct_species=3,
        species_seed=11,
        photo_master_seed=12,
        null_master_seed=13,
    )
    assert np.isclose(result.observed_concentration, 0.0, atol=1e-15)
    assert np.allclose(result.null_concentrations, 0.0, atol=1e-15)
    assert result.null_unit.startswith("complete_photo_colour_vector")


def test_end_to_end_planted_shared_transition_exceeds_small_matched_null_panel():
    result = run_g1_shard(
        _synthetic_pool(n_species=12, photos_per_species=16, planted=True),
        null_indices=list(range(12)),
        n_outer=8,
        species_per_outer=10,
        photos_per_species=12,
        minimum_pool_photos_per_species=14,
        k=3,
        n_lon=18,
        n_sinlat=9,
        kernel_km=900.0,
        minimum_distinct_species=4,
        species_seed=21,
        photo_master_seed=22,
        null_master_seed=23,
    )
    assert np.isfinite(result.observed_concentration)
    assert np.isfinite(result.null_concentrations).all()
    assert result.observed_concentration > float(np.mean(result.null_concentrations))
    assert np.count_nonzero(result.null_concentrations >= result.observed_concentration) <= 2
