import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_atlas_v2 import run_boundary_persistence
from fcp_pipeline.photo_first_h1_fast import (
    evaluate_label_codes,
    permute_label_codes_within_species,
    persistence_null_test_cached,
    prepare_h1_plan,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


def synthetic_table(*, include_uncertain=True):
    grid = EqualAreaGrid(n_lon=6, n_sinlat=3)
    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    rows = []
    for species_index in range(8):
        species = f"species_{species_index:02d}"
        for cid, lat, lon in zip(cell_id, latitude, longitude, strict=True):
            col = int(cid) % grid.n_lon
            for repeat in range(4):
                if include_uncertain and repeat == 0 and (int(cid) + species_index) % 5 == 0:
                    morph = "mixed_uncertain"
                elif col < 2:
                    morph = "red_pink"
                elif col < 4:
                    morph = "white"
                else:
                    morph = "blue_purple"
                rows.append(
                    {
                        "species": species,
                        "latitude": float(lat),
                        "longitude": float(lon),
                        "morph": morph,
                    }
                )
    return grid, pd.DataFrame(rows)


def test_cached_observed_matches_frozen_v2_transition_counts_and_concentration():
    grid, photos = synthetic_table()
    kwargs = dict(
        grid=grid,
        target_n=180,
        n_replicates=16,
        species_cap_per_cell=2,
        min_photos_per_cell=4,
        transition_quantile=0.90,
        random_seed=20260903,
    )
    slow = run_boundary_persistence(photos, **kwargs)
    plan = prepare_h1_plan(photos, **kwargs)
    fast = evaluate_label_codes(plan, plan.observed_label_codes)

    pd.testing.assert_series_equal(
        slow.edge_table["edge_id"], fast.edge_table["edge_id"], check_names=False
    )
    np.testing.assert_array_equal(
        slow.edge_table["opportunities"].to_numpy(),
        fast.edge_table["opportunities"].to_numpy(),
    )
    np.testing.assert_array_equal(
        slow.edge_table["transition_count"].to_numpy(),
        fast.edge_table["transition_count"].to_numpy(),
    )
    np.testing.assert_allclose(
        slow.edge_table["persistence"].to_numpy(dtype=float),
        fast.edge_table["persistence"].to_numpy(dtype=float),
        equal_nan=True,
        rtol=0,
        atol=0,
    )
    assert np.isclose(slow.concentration, fast.concentration, rtol=0, atol=1e-15)
    assert np.isclose(slow.transition_rate, fast.transition_rate, rtol=0, atol=1e-15)


def test_vectorized_species_permutation_preserves_species_marginals_and_missing_mask():
    grid, photos = synthetic_table()
    plan = prepare_h1_plan(
        photos,
        grid=grid,
        target_n=120,
        n_replicates=4,
        species_cap_per_cell=2,
        min_photos_per_cell=2,
        random_seed=7,
    )
    permuted = permute_label_codes_within_species(
        plan, rng=np.random.default_rng(12345)
    )

    np.testing.assert_array_equal(
        plan.observed_label_codes < 0,
        permuted < 0,
    )
    for species_code in np.unique(plan.species_codes):
        idx = np.flatnonzero(
            (plan.species_codes == species_code) & (plan.observed_label_codes >= 0)
        )
        before = np.bincount(
            plan.observed_label_codes[idx], minlength=len(plan.morph_levels)
        )
        after = np.bincount(permuted[idx], minlength=len(plan.morph_levels))
        np.testing.assert_array_equal(before, after)


def test_prepared_plan_reuses_fixed_sampling_and_opportunity_denominators():
    grid, photos = synthetic_table()
    plan = prepare_h1_plan(
        photos,
        grid=grid,
        target_n=160,
        n_replicates=10,
        species_cap_per_cell=2,
        min_photos_per_cell=3,
        random_seed=29,
    )
    assert plan.sampled_indices.shape == (10, 160)
    assert plan.evaluable_edges.shape[0] == 10
    np.testing.assert_array_equal(
        plan.opportunities,
        plan.evaluable_edges.sum(axis=0),
    )
    assert np.isfinite(plan.tie_breaks[plan.evaluable_edges]).all()
    assert np.isnan(plan.tie_breaks[~plan.evaluable_edges]).all()


def test_cached_null_is_deterministic_and_returns_valid_monte_carlo_p():
    grid, photos = synthetic_table()
    kwargs = dict(
        grid=grid,
        target_n=180,
        n_replicates=12,
        species_cap_per_cell=2,
        min_photos_per_cell=4,
        transition_quantile=0.90,
        n_permutations=19,
        sampling_seed=20260903,
        permutation_seed=20260904,
    )
    first_observed, first_null, first_p = persistence_null_test_cached(photos, **kwargs)
    second_observed, second_null, second_p = persistence_null_test_cached(photos, **kwargs)
    np.testing.assert_array_equal(first_null, second_null)
    assert first_observed.concentration == second_observed.concentration
    assert first_p == second_p
    assert 0.0 < first_p <= 1.0
    assert len(first_null) == 19
