import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_h1_null_maps import persistence_null_maps_cached
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


def synthetic_photos():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    rows = []
    for species_index in range(6):
        for cid, lat, lon in zip(cell_id, latitude, longitude, strict=True):
            for repeat in range(4):
                if repeat == 0 and (int(cid) + species_index) % 4 == 0:
                    morph = "mixed_uncertain"
                elif int(cid) % grid.n_lon < 2:
                    morph = "red_pink"
                else:
                    morph = "blue_purple"
                rows.append(
                    {
                        "species": f"species_{species_index}",
                        "latitude": float(lat),
                        "longitude": float(lon),
                        "morph": morph,
                    }
                )
    return grid, pd.DataFrame(rows)


def test_null_maps_share_observed_edge_order_and_fixed_opportunity_support():
    grid, photos = synthetic_photos()
    result = persistence_null_maps_cached(
        photos,
        grid=grid,
        target_n=72,
        n_replicates=10,
        species_cap_per_cell=2,
        min_photos_per_cell=3,
        n_permutations=11,
        sampling_seed=20260903,
        permutation_seed=20260904,
    )
    assert result.null_persistence.shape == (11, len(result.edge_ids))
    assert result.null_concentrations.shape == (11,)
    assert tuple(result.observed.edge_table.edge_id.astype(str)) == result.edge_ids
    supported = result.observed.edge_table.opportunities.to_numpy(dtype=int) > 0
    assert np.isfinite(result.null_persistence[:, supported]).all()
    assert 0.0 < result.p_upper <= 1.0


def test_null_map_execution_is_deterministic_for_frozen_seeds():
    grid, photos = synthetic_photos()
    kwargs = dict(
        grid=grid,
        target_n=72,
        n_replicates=8,
        species_cap_per_cell=2,
        min_photos_per_cell=3,
        n_permutations=7,
        sampling_seed=101,
        permutation_seed=202,
    )
    first = persistence_null_maps_cached(photos, **kwargs)
    second = persistence_null_maps_cached(photos, **kwargs)
    np.testing.assert_array_equal(first.null_concentrations, second.null_concentrations)
    np.testing.assert_allclose(first.null_persistence, second.null_persistence, equal_nan=True)
    assert first.edge_ids == second.edge_ids
    assert first.p_upper == second.p_upper
