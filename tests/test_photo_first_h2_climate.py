import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_atlas import adjacent_grid_edges
from fcp_pipeline.photo_first_h2_climate import (
    aggregate_climate_to_h1_grid,
    build_edge_climate_contrasts,
    climate_concordance_test,
    holm_adjust,
    weighted_pearson,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


def climate_source_for_grid(grid):
    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    rows = []
    for cid, lat, lon in zip(cell_id, latitude, longitude, strict=True):
        col = int(cid) % grid.n_lon
        row = int(cid) // grid.n_lon
        rows.append(
            {
                "latitude": float(lat),
                "longitude": float(lon),
                "bio1": 5.0 + 3.0 * col + 0.5 * row,
                "bio4": 100.0 + 10.0 * row + col,
                "bio12": 500.0 + 100.0 * col + 20.0 * row,
                "bio15": 30.0 + 2.0 * row + 0.25 * col,
                "realm": "realm_west" if col < grid.n_lon // 2 else "realm_east",
                "biome": "biome_a" if row == 0 else "biome_b",
            }
        )
    return pd.DataFrame(rows)


def test_climate_aggregation_is_independent_equal_area_cell_summary():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    source = climate_source_for_grid(grid)
    # Add a second equal-area source centroid to the same H1 cell; the H1 value
    # must be the arithmetic mean and not use photo density or colour.
    extra = source.iloc[[0]].copy()
    extra["bio1"] = extra["bio1"] + 4.0
    extra["realm"] = "realm_west"
    extra["biome"] = "biome_a"
    combined = pd.concat([source, extra], ignore_index=True)
    climate = aggregate_climate_to_h1_grid(combined, grid=grid)
    first = climate.set_index("cell_id").loc[0]
    assert int(first.source_cell_n) == 2
    assert np.isclose(float(first.bio1), float(source.iloc[0].bio1) + 2.0)
    assert first.realm == "realm_west"
    assert first.biome == "biome_a"
    assert climate[["z_bio1", "z_bio4", "z_bio12", "z_bio15"]].notna().all().all()


def test_edge_climate_distance_is_zero_for_equal_cells_and_positive_for_contrast():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    climate = aggregate_climate_to_h1_grid(climate_source_for_grid(grid), grid=grid)
    # Force adjacent cells 0 and 1 to have identical standardized climate.
    climate = climate.copy()
    for variable in ("z_bio1", "z_bio4", "z_bio12", "z_bio15"):
        climate.loc[climate.cell_id == 1, variable] = climate.loc[
            climate.cell_id == 0, variable
        ].iloc[0]
    edges = build_edge_climate_contrasts(climate, grid=grid).set_index("edge_id")
    assert np.isclose(edges.loc["0:1", "multivariate_climate_distance"], 0.0)
    assert edges["multivariate_climate_distance"].dropna().max() > 0.0


def test_weighted_pearson_has_expected_sign_and_unit_endpoint():
    assert np.isclose(weighted_pearson([0, 1, 2], [0, 2, 4], [1, 2, 3]), 1.0)
    assert np.isclose(weighted_pearson([0, 1, 2], [4, 2, 0], [1, 2, 3]), -1.0)


def test_h2_uses_fixed_h1_edge_null_and_detects_planted_climate_alignment():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    climate = aggregate_climate_to_h1_grid(climate_source_for_grid(grid), grid=grid)
    climate_edges = build_edge_climate_contrasts(climate, grid=grid)
    edge_ids = [f"{int(a)}:{int(b)}" for a, b in adjacent_grid_edges(grid)]
    distances = climate_edges.set_index("edge_id").loc[
        edge_ids, "multivariate_climate_distance"
    ].to_numpy(dtype=float)
    observed_persistence = (distances - distances.min()) / (
        distances.max() - distances.min()
    )
    h1 = pd.DataFrame(
        {
            "edge_id": edge_ids,
            "opportunities": np.full(len(edge_ids), 100),
            "transition_count": np.rint(observed_persistence * 100).astype(int),
            "persistence": observed_persistence,
        }
    )
    rng = np.random.default_rng(20260903)
    null = rng.random((39, len(edge_ids)))
    result = climate_concordance_test(
        h1,
        null,
        edge_ids,
        climate_edges,
        minimum_supported_edges=8,
    )
    assert result.supported_edges == len(edge_ids)
    assert result.statistic > 0.99
    assert result.p_upper <= 0.05


def test_h2_within_biome_subset_is_predeclared_and_uses_same_null():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    climate = aggregate_climate_to_h1_grid(climate_source_for_grid(grid), grid=grid)
    climate_edges = build_edge_climate_contrasts(climate, grid=grid)
    edge_ids = climate_edges.edge_id.astype(str).tolist()
    h1 = pd.DataFrame(
        {
            "edge_id": edge_ids,
            "opportunities": np.full(len(edge_ids), 50),
            "persistence": np.linspace(0.0, 1.0, len(edge_ids)),
        }
    )
    null = np.tile(np.linspace(1.0, 0.0, len(edge_ids)), (9, 1))
    result = climate_concordance_test(
        h1,
        null,
        edge_ids,
        climate_edges,
        subset="within_biome",
        minimum_supported_edges=2,
    )
    assert result.subset == "within_biome"
    assert result.supported_edges >= 2


def test_holm_adjustment_is_monotone_in_sorted_p_values():
    raw = np.array([0.01, 0.04, 0.03, 0.20])
    adjusted = holm_adjust(raw)
    np.testing.assert_allclose(adjusted, [0.04, 0.09, 0.09, 0.20])
