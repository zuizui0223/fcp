import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.photo_first_atlas_v2 import (
    BIOLOGICAL_MORPH_LEVELS,
    adjacent_grid_edges,
    cell_first_species_capped_sample,
    coarse_morph_from_palette,
    jensen_shannon_divergence,
    prepare_photo_grid,
    replicate_edge_table,
    run_boundary_persistence,
    species_conditioned_morph_permutation,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid


def test_coarse_palette_mapping_is_species_independent_and_keeps_ambiguity():
    red = {
        "blue": 0.00,
        "bronze": 0.02,
        "magenta": 0.25,
        "orange": 0.01,
        "pink": 0.20,
        "purple": 0.02,
        "red": 0.45,
        "white": 0.04,
        "yellow": 0.01,
    }
    assert coarse_morph_from_palette(red) == "red_pink"

    ambiguous = dict(red)
    ambiguous.update(
        {"magenta": 0.20, "pink": 0.10, "red": 0.10, "white": 0.35, "yellow": 0.15}
    )
    assert coarse_morph_from_palette(ambiguous) == "mixed_uncertain"


def test_jensen_shannon_divergence_has_expected_endpoints():
    assert jensen_shannon_divergence([1, 0, 0], [1, 0, 0]) == 0.0
    assert np.isclose(jensen_shannon_divergence([1, 0], [0, 1]), 1.0)


def test_grid_edges_wrap_longitude_without_duplicate_edges():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    edges = adjacent_grid_edges(grid)
    edge_set = {tuple(row) for row in edges.tolist()}
    assert (0, 3) in edge_set
    assert len(edge_set) == len(edges)
    assert len(edges) == 12


def test_cell_first_sampler_enforces_species_cap_and_target():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    rows = []
    for species in ("a", "b", "c"):
        for repeat in range(5):
            rows.append(
                {
                    "species": species,
                    "latitude": -20.0,
                    "longitude": -150.0,
                    "morph": "red_pink" if repeat % 2 else "white",
                }
            )
    photos = prepare_photo_grid(pd.DataFrame(rows), grid=grid)
    sampled = cell_first_species_capped_sample(
        photos,
        target_n=6,
        species_cap_per_cell=2,
        rng=np.random.default_rng(7),
    )
    assert len(sampled) == 6
    counts = sampled.groupby(["cell_id", "species"]).size()
    assert int(counts.max()) <= 2


def test_species_conditioned_permutation_preserves_each_species_biological_marginal():
    photos = pd.DataFrame(
        {
            "species": ["a", "a", "a", "b", "b", "b"],
            "latitude": [0, 1, 2, 3, 4, 5],
            "longitude": [0, 1, 2, 3, 4, 5],
            "morph": [
                "red_pink",
                "red_pink",
                "white",
                "blue_purple",
                "yellow_orange",
                "yellow_orange",
            ],
        }
    )
    permuted = species_conditioned_morph_permutation(
        photos, rng=np.random.default_rng(123)
    )
    for species in ("a", "b"):
        before = sorted(photos.loc[photos.species == species, "morph"].tolist())
        after = sorted(permuted.loc[permuted.species == species, "morph"].tolist())
        assert before == after


def test_mixed_uncertain_positions_are_fixed_under_species_conditioned_null():
    photos = pd.DataFrame(
        {
            "species": ["a"] * 6 + ["b"] * 6,
            "latitude": np.arange(12),
            "longitude": np.arange(12),
            "morph": [
                "mixed_uncertain",
                "red_pink",
                "white",
                "red_pink",
                "white",
                "mixed_uncertain",
                "blue_purple",
                "mixed_uncertain",
                "yellow_orange",
                "blue_purple",
                "yellow_orange",
                "mixed_uncertain",
            ],
        }
    )
    uncertain_before = photos.morph.eq("mixed_uncertain").to_numpy()
    permuted = species_conditioned_morph_permutation(
        photos, rng=np.random.default_rng(17)
    )
    uncertain_after = permuted.morph.eq("mixed_uncertain").to_numpy()
    np.testing.assert_array_equal(uncertain_before, uncertain_after)


def test_mixed_uncertain_is_not_a_fifth_biological_state_or_cell_support():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    rows = []
    # Cell 0 has six raw photos but only one classifiable photo. It must fail a
    # minimum classifiable-photo threshold of two.
    for morph in [
        "red_pink",
        "mixed_uncertain",
        "mixed_uncertain",
        "mixed_uncertain",
        "mixed_uncertain",
        "mixed_uncertain",
    ]:
        rows.append({"species": "a", "latitude": -20.0, "longitude": -150.0, "morph": morph})
    for morph in ["blue_purple", "blue_purple", "white", "white"]:
        rows.append({"species": "b", "latitude": -20.0, "longitude": -60.0, "morph": morph})
    sampled = prepare_photo_grid(pd.DataFrame(rows), grid=grid)
    table = replicate_edge_table(
        sampled,
        grid=grid,
        morph_levels=BIOLOGICAL_MORPH_LEVELS,
        min_photos_per_cell=2,
        transition_quantile=0.90,
        rng=np.random.default_rng(2),
    )
    row = table[((table.cell_i == 0) & (table.cell_j == 1)) | ((table.cell_i == 1) & (table.cell_j == 0))]
    assert len(row) == 1
    assert int(row.iloc[0].raw_n_i) + int(row.iloc[0].raw_n_j) == 10
    assert bool(row.iloc[0].evaluable) is False


def test_persistence_uses_evaluable_replicates_as_edge_denominator():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    rows = []
    coordinates = [(-20.0, -150.0), (-20.0, -60.0), (20.0, -150.0)]
    for species in ("a", "b"):
        for cell_index, (lat, lon) in enumerate(coordinates):
            for repeat in range(6):
                rows.append(
                    {
                        "species": species,
                        "latitude": lat,
                        "longitude": lon,
                        "morph": "red_pink" if cell_index == 0 else "white",
                    }
                )
    result = run_boundary_persistence(
        pd.DataFrame(rows),
        grid=grid,
        target_n=8,
        n_replicates=20,
        species_cap_per_cell=2,
        min_photos_per_cell=2,
        transition_quantile=0.90,
        random_seed=19,
    )
    supported = result.edge_table[result.edge_table.opportunities > 0]
    assert len(supported) > 0
    assert (supported.opportunities <= 20).all()
    expected = supported.transition_count / supported.opportunities
    np.testing.assert_allclose(supported.persistence, expected)


def test_replicate_sampling_is_deterministic_for_fixed_seed():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    rows = []
    for species in ("a", "b", "c"):
        for lon in (-150.0, -60.0, 30.0, 120.0):
            for repeat in range(3):
                rows.append(
                    {
                        "species": species,
                        "latitude": -20.0 if repeat % 2 == 0 else 20.0,
                        "longitude": lon,
                        "morph": "red_pink" if lon < 0 else "blue_purple",
                    }
                )
    kwargs = dict(
        grid=grid,
        target_n=20,
        n_replicates=10,
        species_cap_per_cell=2,
        min_photos_per_cell=2,
        transition_quantile=0.90,
        random_seed=222,
    )
    first = run_boundary_persistence(pd.DataFrame(rows), **kwargs)
    second = run_boundary_persistence(pd.DataFrame(rows), **kwargs)
    pd.testing.assert_frame_equal(first.edge_table, second.edge_table)
    assert first.concentration == second.concentration


def test_fixed_replicate_size_fails_closed_when_species_capped_capacity_is_too_small():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    photos = pd.DataFrame(
        {
            "species": ["a"] * 8,
            "latitude": [-20.0] * 8,
            "longitude": [-150.0] * 8,
            "morph": ["red_pink", "white"] * 4,
        }
    )
    with pytest.raises(ValueError, match="not_evaluable_fixed_replicate_size"):
        run_boundary_persistence(
            photos,
            grid=grid,
            target_n=3,
            n_replicates=2,
            species_cap_per_cell=2,
            min_photos_per_cell=1,
        )
