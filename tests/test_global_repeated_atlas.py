from __future__ import annotations

import numpy as np
import pytest

from fcp_pipeline.global_repeated_atlas import (
    build_repeated_atlas_schedule,
    consensus_field,
    null_source_photo_ids,
    odd_even_consensus,
    running_consensus,
    schedule_audit,
)


def _photo_pool(n_species: int = 12, photos_per_species: int = 12):
    photo_ids = []
    species = []
    for s in range(n_species):
        label = f"species_{s:03d}"
        for p in range(photos_per_species):
            photo_ids.append(1_000_000 + s * 1000 + p)
            species.append(label)
    return np.asarray(photo_ids, dtype=np.int64), np.asarray(species, dtype=object)


def test_repeated_schedule_is_input_order_invariant_and_balanced():
    photo_ids, species = _photo_pool(n_species=12, photos_per_species=12)
    kwargs = dict(
        n_outer=20,
        species_per_outer=7,
        photos_per_species=5,
        minimum_pool_photos_per_species=10,
        species_seed=123,
        photo_master_seed=456,
    )
    a = build_repeated_atlas_schedule(photo_ids, species, **kwargs)
    rng = np.random.default_rng(99)
    order = rng.permutation(len(photo_ids))
    b = build_repeated_atlas_schedule(photo_ids[order], species[order], **kwargs)

    assert a.species_labels == b.species_labels
    assert np.array_equal(a.outer_species, b.outer_species)
    assert np.array_equal(a.outer_photo_ids, b.outer_photo_ids)
    assert np.array_equal(a.species_inclusion_counts, b.species_inclusion_counts)

    audit = schedule_audit(a)
    assert audit["species_max_inclusion_imbalance"] <= 1
    assert audit["photo_max_inclusion_imbalance_within_species"] <= 1
    assert audit["total_species_inclusions"] == 20 * 7

    sorted_draws = np.sort(a.outer_photo_ids, axis=2)
    assert not np.any(np.diff(sorted_draws, axis=2) == 0)


def test_schedule_fails_closed_if_species_pool_below_frozen_minimum():
    photo_ids, species = _photo_pool(n_species=5, photos_per_species=9)
    with pytest.raises(ValueError, match="fail minimum_pool"):
        build_repeated_atlas_schedule(
            photo_ids,
            species,
            n_outer=5,
            species_per_outer=4,
            photos_per_species=5,
            minimum_pool_photos_per_species=10,
        )


def test_null_mapping_is_row_order_invariant_and_species_conditioned():
    photo_ids, species = _photo_pool(n_species=6, photos_per_species=8)
    mapping_a = null_source_photo_ids(
        photo_ids,
        species,
        permutation_index=17,
        master_seed=789,
    )
    rng = np.random.default_rng(5)
    order = rng.permutation(len(photo_ids))
    mapping_b_shuffled = null_source_photo_ids(
        photo_ids[order],
        species[order],
        permutation_index=17,
        master_seed=789,
    )
    recovered = np.empty_like(mapping_b_shuffled)
    recovered[order] = mapping_b_shuffled
    assert np.array_equal(mapping_a, recovered)

    species_by_id = {int(pid): str(sp) for pid, sp in zip(photo_ids, species)}
    for target, source in zip(photo_ids, mapping_a):
        assert species_by_id[int(target)] == species_by_id[int(source)]

    for label in np.unique(species):
        idx = np.flatnonzero(species == label)
        assert sorted(mapping_a[idx].tolist()) == sorted(photo_ids[idx].tolist())


def test_null_mapping_changes_by_permutation_index_without_changing_pool():
    photo_ids, species = _photo_pool(n_species=4, photos_per_species=10)
    a = null_source_photo_ids(photo_ids, species, permutation_index=1, master_seed=101)
    b = null_source_photo_ids(photo_ids, species, permutation_index=2, master_seed=101)
    assert not np.array_equal(a, b)
    assert sorted(a.tolist()) == sorted(photo_ids.tolist())
    assert sorted(b.tolist()) == sorted(photo_ids.tolist())


def test_consensus_equals_explicit_pooled_numerator_over_opportunity():
    fields = np.array(
        [
            [1.0, 10.0, np.nan],
            [3.0, 20.0, 5.0],
            [5.0, np.nan, 9.0],
        ]
    )
    opportunities = np.array(
        [
            [1.0, 1.0, 0.0],
            [3.0, 2.0, 4.0],
            [6.0, 0.0, 1.0],
        ]
    )
    result = consensus_field(fields, opportunities)
    expected = np.array(
        [
            (1.0 * 1.0 + 3.0 * 3.0 + 5.0 * 6.0) / 10.0,
            (10.0 * 1.0 + 20.0 * 2.0) / 3.0,
            (5.0 * 4.0 + 9.0 * 1.0) / 5.0,
        ]
    )
    assert np.allclose(result.field, expected)
    assert np.allclose(result.aggregate_opportunity, [10.0, 3.0, 5.0])
    assert np.array_equal(result.evaluable_outer_counts, [3, 2, 2])
    assert result.concentration >= 0.0


def test_unsupported_cells_remain_nan_not_biological_zero():
    fields = np.array([[1.0, np.nan], [3.0, np.nan]])
    opportunities = np.array([[1.0, 0.0], [1.0, 0.0]])
    result = consensus_field(fields, opportunities)
    assert result.field[0] == pytest.approx(2.0)
    assert np.isnan(result.field[1])
    assert result.aggregate_opportunity[1] == 0.0
    assert result.evaluable_outer_counts[1] == 0


def test_running_and_odd_even_consensus_use_same_aggregation_rule():
    fields = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]])
    opportunities = np.ones_like(fields)
    running = running_consensus(fields, opportunities, checkpoints=(2, 4))
    assert np.allclose(running[2].field, [2.0, 3.0])
    assert np.allclose(running[4].field, [4.0, 5.0])

    odd, even = odd_even_consensus(fields, opportunities)
    assert np.allclose(odd.field, [3.0, 4.0])  # rows 1 and 3 in human numbering
    assert np.allclose(even.field, [5.0, 6.0]) # rows 2 and 4 in human numbering
