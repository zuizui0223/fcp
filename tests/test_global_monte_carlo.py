import numpy as np
import pytest

from fcp_pipeline.global_monte_carlo import (
    balanced_random_schedule,
    draw_labels,
    opportunity_normalized_field,
    pair_count,
    schedule_audit,
)


def test_primary_scale_schedule_is_fixed_size_unique_and_exactly_balanced():
    species = [f"species_{i:04d}" for i in range(1000)]
    schedule = balanced_random_schedule(
        species,
        n_replicates=200,
        items_per_replicate=250,
        seed=20260904,
    )
    assert schedule.draws.shape == (200, 250)
    assert all(len(set(row.tolist())) == 250 for row in schedule.draws)
    # 200 * 250 / 1000 = exactly 50 inclusions per species.
    assert np.all(schedule.inclusion_counts == 50)
    assert schedule.max_inclusion_imbalance == 0


def test_nondivisible_pool_still_differs_by_at_most_one_inclusion():
    species = [f"species_{i:03d}" for i in range(317)]
    schedule = balanced_random_schedule(
        species,
        n_replicates=200,
        items_per_replicate=250,
        seed=17,
    )
    assert schedule.max_inclusion_imbalance <= 1
    assert int(schedule.inclusion_counts.sum()) == 200 * 250
    assert all(len(set(row.tolist())) == 250 for row in schedule.draws)


def test_same_seed_is_reproducible_and_new_seed_changes_schedule():
    items = [f"p{i}" for i in range(47)]
    first = balanced_random_schedule(items, n_replicates=25, items_per_replicate=20, seed=5)
    second = balanced_random_schedule(items, n_replicates=25, items_per_replicate=20, seed=5)
    third = balanced_random_schedule(items, n_replicates=25, items_per_replicate=20, seed=6)
    assert np.array_equal(first.draws, second.draws)
    assert not np.array_equal(first.draws, third.draws)
    assert len(draw_labels(first, 0)) == 20


def test_schedule_refuses_duplicate_items_or_oversized_draw():
    with pytest.raises(ValueError):
        balanced_random_schedule(["a", "a"], n_replicates=2, items_per_replicate=1, seed=1)
    with pytest.raises(ValueError):
        balanced_random_schedule(["a", "b"], n_replicates=2, items_per_replicate=3, seed=1)


def test_pair_count_matches_twenty_photo_design():
    assert pair_count(20) == 190
    assert pair_count(0) == 0
    with pytest.raises(ValueError):
        pair_count(-1)


def test_opportunity_normalization_preserves_not_evaluable_cells_as_nan():
    field = opportunity_normalized_field(
        [4.0, 2.0, 0.0, 5.0],
        [2.0, 0.0, 1.0, 5.0],
        minimum_opportunity=0.0,
    )
    assert field[0] == pytest.approx(2.0)
    assert np.isnan(field[1])
    assert field[2] == pytest.approx(0.0)
    assert field[3] == pytest.approx(1.0)


def test_schedule_audit_is_json_serializable_shape():
    schedule = balanced_random_schedule(
        [f"s{i}" for i in range(300)],
        n_replicates=200,
        items_per_replicate=250,
        seed=3,
    )
    audit = schedule_audit(schedule)
    assert audit["eligible_items"] == 300
    assert audit["n_replicates"] == 200
    assert audit["items_per_replicate"] == 250
    assert audit["total_item_inclusions"] == 50000
    assert audit["max_inclusion_imbalance"] <= 1
