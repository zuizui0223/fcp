from scripts.data.build_inaturalist_photo_development_sample import (
    assign_partition,
    parse_csv,
    round_robin_sample,
    stable_rank,
)


def candidate(obs_id, observer, cell_x, cell_y, week):
    return {
        "id": obs_id,
        "user_id": observer,
        "grid_x": cell_x,
        "grid_y": cell_y,
        "flowering_week": week,
    }


def test_round_robin_deduplicates_observer_cell_week_and_spreads_cells():
    rows = [
        candidate(1, "A", 1, 1, 10),
        candidate(2, "A", 1, 1, 10),
        candidate(3, "B", 1, 1, 11),
        candidate(4, "C", 2, 2, 10),
        candidate(5, "D", 3, 3, 10),
        candidate(6, "E", 1, 1, 12),
    ]
    selected = round_robin_sample(rows, "Species alpha", "seed", target=4)
    keys = {
        (row["user_id"], row["grid_x"], row["grid_y"], row["flowering_week"])
        for row in selected
    }
    assert len(keys) == len(selected) == 4
    assert len({(row["grid_x"], row["grid_y"]) for row in selected[:3]}) == 3


def test_sample_and_partition_are_deterministic_and_exact_40_60():
    rows = [candidate(i, f"observer-{i}", i % 7, i % 5, i % 20 + 1) for i in range(200)]
    first = round_robin_sample(rows, "Species alpha", "seed", target=100)
    second = round_robin_sample(reversed(rows), "Species alpha", "seed", target=100)
    assert [row["id"] for row in first] == [row["id"] for row in second]
    assign_partition(first, "Species alpha", "seed")
    assert sum(row["annotation_partition"] == "development_40" for row in first) == 40
    assert sum(row["annotation_partition"] == "locked_60" for row in first) == 60


def test_stable_rank_changes_with_seed_and_identity():
    assert stable_rank("a", 1) == stable_rank("a", 1)
    assert stable_rank("a", 1) != stable_rank("a", 2)
    assert stable_rank("a", 1) != stable_rank("b", 1)


def test_csv_parser_preserves_quoted_embedded_newlines():
    payload = b'id,description\n1,"line one\nline two"\n'
    assert parse_csv(payload) == [{"id": "1", "description": "line oneline two"}]
