from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from scripts.data.extract_jbi_ch1_direct_florence_colour import (
    read_frozen_shard,
    select_largest_valid_box,
    stable_blind_id,
)


def test_direct_box_selection_is_geometry_only_and_deterministic() -> None:
    parsed = {
        "<OPEN_VOCABULARY_DETECTION>": {
            "bboxes": [
                [10.0, 10.0, 50.0, 50.0],
                [5.0, 5.0, 95.0, 80.0],
                [10.0, 10.0, 50.0, 50.0],
            ],
            "labels": ["anything", "anything", "anything"],
            "scores": [0.99, 0.01, 1.0],
        }
    }
    box, path, raw = select_largest_valid_box(parsed, image_width=100, image_height=100)
    assert box == (5, 5, 95, 80)
    assert "bboxes[1]" in path
    assert len(raw) == 3


def test_invalid_and_out_of_frame_boxes_do_not_create_whole_image_fallback() -> None:
    parsed = {"bboxes": [[0, 0, 0, 0], [200, 200, 201, 201]]}
    try:
        select_largest_valid_box(parsed, image_width=100, image_height=100)
    except ValueError as exc:
        assert "no valid flower box" in str(exc)
    else:
        raise AssertionError("invalid Florence output must fail closed")


def test_blind_id_depends_only_on_frozen_identity() -> None:
    first = stable_blind_id("Species a", "evaluation", "123")
    second = stable_blind_id("Species a", "evaluation", "123")
    assert first == second
    assert first != stable_blind_id("Species a", "calibration", "123")
    assert first != stable_blind_id("Species b", "evaluation", "123")


def test_frozen_sharding_uses_only_split_rank_and_identity(tmp_path: Path) -> None:
    path = tmp_path / "split.csv"
    fields = ["species", "split", "photo_id", "photo_url", "split_rank_hash", "latitude", "longitude"]
    rows = []
    for index in range(80):
        rows.append(
            {
                "species": "Species a",
                "split": "calibration",
                "photo_id": str(1000 + index),
                "photo_url": f"https://example.invalid/{index}.jpg",
                "split_rank_hash": f"{79-index:03d}",
                "latitude": str(index),
                "longitude": str(-index),
            }
        )
    for index in range(120):
        rows.append(
            {
                "species": "Species a",
                "split": "evaluation",
                "photo_id": str(2000 + index),
                "photo_url": f"https://example.invalid/e{index}.jpg",
                "split_rank_hash": f"{119-index:03d}",
                "latitude": str(index),
                "longitude": str(-index),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    calibration = read_frozen_shard(
        path,
        species="Species a",
        split="calibration",
        shard_index=0,
        shard_count=4,
    )
    evaluation = read_frozen_shard(
        path,
        species="Species a",
        split="evaluation",
        shard_index=5,
        shard_count=6,
    )
    assert len(calibration) == 20
    assert len(evaluation) == 20
    assert [row["split_rank_hash"] for row in calibration] == [f"{i:03d}" for i in range(20)]
    assert [row["split_rank_hash"] for row in evaluation] == [f"{i:03d}" for i in range(100, 120)]
    # Coordinates vary strongly but do not enter ordering.
    assert np.std([float(row["latitude"]) for row in calibration]) > 0
