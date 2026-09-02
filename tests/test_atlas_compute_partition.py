from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_compute_partition import (
    ComputePartitionError,
    PARTITIONS_PER_SEMANTIC_SHARD,
    SEMANTIC_SHARDS,
    compute_partition,
    compute_partition_coordinates,
    select_compute_partition,
    validate_compute_partition_contract,
    validate_partition_coverage,
)

CONTRACT = Path("docs/supporting/jbi_atlas_compute_partition_amendment_v1.json")


def test_compute_partition_contract_is_pre_pixel_and_fixed() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_compute_partition_contract(contract)
    assert contract["partition"]["total_compute_partitions"] == 256
    assert contract["partition"]["assignment_uses"] == ["measurement_id only"]


def test_assignment_is_deterministic_and_bounded() -> None:
    ids = [f"FCPM-{index:06d}" for index in range(1000)]
    first = [compute_partition(value) for value in ids]
    second = [compute_partition(value) for value in ids]
    assert first == second
    assert all(0 <= value < PARTITIONS_PER_SEMANTIC_SHARD for value in first)
    coords = [compute_partition_coordinates(value) for value in ids]
    assert all(0 <= semantic < SEMANTIC_SHARDS and 0 <= partition < PARTITIONS_PER_SEMANTIC_SHARD for semantic, partition in coords)


def test_select_compute_partition_uses_measurement_id_only() -> None:
    rows = [
        {
            "measurement_id": f"FCPM-{index:04d}",
            "species": f"Species {index % 3}",
            "latitude": index * 0.1,
            "colour": index,
        }
        for index in range(100)
    ]
    target = compute_partition_coordinates(rows[17]["measurement_id"])
    selected = select_compute_partition(
        rows,
        semantic_shard_index=target[0],
        compute_partition_index=target[1],
    )
    expected = [row for row in rows if compute_partition_coordinates(row["measurement_id"]) == target]
    assert selected == expected


def test_exact_256_partition_union_reassembles_frozen_denominator() -> None:
    ids = [f"FCPM-{index:06d}" for index in range(4096)]
    partitions = {
        (semantic, partition): []
        for semantic in range(SEMANTIC_SHARDS)
        for partition in range(PARTITIONS_PER_SEMANTIC_SHARD)
    }
    for measurement_id in ids:
        partitions[compute_partition_coordinates(measurement_id)].append({"measurement_id": measurement_id})
    receipt = validate_partition_coverage(ids, partitions)
    assert receipt["status"] == "pass_exact_256_compute_partition_coverage"
    assert receipt["measurement_ids"] == len(ids)

    broken = {key: list(rows) for key, rows in partitions.items()}
    source_key = next(key for key, rows in broken.items() if rows)
    removed = broken[source_key].pop()
    with pytest.raises(ComputePartitionError, match="union differs"):
        validate_partition_coverage(ids, broken)
    broken[source_key].append(removed)

    wrong_key = next(key for key in broken if key != source_key)
    broken[wrong_key].append(dict(removed))
    with pytest.raises(ComputePartitionError):
        validate_partition_coverage(ids, broken)


def test_missing_partition_key_fails_closed() -> None:
    ids = ["FCPM-a"]
    partitions = {
        (semantic, partition): []
        for semantic in range(SEMANTIC_SHARDS)
        for partition in range(PARTITIONS_PER_SEMANTIC_SHARD)
    }
    partitions[compute_partition_coordinates(ids[0])].append({"measurement_id": ids[0]})
    partitions.pop((0, 0))
    with pytest.raises(ComputePartitionError, match="partition key set changed"):
        validate_partition_coverage(ids, partitions)
