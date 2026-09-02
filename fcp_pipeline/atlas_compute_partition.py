"""Deterministic colour-blind compute partitioning for terminal atlas v5."""

from __future__ import annotations

import hashlib
from typing import Any, Iterable, Mapping, Sequence

from .atlas_measurement import measurement_shard

PROTOCOL = "jbi-atlas-compute-partition-amendment-v1"
SEMANTIC_SHARDS = 16
PARTITIONS_PER_SEMANTIC_SHARD = 16
TOTAL_COMPUTE_PARTITIONS = 256
PARTITION_SALT = "fcp-atlas-v5-compute-partition"


class ComputePartitionError(ValueError):
    pass


def validate_compute_partition_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != PROTOCOL:
        raise ComputePartitionError("unexpected compute partition protocol")
    if contract.get("status") != "prospectively_frozen_before_any_terminal_scaleout_candidate_pixel":
        raise ComputePartitionError("compute partition contract was not frozen pre-pixel")
    trigger = contract.get("trigger", {})
    if trigger.get("candidate_image_pixels_opened") is not False or trigger.get("terminal_scaleout_colour_measured") is not False:
        raise ComputePartitionError("compute partition freeze contains opened terminal outcomes")
    parent = contract.get("immutable_parent", {})
    if (
        parent.get("semantic_acquisition_shards") != SEMANTIC_SHARDS
        or parent.get("semantic_measurement_shards") != SEMANTIC_SHARDS
        or parent.get("torch_threads_per_measurement_worker") != 2
        or parent.get("maximum_concurrent_workers") != 8
    ):
        raise ComputePartitionError("compute partition parent execution rules changed")
    partition = contract.get("partition", {})
    if (
        partition.get("partitions_per_semantic_shard") != PARTITIONS_PER_SEMANTIC_SHARD
        or partition.get("total_compute_partitions") != TOTAL_COMPUTE_PARTITIONS
        or partition.get("assignment_uses") != ["measurement_id only"]
        or partition.get("same_partition_for_acquisition_and_measurement") is not True
        or partition.get("all_256_partitions_required") is not True
        or partition.get("early_stopping") is not False
        or partition.get("favourable_partition_selection_forbidden") is not True
    ):
        raise ComputePartitionError("compute partition semantics changed")
    reassembly = contract.get("reassembly", {})
    if (
        reassembly.get("target_semantic_shards") != SEMANTIC_SHARDS
        or reassembly.get("required_compute_partitions_per_semantic_shard") != PARTITIONS_PER_SEMANTIC_SHARD
        or reassembly.get("duplicate_measurement_id_rule") != "fail closed"
        or reassembly.get("missing_measurement_id_rule") != "fail closed"
        or reassembly.get("completeness_gate_unchanged") is not True
    ):
        raise ComputePartitionError("compute partition reassembly rule changed")


def compute_partition(measurement_id: str, partitions_per_shard: int = PARTITIONS_PER_SEMANTIC_SHARD) -> int:
    """Assign one blinded measurement ID to a deterministic subpartition."""
    if partitions_per_shard != PARTITIONS_PER_SEMANTIC_SHARD:
        raise ComputePartitionError("terminal compute partition count is frozen at 16")
    measurement_id = str(measurement_id)
    if not measurement_id:
        raise ComputePartitionError("measurement_id cannot be empty")
    digest = hashlib.sha256(f"{PARTITION_SALT}\x1f{measurement_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False) % partitions_per_shard


def compute_partition_coordinates(measurement_id: str) -> tuple[int, int]:
    return (
        measurement_shard(str(measurement_id), SEMANTIC_SHARDS),
        compute_partition(str(measurement_id)),
    )


def select_compute_partition(
    rows: Sequence[Mapping[str, Any]], *, semantic_shard_index: int, compute_partition_index: int
) -> list[dict[str, Any]]:
    if not 0 <= semantic_shard_index < SEMANTIC_SHARDS:
        raise ComputePartitionError("semantic_shard_index must lie in 0..15")
    if not 0 <= compute_partition_index < PARTITIONS_PER_SEMANTIC_SHARD:
        raise ComputePartitionError("compute_partition_index must lie in 0..15")
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in rows:
        measurement_id = str(raw.get("measurement_id") or "")
        if not measurement_id:
            raise ComputePartitionError("row lacks measurement_id")
        if measurement_id in seen:
            raise ComputePartitionError("input contains duplicate measurement_id")
        seen.add(measurement_id)
        if compute_partition_coordinates(measurement_id) == (semantic_shard_index, compute_partition_index):
            selected.append(dict(raw))
    return selected


def validate_partition_coverage(
    expected_measurement_ids: Iterable[str],
    partition_rows: Mapping[tuple[int, int], Sequence[Mapping[str, Any]]],
) -> dict[str, Any]:
    """Require all 256 partitions and exact union membership before reassembly."""
    expected = {str(value) for value in expected_measurement_ids}
    if not expected or "" in expected:
        raise ComputePartitionError("expected measurement IDs must be non-empty")
    required_keys = {
        (semantic, partition)
        for semantic in range(SEMANTIC_SHARDS)
        for partition in range(PARTITIONS_PER_SEMANTIC_SHARD)
    }
    if set(partition_rows) != required_keys:
        missing = sorted(required_keys - set(partition_rows))
        extra = sorted(set(partition_rows) - required_keys)
        raise ComputePartitionError(f"partition key set changed; missing={missing[:5]} extra={extra[:5]}")

    observed: set[str] = set()
    counts: dict[str, int] = {}
    by_semantic: dict[int, int] = {index: 0 for index in range(SEMANTIC_SHARDS)}
    for key in sorted(required_keys):
        semantic, partition = key
        local_seen: set[str] = set()
        for row in partition_rows[key]:
            measurement_id = str(row.get("measurement_id") or "")
            if not measurement_id:
                raise ComputePartitionError(f"partition {key} contains an empty measurement_id")
            if measurement_id in local_seen:
                raise ComputePartitionError(f"partition {key} repeats measurement_id {measurement_id}")
            local_seen.add(measurement_id)
            if compute_partition_coordinates(measurement_id) != key:
                raise ComputePartitionError(f"measurement_id {measurement_id} is stored in the wrong compute partition")
            counts[measurement_id] = counts.get(measurement_id, 0) + 1
            observed.add(measurement_id)
            by_semantic[semantic] += 1
    duplicates = sorted(mid for mid, count in counts.items() if count != 1)
    if duplicates:
        raise ComputePartitionError(f"measurement IDs appear in multiple partitions: {duplicates[:10]}")
    missing_ids = sorted(expected - observed)
    extra_ids = sorted(observed - expected)
    if missing_ids or extra_ids:
        raise ComputePartitionError(f"compute partition union differs from frozen denominator; missing={missing_ids[:10]} extra={extra_ids[:10]}")

    expected_by_semantic = {
        semantic: sum(measurement_shard(mid, SEMANTIC_SHARDS) == semantic for mid in expected)
        for semantic in range(SEMANTIC_SHARDS)
    }
    if by_semantic != expected_by_semantic:
        raise ComputePartitionError("semantic shard membership changed during compute partitioning")
    return {
        "status": "pass_exact_256_compute_partition_coverage",
        "measurement_ids": len(expected),
        "semantic_shards": SEMANTIC_SHARDS,
        "compute_partitions": TOTAL_COMPUTE_PARTITIONS,
        "rows_by_semantic_shard": by_semantic,
    }


__all__ = [
    "ComputePartitionError",
    "PARTITIONS_PER_SEMANTIC_SHARD",
    "SEMANTIC_SHARDS",
    "TOTAL_COMPUTE_PARTITIONS",
    "compute_partition",
    "compute_partition_coordinates",
    "select_compute_partition",
    "validate_compute_partition_contract",
    "validate_partition_coverage",
]
