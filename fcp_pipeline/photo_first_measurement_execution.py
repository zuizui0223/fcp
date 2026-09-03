"""Pre-pixel firewall and complete-result reassembly for the fresh photo-first atlas.

This module contains no network or image decoding. It validates the one-shot
candidate metadata gate, creates a location-blind worker interface, assigns
stable compute partitions, and refuses the coordinate-colour join until every
frozen candidate measurement ID has exactly one terminal result.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd


EXECUTION_PROTOCOL = "random-photo-first-measurement-execution-v1"
MEASUREMENT_PROTOCOL = "random-photo-first-measurement-v1"
CANDIDATE_PROTOCOL = "random-photo-first-candidate-pool-v1"
WORKER_FIELDS = ("measurement_id", "image_filename", "photo_license")
ACQUISITION_FIELDS = (
    "measurement_id",
    "image_filename",
    "photo_url_large",
    "photo_license",
)
BIOLOGICAL_MORPHS = frozenset(
    {"white", "yellow_orange", "red_pink", "blue_purple"}
)
ALL_MORPHS = BIOLOGICAL_MORPHS | {"mixed_uncertain"}
TERMINAL_STATUSES = frozenset(
    {
        "classified_four_state_morph",
        "not_evaluable_insufficient_flower_pixels",
        "not_evaluable_no_biological_palette_mass",
        "not_evaluable_ambiguous_palette_composition",
        "not_evaluable_roi_or_flip_gate",
        "image_acquisition_failed",
    }
)


@dataclass(frozen=True)
class MeasurementFirewall:
    worker_manifest: pd.DataFrame
    acquisition_key: pd.DataFrame
    metadata_join_key: pd.DataFrame
    manifest: dict[str, Any]


@dataclass(frozen=True)
class MeasurementReassembly:
    joined_photos: pd.DataFrame
    result_manifest: dict[str, Any]


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def canonical_candidate_sha256(candidate: pd.DataFrame) -> str:
    required = {"cell_id", "observation_id", "photo_id"}
    missing = sorted(required.difference(candidate.columns))
    if missing:
        raise ValueError(f"candidate table missing hash columns: {missing}")
    canonical = candidate.sort_values(
        ["cell_id", "observation_id", "photo_id"]
    ).to_csv(index=False, lineterminator="\n")
    return _sha256_text(canonical)


def _blind(salt: str, label: str, value: object, *, length: int = 24) -> str:
    payload = f"{salt}\x1f{label}\x1f{value}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest().upper()[:length]


def measurement_id(photo_id: object, *, salt: str) -> str:
    return "FCPR-" + _blind(salt, "photo", photo_id)


def _stable_mod(prefix: str, measurement: str, modulus: int) -> int:
    if int(modulus) < 1:
        raise ValueError("partition modulus must be positive")
    digest = hashlib.sha256(
        f"{prefix}\x1f{measurement}".encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big") % int(modulus)


def semantic_shard(measurement: str, shard_count: int = 32) -> int:
    return _stable_mod("fcp-random-photo-first-semantic", measurement, shard_count)


def compute_partition(measurement: str, partition_count: int = 4) -> int:
    return _stable_mod("fcp-random-photo-first-compute", measurement, partition_count)


def validate_execution_contract(
    measurement_contract: Mapping[str, Any],
    execution_contract: Mapping[str, Any],
) -> None:
    if measurement_contract.get("protocol") != MEASUREMENT_PROTOCOL:
        raise ValueError("unexpected photo-first measurement contract")
    if measurement_contract.get("status") != "frozen_before_any_fresh_candidate_image_pixel":
        raise ValueError("measurement contract is not frozen pre-pixel")
    if execution_contract.get("protocol") != EXECUTION_PROTOCOL:
        raise ValueError("unexpected measurement execution protocol")
    if execution_contract.get("status") != (
        "prospectively_frozen_before_any_fresh_candidate_image_pixel"
    ):
        raise ValueError("measurement execution is not frozen pre-pixel")
    if execution_contract.get("parent_measurement_protocol") != MEASUREMENT_PROTOCOL:
        raise ValueError("measurement execution parent protocol drifted")
    if any(
        value is not False
        for value in execution_contract.get("outcome_firewall", {}).values()
    ):
        raise ValueError("measurement execution contract contains opened outcomes")

    blinding = execution_contract.get("blinding", {})
    if tuple(blinding.get("worker_allowed_fields", ())) != WORKER_FIELDS:
        raise ValueError("location-blind worker interface changed")
    if tuple(blinding.get("sealed_acquisition_key_fields", ())) != ACQUISITION_FIELDS:
        raise ValueError("sealed acquisition interface changed")
    if not str(blinding.get("salt") or ""):
        raise ValueError("measurement blinding salt is missing")

    partition = execution_contract.get("partitioning", {})
    if (
        int(partition.get("semantic_shards", 0)) != 32
        or int(partition.get("compute_partitions_per_semantic_shard", 0)) != 4
        or int(partition.get("total_compute_partitions", 0)) != 128
        or int(partition.get("maximum_concurrent_compute_partitions", 0)) != 8
        or partition.get("all_partitions_required") is not True
        or partition.get("early_stopping") is not False
    ):
        raise ValueError("measurement partitioning changed")

    measurement = execution_contract.get("measurement", {})
    if (
        measurement.get("all_candidate_rows_reach_exactly_one_terminal_result")
        is not True
        or measurement.get("terminal_result_contains_coordinates") is not False
        or measurement.get("terminal_result_contains_species") is not False
        or measurement.get("persist_flower_mask") is not False
        or measurement.get("persist_image_pixels") is not False
    ):
        raise ValueError("measurement result firewall changed")


def validate_candidate_pixel_gate(
    candidate: pd.DataFrame,
    candidate_manifest: Mapping[str, Any],
    measurement_contract: Mapping[str, Any],
    execution_contract: Mapping[str, Any],
) -> None:
    """Fail closed unless the exact one-shot metadata pool authorizes pixels."""

    validate_execution_contract(measurement_contract, execution_contract)
    if candidate_manifest.get("protocol") != CANDIDATE_PROTOCOL:
        raise ValueError("unexpected candidate-pool protocol")
    if candidate_manifest.get("status") != (
        "metadata_pool_frozen_before_candidate_image_pixels"
    ):
        raise ValueError("candidate pool has not reached the frozen metadata state")
    firewall = candidate_manifest.get("outcome_firewall", {})
    if (
        firewall.get("candidate_image_pixels_opened") is not False
        or firewall.get("morph_used_for_selection") is not False
        or firewall.get("species_list_fixed_or_targeted") is not False
        or firewall.get("legacy_pr21_terminal_records_used") is not False
    ):
        raise ValueError("candidate-pool outcome firewall failed")
    gate = candidate_manifest.get("premeasurement_h1_gate", {})
    if gate.get("pixels_may_open") is not True:
        raise RuntimeError(
            "not_evaluable_candidate_sampling_capacity_before_pixels"
        )
    expected_count = int(candidate_manifest.get("counts", {}).get("observations", -1))
    if expected_count != len(candidate):
        raise ValueError("candidate row count differs from frozen manifest")
    if candidate["observation_id"].nunique() != len(candidate):
        raise ValueError("candidate observation IDs are not unique")
    if candidate["photo_id"].nunique() != len(candidate):
        raise ValueError("candidate photo IDs are not unique")
    observed_sha = canonical_candidate_sha256(candidate)
    if observed_sha != str(candidate_manifest.get("candidate_table_sha256") or ""):
        raise ValueError("candidate table SHA-256 differs from frozen manifest")

    parent_gate = measurement_contract.get("pixel_opening_gate", {})
    if (
        parent_gate.get("required_gate_value") is not True
        or parent_gate.get("pr21_terminal_records_are_input") is not False
    ):
        raise ValueError("measurement parent pixel-opening rule drifted")
    if execution_contract.get("pixel_gate", {}).get(
        "legacy_pr21_terminal_records_permitted"
    ) is not False:
        raise ValueError("execution contract permits forbidden terminal records")


def build_measurement_firewall(
    candidate: pd.DataFrame,
    candidate_manifest: Mapping[str, Any],
    measurement_contract: Mapping[str, Any],
    execution_contract: Mapping[str, Any],
) -> MeasurementFirewall:
    validate_candidate_pixel_gate(
        candidate,
        candidate_manifest,
        measurement_contract,
        execution_contract,
    )
    required = {
        "photo_id",
        "photo_url_large",
        "photo_license",
        "species",
        "latitude",
        "longitude",
    }
    missing = sorted(required.difference(candidate.columns))
    if missing:
        raise ValueError(f"candidate table missing measurement fields: {missing}")
    salt = str(execution_contract["blinding"]["salt"])
    metadata = candidate.reset_index(drop=True).copy()
    metadata.insert(
        0,
        "measurement_id",
        [measurement_id(value, salt=salt) for value in metadata["photo_id"]],
    )
    if metadata["measurement_id"].nunique() != len(metadata):
        raise ValueError("blinded measurement IDs are not unique")
    worker = pd.DataFrame(
        {
            "measurement_id": metadata["measurement_id"].astype(str),
            "image_filename": metadata["measurement_id"].astype(str) + ".jpg",
            "photo_license": metadata["photo_license"].astype(str),
        }
    )
    acquisition = pd.DataFrame(
        {
            "measurement_id": worker["measurement_id"],
            "image_filename": worker["image_filename"],
            "photo_url_large": metadata["photo_url_large"].astype(str),
            "photo_license": worker["photo_license"],
        }
    )
    metadata_join = metadata.drop(columns=["photo_url_large"]).copy()
    manifest = {
        "protocol": EXECUTION_PROTOCOL,
        "status": "measurement_firewall_frozen_before_candidate_pixels",
        "candidate_table_sha256": canonical_candidate_sha256(candidate),
        "candidate_rows": int(len(candidate)),
        "measurement_ids": int(worker["measurement_id"].nunique()),
        "worker_fields": list(worker.columns),
        "acquisition_fields": list(acquisition.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "coordinate_colour_join_opened": False,
        "candidate_pixels_opened": False,
    }
    return MeasurementFirewall(
        worker_manifest=worker,
        acquisition_key=acquisition,
        metadata_join_key=metadata_join,
        manifest=manifest,
    )


def select_measurement_partition(
    worker_manifest: pd.DataFrame,
    acquisition_key: pd.DataFrame,
    *,
    semantic_shard_index: int,
    compute_partition_index: int,
    semantic_shards: int = 32,
    compute_partitions_per_shard: int = 4,
) -> pd.DataFrame:
    if tuple(worker_manifest.columns) != WORKER_FIELDS:
        raise ValueError("worker manifest fields changed or leaked")
    if tuple(acquisition_key.columns) != ACQUISITION_FIELDS:
        raise ValueError("acquisition key fields changed")
    if not (0 <= int(semantic_shard_index) < int(semantic_shards)):
        raise ValueError("semantic shard index is outside the frozen range")
    if not (
        0
        <= int(compute_partition_index)
        < int(compute_partitions_per_shard)
    ):
        raise ValueError("compute partition index is outside the frozen range")
    if worker_manifest["measurement_id"].nunique() != len(worker_manifest):
        raise ValueError("worker manifest measurement IDs are not unique")
    if acquisition_key["measurement_id"].nunique() != len(acquisition_key):
        raise ValueError("acquisition measurement IDs are not unique")
    if set(worker_manifest["measurement_id"]) != set(acquisition_key["measurement_id"]):
        raise ValueError("worker and acquisition keys have different denominators")

    merged = worker_manifest.merge(
        acquisition_key,
        on=["measurement_id", "image_filename", "photo_license"],
        how="inner",
        validate="one_to_one",
    )
    keep = []
    for value in merged["measurement_id"].astype(str):
        keep.append(
            semantic_shard(value, semantic_shards) == int(semantic_shard_index)
            and compute_partition(value, compute_partitions_per_shard)
            == int(compute_partition_index)
        )
    return merged.loc[np.asarray(keep, dtype=bool)].reset_index(drop=True)


def validate_terminal_partition_results(
    results: pd.DataFrame,
    expected_measurement_ids: Sequence[str],
) -> None:
    expected = tuple(str(value) for value in expected_measurement_ids)
    if len(set(expected)) != len(expected):
        raise ValueError("partition expected IDs are not unique")
    if "measurement_id" not in results.columns:
        raise ValueError("terminal results lack measurement_id")
    observed = results["measurement_id"].astype(str)
    if observed.nunique() != len(results):
        raise ValueError("terminal partition results contain duplicate IDs")
    if set(observed) != set(expected) or len(results) != len(expected):
        raise ValueError("terminal partition result coverage is incomplete")
    forbidden = (
        "species",
        "taxon",
        "observation",
        "photo_id",
        "photo_url",
        "latitude",
        "longitude",
        "observer",
        "observed_on",
        "climate",
        "pollinator",
    )
    leaked = [
        str(column)
        for column in results.columns
        if any(token in str(column).casefold() for token in forbidden)
    ]
    if leaked:
        raise ValueError(f"terminal measurement result leaked metadata: {leaked}")
    for row in results.itertuples(index=False):
        morph = str(getattr(row, "morph", ""))
        status = str(getattr(row, "measurement_status", ""))
        if morph not in ALL_MORPHS or status not in TERMINAL_STATUSES:
            raise ValueError("terminal result contains unknown morph/status")
        if status == "classified_four_state_morph" and morph not in BIOLOGICAL_MORPHS:
            raise ValueError("classified state is not one of the four biological morphs")
        if status != "classified_four_state_morph" and morph != "mixed_uncertain":
            raise ValueError("non-classified terminal state must be mixed_uncertain")


def reassemble_complete_measurement(
    partition_results: Sequence[pd.DataFrame],
    worker_manifest: pd.DataFrame,
    metadata_join_key: pd.DataFrame,
    *,
    expected_partition_receipts: int = 128,
) -> MeasurementReassembly:
    """Open the metadata join only after exact complete terminal measurement."""

    if len(partition_results) != int(expected_partition_receipts):
        raise ValueError(
            "not_evaluable_incomplete_measurement_partitions: "
            f"{len(partition_results)} != {int(expected_partition_receipts)}"
        )
    expected_ids = worker_manifest["measurement_id"].astype(str)
    if expected_ids.nunique() != len(worker_manifest):
        raise ValueError("worker denominator contains duplicate measurement IDs")
    combined = pd.concat(partition_results, ignore_index=True)
    if combined["measurement_id"].astype(str).nunique() != len(combined):
        raise ValueError("reassembled terminal results contain duplicate IDs")
    if len(combined) != len(worker_manifest) or set(combined["measurement_id"].astype(str)) != set(expected_ids):
        raise ValueError("not_evaluable_incomplete_terminal_measurement_coverage")
    if metadata_join_key["measurement_id"].astype(str).nunique() != len(metadata_join_key):
        raise ValueError("sealed metadata join key contains duplicate IDs")
    if set(metadata_join_key["measurement_id"].astype(str)) != set(expected_ids):
        raise ValueError("sealed metadata join key does not match worker denominator")

    joined = metadata_join_key.merge(
        combined,
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != len(worker_manifest):
        raise ValueError("coordinate-colour join lost frozen candidate rows")
    counts = joined["measurement_status"].astype(str).value_counts().to_dict()
    morph_counts = joined["morph"].astype(str).value_counts().to_dict()
    manifest = {
        "protocol": EXECUTION_PROTOCOL,
        "status": "complete_fresh_location_blind_measurement_and_join",
        "frozen_candidate_rows": int(len(worker_manifest)),
        "terminal_result_rows": int(len(combined)),
        "joined_rows": int(len(joined)),
        "partition_receipts": int(len(partition_results)),
        "coordinate_colour_join_opened_after_complete_measurement": True,
        "terminal_status_counts": {str(k): int(v) for k, v in counts.items()},
        "morph_counts": {str(k): int(v) for k, v in morph_counts.items()},
        "classified_rows": int(joined["morph"].isin(BIOLOGICAL_MORPHS).sum()),
        "mixed_uncertain_rows": int(joined["morph"].eq("mixed_uncertain").sum()),
        "legacy_pr21_terminal_records_used": False,
        "h1_run": False,
        "h2_run": False,
    }
    return MeasurementReassembly(joined_photos=joined, result_manifest=manifest)


__all__ = [
    "ACQUISITION_FIELDS",
    "ALL_MORPHS",
    "BIOLOGICAL_MORPHS",
    "EXECUTION_PROTOCOL",
    "MeasurementFirewall",
    "MeasurementReassembly",
    "TERMINAL_STATUSES",
    "WORKER_FIELDS",
    "build_measurement_firewall",
    "canonical_candidate_sha256",
    "compute_partition",
    "measurement_id",
    "reassemble_complete_measurement",
    "select_measurement_partition",
    "semantic_shard",
    "validate_candidate_pixel_gate",
    "validate_execution_contract",
    "validate_terminal_partition_results",
]
