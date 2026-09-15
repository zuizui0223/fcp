"""Fail-closed chronology guard for the prospective P500 H2 confirmation.

This module does not acquire images, measure colour, or calculate H2.  Its only
job is to make the execution chronology explicit and machine-checkable from a
new prospective freeze forward.  It deliberately makes no retrospective claim
about whether pixels were accessed before this epoch was created.
"""

from __future__ import annotations

import copy
import hashlib
import json
from typing import Any, Mapping


SCHEMA_VERSION = "p500_execution_epoch_v1"

STATES = (
    "PREOPEN_FROZEN",
    "DISPATCH_RECORDED",
    "MEASUREMENT_COMPLETE",
    "INFERENCE_OPENED",
)

ALLOWED_TRANSITIONS = {
    "PREOPEN_FROZEN": "DISPATCH_RECORDED",
    "DISPATCH_RECORDED": "MEASUREMENT_COMPLETE",
    "MEASUREMENT_COMPLETE": "INFERENCE_OPENED",
}

EXPECTED_FREEZE = {
    "schema_version": SCHEMA_VERSION,
    "epoch_id": "P500_H2_PROSPECTIVE_20260915_E1",
    "state": "PREOPEN_FROZEN",
    "base_branch": "analysis/polymorphism-42111-h1-h2-gates-20260912",
    "base_sha": "26f7a792f93fb53bc6e8c9a64efb45ccb548b9ea",
    "protocol_path": "docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260915.md",
    "protocol_blob_sha": "d1f81202a7a8eedf8625b9c9db58845c0c96678a",
    "authorization_path": "docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_AUTHORIZATION_20260915.json",
    "authorization_blob_sha": "3310b4d573a2dd321f9684250892cee3b1982954",
    "P500_sha256": "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4",
    "metadata_sha256": "a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595",
    "full100_species_sha256": "568a64e692dc21f4a5f76c5eadb615a63b07d71f01bc9cd91e8496ccd0fa3f9e",
    "measurement_machine_source_commit": "9fae6ccdf684a46026f72ba12e98de2c5c54bf2a",
    "frozen_species": 499,
    "frozen_rows": 49900,
    "target_rows_per_species": 100,
    "terminal_partitions": 256,
    "minimum_classifiable_photos_per_species": 40,
    "minimum_measurement_evaluable_species": 250,
    "minimum_primary_H2_vector_species": 20,
    "H2_target": "fixed_q_white_W_structured_null",
    "axis_refit_allowed": False,
    "primary_threshold": 0.10,
    "strict_sensitivity_threshold": 0.20,
    "structured_null_replicates": 999,
    "primary_seed": 20260915,
    "strict_seed": 20261015,
    "species_replacement_allowed": False,
    "row_replacement_allowed": False,
    "target_relaxation_allowed": False,
    "rerun_after_outcome_allowed": False,
    "historical_nonaccess_claimed": False,
    "provenance_scope": "prospective_from_this_receipt_forward_no_retrospective_claim",
}

PREOPEN_FALSE_FLAGS = (
    "measurement_dispatched",
    "image_requests_started",
    "image_bytes_opened",
    "flower_colour_opened",
    "palette_opened",
    "D_opened",
    "H2_W_opened",
    "metadata_colour_join_opened",
)


def canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    """Return a deterministic SHA-256 over a JSON object."""
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _require_bool(value: Any, label: str) -> None:
    if type(value) is not bool:  # bool is intentionally stricter than truthiness.
        raise RuntimeError(f"{label} must be boolean")


def validate_preopen_freeze(receipt: Mapping[str, Any]) -> None:
    """Validate the immutable pre-opening epoch receipt.

    Passing this function means only that the *new epoch receipt* is internally
    consistent.  It does not prove any historical non-access before the receipt.
    """
    for key, expected in EXPECTED_FREEZE.items():
        if receipt.get(key) != expected:
            raise RuntimeError(f"frozen field mismatch: {key}")

    if receipt.get("history") != []:
        raise RuntimeError("preopen history must be empty")

    flags = receipt.get("opening_flags")
    if not isinstance(flags, dict):
        raise RuntimeError("opening_flags missing")
    for flag in PREOPEN_FALSE_FLAGS:
        if flag not in flags:
            raise RuntimeError(f"opening flag missing: {flag}")
        _require_bool(flags[flag], flag)
        if flags[flag]:
            raise RuntimeError(f"preopen flag already opened: {flag}")

    if set(flags) != set(PREOPEN_FALSE_FLAGS):
        raise RuntimeError("unexpected opening flag in preopen receipt")


def _validate_common_event(event: Mapping[str, Any], *, kind: str) -> None:
    if event.get("kind") != kind:
        raise RuntimeError(f"event kind must be {kind}")
    event_id = event.get("event_id")
    if not isinstance(event_id, str) or not event_id.strip():
        raise RuntimeError("event_id must be a non-empty string")
    timestamp = event.get("timestamp_utc")
    if not isinstance(timestamp, str) or not timestamp.endswith("Z"):
        raise RuntimeError("timestamp_utc must be an explicit UTC Z timestamp")


def advance_epoch(
    receipt: Mapping[str, Any], next_state: str, event: Mapping[str, Any]
) -> dict[str, Any]:
    """Advance exactly one chronology state after validating its evidence.

    The caller is responsible for durably saving the returned receipt before
    performing any action that belongs to the new state.
    """
    current = receipt.get("state")
    if current not in STATES:
        raise RuntimeError("unknown current state")
    expected_next = ALLOWED_TRANSITIONS.get(current)
    if expected_next is None:
        raise RuntimeError("terminal epoch cannot advance")
    if next_state != expected_next:
        raise RuntimeError(f"illegal state transition: {current} -> {next_state}")

    prior_history = receipt.get("history")
    if not isinstance(prior_history, list):
        raise RuntimeError("history must be a list")
    seen_ids = {item.get("event_id") for item in prior_history if isinstance(item, dict)}
    if event.get("event_id") in seen_ids:
        raise RuntimeError("duplicate event_id")

    prior_digest = canonical_json_sha256(receipt)
    if event.get("prior_receipt_sha256") != prior_digest:
        raise RuntimeError("prior receipt digest mismatch")

    if current == "PREOPEN_FROZEN":
        validate_preopen_freeze(receipt)
        _validate_common_event(event, kind="first_dispatch_record")
        for field in (
            "image_requests_started_before_record",
            "image_bytes_opened_before_record",
            "measurement_started_before_record",
            "H2_opened_before_record",
        ):
            _require_bool(event.get(field), field)
            if event[field]:
                raise RuntimeError(f"dispatch chronology violation: {field}")
        if event.get("frozen_rows") != EXPECTED_FREEZE["frozen_rows"]:
            raise RuntimeError("dispatch frozen row count mismatch")
        if event.get("frozen_species") != EXPECTED_FREEZE["frozen_species"]:
            raise RuntimeError("dispatch frozen species count mismatch")

    elif current == "DISPATCH_RECORDED":
        _validate_common_event(event, kind="measurement_completion_record")
        required_counts = {
            "terminal_rows": EXPECTED_FREEZE["frozen_rows"],
            "unique_measurement_ids": EXPECTED_FREEZE["frozen_rows"],
            "terminal_partitions": EXPECTED_FREEZE["terminal_partitions"],
            "duplicate_measurement_ids": 0,
            "missing_measurement_ids": 0,
        }
        for key, expected in required_counts.items():
            if event.get(key) != expected:
                raise RuntimeError(f"measurement completion mismatch: {key}")
        _require_bool(event.get("early_stopping_used"), "early_stopping_used")
        if event["early_stopping_used"]:
            raise RuntimeError("early stopping is forbidden")
        _require_bool(event.get("replacement_used"), "replacement_used")
        if event["replacement_used"]:
            raise RuntimeError("replacement is forbidden")

    elif current == "MEASUREMENT_COMPLETE":
        _validate_common_event(event, kind="inference_opening_record")
        _require_bool(event.get("complete_terminal_census_verified"), "complete_terminal_census_verified")
        if not event["complete_terminal_census_verified"]:
            raise RuntimeError("inference cannot open before complete census verification")
        _require_bool(event.get("metadata_colour_join_opened_before_record"), "metadata_colour_join_opened_before_record")
        if event["metadata_colour_join_opened_before_record"]:
            raise RuntimeError("join was opened before inference-opening record")
        _require_bool(event.get("H2_opened_before_record"), "H2_opened_before_record")
        if event["H2_opened_before_record"]:
            raise RuntimeError("H2 was opened before inference-opening record")

    updated = copy.deepcopy(dict(receipt))
    updated["state"] = next_state
    updated["history"] = prior_history + [dict(event)]
    return updated
