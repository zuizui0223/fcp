"""Forward-only execution gate for the prospective P500 H2 confirmation.

This module does not open images or compute flower-colour outcomes.  It validates
frozen metadata/authorization inputs, creates a forward execution chronology,
and enforces the permitted stage order for the future P500 run.
"""

from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
from typing import Any


P500_SHA256 = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"
METADATA_SHA256 = "a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595"
FULL100_SHA256 = "568a64e692dc21f4a5f76c5eadb615a63b07d71f01bc9cd91e8496ccd0fa3f9e"
MEASUREMENT_SOURCE_COMMIT = "9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
EXPECTED_SPECIES = 499
EXPECTED_ROWS = 49_900
EXPECTED_TERMINAL_PARTITIONS = 256
MIN_MEASUREMENT_EVALUABLE_SPECIES = 250
MIN_H2_VECTOR_SPECIES = 20

STAGES = (
    "PREOPENING",
    "ACQUISITION_COMPLETE",
    "MEASUREMENT_COMPLETE",
    "REASSEMBLY_COMPLETE",
    "SUPPORT_GATE_COMPLETE",
    "H2_COMPLETE",
)


def _require_bool(value: Any, label: str) -> None:
    if type(value) is not bool:
        raise RuntimeError(f"{label} must be boolean")


def validate_candidate_metadata(candidate: Mapping[str, Any]) -> None:
    if candidate.get("analysis") != "polymorphism_h2_p500_candidate_metadata_acquisition":
        raise RuntimeError("unexpected candidate metadata analysis")
    if candidate.get("status") != "complete_fresh_metadata_draw_before_pixels":
        raise RuntimeError("candidate metadata was not frozen before pixels")
    if candidate.get("decision", {}).get("verdict") != "P500_CANDIDATE_METADATA_GATE_PASS":
        raise RuntimeError("candidate metadata gate did not pass")
    if candidate.get("P500", {}).get("sha256") != P500_SHA256:
        raise RuntimeError("P500 identity mismatch")
    if candidate.get("capacity", {}).get("full100_species") != EXPECTED_SPECIES:
        raise RuntimeError("full100 species count mismatch")
    if candidate.get("capacity", {}).get("selected_metadata_rows_full100_species") != EXPECTED_ROWS:
        raise RuntimeError("full100 row count mismatch")
    files = candidate.get("files", {})
    if files.get("metadata_sha256") != METADATA_SHA256:
        raise RuntimeError("candidate metadata hash mismatch")
    if files.get("full100_species_sha256") != FULL100_SHA256:
        raise RuntimeError("full100 species hash mismatch")
    if candidate.get("transport", {}).get("request_errors") != 0:
        raise RuntimeError("candidate metadata contains request errors")
    decision = candidate.get("decision", {})
    if decision.get("no_species_replacement") is not True:
        raise RuntimeError("species replacement prohibition missing")
    if decision.get("no_target_relaxation") is not True:
        raise RuntimeError("target relaxation prohibition missing")
    firewall = candidate.get("outcome_firewall", {})
    required_false = (
        "image_pixels_opened",
        "flower_colour_opened",
        "palette_opened",
        "D_opened",
        "H2_W_opened",
        "morph_opened",
    )
    for key in required_false:
        if firewall.get(key) is not False:
            raise RuntimeError(f"candidate metadata firewall not closed: {key}")


def validate_authorization(auth: Mapping[str, Any]) -> None:
    if auth.get("status") != "authorize_exactly_one_prospective_h2_p500_location_blind_measurement":
        raise RuntimeError("unexpected authorization status")
    if auth.get("P500_sha256") != P500_SHA256:
        raise RuntimeError("authorization P500 hash mismatch")
    if auth.get("metadata_sha256") != METADATA_SHA256:
        raise RuntimeError("authorization metadata hash mismatch")
    if auth.get("full100_species_sha256") != FULL100_SHA256:
        raise RuntimeError("authorization full100 hash mismatch")
    if auth.get("frozen_species") != EXPECTED_SPECIES or auth.get("frozen_rows") != EXPECTED_ROWS:
        raise RuntimeError("authorization denominator mismatch")
    if auth.get("species_replacement_allowed") is not False:
        raise RuntimeError("authorization permits species replacement")
    if auth.get("row_replacement_allowed") is not False:
        raise RuntimeError("authorization permits row replacement")
    if auth.get("target_relaxation_allowed") is not False:
        raise RuntimeError("authorization permits target relaxation")
    if auth.get("measurement_machine_source_commit") != MEASUREMENT_SOURCE_COMMIT:
        raise RuntimeError("measurement source commit mismatch")
    if auth.get("H2_target") != "fixed_q_white_W_structured_null":
        raise RuntimeError("H2 target mismatch")
    if auth.get("axis_refit_allowed") is not False:
        raise RuntimeError("axis refit must remain forbidden")
    if auth.get("primary_threshold") != 0.10 or auth.get("strict_sensitivity_threshold") != 0.20:
        raise RuntimeError("H2 threshold mismatch")
    if auth.get("structured_null_replicates") != 999:
        raise RuntimeError("structured-null replicate mismatch")
    if auth.get("minimum_classifiable_photos_per_species") != 40:
        raise RuntimeError("classifiable-photo threshold mismatch")
    if auth.get("minimum_measurement_evaluable_species") != MIN_MEASUREMENT_EVALUABLE_SPECIES:
        raise RuntimeError("measurement-support threshold mismatch")
    if auth.get("minimum_primary_H2_vector_species") != MIN_H2_VECTOR_SPECIES:
        raise RuntimeError("H2 vector threshold mismatch")


def build_start_record(*, branch: str, head_sha: str, known_output_paths_absent: bool,
                       working_tree_inputs_verified: bool) -> dict[str, Any]:
    """Create a forward chronology receipt immediately before first pixel opening.

    The two booleans must be established by the execution workflow itself.  They
    are deliberately not inferred from historical receipts.
    """
    _require_bool(known_output_paths_absent, "known_output_paths_absent")
    _require_bool(working_tree_inputs_verified, "working_tree_inputs_verified")
    if not head_sha or len(head_sha) != 40:
        raise RuntimeError("head_sha must be a 40-character commit SHA")
    now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    authorized = known_output_paths_absent and working_tree_inputs_verified
    return {
        "schema": "p500_prospective_execution_start_v1",
        "stage": "PREOPENING",
        "created_at_utc": now,
        "branch": branch,
        "head_sha": head_sha,
        "P500_sha256": P500_SHA256,
        "metadata_sha256": METADATA_SHA256,
        "full100_species_sha256": FULL100_SHA256,
        "measurement_source_commit": MEASUREMENT_SOURCE_COMMIT,
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "image_pixels_opened_before_start": False,
        "flower_colour_opened_before_start": False,
        "H2_opened_before_start": False,
        "known_output_paths_absent": known_output_paths_absent,
        "working_tree_inputs_verified": working_tree_inputs_verified,
        "opening_authorized": authorized,
        "claim_scope": "forward chronology only; no retrospective certification of unrecorded external access",
    }


def validate_start_record(record: Mapping[str, Any], *, expected_head_sha: str) -> None:
    if record.get("schema") != "p500_prospective_execution_start_v1":
        raise RuntimeError("unexpected start-record schema")
    if record.get("stage") != "PREOPENING":
        raise RuntimeError("start record must be PREOPENING")
    if record.get("head_sha") != expected_head_sha:
        raise RuntimeError("start record head SHA mismatch")
    if record.get("P500_sha256") != P500_SHA256:
        raise RuntimeError("start record P500 mismatch")
    if record.get("metadata_sha256") != METADATA_SHA256:
        raise RuntimeError("start record metadata mismatch")
    if record.get("full100_species_sha256") != FULL100_SHA256:
        raise RuntimeError("start record full100 mismatch")
    if record.get("measurement_source_commit") != MEASUREMENT_SOURCE_COMMIT:
        raise RuntimeError("start record measurement source mismatch")
    for key in (
        "image_pixels_opened_before_start",
        "flower_colour_opened_before_start",
        "H2_opened_before_start",
    ):
        if record.get(key) is not False:
            raise RuntimeError(f"preopening outcome flag not false: {key}")
    if record.get("known_output_paths_absent") is not True:
        raise RuntimeError("known P500 output paths were not verified absent")
    if record.get("working_tree_inputs_verified") is not True:
        raise RuntimeError("working-tree inputs were not verified")
    if record.get("opening_authorized") is not True:
        raise RuntimeError("start record does not authorize opening")
    raw_time = record.get("created_at_utc")
    if not isinstance(raw_time, str) or not raw_time.endswith("Z"):
        raise RuntimeError("invalid start-record timestamp")
    datetime.fromisoformat(raw_time[:-1] + "+00:00")


def validate_stage_transition(previous: Mapping[str, Any], current: Mapping[str, Any]) -> None:
    """Fail closed unless execution advances exactly one frozen stage."""
    prev = previous.get("stage")
    cur = current.get("stage")
    if prev not in STAGES or cur not in STAGES:
        raise RuntimeError("unknown P500 execution stage")
    if STAGES.index(cur) != STAGES.index(prev) + 1:
        raise RuntimeError(f"illegal P500 stage transition: {prev} -> {cur}")
    if current.get("species") != EXPECTED_SPECIES or current.get("rows") != EXPECTED_ROWS:
        raise RuntimeError("P500 denominator changed across stages")
    if current.get("replacement_rows", 0) != 0 or current.get("replacement_species", 0) != 0:
        raise RuntimeError("replacement is forbidden")

    if cur == "ACQUISITION_COMPLETE":
        if current.get("terminal_acquisition_rows") != EXPECTED_ROWS:
            raise RuntimeError("acquisition is not complete for all frozen rows")
    elif cur == "MEASUREMENT_COMPLETE":
        if current.get("terminal_measurement_rows") != EXPECTED_ROWS:
            raise RuntimeError("measurement is not complete for all frozen rows")
        if current.get("terminal_partitions") != EXPECTED_TERMINAL_PARTITIONS:
            raise RuntimeError("terminal partition census mismatch")
    elif cur == "REASSEMBLY_COMPLETE":
        if current.get("unique_measurement_ids") != EXPECTED_ROWS:
            raise RuntimeError("reassembly measurement-ID census mismatch")
        if current.get("duplicate_measurement_ids") != 0:
            raise RuntimeError("duplicate measurement IDs in reassembly")
    elif cur == "SUPPORT_GATE_COMPLETE":
        n = current.get("measurement_evaluable_species")
        if type(n) is not int or not 0 <= n <= EXPECTED_SPECIES:
            raise RuntimeError("invalid measurement-evaluable species count")
        decision = current.get("support_decision")
        expected = "PASS" if n >= MIN_MEASUREMENT_EVALUABLE_SPECIES else "NOT_EVALUABLE"
        if decision != expected:
            raise RuntimeError("measurement-support decision mismatch")
    elif cur == "H2_COMPLETE":
        if previous.get("support_decision") != "PASS":
            raise RuntimeError("H2 cannot open after a failed support gate")
        nvec = current.get("primary_h2_vector_species")
        if type(nvec) is not int or nvec < 0:
            raise RuntimeError("invalid H2 vector count")
        if nvec < MIN_H2_VECTOR_SPECIES:
            if current.get("h2_decision") != "NOT_EVALUABLE_VECTOR_SUPPORT":
                raise RuntimeError("H2 vector-support decision mismatch")
        else:
            p = current.get("primary_structured_null_p")
            if not isinstance(p, (int, float)) or not 0 <= float(p) <= 1:
                raise RuntimeError("invalid H2 p-value")
            expected = "CONFIRMED" if float(p) < 0.05 else "NOT_CONFIRMED"
            if current.get("h2_decision") != expected:
                raise RuntimeError("H2 decision mismatch")


__all__ = [
    "EXPECTED_ROWS",
    "EXPECTED_SPECIES",
    "EXPECTED_TERMINAL_PARTITIONS",
    "FULL100_SHA256",
    "MEASUREMENT_SOURCE_COMMIT",
    "METADATA_SHA256",
    "MIN_H2_VECTOR_SPECIES",
    "MIN_MEASUREMENT_EVALUABLE_SPECIES",
    "P500_SHA256",
    "STAGES",
    "build_start_record",
    "validate_authorization",
    "validate_candidate_metadata",
    "validate_stage_transition",
    "validate_start_record",
]
