import copy
import json
from pathlib import Path

import pytest

from fcp_pipeline.p500_prospective_execution_gate import (
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    EXPECTED_TERMINAL_PARTITIONS,
    build_start_record,
    validate_authorization,
    validate_candidate_metadata,
    validate_stage_transition,
    validate_start_record,
)

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "results/polymorphism_h2_p500_candidate_metadata_20260913/result.json"
AUTH = ROOT / "docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_AUTHORIZATION_20260915.json"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _base_record(stage):
    return {
        "stage": stage,
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
    }


def test_real_frozen_candidate_and_authorization_are_accepted():
    validate_candidate_metadata(_load(CANDIDATE))
    validate_authorization(_load(AUTH))


def test_candidate_firewall_change_is_rejected():
    candidate = _load(CANDIDATE)
    candidate["outcome_firewall"]["image_pixels_opened"] = True
    with pytest.raises(RuntimeError, match="firewall"):
        validate_candidate_metadata(candidate)


def test_authorization_estimand_change_is_rejected():
    auth = _load(AUTH)
    auth["axis_refit_allowed"] = True
    with pytest.raises(RuntimeError, match="axis refit"):
        validate_authorization(auth)


def test_forward_start_record_can_authorize_only_after_runtime_checks():
    sha = "a" * 40
    record = build_start_record(
        branch="analysis/p500-prospective-execution-gate-20260915",
        head_sha=sha,
        known_output_paths_absent=True,
        working_tree_inputs_verified=True,
    )
    assert record["opening_authorized"] is True
    validate_start_record(record, expected_head_sha=sha)


def test_forward_start_record_fails_closed_without_runtime_checks():
    sha = "b" * 40
    record = build_start_record(
        branch="analysis/p500-prospective-execution-gate-20260915",
        head_sha=sha,
        known_output_paths_absent=False,
        working_tree_inputs_verified=True,
    )
    assert record["opening_authorized"] is False
    with pytest.raises(RuntimeError, match="output paths"):
        validate_start_record(record, expected_head_sha=sha)


def test_stage_ladder_accepts_complete_nominal_path():
    pre = _base_record("PREOPENING")
    acq = _base_record("ACQUISITION_COMPLETE")
    acq["terminal_acquisition_rows"] = EXPECTED_ROWS
    validate_stage_transition(pre, acq)

    meas = _base_record("MEASUREMENT_COMPLETE")
    meas["terminal_measurement_rows"] = EXPECTED_ROWS
    meas["terminal_partitions"] = EXPECTED_TERMINAL_PARTITIONS
    validate_stage_transition(acq, meas)

    reasm = _base_record("REASSEMBLY_COMPLETE")
    reasm["unique_measurement_ids"] = EXPECTED_ROWS
    reasm["duplicate_measurement_ids"] = 0
    validate_stage_transition(meas, reasm)

    support = _base_record("SUPPORT_GATE_COMPLETE")
    support["measurement_evaluable_species"] = 300
    support["support_decision"] = "PASS"
    validate_stage_transition(reasm, support)

    h2 = _base_record("H2_COMPLETE")
    h2["primary_h2_vector_species"] = 40
    h2["primary_structured_null_p"] = 0.031
    h2["h2_decision"] = "CONFIRMED"
    validate_stage_transition(support, h2)


def test_stage_skipping_is_rejected():
    pre = _base_record("PREOPENING")
    meas = _base_record("MEASUREMENT_COMPLETE")
    meas["terminal_measurement_rows"] = EXPECTED_ROWS
    meas["terminal_partitions"] = EXPECTED_TERMINAL_PARTITIONS
    with pytest.raises(RuntimeError, match="illegal"):
        validate_stage_transition(pre, meas)


def test_incomplete_measurement_denominator_is_rejected():
    acq = _base_record("ACQUISITION_COMPLETE")
    meas = _base_record("MEASUREMENT_COMPLETE")
    meas["terminal_measurement_rows"] = EXPECTED_ROWS - 1
    meas["terminal_partitions"] = EXPECTED_TERMINAL_PARTITIONS
    with pytest.raises(RuntimeError, match="all frozen rows"):
        validate_stage_transition(acq, meas)


def test_support_failure_stops_before_h2():
    reasm = _base_record("REASSEMBLY_COMPLETE")
    support = _base_record("SUPPORT_GATE_COMPLETE")
    support["measurement_evaluable_species"] = 249
    support["support_decision"] = "NOT_EVALUABLE"
    validate_stage_transition(reasm, support)

    h2 = _base_record("H2_COMPLETE")
    h2["primary_h2_vector_species"] = 30
    h2["primary_structured_null_p"] = 0.01
    h2["h2_decision"] = "CONFIRMED"
    with pytest.raises(RuntimeError, match="failed support gate"):
        validate_stage_transition(support, h2)


def test_vector_support_is_not_rescued_by_pvalue():
    support = _base_record("SUPPORT_GATE_COMPLETE")
    support["measurement_evaluable_species"] = 300
    support["support_decision"] = "PASS"
    h2 = _base_record("H2_COMPLETE")
    h2["primary_h2_vector_species"] = 19
    h2["primary_structured_null_p"] = 0.001
    h2["h2_decision"] = "CONFIRMED"
    with pytest.raises(RuntimeError, match="vector-support"):
        validate_stage_transition(support, h2)


def test_replacement_is_never_allowed():
    pre = _base_record("PREOPENING")
    acq = _base_record("ACQUISITION_COMPLETE")
    acq["terminal_acquisition_rows"] = EXPECTED_ROWS
    acq["replacement_rows"] = 1
    with pytest.raises(RuntimeError, match="replacement"):
        validate_stage_transition(pre, acq)
