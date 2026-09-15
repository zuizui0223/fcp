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


def _base(stage):
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


def test_forward_start_record_requires_runtime_checks():
    sha = "a" * 40
    record = build_start_record(
        branch="analysis/p500-prospective-execution-gate-20260915",
        head_sha=sha,
        known_output_paths_absent=True,
        working_tree_inputs_verified=True,
    )
    validate_start_record(record, expected_head_sha=sha)
    assert record["opening_authorized"] is True

    blocked = build_start_record(
        branch="analysis/p500-prospective-execution-gate-20260915",
        head_sha=sha,
        known_output_paths_absent=False,
        working_tree_inputs_verified=True,
    )
    assert blocked["opening_authorized"] is False
    with pytest.raises(RuntimeError, match="output paths"):
        validate_start_record(blocked, expected_head_sha=sha)


def test_stage_ladder_matches_partition_local_ephemeral_execution():
    pre = _base("PREOPENING")
    firewall = _base("FIREWALL_FROZEN")
    firewall.update(
        measurement_ids=EXPECTED_ROWS,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        candidate_pixels_opened=False,
    )
    validate_stage_transition(pre, firewall)

    measured = _base("PARTITION_MEASUREMENT_COMPLETE")
    measured.update(
        terminal_acquisition_rows=EXPECTED_ROWS,
        terminal_measurement_rows=EXPECTED_ROWS,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        persisted_image_pixels=False,
    )
    validate_stage_transition(firewall, measured)

    reasm = _base("REASSEMBLY_COMPLETE")
    reasm.update(unique_measurement_ids=EXPECTED_ROWS, duplicate_measurement_ids=0)
    validate_stage_transition(measured, reasm)

    support = _base("SUPPORT_GATE_COMPLETE")
    support.update(measurement_evaluable_species=300, support_decision="PASS")
    validate_stage_transition(reasm, support)

    h2 = _base("H2_COMPLETE")
    h2.update(
        primary_h2_vector_species=40,
        primary_structured_null_p=0.031,
        h2_decision="CONFIRMED",
    )
    validate_stage_transition(support, h2)


def test_stage_skipping_is_rejected():
    pre = _base("PREOPENING")
    measured = _base("PARTITION_MEASUREMENT_COMPLETE")
    measured.update(
        terminal_acquisition_rows=EXPECTED_ROWS,
        terminal_measurement_rows=EXPECTED_ROWS,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        persisted_image_pixels=False,
    )
    with pytest.raises(RuntimeError, match="illegal"):
        validate_stage_transition(pre, measured)


def test_incomplete_partition_measurement_is_rejected():
    firewall = _base("FIREWALL_FROZEN")
    measured = _base("PARTITION_MEASUREMENT_COMPLETE")
    measured.update(
        terminal_acquisition_rows=EXPECTED_ROWS,
        terminal_measurement_rows=EXPECTED_ROWS - 1,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        persisted_image_pixels=False,
    )
    with pytest.raises(RuntimeError, match="measurement coverage"):
        validate_stage_transition(firewall, measured)


def test_persisted_pixels_are_rejected():
    firewall = _base("FIREWALL_FROZEN")
    measured = _base("PARTITION_MEASUREMENT_COMPLETE")
    measured.update(
        terminal_acquisition_rows=EXPECTED_ROWS,
        terminal_measurement_rows=EXPECTED_ROWS,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        persisted_image_pixels=True,
    )
    with pytest.raises(RuntimeError, match="persisted"):
        validate_stage_transition(firewall, measured)


def test_support_failure_stops_before_h2():
    reasm = _base("REASSEMBLY_COMPLETE")
    support = _base("SUPPORT_GATE_COMPLETE")
    support.update(measurement_evaluable_species=249, support_decision="NOT_EVALUABLE")
    validate_stage_transition(reasm, support)

    h2 = _base("H2_COMPLETE")
    h2.update(primary_h2_vector_species=30, primary_structured_null_p=0.01, h2_decision="CONFIRMED")
    with pytest.raises(RuntimeError, match="failed support gate"):
        validate_stage_transition(support, h2)


def test_vector_support_is_not_rescued_by_pvalue():
    support = _base("SUPPORT_GATE_COMPLETE")
    support.update(measurement_evaluable_species=300, support_decision="PASS")
    h2 = _base("H2_COMPLETE")
    h2.update(primary_h2_vector_species=19, primary_structured_null_p=0.001, h2_decision="CONFIRMED")
    with pytest.raises(RuntimeError, match="vector-support"):
        validate_stage_transition(support, h2)


def test_replacement_is_never_allowed():
    pre = _base("PREOPENING")
    firewall = _base("FIREWALL_FROZEN")
    firewall.update(
        measurement_ids=EXPECTED_ROWS,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        candidate_pixels_opened=False,
        replacement_rows=1,
    )
    with pytest.raises(RuntimeError, match="replacement"):
        validate_stage_transition(pre, firewall)
