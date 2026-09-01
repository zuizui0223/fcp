from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_inference_cascade import (
    BRANCHES,
    STOP_CONFIRMATORY,
    classify_shared_transition,
    classify_spatial_organization,
    load_contract,
    next_step,
    validate_decision_ledger_rows,
    walk_confirmatory_cascade,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"


def _contract():
    return load_contract(CONTRACT_PATH)


def _git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def test_contract_freezes_terminal_dimensions_and_branch_order():
    contract = _contract()
    sampling = contract["terminal_scaleout_sampling"]
    assert sampling["species_total"] == 200
    assert sampling["photo_target_total"] == 60000
    assert sampling["cohort_count"] == 8
    assert sampling["species_per_cohort"] == 25
    assert tuple(contract["confirmatory_sequence"]) == BRANCHES


def test_immutable_legacy_blob_pins_match_repo():
    contract = _contract()
    for relative, expected in contract["immutable_git_blobs"].items():
        assert _git_blob_sha(ROOT / relative) == expected, relative


def test_not_evaluable_never_advances_confirmatory_branch():
    contract = _contract()
    for branch in BRANCHES[:-1]:
        assert next_step(branch, "not_evaluable", contract) == STOP_CONFIRMATORY


def test_spatial_support_requires_pooled_signal_and_six_of_eight_directions():
    contract = _contract()
    assert (
        classify_spatial_organization(
            pooled_p_value=0.01,
            cohort_directions=[1, 1, 1, 1, 1, 1, -1, -1],
            contract=contract,
        )
        == "supported"
    )
    assert (
        classify_spatial_organization(
            pooled_p_value=0.01,
            cohort_directions=[1, 1, 1, 1, 1, -1, -1, -1],
            contract=contract,
        )
        == "unsupported"
    )
    assert (
        classify_spatial_organization(
            pooled_p_value=0.20,
            cohort_directions=[1] * 8,
            contract=contract,
        )
        == "unsupported"
    )


def test_missing_cohort_makes_spatial_branch_not_evaluable():
    contract = _contract()
    assert (
        classify_spatial_organization(
            pooled_p_value=0.01,
            cohort_directions=[1, 1, 1, 1, 1, 1, 1, None],
            contract=contract,
        )
        == "not_evaluable"
    )


def test_shared_transition_requires_preimage_qualification():
    contract = _contract()
    assert (
        classify_shared_transition(
            pooled_p_value=0.001,
            cohort_directions=[1] * 8,
            qualification_passed=False,
            contract=contract,
        )
        == "not_evaluable"
    )
    assert (
        classify_shared_transition(
            pooled_p_value=0.001,
            cohort_directions=[1] * 8,
            qualification_passed=True,
            contract=contract,
        )
        == "supported"
    )


def test_spatial_unsupported_skips_shared_and_runs_environment():
    contract = _contract()
    decisions = walk_confirmatory_cascade(
        {
            "species_conditioned_spatial_organization": "unsupported",
            "environmental_concordance": "unsupported",
            "pollinator_biogeographic_concordance": "supported",
        },
        contract,
    )
    assert [item.branch for item in decisions] == [
        "species_conditioned_spatial_organization",
        "environmental_concordance",
        "pollinator_biogeographic_concordance",
    ]


def test_spatial_supported_then_shared_unsupported_runs_environment():
    contract = _contract()
    decisions = walk_confirmatory_cascade(
        {
            "species_conditioned_spatial_organization": "supported",
            "shared_transition": "unsupported",
            "environmental_concordance": "supported",
        },
        contract,
    )
    assert [item.branch for item in decisions] == [
        "species_conditioned_spatial_organization",
        "shared_transition",
        "environmental_concordance",
    ]


def test_not_evaluable_rejects_posthoc_downstream_confirmatory_outcome():
    contract = _contract()
    with pytest.raises(ValueError, match="not confirmatorily reached"):
        walk_confirmatory_cascade(
            {
                "species_conditioned_spatial_organization": "not_evaluable",
                "environmental_concordance": "supported",
            },
            contract,
        )


def test_supported_shared_transition_stops_before_environment():
    contract = _contract()
    with pytest.raises(ValueError, match="not confirmatorily reached"):
        walk_confirmatory_cascade(
            {
                "species_conditioned_spatial_organization": "supported",
                "shared_transition": "supported",
                "environmental_concordance": "supported",
            },
            contract,
        )


def test_decision_ledger_allows_diagnostic_only_after_not_evaluable_but_not_promotion():
    contract = _contract()
    columns = contract["required_decision_ledger_columns"]

    def row(branch, outcome, state, reason):
        value = {column: "frozen" for column in columns}
        value.update(
            branch=branch,
            outcome=outcome,
            execution_state=state,
            advance_reason=reason,
        )
        return value

    rows = [
        row(
            "species_conditioned_spatial_organization",
            "not_evaluable",
            "stopped_not_evaluable",
            "geometry gate failed",
        ),
        row(
            "environmental_concordance",
            "supported",
            "diagnostic_only",
            "descriptive diagnostic after stopped confirmatory cascade",
        ),
    ]
    validate_decision_ledger_rows(rows, contract)

    rows[1]["advance_reason"] = "confirmatory promotion"
    with pytest.raises(ValueError, match="diagnostic-only"):
        validate_decision_ledger_rows(rows, contract)
