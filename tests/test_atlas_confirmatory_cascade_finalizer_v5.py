from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.data.finalize_jbi_atlas_confirmatory_cascade_v5 import (
    BOMBUS_NE,
    build_final_decision,
)


INFERENCE = Path("docs/supporting/jbi_image_first_atlas_inference_contract_v5.json")


def _inference() -> dict:
    return json.loads(INFERENCE.read_text(encoding="utf-8"))


def _gate(status: str) -> dict:
    passed = status == "pass_scaleout_measurement_completeness"
    return {
        "status": status,
        "frozen_measurements": 60000,
        "coordinates_opened": False,
        "coordinate_join_permitted": passed,
    }


def _primary(branches: list[tuple[str, str]], next_branch: str, colour_hash: str = "abc") -> dict:
    return {
        "protocol": "jbi-atlas-real-colour-inference-amendment-v5-v1",
        "status": "complete_reached_primary_spatial_branches_v5",
        "coordinate_colour_join_performed": True,
        "same_standardized_colour_field_required_for_all_downstream_branches": True,
        "standardized_colour_field_sha256": colour_hash,
        "branches": [{"branch": branch, "outcome": outcome} for branch, outcome in branches],
        "next_confirmatory_branch": next_branch,
    }


def _environment(outcome: str, colour_hash: str = "abc") -> dict:
    return {
        "protocol": "jbi-atlas-environmental-concordance-v5-v1",
        "status": "complete_environmental_concordance_v5",
        "branch": "environmental_concordance",
        "outcome": outcome,
        "standardized_flower_field_remeasured": False,
        "standardized_flower_field_restandardized": False,
        "standardized_flower_colour_field_sha256": colour_hash,
    }


def test_measurement_not_evaluable_stops_without_coordinate_join() -> None:
    result = build_final_decision(
        measurement_gate=_gate("not_evaluable_scaleout_measurement_completeness"),
        inference=_inference(),
    )
    assert result["coordinate_colour_join_performed"] is False
    assert result["final_confirmatory_state"] == "STOP_CONFIRMATORY"
    assert result["branch_ledger"] == [
        {
            "branch": "species_conditioned_spatial_organization",
            "outcome": "not_evaluable",
            "next_step": "STOP_CONFIRMATORY",
        }
    ]


def test_measurement_not_evaluable_rejects_downstream_result() -> None:
    with pytest.raises(ValueError, match="downstream evidence"):
        build_final_decision(
            measurement_gate=_gate("not_evaluable_scaleout_measurement_completeness"),
            inference=_inference(),
            primary=_primary(
                [("species_conditioned_spatial_organization", "not_evaluable")],
                "STOP_CONFIRMATORY",
            ),
        )


def test_spatial_supported_shared_supported_stops() -> None:
    result = build_final_decision(
        measurement_gate=_gate("pass_scaleout_measurement_completeness"),
        inference=_inference(),
        primary=_primary(
            [
                ("species_conditioned_spatial_organization", "supported"),
                ("shared_transition", "supported"),
            ],
            "STOP_CONFIRMATORY",
        ),
    )
    assert result["supported_branches"] == [
        "species_conditioned_spatial_organization",
        "shared_transition",
    ]
    assert result["final_confirmatory_state"] == "STOP_CONFIRMATORY"


def test_spatial_unsupported_environment_supported_stops() -> None:
    result = build_final_decision(
        measurement_gate=_gate("pass_scaleout_measurement_completeness"),
        inference=_inference(),
        primary=_primary(
            [("species_conditioned_spatial_organization", "unsupported")],
            "environmental_concordance",
        ),
        environmental=_environment("supported"),
    )
    assert [row["branch"] for row in result["branch_ledger"]] == [
        "species_conditioned_spatial_organization",
        "environmental_concordance",
    ]
    assert result["supported_branches"] == ["environmental_concordance"]
    assert result["final_confirmatory_state"] == "STOP_CONFIRMATORY"


def test_shared_unsupported_environment_unsupported_reaches_binding_pollinator_ne() -> None:
    result = build_final_decision(
        measurement_gate=_gate("pass_scaleout_measurement_completeness"),
        inference=_inference(),
        primary=_primary(
            [
                ("species_conditioned_spatial_organization", "supported"),
                ("shared_transition", "unsupported"),
            ],
            "environmental_concordance",
        ),
        environmental=_environment("unsupported"),
        bombus_gate={"status": BOMBUS_NE},
    )
    assert [row["outcome"] for row in result["branch_ledger"]] == [
        "supported",
        "unsupported",
        "unsupported",
        "not_evaluable",
    ]
    assert result["final_confirmatory_state"] == "TERMINAL"
    assert result["post_colour_pollinator_source_substitution_used"] is False


def test_environment_cannot_be_supplied_when_shared_supported() -> None:
    with pytest.raises(ValueError, match="not reached"):
        build_final_decision(
            measurement_gate=_gate("pass_scaleout_measurement_completeness"),
            inference=_inference(),
            primary=_primary(
                [
                    ("species_conditioned_spatial_organization", "supported"),
                    ("shared_transition", "supported"),
                ],
                "STOP_CONFIRMATORY",
            ),
            environmental=_environment("supported"),
        )


def test_environment_must_reuse_exact_flower_hash() -> None:
    with pytest.raises(ValueError, match="frozen v5 colour field"):
        build_final_decision(
            measurement_gate=_gate("pass_scaleout_measurement_completeness"),
            inference=_inference(),
            primary=_primary(
                [("species_conditioned_spatial_organization", "unsupported")],
                "environmental_concordance",
                colour_hash="one",
            ),
            environmental=_environment("supported", colour_hash="two"),
        )
