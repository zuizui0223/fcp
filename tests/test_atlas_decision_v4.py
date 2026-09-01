from __future__ import annotations

import json
from pathlib import Path

from fcp_pipeline.atlas_decision_v4 import strict_cascade_decisions, validate_v4_contract

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v4.json"


def test_committed_v4_contract_freezes_legacy_sampling_roi_and_branch_order() -> None:
    validate_v4_contract(json.loads(CONTRACT.read_text(encoding="utf-8")), root=ROOT)


def test_not_evaluable_spatial_branch_stops_confirmatory_fallback() -> None:
    result = strict_cascade_decisions({
        "species_conditioned_shared_transition": {"outcome": "not_evaluable", "reason_code": "signal_recovery_failed"},
        "environmental_boundary_concordance": {"outcome": "supported", "reason_code": "diagnostic_environment"},
    })
    assert result["terminal_outcome"] == "not_evaluable"
    assert result["promoted_conclusion"] is None
    assert result["confirmatory_results"]["environmental_boundary_concordance"]["outcome"] == "not_evaluable"
    assert result["diagnostic_results"]["environmental_boundary_concordance"]["outcome"] == "supported"
    assert result["diagnostic_results"]["environmental_boundary_concordance"]["confirmatory_status"] == "diagnostic_only"


def test_only_unsupported_authorizes_next_branch() -> None:
    result = strict_cascade_decisions({
        "species_conditioned_shared_transition": {"outcome": "unsupported"},
        "environmental_boundary_concordance": {"outcome": "supported"},
        "pollinator_biogeographic_concordance": {"outcome": "supported"},
    })
    assert result["terminal_outcome"] == "supported"
    assert result["promoted_conclusion"] == "environmental_boundary_concordance"
    assert result["confirmatory_results"]["species_conditioned_shared_transition"]["outcome"] == "unsupported"
    assert result["confirmatory_results"]["pollinator_biogeographic_concordance"]["outcome"] == "not_evaluable"
    assert result["diagnostic_results"]["pollinator_biogeographic_concordance"]["confirmatory_status"] == "diagnostic_only"


def test_pollinator_is_reached_only_after_two_unsupported_branches() -> None:
    result = strict_cascade_decisions({
        "species_conditioned_shared_transition": {"outcome": "unsupported"},
        "environmental_boundary_concordance": {"outcome": "unsupported"},
        "pollinator_biogeographic_concordance": {"outcome": "supported"},
    })
    assert result["terminal_outcome"] == "supported"
    assert result["promoted_conclusion"] == "pollinator_biogeographic_concordance"


def test_complete_evaluable_negative_tree_is_unsupported() -> None:
    result = strict_cascade_decisions({
        "species_conditioned_shared_transition": {"outcome": "unsupported"},
        "environmental_boundary_concordance": {"outcome": "unsupported"},
        "pollinator_biogeographic_concordance": {"outcome": "unsupported"},
    })
    assert result["terminal_outcome"] == "unsupported"
    assert result["promoted_conclusion"] is None
