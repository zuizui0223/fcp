import importlib.util
from pathlib import Path

import pytest


PATH = Path(__file__).parents[1] / "scripts" / "data" / "screen_jbi_ch1_calibration_species.py"
SPEC = importlib.util.spec_from_file_location("semantic_screen", PATH)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def base_result():
    return {
        "blind_id": "abc123",
        "flower_visibility": "evaluable",
        "visibility_failure_code": "",
        "flower_condition": "fresh",
        "flower_region": "single_target_clear",
        "within_photo_flower_consistency": "single_flower",
        "segmentation_feasibility": "feasible",
        "candidate_state": "yellow",
        "candidate_state_confidence": 0.9,
        "apparent_petals_colour_terms": ["yellow"],
        "notes": "clear fresh flower",
    }


def test_valid_fresh_candidate_state_passes():
    mod.validate(base_result(), "abc123", {"yellow", "orange", "unresolved"})


def test_state_outside_frozen_codebook_is_rejected():
    row = base_result()
    row["candidate_state"] = "green"
    with pytest.raises(ValueError, match="frozen codebook"):
        mod.validate(row, "abc123", {"yellow", "orange", "unresolved"})


def test_senescent_flower_must_be_unresolved():
    row = base_result()
    row["flower_condition"] = "senescent"
    with pytest.raises(ValueError, match="requires candidate_state=unresolved"):
        mod.validate(row, "abc123", {"yellow", "orange", "unresolved"})
    row["candidate_state"] = "unresolved"
    mod.validate(row, "abc123", {"yellow", "orange", "unresolved"})


def test_discordant_multiple_flowers_must_be_unresolved():
    row = base_result()
    row["flower_region"] = "multiple_flowers_clear"
    row["within_photo_flower_consistency"] = "variable_between_flowers"
    with pytest.raises(ValueError, match="requires candidate_state=unresolved"):
        mod.validate(row, "abc123", {"yellow", "orange", "unresolved"})


def test_uncertain_segmentation_must_be_unresolved():
    row = base_result()
    row["segmentation_feasibility"] = "uncertain"
    with pytest.raises(ValueError, match="requires candidate_state=unresolved"):
        mod.validate(row, "abc123", {"yellow", "orange", "unresolved"})
