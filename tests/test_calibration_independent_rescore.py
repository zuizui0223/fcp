import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
PATH = ROOT / "scripts" / "data" / "rescore_jbi_ch1_calibration_states.py"
SPEC = importlib.util.spec_from_file_location("independent_rescore", PATH)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def result(state="yellow"):
    return {
        "blind_id": "abc",
        "flower_visibility": "evaluable",
        "visibility_failure_code": "",
        "flower_condition": "fresh",
        "flower_region": "single_target_clear",
        "within_photo_flower_consistency": "single_flower",
        "segmentation_feasibility": "feasible",
        "candidate_state": state,
        "candidate_state_confidence": 0.51,
        "notes": "independent review",
    }


def test_valid_independent_state_passes_without_confidence_threshold():
    mod.validate(result(), "abc", {"yellow", "orange", "unresolved"})


def test_state_outside_frozen_codebook_fails():
    with pytest.raises(ValueError, match="outside frozen codebook"):
        mod.validate(result("green"), "abc", {"yellow", "orange", "unresolved"})


def test_nonfresh_must_be_unresolved():
    row = result()
    row["flower_condition"] = "senescent"
    with pytest.raises(ValueError, match="requires unresolved"):
        mod.validate(row, "abc", {"yellow", "orange", "unresolved"})


def test_variable_multiple_flowers_must_be_unresolved():
    row = result()
    row["flower_region"] = "multiple_flowers_clear"
    row["within_photo_flower_consistency"] = "variable_between_flowers"
    with pytest.raises(ValueError, match="requires unresolved"):
        mod.validate(row, "abc", {"yellow", "orange", "unresolved"})
