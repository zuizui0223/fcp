import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[3]
BENCHMARK = Path(__file__).resolve().parents[1] / "benchmarks" / "performance_surface.py"
RESULT_DIR = ROOT / "results" / "disttrait_performance_surface_v0_3_20260919"
CELLS = RESULT_DIR / "cells.csv"
RECEIPT = RESULT_DIR / "result.json"


def _load_benchmark():
    spec = importlib.util.spec_from_file_location("disttrait_performance_surface", BENCHMARK)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_performance_surface_matches_frozen_receipt() -> None:
    module = _load_benchmark()
    observed = module.run_surface()

    expected_cells = pd.read_csv(CELLS)
    observed_cells = pd.DataFrame(observed["cells"])
    pd.testing.assert_frame_equal(
        observed_cells,
        expected_cells,
        check_exact=False,
        atol=1e-12,
        rtol=0,
        check_dtype=False,
    )

    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    for key, expected in receipt["summary"].items():
        actual = observed["summary"][key]
        if isinstance(expected, list):
            assert actual == pytest.approx(expected, abs=1e-12, rel=0)
        else:
            assert actual == pytest.approx(expected, abs=1e-12, rel=0)

    assert observed["summary"]["max_species_conditioned_null_false_positive_fraction"] <= 0.05
    assert observed["summary"]["min_naive_null_false_positive_fraction"] == 1.0
    assert observed["summary"]["effect_0_8_conditioned_detection_range"] == pytest.approx(
        [0.225, 0.675], abs=1e-12, rel=0
    )
    assert observed["summary"]["effect_1_6_conditioned_detection_range"] == pytest.approx(
        [0.75, 1.0], abs=1e-12, rel=0
    )
    assert observed["summary"]["effect_2_4_conditioned_detection_range"] == pytest.approx(
        [1.0, 1.0], abs=1e-12, rel=0
    )
