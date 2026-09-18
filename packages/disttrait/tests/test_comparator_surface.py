import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[3]
BENCHMARK = Path(__file__).resolve().parents[1] / "benchmarks" / "comparator_surface.py"
RESULT_DIR = ROOT / "results" / "disttrait_comparator_surface_v0_4_20260919"
CELLS = RESULT_DIR / "cells.csv"
RECEIPT = RESULT_DIR / "result.json"


def _load_benchmark():
    spec = importlib.util.spec_from_file_location("disttrait_comparator_surface", BENCHMARK)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_comparator_surface_matches_frozen_receipt() -> None:
    module = _load_benchmark()
    observed = module.run_benchmark()

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

    assert observed["summary"]["max_equal_matched_null_fpr"] <= 0.05
    assert observed["summary"]["max_pair_weighted_matched_null_fpr"] <= 0.05
    assert observed["summary"]["max_species_rho_ttest_null_fpr"] <= 0.05
    assert observed["summary"]["min_naive_null_fpr"] == 1.0

    # Under weak signal plus 4x imbalance, equal weighting is never worse than
    # the pair-weighted matched-null range in this frozen benchmark.
    assert observed["summary"]["effect_0_8_imbalance4_equal_detection_range"][0] >= (
        observed["summary"]["effect_0_8_imbalance4_pair_weighted_detection_range"][0]
    )
    assert observed["summary"]["effect_0_8_imbalance4_equal_detection_range"][1] >= (
        observed["summary"]["effect_0_8_imbalance4_pair_weighted_detection_range"][1]
    )
