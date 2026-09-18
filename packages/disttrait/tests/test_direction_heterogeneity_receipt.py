import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "results" / "disttrait_direction_heterogeneity_v0_7_20260919"
RESULT = RESULT_DIR / "result.json"
CELLS = RESULT_DIR / "cells.csv"


def test_direction_heterogeneity_receipt_contract() -> None:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    cells = pd.read_csv(CELLS)

    assert result["schema"] == "disttrait_direction_heterogeneity_surface_v1"
    assert len(cells) == 24
    assert set(cells["effect_size"]) == {0.0, 0.4, 0.8, 1.2}
    assert set(cells["reversal_fraction"]) == {0.0, 0.25, 0.5}
    assert set(cells["missing_fraction"]) == {0.0, 0.5}
    assert set(cells["worlds"]) == {40}

    summary = result["summary"]
    assert summary["max_equal_matched_null_fpr"] <= 0.05
    assert summary["min_naive_null_fpr"] == 1.0

    # With shared direction, the correctly aligned common-slope model is stronger
    # at the weak effect.
    assert summary["effect_0_4_reversal0_common_slope_detection_range"] == [1.0, 1.0]
    assert summary["effect_0_4_reversal0_equal_detection_range"][1] <= 0.8

    # With 50% direction reversal, the common signed slope cancels but the
    # direction-invariant matched-null statistic retains meaningful power.
    assert summary["effect_0_4_reversal0_5_common_slope_detection_range"][1] <= 0.10
    assert summary["effect_0_4_reversal0_5_equal_detection_range"][0] >= 0.475
    assert summary["effect_0_8_reversal0_5_common_slope_detection_range"][1] <= 0.125
    assert summary["effect_0_8_reversal0_5_equal_detection_range"] == [1.0, 1.0]

    assert "estimand alignment" in result["interpretation"]["supported"]
