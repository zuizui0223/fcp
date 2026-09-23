from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "reassemble_fcp_v2_biological_pass_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_bio_reassemble", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_condition_denominators_are_frozen() -> None:
    m = load_module()
    assert m.NONHEAVY_CONDITIONS == 6
    assert m.HEAVY_CONDITIONS == 18
    assert m.EXPECTED_ROWS == 40_000
    assert m.EXPECTED_PARTITIONS == 256
    assert m.EXPECTED_HEAVY == 8_000


def test_bool_parser_is_stable() -> None:
    m = load_module()
    x = pd.Series([True, False, "true", "1", "yes", "no", ""])
    observed = m._bool_series(x).tolist()
    assert observed == [True, False, True, True, True, False, False]
