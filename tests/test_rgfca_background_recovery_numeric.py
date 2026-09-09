"""Replay the actual CSV parser failure without loading the image model."""
import ast
from pathlib import Path

import numpy as np
import pytest

SOURCE = Path(__file__).resolve().parents[1] / "scripts/analysis/recover_global_rgfca_background_control_partition.py"


def parser_from_source():
    # Execute the real pure helper without importing Torch/ONNX/image models.
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "_expected_integer")
    namespace = {"np": np}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(SOURCE), "exec"), namespace)
    return namespace["_expected_integer"]


@pytest.mark.parametrize("serialized,expected", [("1818.0", 1818), ("56230.0", 56230), ("0.0", 0), ("100", 100)])
def test_real_integral_csv_values_preserved(serialized, expected):
    parse = parser_from_source()
    assert parse(serialized, label="frozen pixel count") == expected
    if "." in serialized:
        # This is precisely the old call-site failure, not a synthetic ROI test.
        with pytest.raises(ValueError, match="invalid literal for int"):
            int(serialized)


@pytest.mark.parametrize("serialized", ["1818.5", "-1", "nan", "inf", "", "garbage"])
def test_invalid_values_never_rounded_or_imputed(serialized):
    with pytest.raises((RuntimeError, ValueError)):
        parser_from_source()(serialized, label="frozen pixel count")


def test_all_three_frozen_reproduction_counts_use_guard():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name) and node.func.id == "_expected_integer"]
    assert len(calls) == 3  # flower pixels, background pixels, each flower palette count
