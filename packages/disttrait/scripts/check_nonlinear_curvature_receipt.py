#!/usr/bin/env python3
"""Verify the v0.11 nonlinear curvature benchmark receipt."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


SUMMARY_KEYS = (
    "max_equal_matched_null_fpr",
    "max_linear_null_fpr",
    "max_quadratic_null_fpr",
    "min_naive_null_fpr",
    "effect_0_4_reversal0_equal_detection_range",
    "effect_0_4_reversal0_linear_detection_range",
    "effect_0_4_reversal0_quadratic_detection_range",
    "effect_0_4_reversal0_5_equal_detection_range",
    "effect_0_4_reversal0_5_quadratic_detection_range",
    "effect_0_8_reversal0_5_equal_detection_range",
    "effect_0_8_reversal0_5_quadratic_detection_range",
)


def _close(a, b, tol: float = 1e-12) -> bool:
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(
            _close(x, y, tol) for x, y in zip(a, b, strict=True)
        )
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    args = parser.parse_args()

    observed = json.loads(args.observed.read_text(encoding="utf-8"))
    expected = json.loads(args.expected.read_text(encoding="utf-8"))
    expected_cells = pd.read_csv(args.cells)
    observed_cells = pd.DataFrame(observed["cells"])

    if len(observed_cells) != 24 or len(expected_cells) != 24:
        raise RuntimeError("expected exactly 24 nonlinear benchmark cells")

    pd.testing.assert_frame_equal(
        observed_cells[expected_cells.columns],
        expected_cells,
        check_exact=False,
        atol=1e-12,
        rtol=0,
        check_dtype=False,
    )

    for key in SUMMARY_KEYS:
        a = observed["summary"][key]
        b = expected["summary"][key]
        if not _close(a, b):
            raise RuntimeError(f"{key} drift: observed={a!r}, expected={b!r}")

    shared_quad = observed["summary"]["effect_0_4_reversal0_quadratic_detection_range"]
    shared_linear = observed["summary"]["effect_0_4_reversal0_linear_detection_range"]
    reversed_equal = observed["summary"]["effect_0_8_reversal0_5_equal_detection_range"]
    reversed_quad = observed["summary"]["effect_0_8_reversal0_5_quadratic_detection_range"]

    if min(shared_quad) < 0.95:
        raise RuntimeError("shared-curvature quadratic power weakened below contract")
    if max(shared_linear) > 0.05:
        raise RuntimeError("common linear model unexpectedly detects weak quadratic signal")
    if min(reversed_equal) < 0.80:
        raise RuntimeError("matched-null power under curvature reversal weakened below contract")
    if max(reversed_quad) > 0.10:
        raise RuntimeError("common quadratic model no longer cancels under 50% reversal")

    print(json.dumps({
        "status": "PASS",
        "cells": 24,
        "schema": expected["schema"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
