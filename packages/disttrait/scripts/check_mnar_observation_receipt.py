#!/usr/bin/env python3
"""Verify the v0.10 MNAR observation benchmark against its frozen receipt."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


SUMMARY_KEYS = (
    "minimum_evaluable_worlds",
    "max_equal_rejection_nonjoint",
    "max_fixed_logit_rejection_nonjoint",
    "joint_equal_rejection_range",
    "joint_fixed_logit_rejection_range",
    "joint_observed_abs_state_position_rho_range",
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

    if len(observed_cells) != 8 or len(expected_cells) != 8:
        raise RuntimeError("expected exactly 8 MNAR benchmark cells")

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

    if int(observed["summary"]["minimum_evaluable_worlds"]) != 40:
        raise RuntimeError("all benchmark cells must retain 40 evaluable worlds")
    if float(observed["summary"]["max_equal_rejection_nonjoint"]) > 0.05 + 1e-12:
        raise RuntimeError("nonjoint matched-null rejection exceeds frozen 0.05 ceiling")
    if float(observed["summary"]["joint_equal_rejection_range"][0]) < 0.80 - 1e-12:
        raise RuntimeError("joint-selection stress effect weakened below frozen contract")

    print(json.dumps({
        "status": "PASS",
        "cells": 8,
        "schema": expected["schema"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
