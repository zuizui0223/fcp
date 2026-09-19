#!/usr/bin/env python3
"""Verify the v0.12 multivariate orientation benchmark receipt."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


SUMMARY_KEYS = (
    "max_equal_matched_null_fpr",
    "max_common_vector_null_fpr",
    "min_naive_null_fpr",
    "effect_0_4_spread0_equal_detection_range",
    "effect_0_4_spread0_common_detection_range",
    "effect_0_4_spread1_equal_detection_range",
    "effect_0_4_spread1_common_detection_range",
    "effect_0_8_spread1_equal_detection_range",
    "effect_0_8_spread1_common_detection_range",
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
        raise RuntimeError("expected exactly 24 multivariate benchmark cells")

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

    if observed["summary"]["max_equal_matched_null_fpr"] > 0.05:
        raise RuntimeError("matched-null multivariate null rejection exceeds contract")
    if observed["summary"]["max_common_vector_null_fpr"] > 0.05:
        raise RuntimeError("common-vector null rejection exceeds contract")
    if observed["summary"]["min_naive_null_fpr"] != 1.0:
        raise RuntimeError("pooled confounding benchmark weakened below contract")

    shared_equal = observed["summary"]["effect_0_4_spread0_equal_detection_range"]
    shared_common = observed["summary"]["effect_0_4_spread0_common_detection_range"]
    dispersed_equal = observed["summary"]["effect_0_4_spread1_equal_detection_range"]
    dispersed_common = observed["summary"]["effect_0_4_spread1_common_detection_range"]
    strong_dispersed_equal = observed["summary"]["effect_0_8_spread1_equal_detection_range"]
    strong_dispersed_common = observed["summary"]["effect_0_8_spread1_common_detection_range"]

    if min(shared_equal) < 0.95 or min(shared_common) < 0.95:
        raise RuntimeError("shared-orientation power weakened below contract")
    if min(dispersed_equal) < 0.90:
        raise RuntimeError("multivariate matched-null loses weak dispersed-orientation power")
    if max(dispersed_common) > 0.10:
        raise RuntimeError("common multivariate vector no longer cancels under full spread")
    if min(strong_dispersed_equal) < 0.99:
        raise RuntimeError("strong dispersed-orientation matched-null power weakened")
    if max(strong_dispersed_common) > 0.30:
        raise RuntimeError("common vector unexpectedly robust under full orientation spread")

    print(json.dumps({
        "status": "PASS",
        "cells": 24,
        "schema": expected["schema"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
