#!/usr/bin/env python3
"""Verify the frozen ShareTrait Gammarus empirical transport receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _close(a: float, b: float, tol: float = 1e-12) -> None:
    if not np.isclose(float(a), float(b), rtol=0.0, atol=tol):
        raise RuntimeError(f"value mismatch: {a} != {b}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--observed", type=Path, required=True)
    p.add_argument("--expected", type=Path, required=True)
    args = p.parse_args()

    observed = json.loads(args.observed.read_text(encoding="utf-8"))
    expected = json.loads(args.expected.read_text(encoding="utf-8"))

    if observed["schema"] != expected["schema"]:
        raise RuntimeError("schema drift")
    for key in ("rows", "individuals", "populations", "sites"):
        if int(observed["fixture"][key]) != int(expected["fixture"][key]):
            raise RuntimeError(f"fixture denominator drift: {key}")

    for section in ("primary_log_rate", "raw_rate_sensitivity"):
        for key in (
            "observed_rho",
            "p_upper",
            "null_mean",
            "null_q025",
            "null_q975",
        ):
            _close(observed[section][key], expected[section][key])
        if int(observed[section]["permutations"]) != int(expected[section]["permutations"]):
            raise RuntimeError(f"permutation count drift: {section}")

    expected_counts = {
        row["population_id"]: int(row["n"])
        for row in expected["population_descriptives"]
    }
    observed_counts = {
        row["population_id"]: int(row["n"])
        for row in observed["population_descriptives"]
    }
    if observed_counts != expected_counts:
        raise RuntimeError("population denominator drift")

    # The frozen result is intentionally a null empirical transport.
    if float(observed["primary_log_rate"]["p_upper"]) < 0.05:
        raise RuntimeError("primary result changed from frozen non-support")
    if float(observed["raw_rate_sensitivity"]["p_upper"]) < 0.05:
        raise RuntimeError("raw sensitivity changed from frozen non-support")

    print(json.dumps({
        "status": "PASS",
        "primary_rho": observed["primary_log_rate"]["observed_rho"],
        "primary_p": observed["primary_log_rate"]["p_upper"],
        "raw_rho": observed["raw_rate_sensitivity"]["observed_rho"],
        "raw_p": observed["raw_rate_sensitivity"]["p_upper"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
