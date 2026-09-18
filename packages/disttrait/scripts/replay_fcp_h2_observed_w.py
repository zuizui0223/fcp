#!/usr/bin/env python3
"""Replay frozen FCP observed white-axis W values from legacy delta-vector artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from disttrait import alignment_statistic, one_vs_rest_contrast

PALETTE = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
DELTA_COLUMNS = [f"delta_{x}" for x in PALETTE]


def observed_w(path: Path) -> tuple[int, float]:
    frame = pd.read_csv(path)
    missing = sorted(set(DELTA_COLUMNS) - set(frame.columns))
    if missing:
        raise RuntimeError(f"{path}: missing delta columns {missing}")
    axes = frame[DELTA_COLUMNS].to_numpy(float)
    if np.any(~np.isfinite(axes)):
        raise RuntimeError(f"{path}: non-finite delta values")
    q = one_vs_rest_contrast(len(PALETTE), focal_index=0)
    return int(len(frame)), alignment_statistic(axes, q)


def check(label: str, path: Path, expected: dict, tol: float) -> dict:
    n, w = observed_w(path)
    tier, cohort = label.split(":", 1)
    target = expected["thresholds"][tier][cohort]
    expected_n = int(target["species"])
    expected_w = float(target["observed_mean_squared_white_axis_alignment"])
    if n != expected_n:
        raise RuntimeError(f"{label}: species={n}, expected={expected_n}")
    if not np.isclose(w, expected_w, rtol=0.0, atol=tol):
        raise RuntimeError(f"{label}: W={w}, expected={expected_w}")
    return {"species": n, "W": w}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--primary-discovery", type=Path, required=True)
    p.add_argument("--primary-reserve", type=Path, required=True)
    p.add_argument("--strict-discovery", type=Path, required=True)
    p.add_argument("--strict-reserve", type=Path, required=True)
    p.add_argument("--expected", type=Path, required=True)
    p.add_argument("--tolerance", type=float, default=1e-12)
    args = p.parse_args()

    expected = json.loads(args.expected.read_text(encoding="utf-8"))
    out = {
        "primary_0_10:discovery": check(
            "primary_0_10:discovery", args.primary_discovery, expected, args.tolerance
        ),
        "primary_0_10:reserve": check(
            "primary_0_10:reserve", args.primary_reserve, expected, args.tolerance
        ),
        "strict_0_20:discovery": check(
            "strict_0_20:discovery", args.strict_discovery, expected, args.tolerance
        ),
        "strict_0_20:reserve": check(
            "strict_0_20:reserve", args.strict_reserve, expected, args.tolerance
        ),
    }
    print(json.dumps({"status": "PASS", "results": out}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
