#!/usr/bin/env python3
"""Replay the frozen FCP H1 observer-disjoint reliability with disttrait.

This script is intentionally artifact-backed: the clean submission branch does
not carry the original 100,000-row discovery/reserve measured tables.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from disttrait import observer_disjoint_reliability

STATES = ["white", "yellow_orange", "red_pink", "blue_purple"]


def _assert_close(label: str, observed: float, expected: float, tol: float) -> None:
    if not np.isclose(observed, expected, rtol=0.0, atol=tol, equal_nan=True):
        raise RuntimeError(f"{label}: observed={observed!r}, expected={expected!r}")


def replay(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, low_memory=False)
    required = {"species", "observer_id", "global_classifiable", "morph"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"{path}: missing columns {missing}")
    result = observer_disjoint_reliability(
        frame,
        species_col="species",
        observer_col="observer_id",
        state_col="morph",
        classifiable_col="global_classifiable",
        states=STATES,
        n_partitions=200,
        base_seed=20260913,
        min_full_classifiable=40,
        min_half_classifiable=20,
        minimum_distinct_observers=2,
    )
    return result.summary


def check_cohort(
    cohort: str,
    observed: dict[str, float],
    expected: dict,
    *,
    tol: float,
) -> None:
    target = expected[cohort]["primary20"]
    mapping = {
        "paired_n_median": "paired_n_median",
        "rho_median": "rho_median",
        "rho_q05": "rho_q05",
        "rho_q95": "rho_q95",
        "ccc_median": "ccc_median",
        "spearman_brown_median": "spearman_brown_median",
        "mae_median": "mae_median",
    }
    for ours, theirs in mapping.items():
        _assert_close(f"{cohort}.{ours}", observed[ours], float(target[theirs]), tol)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--discovery", type=Path, required=True)
    p.add_argument("--reserve", type=Path, required=True)
    p.add_argument("--expected", type=Path, required=True)
    p.add_argument("--tolerance", type=float, default=1e-12)
    args = p.parse_args()

    expected = json.loads(args.expected.read_text(encoding="utf-8"))
    discovery = replay(args.discovery)
    reserve = replay(args.reserve)
    check_cohort("discovery", discovery, expected, tol=args.tolerance)
    check_cohort("reserve", reserve, expected, tol=args.tolerance)

    print(
        json.dumps(
            {
                "status": "PASS",
                "tolerance": args.tolerance,
                "discovery": discovery,
                "reserve": reserve,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
