#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))

import run_polymorphism_white_axis_targeted_test_20260912 as legacy  # noqa: E402
import run_polymorphism_h2_p500_prospective_white_axis_20260915 as prospective  # noqa: E402

DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
VECTOR_DIR = ROOT / "results" / "polymorphism_delta_geometry_validation_20260912"
FROZEN_DIR = ROOT / "results" / "polymorphism_white_axis_targeted_test_20260912"
FROZEN_RESULT = FROZEN_DIR / "result.json"
LABELS = ("primary_0_10", "strict_0_20")
LEGACY_SEEDS = {
    ("primary_0_10", "discovery"): 20260912,
    ("primary_0_10", "reserve"): 20260913,
    ("strict_0_20", "discovery"): 20261012,
    ("strict_0_20", "reserve"): 20261013,
}
ATOL = 1e-15


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Qualify the prospective P500 H2 executor by exact replay of the frozen 2026-09-12 "
            "white-axis result using already-opened discovery/reserve cohorts only."
        )
    )
    p.add_argument("--output", type=Path, default=None)
    return p.parse_args()


def _assert_close(name: str, got: float, expected: float) -> None:
    if not np.isclose(float(got), float(expected), rtol=0.0, atol=ATOL):
        raise RuntimeError(f"{name} mismatch: got={got!r} expected={expected!r}")


def _check_case(
    *,
    label: str,
    cohort: str,
    work: pd.DataFrame,
    vectors: pd.DataFrame,
    frozen: dict,
) -> dict:
    q_legacy = legacy.white_contrast()
    q_prospective = prospective.white_contrast()
    if not np.array_equal(q_legacy, q_prospective):
        raise RuntimeError("prospective q_white differs from frozen legacy q_white")

    seed = LEGACY_SEEDS[(label, cohort)]
    got, null = prospective.run_threshold(work, vectors, q_prospective, seed)
    if got.get("evaluable") is not True or null is None:
        raise RuntimeError(f"legacy replay unexpectedly not evaluable for {label}/{cohort}")

    expected = frozen["thresholds"][label][cohort]
    if int(got["species"]) != int(expected["species"]):
        raise RuntimeError(
            f"{label}/{cohort} species mismatch: {got['species']} != {expected['species']}"
        )
    _assert_close(
        f"{label}/{cohort} W",
        got["observed_W"],
        expected["observed_mean_squared_white_axis_alignment"],
    )
    _assert_close(
        f"{label}/{cohort} p",
        got["structured_null_upper_p"],
        expected["structured_null_upper_p"],
    )
    _assert_close(
        f"{label}/{cohort} observed-minus-null-median",
        got["observed_minus_null_median"],
        expected["observed_minus_null_median"],
    )
    _assert_close(
        f"{label}/{cohort} observed-to-null-median-ratio",
        got["observed_to_null_median_ratio"],
        expected["observed_to_null_median_ratio"],
    )
    for key, value in expected["structured_null_summary"].items():
        _assert_close(
            f"{label}/{cohort} null-summary/{key}",
            got["structured_null_summary"][key],
            value,
        )
    if bool(got["pass"]) != bool(expected["pass"]):
        raise RuntimeError(f"{label}/{cohort} pass/fail differs from frozen result")

    frozen_null_path = FROZEN_DIR / f"{label}_{cohort}_structured_null.csv"
    frozen_null = pd.read_csv(frozen_null_path)["mean_squared_white_axis_alignment"].to_numpy(float)
    if len(frozen_null) != len(null):
        raise RuntimeError(
            f"{label}/{cohort} null length mismatch: {len(null)} != {len(frozen_null)}"
        )
    if not np.array_equal(null, frozen_null):
        max_abs = float(np.max(np.abs(null - frozen_null)))
        raise RuntimeError(
            f"{label}/{cohort} structured-null replay is not exact; max_abs_diff={max_abs:.17g}"
        )

    return {
        "label": label,
        "cohort": cohort,
        "seed": seed,
        "species": int(got["species"]),
        "observed_W": float(got["observed_W"]),
        "structured_null_upper_p": float(got["structured_null_upper_p"]),
        "null_replicates": int(len(null)),
        "exact_null_vector_match": True,
        "frozen_numeric_summary_match": True,
        "pass": bool(got["pass"]),
    }


def main() -> int:
    args = parse_args()
    frozen = json.loads(FROZEN_RESULT.read_text(encoding="utf-8"))

    discovery = prospective.load_classifiable(DISCOVERY)
    reserve = prospective.load_classifiable(RESERVE)
    work_by_cohort = {"discovery": discovery, "reserve": reserve}

    cases: list[dict] = []
    for label in LABELS:
        for cohort in ("discovery", "reserve"):
            vectors = pd.read_csv(VECTOR_DIR / f"{label}_{cohort}_delta_vectors.csv")
            cases.append(
                _check_case(
                    label=label,
                    cohort=cohort,
                    work=work_by_cohort[cohort],
                    vectors=vectors,
                    frozen=frozen,
                )
            )

    result = {
        "analysis": "polymorphism_h2_p500_executor_legacy_replay_qualification",
        "date_jst": "2026-09-15",
        "status": "executor_qualified_by_exact_frozen_legacy_replay",
        "qualification_pass": True,
        "prospective_executor": "scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py",
        "frozen_reference": "results/polymorphism_white_axis_targeted_test_20260912/result.json",
        "cases": cases,
        "cases_passed": len(cases),
        "cases_required": 4,
        "candidate_pixels_opened": False,
        "uses_existing_opened_cohorts_only": True,
        "claim_boundary": (
            "This qualification establishes implementation continuity of q_white, W, and the structured-null "
            "kernel by exact replay of already-opened legacy cohorts. It does not open P500 pixels and is not "
            "prospective H2 evidence."
        ),
    }

    text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
