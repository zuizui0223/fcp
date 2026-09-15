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
FROZEN_RESULT = ROOT / "results" / "polymorphism_white_axis_targeted_test_20260912" / "result.json"
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
            "Qualify the prospective P500 H2 executor by replaying the frozen 2026-09-12 white-axis "
            "result on already-opened discovery/reserve cohorts. Delta vectors and null draws are rebuilt "
            "in memory because the historical CSV intermediates were not committed."
        )
    )
    p.add_argument("--output", type=Path, default=None)
    return p.parse_args()


def _assert_close(name: str, got: float, expected: float) -> None:
    if not np.isclose(float(got), float(expected), rtol=0.0, atol=ATOL):
        raise RuntimeError(f"{name} mismatch: got={got!r} expected={expected!r}")


def _check_frozen_summary(*, label: str, cohort: str, got: dict, expected: dict) -> None:
    if int(got["species"]) != int(expected["species"]):
        raise RuntimeError(
            f"{label}/{cohort} species mismatch: {got['species']} != {expected['species']}"
        )
    _assert_close(
        f"{label}/{cohort} W",
        got["observed_mean_squared_white_axis_alignment"],
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


def _check_case(
    *,
    label: str,
    cohort: str,
    work: pd.DataFrame,
    frozen: dict,
) -> dict:
    q_legacy = legacy.white_contrast()
    q_prospective = prospective.white_contrast()
    if not np.array_equal(q_legacy, q_prospective):
        raise RuntimeError("prospective q_white differs from frozen legacy q_white")

    threshold = prospective.THRESHOLDS[label]
    vectors, vector_matrix, gate_counts = prospective.species_delta_vectors(work, threshold)
    if len(vectors) != len(vector_matrix):
        raise RuntimeError(f"{label}/{cohort} rebuilt vector table/matrix denominator mismatch")

    seed = LEGACY_SEEDS[(label, cohort)]
    legacy_got, legacy_null = legacy.run_cohort(work, vectors, q_legacy, seed)
    prospective_got, prospective_null = prospective.run_threshold(
        work, vectors, q_prospective, seed
    )
    if prospective_got.get("evaluable") is not True or prospective_null is None:
        raise RuntimeError(f"legacy replay unexpectedly not evaluable for {label}/{cohort}")

    # First anchor the rebuilt historical computation to the committed frozen result.json.
    expected = frozen["thresholds"][label][cohort]
    _check_frozen_summary(label=label, cohort=cohort, got=legacy_got, expected=expected)

    # Then require the prospective executor to be numerically identical to that legacy computation.
    if int(prospective_got["species"]) != int(legacy_got["species"]):
        raise RuntimeError(f"{label}/{cohort} prospective species denominator differs from legacy")
    _assert_close(
        f"{label}/{cohort} prospective-vs-legacy W",
        prospective_got["observed_W"],
        legacy_got["observed_mean_squared_white_axis_alignment"],
    )
    _assert_close(
        f"{label}/{cohort} prospective-vs-legacy p",
        prospective_got["structured_null_upper_p"],
        legacy_got["structured_null_upper_p"],
    )
    _assert_close(
        f"{label}/{cohort} prospective-vs-legacy observed-minus-null-median",
        prospective_got["observed_minus_null_median"],
        legacy_got["observed_minus_null_median"],
    )
    _assert_close(
        f"{label}/{cohort} prospective-vs-legacy observed-to-null-median-ratio",
        prospective_got["observed_to_null_median_ratio"],
        legacy_got["observed_to_null_median_ratio"],
    )
    for key, value in legacy_got["structured_null_summary"].items():
        _assert_close(
            f"{label}/{cohort} prospective-vs-legacy null-summary/{key}",
            prospective_got["structured_null_summary"][key],
            value,
        )
    if bool(prospective_got["pass"]) != bool(legacy_got["pass"]):
        raise RuntimeError(f"{label}/{cohort} prospective pass/fail differs from legacy")

    if len(legacy_null) != len(prospective_null):
        raise RuntimeError(
            f"{label}/{cohort} null length mismatch: {len(prospective_null)} != {len(legacy_null)}"
        )
    if not np.array_equal(prospective_null, legacy_null):
        max_abs = float(np.max(np.abs(prospective_null - legacy_null)))
        raise RuntimeError(
            f"{label}/{cohort} prospective structured-null vector differs from legacy; "
            f"max_abs_diff={max_abs:.17g}"
        )

    return {
        "label": label,
        "cohort": cohort,
        "threshold": float(threshold),
        "seed": seed,
        "species": int(prospective_got["species"]),
        "observed_W": float(prospective_got["observed_W"]),
        "structured_null_upper_p": float(prospective_got["structured_null_upper_p"]),
        "null_replicates": int(len(prospective_null)),
        "construction_gate_counts": {str(k): int(v) for k, v in gate_counts.items()},
        "exact_legacy_prospective_null_vector_match": True,
        "frozen_result_json_numeric_match": True,
        "pass": bool(prospective_got["pass"]),
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
            cases.append(
                _check_case(
                    label=label,
                    cohort=cohort,
                    work=work_by_cohort[cohort],
                    frozen=frozen,
                )
            )

    result = {
        "analysis": "polymorphism_h2_p500_executor_legacy_replay_qualification",
        "date_jst": "2026-09-15",
        "status": "executor_qualified_by_frozen_summary_and_exact_legacy_kernel_replay",
        "qualification_pass": True,
        "prospective_executor": "scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py",
        "legacy_executor": "scripts/analysis/run_polymorphism_white_axis_targeted_test_20260912.py",
        "frozen_reference": "results/polymorphism_white_axis_targeted_test_20260912/result.json",
        "historical_intermediate_csvs_committed": False,
        "replay_method": (
            "rebuild Delta vectors from the same already-opened discovery/reserve source tables using the "
            "legacy construction; require the rebuilt legacy summary to match frozen result.json; then require "
            "prospective and legacy executors to return exactly identical 999-draw structured-null vectors"
        ),
        "cases": cases,
        "cases_passed": len(cases),
        "cases_required": 4,
        "candidate_pixels_opened": False,
        "uses_existing_opened_cohorts_only": True,
        "claim_boundary": (
            "This qualification establishes implementation continuity of q_white, W, and the structured-null "
            "kernel. The historical vector/null CSV intermediates were not committed, so frozen-result continuity "
            "is established through exact reconstruction of the legacy computation plus committed result.json, "
            "not by comparison to unavailable historical CSV bytes. It does not open P500 pixels and is not "
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
