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

from fcp_pipeline.p500_prospective_execution_gate import (  # noqa: E402
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    MIN_H2_VECTOR_SPECIES,
    MIN_MEASUREMENT_EVALUABLE_SPECIES,
    validate_stage_transition,
)
from run_polymorphism_delta_geometry_validation_20260912 import (  # noqa: E402
    BIO,
    load_classifiable,
    nclass40_count,
    species_delta_vectors,
)
from run_polymorphism_delta_geometry_structured_null_20260912 import (  # noqa: E402
    prepare_selected,
    permute_within_coarse_morph,
    delta_units_from_compositions,
)

PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260915.md"
N_NULL = 999
MIN_VECTOR_SPECIES = MIN_H2_VECTOR_SPECIES
THRESHOLDS = {"primary_0_10": 0.10, "strict_0_20": 0.20}
SEEDS = {"primary_0_10": 20260915, "strict_0_20": 20261015}
EPS = 1e-12


def white_contrast() -> np.ndarray:
    q = np.array([1.0] + [-1.0 / 8.0] * 8, dtype=float)
    q /= np.linalg.norm(q)
    if not np.isclose(float(q.sum()), 0.0, atol=1e-12):
        raise RuntimeError("fixed q_white is not in the zero-sum subspace")
    return q


def units_from_vectors(table: pd.DataFrame) -> np.ndarray:
    x = table[[f"delta_{c}" for c in BIO]].to_numpy(float)
    norms = np.linalg.norm(x, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= EPS):
        raise RuntimeError("prospective H2 table contains invalid Delta vectors")
    return x / norms[:, None]


def statistic(units: np.ndarray, q: np.ndarray) -> float:
    return float(np.mean((units @ q) ** 2))


def mc_upper(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def qsummary(x: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)),
        "q025": float(np.quantile(x, 0.025)),
        "q50": float(np.quantile(x, 0.50)),
        "q975": float(np.quantile(x, 0.975)),
    }


def run_threshold(
    work: pd.DataFrame,
    table: pd.DataFrame,
    q: np.ndarray,
    seed: int,
) -> tuple[dict, np.ndarray | None]:
    n_species = int(len(table))
    if n_species < MIN_VECTOR_SPECIES:
        return {
            "evaluable": False,
            "species": n_species,
            "minimum_vector_species": MIN_VECTOR_SPECIES,
            "reason": "insufficient_frozen_h2_vector_support",
            "structured_null_upper_p": None,
            "pass": None,
        }, None

    observed_units = units_from_vectors(table)
    observed = statistic(observed_units, q)
    prepared = prepare_selected(work, table["species"].astype(str).tolist())
    rng = np.random.default_rng(seed)
    null = np.empty(N_NULL, dtype=float)
    for b in range(N_NULL):
        permuted = permute_within_coarse_morph(prepared, rng)
        units = delta_units_from_compositions(permuted, prepared["species_indices"])
        null[b] = statistic(units, q)

    p = mc_upper(observed, null)
    summary = qsummary(null)
    return {
        "evaluable": True,
        "species": n_species,
        "minimum_vector_species": MIN_VECTOR_SPECIES,
        "observed_W": observed,
        "structured_null_replicates": N_NULL,
        "structured_null_upper_p": p,
        "structured_null_summary": summary,
        "observed_minus_null_median": float(observed - summary["q50"]),
        "observed_to_null_median_ratio": float(observed / summary["q50"]),
        "isotropic_8d_expectation": 0.125,
        "pass": bool(p < 0.05),
    }, null


def measurement_support_from_receipt(measurement: dict) -> dict:
    if measurement.get("schema") != "p500_prospective_measurement_result_v1":
        raise RuntimeError("unexpected prospective measurement receipt schema")
    if measurement.get("status") != "complete_p500_location_blind_measurement_and_support_gate":
        raise RuntimeError("prospective measurement receipt is incomplete")
    if measurement.get("stage") != "SUPPORT_GATE_COMPLETE":
        raise RuntimeError("prospective measurement receipt did not reach SUPPORT_GATE_COMPLETE")
    if measurement.get("species") != EXPECTED_SPECIES or measurement.get("rows") != EXPECTED_ROWS:
        raise RuntimeError("prospective measurement receipt denominator mismatch")
    if measurement.get("H2_opened") is not False:
        raise RuntimeError("measurement receipt says H2 was already opened")
    if measurement.get("replacement_rows", 0) != 0 or measurement.get("replacement_species", 0) != 0:
        raise RuntimeError("measurement receipt contains forbidden replacement")
    if measurement.get("persisted_image_pixels") is not False:
        raise RuntimeError("measurement receipt says image pixels persisted")

    n = measurement.get("measurement_evaluable_species")
    if type(n) is not int or not 0 <= n <= EXPECTED_SPECIES:
        raise RuntimeError("invalid measurement-evaluable species count")
    expected = "PASS" if n >= MIN_MEASUREMENT_EVALUABLE_SPECIES else "NOT_EVALUABLE"
    if measurement.get("support_decision") != expected:
        raise RuntimeError("measurement-support decision mismatch")

    support_stage = measurement.get("stage_receipts", {}).get("support")
    if not isinstance(support_stage, dict):
        raise RuntimeError("measurement receipt lacks frozen support-stage record")
    if support_stage.get("stage") != "SUPPORT_GATE_COMPLETE":
        raise RuntimeError("embedded support-stage record is invalid")
    if support_stage.get("measurement_evaluable_species") != n:
        raise RuntimeError("embedded support-stage count differs from measurement receipt")
    if support_stage.get("support_decision") != expected:
        raise RuntimeError("embedded support-stage decision differs from measurement receipt")

    return {
        "pass": expected == "PASS",
        "decision": expected,
        "evaluable_species": n,
        "minimum_evaluable_species": MIN_MEASUREMENT_EVALUABLE_SPECIES,
        "support_stage": support_stage,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--measured-csv", type=Path, required=True)
    p.add_argument("--measurement-result", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def write_result(out: Path, result: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


def repo_display_path(path: Path) -> str:
    """Return a stable repo-relative display path when possible.

    CLI arguments may arrive as either repository-relative or absolute Paths.  This
    helper is presentation-only and must not change any scientific input or output.
    """
    if not path.is_absolute():
        return str(path)
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def main() -> int:
    args = parse_args()
    measurement = json.loads(args.measurement_result.read_text(encoding="utf-8"))
    support = measurement_support_from_receipt(measurement)

    if support["pass"] is not True:
        result = {
            "analysis": "polymorphism_h2_p500_prospective_white_axis",
            "date_jst": "2026-09-15",
            "status": "untouched_prospective_target_not_evaluable",
            "protocol": str(PROTOCOL.relative_to(ROOT)),
            "decision": {
                "verdict": "H2_PROSPECTIVE_NOT_EVALUABLE_MEASUREMENT_SUPPORT",
                "primary_confirmed": False,
            },
            "measurement_support": {k: v for k, v in support.items() if k != "support_stage"},
            "H2_W_computed": False,
            "hard_nonclaims": [
                "measurement-support failure is not evidence against polymorphism",
                "measurement-support failure is not evidence against the white/non-white axis",
            ],
        }
        write_result(args.output_dir, result)
        return 0

    work = load_classifiable(args.measured_csv)
    n40 = nclass40_count(work)
    if n40 != int(support["evaluable_species"]):
        raise RuntimeError(
            f"n_classifiable>=40 count differs between measurement receipt and H2 loader: {n40}"
        )

    q = white_contrast()
    threshold_results: dict[str, dict] = {}
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for label, threshold in THRESHOLDS.items():
        table, vectors, gate_counts = species_delta_vectors(work, threshold)
        if len(table) != len(vectors):
            raise RuntimeError("prospective vector table and matrix have different denominators")
        table.to_csv(args.output_dir / f"{label}_delta_vectors.csv", index=False, lineterminator="\n")
        test, null = run_threshold(work, table, q, SEEDS[label])
        test["threshold"] = threshold
        test["construction_gate_counts"] = {str(k): int(v) for k, v in gate_counts.items()}
        threshold_results[label] = test
        if null is not None:
            pd.DataFrame({"W": null}).to_csv(
                args.output_dir / f"{label}_structured_null.csv", index=False, lineterminator="\n"
            )

    primary = threshold_results["primary_0_10"]
    strict = threshold_results["strict_0_20"]
    if not primary["evaluable"]:
        verdict = "H2_PROSPECTIVE_NOT_EVALUABLE_VECTOR_SUPPORT"
        primary_confirmed = False
        h2_decision = "NOT_EVALUABLE_VECTOR_SUPPORT"
        primary_p = None
    elif primary["pass"]:
        verdict = "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
        primary_confirmed = True
        h2_decision = "CONFIRMED"
        primary_p = float(primary["structured_null_upper_p"])
    else:
        verdict = "H2_PROSPECTIVE_WHITE_AXIS_NOT_CONFIRMED"
        primary_confirmed = False
        h2_decision = "NOT_CONFIRMED"
        primary_p = float(primary["structured_null_upper_p"])

    h2_stage = {
        "stage": "H2_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "primary_h2_vector_species": int(primary["species"]),
        "primary_structured_null_p": primary_p,
        "h2_decision": h2_decision,
    }
    validate_stage_transition(support["support_stage"], h2_stage)

    result = {
        "analysis": "polymorphism_h2_p500_prospective_white_axis",
        "date_jst": "2026-09-15",
        "status": "untouched_prospective_test_of_previously_frozen_axis",
        "stage": "H2_COMPLETE",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "measurement_result": repo_display_path(args.measurement_result),
        "measurement_support": {k: v for k, v in support.items() if k != "support_stage"},
        "fixed_axis_palette_loadings": {BIO[i]: float(q[i]) for i in range(len(BIO))},
        "statistic": "W = mean_i (u_i dot q_white)^2",
        "structured_null": (
            "within prospective cohort and frozen coarse morph, permute normalized nine-colour "
            "rows across already selected H2 species while preserving species x coarse-morph row counts; "
            "refit unchanged label-free Hellinger two-means"
        ),
        "thresholds": threshold_results,
        "stage_receipt": h2_stage,
        "decision": {
            "primary_threshold": 0.10,
            "primary_evaluable": bool(primary["evaluable"]),
            "primary_confirmed": bool(primary_confirmed),
            "strict_sensitivity_evaluable": bool(strict["evaluable"]),
            "strict_sensitivity_pass": strict["pass"],
            "verdict": verdict,
        },
        "H2_W_computed": bool(primary["evaluable"]),
        "hard_nonclaims": [
            "does not estimate global flower-colour polymorphism prevalence",
            "does not establish pigment chemistry or a pigment-loss mechanism",
            "does not identify evolutionary direction of white/non-white transitions",
            "does not establish pollinator, climate, or other adaptive causation",
            "does not support a recurrent non-white hue axis",
        ],
    }
    write_result(args.output_dir, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
