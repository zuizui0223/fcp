#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))

from run_polymorphism_delta_geometry_validation_20260912 import (  # noqa: E402
    BIO,
    THRESHOLDS,
    load_classifiable,
)
from run_polymorphism_delta_geometry_structured_null_20260912 import (  # noqa: E402
    prepare_selected,
    permute_within_coarse_morph,
    delta_units_from_compositions,
)

DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
H2_OUT = ROOT / "results" / "polymorphism_delta_geometry_validation_20260912"
OUT = ROOT / "results" / "polymorphism_white_axis_targeted_test_20260912"
N_NULL = 999
SEED = 20260912
EPS = 1e-12


def white_contrast() -> np.ndarray:
    q = np.array([1.0] + [-1.0 / 8.0] * 8, dtype=float)
    q /= np.linalg.norm(q)
    if not np.isclose(float(q.sum()), 0.0, atol=1e-12):
        raise RuntimeError("white contrast not zero-sum")
    return q


def units_from_table(table: pd.DataFrame) -> np.ndarray:
    x = table[[f"delta_{c}" for c in BIO]].to_numpy(float)
    n = np.linalg.norm(x, axis=1)
    if np.any(~np.isfinite(n)) or np.any(n <= EPS):
        raise RuntimeError("invalid observed Delta vector")
    return x / n[:, None]


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


def run_cohort(
    work: pd.DataFrame,
    vector_table: pd.DataFrame,
    q: np.ndarray,
    seed: int,
) -> tuple[dict, np.ndarray]:
    observed_units = units_from_table(vector_table)
    observed = statistic(observed_units, q)
    prepared = prepare_selected(work, vector_table["species"].astype(str).tolist())
    rng = np.random.default_rng(seed)
    null = np.empty(N_NULL, dtype=float)
    for b in range(N_NULL):
        permuted = permute_within_coarse_morph(prepared, rng)
        units = delta_units_from_compositions(permuted, prepared["species_indices"])
        null[b] = statistic(units, q)

    p = mc_upper(observed, null)
    summary = qsummary(null)
    return {
        "species": int(len(vector_table)),
        "observed_mean_squared_white_axis_alignment": observed,
        "structured_null_upper_p": p,
        "structured_null_summary": summary,
        "observed_minus_null_median": float(observed - summary["q50"]),
        "observed_to_null_median_ratio": float(observed / summary["q50"]),
        "isotropic_8d_expectation": 0.125,
        "pass": bool(p < 0.05),
    }, null


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    discovery = load_classifiable(DISCOVERY)
    reserve = load_classifiable(RESERVE)
    q = white_contrast()

    threshold_results: dict[str, dict] = {}
    for j, label in enumerate(THRESHOLDS):
        dtab = pd.read_csv(H2_OUT / f"{label}_discovery_delta_vectors.csv")
        rtab = pd.read_csv(H2_OUT / f"{label}_reserve_delta_vectors.csv")
        dres, dnull = run_cohort(discovery, dtab, q, SEED + 100 * j)
        rres, rnull = run_cohort(reserve, rtab, q, SEED + 100 * j + 1)
        threshold_results[label] = {
            "discovery": dres,
            "reserve": rres,
            "both_cohorts_pass": bool(dres["pass"] and rres["pass"]),
        }
        pd.DataFrame({"mean_squared_white_axis_alignment": dnull}).to_csv(
            OUT / f"{label}_discovery_structured_null.csv", index=False
        )
        pd.DataFrame({"mean_squared_white_axis_alignment": rnull}).to_csv(
            OUT / f"{label}_reserve_structured_null.csv", index=False
        )

    primary = threshold_results["primary_0_10"]["both_cohorts_pass"]
    strict = threshold_results["strict_0_20"]["both_cohorts_pass"]
    if primary and strict:
        verdict = "WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT"
    elif primary:
        verdict = "WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_ONLY"
    else:
        verdict = "WHITE_AXIS_TARGETED_NOT_SUPPORTED"

    result = {
        "analysis": "polymorphism_white_axis_targeted_test",
        "date_jst": "2026-09-12",
        "status": "targeted_post_audit_decomposition_not_untouched_confirmatory_test",
        "freeze": "docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md",
        "fixed_axis_palette_loadings": {BIO[i]: float(q[i]) for i in range(len(BIO))},
        "statistic": "mean_i (u_i dot q_white)^2",
        "structured_null_replicates": N_NULL,
        "thresholds": threshold_results,
        "decision": {
            "primary_both_cohorts_pass": bool(primary),
            "strict_both_cohorts_pass": bool(strict),
            "verdict": verdict,
        },
        "claim_boundary": (
            "Positive existing-cohort evidence is retrospective because the fixed achromatic/chromatic "
            "candidate was motivated by already opened H2 geometry. It may be frozen as the prospective "
            "axis target for the global high-depth expansion. No pigment pathway, ecological mechanism, "
            "or adaptive interpretation is identified here."
        ),
    }
    (OUT / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# Targeted H2 achromatic/chromatic axis test",
        "",
        f"**Verdict: `{verdict}`**",
        "",
        "This is a targeted post-audit decomposition, not an untouched confirmatory test.",
        "",
    ]
    for label in THRESHOLDS:
        x = threshold_results[label]
        lines += [
            f"## {label}",
            "",
            f"- discovery N: **{x['discovery']['species']}**, W: **{x['discovery']['observed_mean_squared_white_axis_alignment']:.6f}**, structured p: **{x['discovery']['structured_null_upper_p']:.6g}**",
            f"- reserve N: **{x['reserve']['species']}**, W: **{x['reserve']['observed_mean_squared_white_axis_alignment']:.6f}**, structured p: **{x['reserve']['structured_null_upper_p']:.6g}**",
            f"- both cohorts pass: **{x['both_cohorts_pass']}**",
            "",
        ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
