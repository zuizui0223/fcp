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
    FRACTIONS,
    MORPHS,
    OUT as H2_OUT,
    THRESHOLDS,
    deterministic_two_means,
    load_classifiable,
)

DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
OUT = ROOT / "results" / "polymorphism_delta_geometry_structured_null_20260912"
N_NULL = 999
SEED = 20260912
EPS = 1e-12


def prepare_selected(work: pd.DataFrame, selected_species: list[str]) -> dict:
    selected = set(map(str, selected_species))
    sub = work.loc[work["species"].astype(str).isin(selected)].copy().reset_index(drop=True)
    found = set(sub["species"].astype(str).unique())
    missing = sorted(selected - found)
    if missing:
        raise RuntimeError(f"selected species missing from measured table: {missing[:5]}")

    compositions = sub[FRACTIONS].to_numpy(float)
    morphs = sub["morph"].astype(str).to_numpy()
    species = sub["species"].astype(str).to_numpy()
    species_indices = {
        sp: np.flatnonzero(species == sp) for sp in sorted(selected)
    }
    morph_indices = {
        morph: np.flatnonzero(morphs == morph) for morph in MORPHS
    }
    return {
        "compositions": compositions,
        "species_indices": species_indices,
        "morph_indices": morph_indices,
    }


def delta_units_from_compositions(compositions: np.ndarray, species_indices: dict[str, np.ndarray]) -> np.ndarray:
    units: list[np.ndarray] = []
    for sp in sorted(species_indices):
        idx = species_indices[sp]
        p = compositions[idx]
        labels, counts, _ = deterministic_two_means(np.sqrt(p))
        if counts[0] > counts[1]:
            major, minor = 0, 1
        elif counts[1] > counts[0]:
            major, minor = 1, 0
        else:
            c0 = p[labels == 0].mean(axis=0)
            c1 = p[labels == 1].mean(axis=0)
            major, minor = (0, 1) if tuple(c0.tolist()) <= tuple(c1.tolist()) else (1, 0)
        mu1 = p[labels == major].mean(axis=0)
        mu2 = p[labels == minor].mean(axis=0)
        delta = mu2 - mu1
        if not np.isclose(float(delta.sum()), 0.0, atol=1e-10):
            raise RuntimeError(f"null Delta outside zero-sum subspace for {sp}")
        norm = float(np.linalg.norm(delta))
        if not np.isfinite(norm) or norm <= EPS:
            raise RuntimeError(f"degenerate null Delta for {sp}")
        units.append(delta / norm)
    return np.vstack(units)


def permute_within_coarse_morph(prepared: dict, rng: np.random.Generator) -> np.ndarray:
    source = prepared["compositions"]
    permuted = source.copy()
    for morph in MORPHS:
        idx = prepared["morph_indices"][morph]
        if len(idx) > 1:
            permuted[idx] = source[idx[rng.permutation(len(idx))]]
    return permuted


def lambda1(units9: np.ndarray) -> float:
    matrix = units9.T @ units9 / len(units9)
    return float(np.linalg.eigvalsh(matrix)[-1])


def mean_squared_projection(units9: np.ndarray, axis9: np.ndarray) -> float:
    axis9 = np.asarray(axis9, dtype=float)
    axis9 = axis9 / np.linalg.norm(axis9)
    return float(np.mean((units9 @ axis9) ** 2))


def mc_upper(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def qsummary(x: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)),
        "q025": float(np.quantile(x, 0.025)),
        "q50": float(np.quantile(x, 0.5)),
        "q975": float(np.quantile(x, 0.975)),
    }


def run_threshold(
    label: str,
    discovery: pd.DataFrame,
    reserve: pd.DataFrame,
    h2_result: dict,
    seed_offset: int,
) -> tuple[dict, np.ndarray, np.ndarray]:
    d_vectors = pd.read_csv(H2_OUT / f"{label}_discovery_delta_vectors.csv")
    r_vectors = pd.read_csv(H2_OUT / f"{label}_reserve_delta_vectors.csv")
    if d_vectors.empty or r_vectors.empty:
        raise RuntimeError(f"empty H2 selected set for {label}")

    dprep = prepare_selected(discovery, d_vectors["species"].astype(str).tolist())
    rprep = prepare_selected(reserve, r_vectors["species"].astype(str).tolist())

    observed_lambda = float(h2_result["thresholds"][label]["discovery"]["lambda1"])
    observed_transport = float(
        h2_result["thresholds"][label]["transport"][
            "reserve_mean_squared_projection_on_frozen_discovery_axis"
        ]
    )
    axis9 = np.array(
        [
            h2_result["thresholds"][label]["transport"]["discovery_axis_palette_loadings"][c]
            for c in BIO
        ],
        dtype=float,
    )
    if not np.isclose(float(axis9.sum()), 0.0, atol=1e-9):
        raise RuntimeError("frozen discovery axis is not zero-sum")

    rng_d = np.random.default_rng(SEED + seed_offset)
    rng_r = np.random.default_rng(SEED + seed_offset + 1)
    null_lambda = np.empty(N_NULL, dtype=float)
    null_transport = np.empty(N_NULL, dtype=float)

    for b in range(N_NULL):
        dperm = permute_within_coarse_morph(dprep, rng_d)
        dunits = delta_units_from_compositions(dperm, dprep["species_indices"])
        null_lambda[b] = lambda1(dunits)

        rperm = permute_within_coarse_morph(rprep, rng_r)
        runits = delta_units_from_compositions(rperm, rprep["species_indices"])
        null_transport[b] = mean_squared_projection(runits, axis9)

    lambda_p = mc_upper(observed_lambda, null_lambda)
    transport_p = mc_upper(observed_transport, null_transport)
    result = {
        "conditional_species_sets_fixed_from_observed_H2_gate": True,
        "discovery_H2_species": int(len(d_vectors)),
        "reserve_H2_species": int(len(r_vectors)),
        "null_definition": (
            "Within each cohort and each frozen coarse morph state, permute normalized 9-colour "
            "palette rows across the already-selected H2 species while preserving every species x "
            "coarse-morph row count; then refit label-free Hellinger two-means. This preserves the "
            "global coarse-state-to-palette mapping and H1 state composition but destroys "
            "species-specific continuous palette structure."
        ),
        "null_replicates": N_NULL,
        "discovery_axis_concentration": {
            "observed_lambda1": observed_lambda,
            "structured_null_upper_p": lambda_p,
            "null_summary": qsummary(null_lambda),
            "pass": bool(lambda_p < 0.05),
        },
        "reserve_frozen_axis_transport": {
            "observed_mean_squared_projection": observed_transport,
            "structured_null_upper_p": transport_p,
            "null_summary": qsummary(null_transport),
            "pass": bool(transport_p < 0.05),
        },
        "decision": {
            "structured_null_discovery_pass": bool(lambda_p < 0.05),
            "structured_null_transport_pass": bool(transport_p < 0.05),
            "structured_null_H2_pass": bool(lambda_p < 0.05 and transport_p < 0.05),
        },
    }
    return result, null_lambda, null_transport


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    h2_result = json.loads((H2_OUT / "result.json").read_text(encoding="utf-8"))
    discovery = load_classifiable(DISCOVERY)
    reserve = load_classifiable(RESERVE)

    threshold_results: dict[str, dict] = {}
    for j, label in enumerate(THRESHOLDS):
        result, nl, nt = run_threshold(label, discovery, reserve, h2_result, 100 * j)
        threshold_results[label] = result
        pd.DataFrame({"lambda1": nl}).to_csv(OUT / f"{label}_structured_lambda_null.csv", index=False)
        pd.DataFrame({"mean_squared_projection": nt}).to_csv(
            OUT / f"{label}_structured_transport_null.csv", index=False
        )

    primary = threshold_results["primary_0_10"]["decision"]["structured_null_H2_pass"]
    strict = threshold_results["strict_0_20"]["decision"]["structured_null_H2_pass"]
    if primary and strict:
        verdict = "H2_SURVIVES_COARSE_STATE_PRESERVING_STRUCTURED_NULL"
    elif primary:
        verdict = "H2_PRIMARY_ONLY_STRUCTURED_NULL_STRICT_FAILS"
    else:
        verdict = "H2_EXPLAINED_BY_COARSE_STATE_PRESERVING_CONSTRUCTION_NULL"

    out = {
        "analysis": "polymorphism_delta_geometry_coarse_state_preserving_structured_null",
        "date_jst": "2026-09-12",
        "source_H2_result": str(H2_OUT.relative_to(ROOT) / "result.json"),
        "null_replicates": N_NULL,
        "thresholds": threshold_results,
        "decision": {
            "primary_pass": bool(primary),
            "strict_pass": bool(strict),
            "verdict": verdict,
        },
        "interpretation_boundary": (
            "This is a construction-control test conditional on the observed H2-selected species. "
            "Failure means that the recurrent axis can be reproduced after preserving the frozen "
            "coarse-state composition and its global palette mapping, so it cannot be presented as "
            "species-specific continuous geometry beyond state construction. Passing does not by "
            "itself establish an ecological mechanism."
        ),
    }
    (OUT / "result.json").write_text(
        json.dumps(out, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# H2 coarse-state-preserving structured null",
        "",
        f"**Verdict: `{verdict}`**",
        "",
    ]
    for label in THRESHOLDS:
        x = threshold_results[label]
        lines += [
            f"## {label}",
            "",
            f"- discovery observed lambda1: **{x['discovery_axis_concentration']['observed_lambda1']:.6f}**",
            f"- structured-null p: **{x['discovery_axis_concentration']['structured_null_upper_p']:.6g}**",
            f"- reserve frozen-axis projection: **{x['reserve_frozen_axis_transport']['observed_mean_squared_projection']:.6f}**",
            f"- structured-null transport p: **{x['reserve_frozen_axis_transport']['structured_null_upper_p']:.6g}**",
            "",
        ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
