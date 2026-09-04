#!/usr/bin/env python3
"""Recover the predeclared regional H1 sensitivity from frozen H1 null maps.

This is a post-outcome *execution* of an analysis scope that was prospectively
registered in random_photo_first_boundary_persistence_contract_v1.json before the
fresh candidate pool was opened: report global, continent-within and biome-within.

No H1 sampling, transition rule, morph classification, seed, permutation or null
is rerun or changed here. The script only restricts the already completed observed
edge-persistence vector and the exact matched 999 null persistence maps to fixed
edge subsets. It first replays the global concentration and p-value exactly; any
mismatch fails closed.

The frozen environment frame contains biome and realm labels but no exact
continent field. Therefore continent-within is reported as not evaluable rather
than being silently replaced by realm or by an outcome-postdated geography source.
Realm is retained only as a separately labelled diagnostic because it was already
frozen for the H2 sensitivity analysis.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"
H1_CONTRACT = ROOT / "docs/supporting/random_photo_first_boundary_persistence_contract_v1.json"
H2_CONTRACT = ROOT / "docs/supporting/random_photo_first_h2_climate_contract_v1.json"
H1_EDGES = ROOT / "data/derived/random_photo_first_h1_edges_v1.csv"
H1_NULL_MAPS = ROOT / "data/derived/random_photo_first_h1_null_maps_v1.npz"
H2_CLIMATE_EDGES = ROOT / "data/derived/random_photo_first_h2_climate_edges_250km_v1.csv"
OUT_JSON = ROOT / "docs/supporting/random_photo_first_h1_region_audit_v1.json"
OUT_CSV = ROOT / "data/derived/random_photo_first_h1_region_audit_v1.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def concentration_from_persistence(
    persistence: np.ndarray,
    opportunities: np.ndarray,
    subset_mask: np.ndarray,
) -> tuple[float, float, int]:
    p = np.asarray(persistence, dtype=float)
    w = np.asarray(opportunities, dtype=float)
    mask = np.asarray(subset_mask, dtype=bool)
    if p.shape != w.shape or p.shape != mask.shape:
        raise ValueError("persistence, opportunities and subset mask must align")
    supported = mask & (w > 0) & np.isfinite(p)
    n_supported = int(np.count_nonzero(supported))
    if n_supported == 0:
        raise ValueError("regional subset has no supported H1 edges")
    total_w = float(w[supported].sum())
    if total_w <= 0.0:
        raise ValueError("regional subset has zero total opportunity")
    rate = float(np.sum(w[supported] * p[supported]) / total_w)
    concentration = float(np.sum(w[supported] * (p[supported] - rate) ** 2) / total_w)
    return concentration, rate, n_supported


def regional_test(
    observed_persistence: np.ndarray,
    null_persistence: np.ndarray,
    opportunities: np.ndarray,
    subset_mask: np.ndarray,
) -> dict[str, Any]:
    observed_concentration, transition_rate, n_supported = concentration_from_persistence(
        observed_persistence, opportunities, subset_mask
    )
    null = np.empty(null_persistence.shape[0], dtype=float)
    for i in range(null_persistence.shape[0]):
        null[i], _, null_supported = concentration_from_persistence(
            null_persistence[i], opportunities, subset_mask
        )
        if null_supported != n_supported:
            raise RuntimeError("regional supported-edge denominator drifted across matched null maps")
    if not np.isfinite(null).all():
        raise ValueError("regional null concentration contains non-finite values")
    p_upper = float((1 + np.count_nonzero(null >= observed_concentration)) / (len(null) + 1))
    return {
        "status": "complete_evaluable",
        "observed_concentration": observed_concentration,
        "realized_transition_rate": transition_rate,
        "p_upper": p_upper,
        "null_mean": float(np.mean(null)),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "supported_edges": n_supported,
        "permutations": int(len(null)),
    }


def main() -> int:
    h1 = load_json(H1_RESULT)
    h1_contract = load_json(H1_CONTRACT)
    h2_contract = load_json(H2_CONTRACT)
    if h1.get("status") != "complete_h1_evaluable":
        raise RuntimeError("completed evaluable H1 is required")
    if h1_contract.get("predeclared_sensitivities", {}).get(
        "report_global_continent_within_and_biome_within"
    ) is not True:
        raise RuntimeError("regional H1 sensitivity was not prospectively registered")
    if h2_contract.get("status") != "prospectively_frozen_before_h1_biological_outcome_and_before_any_h2_join":
        raise RuntimeError("biome/realm source was not frozen before H1 outcome")

    expected_hashes = h1["file_sha256"]
    if sha256(H1_EDGES) != expected_hashes["edge_table"]:
        raise RuntimeError("H1 edge table hash drifted")
    if sha256(H1_NULL_MAPS) != expected_hashes["null_persistence_maps"]:
        raise RuntimeError("H1 null-map hash drifted")

    edges = pd.read_csv(H1_EDGES)
    climate_edges = pd.read_csv(H2_CLIMATE_EDGES)
    required_h1 = {"edge_id", "opportunities", "persistence"}
    required_region = {"edge_id", "within_biome", "within_realm"}
    if not required_h1.issubset(edges.columns):
        raise RuntimeError("H1 edge table lacks required fields")
    if not required_region.issubset(climate_edges.columns):
        raise RuntimeError("frozen region edge table lacks biome/realm fields")

    with np.load(H1_NULL_MAPS, allow_pickle=False) as archive:
        edge_ids = archive["edge_ids"].astype(str)
        null_persistence = archive["null_persistence"].astype(float)
        null_concentrations = archive["null_concentrations"].astype(float)
    observed_ids = edges["edge_id"].astype(str).to_numpy()
    if not np.array_equal(edge_ids, observed_ids):
        raise RuntimeError("H1 null edge order drifted from observed edge table")
    if null_persistence.shape != (999, len(edges)) or len(null_concentrations) != 999:
        raise RuntimeError("expected exact matched 999 H1 null maps")

    joined = edges[["edge_id", "opportunities", "persistence"]].merge(
        climate_edges[["edge_id", "within_biome", "within_realm"]],
        on="edge_id",
        how="left",
        validate="one_to_one",
        sort=False,
    )
    if len(joined) != len(edges) or joined[["within_biome", "within_realm"]].isna().any().any():
        raise RuntimeError("regional labels do not cover the frozen H1 edge graph exactly")
    if not np.array_equal(joined["edge_id"].astype(str).to_numpy(), observed_ids):
        raise RuntimeError("regional join changed H1 edge order")

    observed = joined["persistence"].to_numpy(dtype=float)
    opportunities = joined["opportunities"].to_numpy(dtype=float)
    global_mask = np.ones(len(joined), dtype=bool)
    global_replay = regional_test(observed, null_persistence, opportunities, global_mask)

    # Fail closed unless the subset implementation exactly recovers the completed H1.
    if not np.isclose(
        global_replay["observed_concentration"],
        float(h1["primary"]["observed_concentration"]),
        rtol=0.0,
        atol=1e-12,
    ):
        raise RuntimeError("regional statistic does not exactly replay H1 global concentration")
    if not np.isclose(global_replay["p_upper"], float(h1["primary"]["p_upper"]), rtol=0.0, atol=1e-12):
        raise RuntimeError("regional null does not exactly replay H1 global p-value")
    if not np.allclose(null_concentrations, [
        concentration_from_persistence(row, opportunities, global_mask)[0]
        for row in null_persistence
    ], rtol=0.0, atol=1e-12):
        raise RuntimeError("stored H1 null concentrations do not match replayed null maps")

    within_biome = regional_test(
        observed,
        null_persistence,
        opportunities,
        joined["within_biome"].astype(bool).to_numpy(),
    )
    within_realm = regional_test(
        observed,
        null_persistence,
        opportunities,
        joined["within_realm"].astype(bool).to_numpy(),
    )
    alpha = float(h1["alpha"])
    within_biome["support"] = bool(within_biome["p_upper"] < alpha)
    within_realm["support"] = bool(within_realm["p_upper"] < alpha)

    continent = {
        "status": "not_evaluable_missing_prospectively_frozen_continent_labels",
        "reason": (
            "The prospectively frozen region frame contains deterministic modal biome and realm labels, "
            "but no exact continent field. Realm is not silently substituted for continent, and no "
            "outcome-postdated geography source is introduced."
        ),
        "support": False,
    }

    summary = pd.DataFrame(
        [
            {"subset": "global_replay", **global_replay, "support": bool(global_replay["p_upper"] < alpha)},
            {"subset": "within_biome_predeclared", **within_biome},
            {"subset": "within_realm_existing_diagnostic", **within_realm},
            {"subset": "within_continent_predeclared", **continent},
        ]
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_CSV, index=False, lineterminator="\n")

    payload = {
        "protocol": "random-photo-first-h1-region-audit-v1",
        "status": "complete_predeclared_regional_sensitivity_recovery",
        "role": "sensitivity_only_primary_global_h1_decision_unchanged",
        "primary_h1_decision": h1["decision"],
        "primary_h1_p_upper": float(h1["primary"]["p_upper"]),
        "alpha": alpha,
        "global_replay": global_replay,
        "within_biome_predeclared": within_biome,
        "within_continent_predeclared": continent,
        "within_realm_existing_diagnostic": within_realm,
        "interpretation_rule": (
            "Regional sensitivity cannot replace or rescue the prospectively primary global H1 result. "
            "A regional p-value below alpha is reported only as scale/subset-specific sensitivity evidence."
        ),
        "lineage": {
            "h1_edge_table_sha256": sha256(H1_EDGES),
            "h1_null_maps_sha256": sha256(H1_NULL_MAPS),
            "h1_permutations": int(h1["primary"]["permutations"]),
            "h1_sampling_seed": int(h1["primary"]["sampling_seed"]),
            "h1_permutation_seed": int(h1["primary"]["permutation_seed"]),
            "biome_realm_source_frozen_before_h1_outcome": True,
        },
        "files": {"summary": str(OUT_CSV.relative_to(ROOT))},
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
