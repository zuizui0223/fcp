#!/usr/bin/env python3
"""Metadata-only recovery of the frozen five-block environmental Holm finalizer.

The original finalizer failed after loading the five raw block results because
terrain_structure uses the legacy reporting key
`terrain_edge_occurrence_coverage_fraction` rather than the generic
`block_edge_occurrence_coverage_fraction`. This recovery changes only that
reporting-field lookup. Family, raw p-values, Holm algorithm, effect-direction
gate, alpha, and all biological inputs are unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "docs/supporting/global_rgfca_expanded_environmental_process_panel_contract_v1.json"
ORDER = ROOT / "docs/supporting/global_rgfca_expanded_environmental_panel_execution_v1.json"
EXPECTED = [
    "terrain_structure",
    "thermal_regime",
    "water_balance",
    "atmospheric_energy_dryness",
    "edaphic_regime",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def holm(raw: dict[str, float]) -> dict[str, float]:
    if list(raw) != EXPECTED:
        raise RuntimeError("Holm family must contain exactly the five frozen blocks in frozen order")
    vals = {k: float(v) for k, v in raw.items()}
    if any((not np.isfinite(v)) or v < 0.0 or v > 1.0 for v in vals.values()):
        raise RuntimeError("all raw p-values must be finite in [0,1]")
    ordered = sorted(vals.items(), key=lambda kv: kv[1])
    m = len(ordered)
    running = 0.0
    adjusted: dict[str, float] = {}
    for i, (name, p) in enumerate(ordered):
        running = max(running, min(1.0, (m - i) * p))
        adjusted[name] = float(min(1.0, running))
    return {name: adjusted[name] for name in EXPECTED}


def coverage_fraction(name: str, r: dict) -> float:
    if "block_edge_occurrence_coverage_fraction" in r:
        return float(r["block_edge_occurrence_coverage_fraction"])
    if name == "terrain_structure" and "terrain_edge_occurrence_coverage_fraction" in r:
        return float(r["terrain_edge_occurrence_coverage_fraction"])
    raise RuntimeError(f"missing frozen coverage field for {name}")


def main() -> int:
    ap = argparse.ArgumentParser()
    for name in EXPECTED:
        ap.add_argument(f"--{name.replace('_','-')}", dest=name, type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    panel = json.loads(PANEL.read_text())
    order = json.loads(ORDER.read_text())
    if panel.get("status") != "postoutcome_secondary_panel_frozen_before_any_expanded_environmental_panel_result_or_new_source_payload_is_opened":
        raise RuntimeError("panel contract drift")
    if order.get("status") != "technical_execution_order_frozen_before_any_expanded_panel_block_colour_alignment_result":
        raise RuntimeError("panel execution-order drift")
    if order.get("execution_order") != EXPECTED:
        raise RuntimeError("frozen block order drift")
    if panel["primary_edge_test"].get("multiplicity") != "Holm across all evaluable process blocks in this five-block panel":
        raise RuntimeError("multiplicity contract drift")

    results: dict[str, dict] = {}
    paths: dict[str, Path] = {name: getattr(args, name) for name in EXPECTED}
    for name in EXPECTED:
        p = paths[name]
        r = json.loads(p.read_text())
        if r.get("block") != name:
            raise RuntimeError(f"block identity mismatch for {name}")
        if r.get("null_permutations") != 999 or r.get("null_indices_exact_0_998") is not True:
            raise RuntimeError(f"incomplete null execution for {name}")
        if r.get("final_panel_support_decision_available") is not False:
            raise RuntimeError(f"raw block {name} opened final support early")
        if r.get("p_holm") is not None or r.get("panel_supported") is not None:
            raise RuntimeError(f"raw block {name} already contains a final family decision")
        if r.get("individual_variable_decomposition_open") is not False:
            raise RuntimeError(f"raw block {name} opened variable decomposition early")
        if int(r.get("n_evaluable_species", 0)) < int(panel["primary_edge_test"]["minimum_evaluable_species"]):
            raise RuntimeError(f"frozen block {name} is not evaluable")
        coverage_fraction(name, r)
        results[name] = r

    raw = {name: float(results[name]["raw_p_upper"]) for name in EXPECTED}
    adj = holm(raw)
    alpha = float(panel["primary_edge_test"]["alpha"])
    block_rows = []
    supported_blocks: list[str] = []
    for name in EXPECTED:
        r = results[name]
        mean_rho = float(r["observed_mean_species_partial_rho"])
        supported = bool(mean_rho > 0.0 and adj[name] < alpha)
        if supported:
            supported_blocks.append(name)
        block_rows.append({
            "block": name,
            "observed_mean_species_partial_rho": mean_rho,
            "observed_median_species_partial_rho": float(r["observed_median_species_partial_rho"]),
            "observed_positive_species_fraction": float(r["observed_positive_species_fraction"]),
            "n_evaluable_species": int(r["n_evaluable_species"]),
            "block_edge_occurrence_coverage_fraction": coverage_fraction(name, r),
            "raw_p_upper": raw[name],
            "p_holm": adj[name],
            "supported_after_five_block_holm": supported,
            "individual_variable_decomposition_open": supported,
            "raw_result_sha256": sha256_file(paths[name]),
        })

    payload = {
        "protocol": "global-rgfca-expanded-environmental-process-panel-final-v1",
        "status": "complete_fixed_five_block_environmental_panel_holm_technical_recovery",
        "technical_recovery": {
            "original_finalization_run_id": 34021682120,
            "original_failure": "KeyError: block_edge_occurrence_coverage_fraction on terrain_structure reporting metadata",
            "inferential_change": False,
            "only_change": "terrain coverage reporting key fallback to terrain_edge_occurrence_coverage_fraction",
        },
        "inferential_role": panel["inferential_role"],
        "family_size": 5,
        "family": EXPECTED,
        "multiplicity": "Holm across exactly the five frozen parent process blocks",
        "alpha": alpha,
        "blocks": block_rows,
        "supported_blocks": supported_blocks,
        "any_environmental_process_block_supported": bool(supported_blocks),
        "individual_variable_decomposition_open_for": supported_blocks,
        "unsupported_parent_blocks_may_not_be_decomposed_for_significance": True,
        "cannot_rescue_or_reclassify_primary_G1": True,
        "cannot_reclassify_species_disjoint_commonness": True,
        "causal_or_local_adaptation_language_allowed": False,
        "lineage": {
            "panel_contract_sha256": sha256_file(PANEL),
            "panel_execution_sha256": sha256_file(ORDER),
            "raw_block_result_sha256": {name: sha256_file(paths[name]) for name in EXPECTED},
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
