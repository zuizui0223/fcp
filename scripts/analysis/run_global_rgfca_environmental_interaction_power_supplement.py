#!/usr/bin/env python3
"""Exact joint-null design-power supplement for the frozen 10-interaction family."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_environmental_interaction_power_supplement_contract_v1.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def holm_adjust_matrix(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = p.shape[1]
    order = np.argsort(p, axis=1)
    sorted_p = np.take_along_axis(p, order, axis=1)
    adjusted_sorted = np.maximum.accumulate(sorted_p * np.arange(m, 0, -1, dtype=float)[None, :], axis=1)
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)
    out = np.empty_like(adjusted_sorted)
    rows = np.arange(p.shape[0])[:, None]
    out[rows, order] = adjusted_sorted
    return out


def empirical_two_sided_p(stat: np.ndarray, sorted_abs_null: list[np.ndarray]) -> np.ndarray:
    stat = np.asarray(stat, dtype=float)
    out = np.empty_like(stat)
    n = len(sorted_abs_null[0])
    for j, ref in enumerate(sorted_abs_null):
        count_ge = n - np.searchsorted(ref, np.abs(stat[:, j]), side="left")
        out[:, j] = (1.0 + count_ge) / float(n + 1)
    return out


def main() -> int:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_before_exact_interaction_family_power_result":
        raise RuntimeError("interaction power supplement contract drift")
    source = contract["source"]
    files = sorted(args.artifact_root.rglob(str(source["null_filename"])))
    if len(files) != 1:
        raise RuntimeError(f"expected one exact null file, found {files}")
    path = files[0]
    digest = sha256_file(path)
    if digest != str(source["null_sha256"]):
        raise RuntimeError(f"null SHA mismatch: {digest}")
    df = pd.read_csv(path)
    if "null_index" not in df.columns:
        raise RuntimeError("exact null file lacks null_index")
    if len(df) != int(source["expected_null_rows"]):
        raise RuntimeError(f"null row count drift: {len(df)}")
    if set(df["null_index"].astype(int)) != set(range(999)):
        raise RuntimeError("null indices are not exactly 0..998")
    features = [c for c in df.columns if c != "null_index"]
    if len(features) != int(source["expected_family_size"]):
        raise RuntimeError(f"interaction family size drift: {len(features)}")
    null = df[features].to_numpy(float)
    if not np.isfinite(null).all():
        raise RuntimeError("non-finite exact null statistic")
    sorted_abs = [np.sort(np.abs(null[:, j])) for j in range(null.shape[1])]
    corr = pd.DataFrame(null, columns=features).corr(method="pearson")

    cfg = contract["simulation"]
    reps = int(cfg["replicates"])
    effects = [float(x) for x in cfg["effect_grid"]]
    rng = np.random.default_rng(int(cfg["seed"]))
    rows = []
    for target_index, target in enumerate(features):
        sampled_index = rng.integers(0, len(null), size=reps)
        base = null[sampled_index].copy()
        for effect in effects:
            stat = base.copy()
            stat[:, target_index] += effect
            p_raw = empirical_two_sided_p(stat, sorted_abs)
            p_holm = holm_adjust_matrix(p_raw)
            one = p_raw[:, target_index] < 0.05
            target_holm = p_holm[:, target_index] < 0.05
            any_holm = np.any(p_holm < 0.05, axis=1)
            rows.append({
                "target_interaction": target,
                "true_target_effect": effect,
                "replicates": reps,
                "target_one_test_power": float(np.mean(one)),
                "target_holm_power": float(np.mean(target_holm)),
                "family_any_discovery_probability": float(np.mean(any_holm)),
                "holm_power_cost": float(np.mean(one) - np.mean(target_holm)),
            })
    power = pd.DataFrame(rows)
    mde = []
    for target, group in power.groupby("target_interaction", sort=False):
        g = group.sort_values("true_target_effect")
        hit = g[g["target_holm_power"] >= 0.8]
        mde.append({
            "target_interaction": target,
            "minimum_grid_effect_for_holm_power_ge_0_8": None if hit.empty else float(hit.iloc[0]["true_target_effect"]),
            "maximum_grid_power": float(g["target_holm_power"].max()),
        })
    mde_df = pd.DataFrame(mde)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    power.to_csv(args.output_dir / "global_rgfca_environmental_interaction_power_supplement_v1.csv", index=False)
    corr.rename_axis(index="interaction_a", columns="interaction_b").stack().rename("pearson_r").reset_index().to_csv(
        args.output_dir / "global_rgfca_environmental_interaction_null_correlations_v1.csv", index=False
    )
    mde_df.to_csv(args.output_dir / "global_rgfca_environmental_interaction_power_mde_v1.csv", index=False)
    observed = contract.get("observed_reference_only_after_simulation", {})
    observed_positions = {}
    for target, value in observed.items():
        if target == "rule" or target not in features:
            continue
        g = power[power["target_interaction"] == target].sort_values("true_target_effect")
        below = g[g["true_target_effect"] <= float(value)]
        above = g[g["true_target_effect"] >= float(value)]
        observed_positions[target] = {
            "observed_effect": float(value),
            "nearest_lower_grid_effect": None if below.empty else float(below.iloc[-1]["true_target_effect"]),
            "nearest_lower_holm_power": None if below.empty else float(below.iloc[-1]["target_holm_power"]),
            "nearest_upper_grid_effect": None if above.empty else float(above.iloc[0]["true_target_effect"]),
            "nearest_upper_holm_power": None if above.empty else float(above.iloc[0]["target_holm_power"]),
        }
    result = {
        "protocol": contract["protocol"],
        "status": "complete_exact_joint_null_interaction_family_power_supplement",
        "family_size": len(features),
        "features": features,
        "null_rows": len(null),
        "null_source_sha256": digest,
        "simulation_replicates_per_target": reps,
        "effect_grid": effects,
        "mde_grid": mde,
        "observed_reference_positions": observed_positions,
        "decision_guard": {
            "observed_interaction_family_reclassified": False,
            "variable_level_decomposition_opened": False,
            "family_of_one_retest_used": False,
        },
    }
    (args.output_dir / "global_rgfca_environmental_interaction_power_supplement_result_v1.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
