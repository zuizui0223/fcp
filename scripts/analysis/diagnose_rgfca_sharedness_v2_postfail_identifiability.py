#!/usr/bin/env python3
"""Post-fail identifiability audit for frozen RGFCA sharedness-v2 qualification.

This is a descriptive diagnostic only. It does not alter or rerun the frozen
qualification decision, tune any statistic/grid/sample size/power floor, open
image pixels, or authorize empirical measurement. It asks whether the shared
boundary injected by the frozen positive generator is actually observable on
the retained real metadata geometry for the species assigned to share it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_rgfca_sharedness_v2_synthetic_qualification_transport as transport
import run_rgfca_sharedness_v2_synthetic_qualification as base
from benchmark_rgfca_sharedness_v2_exact_kernel import canonical_axes, THRESHOLDS

ROOT = Path(__file__).resolve().parents[2]
QUAL = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_qualification_result_v1.json"
OUT_DEFAULT = ROOT / "docs/supporting/rgfca_sharedness_v2_postfail_identifiability_v1.json"
REPS = 250
HARD = (
    {"amplitude": 1.0, "shared_fraction": 1.0, "retention": 0.4},
    {"amplitude": 1.0, "shared_fraction": 0.5, "retention": 0.4},
    {"amplitude": 2.0, "shared_fraction": 0.5, "retention": 0.4},
)
THRESHOLD_SDS = (0.0, 0.25)


def generator_truth(root: int, arm: dict, rep: int):
    """Reproduce only the pre-colour truth draws of the frozen positive generator."""
    rng = np.random.default_rng(base.seed_for(root, "world", "evaluation", base.arm_id(arm), rep))
    normals = base.unit_vectors(rng, (150, 3))
    sd = float(arm["threshold_sd"])
    thresholds = rng.normal(0.0, sd, 150) if sd > 0 else np.zeros(150)
    common = base.unit_vectors(rng, (1, 3))[0]
    chosen = rng.permutation(150)[: round(150 * float(arm["shared_fraction"]))]
    normals[chosen] = common
    return common, thresholds, np.asarray(chosen, dtype=int)


def best_partition_agreement(x: np.ndarray, true_side: np.ndarray) -> float:
    """Best agreement with the frozen 96x5 candidate partition grid; descriptive only."""
    proj = x @ canonical_axes().T
    best = 0.5
    for t in THRESHOLDS:
        cand = proj > t
        for k in range(cand.shape[1]):
            a = float(np.mean(cand[:, k] == true_side))
            a = max(a, 1.0 - a)  # antipodal/sign invariance
            best = max(best, a)
    return best


def summarize(rows: pd.DataFrame, ids: np.ndarray) -> dict:
    z = rows[rows["species_id"].isin(set(map(int, ids)))]
    if z.empty:
        return {"shared_species_instances": 0}
    return {
        "shared_species_instances": int(len(z)),
        "straddle_fraction": float(z["straddles_true_boundary"].mean()),
        "balanced_10pct_fraction": float((z["minority_side_fraction"] >= 0.10).mean()),
        "balanced_20pct_fraction": float((z["minority_side_fraction"] >= 0.20).mean()),
        "median_minority_side_fraction": float(z["minority_side_fraction"].median()),
        "median_latent_dynamic_range": float(z["latent_dynamic_range"].median()),
        "median_best_frozen_grid_partition_agreement": float(z["best_grid_partition_agreement"].median()),
        "q10_best_frozen_grid_partition_agreement": float(z["best_grid_partition_agreement"].quantile(0.10)),
    }


def main(output: Path) -> int:
    q = json.loads(QUAL.read_text())
    if q.get("qualification_pass") is not False or q.get("image_acquisition_permitted") is not False:
        raise RuntimeError("post-fail audit requires immutable synthetic FAIL with acquisition blocked")
    if q.get("observed_flower_colour_opened") or q.get("observed_background_colour_opened") or q.get("image_pixels_opened"):
        raise RuntimeError("empirical firewall drift")

    _, mapping, _, xyz, train_ids, test_ids = transport.load_transport_inputs()
    root = int(mapping["execution"]["deterministic_seed_root"])
    species_rows = []
    world_rows = []

    for hard_index, h in enumerate(HARD):
        for sd in THRESHOLD_SDS:
            arm = {
                "id": "partially_shared_boundaries",
                "kind": "positive",
                "amplitude": float(h["amplitude"]),
                "shared_fraction": float(h["shared_fraction"]),
                "threshold_sd": float(sd),
                "retention": float(h["retention"]),
            }
            for rep in range(REPS):
                common, thresholds, chosen = generator_truth(root, arm, rep)
                idx = base.retention_indices(root, "evaluation", rep, float(h["retention"]), xyz, False, 0.0)
                wr = []
                for s in chosen:
                    take = idx[s]
                    x = xyz[s, take]
                    p = x @ common
                    side = p > thresholds[s]
                    frac = float(side.mean())
                    minority = min(frac, 1.0 - frac)
                    latent = np.tanh((p - thresholds[s]) / 0.15)
                    r = {
                        "hard_index": hard_index,
                        "amplitude": float(h["amplitude"]),
                        "shared_fraction": float(h["shared_fraction"]),
                        "retention": float(h["retention"]),
                        "threshold_sd": float(sd),
                        "replicate": rep,
                        "species_id": int(s),
                        "role": "training" if int(s) in set(map(int, train_ids)) else "evaluation",
                        "straddles_true_boundary": bool(side.any() and (~side).any()),
                        "minority_side_fraction": minority,
                        "latent_dynamic_range": float(latent.max() - latent.min()),
                        "best_grid_partition_agreement": best_partition_agreement(x, side),
                    }
                    species_rows.append(r)
                    wr.append(r)
                wdf = pd.DataFrame(wr)
                world_rows.append({
                    "hard_index": hard_index,
                    "amplitude": float(h["amplitude"]),
                    "shared_fraction": float(h["shared_fraction"]),
                    "retention": float(h["retention"]),
                    "threshold_sd": float(sd),
                    "replicate": rep,
                    "shared_species": int(len(chosen)),
                    "straddle_fraction": float(wdf["straddles_true_boundary"].mean()),
                    "balanced_10pct_fraction": float((wdf["minority_side_fraction"] >= 0.10).mean()),
                    "balanced_20pct_fraction": float((wdf["minority_side_fraction"] >= 0.20).mean()),
                    "median_latent_dynamic_range": float(wdf["latent_dynamic_range"].median()),
                    "median_best_grid_partition_agreement": float(wdf["best_grid_partition_agreement"].median()),
                })

    sdf = pd.DataFrame(species_rows)
    wdf = pd.DataFrame(world_rows)
    scenarios = []
    for (hi, sd), z in sdf.groupby(["hard_index", "threshold_sd"], sort=True):
        h = HARD[int(hi)]
        wz = wdf[(wdf["hard_index"] == hi) & (wdf["threshold_sd"] == sd)]
        scenarios.append({
            "hard_index": int(hi),
            "amplitude": float(h["amplitude"]),
            "shared_fraction": float(h["shared_fraction"]),
            "retention": float(h["retention"]),
            "threshold_sd": float(sd),
            "replicates": int(wz["replicate"].nunique()),
            "all_shared_species": summarize(z, np.arange(150)),
            "training_shared_species": summarize(z, train_ids),
            "evaluation_shared_species": summarize(z, test_ids),
            "world_level": {
                "median_straddle_fraction": float(wz["straddle_fraction"].median()),
                "q10_straddle_fraction": float(wz["straddle_fraction"].quantile(0.10)),
                "median_balanced_10pct_fraction": float(wz["balanced_10pct_fraction"].median()),
                "median_balanced_20pct_fraction": float(wz["balanced_20pct_fraction"].median()),
                "median_latent_dynamic_range": float(wz["median_latent_dynamic_range"].median()),
                "median_best_frozen_grid_partition_agreement": float(wz["median_best_grid_partition_agreement"].median()),
            },
        })

    payload = {
        "status": "complete_postfail_identifiability_audit",
        "protocol": "rgfca-sharedness-v2-postfail-identifiability-descriptive-v1",
        "parent_qualification_status": q["status"],
        "parent_qualification_pass": False,
        "parent_global_threshold": q["global_threshold"],
        "qualification_reopened": False,
        "qualification_decision_changed": False,
        "parameter_search_performed": False,
        "image_acquisition_permitted": False,
        "observed_flower_colour_opened": False,
        "observed_background_colour_opened": False,
        "image_pixels_opened": False,
        "diagnostic_question": "Does the frozen positive generator's true shared boundary have within-species observational support on the retained real metadata geometry, and can the frozen 96x5 grid geometrically approximate its partition?",
        "interpretation_rule": "Descriptive identifiability/support diagnosis only. It cannot convert FAIL to PASS or authorize empirical acquisition.",
        "worlds_audited": int(len(wdf)),
        "shared_species_instances_audited": int(len(sdf)),
        "scenarios": scenarios,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()
    raise SystemExit(main(args.output))
