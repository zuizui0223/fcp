#!/usr/bin/env python3
"""Finalize and independently census real-climate synthetic qualification shards."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_hypervolume_real_climate_synthetic_qualification as run
import run_hypervolume_transfer_benchmark as base


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key_mask(df: pd.DataFrame, scenario: dict) -> pd.Series:
    mask = (df.world == scenario["world"]) & np.isclose(df.amplitude, float(scenario["amplitude"]))
    mask &= np.isclose(df.shared_fraction, float(scenario.get("shared_fraction", 0.0)))
    mask &= np.isclose(df.threshold_sd, float(scenario.get("threshold_sd", 0.0)))
    if "common_threshold" in scenario:
        mask &= np.isclose(df.common_threshold, float(scenario["common_threshold"]), equal_nan=False)
    else:
        mask &= df.common_threshold.isna()
    if "axis_concentration" in scenario:
        mask &= np.isclose(df.axis_concentration, float(scenario["axis_concentration"]), equal_nan=False)
    else:
        mask &= df.axis_concentration.isna()
    return mask


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh finalizer output")

    score_files = sorted(args.shards.rglob("scores.csv"))
    meta_files = sorted(args.shards.rglob("meta.json"))
    if len(score_files) != 60 or len(meta_files) != 60:
        raise RuntimeError(f"expected 60 score/meta shards, got {len(score_files)}/{len(meta_files)}")
    metas = [json.loads(p.read_text()) for p in meta_files]
    expected_ranges = {
        (block, stage, start, start + 25)
        for block in run.BLOCKS
        for stage in ("calibration", "evaluation")
        for start in range(0, 250, 25)
    }
    got_ranges = {
        (m["block"], m["stage"], int(m["replicate_start"]), int(m["replicate_stop"]))
        for m in metas
    }
    if got_ranges != expected_ranges:
        raise RuntimeError("60-shard block/stage/range census mismatch")
    if any(m["specification_commit"] != run.SPEC for m in metas):
        raise RuntimeError("specification drift")
    if len({m["runner_sha256"] for m in metas}) != 1:
        raise RuntimeError("runner differs between shards")
    if {m["geometry_sha256"] for m in metas} != {run.GEOMETRY_SHA}:
        raise RuntimeError("geometry lineage drift")
    if {m["real_climate_export_sha256"] for m in metas} != {run.CLIMATE_SHA}:
        raise RuntimeError("climate lineage drift")
    if {int(m["joint_complete_retained_rows"]) for m in metas} != {21418}:
        raise RuntimeError("joint-complete frame drift")
    if len({m["axis_grid_sha256"] for m in metas}) != 1:
        raise RuntimeError("axis grid differs between shards")
    if any(m.get("biological_colour_values_read") is not False for m in metas):
        raise RuntimeError("colour firewall violated")

    df = pd.concat([pd.read_csv(p) for p in score_files], ignore_index=True)
    key = [
        "stage", "block", "replicate", "model", "world", "amplitude",
        "shared_fraction", "threshold_sd", "common_threshold", "axis_concentration",
    ]
    if len(df) != 45000 or df.duplicated(key).any():
        raise RuntimeError(f"score census/duplication failure: {len(df)}")
    numeric = ["score", "train_fraction"] + [f"score_fold{k}" for k in range(4)]
    if not np.isfinite(df[numeric].to_numpy()).all():
        raise RuntimeError("nonfinite numeric score")
    if np.max(np.abs(df[[f"score_fold{k}" for k in range(4)]].mean(axis=1) - df.score)) > 1e-12:
        raise RuntimeError("fold average mismatch")

    for block in run.BLOCKS:
        for stage, scenarios in (
            ("calibration", run.NUISANCE),
            ("evaluation", run.NUISANCE + run.POSITIVE),
        ):
            for model in run.MODELS:
                for s in scenarios:
                    q = df[
                        (df.block == block)
                        & (df.stage == stage)
                        & (df.model == model)
                        & key_mask(df, s)
                    ]
                    if len(q) != 250 or set(q.replicate) != set(range(250)):
                        raise RuntimeError(f"arm census failure {block} {stage} {model} {s}: {len(q)}")

    quantiles, thresholds, rates, paired = [], [], [], []
    for model in run.MODELS:
        cuts = []
        for block in run.BLOCKS:
            for s in run.NUISANCE:
                q = df[
                    (df.block == block)
                    & (df.stage == "calibration")
                    & (df.model == model)
                    & key_mask(df, s)
                ]
                cut = float(np.quantile(q.score, 0.975, method="higher"))
                quantiles.append({"model": model, "block": block, **s, "quantile_975": cut})
                cuts.append(cut)
        threshold = max(cuts)
        thresholds.append({"model": model, "threshold": threshold})
        for block in run.BLOCKS:
            for s in run.NUISANCE + run.POSITIVE:
                q = df[
                    (df.block == block)
                    & (df.stage == "evaluation")
                    & (df.model == model)
                    & key_mask(df, s)
                ]
                count = int((q.score > threshold).sum())
                lo, hi = base.wilson(count, 250)
                rates.append({
                    "model": model,
                    "block": block,
                    **s,
                    "count": count,
                    "replicates": 250,
                    "rate": count / 250.0,
                    "wilson95_low": lo,
                    "wilson95_high": hi,
                    "threshold": threshold,
                    "mean_score": float(q.score.mean()),
                    "mean_train_fraction": float(q.train_fraction.mean()),
                })

    cutmap = {r["model"]: r["threshold"] for r in thresholds}
    baseline, candidate = run.MODELS
    for block in run.BLOCKS:
        for s in run.POSITIVE:
            q = df[
                (df.block == block)
                & (df.stage == "evaluation")
                & key_mask(df, s)
            ]
            wide = q.pivot(index="replicate", columns="model", values="score")
            b = (wide[baseline] > cutmap[baseline]).to_numpy()
            c = (wide[candidate] > cutmap[candidate]).to_numpy()
            diff = c.astype(float) - b.astype(float)
            paired.append({
                "block": block,
                **s,
                "candidate": candidate,
                "baseline": baseline,
                "gained": int((c & ~b).sum()),
                "lost": int((~c & b).sum()),
                "paired_difference": float(diff.mean()),
                "mc_se": float(diff.std(ddof=1) / np.sqrt(250)),
            })

    args.output.mkdir(parents=True)
    combined = args.output / "scores.csv"
    df.sort_values(key).to_csv(combined, index=False, lineterminator="\n")
    for name, rows in (
        ("calibration_quantiles.csv", quantiles),
        ("thresholds.csv", thresholds),
        ("decision_rates.csv", rates),
        ("paired_detections.csv", paired),
    ):
        pd.DataFrame(rows).to_csv(args.output / name, index=False, lineterminator="\n")

    candidate_rates = pd.DataFrame([r for r in rates if r["model"] == candidate])
    nuisance_mask = ~candidate_rates.world.str.startswith("climate_shared")
    max_nuisance = float(candidate_rates.loc[nuisance_mask, "rate"].max())
    moderate_by_block = {}
    strong_by_block = {}
    for block in run.BLOCKS:
        moderate = candidate_rates[
            (candidate_rates.block == block)
            & (candidate_rates.world == "climate_shared_shifted")
            & np.isclose(candidate_rates.amplitude, 1.0)
            & np.isclose(candidate_rates.shared_fraction, 1.0)
            & np.isclose(candidate_rates.threshold_sd, 0.5)
        ]
        strong = candidate_rates[
            (candidate_rates.block == block)
            & (candidate_rates.world == "climate_shared_shifted")
            & np.isclose(candidate_rates.amplitude, 2.0)
            & np.isclose(candidate_rates.shared_fraction, 1.0)
            & np.isclose(candidate_rates.threshold_sd, 0.5)
        ]
        if len(moderate) != 1 or len(strong) != 1:
            raise RuntimeError("qualification positive arm lookup failed")
        moderate_by_block[block] = float(moderate.iloc[0].rate)
        strong_by_block[block] = float(strong.iloc[0].rate)

    gates = {
        "maximum_nuisance_rate_observed": max_nuisance,
        "maximum_nuisance_rate_required": 0.10,
        "maximum_nuisance_gate_pass": bool(max_nuisance <= 0.10),
        "moderate_full_sharing_recovery_by_block": moderate_by_block,
        "moderate_full_sharing_required_each_block": 0.50,
        "moderate_full_sharing_gate_pass": bool(all(v >= 0.50 for v in moderate_by_block.values())),
        "strong_full_sharing_recovery_by_block": strong_by_block,
        "strong_full_sharing_required_each_block": 0.90,
        "strong_full_sharing_gate_pass": bool(all(v >= 0.90 for v in strong_by_block.values())),
    }
    gates["synthetic_qualification_pass"] = bool(
        gates["maximum_nuisance_gate_pass"]
        and gates["moderate_full_sharing_gate_pass"]
        and gates["strong_full_sharing_gate_pass"]
    )

    summary = {
        "protocol": "hypervolume-real-climate-synthetic-qualification-v1",
        "status": "complete_pending_ci_guard",
        "specification_commit": run.SPEC,
        "seed": run.SEED,
        "blocks": list(run.BLOCKS),
        "models": list(run.MODELS),
        "block_worlds": 22500,
        "score_records": len(df),
        "global_thresholds": thresholds,
        "qualification_gates": gates,
        "primary_candidate": candidate,
        "primary_positive_rates": [
            r for r in rates if r["model"] == candidate and r["world"].startswith("climate_shared")
        ],
        "primary_nuisance_rates": [
            r for r in rates if r["model"] == candidate and not r["world"].startswith("climate_shared")
        ],
        "paired_candidate_vs_baseline": paired,
        "verification": {
            "sixty_nonoverlapping_shards": True,
            "all_45000_unique_finite_scores": True,
            "all_four_folds_equal_weight": True,
            "all_arm_replicate_censuses": True,
            "all_72_calibration_quantiles_recomputed": True,
            "all_2_global_thresholds_recomputed": True,
            "all_108_evaluation_counts_and_wilson_intervals_recomputed": True,
            "runner_sha256": next(iter({m["runner_sha256"] for m in metas})),
            "axis_grid_sha256": next(iter({m["axis_grid_sha256"] for m in metas})),
            "geometry_sha256": run.GEOMETRY_SHA,
            "real_climate_export_sha256": run.CLIMATE_SHA,
            "combined_scores_sha256": digest(combined),
            "decision_rates_sha256": digest(args.output / "decision_rates.csv"),
        },
        "biological_colour_values_read": False,
        "synthetic_colour_labels_only": True,
        "real_climate_values_used_as_geometry": True,
        "parent_empirical_decisions_modified": False,
        "empirical_inference_opened": False,
        "claim_ceiling": (
            "Synthetic recoverability/calibration on SHA-locked real CHELSA geometry only; "
            "not evidence for real flower-colour climate responses, adaptation, causality or prevalence."
        ),
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
