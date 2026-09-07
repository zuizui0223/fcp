#!/usr/bin/env python3
"""Combine frozen threshold-robustness shards and verify all recorded decisions."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_hypervolume_threshold_robustness as run
import run_hypervolume_transfer_benchmark as base


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key_mask(df: pd.DataFrame, s: dict) -> pd.Series:
    mask = (df.world == s["world"]) & np.isclose(df.amplitude, float(s["amplitude"]))
    mask &= np.isclose(df.shared_fraction, float(s.get("shared_fraction", 0.0)))
    mask &= np.isclose(df.threshold_sd, float(s.get("threshold_sd", 0.0)))
    if "common_threshold" in s:
        mask &= np.isclose(df.common_threshold, float(s["common_threshold"]), equal_nan=False)
    else:
        mask &= df.common_threshold.isna()
    if "axis_concentration" in s:
        mask &= np.isclose(df.axis_concentration, float(s["axis_concentration"]), equal_nan=False)
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
    if len(score_files) != 20 or len(meta_files) != 20:
        raise RuntimeError(f"expected 20 score/meta shards, got {len(score_files)}/{len(meta_files)}")
    metas = [json.loads(p.read_text()) for p in meta_files]
    expected_ranges = {(stage, start, start + 25) for stage in ("calibration", "evaluation") for start in range(0, 250, 25)}
    got_ranges = {(m["stage"], int(m["replicate_start"]), int(m["replicate_stop"])) for m in metas}
    if got_ranges != expected_ranges:
        raise RuntimeError("shard range census mismatch")
    if any(m["specification_commit"] != run.SPEC for m in metas):
        raise RuntimeError("specification drift")
    if len({m["runner_sha256"] for m in metas}) != 1:
        raise RuntimeError("runner hash differs between shards")
    if len({m["geometry_sha256"] for m in metas}) != 1:
        raise RuntimeError("geometry hash differs between shards")

    df = pd.concat([pd.read_csv(p) for p in score_files], ignore_index=True)
    key = ["stage", "replicate", "model", "world", "amplitude", "shared_fraction",
           "threshold_sd", "common_threshold", "axis_concentration"]
    if len(df) != 22500 or df.duplicated(key).any():
        raise RuntimeError(f"score census/duplication failure: {len(df)}")
    if not np.isfinite(df[["score", "train_fraction"] + [f"score_fold{k}" for k in range(4)]].to_numpy()).all():
        raise RuntimeError("nonfinite numeric score")
    if np.max(np.abs(df[[f"score_fold{k}" for k in range(4)]].mean(axis=1) - df.score)) > 1e-12:
        raise RuntimeError("fold average mismatch")

    # Exact arm/repetition census.
    for stage, scenarios in (("calibration", run.NUISANCE), ("evaluation", run.NUISANCE + run.POSITIVE)):
        for model in run.MODELS:
            for s in scenarios:
                q = df[(df.stage == stage) & (df.model == model) & key_mask(df, s)]
                if len(q) != 250 or set(q.replicate) != set(range(250)):
                    raise RuntimeError(f"arm census failure {stage} {model} {s}: {len(q)}")

    quantiles, thresholds, rates, paired = [], [], [], []
    for model in run.MODELS:
        cuts = []
        for s in run.NUISANCE:
            q = df[(df.stage == "calibration") & (df.model == model) & key_mask(df, s)]
            cut = float(np.quantile(q.score, 0.975, method="higher"))
            quantiles.append({"model": model, **s, "quantile_975": cut})
            cuts.append(cut)
        threshold = max(cuts)
        thresholds.append(dict(model=model, threshold=threshold))
        for s in run.NUISANCE + run.POSITIVE:
            q = df[(df.stage == "evaluation") & (df.model == model) & key_mask(df, s)]
            count = int((q.score > threshold).sum())
            lo, hi = base.wilson(count, 250)
            rates.append({"model": model, **s, "count": count, "replicates": 250,
                          "rate": count / 250.0, "wilson95_low": lo, "wilson95_high": hi,
                          "threshold": threshold, "mean_score": float(q.score.mean()),
                          "mean_train_fraction": float(q.train_fraction.mean())})

    cutmap = {r["model"]: r["threshold"] for r in thresholds}
    baseline = run.MODELS[0]
    for s in run.POSITIVE:
        q = df[(df.stage == "evaluation") & key_mask(df, s)]
        wide = q.pivot(index="replicate", columns="model", values="score")
        b = (wide[baseline] > cutmap[baseline]).to_numpy()
        for candidate in run.MODELS[1:]:
            c = (wide[candidate] > cutmap[candidate]).to_numpy()
            diff = c.astype(float) - b.astype(float)
            paired.append({**s, "candidate": candidate, "baseline": baseline,
                           "gained": int((c & ~b).sum()), "lost": int((~c & b).sum()),
                           "paired_difference": float(diff.mean()),
                           "mc_se": float(diff.std(ddof=1) / np.sqrt(250))})

    args.output.mkdir(parents=True)
    combined = args.output / "scores.csv"
    df.sort_values(key).to_csv(combined, index=False, lineterminator="\n")
    for name, rows in (("calibration_quantiles.csv", quantiles), ("thresholds.csv", thresholds),
                       ("decision_rates.csv", rates), ("paired_detections.csv", paired)):
        pd.DataFrame(rows).to_csv(args.output / name, index=False, lineterminator="\n")

    nuisance_rates = [r for r in rates if any(all(r.get(k) == v for k, v in s.items()) for s in run.NUISANCE)]
    maximum_nuisance = max(nuisance_rates, key=lambda r: r["rate"])
    primary = [r for r in rates if r["model"] == "linear_logit_joint_train_intercept_marginal" and
               r["world"].startswith("environmental_shared")]
    summary = {
        "protocol": "hypervolume-threshold-robustness-v1",
        "status": "complete_pending_external_ci_confirmation",
        "specification_commit": run.SPEC,
        "seed": run.SEED,
        "unique_worlds": 7500,
        "calibration_worlds": 3000,
        "evaluation_worlds": 4500,
        "score_records": len(df),
        "models": list(run.MODELS),
        "primary_candidate": "linear_logit_joint_train_intercept_marginal",
        "primary_positive_rates": primary,
        "maximum_nuisance_rejection": maximum_nuisance,
        "paired_candidate_vs_baseline": [x for x in paired if x["candidate"] == "linear_logit_joint_train_intercept_marginal"],
        "verification": {
            "twenty_nonoverlapping_shards": True,
            "all_22500_unique_finite_scores": True,
            "all_four_folds_equal_weight": True,
            "all_arm_replicate_censuses": True,
            "all_36_calibration_quantiles_recomputed": True,
            "all_3_thresholds_recomputed": True,
            "all_54_evaluation_counts_and_wilson_intervals_recomputed": True,
            "runner_sha256": next(iter({m["runner_sha256"] for m in metas})),
            "geometry_sha256": next(iter({m["geometry_sha256"] for m in metas})),
            "combined_scores_sha256": digest(combined),
            "decision_rates_sha256": digest(args.output / "decision_rates.csv"),
        },
        "biological_colour_values_read": False,
        "real_climate_used": False,
        "parent_decisions_modified": False,
        "empirical_inference_opened": False,
        "claim_ceiling": "Synthetic threshold/richer-null robustness on fixed actual-coordinate support only; not real flower-colour or climate evidence and not a universally certified <=5% test."
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
