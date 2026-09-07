#!/usr/bin/env python3
"""Finalize post-failure real-climate identifiability diagnosis."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_hypervolume_real_climate_identifiability_diagnostic as diag


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scenario_mask(df: pd.DataFrame, scenario: dict) -> pd.Series:
    return diag.scenario_mask(df, scenario)


def qstats(x: pd.Series) -> dict[str, float]:
    a = x.to_numpy(float)
    return {
        "mean": float(np.mean(a)),
        "median": float(np.median(a)),
        "p10": float(np.quantile(a, 0.10)),
        "p90": float(np.quantile(a, 0.90)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=Path, required=True)
    ap.add_argument("--parent-rates", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh diagnostic finalizer output")

    row_files = sorted(args.shards.rglob("diagnostic_rows.csv"))
    meta_files = sorted(args.shards.rglob("meta.json"))
    if len(row_files) != 30 or len(meta_files) != 30:
        raise RuntimeError(f"expected 30 diagnostic row/meta shards, got {len(row_files)}/{len(meta_files)}")
    metas = [json.loads(p.read_text()) for p in meta_files]
    expected = {
        (block, start, start + 25)
        for block in diag.parent.BLOCKS
        for start in range(0, 250, 25)
    }
    got = {(m["block"], int(m["replicate_start"]), int(m["replicate_stop"])) for m in metas}
    if got != expected:
        raise RuntimeError("diagnostic shard census mismatch")
    if {m["specification_commit"] for m in metas} != {diag.SPECIFICATION_COMMIT}:
        raise RuntimeError("diagnostic specification drift")
    if {m["geometry_sha256"] for m in metas} != {diag.parent.GEOMETRY_SHA}:
        raise RuntimeError("geometry lineage drift")
    if {m["real_climate_export_sha256"] for m in metas} != {diag.parent.CLIMATE_SHA}:
        raise RuntimeError("climate lineage drift")
    if len({m["runner_sha256"] for m in metas}) != 1:
        raise RuntimeError("runner differs across diagnostic shards")
    if any(m.get("biological_colour_values_read") is not False for m in metas):
        raise RuntimeError("colour firewall violated")
    if any(m.get("parent_failure_reclassified") is not False for m in metas):
        raise RuntimeError("parent failure guard violated")

    df = pd.concat([pd.read_csv(p) for p in row_files], ignore_index=True)
    key = ["block", "replicate", "world", "amplitude", "shared_fraction", "threshold_sd", "common_threshold"]
    if len(df) != 3000 or df.duplicated(key).any():
        raise RuntimeError(f"diagnostic row census/duplication failure: {len(df)}")
    numeric = [
        "parent_score", "augmented_true_axis_score", "augmented_minus_parent_score",
        "augmented_mean_sharing_fraction", "privileged_true_axis_posterior_mass",
        "true_axis_vector_coverage", "true_axis_projection_corr_mean",
        "true_axis_projection_corr_min_fold", "train_transition_exposure",
        "evaluation_transition_exposure", "train_colour_variability",
        "evaluation_colour_variability", "train_true_projection_sd",
        "evaluation_true_projection_sd", "parent_score_reconstruction_abs_error",
    ] + [f"true_axis_projection_corr_fold{i}" for i in range(4)]
    if not np.isfinite(df[numeric].to_numpy()).all():
        raise RuntimeError("nonfinite diagnostic metric")
    max_reconstruction = float(df.parent_score_reconstruction_abs_error.max())
    if max_reconstruction > 2e-12:
        raise RuntimeError("parent score reconstruction exceeds frozen tolerance")

    parent_rates = pd.read_csv(args.parent_rates)
    for block in diag.parent.BLOCKS:
        for scenario in diag.SCENARIOS:
            q = df[(df.block == block) & scenario_mask(df, scenario)]
            if len(q) != 250 or set(q.replicate) != set(range(250)):
                raise RuntimeError(f"diagnostic scenario census failure {block} {scenario}")
            observed = int(q.parent_detected.astype(bool).sum())
            p = parent_rates[
                (parent_rates.model == diag.PARENT_CANDIDATE)
                & (parent_rates.block == block)
                & scenario_mask(parent_rates, scenario)
            ]
            if len(p) != 1 or observed != int(p.iloc[0]["count"]):
                raise RuntimeError(f"parent decision reconstruction mismatch {block} {scenario}: {observed}")

    summary_rows = []
    metric_names = [
        "parent_score",
        "augmented_true_axis_score",
        "augmented_minus_parent_score",
        "privileged_true_axis_posterior_mass",
        "true_axis_vector_coverage",
        "true_axis_projection_corr_mean",
        "true_axis_projection_corr_min_fold",
        "train_transition_exposure",
        "evaluation_transition_exposure",
        "train_colour_variability",
        "evaluation_colour_variability",
        "train_true_projection_sd",
        "evaluation_true_projection_sd",
    ]
    for block in diag.parent.BLOCKS:
        for scenario in diag.SCENARIOS:
            q = df[(df.block == block) & scenario_mask(df, scenario)].copy()
            row = {
                "block": block,
                **scenario,
                "replicates": len(q),
                "parent_detection_rate": float(q.parent_detected.astype(bool).mean()),
            }
            for metric in metric_names:
                stats = qstats(q[metric])
                row.update({f"{metric}_{k}": v for k, v in stats.items()})
            summary_rows.append(row)

    assoc_rows = []
    assoc_metrics = [
        "augmented_minus_parent_score",
        "true_axis_vector_coverage",
        "true_axis_projection_corr_mean",
        "train_transition_exposure",
        "evaluation_transition_exposure",
        "train_colour_variability",
        "evaluation_colour_variability",
        "train_true_projection_sd",
        "evaluation_true_projection_sd",
    ]
    for block in diag.parent.BLOCKS:
        qb = df[df.block == block]
        for scenario_label, qs in [("all_four_full_sharing_scenarios", qb)] + [
            (diag.parent.scenario_id(s), qb[scenario_mask(qb, s)]) for s in diag.SCENARIOS
        ]:
            for metric in assoc_metrics:
                rho = float(qs[["parent_score", metric]].corr(method="spearman").iloc[0, 1])
                assoc_rows.append({
                    "block": block,
                    "scenario": scenario_label,
                    "metric": metric,
                    "spearman_parent_score": rho,
                    "n": len(qs),
                })

    args.output.mkdir(parents=True)
    rows_path = args.output / "diagnostic_rows.csv"
    summary_path = args.output / "diagnostic_summary.csv"
    assoc_path = args.output / "diagnostic_associations.csv"
    df.sort_values(key).to_csv(rows_path, index=False, lineterminator="\n")
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False, lineterminator="\n")
    pd.DataFrame(assoc_rows).to_csv(assoc_path, index=False, lineterminator="\n")

    compact = {
        "protocol": "hypervolume-real-climate-identifiability-diagnostic-v1",
        "status": "complete_pending_ci_guard",
        "specification_commit": diag.SPECIFICATION_COMMIT,
        "parent_qualification_decision": "FAIL",
        "parent_qualification_reclassified": False,
        "rows": len(df),
        "blocks": list(diag.parent.BLOCKS),
        "scenarios_per_block": len(diag.SCENARIOS),
        "maximum_parent_score_reconstruction_abs_error": max_reconstruction,
        "summaries": summary_rows,
        "verification": {
            "thirty_nonoverlapping_shards": True,
            "all_3000_rows_unique_finite": True,
            "all_parent_scores_reconstructed": True,
            "all_parent_decision_counts_reconstructed": True,
            "runner_sha256": next(iter({m["runner_sha256"] for m in metas})),
            "geometry_sha256": diag.parent.GEOMETRY_SHA,
            "real_climate_export_sha256": diag.parent.CLIMATE_SHA,
            "diagnostic_rows_sha256": digest(rows_path),
            "diagnostic_summary_sha256": digest(summary_path),
            "diagnostic_associations_sha256": digest(assoc_path),
        },
        "biological_colour_values_read": False,
        "privileged_true_axis_is_empirically_available": False,
        "privileged_axis_score_used_for_inference": False,
        "empirical_inference_opened": False,
        "claim_ceiling": (
            "Post-failure synthetic identifiability diagnosis only; cannot rescue the failed Stage-6 qualification "
            "and is not evidence for real flower-colour responses to any climate block."
        ),
    }
    (args.output / "result.json").write_text(json.dumps(compact, indent=2) + "\n")
    print(json.dumps(compact, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
