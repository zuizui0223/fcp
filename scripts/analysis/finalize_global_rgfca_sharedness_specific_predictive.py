#!/usr/bin/env python3
"""Finalize frozen sharedness-specific predictive shards and verify every decision."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_global_rgfca_sharedness_specific_predictive as run
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


def summarize_inclusions(meta_files: list[Path], output: Path) -> dict:
    species_parts = []
    photo_parts = []
    for meta_path in meta_files:
        meta = json.loads(meta_path.read_text())
        stage = meta["stage"]
        sp = pd.read_csv(meta_path.with_name("species_inclusions.csv"))
        ph = pd.read_csv(meta_path.with_name("photo_inclusions.csv"))
        sp["stage"] = stage
        ph["stage"] = stage
        species_parts.append(sp)
        photo_parts.append(ph)

    species = pd.concat(species_parts, ignore_index=True)
    agg_species = (
        species.groupby(["stage", "species"], as_index=False)[["training_inclusions", "evaluation_inclusions"]]
        .sum()
        .sort_values(["stage", "species"], kind="stable")
    )
    for stage in ("calibration", "evaluation"):
        q = agg_species[agg_species.stage == stage]
        if len(q) != 369:
            raise RuntimeError(f"species inclusion census failure for {stage}: {len(q)}")
        if int(q.training_inclusions.sum()) != run.REPS * run.N_TRAIN:
            raise RuntimeError(f"training species-slot total mismatch for {stage}")
        if int(q.evaluation_inclusions.sum()) != run.REPS * run.N_TEST:
            raise RuntimeError(f"evaluation species-slot total mismatch for {stage}")
        if int((q.training_inclusions > 0).sum()) != 184 or int((q.evaluation_inclusions > 0).sum()) != 185:
            raise RuntimeError(f"fixed 184/185 split not recovered for {stage}")
        if ((q.training_inclusions > 0) & (q.evaluation_inclusions > 0)).any():
            raise RuntimeError(f"training/evaluation species overlap in {stage}")

    photos = pd.concat(photo_parts, ignore_index=True)
    agg_photos = (
        photos.groupby(["stage", "photo_id"], as_index=False).inclusions.sum()
        .sort_values(["stage", "photo_id"], kind="stable")
    )
    for stage in ("calibration", "evaluation"):
        total = int(agg_photos.loc[agg_photos.stage == stage, "inclusions"].sum())
        expected = run.REPS * (run.N_TRAIN + run.N_TEST) * run.N_PHOTOS
        if total != expected:
            raise RuntimeError(f"photo-slot total mismatch for {stage}: {total} != {expected}")

    agg_species.to_csv(output / "species_inclusion_census.csv", index=False, lineterminator="\n")
    agg_photos.to_csv(output / "photo_inclusion_census.csv", index=False, lineterminator="\n")
    return {
        "calibration_training_species_inclusion_min": int(
            agg_species[(agg_species.stage == "calibration") & (agg_species.training_inclusions > 0)].training_inclusions.min()
        ),
        "calibration_training_species_inclusion_max": int(
            agg_species[(agg_species.stage == "calibration") & (agg_species.training_inclusions > 0)].training_inclusions.max()
        ),
        "calibration_evaluation_species_inclusion_min": int(
            agg_species[(agg_species.stage == "calibration") & (agg_species.evaluation_inclusions > 0)].evaluation_inclusions.min()
        ),
        "calibration_evaluation_species_inclusion_max": int(
            agg_species[(agg_species.stage == "calibration") & (agg_species.evaluation_inclusions > 0)].evaluation_inclusions.max()
        ),
        "evaluation_training_species_inclusion_min": int(
            agg_species[(agg_species.stage == "evaluation") & (agg_species.training_inclusions > 0)].training_inclusions.min()
        ),
        "evaluation_training_species_inclusion_max": int(
            agg_species[(agg_species.stage == "evaluation") & (agg_species.training_inclusions > 0)].training_inclusions.max()
        ),
        "evaluation_evaluation_species_inclusion_min": int(
            agg_species[(agg_species.stage == "evaluation") & (agg_species.evaluation_inclusions > 0)].evaluation_inclusions.min()
        ),
        "evaluation_evaluation_species_inclusion_max": int(
            agg_species[(agg_species.stage == "evaluation") & (agg_species.evaluation_inclusions > 0)].evaluation_inclusions.max()
        ),
        "calibration_unique_photos_used": int((agg_photos.stage == "calibration").sum()),
        "evaluation_unique_photos_used": int((agg_photos.stage == "evaluation").sum()),
    }


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
    expected_ranges = {
        (stage, start, start + 25)
        for stage in ("calibration", "evaluation")
        for start in range(0, run.REPS, 25)
    }
    got_ranges = {
        (m["stage"], int(m["replicate_start"]), int(m["replicate_stop"])) for m in metas
    }
    if got_ranges != expected_ranges:
        raise RuntimeError("shard range census mismatch")
    if any(m["specification_commit"] != run.SPECIFICATION_COMMIT for m in metas):
        raise RuntimeError("specification drift")
    if any(m["geometry_sha256"] != run.GEOMETRY_SHA for m in metas):
        raise RuntimeError("geometry checksum drift across shards")
    if len({m["runner_sha256"] for m in metas}) != 1 or len({m["axis_grid_sha256"] for m in metas}) != 1:
        raise RuntimeError("runner or axis-grid hash differs between shards")
    if any(m.get("biological_colour_values_read") is not False for m in metas):
        raise RuntimeError("biological-colour firewall failure")
    if any(m.get("structured_independent_boundaries_in_null") is not True for m in metas):
        raise RuntimeError("structured independent-boundary null guard failure")
    if any(m.get("parent_G1_reclassified") is not False for m in metas):
        raise RuntimeError("parent G1 decision guard failure")

    df = pd.concat([pd.read_csv(p) for p in score_files], ignore_index=True)
    key = [
        "stage", "replicate", "world", "amplitude", "shared_fraction",
        "threshold_sd", "common_threshold", "axis_concentration",
    ]
    if len(df) != 7000 or df.duplicated(key).any():
        raise RuntimeError(f"score census/duplication failure: {len(df)}")
    numeric = [
        "score", "posterior_mean_shared_fraction",
        "posterior_axis_entropy_normalized", "posterior_max_axis_mass",
    ]
    if not np.isfinite(df[numeric].to_numpy()).all():
        raise RuntimeError("nonfinite score/posterior diagnostic")
    if not df.posterior_mean_shared_fraction.between(0.0, 1.0).all():
        raise RuntimeError("posterior mean shared fraction outside [0,1]")
    if not df.posterior_axis_entropy_normalized.between(-1e-12, 1.0 + 1e-12).all():
        raise RuntimeError("normalized axis entropy outside [0,1]")
    if not df.posterior_max_axis_mass.between(0.0, 1.0).all():
        raise RuntimeError("posterior max axis mass outside [0,1]")

    for stage, scenarios in (
        ("calibration", run.NUISANCE),
        ("evaluation", run.NUISANCE + run.POSITIVE),
    ):
        for scenario in scenarios:
            q = df[(df.stage == stage) & key_mask(df, scenario)]
            if len(q) != run.REPS or set(q.replicate) != set(range(run.REPS)):
                raise RuntimeError(f"arm census failure {stage} {scenario}: {len(q)}")

    calibration = []
    cuts = []
    for scenario in run.NUISANCE:
        q = df[(df.stage == "calibration") & key_mask(df, scenario)]
        cut = float(np.quantile(q.score, 0.975, method="higher"))
        calibration.append({**scenario, "quantile_975": cut})
        cuts.append(cut)
    threshold = float(max(cuts))

    rates = []
    for scenario in run.NUISANCE + run.POSITIVE:
        q = df[(df.stage == "evaluation") & key_mask(df, scenario)]
        count = int((q.score > threshold).sum())
        lo, hi = base.wilson(count, run.REPS)
        rates.append({
            **scenario,
            "count": count,
            "replicates": run.REPS,
            "rate": count / run.REPS,
            "wilson95_low": lo,
            "wilson95_high": hi,
            "threshold": threshold,
            "mean_score": float(q.score.mean()),
            "median_score": float(q.score.median()),
            "mean_posterior_shared_fraction": float(q.posterior_mean_shared_fraction.mean()),
            "median_posterior_shared_fraction": float(q.posterior_mean_shared_fraction.median()),
            "mean_axis_entropy_normalized": float(q.posterior_axis_entropy_normalized.mean()),
            "mean_max_axis_mass": float(q.posterior_max_axis_mass.mean()),
        })

    nuisance_rates = rates[: len(run.NUISANCE)]
    positive_rates = rates[len(run.NUISANCE):]
    worst_nuisance = max(nuisance_rates, key=lambda x: x["rate"])

    def positive_rate(amplitude: float, shared_fraction: float, threshold_sd: float = 0.25) -> float:
        matches = [
            r for r in positive_rates
            if r["world"] == "partially_shared_boundaries"
            and np.isclose(r["amplitude"], amplitude)
            and np.isclose(r["shared_fraction"], shared_fraction)
            and np.isclose(r["threshold_sd"], threshold_sd)
        ]
        if len(matches) != 1:
            raise RuntimeError(f"positive gate arm missing/duplicated: amp={amplitude}, share={shared_fraction}")
        return float(matches[0]["rate"])

    gates = {
        "maximum_evaluation_nuisance_rejection_rate_observed": float(worst_nuisance["rate"]),
        "maximum_evaluation_nuisance_rejection_rate_required": 0.10,
        "maximum_evaluation_nuisance_rejection_gate_pass": bool(worst_nuisance["rate"] <= 0.10),
        "amp1_shared_fraction_1_recovery": positive_rate(1.0, 1.0),
        "amp1_shared_fraction_1_required": 0.80,
        "amp1_shared_fraction_1_gate_pass": bool(positive_rate(1.0, 1.0) >= 0.80),
        "amp1_shared_fraction_0_5_recovery": positive_rate(1.0, 0.5),
        "amp1_shared_fraction_0_5_required": 0.50,
        "amp1_shared_fraction_0_5_gate_pass": bool(positive_rate(1.0, 0.5) >= 0.50),
        "amp2_shared_fraction_0_5_recovery": positive_rate(2.0, 0.5),
        "amp2_shared_fraction_0_5_required": 0.80,
        "amp2_shared_fraction_0_5_gate_pass": bool(positive_rate(2.0, 0.5) >= 0.80),
    }
    gates["synthetic_qualification_pass"] = bool(
        gates["maximum_evaluation_nuisance_rejection_gate_pass"]
        and gates["amp1_shared_fraction_1_gate_pass"]
        and gates["amp1_shared_fraction_0_5_gate_pass"]
        and gates["amp2_shared_fraction_0_5_gate_pass"]
    )

    args.output.mkdir(parents=True)
    combined = args.output / "scores.csv"
    df.sort_values(key, kind="stable").to_csv(combined, index=False, lineterminator="\n")
    pd.DataFrame(calibration).to_csv(args.output / "calibration_quantiles.csv", index=False, lineterminator="\n")
    pd.DataFrame([{"threshold": threshold}]).to_csv(args.output / "threshold.csv", index=False, lineterminator="\n")
    pd.DataFrame(rates).to_csv(args.output / "decision_rates.csv", index=False, lineterminator="\n")
    inclusion_summary = summarize_inclusions(meta_files, args.output)

    summary = {
        "protocol": "global-rgfca-sharedness-specific-predictive-v1",
        "status": "complete_pending_external_ci_confirmation",
        "specification_commit": run.SPECIFICATION_COMMIT,
        "seed": run.SEED,
        "axis_seed": run.AXIS_SEED,
        "eligible_species": 369,
        "training_species_pool": 184,
        "evaluation_species_pool": 185,
        "training_species_per_replicate": run.N_TRAIN,
        "evaluation_species_per_replicate": run.N_TEST,
        "photos_per_species": run.N_PHOTOS,
        "calibration_worlds": 2250,
        "evaluation_nuisance_worlds": 2250,
        "evaluation_positive_worlds": 2500,
        "total_worlds": len(df),
        "global_threshold": threshold,
        "qualification_gates": gates,
        "worst_evaluation_nuisance": worst_nuisance,
        "positive_rates": positive_rates,
        "inclusion_summary": inclusion_summary,
        "verification": {
            "twenty_nonoverlapping_shards": True,
            "all_7000_unique_finite_scores": True,
            "all_9_calibration_quantiles_recomputed": True,
            "global_threshold_recomputed_as_max_null_q975": True,
            "all_19_evaluation_arm_counts_and_wilson_intervals_recomputed": True,
            "fixed_184_185_species_split_recovered": True,
            "species_and_photo_slot_censuses_verified": True,
            "runner_sha256": next(iter({m["runner_sha256"] for m in metas})),
            "axis_grid_sha256": next(iter({m["axis_grid_sha256"] for m in metas})),
            "geometry_sha256": run.GEOMETRY_SHA,
            "combined_scores_sha256": digest(combined),
            "decision_rates_sha256": digest(args.output / "decision_rates.csv"),
        },
        "biological_colour_values_read": False,
        "real_climate_values_used": False,
        "structured_independent_boundaries_in_null": True,
        "parent_G1_reclassified": False,
        "six_species_or_34species_used": False,
        "claim_ceiling": "Synthetic qualification of a species-disjoint cross-species geographic sharedness statistic on the 369-species RGFCA geometry only. It does not establish real shared flower-colour boundaries and cannot rescue G1 or any prior empirical result.",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
