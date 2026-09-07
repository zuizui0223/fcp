#!/usr/bin/env python3
"""Diagnose the failed real-climate synthetic qualification without rescuing it.

Uses the exact Stage-6 evaluation worlds and SHA-locked real-climate geometry.
Observed flower-colour values are forbidden. The simulator true shared axis is used
only as a privileged diagnostic augmentation and never receives an inferential gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

import run_hypervolume_real_climate_synthetic_qualification as parent
import run_hypervolume_transfer_benchmark as base

SPECIFICATION_COMMIT = "822fcc03b9b0db1265a07ac398f68c626244d9aa"
PARENT_CANDIDATE = "linear_logit_joint_train_intercept_marginal"
PARENT_SCORE_THRESHOLD = 0.0078392153007344
SCENARIOS = (
    dict(world="climate_shared_shifted", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.5),
    dict(world="climate_shared_shifted", amplitude=1.0, shared_fraction=1.0, threshold_sd=1.0),
    dict(world="climate_shared_shifted", amplitude=2.0, shared_fraction=1.0, threshold_sd=0.5),
    dict(world="climate_shared_common_offset", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.0, common_threshold=0.5),
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scenario_mask(df: pd.DataFrame, scenario: dict) -> pd.Series:
    m = (df.world == scenario["world"]) & np.isclose(df.amplitude, float(scenario["amplitude"]))
    m &= np.isclose(df.shared_fraction, float(scenario.get("shared_fraction", 0.0)))
    m &= np.isclose(df.threshold_sd, float(scenario.get("threshold_sd", 0.0)))
    if "common_threshold" in scenario:
        m &= np.isclose(df.common_threshold, float(scenario["common_threshold"]), equal_nan=False)
    else:
        m &= df.common_threshold.isna()
    return m


def shared_params_and_labels(
    geo: parent.ClimateGeometry,
    block: str,
    scenario: dict,
    replicate: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact reconstruction of parent labels plus privileged true axis/thresholds."""
    rng = np.random.default_rng(parent.seed_for("labels", block, "evaluation", parent.scenario_id(scenario), replicate))
    amplitude = float(scenario["amplitude"])
    fraction = float(scenario.get("shared_fraction", 0.0))
    if fraction != 1.0 or not scenario["world"].startswith("climate_shared"):
        raise ValueError("diagnostic is frozen to full-shared climate worlds")
    x = geo.climate[block]

    ne = base.unit_vectors(rng, (369, 4))
    _ng = base.unit_vectors(rng, (369, 3))  # consumed exactly as in parent generator
    common = base.unit_vectors(rng, (1, 4))[0]
    chosen = rng.permutation(369)[: round(369 * fraction)]
    ne[chosen] = common
    projection = np.einsum("nd,nd->n", x, ne[geo.sid])

    if "common_threshold" in scenario:
        threshold = np.full(369, float(scenario["common_threshold"]), dtype=float)
    elif float(scenario.get("threshold_sd", 0.0)) > 0:
        threshold = rng.normal(0.0, float(scenario["threshold_sd"]), 369)
    else:
        threshold = np.zeros(369, dtype=float)
    signs = rng.choice([-1.0, 1.0], 369)
    latent = amplitude * signs[geo.sid] * np.tanh((projection - threshold[geo.sid]) / 0.5)
    labels = latent + rng.normal(size=len(geo.df)) > 0

    reference = parent.labels_for(geo, block, "evaluation", scenario, replicate)
    if not np.array_equal(labels, reference):
        raise RuntimeError("parent synthetic-label reconstruction mismatch")
    return labels, common, threshold


def likelihood_features_axes(x: np.ndarray, axes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    axes = np.asarray(axes, dtype=float)
    if axes.ndim != 2 or axes.shape[1] != 4 or not np.isfinite(axes).all():
        raise ValueError("finite Kx4 axes required")
    norm = np.linalg.norm(axes, axis=1)
    if np.any(norm <= 1e-12):
        raise ValueError("zero diagnostic axis")
    axes = axes / norm[:, None]
    projection = np.einsum("fsnd,kd->fskn", x, axes, optimize=True)
    shifted = projection[:, :, :, None, None, :] - parent.OFFSETS[None, None, None, None, :, None]
    z = shifted * parent.SLOPES[None, None, None, :, None, None]
    zero = (-np.logaddexp(0.0, z)).sum(axis=-1)
    sum_z = z.sum(axis=-1)
    return z, zero, sum_z


def candidate_ll(y: np.ndarray, features: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    _baseline, candidate = parent.axis_log_likelihood_pair(y, features)
    return candidate


def score_generic(ll: np.ndarray, privileged_axis_index: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parent score algebra generalized from 32 to arbitrary candidate-axis count."""
    if ll.ndim != 4 or not np.isfinite(ll).all():
        raise ValueError("finite world,fold,species,axis log likelihoods required")
    k = ll.shape[-1]
    lr = ll - logsumexp(ll, axis=-1, keepdims=True) + np.log(k)
    tr, te = lr[..., : parent.N_TRAIN, :], lr[..., parent.N_TRAIN :, :]

    lf = np.full(len(parent.FRACTIONS), -np.inf)
    lnf = np.full(len(parent.FRACTIONS), -np.inf)
    np.log(parent.FRACTIONS, out=lf, where=parent.FRACTIONS > 0)
    np.log(1.0 - parent.FRACTIONS, out=lnf, where=parent.FRACTIONS < 1)
    mix = np.logaddexp(lnf[:, None], lf[:, None] + tr[..., :, None, :])
    lp = mix.sum(axis=-3)
    posterior = np.exp(lp - logsumexp(lp, axis=(-2, -1), keepdims=True))

    weighted_axis = np.einsum("...fk,f->...k", posterior, parent.FRACTIONS)
    independent = np.einsum("...fk,f->...", posterior, 1.0 - parent.FRACTIONS)
    ratio = independent[..., None] + np.einsum("...k,...tk->...t", weighted_axis, np.exp(te))
    gain = np.log(np.maximum(ratio, np.finfo(float).tiny)).mean(axis=-1) / parent.N_PHOTOS
    mean_f = np.einsum("...fk,f->...", posterior, parent.FRACTIONS)
    if privileged_axis_index is None:
        axis_mass = np.full(gain.shape, np.nan)
    else:
        axis_mass = posterior[..., :, privileged_axis_index].sum(axis=-1)
    return gain, mean_f, axis_mass


def max_projection_corr(x_fold: np.ndarray, true_axis: np.ndarray) -> float:
    flat = x_fold.reshape(-1, 4)
    true_p = flat @ true_axis
    cand_p = flat @ parent.AXES.T
    true_c = true_p - true_p.mean()
    cand_c = cand_p - cand_p.mean(axis=0, keepdims=True)
    true_sd = float(np.sqrt(np.mean(true_c * true_c)))
    cand_sd = np.sqrt(np.mean(cand_c * cand_c, axis=0))
    good = cand_sd > 1e-12
    if true_sd <= 1e-12 or not np.any(good):
        return 0.0
    corr = np.zeros(cand_p.shape[1], dtype=float)
    corr[good] = np.mean(cand_c[:, good] * true_c[:, None], axis=0) / (cand_sd[good] * true_sd)
    return float(np.max(np.abs(corr)))


def species_metrics(
    geo: parent.ClimateGeometry,
    idx: np.ndarray,
    x: np.ndarray,
    labels_selected: np.ndarray,
    true_axis: np.ndarray,
    threshold: np.ndarray,
) -> dict[str, float]:
    true_projection = np.einsum("fsnd,d->fsn", x, true_axis)
    sids = geo.sid[idx[:, :, 0]]
    thresholds = threshold[sids]
    straddle = (true_projection.min(axis=-1) <= thresholds) & (true_projection.max(axis=-1) >= thresholds)
    variable = labels_selected.any(axis=-1) & (~labels_selected).any(axis=-1)
    proj_sd = np.std(true_projection, axis=-1, ddof=0)
    result = {}
    for role, sl in (("train", slice(0, parent.N_TRAIN)), ("evaluation", slice(parent.N_TRAIN, None))):
        result[f"{role}_transition_exposure"] = float(straddle[:, sl].mean())
        result[f"{role}_colour_variability"] = float(variable[:, sl].mean())
        result[f"{role}_true_projection_sd"] = float(proj_sd[:, sl].mean())
    return result


def lookup_parent_row(parent_scores: pd.DataFrame, block: str, scenario: dict, replicate: int) -> pd.Series:
    q = parent_scores[
        (parent_scores.stage == "evaluation")
        & (parent_scores.block == block)
        & (parent_scores.model == PARENT_CANDIDATE)
        & (parent_scores.replicate == replicate)
        & scenario_mask(parent_scores, scenario)
    ]
    if len(q) != 1:
        raise RuntimeError(f"parent score lookup failed: {block} {scenario} {replicate}: {len(q)}")
    return q.iloc[0]


def run_shard(
    geometry_path: Path,
    climate_path: Path,
    parent_scores_path: Path,
    block: str,
    start: int,
    stop: int,
    output: Path,
) -> None:
    if block not in parent.BLOCKS or not (0 <= start < stop <= 250):
        raise ValueError("invalid diagnostic shard")
    if output.exists():
        raise FileExistsError("use a fresh diagnostic output")
    geo = parent.ClimateGeometry(geometry_path, climate_path)
    parent_scores = pd.read_csv(parent_scores_path)
    if len(parent_scores) != 45000:
        raise RuntimeError("parent full-score census drift")

    records: list[dict[str, object]] = []
    max_discrepancy = 0.0
    for replicate in range(start, stop):
        idx = geo.schedule("evaluation", replicate)
        x = geo.climate[block][idx]
        fixed_features = parent.likelihood_features(x)
        for scenario in SCENARIOS:
            labels, true_axis, threshold = shared_params_and_labels(geo, block, scenario, replicate)
            y_selected = labels[idx]
            y = y_selected[None, ...]

            _bll, fixed_ll = parent.axis_log_likelihood_pair(y, fixed_features)
            fixed_score_fold, _fixed_f = parent.score_from_ll(fixed_ll)
            fixed_score = float(np.mean(fixed_score_fold[0]))
            parent_row = lookup_parent_row(parent_scores, block, scenario, replicate)
            discrepancy = abs(fixed_score - float(parent_row.score))
            fold_discrepancy = max(
                abs(float(fixed_score_fold[0, f]) - float(parent_row[f"score_fold{f}"])) for f in range(4)
            )
            max_discrepancy = max(max_discrepancy, discrepancy, fold_discrepancy)
            if discrepancy > 2e-12 or fold_discrepancy > 2e-12:
                raise RuntimeError(f"parent score reconstruction mismatch: {discrepancy} {fold_discrepancy}")

            true_features = likelihood_features_axes(x, true_axis[None, :])
            true_ll = candidate_ll(y, true_features)
            augmented_ll = np.concatenate([fixed_ll, true_ll], axis=-1)
            augmented_score_fold, augmented_f, true_axis_mass_fold = score_generic(
                augmented_ll, privileged_axis_index=augmented_ll.shape[-1] - 1
            )
            augmented_score = float(np.mean(augmented_score_fold[0]))

            vector_coverage = float(np.max(np.abs(parent.AXES @ true_axis)))
            corr_by_fold = [max_projection_corr(x[f], true_axis) for f in range(4)]
            metrics = species_metrics(geo, idx, x, y_selected, true_axis, threshold)
            record: dict[str, object] = {
                "block": block,
                "replicate": replicate,
                "world": scenario["world"],
                "amplitude": float(scenario["amplitude"]),
                "shared_fraction": float(scenario["shared_fraction"]),
                "threshold_sd": float(scenario.get("threshold_sd", 0.0)),
                "common_threshold": float(scenario.get("common_threshold", np.nan)),
                "parent_score": fixed_score,
                "parent_detected": bool(fixed_score > PARENT_SCORE_THRESHOLD),
                "augmented_true_axis_score": augmented_score,
                "augmented_minus_parent_score": augmented_score - fixed_score,
                "augmented_mean_sharing_fraction": float(np.mean(augmented_f[0])),
                "privileged_true_axis_posterior_mass": float(np.mean(true_axis_mass_fold[0])),
                "true_axis_vector_coverage": vector_coverage,
                "true_axis_projection_corr_mean": float(np.mean(corr_by_fold)),
                "true_axis_projection_corr_min_fold": float(np.min(corr_by_fold)),
                "parent_score_reconstruction_abs_error": discrepancy,
            }
            record.update({f"true_axis_projection_corr_fold{f}": corr_by_fold[f] for f in range(4)})
            record.update(metrics)
            records.append(record)
        print(json.dumps({"block": block, "completed_replicate": replicate + 1, "stop": stop}), flush=True)

    output.mkdir(parents=True)
    table = pd.DataFrame(records)
    table.to_csv(output / "diagnostic_rows.csv", index=False, lineterminator="\n")
    meta = {
        "protocol": "hypervolume-real-climate-identifiability-diagnostic-v1",
        "specification_commit": SPECIFICATION_COMMIT,
        "block": block,
        "replicate_start": start,
        "replicate_stop": stop,
        "rows": len(table),
        "maximum_parent_score_reconstruction_abs_error": max_discrepancy,
        "geometry_sha256": digest(geometry_path),
        "real_climate_export_sha256": digest(climate_path),
        "parent_scores_sha256": digest(parent_scores_path),
        "runner_sha256": digest(Path(__file__)),
        "rows_sha256": digest(output / "diagnostic_rows.csv"),
        "biological_colour_values_read": False,
        "parent_failure_reclassified": False,
        "privileged_true_axis_used_for_inference": False,
    }
    (output / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--climate", type=Path, required=True)
    ap.add_argument("--parent-scores", type=Path, required=True)
    ap.add_argument("--block", choices=list(parent.BLOCKS), required=True)
    ap.add_argument("--rep-start", type=int, required=True)
    ap.add_argument("--rep-stop", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    run_shard(args.geometry, args.climate, args.parent_scores, args.block, args.rep_start, args.rep_stop, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
