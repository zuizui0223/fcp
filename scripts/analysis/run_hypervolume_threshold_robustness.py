#!/usr/bin/env python3
"""Threshold-shift and richer-null robustness on frozen actual coordinate support.

The final specification was frozen at commit 0e9b50009b4706c5c6ad5360988e1c9953487ee1
before any synthetic outcome from this stage existed. This runner never reads
biological flower-colour values or real climate rasters.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import log_ndtr, logsumexp

import run_hypervolume_actual_support_blocked as support
import run_hypervolume_pooled_response as pooled
import run_hypervolume_transfer_benchmark as base

SEED = 2026090705
REPS = 250
SPEC = "0e9b50009b4706c5c6ad5360988e1c9953487ee1"
K = pooled.K
AXES = pooled.AXES
FRACTIONS = pooled.FRACTIONS
OFFSETS = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)
SLOPES = {
    "linear_logit": np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0]),
    "saturating_probit": np.array([0.0, 0.5, 1.0, 2.0]),
}
MODELS = (
    "linear_logit_joint_train",
    "linear_logit_joint_train_intercept_marginal",
    "saturating_probit_joint_train_intercept_marginal",
)
PARENT_HASHES = {
    "run_hypervolume_pooled_response.py": "e87f8343fdd84f64c7245dff0399fb6b33ec8040f5290ffe5a8fdf1a1158581d",
    "run_hypervolume_actual_support_blocked.py": "dbedd48765816f6a2c8678585c62c52add32bdbd75a3fab59a66f06d8084891e",
    "run_hypervolume_transfer_benchmark.py": "b5312f75128992f99b45898a6e5190aab1a86e18ad59573d26ca5fae66d8219b",
}

NUISANCE = (
    dict(world="no_structure", amplitude=0.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="environmental_independent", amplitude=0.5, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="environmental_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="environmental_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="geographic_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="mixed_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="environmental_independent_shifted", amplitude=0.5, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="environmental_independent_shifted", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="environmental_independent_shifted", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="support_aligned_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="directionally_clustered_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5, axis_concentration=2.0),
    dict(world="nonlinear_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5),
)
POSITIVE = (
    dict(world="environmental_shared_shifted", amplitude=1.0, shared_fraction=0.25, threshold_sd=0.5),
    dict(world="environmental_shared_shifted", amplitude=1.0, shared_fraction=0.5, threshold_sd=0.5),
    dict(world="environmental_shared_shifted", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.5),
    dict(world="environmental_shared_shifted", amplitude=1.0, shared_fraction=1.0, threshold_sd=1.0),
    dict(world="environmental_shared_shifted", amplitude=2.0, shared_fraction=1.0, threshold_sd=0.5),
    dict(world="environmental_shared_common_offset", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.0, common_threshold=0.5),
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scenario_id(s: dict) -> str:
    return "|".join(f"{k}={s[k]}" for k in sorted(s))


def seed_for(*parts: object) -> int:
    text = "|".join(map(str, (SEED,) + parts))
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def schedule(geo: support.Geometry, stage: str, replicate: int) -> np.ndarray:
    rng = np.random.default_rng(seed_for("schedule", stage, replicate))
    priority = rng.random(len(geo.df))
    result = np.empty((4, 40, 20), dtype=np.int32)
    for fold in range(4):
        chosen = np.r_[rng.choice(geo.training_eligible[fold], 20, replace=False),
                       rng.choice(geo.test_eligible[fold], 20, replace=False)]
        for j, sid in enumerate(chosen):
            pool = geo.pools[(1, fold, int(sid))]
            result[fold, j] = pool[np.argsort(priority[pool], kind="stable")[:20]]
    return result


def support_aligned_axes(geo: support.Geometry) -> np.ndarray:
    total = np.zeros((369, 2), dtype=float)
    count = np.zeros(369, dtype=float)
    np.add.at(total, geo.sid, geo.e)
    np.add.at(count, geo.sid, 1.0)
    mean = total / np.maximum(count[:, None], 1.0)
    norm = np.linalg.norm(mean, axis=1)
    axes = np.zeros_like(mean)
    good = norm > 1e-12
    axes[good] = mean[good] / norm[good, None]
    # Deterministic fallback for the vanishingly small mean-support case.
    angle = (np.arange(369) * np.pi * (3.0 - np.sqrt(5.0))) % (2.0 * np.pi)
    axes[~good] = np.stack([np.cos(angle[~good]), np.sin(angle[~good])], axis=1)
    return axes


def labels_for(geo: support.Geometry, stage: str, scenario: dict, replicate: int,
               aligned_axes: np.ndarray) -> np.ndarray:
    """Generate one complete synthetic world; thresholds/directions are species-level."""
    rng = np.random.default_rng(seed_for("labels", stage, scenario_id(scenario), replicate))
    world = scenario["world"]
    amplitude = float(scenario["amplitude"])
    fraction = float(scenario.get("shared_fraction", 0.0))
    threshold_sd = float(scenario.get("threshold_sd", 0.0))

    ne = base.unit_vectors(rng, (369, 2))
    ng = base.unit_vectors(rng, (369, 3))
    nonlinear = False

    if world.startswith("environmental_shared"):
        common = base.unit_vectors(rng, (1, 2))[0]
        chosen = rng.permutation(369)[:round(369 * fraction)]
        ne[chosen] = common
    elif world == "support_aligned_independent":
        ne = aligned_axes.copy()
    elif world == "directionally_clustered_independent":
        mu = rng.uniform(-np.pi, np.pi)
        kappa = float(scenario.get("axis_concentration", 2.0))
        theta = rng.vonmises(mu, kappa, 369)
        ne = np.stack([np.cos(theta), np.sin(theta)], axis=1)
    elif world == "nonlinear_independent":
        nonlinear = True
    elif world not in {
        "no_structure", "environmental_independent", "environmental_independent_shifted",
        "geographic_independent", "mixed_independent",
    }:
        raise ValueError("unknown scenario " + world)

    pe = np.einsum("nd,nd->n", geo.e, ne[geo.sid])
    pg = np.einsum("nd,nd->n", geo.g, ng[geo.sid])
    if nonlinear:
        perp = np.stack([-ne[:, 1], ne[:, 0]], axis=1)
        qe = np.einsum("nd,nd->n", geo.e, perp[geo.sid])
        projection = pe + 0.75 * (qe * qe - pe * pe)
    elif world == "geographic_independent":
        projection = pg
    elif world == "mixed_independent":
        use_env = rng.random(369) < 0.5
        projection = np.where(use_env[geo.sid], pe, pg)
    elif world == "no_structure":
        projection = np.zeros(len(geo.df), dtype=float)
    else:
        projection = pe

    if "common_threshold" in scenario:
        threshold = np.full(369, float(scenario["common_threshold"]))
    elif threshold_sd > 0:
        threshold = rng.normal(0.0, threshold_sd, 369)
    else:
        threshold = np.zeros(369)
    signs = rng.choice([-1.0, 1.0], 369)
    latent = amplitude * signs[geo.sid] * np.tanh((projection - threshold[geo.sid]) / 0.5)
    return latent + rng.normal(size=len(geo.df)) > 0


def axis_log_likelihood_intercept(y: np.ndarray, x: np.ndarray, link: str) -> np.ndarray:
    """Integrate sign, slope and species-specific threshold for every candidate axis.

    y shape (fold,species,photo), x shape (fold,species,photo,2).
    Returns (fold,species,axis).
    """
    if link not in SLOPES:
        raise ValueError("unknown link")
    if y.ndim != 3 or x.shape != y.shape + (2,):
        raise ValueError("label/coordinate shape mismatch")
    projection = np.einsum("fsnd,kd->fskn", x, AXES)
    shifted = projection[:, :, :, None, None, :] - OFFSETS[None, None, None, None, :, None]
    slope = SLOPES[link][None, None, None, :, None, None]
    if link == "saturating_probit":
        latent = np.tanh(shifted / 0.5) * slope
        lp, lq = log_ndtr(latent), log_ndtr(-latent)
    else:
        latent = shifted * slope
        lp = -np.logaddexp(0.0, -latent)
        lq = -np.logaddexp(0.0, latent)
    yy = y[:, :, None, None, None, :].astype(float)
    pos = (yy * lp + (1.0 - yy) * lq).sum(axis=-1)
    neg = (yy * lq + (1.0 - yy) * lp).sum(axis=-1)
    sign_integrated = np.logaddexp(pos, neg) - np.log(2.0)
    return logsumexp(sign_integrated, axis=(-2, -1)) - np.log(len(OFFSETS) * len(SLOPES[link]))


def baseline_score(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    features = pooled.likelihood_features(x, "linear_logit")
    ll = pooled.axis_log_likelihood(y[None, ...], features)
    lr = pooled.normalized_log_evidence(ll)
    scored, mean_f = pooled.marginal_scores(lr)
    return scored["joint_train"][0], mean_f[0]


def intercept_score(y: np.ndarray, x: np.ndarray, link: str) -> tuple[np.ndarray, np.ndarray]:
    ll = axis_log_likelihood_intercept(y, x, link)
    lr = pooled.normalized_log_evidence(ll[None, ...])
    scored, mean_f = pooled.marginal_scores(lr)
    return scored["joint_train"][0], mean_f[0]


def run_shard(geometry: Path, output: Path, stage: str, start: int, stop: int) -> None:
    if output.exists():
        raise FileExistsError("use a fresh shard output")
    if stage not in {"calibration", "evaluation"} or not (0 <= start < stop <= REPS):
        raise ValueError("invalid shard")
    for name, expected in PARENT_HASHES.items():
        path = Path(__file__).with_name(name)
        if digest(path) != expected:
            raise RuntimeError("parent code drift: " + name)
    geo = support.Geometry(geometry)
    aligned = support_aligned_axes(geo)
    scenarios = NUISANCE + (POSITIVE if stage == "evaluation" else ())
    records = []
    for r in range(start, stop):
        idx = schedule(geo, stage, r)
        x = geo.e[idx]
        for s in scenarios:
            y = labels_for(geo, stage, s, r, aligned)[idx]
            bscore, bf = baseline_score(y, x)
            lscore, lf = intercept_score(y, x, "linear_logit")
            pscore, pf = intercept_score(y, x, "saturating_probit")
            model_values = {
                "linear_logit_joint_train": (bscore, bf),
                "linear_logit_joint_train_intercept_marginal": (lscore, lf),
                "saturating_probit_joint_train_intercept_marginal": (pscore, pf),
            }
            for model, (score, mean_f) in model_values.items():
                row = dict(stage=stage, replicate=r, model=model,
                           world=s["world"], amplitude=float(s["amplitude"]),
                           shared_fraction=float(s.get("shared_fraction", 0.0)),
                           threshold_sd=float(s.get("threshold_sd", 0.0)),
                           common_threshold=float(s.get("common_threshold", np.nan)),
                           axis_concentration=float(s.get("axis_concentration", np.nan)),
                           score=float(np.mean(score)), train_fraction=float(np.mean(mean_f)))
                row.update({f"score_fold{k}": float(score[k]) for k in range(4)})
                records.append(row)
        print(json.dumps({"stage": stage, "replicate": r + 1, "stop": stop, "rows": len(records)}), flush=True)
    output.mkdir(parents=True)
    table = pd.DataFrame(records)
    table.to_csv(output / "scores.csv", index=False, lineterminator="\n")
    meta = {
        "protocol": "hypervolume-threshold-robustness-v1",
        "specification_commit": SPEC,
        "stage": stage,
        "replicate_start": start,
        "replicate_stop": stop,
        "score_rows": len(table),
        "models": list(MODELS),
        "biological_colour_values_read": False,
        "real_climate_used": False,
        "parent_decisions_modified": False,
        "runner_sha256": digest(Path(__file__)),
        "geometry_sha256": digest(geometry),
        "scores_sha256": digest(output / "scores.csv"),
    }
    (output / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--stage", choices=["calibration", "evaluation"], required=True)
    ap.add_argument("--rep-start", type=int, required=True)
    ap.add_argument("--rep-stop", type=int, required=True)
    args = ap.parse_args()
    run_shard(args.geometry, args.output, args.stage, args.rep_start, args.rep_stop)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
