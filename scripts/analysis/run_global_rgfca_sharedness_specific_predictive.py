#!/usr/bin/env python3
"""Species-disjoint predictive test of cross-species geographic boundary sharing.

The structured null explicitly allows every species to possess its own strong
spatial boundary. No biological flower-colour value is read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

import run_hypervolume_actual_support_blocked as support
import run_hypervolume_transfer_benchmark as base

SPECIFICATION_COMMIT = "34fe3093cd17c6d20fb8044259984314761cd210"
SEED = 2026090708
AXIS_SEED = 202609070801
REPS = 250
N_TRAIN = 100
N_TEST = 100
N_PHOTOS = 20
GEOMETRY_SHA = "5fa0c45a3f54f7441e2e372fa4ff4e813548c4cc44896a12a5ba8e97fb5af743"
FRACTIONS = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 1.0], dtype=float)
SLOPES = np.array([0.0, 1.0, 2.0, 4.0, 8.0], dtype=float)
OFFSETS = np.array([-0.5, -0.25, 0.0, 0.25, 0.5], dtype=float)

NUISANCE = (
    dict(world="no_structure", amplitude=0.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="independent_boundaries", amplitude=0.5, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="independent_boundaries", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="independent_boundaries", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="independent_shifted_boundaries", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.25),
    dict(world="independent_shifted_boundaries", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.25),
    dict(world="directionally_clustered_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.25, axis_concentration=2.0),
    dict(world="directionally_clustered_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.25, axis_concentration=2.0),
    dict(world="bimodal_clustered_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.25, axis_concentration=3.0),
)
POSITIVE = (
    dict(world="partially_shared_boundaries", amplitude=0.5, shared_fraction=0.5, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=1.0, shared_fraction=0.1, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=1.0, shared_fraction=0.25, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=1.0, shared_fraction=0.5, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=1.0, shared_fraction=0.75, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=2.0, shared_fraction=0.25, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=2.0, shared_fraction=0.5, threshold_sd=0.25),
    dict(world="partially_shared_boundaries", amplitude=2.0, shared_fraction=1.0, threshold_sd=0.25),
    dict(world="shared_common_offset", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.0, common_threshold=0.25),
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_for(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("|".join(map(str, (SEED,) + parts)).encode()).digest()[:8], "little"
    )


def scenario_id(s: dict) -> str:
    return "|".join(f"{k}={s[k]}" for k in sorted(s))


def canonical_axis(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).copy()
    for j in range(v.shape[-1]):
        flip = np.isclose(v[..., :j], 0.0).all(axis=-1) & (v[..., j] < 0) if j else (v[..., j] < 0)
        v[flip] *= -1.0
        decided = ~np.isclose(v[..., j], 0.0)
        if np.all(decided):
            break
    return v


def fixed_axes() -> np.ndarray:
    rng = np.random.default_rng(AXIS_SEED)
    axes = base.unit_vectors(rng, (96, 3))
    axes = canonical_axis(axes)
    if axes.shape != (96, 3) or not np.allclose(np.linalg.norm(axes, axis=1), 1.0, atol=1e-12):
        raise RuntimeError("fixed geographic axis grid failed")
    return axes


AXES = fixed_axes()
K = len(AXES)


class Geometry:
    def __init__(self, path: Path):
        if digest(path) != GEOMETRY_SHA:
            raise ValueError("geometry checksum drift")
        parent = support.Geometry(path)
        self.df = parent.df.copy()
        self.species = parent.species
        self.sid = parent.sid.copy()
        self.xyz = parent.xyz.copy()
        self.mask = parent.mask.copy()
        self.train_sid = set(parent.train_sid)
        self.test_sid = set(parent.test_sid)
        if len(self.species) != 369 or int(self.mask.sum()) != 21424:
            raise ValueError("eligible geometry drift")
        self.pools = {
            i: np.flatnonzero((self.sid == i) & self.mask).astype(np.int32)
            for i in range(369)
        }
        if min(map(len, self.pools.values())) < 40:
            raise ValueError("frozen >=40 classifiable gate drift")
        self.train_ids = np.array(sorted(self.train_sid), dtype=np.int32)
        self.test_ids = np.array(sorted(self.test_sid), dtype=np.int32)
        if len(self.train_ids) != 184 or len(self.test_ids) != 185 or set(self.train_ids) & set(self.test_ids):
            raise ValueError("fixed species split drift")

    def schedule(self, stage: str, replicate: int) -> np.ndarray:
        rng = np.random.default_rng(seed_for("schedule", stage, replicate))
        chosen = np.r_[
            rng.choice(self.train_ids, N_TRAIN, replace=False),
            rng.choice(self.test_ids, N_TEST, replace=False),
        ]
        priority = rng.random(len(self.df))
        out = np.empty((N_TRAIN + N_TEST, N_PHOTOS), dtype=np.int32)
        for j, sid in enumerate(chosen):
            pool = self.pools[int(sid)]
            out[j] = pool[np.argsort(priority[pool], kind="stable")[:N_PHOTOS]]
        return out


def clustered_axes(rng: np.random.Generator, n: int, concentration: float) -> np.ndarray:
    center = base.unit_vectors(rng, (1, 3))[0]
    raw = concentration * center[None, :] + rng.normal(size=(n, 3))
    norm = np.linalg.norm(raw, axis=1)
    if np.any(norm <= 1e-12):
        raise RuntimeError("clustered geographic axis failure")
    return raw / norm[:, None]


def bimodal_axes(rng: np.random.Generator, n: int, concentration: float) -> np.ndarray:
    c1 = base.unit_vectors(rng, (1, 3))[0]
    raw2 = rng.normal(size=3)
    raw2 -= raw2.dot(c1) * c1
    if np.linalg.norm(raw2) <= 1e-12:
        raw2 = np.roll(c1, 1) - np.roll(c1, 1).dot(c1) * c1
    c2 = raw2 / np.linalg.norm(raw2)
    group = rng.integers(0, 2, size=n)
    center = np.where(group[:, None] == 0, c1[None, :], c2[None, :])
    raw = concentration * center + rng.normal(size=(n, 3))
    norm = np.linalg.norm(raw, axis=1)
    if np.any(norm <= 1e-12):
        raise RuntimeError("bimodal geographic axis failure")
    return raw / norm[:, None]


def generator_params(geo: Geometry, stage: str, scenario: dict, replicate: int):
    rng = np.random.default_rng(seed_for("labels", stage, scenario_id(scenario), replicate))
    world = scenario["world"]
    normals = base.unit_vectors(rng, (369, 3))
    if world == "directionally_clustered_independent":
        normals = clustered_axes(rng, 369, float(scenario["axis_concentration"]))
    elif world == "bimodal_clustered_independent":
        normals = bimodal_axes(rng, 369, float(scenario["axis_concentration"]))
    elif world in {"partially_shared_boundaries", "shared_common_offset"}:
        common = base.unit_vectors(rng, (1, 3))[0]
        chosen = rng.permutation(369)[: round(369 * float(scenario["shared_fraction"]))]
        normals[chosen] = common
    elif world not in {
        "no_structure", "independent_boundaries", "independent_shifted_boundaries"
    }:
        raise ValueError("unknown generator world " + world)

    if "common_threshold" in scenario:
        threshold = np.full(369, float(scenario["common_threshold"]), dtype=float)
    elif float(scenario.get("threshold_sd", 0.0)) > 0:
        threshold = rng.normal(0.0, float(scenario["threshold_sd"]), 369)
    else:
        threshold = np.zeros(369, dtype=float)
    sign = rng.choice([-1.0, 1.0], 369)
    return rng, normals, threshold, sign


def labels_for(geo: Geometry, stage: str, scenario: dict, replicate: int) -> np.ndarray:
    rng, normals, threshold, sign = generator_params(geo, stage, scenario, replicate)
    if scenario["world"] == "no_structure":
        projection = np.zeros(len(geo.df), dtype=float)
    else:
        projection = np.einsum("nd,nd->n", geo.xyz, normals[geo.sid])
    latent = float(scenario["amplitude"]) * sign[geo.sid] * np.tanh(
        (projection - threshold[geo.sid]) / 0.15
    )
    return latent + rng.normal(size=len(geo.df)) > 0


def precompute_features(x: np.ndarray):
    if x.shape != (N_TRAIN + N_TEST, N_PHOTOS, 3) or not np.isfinite(x).all():
        raise ValueError("expected finite 200x20x3 selected coordinates")
    projection = np.einsum("snd,kd->skn", x, AXES, optimize=True)
    sum_projection = projection.sum(axis=-1)
    zero = np.empty((len(x), K, len(SLOPES), len(OFFSETS)), dtype=float)
    sum_z = np.empty_like(zero)
    for li, slope in enumerate(SLOPES):
        for oi, offset in enumerate(OFFSETS):
            z = slope * (projection - offset)
            zero[:, :, li, oi] = (-np.logaddexp(0.0, z)).sum(axis=-1)
            sum_z[:, :, li, oi] = slope * (sum_projection - N_PHOTOS * offset)
    return projection, zero, sum_z


def axis_log_likelihood(y: np.ndarray, features) -> np.ndarray:
    """Return world x species x axis log likelihood marginalized over sign/slope/offset."""
    projection, zero, sum_z = features
    if y.ndim != 3 or y.shape[1:] != (N_TRAIN + N_TEST, N_PHOTOS):
        raise ValueError("unexpected synthetic-label batch")
    yp = np.einsum("wsn,skn->wsk", y.astype(float), projection, optimize=True)
    n1 = y.sum(axis=-1).astype(float)
    signed_terms = []
    for li, slope in enumerate(SLOPES):
        offset_terms = []
        for oi, offset in enumerate(OFFSETS):
            yz = slope * (yp - n1[..., None] * offset)
            pos = zero[None, :, :, li, oi] + yz
            neg = zero[None, :, :, li, oi] + sum_z[None, :, :, li, oi] - yz
            offset_terms.append(np.logaddexp(pos, neg) - np.log(2.0))
        signed_terms.append(np.stack(offset_terms, axis=-1))
    stacked = np.stack(signed_terms, axis=-2)
    return logsumexp(stacked, axis=(-2, -1)) - np.log(len(SLOPES) * len(OFFSETS))


def score_from_ll(ll: np.ndarray):
    if ll.ndim != 3 or ll.shape[1:] != (N_TRAIN + N_TEST, K) or not np.isfinite(ll).all():
        raise ValueError("finite world x 200 species x K likelihood required")
    lr = ll - logsumexp(ll, axis=-1, keepdims=True) + np.log(K)
    tr, te = lr[:, :N_TRAIN], lr[:, N_TRAIN:]
    lf = np.full(len(FRACTIONS), -np.inf)
    lnf = np.full(len(FRACTIONS), -np.inf)
    np.log(FRACTIONS, out=lf, where=FRACTIONS > 0)
    np.log(1.0 - FRACTIONS, out=lnf, where=FRACTIONS < 1)
    mix = np.logaddexp(lnf[:, None, None], lf[:, None, None] + tr[:, :, None, :])
    lp = mix.sum(axis=1)
    posterior = np.exp(lp - logsumexp(lp, axis=(-2, -1), keepdims=True))
    weighted_axis = np.einsum("wfk,f->wk", posterior, FRACTIONS)
    independent = np.einsum("wfk,f->w", posterior, 1.0 - FRACTIONS)
    ratio = independent[:, None] + np.einsum("wk,wtk->wt", weighted_axis, np.exp(te))
    score = np.log(np.maximum(ratio, np.finfo(float).tiny)).mean(axis=1) / N_PHOTOS
    mean_f = np.einsum("wfk,f->w", posterior, FRACTIONS)
    axis_mass = posterior.sum(axis=1)
    entropy = -np.sum(axis_mass * np.log(np.maximum(axis_mass, np.finfo(float).tiny)), axis=1) / np.log(K)
    max_axis_mass = axis_mass.max(axis=1)
    return score, mean_f, entropy, max_axis_mass


def run_shard(geometry_path: Path, output: Path, stage: str, start: int, stop: int) -> None:
    if stage not in {"calibration", "evaluation"} or not (0 <= start < stop <= REPS):
        raise ValueError("invalid shard")
    if output.exists():
        raise FileExistsError("use a fresh output")
    geo = Geometry(geometry_path)
    scenarios = NUISANCE + (POSITIVE if stage == "evaluation" else ())
    records = []
    inclusion_train = np.zeros(369, dtype=int)
    inclusion_eval = np.zeros(369, dtype=int)
    photo_counts = np.zeros(len(geo.df), dtype=int)

    for replicate in range(start, stop):
        idx = geo.schedule(stage, replicate)
        selected_sid = geo.sid[idx[:, 0]]
        inclusion_train[selected_sid[:N_TRAIN]] += 1
        inclusion_eval[selected_sid[N_TRAIN:]] += 1
        np.add.at(photo_counts, idx.ravel(), 1)
        features = precompute_features(geo.xyz[idx])
        labels = np.stack([labels_for(geo, stage, s, replicate)[idx] for s in scenarios])
        ll = axis_log_likelihood(labels, features)
        score, mean_f, entropy, max_axis_mass = score_from_ll(ll)
        for si, scenario in enumerate(scenarios):
            records.append({
                "stage": stage,
                "replicate": replicate,
                "world": scenario["world"],
                "amplitude": float(scenario["amplitude"]),
                "shared_fraction": float(scenario.get("shared_fraction", 0.0)),
                "threshold_sd": float(scenario.get("threshold_sd", 0.0)),
                "common_threshold": float(scenario.get("common_threshold", np.nan)),
                "axis_concentration": float(scenario.get("axis_concentration", np.nan)),
                "score": float(score[si]),
                "posterior_mean_shared_fraction": float(mean_f[si]),
                "posterior_axis_entropy_normalized": float(entropy[si]),
                "posterior_max_axis_mass": float(max_axis_mass[si]),
            })
        print(json.dumps({"stage": stage, "completed_replicate": replicate + 1, "stop": stop}), flush=True)

    output.mkdir(parents=True)
    table = pd.DataFrame(records)
    table.to_csv(output / "scores.csv", index=False, lineterminator="\n")
    pd.DataFrame({
        "species": geo.species,
        "training_inclusions": inclusion_train,
        "evaluation_inclusions": inclusion_eval,
    }).to_csv(output / "species_inclusions.csv", index=False, lineterminator="\n")
    used = np.flatnonzero(photo_counts)
    pd.DataFrame({
        "photo_id": geo.df.photo_id.iloc[used].astype(str).to_numpy(),
        "inclusions": photo_counts[used],
    }).to_csv(output / "photo_inclusions.csv", index=False, lineterminator="\n")
    meta = {
        "protocol": "global-rgfca-sharedness-specific-predictive-v1",
        "specification_commit": SPECIFICATION_COMMIT,
        "stage": stage,
        "replicate_start": start,
        "replicate_stop": stop,
        "score_rows": len(table),
        "geometry_sha256": digest(geometry_path),
        "axis_grid_sha256": hashlib.sha256(AXES.tobytes()).hexdigest(),
        "runner_sha256": digest(Path(__file__)),
        "scores_sha256": digest(output / "scores.csv"),
        "biological_colour_values_read": False,
        "structured_independent_boundaries_in_null": True,
        "parent_G1_reclassified": False,
    }
    (output / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--stage", choices=["calibration", "evaluation"], required=True)
    ap.add_argument("--rep-start", type=int, required=True)
    ap.add_argument("--rep-stop", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    run_shard(args.geometry, args.output, args.stage, args.rep_start, args.rep_stop)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
