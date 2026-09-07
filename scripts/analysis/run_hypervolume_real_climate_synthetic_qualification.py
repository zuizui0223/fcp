#!/usr/bin/env python3
"""Synthetic response-sharing qualification on SHA-locked real CHELSA geometry.

Specification commit da9c7a115bdba3bd7c133bf61f8eaf71c4d387fc precedes every
synthetic outcome from this stage. Biological flower-colour values are forbidden.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.spatial import cKDTree

import run_hypervolume_actual_support_blocked as support
import run_hypervolume_transfer_benchmark as base

SPEC = "da9c7a115bdba3bd7c133bf61f8eaf71c4d387fc"
SEED = 2026090707
AXIS_SEED = 202609070701
REPS = 250
N_TRAIN = 20
N_TEST = 20
N_PHOTOS = 20
BUFFER_KM = 500.0
CLIMATE_SHA = "33ab7c56df9bcd5a1f88ae97346ab7108d4a12cdffb4168359448b308ef5a4eb"
GEOMETRY_SHA = "5fa0c45a3f54f7441e2e372fa4ff4e813548c4cc44896a12a5ba8e97fb5af743"
BLOCKS = {
    "thermal_regime": ("env_bio01", "env_bio04", "env_bio05", "env_bio06"),
    "water_balance": ("env_bio12", "env_bio14", "env_bio15", "env_cmi_mean"),
    "atmospheric_energy_dryness": ("env_vpd_mean", "env_rsds_mean", "env_gdd5", "env_sfcWind_mean"),
}
ALL_ENV = tuple(v for values in BLOCKS.values() for v in values)
FRACTIONS = np.array([0.0, 0.25, 0.5, 0.75, 1.0], dtype=float)
SLOPES = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0], dtype=float)
OFFSETS = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)
MODELS = (
    "linear_logit_joint_train_fixed_zero_threshold",
    "linear_logit_joint_train_intercept_marginal",
)
NUISANCE = (
    dict(world="no_structure", amplitude=0.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="climate_independent", amplitude=0.5, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="climate_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="climate_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="geographic_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="mixed_independent", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.0),
    dict(world="climate_independent_shifted", amplitude=0.5, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="climate_independent_shifted", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="climate_independent_shifted", amplitude=2.0, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="support_aligned_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5),
    dict(world="directionally_clustered_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5, axis_concentration=2.0),
    dict(world="nonlinear_independent", amplitude=1.0, shared_fraction=0.0, threshold_sd=0.5),
)
POSITIVE = (
    dict(world="climate_shared_shifted", amplitude=1.0, shared_fraction=0.25, threshold_sd=0.5),
    dict(world="climate_shared_shifted", amplitude=1.0, shared_fraction=0.5, threshold_sd=0.5),
    dict(world="climate_shared_shifted", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.5),
    dict(world="climate_shared_shifted", amplitude=1.0, shared_fraction=1.0, threshold_sd=1.0),
    dict(world="climate_shared_shifted", amplitude=2.0, shared_fraction=1.0, threshold_sd=0.5),
    dict(world="climate_shared_common_offset", amplitude=1.0, shared_fraction=1.0, threshold_sd=0.0, common_threshold=0.5),
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_for(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("|".join(map(str, (SEED,) + parts)).encode()).digest()[:8],
        "little",
    )


def scenario_id(s: dict) -> str:
    return "|".join(f"{k}={s[k]}" for k in sorted(s))


def fixed_axes() -> np.ndarray:
    rng = np.random.default_rng(AXIS_SEED)
    axes = base.unit_vectors(rng, (32, 4))
    if axes.shape != (32, 4) or not np.allclose(np.linalg.norm(axes, axis=1), 1.0):
        raise RuntimeError("fixed R4 axis grid construction failed")
    return axes


AXES = fixed_axes()
K = len(AXES)


class ClimateGeometry:
    """Exact prior geometry with one joint 12-variable complete photo frame."""

    def __init__(self, geometry_path: Path, climate_path: Path):
        if digest(geometry_path) != GEOMETRY_SHA:
            raise ValueError("geometry checksum drift")
        if digest(climate_path) != CLIMATE_SHA:
            raise ValueError("real-climate export checksum drift")

        raw_g = pd.read_csv(geometry_path)
        raw_c = pd.read_csv(climate_path)
        geom_cols = ["photo_id", "species", "latitude", "longitude", "global_classifiable"]
        if list(raw_c.columns[:5]) != geom_cols or list(raw_g.columns) != geom_cols:
            raise ValueError("geometry/climate leading columns drift")
        if len(raw_g) != 50000 or len(raw_c) != 50000:
            raise ValueError("expected 50000 source rows")
        for col in geom_cols:
            if col in {"latitude", "longitude"}:
                if not np.array_equal(raw_g[col].to_numpy(float), raw_c[col].to_numpy(float)):
                    raise ValueError("coordinate order/value drift: " + col)
            else:
                if not np.array_equal(raw_g[col].astype(str).to_numpy(), raw_c[col].astype(str).to_numpy()):
                    raise ValueError("identity/mask order drift: " + col)
        if tuple(raw_c.columns[5:]) != ALL_ENV:
            raise ValueError("real-climate variable order/family drift")

        base_geo = support.Geometry(geometry_path)
        c = raw_c[raw_c.species.astype(str).isin(base_geo.species)].copy()
        c = c.sort_values(["species", "photo_id"], kind="stable").reset_index(drop=True)
        if not np.array_equal(c.photo_id.astype(str).to_numpy(), base_geo.df.photo_id.astype(str).to_numpy()):
            raise ValueError("eligible climate/geometry photo alignment drift")

        env = c.loc[:, ALL_ENV].to_numpy(float)
        joint = np.isfinite(env).all(axis=1)
        mask = base_geo.mask & joint
        if int(mask.sum()) != 21418:
            raise ValueError(f"joint 12-variable retained frame drift: {int(mask.sum())}")

        self.df = base_geo.df.copy()
        self.species = base_geo.species
        self.sid = base_geo.sid.copy()
        self.xyz = base_geo.xyz.copy()
        self.g = base_geo.g.copy()
        self.sector = base_geo.sector.copy()
        self.train_sid = set(base_geo.train_sid)
        self.test_sid = set(base_geo.test_sid)
        self.mask = mask
        self.climate = {block: c.loc[:, cols].to_numpy(float) for block, cols in BLOCKS.items()}
        if any(not np.isfinite(x[self.mask]).all() for x in self.climate.values()):
            raise ValueError("joint frame contains nonfinite block coordinates")

        self.pools = {}
        self.training_eligible = []
        self.test_eligible = []
        self.capacity = []
        for fold in range(4):
            inside = self.sector == fold
            dist = support.chord_to_km(cKDTree(self.xyz[inside]).query(self.xyz)[0])
            away = (~inside) & (dist >= BUFFER_KM)
            tr, te = [], []
            for i in range(369):
                sm = self.sid == i
                ntr = int(np.count_nonzero(sm & self.mask & away))
                nte = int(np.count_nonzero(sm & self.mask & inside))
                if i in self.train_sid and ntr >= N_PHOTOS:
                    tr.append(i)
                if i in self.test_sid and nte >= N_PHOTOS:
                    te.append(i)
                loc = away if i in self.train_sid else inside
                self.pools[(fold, i)] = np.flatnonzero(sm & self.mask & loc).astype(np.int32)
            if len(tr) < N_TRAIN or len(te) < N_TEST:
                raise ValueError(f"not evaluable fold {fold}: training={len(tr)} evaluation={len(te)}")
            self.training_eligible.append(np.asarray(tr, dtype=np.int32))
            self.test_eligible.append(np.asarray(te, dtype=np.int32))
            usable = away & self.mask & np.isin(self.sid, tr)
            self.capacity.append({
                "fold": fold,
                "training_species": len(tr),
                "evaluation_species": len(te),
                "joint_complete_training_photos": int(np.count_nonzero(usable)),
                "joint_complete_evaluation_photos": int(np.count_nonzero(inside & self.mask & np.isin(self.sid, te))),
                "minimum_training_to_sector_photo_km": float(dist[usable].min()),
            })
        if [x["training_species"] for x in self.capacity] != [133, 135, 141, 165]:
            raise ValueError("joint climate training capacity drift")
        if [x["evaluation_species"] for x in self.capacity] != [60, 62, 59, 26]:
            raise ValueError("joint climate evaluation capacity drift")

        self.support_axes = {}
        for block, x in self.climate.items():
            total = np.zeros((369, 4), dtype=float)
            count = np.zeros(369, dtype=float)
            np.add.at(total, self.sid[self.mask], x[self.mask])
            np.add.at(count, self.sid[self.mask], 1.0)
            mean = total / np.maximum(count[:, None], 1.0)
            norm = np.linalg.norm(mean, axis=1)
            axes = np.zeros_like(mean)
            good = norm > 1e-12
            axes[good] = mean[good] / norm[good, None]
            axes[~good] = AXES[np.arange(np.count_nonzero(~good)) % K]
            self.support_axes[block] = axes

    def schedule(self, stage: str, replicate: int) -> np.ndarray:
        rng = np.random.default_rng(seed_for("schedule", stage, replicate))
        priority = rng.random(len(self.df))
        out = np.empty((4, N_TRAIN + N_TEST, N_PHOTOS), dtype=np.int32)
        for fold in range(4):
            chosen = np.r_[
                rng.choice(self.training_eligible[fold], N_TRAIN, replace=False),
                rng.choice(self.test_eligible[fold], N_TEST, replace=False),
            ]
            for j, sid in enumerate(chosen):
                pool = self.pools[(fold, int(sid))]
                out[fold, j] = pool[np.argsort(priority[pool], kind="stable")[:N_PHOTOS]]
        return out


def clustered_axes(rng: np.random.Generator, n: int, concentration: float) -> np.ndarray:
    center = base.unit_vectors(rng, (1, 4))[0]
    raw = concentration * center[None, :] + rng.normal(size=(n, 4))
    norm = np.linalg.norm(raw, axis=1)
    if np.any(norm <= 1e-12):
        raise RuntimeError("clustered axis normalization failed")
    return raw / norm[:, None]


def perpendicular_axes(rng: np.random.Generator, axes: np.ndarray) -> np.ndarray:
    raw = rng.normal(size=axes.shape)
    raw -= (raw * axes).sum(axis=1)[:, None] * axes
    norm = np.linalg.norm(raw, axis=1)
    bad = norm <= 1e-12
    if np.any(bad):
        raw[bad] = np.roll(axes[bad], 1, axis=1)
        raw[bad] -= (raw[bad] * axes[bad]).sum(axis=1)[:, None] * axes[bad]
        norm = np.linalg.norm(raw, axis=1)
    if np.any(norm <= 1e-12):
        raise RuntimeError("perpendicular axis construction failed")
    return raw / norm[:, None]


def labels_for(geo: ClimateGeometry, block: str, stage: str, scenario: dict, replicate: int) -> np.ndarray:
    rng = np.random.default_rng(seed_for("labels", block, stage, scenario_id(scenario), replicate))
    world = scenario["world"]
    amplitude = float(scenario["amplitude"])
    fraction = float(scenario.get("shared_fraction", 0.0))
    threshold_sd = float(scenario.get("threshold_sd", 0.0))
    x = geo.climate[block]

    ne = base.unit_vectors(rng, (369, 4))
    ng = base.unit_vectors(rng, (369, 3))
    nonlinear = False
    if world.startswith("climate_shared"):
        common = base.unit_vectors(rng, (1, 4))[0]
        chosen = rng.permutation(369)[: round(369 * fraction)]
        ne[chosen] = common
    elif world == "support_aligned_independent":
        ne = geo.support_axes[block].copy()
    elif world == "directionally_clustered_independent":
        ne = clustered_axes(rng, 369, float(scenario.get("axis_concentration", 2.0)))
    elif world == "nonlinear_independent":
        nonlinear = True
    elif world not in {
        "no_structure",
        "climate_independent",
        "climate_independent_shifted",
        "geographic_independent",
        "mixed_independent",
    }:
        raise ValueError("unknown world " + world)

    pc = np.einsum("nd,nd->n", x, ne[geo.sid])
    pg = np.einsum("nd,nd->n", geo.g, ng[geo.sid])
    if nonlinear:
        perp = perpendicular_axes(rng, ne)
        qc = np.einsum("nd,nd->n", x, perp[geo.sid])
        projection = pc + 0.75 * (qc * qc - pc * pc)
    elif world == "geographic_independent":
        projection = pg
    elif world == "mixed_independent":
        use_climate = rng.random(369) < 0.5
        projection = np.where(use_climate[geo.sid], pc, pg)
    elif world == "no_structure":
        projection = np.zeros(len(geo.df), dtype=float)
    else:
        projection = pc

    if "common_threshold" in scenario:
        threshold = np.full(369, float(scenario["common_threshold"]))
    elif threshold_sd > 0:
        threshold = rng.normal(0.0, threshold_sd, 369)
    else:
        threshold = np.zeros(369)
    signs = rng.choice([-1.0, 1.0], 369)
    latent = amplitude * signs[geo.sid] * np.tanh((projection - threshold[geo.sid]) / 0.5)
    return latent + rng.normal(size=len(geo.df)) > 0


def likelihood_features(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute logistic terms once per block/replicate schedule."""
    if x.ndim != 4 or x.shape[-1] != 4 or not np.isfinite(x).all():
        raise ValueError("finite fold,species,photo,R4 coordinates required")
    projection = np.einsum("fsnd,kd->fskn", x, AXES, optimize=True)
    shifted = projection[:, :, :, None, None, :] - OFFSETS[None, None, None, None, :, None]
    z = shifted * SLOPES[None, None, None, :, None, None]
    zero = (-np.logaddexp(0.0, z)).sum(axis=-1)
    sum_z = z.sum(axis=-1)
    return z, zero, sum_z


def axis_log_likelihood_pair(
    y: np.ndarray,
    features: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """One expensive y×feature contraction yields both frozen models."""
    z, zero, sum_z = features
    if y.ndim != 4 or y.shape[1:3] != z.shape[:2] or y.shape[-1] != z.shape[-1]:
        raise ValueError("batch label/feature shape mismatch")
    yz = np.einsum("wfsn,fsklon->wfsklo", y.astype(float), z, optimize=True)
    pos = zero[None, ...] + yz
    neg = zero[None, ...] + sum_z[None, ...] - yz
    signed = np.logaddexp(pos, neg) - np.log(2.0)
    candidate = logsumexp(signed, axis=(-2, -1)) - np.log(len(SLOPES) * len(OFFSETS))
    zero_idx = np.flatnonzero(np.isclose(OFFSETS, 0.0))
    if len(zero_idx) != 1:
        raise RuntimeError("zero offset missing/duplicated")
    baseline = logsumexp(signed[..., int(zero_idx[0])], axis=-1) - np.log(len(SLOPES))
    return baseline, candidate


def normalized_log_evidence(ll: np.ndarray) -> np.ndarray:
    if ll.ndim != 4 or ll.shape[-1] != K or not np.isfinite(ll).all():
        raise ValueError("finite world,fold,species,axis likelihoods required")
    return ll - logsumexp(ll, axis=-1, keepdims=True) + np.log(K)


def train_posterior(train_lr: np.ndarray) -> np.ndarray:
    if train_lr.ndim != 4 or train_lr.shape[-1] != K or train_lr.shape[-2] != N_TRAIN:
        raise ValueError("unexpected training likelihood-ratio shape")
    lf = np.full(len(FRACTIONS), -np.inf)
    lnf = np.full(len(FRACTIONS), -np.inf)
    np.log(FRACTIONS, out=lf, where=FRACTIONS > 0)
    np.log(1.0 - FRACTIONS, out=lnf, where=FRACTIONS < 1)
    mix = np.logaddexp(lnf[:, None], lf[:, None] + train_lr[..., :, None, :])
    lp = mix.sum(axis=-3)
    return np.exp(lp - logsumexp(lp, axis=(-2, -1), keepdims=True))


def score_from_ll(ll: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lr = normalized_log_evidence(ll)
    tr, te = lr[..., :N_TRAIN, :], lr[..., N_TRAIN:, :]
    posterior = train_posterior(tr)
    weighted_axis = np.einsum("...fk,f->...k", posterior, FRACTIONS)
    independent = np.einsum("...fk,f->...", posterior, 1.0 - FRACTIONS)
    ratio = independent[..., None] + np.einsum("...k,...tk->...t", weighted_axis, np.exp(te))
    gain = np.log(np.maximum(ratio, np.finfo(float).tiny)).mean(axis=-1) / N_PHOTOS
    mean_f = np.einsum("...fk,f->...", posterior, FRACTIONS)
    return gain, mean_f


def run_shard(geometry_path: Path, climate_path: Path, block: str, output: Path,
              stage: str, start: int, stop: int) -> None:
    if block not in BLOCKS:
        raise ValueError("unknown climate block")
    if stage not in {"calibration", "evaluation"} or not (0 <= start < stop <= REPS):
        raise ValueError("invalid shard")
    if output.exists():
        raise FileExistsError("use a fresh output")
    geo = ClimateGeometry(geometry_path, climate_path)
    scenarios = NUISANCE + (POSITIVE if stage == "evaluation" else ())
    records = []

    for r in range(start, stop):
        idx = geo.schedule(stage, r)
        x = geo.climate[block][idx]
        features = likelihood_features(x)
        y = np.stack([labels_for(geo, block, stage, s, r)[idx] for s in scenarios])
        bll, cll = axis_log_likelihood_pair(y, features)
        bscore, bf = score_from_ll(bll)
        cscore, cf = score_from_ll(cll)
        for si, s in enumerate(scenarios):
            for model, score, mean_f in (
                (MODELS[0], bscore[si], bf[si]),
                (MODELS[1], cscore[si], cf[si]),
            ):
                row = {
                    "stage": stage,
                    "block": block,
                    "replicate": r,
                    "model": model,
                    "world": s["world"],
                    "amplitude": float(s["amplitude"]),
                    "shared_fraction": float(s.get("shared_fraction", 0.0)),
                    "threshold_sd": float(s.get("threshold_sd", 0.0)),
                    "common_threshold": float(s.get("common_threshold", np.nan)),
                    "axis_concentration": float(s.get("axis_concentration", np.nan)),
                    "score": float(np.mean(score)),
                    "train_fraction": float(np.mean(mean_f)),
                }
                row.update({f"score_fold{k}": float(score[k]) for k in range(4)})
                records.append(row)
        print(json.dumps({"block": block, "stage": stage, "replicate": r + 1, "stop": stop}), flush=True)

    output.mkdir(parents=True)
    table = pd.DataFrame(records)
    table.to_csv(output / "scores.csv", index=False, lineterminator="\n")
    meta = {
        "protocol": "hypervolume-real-climate-synthetic-qualification-v1",
        "specification_commit": SPEC,
        "block": block,
        "stage": stage,
        "replicate_start": start,
        "replicate_stop": stop,
        "score_rows": len(table),
        "models": list(MODELS),
        "axis_grid_sha256": hashlib.sha256(AXES.tobytes()).hexdigest(),
        "geometry_sha256": digest(geometry_path),
        "real_climate_export_sha256": digest(climate_path),
        "runner_sha256": digest(Path(__file__)),
        "scores_sha256": digest(output / "scores.csv"),
        "joint_complete_retained_rows": int(geo.mask.sum()),
        "capacity": geo.capacity,
        "biological_colour_values_read": False,
        "synthetic_colour_labels_only": True,
        "real_climate_values_used_as_geometry": True,
        "parent_empirical_decisions_modified": False,
    }
    (output / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--climate", type=Path, required=True)
    ap.add_argument("--block", choices=list(BLOCKS), required=True)
    ap.add_argument("--stage", choices=["calibration", "evaluation"], required=True)
    ap.add_argument("--rep-start", type=int, required=True)
    ap.add_argument("--rep-stop", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    run_shard(args.geometry, args.climate, args.block, args.output, args.stage, args.rep_start, args.rep_stop)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
