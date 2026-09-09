#!/usr/bin/env python3
"""Frozen, synthetic species-disjoint hypervolume transfer benchmark.

No biological colour data or real climate rasters are read. Gaussian affinity is
not geometric ellipsoid intersection. All rejection rates are method-specific,
not a reclassification of RGFCA. Specification: bb23f879f6d84a27d0c7fdbde430b9d4f4b4f7ad.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import chi2

SPECIFICATION_COMMIT = "bb23f879f6d84a27d0c7fdbde430b9d4f4b4f7ad"
SEED = 2026090701
N_SPECIES, N_TRAIN, N_PHOTOS = 60, 30, 20
AMPLITUDES = (0.5, 1.0, 2.0)
NUISANCE = ("geographic_specific", "environmental_specific", "mixed_specific", "no_structure")
WORLDS = ("geographic_shared", "environmental_shared") + NUISANCE
N_CAL, N_EVAL = 1000, 500
SOURCE_SHA256 = "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4"


def seed_for(*parts: object) -> int:
    text = "|".join(map(str, (SEED,) + parts))
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def unit_vectors(rng: np.random.Generator, shape: tuple[int, ...]) -> np.ndarray:
    a = rng.normal(size=shape)
    return a / np.linalg.norm(a, axis=-1, keepdims=True)


def environment(xyz: np.ndarray, mapping: str = "harmonic") -> np.ndarray:
    if mapping == "confounded":
        return np.sqrt(3.0) * xyz[..., :2]
    if mapping != "harmonic":
        raise ValueError(f"unknown mapping: {mapping}")
    x, y = xyz[..., 0], xyz[..., 1]
    return np.stack([x**4 - 6*x*x*y*y + y**4, 4*x*y*(x*x-y*y)], axis=-1) * np.sqrt(315.0)/8


def generate_world(stage: str, mapping: str, world: str, amplitude: float, replicate: int):
    if stage not in {"calibration", "evaluation", "unit_test"}:
        raise ValueError(stage)
    if world not in WORLDS + ("confounded_shared",):
        raise ValueError(world)
    rng_g = np.random.default_rng(seed_for(stage, "geometry", replicate))
    xyz = unit_vectors(rng_g, (N_SPECIES, N_PHOTOS, 3))
    g, e = np.sqrt(3.0)*xyz, environment(xyz, mapping)
    rng = np.random.default_rng(seed_for(stage, mapping, world, amplitude, replicate))
    ng = unit_vectors(rng, (N_SPECIES, 3))
    ne = unit_vectors(rng, (N_SPECIES, 2))
    if world == "geographic_shared":
        ng[:] = ng[0]
    if world == "environmental_shared":
        ne[:] = ne[0]
    pg, pe = np.einsum("snd,sd->sn", g, ng), np.einsum("snd,sd->sn", e, ne)
    if world.startswith("geographic"):
        projection = pg
    elif world.startswith("environmental"):
        projection = pe
    elif world == "mixed_specific":
        projection = np.where(rng.random((N_SPECIES, 1)) < 0.5, pg, pe)
    elif world == "confounded_shared":
        if mapping != "confounded":
            raise ValueError("confounded world requires confounded mapping")
        projection = g[..., 0]
    else:
        projection = np.zeros_like(pg)
    signs = rng.choice([-1., 1.], size=(N_SPECIES, 1))
    latent = amplitude * signs * np.tanh(projection/0.5) + rng.normal(size=projection.shape)
    return g, e, latent > 0


def shrink_covariance(cov: np.ndarray) -> np.ndarray:
    d = cov.shape[-1]
    target = np.trace(cov, axis1=-2, axis2=-1)[..., None, None] / d
    return 0.8*cov + (0.2*target + 1e-8)*np.eye(d)


def fit_classes(x: np.ndarray, labels: np.ndarray):
    """Input shapes (..., observations, dimension), (..., observations)."""
    if x.shape[:-1] != labels.shape or not np.isfinite(x).all():
        raise ValueError("non-finite coordinates or incompatible labels")
    means, covs, counts = [], [], []
    for k in (False, True):
        w = (labels == k).astype(float)
        n = w.sum(axis=-1)
        mu = (x*w[..., None]).sum(axis=-2) / np.maximum(n[..., None], 1)
        delta = x-mu[..., None, :]
        cov = np.einsum("...ni,...nj,...n->...ij", delta, delta, w)
        cov /= np.maximum(n-1, 1)[..., None, None]
        means.append(mu)
        covs.append(shrink_covariance(cov))
        counts.append(n)
    return np.stack(means, axis=-2), np.stack(covs, axis=-3), np.minimum(*counts) >= 4


def affinity(mu: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Analytic Gaussian Bhattacharyya coefficient; NOT intersection volume."""
    avg = (cov[..., 0, :, :] + cov[..., 1, :, :])/2
    dm = mu[..., 0, :] - mu[..., 1, :]
    quadratic = np.einsum("...i,...ij,...j->...", dm, np.linalg.inv(avg), dm)/8
    ld_avg = np.linalg.slogdet(avg)[1]
    ld = np.linalg.slogdet(cov)[1]
    log_bc = -quadratic - 0.5*ld_avg + 0.25*(ld[..., 0]+ld[..., 1])
    return np.exp(np.minimum(log_bc, 0.0))


def log_total_volume(x: np.ndarray) -> np.ndarray:
    n, d = x.shape[-2:]
    delta = x-x.mean(axis=-2, keepdims=True)
    cov = np.einsum("...ni,...nj->...ij", delta, delta)/(n-1)
    cov = shrink_covariance(cov)
    log_unit_ball = d/2*math.log(math.pi)-math.lgamma(d/2+1)
    return log_unit_ball+d/2*math.log(chi2.ppf(0.95, d))+0.5*np.linalg.slogdet(cov)[1]


def adjusted_rand_binary(truth: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    """Broadcast binary ARI; constant truth/prediction contributes zero information."""
    truth, predicted = np.broadcast_arrays(truth, predicted)
    n = truth.shape[-1]
    if n < 2:
        raise ValueError("ARI requires at least two observations")
    a = np.sum(truth & predicted, axis=-1).astype(float)
    t = truth.sum(axis=-1).astype(float)
    p = predicted.sum(axis=-1).astype(float)
    c2 = lambda z: z*(z-1)/2
    cells = c2(a)+c2(t-a)+c2(p-a)+c2(n-t-p+a)
    row = c2(t)+c2(n-t)
    col = c2(p)+c2(n-p)
    expectation = row*col/c2(n)
    denominator = (row+col)/2-expectation
    out = np.zeros_like(denominator)
    keep = (denominator > 0) & (t > 0) & (t < n) & (p > 0) & (p < n)
    np.divide(cells-expectation, denominator, out=out, where=keep)
    return out


def score_domain(x: np.ndarray, labels: np.ndarray) -> dict[str, np.ndarray]:
    """Batch of worlds: (B, species, photos, dimensions)."""
    mu, cov, valid = fit_classes(x, labels)
    target = x[:, N_TRAIN:]
    discrim = []
    for k in (0, 1):
        m = mu[:, :N_TRAIN, k]
        c = cov[:, :N_TRAIN, k]
        diff = target[:, None] - m[:, :, None, None]
        q = np.einsum("banid,bade,banie->bani", diff, np.linalg.inv(c), diff, optimize=True)
        discrim.append(-0.5*(q+np.linalg.slogdet(c)[1][:, :, None, None]))
    prediction = discrim[1] > discrim[0]
    ari = adjusted_rand_binary(labels[:, None, N_TRAIN:], prediction)
    valid_pairs = valid[:, :N_TRAIN, None] & valid[:, None, N_TRAIN:]
    transfer = np.where(valid_pairs, ari, 0).mean(axis=(1, 2))
    separation = np.where(valid, 1-affinity(mu, cov), 0).mean(axis=1)
    return {"transfer": transfer, "separation": separation,
            "valid_species_fraction": valid.mean(axis=1),
            "log_total_volume_mean": log_total_volume(x).mean(axis=1)}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("no rows")
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def run_arm(output: Path, mapping: str, stage: str, world: str, amplitude: float, batch: int = 10) -> Path:
    reps = N_CAL if stage == "calibration" else N_EVAL
    destination = output / f"{mapping}__{stage}__{world}__a{amplitude:g}.csv"
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite {destination}")
    rows = []
    for start in range(0, reps, batch):
        end = min(start+batch, reps)
        generated = [generate_world(stage, mapping, world, amplitude, r) for r in range(start, end)]
        g, e, y = (np.stack([z[j] for z in generated]) for j in range(3))
        for name, x in (("geographic", g), ("environmental", e)):
            scores = score_domain(x, y)
            for i, r in enumerate(range(start, end)):
                row = {"mapping": mapping, "stage": stage, "world": world,
                       "amplitude": amplitude, "replicate": r, "representation": name}
                row.update({k: float(v[i]) for k, v in scores.items()})
                rows.append(row)
    write_csv(destination, rows)
    print(f"completed {destination.name}: {reps} worlds", flush=True)
    return destination


def wilson(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    q = k/n; a = 1+z*z/n
    center = (q+z*z/(2*n))/a
    radius = z*math.sqrt(q*(1-q)/n+z*z/(4*n*n))/a
    return max(0., center-radius), min(1., center+radius)


def summarize(output: Path) -> dict:
    import pandas as pd
    files = sorted(output.glob("*__*__*__a*.csv"))
    df = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
    thresholds, decisions = [], []
    for mapping in ("harmonic", "confounded"):
        for amp in AMPLITUDES:
            threshold = {}
            for domain in ("geographic", "environmental"):
                for world in NUISANCE:
                    x = df[(df.mapping == mapping) & (df.stage == "calibration") & (df.world == world) & (df.amplitude == amp) & (df.representation == domain)]
                    if len(x) != N_CAL or x.replicate.nunique() != N_CAL:
                        raise RuntimeError(f"incomplete calibration: {mapping}, {amp}, {world}, {domain}")
                    value = float(np.quantile(x.transfer, .975, method="higher"))
                    thresholds.append(dict(mapping=mapping, amplitude=amp, representation=domain, nuisance=world, threshold=value))
                threshold[domain] = max(r["threshold"] for r in thresholds if r["mapping"] == mapping and r["amplitude"] == amp and r["representation"] == domain)
            expected = WORLDS if mapping == "harmonic" else ("confounded_shared",) + NUISANCE
            for world in expected:
                x = df[(df.mapping == mapping) & (df.stage == "evaluation") & (df.world == world) & (df.amplitude == amp)]
                pivot = x.pivot(index="replicate", columns="representation", values="transfer")
                if len(pivot) != N_EVAL or not np.isfinite(pivot.to_numpy()).all():
                    raise RuntimeError(f"incomplete evaluation: {mapping}, {amp}, {world}")
                g, e = pivot.geographic > threshold["geographic"], pivot.environmental > threshold["environmental"]
                cats = {"geography_only": g & ~e, "environment_only": e & ~g, "both": g & e, "neither": ~g & ~e,
                        "any": g | e, "geographic_positive": g, "environmental_positive": e}
                for label, mask in cats.items():
                    k = int(mask.sum()); lo, hi = wilson(k, len(pivot))
                    decisions.append(dict(mapping=mapping, amplitude=amp, world=world, decision=label,
                                          count=k, replicates=len(pivot), rate=k/len(pivot), wilson95_low=lo, wilson95_high=hi))
    write_csv(output/"calibration_thresholds.csv", thresholds)
    write_csv(output/"decision_rates.csv", decisions)
    df.groupby(["mapping", "stage", "world", "amplitude", "representation"], observed=True)[["transfer", "separation", "valid_species_fraction", "log_total_volume_mean"]].mean().to_csv(output/"mean_statistics.csv")
    # Each world uses identical coordinates at a given replicate, so total volume must agree.
    check = df.groupby(["stage", "replicate", "representation", "mapping"]).log_total_volume_mean.agg(["min", "max"])
    max_error = float((check["max"]-check["min"]).abs().max())
    if max_error > 1e-12:
        raise RuntimeError("total volume acquired colour information")
    result = {"protocol": "hypervolume-geo-environment-transfer-benchmark-v1", "status": "complete_synthetic_benchmark",
              "specification_commit": SPECIFICATION_COMMIT, "seed": SEED, "species_per_world": N_SPECIES,
              "training_species": N_TRAIN, "evaluation_species": N_SPECIES-N_TRAIN, "photos_per_species": N_PHOTOS,
              "simulated_worlds": int(len(df)//2), "maximum_colour_invariant_total_volume_error": max_error,
              "numpy_version": np.__version__, "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "biological_colour_values_used": False, "parent_decisions_modified": False,
              "claim_ceiling": "Species-disjoint predictive transfer of low-order Gaussian summaries under specified synthetic worlds. No causal climate/geography identification; no actual-geometry power calibration; no observed niche result.",
              "uncertainty": "Wilson intervals quantify evaluation Monte Carlo error conditional on estimated fixed thresholds; threshold-calibration uncertainty is omitted.",
              "decisions": decisions}
    (output/"summary.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k != "decisions"}, indent=2))
    return result


def geometry_audit(path: Path, output: Path) -> dict:
    import pandas as pd
    audit = json.loads(path.with_name("audit.json").read_text())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if audit["source_sha256"] != SOURCE_SHA256 or digest != audit["export_sha256"]:
        raise RuntimeError("geometry lineage checksum failed")
    cols = ["photo_id", "species", "latitude", "longitude", "global_classifiable"]
    df = pd.read_csv(path, usecols=cols)
    if len(df) != 50000 or df.species.nunique() != 500 or df.photo_id.duplicated().any():
        raise RuntimeError("geometry counts/uniqueness failed")
    if not np.isfinite(df[["latitude", "longitude"]].to_numpy()).all() or not df.latitude.between(-90, 90).all() or not df.longitude.between(-180, 180).all():
        raise RuntimeError("invalid coordinates")
    ok = df.global_classifiable.astype(str).str.lower().isin(["true", "1"])
    if ok.sum() != 25377:
        raise RuntimeError("classifiable count drift")
    counts = df.loc[ok].groupby("species").size()
    eligible = set(counts[counts >= 40].index)
    if len(eligible) != 369:
        raise RuntimeError("eligibility drift")
    rows = []
    for sp, group in df[df.species.isin(eligible)].groupby("species", sort=True):
        for frame, sub in (("complete_same_369", group), ("classifiable_same_369", group.loc[ok.loc[group.index]])):
            lat, lon = np.deg2rad(sub.latitude.to_numpy()), np.deg2rad(sub.longitude.to_numpy())
            xyz = np.stack([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)], axis=-1)
            for domain, x in (("geographic", np.sqrt(3)*xyz), ("synthetic_environmental", environment(xyz))):
                ev = np.maximum(np.linalg.eigvalsh(np.cov(x, rowvar=False)), 0)
                p = ev/ev.sum() if ev.sum() > 0 else np.zeros_like(ev)
                keep = p > 0
                effective_rank = float(np.exp(-np.sum(p[keep]*np.log(p[keep])))) if keep.any() else 0.
                rows.append(dict(species=sp, frame=frame, representation=domain, photos=len(x), effective_rank=effective_rank,
                                 minimum_eigenvalue=float(ev.min()), log_total_volume=float(log_total_volume(x))))
    write_csv(output/"geometry_species.csv", rows)
    r = pd.DataFrame(rows)
    result = {"status": "complete_actual_coordinate_and_mask_audit_only", "source_sha256": SOURCE_SHA256,
              "export_sha256": digest, "rows": len(df), "species": int(df.species.nunique()), "classifiable_rows": int(ok.sum()),
              "eligible_species": len(eligible), "eligible_classifiable_rows": int((ok & df.species.isin(eligible)).sum()),
              "ineligible_species": 500-len(eligible), "biological_colour_values_used": False,
              "real_climate_used": False, "sample_size_confounded_raw_volume_comparison": True,
              "medians": r.groupby(["frame", "representation"])[["photos", "effective_rank", "log_total_volume"]].median().reset_index().to_dict("records")}
    (output/"geometry_audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--mapping", choices=["harmonic", "confounded"], default="harmonic")
    p.add_argument("--stage", choices=["calibration", "evaluation"])
    p.add_argument("--world", choices=list(WORLDS)+["confounded_shared"])
    p.add_argument("--amplitude", type=float, choices=AMPLITUDES)
    p.add_argument("--batch-size", type=int, default=10)
    p.add_argument("--summarize", action="store_true")
    p.add_argument("--geometry", type=Path)
    a = p.parse_args(); a.output_dir.mkdir(parents=True, exist_ok=True)
    if a.geometry:
        geometry_audit(a.geometry, a.output_dir)
    elif a.summarize:
        summarize(a.output_dir)
    elif a.stage and a.world and a.amplitude is not None:
        if a.stage == "calibration" and a.world not in NUISANCE:
            p.error("calibration must use a nuisance world")
        if a.batch_size < 1:
            p.error("batch size must be positive")
        run_arm(a.output_dir, a.mapping, a.stage, a.world, a.amplitude, a.batch_size)
    else:
        p.error("choose --geometry, --summarize, or a complete arm")


if __name__ == "__main__":
    main()
