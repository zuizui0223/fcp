#!/usr/bin/env python3
"""Exact species-sharded RGFCA within-species spatial randomization test.

Postoutcome exploratory extension frozen at contract commit
92db63bae67a91ddb2526eb5da5922ddf2c1e7df. Six-species and 34-species
analyses are not read or used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from fcp_pipeline.global_g3 import _upper_triangle_distances_km, _upper_triangle_jsd
from fcp_pipeline.global_rgfca_engine import COLOUR_COLUMNS

SPECIFICATION_COMMIT = "92db63bae67a91ddb2526eb5da5922ddf2c1e7df"
PROTOCOL = "global-rgfca-within-species-spatial-omnibus-v1"
MASTER_SEED = 202609070901
N_PERM = 999
EXPECTED_ROWS = 21424
EXPECTED_SPECIES = 369
EXPECTED_MEASURED_SHA256 = "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def seed_for(species: str, permutation_index: int) -> int:
    text = f"{MASTER_SEED}|{species}|{permutation_index}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def load_rgfca_pool(measured_path: Path) -> pd.DataFrame:
    if digest(measured_path) != EXPECTED_MEASURED_SHA256:
        raise RuntimeError("measured-table checksum drift")
    frame = pd.read_csv(measured_path)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"missing RGFCA columns: {missing}")
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    counts = pool.groupby("species", observed=True).size()
    eligible = sorted(counts[counts >= 40].index.astype(str))
    pool = pool.loc[pool["species"].astype(str).isin(eligible)].copy()
    pool["species"] = pool["species"].astype(str)
    pool = pool.sort_values(["species", "photo_id"], kind="stable").reset_index(drop=True)
    if len(pool) != EXPECTED_ROWS or pool.species.nunique() != EXPECTED_SPECIES:
        raise RuntimeError("frozen 21,424-row / 369-species frame drift")
    if int(pool.groupby("species", observed=True).size().min()) < 40:
        raise RuntimeError("minimum 40-photo species gate drift")
    if pool.photo_id.duplicated().any():
        raise RuntimeError("photo identity duplication")
    vals = pool[["latitude", "longitude", *COLOUR_COLUMNS]].to_numpy(float)
    if not np.isfinite(vals).all():
        raise RuntimeError("nonfinite coordinate/colour values in frozen RGFCA frame")
    return pool


def _pearson_from_ranks(x_rank: np.ndarray, y_rank: np.ndarray) -> float:
    x = np.asarray(x_rank, dtype=float)
    y = np.asarray(y_rank, dtype=float)
    xc = x - x.mean()
    yc = y - y.mean()
    den = float(np.linalg.norm(xc) * np.linalg.norm(yc))
    if den <= 1e-15:
        return 0.0 if np.ptp(y) <= 1e-15 else float("nan")
    return float(np.dot(xc, yc) / den)


def species_observed_and_null(
    latitude: np.ndarray,
    longitude: np.ndarray,
    colours: np.ndarray,
    *,
    species: str,
    batch_size: int = 64,
) -> tuple[float, np.ndarray, dict[str, float]]:
    n = len(latitude)
    if n < 3 or colours.shape != (n, 4):
        raise ValueError("invalid one-species RGFCA input")
    geo = _upper_triangle_distances_km(latitude, longitude)
    jsd = _upper_triangle_jsd(colours)
    if not np.isfinite(geo).all() or not np.isfinite(jsd).all():
        raise RuntimeError("nonfinite pairwise quantity")
    if np.ptp(geo) <= 1e-12:
        raise RuntimeError(f"species has no geographic-distance information: {species}")
    g_rank = rankdata(geo, method="average")
    c_rank = rankdata(jsd, method="average")
    observed = _pearson_from_ranks(g_rank, c_rank)
    if not np.isfinite(observed):
        raise RuntimeError(f"observed rho is nonfinite: {species}")

    upper = np.triu_indices(n, k=1)
    rank_matrix = np.zeros((n, n), dtype=float)
    rank_matrix[upper] = c_rank
    rank_matrix[(upper[1], upper[0])] = c_rank
    gc = g_rank - g_rank.mean()
    gnorm = float(np.linalg.norm(gc))
    cc = c_rank - c_rank.mean()
    cnorm = float(np.linalg.norm(cc))
    null = np.empty(N_PERM, dtype=float)
    if cnorm <= 1e-15:
        null.fill(0.0)
    else:
        u, v = upper
        for start in range(0, N_PERM, batch_size):
            stop = min(start + batch_size, N_PERM)
            perms = np.stack([
                np.random.default_rng(seed_for(species, p)).permutation(n)
                for p in range(start, stop)
            ])
            values = rank_matrix[perms[:, u], perms[:, v]]
            # A vertex permutation preserves the complete pairwise-rank multiset,
            # so its mean and norm are exactly those of c_rank. Because gc sums
            # to zero, centering values is algebraically unnecessary here.
            null[start:stop] = (values @ gc) / (cnorm * gnorm)
    if not np.isfinite(null).all():
        raise RuntimeError(f"nonfinite null rho: {species}")
    audit = {
        "photos": int(n),
        "pairs": int(len(geo)),
        "geographic_distance_min_km": float(np.min(geo)),
        "geographic_distance_max_km": float(np.max(geo)),
        "colour_jsd_min": float(np.min(jsd)),
        "colour_jsd_max": float(np.max(jsd)),
    }
    return observed, null, audit


def run_shard(measured_path: Path, output: Path, shard_index: int, n_shards: int) -> None:
    if output.exists():
        raise FileExistsError("use a fresh shard output")
    if not (0 <= shard_index < n_shards) or n_shards <= 0:
        raise ValueError("invalid species shard")
    pool = load_rgfca_pool(measured_path)
    species_all = sorted(pool.species.unique())
    chosen = species_all[shard_index::n_shards]
    if not chosen:
        raise RuntimeError("empty species shard")
    rows: list[dict[str, object]] = []
    audits: list[dict[str, object]] = []
    for pos, species in enumerate(chosen):
        q = pool.loc[pool.species == species]
        colours = q[list(COLOUR_COLUMNS)].to_numpy(float)
        observed, null, audit = species_observed_and_null(
            q.latitude.to_numpy(float),
            q.longitude.to_numpy(float),
            colours,
            species=species,
        )
        rows.append({"species": species, "permutation_index": -1, "rho": observed})
        rows.extend(
            {"species": species, "permutation_index": int(p), "rho": float(null[p])}
            for p in range(N_PERM)
        )
        audits.append({"species": species, **audit})
        print(json.dumps({
            "shard": shard_index,
            "species_completed": pos + 1,
            "species_in_shard": len(chosen),
            "species": species,
        }), flush=True)

    output.mkdir(parents=True)
    result = pd.DataFrame(rows)
    audit_frame = pd.DataFrame(audits)
    result.to_csv(output / "species_rho_permutations.csv", index=False, lineterminator="\n")
    audit_frame.to_csv(output / "species_audit.csv", index=False, lineterminator="\n")
    meta = {
        "protocol": PROTOCOL,
        "specification_commit": SPECIFICATION_COMMIT,
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "species_count": len(chosen),
        "species": chosen,
        "permutations": N_PERM,
        "master_seed": MASTER_SEED,
        "measured_table_sha256": digest(measured_path),
        "runner_sha256": digest(Path(__file__)),
        "result_sha256": digest(output / "species_rho_permutations.csv"),
        "six_species_used": False,
        "thirty_four_species_used": False,
        "environmental_variables_used": False,
        "parent_G1_reclassified": False,
    }
    (output / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--measured", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--n-shards", type=int, default=20)
    args = ap.parse_args()
    run_shard(args.measured, args.output, args.shard_index, args.n_shards)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
