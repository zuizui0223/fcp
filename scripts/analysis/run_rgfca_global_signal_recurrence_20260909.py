#!/usr/bin/env python3
"""Hypothesis-free RGFCA flower–background signal recurrence analysis.

Implements the protocol frozen in
`docs/RGFCA_GLOBAL_SIGNAL_RECURRENCE_EXPLORATION_PROTOCOL_20260909.md`.
The script uses only NumPy/Pandas plus the Python standard library. EPSG:6933
coordinates are computed directly with the WGS84 ellipsoidal Lambert
cylindrical equal-area formula (standard parallel 30 deg), avoiding an
additional projection-library dependency.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

PROTOCOL_COMMIT = "958eb931b07d768bef41a6fe08a770bf8cef1793"
PALETTE = [
    "white", "yellow", "orange", "red", "pink", "magenta",
    "purple", "blue", "bronze", "green", "brown", "black",
]
CELL_M = 300_000.0
SEPARATION_M = 500_000.0
MIN_PHOTOS = 3
MIN_OBSERVERS = 2
DEFAULT_SEED = 20260909
DEFAULT_REPS = 200

# WGS84 / NSIDC EASE-Grid 2.0 Global (EPSG:6933), ellipsoidal CEA.
_A = 6378137.0
_F = 1.0 / 298.257223563
_E = math.sqrt(_F * (2.0 - _F))
_PHI1 = math.radians(30.0)
_K0 = math.cos(_PHI1) / math.sqrt(1.0 - _E * _E * math.sin(_PHI1) ** 2)


def epsg6933_xy(longitude_deg: np.ndarray, latitude_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lam = np.radians(np.asarray(longitude_deg, dtype=float))
    phi = np.radians(np.asarray(latitude_deg, dtype=float))
    s = np.sin(phi)
    q = (1.0 - _E * _E) * (
        s / (1.0 - _E * _E * s * s)
        - (1.0 / (2.0 * _E)) * np.log((1.0 - _E * s) / (1.0 + _E * s))
    )
    x = _A * _K0 * lam
    y = _A * q / (2.0 * _K0)
    return x, y


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        type=Path,
        default=Path("data/derived/rgfca_reserve_replication_measured_photos_v1.csv"),
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/rgfca_global_signal_recurrence_20260909"),
    )
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--reps", type=int, default=DEFAULT_REPS)
    return p.parse_args()


def load_measurements(path: Path):
    flower_cols = [f"palette_count_{c}" for c in PALETTE]
    bg_cols = [f"background_palette_count_{c}" for c in PALETTE]
    diff_cols = [f"diff_{c}" for c in PALETTE]
    use = [
        "measurement_id", "species", "latitude", "longitude",
        "positional_accuracy_m", "observer_id", "roi_status", "global_classifiable",
    ] + flower_cols + bg_cols
    df = pd.read_csv(path, usecols=use)
    df["flower_total12"] = df[flower_cols].sum(axis=1)
    df["background_total12"] = df[bg_cols].sum(axis=1)
    keep = (
        df["roi_status"].eq("automated_colour_state_admitted")
        & df["global_classifiable"].fillna(False).astype(bool)
        & df["flower_total12"].gt(0)
        & df["background_total12"].gt(0)
        & df["latitude"].notna()
        & df["longitude"].notna()
        & df["positional_accuracy_m"].fillna(np.inf).le(5000)
    )
    df = df.loc[keep].copy()
    for colour in PALETTE:
        df[f"diff_{colour}"] = (
            df[f"palette_count_{colour}"] / df["flower_total12"]
            - df[f"background_palette_count_{colour}"] / df["background_total12"]
        )
    species_names = np.array(sorted(df["species"].unique()))
    species_map = {s: i for i, s in enumerate(species_names)}
    sp = df["species"].map(species_map).to_numpy(np.int32)
    obs = pd.factorize(df["observer_id"].astype(str), sort=True)[0].astype(np.int32)
    values = df[diff_cols].to_numpy(float)
    x, y = epsg6933_xy(df["longitude"].to_numpy(), df["latitude"].to_numpy())
    return df, sp, obs, values, np.asarray(x), np.asarray(y)


def make_engine(sp, obs, values, x, y):
    cx_offset = 64
    cy_offset = 64
    base = 128
    observer_base = int(obs.max()) + 2

    def groups_for_shift(shift_x: float, shift_y: float):
        cx = np.floor((x + shift_x) / CELL_M).astype(np.int16)
        cy = np.floor((y + shift_y) / CELL_M).astype(np.int16)
        key = ((sp.astype(np.int64) * base + (cx + cx_offset).astype(np.int64)) * base
               + (cy + cy_offset).astype(np.int64))
        unique_key, first, inv, count = np.unique(
            key, return_index=True, return_inverse=True, return_counts=True
        )
        state_observer_key = key * observer_base + obs.astype(np.int64)
        unique_state_observer = np.unique(state_observer_key)
        state_key_from_observer = unique_state_observer // observer_base
        state_index = np.searchsorted(unique_key, state_key_from_observer)
        n_observers = np.bincount(state_index, minlength=len(unique_key)).astype(np.int16)
        eligible = (count >= MIN_PHOTOS) & (n_observers >= MIN_OBSERVERS)
        x_centroid = np.bincount(inv, weights=x, minlength=len(unique_key)) / count
        y_centroid = np.bincount(inv, weights=y, minlength=len(unique_key)) / count
        state_species = sp[first]
        return inv, eligible, x_centroid, y_centroid, state_species

    def state_vectors(inv, eligible, rng=None, bootstrap=False):
        order = np.argsort(inv, kind="stable")
        sorted_inv = inv[order]
        starts = np.r_[0, np.flatnonzero(np.diff(sorted_inv)) + 1]
        ends = np.r_[starts[1:], len(order)]
        state_ids = []
        vectors = []
        for start, end in zip(starts, ends):
            state_id = int(sorted_inv[start])
            if not eligible[state_id]:
                continue
            idx = order[start:end]
            state_obs = obs[idx]
            unique_obs = np.unique(state_obs)
            if bootstrap:
                drawn_obs = rng.choice(unique_obs, size=len(unique_obs), replace=True)
                u = rng.random((len(unique_obs), len(idx)))
                u[drawn_obs[:, None] != state_obs[None, :]] = np.inf
                chosen = idx[np.argmin(u, axis=1)]
                vector = values[chosen].mean(axis=0)
            else:
                observer_means = [values[idx[state_obs == o]].mean(axis=0) for o in unique_obs]
                vector = np.vstack(observer_means).mean(axis=0)
            state_ids.append(state_id)
            vectors.append(vector)
        return np.asarray(state_ids, np.int32), np.vstack(vectors)

    def decompose(state_ids, matrix, x_centroid, y_centroid, state_species):
        species_id = state_species[state_ids]
        order = np.argsort(species_id, kind="stable")
        species_id = species_id[order]
        state_ids = state_ids[order]
        matrix = matrix[order]
        starts = np.r_[0, np.flatnonzero(np.diff(species_id)) + 1]
        ends = np.r_[starts[1:], len(species_id)]
        covariance = np.zeros((len(PALETTE), len(PALETTE)))
        mad = np.zeros(len(PALETTE))
        n_species = 0
        n_states = 0
        for start, end in zip(starts, ends):
            if end - start < 2:
                continue
            ids = state_ids[start:end]
            xx = x_centroid[ids]
            yy = y_centroid[ids]
            dx = xx[:, None] - xx[None, :]
            dy = yy[:, None] - yy[None, :]
            if np.sqrt(dx * dx + dy * dy).max() < SEPARATION_M:
                continue
            a = matrix[start:end]
            z = a - a.mean(axis=0, keepdims=True)
            k = len(a)
            covariance += (z.T @ z) / k
            mad += np.abs(z).mean(axis=0)
            n_species += 1
            n_states += k
        covariance /= n_species
        mad /= n_species
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]
        explained = eigenvalues / eigenvalues.sum()
        return n_species, n_states, eigenvalues, eigenvectors, explained, mad

    return groups_for_shift, state_vectors, decompose


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    df, sp, obs, values, x, y = load_measurements(args.input)
    groups_for_shift, state_vectors, decompose = make_engine(sp, obs, values, x, y)

    inv, eligible, xc, yc, state_sp = groups_for_shift(0.0, 0.0)
    state_ids, matrix = state_vectors(inv, eligible, bootstrap=False)
    ref_n_species, ref_n_states, _, ref_evec, ref_exp, _ = decompose(
        state_ids, matrix, xc, yc, state_sp
    )
    reference = ref_evec[:, :5]
    pd.DataFrame(reference, index=PALETTE, columns=[f"M{i+1}" for i in range(5)]).to_csv(
        args.outdir / "reference_loadings.csv"
    )
    pd.DataFrame({
        "mode": [f"M{i+1}" for i in range(5)],
        "explained": ref_exp[:5],
    }).to_csv(args.outdir / "reference_variance.csv", index=False)

    permutations = list(itertools.permutations(range(5)))
    rng = np.random.default_rng(args.seed)
    rows = []
    for rep in range(1, args.reps + 1):
        shift_x = float(rng.uniform(0, CELL_M))
        shift_y = float(rng.uniform(0, CELL_M))
        inv, eligible, xc, yc, state_sp = groups_for_shift(shift_x, shift_y)
        state_ids, matrix = state_vectors(inv, eligible, rng=rng, bootstrap=True)
        n_species, n_states, _, evec, explained, mad = decompose(
            state_ids, matrix, xc, yc, state_sp
        )
        first5 = evec[:, :5]
        similarity = np.abs(reference.T @ first5)
        perm = max(
            permutations,
            key=lambda p: sum(similarity[i, p[i]] for i in range(5)),
        )
        matched = np.column_stack([first5[:, perm[i]] for i in range(5)])
        cosines = []
        for i in range(5):
            dot = float(reference[:, i] @ matched[:, i])
            if dot < 0:
                matched[:, i] *= -1
                dot = -dot
            cosines.append(dot)
        singular_values = np.linalg.svd(reference.T @ first5, compute_uv=False)
        row = {
            "rep": rep,
            "shift_x_m": shift_x,
            "shift_y_m": shift_y,
            "n_species": n_species,
            "n_states": n_states,
            "subspace_sv_mean": singular_values.mean(),
            "subspace_sv_min": singular_values.min(),
            **{f"subspace_sv{i+1}": singular_values[i] for i in range(5)},
        }
        for i in range(5):
            mode = f"M{i+1}"
            row[f"{mode}_cosine"] = cosines[i]
            row[f"{mode}_matched_run_rank"] = perm[i] + 1
            row[f"{mode}_explained"] = explained[perm[i]]
            top = np.argsort(np.abs(matched[:, i]))[-2:][::-1]
            row[f"{mode}_top1"] = PALETTE[top[0]]
            row[f"{mode}_top2"] = PALETTE[top[1]]
            for colour, value in zip(PALETTE, matched[:, i]):
                row[f"{mode}_loading_{colour}"] = value
        for colour, value in zip(PALETTE, mad):
            row[f"mad_{colour}"] = value
        rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(args.outdir / "recurrence_realizations.csv", index=False)

    mode_summary = []
    top_frequency = []
    for i in range(5):
        mode = f"M{i+1}"
        cosine = result[f"{mode}_cosine"]
        exp = result[f"{mode}_explained"]
        mode_summary.append({
            "mode": mode,
            "reference_explained": ref_exp[i],
            "cosine_median": cosine.median(),
            "cosine_q025": cosine.quantile(0.025),
            "cosine_q975": cosine.quantile(0.975),
            "p_alignment_ge_0_80": (cosine >= 0.80).mean(),
            "p_alignment_ge_0_90": (cosine >= 0.90).mean(),
            "explained_median": exp.median(),
            "explained_q025": exp.quantile(0.025),
            "explained_q975": exp.quantile(0.975),
        })
        for colour in PALETTE:
            top_frequency.append({
                "mode": mode,
                "colour": colour,
                "top2_frequency": (
                    result[f"{mode}_top1"].eq(colour)
                    | result[f"{mode}_top2"].eq(colour)
                ).mean(),
            })
    pd.DataFrame(mode_summary).to_csv(args.outdir / "mode_recurrence_summary.csv", index=False)
    pd.DataFrame(top_frequency).to_csv(args.outdir / "mode_top2_colour_recurrence.csv", index=False)

    global_summary = pd.DataFrame([{
        "n_reps": args.reps,
        "seed": args.seed,
        "n_eligible_photos": len(df),
        "n_eligible_species": df["species"].nunique(),
        "reference_n_species": ref_n_species,
        "reference_n_states": ref_n_states,
        "species_median": result["n_species"].median(),
        "species_min": result["n_species"].min(),
        "species_max": result["n_species"].max(),
        "states_median": result["n_states"].median(),
        "states_min": result["n_states"].min(),
        "states_max": result["n_states"].max(),
        "subspace_sv_mean_median": result["subspace_sv_mean"].median(),
        "subspace_sv_mean_q025": result["subspace_sv_mean"].quantile(0.025),
        "subspace_sv_mean_q975": result["subspace_sv_mean"].quantile(0.975),
        "subspace_sv_min_median": result["subspace_sv_min"].median(),
        "subspace_sv_min_q025": result["subspace_sv_min"].quantile(0.025),
    }])
    global_summary.to_csv(args.outdir / "global_recurrence_summary.csv", index=False)

    mad_summary = []
    for colour in PALETTE:
        s = result[f"mad_{colour}"]
        mad_summary.append({
            "colour": colour,
            "mean": s.mean(),
            "median": s.median(),
            "q025": s.quantile(0.025),
            "q975": s.quantile(0.975),
        })
    pd.DataFrame(mad_summary).sort_values("median", ascending=False).to_csv(
        args.outdir / "palette_mad_recurrence_summary.csv", index=False
    )

    with (args.outdir / "run_manifest.json").open("w", encoding="utf-8") as f:
        json.dump({
            "protocol_commit": PROTOCOL_COMMIT,
            "implementation": "exact observer/photo bootstrap; analytic EPSG:6933",
            "seed": args.seed,
            "n_reps": args.reps,
            "cell_km": CELL_M / 1000.0,
            "min_photos": MIN_PHOTOS,
            "min_observers": MIN_OBSERVERS,
            "min_state_separation_km": SEPARATION_M / 1000.0,
        }, f, indent=2)


if __name__ == "__main__":
    main()
