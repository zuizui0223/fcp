#!/usr/bin/env python3
"""Run frozen Step 2: are chromatic mixed_uncertain rows morph intermediates?"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_mixed_step2_20260909"
OUT.mkdir(parents=True, exist_ok=True)
DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
BIO = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
FRACTIONS = [f"flower_fraction_{x}" for x in BIO]
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
AMBIG_STATUS = "not_evaluable_ambiguous_palette_composition"
SEED = 20260909
N_COORD_PERM = 2000
N_SPECIES_PERM = 20000


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def normalize_rows(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, float)
    finite = np.isfinite(x).all(axis=1)
    sums = np.nansum(x, axis=1)
    good = finite & (sums > 0)
    out = np.full_like(x, np.nan, dtype=float)
    out[good] = x[good] / sums[good, None]
    return out, good


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return float(stats.spearmanr(a, b).statistic)


def upper_perm_p(a: np.ndarray, b: np.ndarray, *, seed: int, n_perm: int) -> tuple[float, np.ndarray]:
    obs = spearman(a, b)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, float)
    for i in range(n_perm):
        null[i] = spearman(a[rng.permutation(len(a))], b)
    p = float((1 + np.sum(null >= obs)) / (n_perm + 1))
    return p, null


def load_table(path: Path) -> pd.DataFrame:
    use = ["species", "morph", "measurement_status", "global_classifiable"] + FRACTIONS
    df = pd.read_csv(path, usecols=use)
    df["species"] = df["species"].astype(str)
    df["morph"] = df["morph"].fillna("").astype(str)
    df["measurement_status"] = df["measurement_status"].fillna("").astype(str)
    df["global_classifiable_bool"] = bool_series(df["global_classifiable"])
    return df


def reserve_reference(reserve: pd.DataFrame) -> tuple[np.ndarray, float, pd.DataFrame]:
    keep = reserve["global_classifiable_bool"] & reserve["morph"].isin(MORPHS)
    work = reserve.loc[keep].copy()
    x, good = normalize_rows(work[FRACTIONS].to_numpy(float))
    work = work.loc[good].copy()
    x = x[good]
    work["_row"] = np.arange(len(work))

    centroids = []
    for morph in MORPHS:
        idx = np.flatnonzero(work["morph"].to_numpy(str) == morph)
        if len(idx) < 10:
            raise RuntimeError(f"reserve morph {morph} has insufficient centroid support: {len(idx)}")
        centroids.append(x[idx].mean(axis=0))
    c = np.vstack(centroids)
    if not np.allclose(c.sum(axis=1), 1.0, atol=1e-10):
        raise RuntimeError("reserve morph centroids do not sum to one")

    own_dist = np.empty(len(work), float)
    morph_arr = work["morph"].to_numpy(str)
    for j, morph in enumerate(MORPHS):
        idx = np.flatnonzero(morph_arr == morph)
        own_dist[idx] = np.linalg.norm(x[idx] - c[j], axis=1)
    radius95 = float(np.quantile(own_dist, 0.95))
    ref = pd.DataFrame(c, index=MORPHS, columns=BIO)
    ref["n_reference_rows"] = [int(np.sum(morph_arr == m)) for m in MORPHS]
    return c, radius95, ref


def bridge_geometry(x: np.ndarray, centroids: np.ndarray, radius95: float) -> pd.DataFrame:
    x = np.asarray(x, float)
    vertex_dist = np.linalg.norm(x[:, None, :] - centroids[None, :, :], axis=2)
    nearest_vertex = vertex_dist.argmin(axis=1)
    d_vertex = vertex_dist[np.arange(len(x)), nearest_vertex]

    pairs = list(combinations(range(len(MORPHS)), 2))
    all_dist = np.empty((len(x), len(pairs)), float)
    all_t = np.empty((len(x), len(pairs)), float)
    for k, (i, j) in enumerate(pairs):
        a = centroids[i]
        v = centroids[j] - a
        denom = float(v @ v)
        if denom <= 0:
            raise RuntimeError(f"degenerate centroid pair {MORPHS[i]}, {MORPHS[j]}")
        t = ((x - a) @ v) / denom
        t_clip = np.clip(t, 0.0, 1.0)
        proj = a[None, :] + t_clip[:, None] * v[None, :]
        all_t[:, k] = t_clip
        all_dist[:, k] = np.linalg.norm(x - proj, axis=1)

    best = all_dist.argmin(axis=1)
    d_segment = all_dist[np.arange(len(x)), best]
    t = all_t[np.arange(len(x)), best]
    pair_left = np.array([MORPHS[pairs[k][0]] for k in best], dtype=object)
    pair_right = np.array([MORPHS[pairs[k][1]] for k in best], dtype=object)
    interior = (t >= 0.20) & (t <= 0.80) & (d_segment <= radius95)
    return pd.DataFrame({
        "nearest_vertex": [MORPHS[i] for i in nearest_vertex],
        "d_vertex": d_vertex,
        "best_pair_left": pair_left,
        "best_pair_right": pair_right,
        "t": t,
        "d_segment": d_segment,
        "bridge_advantage": d_vertex - d_segment,
        "interior_bridge": interior,
    })


def coordinate_null(x: np.ndarray, centroids: np.ndarray, radius95: float) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    null = np.empty(N_COORD_PERM, float)
    for b in range(N_COORD_PERM):
        perm = rng.permutation(x.shape[1])
        g = bridge_geometry(x[:, perm], centroids, radius95)
        null[b] = float(g["interior_bridge"].mean())
    return null


def reconstruct_species_D(discovery: pd.DataFrame) -> pd.DataFrame:
    keep = discovery["global_classifiable_bool"] & discovery["morph"].isin(MORPHS)
    work = discovery.loc[keep]
    rows = []
    for sp, g in work.groupby("species", sort=True):
        if len(g) < 40:
            continue
        cnt = Counter(g["morph"].astype(str))
        p = np.array([cnt.get(m, 0) / len(g) for m in MORPHS], float)
        rows.append({
            "species": sp,
            "n_classifiable": int(len(g)),
            "D": float(1.0 - np.sum(p * p)),
            "second_fraction": float(np.sort(p)[-2]),
        })
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    if len(out) != 369:
        raise RuntimeError(f"discovery fingerprint species mismatch: {len(out)}")
    if abs(float(out["D"].max()) - 0.707645) > 1e-5:
        raise RuntimeError(f"D max fingerprint mismatch: {out['D'].max()}")
    return out


def main() -> None:
    discovery = load_table(DISCOVERY)
    reserve = load_table(RESERVE)
    centroids, radius95, ref = reserve_reference(reserve)
    species = reconstruct_species_D(discovery)

    mixed = discovery["morph"].eq("mixed_uncertain")
    status_counts = discovery.loc[mixed, "measurement_status"].value_counts(dropna=False).sort_index()
    ambiguous = mixed & discovery["measurement_status"].eq(AMBIG_STATUS)
    amb = discovery.loc[ambiguous].copy()
    x_all, good = normalize_rows(amb[FRACTIONS].to_numpy(float))
    amb = amb.loc[good].copy()
    x = x_all[good]
    if len(amb) == 0:
        raise RuntimeError("no chromatically evaluable ambiguous discovery rows")

    geom = bridge_geometry(x, centroids, radius95)
    geom.index = amb.index
    amb_geom = pd.concat([amb[["species", "morph", "measurement_status"]], geom], axis=1)
    observed_bridge = float(geom["interior_bridge"].mean())
    coord_null = coordinate_null(x, centroids, radius95)
    p_bridge = float((1 + np.sum(coord_null >= observed_bridge)) / (N_COORD_PERM + 1))

    # Species-level rates use all raw measured rows as denominator.
    bridge_index = set(amb_geom.index[amb_geom["interior_bridge"]].tolist())
    metrics = []
    for _, row in species.iterrows():
        sp = row["species"]
        g = discovery.loc[discovery["species"].eq(sp)]
        n_raw = len(g)
        if n_raw == 0:
            raise RuntimeError(f"no raw rows for species {sp}")
        gidx = set(g.index.tolist())
        n_mixed = int(g["morph"].eq("mixed_uncertain").sum())
        n_amb = int((g["morph"].eq("mixed_uncertain") & g["measurement_status"].eq(AMBIG_STATUS)).sum())
        n_bridge = int(len(gidx.intersection(bridge_index)))
        n_technical = int(n_mixed - n_amb)
        metrics.append({
            **row.to_dict(),
            "n_raw": int(n_raw),
            "n_mixed": n_mixed,
            "n_ambiguous_palette": n_amb,
            "n_bridge": n_bridge,
            "n_technical_mixed": n_technical,
            "mixed_rate": n_mixed / n_raw,
            "ambiguous_palette_rate": n_amb / n_raw,
            "bridge_rate": n_bridge / n_raw,
            "technical_mixed_rate": n_technical / n_raw,
        })
    sm = pd.DataFrame(metrics)

    d = sm["D"].to_numpy(float)
    bridge_rate = sm["bridge_rate"].to_numpy(float)
    rho_bridge = spearman(d, bridge_rate)
    p_species, species_null = upper_perm_p(d, bridge_rate, seed=SEED + 1, n_perm=N_SPECIES_PERM)
    diagnostics = {
        "rho_D_mixed_rate": spearman(d, sm["mixed_rate"].to_numpy(float)),
        "rho_D_ambiguous_palette_rate": spearman(d, sm["ambiguous_palette_rate"].to_numpy(float)),
        "rho_D_bridge_rate": rho_bridge,
        "rho_D_technical_mixed_rate": spearman(d, sm["technical_mixed_rate"].to_numpy(float)),
    }

    bridge_counts = (
        amb_geom.loc[amb_geom["interior_bridge"]]
        .assign(pair=lambda z: z["best_pair_left"] + "__" + z["best_pair_right"])
        ["pair"].value_counts().sort_values(ascending=False)
    )

    geometry_pass = bool(p_bridge < 0.05)
    species_pass = bool(rho_bridge > 0 and p_species < 0.05)
    survives = bool(geometry_pass and species_pass)
    if survives:
        verdict = "CHROMATIC_BRIDGES_SUPPORT_CLAIM1_INTERPRETATION"
    elif geometry_pass:
        verdict = "GLOBAL_BRIDGE_GEOMETRY_ONLY_NO_D_LINK"
    else:
        verdict = "MIXED_REMAINS_STRUCTURAL_MISSINGNESS"

    all_mixed_n = int(mixed.sum())
    ambiguous_all_n = int(ambiguous.sum())
    ambiguous_valid_n = int(len(amb))
    result = {
        "analysis": "polymorphism_mixed_step2",
        "date_jst": "2026-09-09",
        "protocol": "docs/POLYMORPHISM_MIXED_STEP2_PROTOCOL_20260909.md",
        "new_image_acquisition": False,
        "discovery_rows": int(len(discovery)),
        "discovery_species_D": int(len(sm)),
        "D_fingerprint": {
            "D_min": float(sm["D"].min()),
            "D_max": float(sm["D"].max()),
            "second_ge_0_10": float((sm["second_fraction"] >= 0.10).mean()),
            "second_ge_0_20": float((sm["second_fraction"] >= 0.20).mean()),
        },
        "mixed_decomposition": {
            "all_mixed_uncertain": all_mixed_n,
            "all_mixed_fraction_of_50000": all_mixed_n / len(discovery),
            "ambiguous_palette_rows": ambiguous_all_n,
            "ambiguous_palette_fraction_of_all_mixed": ambiguous_all_n / all_mixed_n if all_mixed_n else None,
            "ambiguous_palette_valid_geometry_rows": ambiguous_valid_n,
            "technical_or_other_mixed_rows": all_mixed_n - ambiguous_all_n,
            "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        },
        "reserve_reference": {
            "within_morph_radius_q95": radius95,
            "centroids": {MORPHS[i]: {BIO[j]: float(centroids[i, j]) for j in range(len(BIO))} for i in range(len(MORPHS))},
            "n_rows_per_morph": {m: int(ref.loc[m, "n_reference_rows"]) for m in MORPHS},
        },
        "bridge_geometry": {
            "observed_interior_bridge_fraction": observed_bridge,
            "coordinate_permutation_reps": N_COORD_PERM,
            "null_mean": float(np.mean(coord_null)),
            "null_q025": float(np.quantile(coord_null, 0.025)),
            "null_q975": float(np.quantile(coord_null, 0.975)),
            "p_bridge_upper": p_bridge,
            "bridge_pair_counts": {str(k): int(v) for k, v in bridge_counts.items()},
        },
        "species_level": {
            **diagnostics,
            "bridge_rate_upper_label_permutation_p": p_species,
            "raw_rows_per_species_min": int(sm["n_raw"].min()),
            "raw_rows_per_species_max": int(sm["n_raw"].max()),
        },
        "decision": {
            "geometry_pass": geometry_pass,
            "species_D_link_pass": species_pass,
            "intermediate_colour_interpretation_supported": survives,
            "verdict": verdict,
            "D_numeric_underestimation_claim_allowed": False,
        },
    }

    ref.to_csv(OUT / "reserve_morph_centroids.csv")
    amb_geom.to_csv(OUT / "ambiguous_row_geometry.csv", index_label="source_row_index")
    sm.to_csv(OUT / "species_mixed_metrics.csv", index=False)
    pd.DataFrame({"bridge_fraction": coord_null}).to_csv(OUT / "coordinate_permutation_null.csv", index=False)
    pd.DataFrame({"rho": species_null}).to_csv(OUT / "species_bridge_rate_label_null.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Polymorphism mixed-uncertain Step 2 — result",
        "",
        f"**Frozen verdict: `{verdict}`.**",
        "",
        "## Decomposition",
        "",
        f"- all mixed_uncertain: **{all_mixed_n} / {len(discovery)} ({all_mixed_n/len(discovery):.3%})**",
        f"- ambiguous-palette mixed rows: **{ambiguous_all_n} ({ambiguous_all_n/all_mixed_n:.3%} of all mixed)**",
        f"- technical/other mixed rows: **{all_mixed_n-ambiguous_all_n}**",
        "",
        "## Independent reserve-defined bridge geometry",
        "",
        f"- valid ambiguous rows tested: **{ambiguous_valid_n}**",
        f"- reserve within-morph radius q95: **{radius95:.6f}**",
        f"- observed interior-bridge fraction: **{observed_bridge:.6f}**",
        f"- coordinate-permutation null mean: **{np.mean(coord_null):.6f}**",
        f"- null 2.5–97.5%: **[{np.quantile(coord_null,0.025):.6f}, {np.quantile(coord_null,0.975):.6f}]**",
        f"- one-sided p_bridge: **{p_bridge:.6g}**",
        "",
        "## Species-level relation to D",
        "",
        f"- rho(D, all mixed rate): **{diagnostics['rho_D_mixed_rate']:.6f}**",
        f"- rho(D, ambiguous-palette rate): **{diagnostics['rho_D_ambiguous_palette_rate']:.6f}**",
        f"- rho(D, bridge rate): **{rho_bridge:.6f}**, permutation p **{p_species:.6g}**",
        f"- rho(D, technical mixed rate): **{diagnostics['rho_D_technical_mixed_rate']:.6f}**",
        "",
        "## Decision",
        "",
        f"- bridge geometry enriched: **{geometry_pass}**",
        f"- positive D–bridge-rate link: **{species_pass}**",
        f"- intermediate-colour interpretation supported: **{survives}**",
        "- numeric claim that hard-state Simpson D is underestimated: **not licensed by this test**",
        "",
        "No new acquisition, threshold tuning, or reclassification was used.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
