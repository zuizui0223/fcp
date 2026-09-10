#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DISC = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
RES = ROOT / "data/derived/rgfca_reserve_replication_measured_photos_v1.csv"
LIT = ROOT / "data/frozen/frozen_34species_coexistence_segregation_v22.csv"
OUT = ROOT / "results/rgfca_local_cooccurrence_segregation_step7a_20260910"
OUT.mkdir(parents=True, exist_ok=True)
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
RADII = [100.0, 250.0, 500.0]
N_PERM = 999
SEED = 20260910


def bseries(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.lower().isin(["true", "1", "yes"])


def haversine_pairs(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    i, j = np.triu_indices(len(lat), 1)
    la1 = np.radians(lat[i]); la2 = np.radians(lat[j])
    dlo = np.radians(lon[j] - lon[i]); dla = la2 - la1
    a = np.sin(dla / 2.0) ** 2 + np.cos(la1) * np.cos(la2) * np.sin(dlo / 2.0) ** 2
    d = 6371.0088 * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    return i, j, d


def load(path: Path) -> pd.DataFrame:
    use = ["species", "morph", "global_classifiable", "latitude", "longitude", "observer_id", "photo_id"]
    df = pd.read_csv(path, usecols=use)
    keep = bseries(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    df = df.loc[keep].dropna(subset=["latitude", "longitude"]).copy()
    return df


def one_species(g: pd.DataFrame, rng: np.random.Generator) -> dict[str, Any] | None:
    n = len(g)
    if n < 40:
        return None
    counts = g["morph"].value_counts()
    p = np.array([counts.get(m, 0) / n for m in MORPHS], float)
    second = float(np.sort(p)[-2])
    if second < 0.10:
        return None
    D = float(1.0 - np.sum(p * p))
    lat = g["latitude"].to_numpy(float); lon = g["longitude"].to_numpy(float)
    morph = g["morph"].astype(str).to_numpy(); obs = g["observer_id"].astype(str).to_numpy()
    pid = g["photo_id"].astype(str).to_numpy()
    ii, jj, dist = haversine_pairs(lat, lon)
    diff = morph[ii] != morph[jj]
    same = ~diff
    if diff.sum() == 0 or same.sum() == 0:
        return None
    global_cross_rate = float(diff.mean())

    row: dict[str, Any] = {"species": str(g["species"].iloc[0]), "n": n, "D": D, "second_fraction": second}
    for r in RADII:
        close = dist <= r
        od = obs[ii] != obs[jj]
        cod = close & od
        cross = cod & diff
        endpoints = np.unique(np.concatenate([ii[cross], jj[cross]])) if cross.any() else np.array([], int)
        cross_observers = np.unique(np.concatenate([obs[ii[cross]], obs[jj[cross]]])) if cross.any() else np.array([], str)
        close_n = int(cod.sum()); cross_n = int(cross.sum())
        local_rate = float(cross_n / close_n) if close_n else np.nan
        row.update({
            f"close_observer_disjoint_pairs_{int(r)}km": close_n,
            f"cross_morph_close_pairs_{int(r)}km": cross_n,
            f"cross_morph_unique_photos_{int(r)}km": int(len(endpoints)),
            f"cross_morph_unique_observers_{int(r)}km": int(len(cross_observers)),
            f"local_cross_rate_{int(r)}km": local_rate,
            f"local_retention_ratio_{int(r)}km": float(local_rate / global_cross_rate) if close_n and global_cross_rate > 0 else np.nan,
        })
    C = bool(
        row["close_observer_disjoint_pairs_100km"] >= 30
        and row["cross_morph_close_pairs_100km"] >= 5
        and row["cross_morph_unique_photos_100km"] >= 4
        and row["cross_morph_unique_observers_100km"] >= 3
    )

    obs_delta = float(np.median(dist[diff]) - np.median(dist[same]))
    null = np.empty(N_PERM, float)
    for b in range(N_PERM):
        mp = morph[rng.permutation(n)]
        dp = mp[ii] != mp[jj]
        if dp.all() or (~dp).all():
            null[b] = np.nan
        else:
            null[b] = float(np.median(dist[dp]) - np.median(dist[~dp]))
    null = null[np.isfinite(null)]
    p_upper = float((1 + np.sum(null >= obs_delta - 1e-12)) / (len(null) + 1))
    S = bool(obs_delta > 0 and p_upper <= 0.05)
    if C and S:
        state = "cooccurrence_and_segregation"
    elif C:
        state = "local_cooccurrence_only"
    elif S:
        state = "spatial_segregation_only"
    else:
        state = "unresolved_photo_state"
    row.update({
        "global_cross_morph_pair_rate": global_cross_rate,
        "C_star": C,
        "S_star": S,
        "segregation_delta_km": obs_delta,
        "segregation_perm_p": p_upper,
        "organization_state": state,
        "second_ge_0_20": bool(second >= 0.20),
    })
    return row


def analyze(df: pd.DataFrame, name: str, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for _, g in df.groupby("species", sort=True):
        r = one_species(g.reset_index(drop=True), rng)
        if r is not None:
            rows.append(r)
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    out.insert(0, "tranche", name)
    return out


def lit_validate(all_states: pd.DataFrame) -> dict[str, Any]:
    lit = pd.read_csv(LIT)
    merged = all_states.merge(lit, left_on="species", right_on="canonical_name", how="inner")
    if merged.empty:
        return {"overlap_rows": 0}
    c_acc = float((merged["C_star"].astype(int) == merged["C_local_coexistence_documented"].astype(int)).mean())
    s_acc = float((merged["S_star"].astype(int) == merged["S_spatial_segregation_documented"].astype(int)).mean())
    return {
        "overlap_rows": int(len(merged)),
        "overlap_species": int(merged["species"].nunique()),
        "C_binary_accuracy": c_acc,
        "S_binary_accuracy": s_acc,
        "rows": merged[["tranche", "species", "organization_state", "organization_state_y", "C_star", "S_star", "C_local_coexistence_documented", "S_spatial_segregation_documented"]].to_dict("records") if "organization_state_y" in merged.columns else [],
    }


def summary(df: pd.DataFrame) -> dict[str, Any]:
    return {
        "eligible_polymorphic_species": int(len(df)),
        "state_counts": {str(k): int(v) for k, v in df["organization_state"].value_counts().to_dict().items()},
        "C_star_fraction": float(df["C_star"].mean()) if len(df) else None,
        "S_star_fraction": float(df["S_star"].mean()) if len(df) else None,
        "median_local_retention_100km": float(df["local_retention_ratio_100km"].median()) if len(df) else None,
        "median_segregation_delta_km": float(df["segregation_delta_km"].median()) if len(df) else None,
    }


def main() -> None:
    d = analyze(load(DISC), "discovery", SEED + 1)
    r = analyze(load(RES), "reserve", SEED + 2)
    all_states = pd.concat([d, r], ignore_index=True)
    all_states.to_csv(OUT / "photo_organization_states.csv", index=False)
    validation = lit_validate(all_states)
    result = {
        "analysis": "rgfca_local_cooccurrence_segregation_step7a",
        "protocol": "docs/RGFCA_LOCAL_COOCCURRENCE_SEGREGATION_STEP7A_PROTOCOL_20260910.md",
        "discovery": summary(d),
        "reserve": summary(r),
        "species_overlap_discovery_reserve": int(len(set(d.species) & set(r.species))),
        "literature_validation": validation,
        "claim_boundary": "C* is repeated <=100-km photo co-occurrence, not within-population coexistence; S* is photo-label spatial segregation, not adaptation/genetic differentiation.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# RGFCA Step 7A — local co-occurrence and spatial segregation",
        "",
        f"- discovery eligible polymorphic species: **{len(d)}**",
        f"- reserve eligible polymorphic species: **{len(r)}**",
        f"- discovery/reserve species overlap: **{result['species_overlap_discovery_reserve']}**",
        "",
        "## Discovery states",
        "",
    ]
    for k, v in result["discovery"]["state_counts"].items(): lines.append(f"- {k}: **{v}**")
    lines += ["", "## Reserve states", ""]
    for k, v in result["reserve"]["state_counts"].items(): lines.append(f"- {k}: **{v}**")
    lines += [
        "", "## Continuous descriptors", "",
        f"- discovery median local-retention ratio (100 km): **{result['discovery']['median_local_retention_100km']:.4f}**",
        f"- reserve median local-retention ratio (100 km): **{result['reserve']['median_local_retention_100km']:.4f}**",
        f"- discovery median segregation delta: **{result['discovery']['median_segregation_delta_km']:.2f} km**",
        f"- reserve median segregation delta: **{result['reserve']['median_segregation_delta_km']:.2f} km**",
        "", "## Literature validation", "",
        f"- overlap rows: **{validation.get('overlap_rows', 0)}**",
    ]
    if validation.get("overlap_rows", 0):
        lines += [f"- C binary accuracy: **{validation['C_binary_accuracy']:.3f}**", f"- S binary accuracy: **{validation['S_binary_accuracy']:.3f}**"]
    lines += ["", "C* is repeated <=100-km photo co-occurrence, not proof of within-population coexistence. S* is photo-label spatial segregation, not evidence of adaptation or genetic differentiation."]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()
