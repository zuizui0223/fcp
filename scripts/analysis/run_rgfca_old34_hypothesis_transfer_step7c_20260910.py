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
STATES = ROOT / "results/rgfca_local_cooccurrence_segregation_step7a_20260910/photo_organization_states.csv"
OUT = ROOT / "results/rgfca_old34_hypothesis_transfer_step7c_20260910"
OUT.mkdir(parents=True, exist_ok=True)
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
R_EARTH_KM = 6371.0088
N_CROSS_PERM = 20_000
N_WITHIN_PERM = 999
SEED = 20260910


def bseries(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def haversine_pairs(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    i, j = np.triu_indices(len(lat), 1)
    la1 = np.radians(lat[i]); la2 = np.radians(lat[j])
    dla = la2 - la1; dlo = np.radians(lon[j] - lon[i])
    a = np.sin(dla / 2.0) ** 2 + np.cos(la1) * np.cos(la2) * np.sin(dlo / 2.0) ** 2
    d = 2.0 * R_EARTH_KM * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    return i, j, d


def rankdata(x: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(x, float)).rank(method="average").to_numpy(float)


def residualize_rank(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    ry = rankdata(y); rx = rankdata(x)
    X = np.column_stack([np.ones(len(rx)), rx])
    return ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(x); ry = rankdata(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def median_contrast_perm(values: np.ndarray, labels: np.ndarray, *, alternative: str, seed: int) -> dict[str, Any]:
    values = np.asarray(values, float); labels = np.asarray(labels, bool)
    keep = np.isfinite(values)
    values = values[keep]; labels = labels[keep]
    if labels.sum() == 0 or (~labels).sum() == 0:
        return {"n": int(len(values)), "n_positive": int(labels.sum()), "delta": None, "p": 1.0}
    obs = float(np.median(values[labels]) - np.median(values[~labels]))
    rng = np.random.default_rng(seed)
    null = np.empty(N_CROSS_PERM, float)
    for b in range(N_CROSS_PERM):
        lp = labels[rng.permutation(len(labels))]
        null[b] = float(np.median(values[lp]) - np.median(values[~lp]))
    if alternative == "greater":
        p = float((1 + np.sum(null >= obs - 1e-15)) / (N_CROSS_PERM + 1))
    else:
        p = float((1 + np.sum(np.abs(null) >= abs(obs) - 1e-15)) / (N_CROSS_PERM + 1))
    return {
        "n": int(len(values)), "n_positive": int(labels.sum()), "delta": obs, "p": p,
        "null_q025": float(np.quantile(null, .025)), "null_q975": float(np.quantile(null, .975)),
    }


def holm_two(p1: float, p2: float) -> tuple[float, float]:
    p = np.array([p1, p2], float); order = np.argsort(p)
    adj_sorted = np.maximum.accumulate(np.minimum(1.0, p[order] * np.array([2.0, 1.0])))
    out = np.empty(2, float); out[order] = adj_sorted
    return float(out[0]), float(out[1])


def largest_component_fraction(n: int, ii: np.ndarray, jj: np.ndarray, dist: np.ndarray, radius: float) -> float:
    parent = np.arange(n, dtype=int); size = np.ones(n, dtype=int)
    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb: return
        if size[ra] < size[rb]: ra, rb = rb, ra
        parent[rb] = ra; size[ra] += size[rb]
    for a, b in zip(ii[dist <= radius], jj[dist <= radius]):
        union(int(a), int(b))
    counts: dict[int, int] = {}
    for k in range(n):
        r = find(k); counts[r] = counts.get(r, 0) + 1
    return float(max(counts.values()) / n)


def doy365(series: pd.Series) -> np.ndarray:
    d = pd.to_datetime(series, errors="coerce")
    out = np.full(len(d), np.nan, float)
    ok = ~d.isna()
    if ok.any():
        vals = d[ok]
        day = vals.dt.dayofyear.to_numpy(float)
        leap = vals.dt.is_leap_year.to_numpy(bool)
        month = vals.dt.month.to_numpy(int); daym = vals.dt.day.to_numpy(int)
        # Collapse Feb 29 onto day 59 and remove the leap-day offset thereafter.
        day[(leap) & (month == 2) & (daym == 29)] = 59.0
        day[(leap) & ((month > 2))] -= 1.0
        out[np.flatnonzero(ok)] = day
    return out


def circular_day_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    z = np.abs(a - b)
    return np.minimum(z, 365.0 - z)


def load_tranche(path: Path, tranche: str) -> pd.DataFrame:
    use = ["species", "morph", "global_classifiable", "latitude", "longitude", "observer_id", "observed_on"]
    d = pd.read_csv(path, usecols=use)
    d = d.loc[bseries(d["global_classifiable"]) & d["morph"].isin(MORPHS)].copy()
    d = d.dropna(subset=["latitude", "longitude"])
    d["tranche"] = tranche
    return d


def one_species(g: pd.DataFrame, state: pd.Series, rng: np.random.Generator) -> dict[str, Any]:
    lat = g.latitude.to_numpy(float); lon = g.longitude.to_numpy(float)
    morph = g.morph.astype(str).to_numpy(); obs = g.observer_id.fillna("").astype(str).to_numpy()
    ii, jj, dist = haversine_pairs(lat, lon)
    span = float(dist.max()) if len(dist) else 0.0
    frag_parts = [1.0 - largest_component_fraction(len(g), ii, jj, dist, km) for km in (100.0, 250.0, 500.0)]
    frag = float(np.mean(frag_parts))

    days = doy365(g.observed_on)
    valid = np.isfinite(days)
    local = (dist <= 100.0) & (obs[ii] != obs[jj]) & valid[ii] & valid[jj]
    same = local & (morph[ii] == morph[jj]); cross = local & (morph[ii] != morph[jj])
    phen_gate = bool(local.sum() >= 30 and same.sum() >= 5 and cross.sum() >= 5)
    phen_delta = np.nan; phen_p = np.nan
    if phen_gate:
        td = circular_day_distance(days[ii], days[jj])
        phen_delta = float(np.median(td[cross]) - np.median(td[same]))
        null = np.empty(N_WITHIN_PERM, float)
        for b in range(N_WITHIN_PERM):
            mp = morph[rng.permutation(len(morph))]
            sp = local & (mp[ii] == mp[jj]); cp = local & (mp[ii] != mp[jj])
            if sp.sum() < 5 or cp.sum() < 5:
                null[b] = np.nan
            else:
                null[b] = float(np.median(td[cp]) - np.median(td[sp]))
        null = null[np.isfinite(null)]
        phen_p = float((1 + np.sum(np.abs(null) >= abs(phen_delta) - 1e-15)) / (len(null) + 1)) if len(null) else np.nan

    return {
        "tranche": state["tranche"], "species": state["species"], "n_classifiable": int(len(g)),
        "D": float(state["D"]), "C_star": bool(state["C_star"]), "S_star": bool(state["S_star"]),
        "organization_state": state["organization_state"],
        "maximum_pairwise_span_km": span,
        "frag_100km": frag_parts[0], "frag_250km": frag_parts[1], "frag_500km": frag_parts[2], "frag_multiscale": frag,
        "valid_dated_photos": int(valid.sum()), "eligible_local_dated_pairs": int(local.sum()),
        "same_morph_local_dated_pairs": int(same.sum()), "cross_morph_local_dated_pairs": int(cross.sum()),
        "phenology_proxy_gate": phen_gate, "phenology_delta_days": phen_delta, "phenology_within_species_perm_p": phen_p,
    }


def run_one(path: Path, tranche: str, states: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    d = load_tranche(path, tranche)
    ss = states[states.tranche.eq(tranche)].copy().sort_values("species")
    rng = np.random.default_rng(seed)
    rows = []
    groups = {sp: g.reset_index(drop=True) for sp, g in d.groupby("species", sort=False)}
    for _, st in ss.iterrows():
        sp = st.species
        if sp not in groups:
            raise RuntimeError(f"{tranche}: missing photo group for {sp}")
        rows.append(one_species(groups[sp], st, rng))
    m = pd.DataFrame(rows)

    h2 = median_contrast_perm(m.frag_multiscale.to_numpy(float), m.S_star.to_numpy(bool), alternative="greater", seed=seed + 11)
    residual_frag = residualize_rank(m.frag_multiscale.to_numpy(float), np.log1p(m.maximum_pairwise_span_km.to_numpy(float)))
    h2_span = median_contrast_perm(residual_frag, m.S_star.to_numpy(bool), alternative="greater", seed=seed + 12)

    hp = m[m.phenology_proxy_gate].copy()
    h9 = median_contrast_perm(hp.phenology_delta_days.to_numpy(float), hp.C_star.to_numpy(bool), alternative="two-sided", seed=seed + 21)
    h9["eligible_species"] = int(len(hp))
    h9["rho_with_D"] = spearman(hp.phenology_delta_days.to_numpy(float), hp.D.to_numpy(float)) if len(hp) >= 3 else np.nan
    h9["rho_with_S_star"] = spearman(hp.phenology_delta_days.to_numpy(float), hp.S_star.astype(int).to_numpy(float)) if len(hp) >= 3 else np.nan

    h2_holm, h9_holm = holm_two(h2["p"], h9["p"])
    h2["holm_p"] = h2_holm; h9["holm_p"] = h9_holm
    summary = {
        "species": int(len(m)),
        "H2_fragmentation": h2,
        "H2_span_adjusted_robustness": h2_span,
        "H9_photo_date_partitioning": h9,
        "H2_primary_pass": bool(h2["delta"] is not None and h2["delta"] > 0 and h2_holm <= .05),
        "H9_primary_pass": bool(h9["delta"] is not None and h9_holm <= .05),
        "phenology_proxy_gate_fraction": float(m.phenology_proxy_gate.mean()),
    }
    return m, summary


def main() -> None:
    states = pd.read_csv(STATES)
    states["C_star"] = bseries(states["C_star"]); states["S_star"] = bseries(states["S_star"])
    dm, ds = run_one(DISC, "discovery", states, SEED + 1)
    rm, rs = run_one(RES, "reserve", states, SEED + 2)
    metrics = pd.concat([dm, rm], ignore_index=True)
    metrics.to_csv(OUT / "species_transfer_metrics.csv", index=False)

    recurrent = {
        "H2": bool(ds["H2_primary_pass"] and rs["H2_primary_pass"]),
        "H9": bool(ds["H9_primary_pass"] and rs["H9_primary_pass"]),
    }
    result = {
        "analysis": "rgfca_old34_hypothesis_transfer_step7c_photo_only",
        "protocol": "docs/RGFCA_OLD34_HYPOTHESIS_TRANSFER_STEP7C_PROTOCOL_20260910.md",
        "discovery": ds,
        "reserve": rs,
        "two_tranche_recurrent": recurrent,
        "environmental_family_status": "NOT_OPENED_REQUIRES_SEPARATE_SOURCE_FREEZE",
        "claim_boundary": "H2 is sampled-coordinate fragmentation, not biological range fragmentation. H9 is photo-date partitioning, not direct flowering phenology. C*/S* are photo-derived states, so this is a transferred test rather than a confirmatory replication of the historical 34-species response.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = ["# RGFCA Step 7C — historical hypothesis transfer, photo-only family", ""]
    for tr, s in [("Discovery", ds), ("Reserve", rs)]:
        h2=s["H2_fragmentation"]; h2s=s["H2_span_adjusted_robustness"]; h9=s["H9_photo_date_partitioning"]
        lines += [f"## {tr}", "",
                  f"- eligible C*/S* species: **{s['species']}**",
                  f"- H2 fragmentation S*=1 minus S*=0 median contrast: **{h2['delta']:.6f}**, raw p **{h2['p']:.6g}**, Holm p **{h2['holm_p']:.6g}**",
                  f"- H2 span-adjusted robustness contrast: **{h2s['delta']:.6f}**, p **{h2s['p']:.6g}**",
                  f"- H9 eligible species: **{h9['eligible_species']}**",
                  f"- H9 C*=1 minus C*=0 phenology-delta contrast: **{h9['delta']:.3f} days**, raw p **{h9['p']:.6g}**, Holm p **{h9['holm_p']:.6g}**",
                  f"- H9 rho(phenology proxy, D): **{h9['rho_with_D']:.4f}**; rho(proxy, S*): **{h9['rho_with_S_star']:.4f}**", ""]
    lines += ["## Two-tranche recurrence", "", f"- H2: **{recurrent['H2']}**", f"- H9: **{recurrent['H9']}**", "",
              "H0/H1a/H1b/H3a/H3b/H3c/H8 remain unopened until the RGFCA environmental source is frozen. H2 is sampled-coordinate fragmentation; H9 is photo-date partitioning rather than direct flowering phenology."]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()
