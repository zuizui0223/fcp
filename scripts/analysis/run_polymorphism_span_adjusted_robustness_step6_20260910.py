#!/usr/bin/env python3
"""Run frozen post-outcome Step 6 sampled-span robustness diagnostics."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_span_adjusted_robustness_step6_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISCOVERY_SPATIAL = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909" / "species_membership_and_observed_rho.csv"
PANEL = ROOT / "results" / "polymorphism_species_attributes_step4_preflight_20260910" / "covariate_panel_preoutcome.csv"
RESERVE_SPATIAL = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "reserve_species_membership_and_observed_metrics.csv"
RESERVE_GEOMETRY = ROOT / "data" / "frozen" / "rgfca_reserve_geometry_audit_v1.csv"
PROTOCOL = "docs/POLYMORPHISM_SPAN_ADJUSTED_ROBUSTNESS_STEP6_PROTOCOL_20260910.md"

SEED = 20260910 + 600
N_PERM = 20_000


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def rank(x: np.ndarray) -> np.ndarray:
    return stats.rankdata(np.asarray(x, float)).astype(float)


def corr(a: np.ndarray, b: np.ndarray) -> float:
    aa = np.asarray(a, float) - float(np.mean(a))
    bb = np.asarray(b, float) - float(np.mean(b))
    denom = float(np.linalg.norm(aa) * np.linalg.norm(bb))
    if denom <= 0:
        return math.nan
    return float(aa @ bb / denom)


def design_controls(*controls: np.ndarray) -> np.ndarray:
    cols = [np.ones(len(controls[0]), float)]
    for c in controls:
        r = rank(np.asarray(c, float))
        r = (r - r.mean()) / r.std(ddof=0)
        cols.append(r)
    return np.column_stack(cols)


def residualize(v: np.ndarray, Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q, _ = np.linalg.qr(Z, mode="reduced")
    fit = q @ (q.T @ np.asarray(v, float))
    return np.asarray(v, float) - fit, fit


def partial_rank_freedman_lane(
    y: np.ndarray,
    x: np.ndarray,
    controls: list[np.ndarray],
    *,
    seed: int,
) -> dict[str, Any]:
    arrays = [np.asarray(y, float), np.asarray(x, float)] + [np.asarray(c, float) for c in controls]
    keep = np.logical_and.reduce([np.isfinite(a) for a in arrays])
    yy = arrays[0][keep]
    xx = arrays[1][keep]
    cc = [a[keep] for a in arrays[2:]]
    ry = rank(yy)
    rx = rank(xx)
    Z = design_controls(*cc)
    ey, fy = residualize(ry, Z)
    ex, _ = residualize(rx, Z)
    observed = corr(ey, ex)

    q, _ = np.linalg.qr(Z, mode="reduced")
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    batch = 400
    for start in range(0, N_PERM, batch):
        end = min(N_PERM, start + batch)
        perms = np.stack([rng.permutation(len(ey)) for _ in range(end - start)], axis=0)
        ystar = fy[None, :] + ey[perms]
        estar = ystar - (ystar @ q) @ q.T
        numer = estar @ ex
        denom = np.linalg.norm(estar, axis=1) * np.linalg.norm(ex)
        null[start:end] = numer / denom
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(yy)),
        "partial_rho": observed,
        "p_two_sided": p,
        "null_mean": float(np.mean(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def ancova_indicator_freedman_lane(
    y: np.ndarray,
    indicator: np.ndarray,
    controls: list[np.ndarray],
    *,
    seed: int,
) -> dict[str, Any]:
    arrays = [np.asarray(y, float), np.asarray(indicator, float)] + [np.asarray(c, float) for c in controls]
    keep = np.logical_and.reduce([np.isfinite(a) for a in arrays])
    yy = arrays[0][keep]
    ii = arrays[1][keep]
    cc = [a[keep] for a in arrays[2:]]
    Z = design_controls(*cc)
    X = np.column_stack([Z, ii])
    pinv_full = np.linalg.pinv(X)
    beta = np.linalg.lstsq(X, yy, rcond=None)[0]
    observed = float(beta[-1])
    ey, fy = residualize(yy, Z)

    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    batch = 400
    for start in range(0, N_PERM, batch):
        end = min(N_PERM, start + batch)
        perms = np.stack([rng.permutation(len(ey)) for _ in range(end - start)], axis=0)
        ystar = fy[None, :] + ey[perms]
        betas = ystar @ pinv_full.T
        null[start:end] = betas[:, -1]
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(yy)),
        "n_indicator_1": int(np.sum(ii == 1)),
        "n_indicator_0": int(np.sum(ii == 0)),
        "adjusted_raw_spatial_rho_difference": observed,
        "p_two_sided": p,
        "null_mean": float(np.mean(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def permutation_spearman_two_sided(x: np.ndarray, y: np.ndarray, *, seed: int) -> dict[str, Any]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    rx = rank(x[keep])
    ry = rank(y[keep])
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    observed = float(rx @ ry / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = float(rx[rng.permutation(len(rx))] @ ry / denom)
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(keep.sum()),
        "rho": observed,
        "p_two_sided": p,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def repeated_group_pairs(labels: pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    lab = labels.fillna("").astype(str).to_numpy()
    groups: dict[str, np.ndarray] = {}
    for name in sorted(set(lab)):
        if not name:
            continue
        idx = np.flatnonzero(lab == name)
        if len(idx) >= 2:
            groups[name] = idx
    repeated_species = sorted({int(i) for idx in groups.values() for i in idx})
    remap = {old: new for new, old in enumerate(repeated_species)}
    pair_i: list[int] = []
    pair_j: list[int] = []
    weights: list[float] = []
    ng = len(groups)
    for idx in groups.values():
        n_pairs = len(idx) * (len(idx) - 1) // 2
        w = 1.0 / (ng * n_pairs) if ng else 0.0
        for a in range(len(idx) - 1):
            for b in range(a + 1, len(idx)):
                pair_i.append(remap[int(idx[a])])
                pair_j.append(remap[int(idx[b])])
                weights.append(w)
    meta = {
        "repeated_groups": int(ng),
        "species_in_repeated_groups": int(len(repeated_species)),
        "group_sizes": {k: int(len(v)) for k, v in groups.items()},
        "gate_pass": bool(ng >= 10 and len(repeated_species) >= 30),
    }
    return (
        np.asarray(repeated_species, int),
        np.asarray(pair_i, int),
        np.asarray(pair_j, int),
        np.asarray(weights, float),
        meta,
    )


def taxonomic_clustering_test(values: np.ndarray, labels: pd.Series, *, seed: int) -> dict[str, Any]:
    species_idx, pair_i, pair_j, weights, meta = repeated_group_pairs(labels)
    out: dict[str, Any] = dict(meta)
    if not out["gate_pass"]:
        out.update({"tested": False, "W_observed": None, "null_mean_W": None, "clustering_gain": None, "p_lower": 1.0})
        return out
    d = np.asarray(values, float)[species_idx]
    if not np.isfinite(d).all():
        raise RuntimeError("non-finite values inside repeated taxonomic groups")
    if len(weights) == 0 or not np.isclose(weights.sum(), 1.0):
        raise RuntimeError("invalid equal-group clustering weights")
    observed = float(np.sum(np.abs(d[pair_i] - d[pair_j]) * weights))
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for start in range(0, N_PERM, 400):
        end = min(N_PERM, start + 400)
        perms = np.stack([rng.permutation(len(d)) for _ in range(end - start)], axis=0)
        vals = d[perms]
        null[start:end] = np.sum(np.abs(vals[:, pair_i] - vals[:, pair_j]) * weights[None, :], axis=1)
    p = float((1 + np.sum(null <= observed + 1e-15)) / (N_PERM + 1))
    null_mean = float(null.mean())
    gain = float(1.0 - observed / null_mean) if null_mean > 0 else math.nan
    out.update({
        "tested": True,
        "W_observed": observed,
        "null_mean_W": null_mean,
        "clustering_gain": gain,
        "p_lower": p,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    })
    return out


def residualize_raw_on_rank_controls(y: np.ndarray, *controls: np.ndarray) -> np.ndarray:
    Z = design_controls(*controls)
    e, _ = residualize(np.asarray(y, float), Z)
    return e


def holm_adjust(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    sorted_p = p[order]
    raw = np.array([(m - i) * sorted_p[i] for i in range(m)], float)
    adj_sorted = np.minimum(1.0, np.maximum.accumulate(raw))
    out = np.empty(m, float)
    out[order] = adj_sorted
    return [float(v) for v in out]


def main() -> None:
    disc = pd.read_csv(DISCOVERY_SPATIAL)
    panel = pd.read_csv(PANEL)
    res = pd.read_csv(RESERVE_SPATIAL)
    geom = pd.read_csv(RESERVE_GEOMETRY)

    if len(disc) != 369 or disc["species"].nunique() != 369:
        raise RuntimeError("discovery Step 5 fingerprint mismatch")
    if len(panel) != 369 or panel["species"].nunique() != 369:
        raise RuntimeError("Step 4 pre-outcome panel fingerprint mismatch")
    if len(res) != 363 or res["species"].nunique() != 363:
        raise RuntimeError("reserve Step 5b fingerprint mismatch")
    if int(disc["primary_second_ge_0_10"].sum()) != 172:
        raise RuntimeError("discovery >=10% subset fingerprint mismatch")
    if int(res["primary_second_ge_0_10"].sum()) != 151:
        raise RuntimeError("reserve >=10% subset fingerprint mismatch")

    raw_disc_rho = float(stats.spearmanr(disc["D"], disc["spatial_observed_rho"]).statistic)
    raw_res_rho = float(stats.spearmanr(res["D"], res["rho_primary"]).statistic)
    if abs(raw_disc_rho - 0.089213) > 5e-6:
        raise RuntimeError(f"discovery D-spatial fingerprint changed: {raw_disc_rho}")
    if abs(raw_res_rho - 0.101601) > 5e-6:
        raise RuntimeError(f"reserve D-spatial fingerprint changed: {raw_res_rho}")

    discj = disc.merge(
        panel[["species", "genus", "family", "log1p_span_primary", "n_classifiable"]],
        on="species",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_panel"),
    )
    if len(discj) != 369:
        raise RuntimeError("discovery span join lost species")
    if not np.array_equal(discj["n_classifiable"].astype(int), discj["n_classifiable_panel"].astype(int)):
        raise RuntimeError("discovery n_classifiable disagreement")

    geom_use = geom[["inat_taxon_id", "species", "maximum_span_km"]].copy()
    resj = res.merge(
        geom_use,
        on="inat_taxon_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_geometry"),
    )
    if len(resj) != 363:
        raise RuntimeError(f"reserve geometry join lost species: {len(resj)}")
    if not (resj["species"] == resj["species_geometry"]).all():
        raise RuntimeError("reserve geometry species-name disagreement")
    if not resj["maximum_span_km"].gt(0).all():
        raise RuntimeError("reserve maximum span contains non-positive values")
    resj["log1p_span_reserve"] = np.log1p(resj["maximum_span_km"].astype(float))
    resj["genus"] = resj["species"].astype(str).str.split().str[0]

    # ---- Primary Claim-2 continuous span-adjusted family ----
    disc_cont = partial_rank_freedman_lane(
        discj["spatial_observed_rho"].to_numpy(float),
        discj["D"].to_numpy(float),
        [discj["log1p_span_primary"].to_numpy(float), discj["n_classifiable"].to_numpy(float)],
        seed=SEED + 1,
    )
    res_cont = partial_rank_freedman_lane(
        resj["rho_primary"].to_numpy(float),
        resj["D"].to_numpy(float),
        [resj["log1p_span_reserve"].to_numpy(float), resj["n_classifiable"].to_numpy(float)],
        seed=SEED + 2,
    )
    cont_adj = holm_adjust([disc_cont["p_two_sided"], res_cont["p_two_sided"]])
    disc_cont["holm_p_two_tranches"] = cont_adj[0]
    res_cont["holm_p_two_tranches"] = cont_adj[1]
    spatial_span_robust = bool(
        disc_cont["partial_rho"] > 0
        and res_cont["partial_rho"] > 0
        and disc_cont["holm_p_two_tranches"] < 0.05
        and res_cont["holm_p_two_tranches"] < 0.05
    )

    # ---- Fixed >=10% threshold sensitivity ----
    disc_thr = ancova_indicator_freedman_lane(
        discj["spatial_observed_rho"].to_numpy(float),
        discj["primary_second_ge_0_10"].astype(int).to_numpy(),
        [discj["log1p_span_primary"].to_numpy(float), discj["n_classifiable"].to_numpy(float)],
        seed=SEED + 3,
    )
    res_thr = ancova_indicator_freedman_lane(
        resj["rho_primary"].to_numpy(float),
        resj["primary_second_ge_0_10"].astype(int).to_numpy(),
        [resj["log1p_span_reserve"].to_numpy(float), resj["n_classifiable"].to_numpy(float)],
        seed=SEED + 4,
    )
    thr_adj = holm_adjust([disc_thr["p_two_sided"], res_thr["p_two_sided"]])
    disc_thr["holm_p_two_tranches"] = thr_adj[0]
    res_thr["holm_p_two_tranches"] = thr_adj[1]
    threshold_span_robust = bool(
        disc_thr["adjusted_raw_spatial_rho_difference"] > 0
        and res_thr["adjusted_raw_spatial_rho_difference"] > 0
        and disc_thr["holm_p_two_tranches"] < 0.05
        and res_thr["holm_p_two_tranches"] < 0.05
    )

    # ---- Reserve-photo recurrence of the two Step-4 supported attributes ----
    reserve_span_D = permutation_spearman_two_sided(
        resj["D"].to_numpy(float), resj["log1p_span_reserve"].to_numpy(float), seed=SEED + 5
    )
    reserve_genus_raw = taxonomic_clustering_test(
        resj["D"].to_numpy(float), resj["genus"], seed=SEED + 6
    )
    recurrence_adj = holm_adjust([reserve_span_D["p_two_sided"], reserve_genus_raw["p_lower"]])
    reserve_span_D["holm_p_two_attribute_recurrence"] = recurrence_adj[0]
    reserve_genus_raw["holm_p_two_attribute_recurrence"] = recurrence_adj[1]
    recurrent_attributes = []
    if reserve_span_D["rho"] > 0 and reserve_span_D["holm_p_two_attribute_recurrence"] < 0.05:
        recurrent_attributes.append("sampled_geographic_span")
    if (
        reserve_genus_raw.get("tested")
        and reserve_genus_raw.get("clustering_gain") is not None
        and reserve_genus_raw["clustering_gain"] > 0
        and reserve_genus_raw["holm_p_two_attribute_recurrence"] < 0.05
    ):
        recurrent_attributes.append("genus_taxonomic_clustering")

    # ---- Genus clustering after span+n removal ----
    disc_D = discj["D"].to_numpy(float)
    disc_n = discj["n_classifiable"].to_numpy(float)
    disc_Dunb = disc_D * disc_n / (disc_n - 1.0)
    disc_resid = residualize_raw_on_rank_controls(
        disc_D, discj["log1p_span_primary"].to_numpy(float), disc_n
    )
    disc_resid_unb = residualize_raw_on_rank_controls(
        disc_Dunb, discj["log1p_span_primary"].to_numpy(float), disc_n
    )
    disc_genus_adj = taxonomic_clustering_test(disc_resid, discj["genus"], seed=SEED + 7)
    disc_genus_adj_unb = taxonomic_clustering_test(disc_resid_unb, discj["genus"], seed=SEED + 8)

    res_D = resj["D"].to_numpy(float)
    res_n = resj["n_classifiable"].to_numpy(float)
    res_Dunb = res_D * res_n / (res_n - 1.0)
    res_resid = residualize_raw_on_rank_controls(
        res_D, resj["log1p_span_reserve"].to_numpy(float), res_n
    )
    res_resid_unb = residualize_raw_on_rank_controls(
        res_Dunb, resj["log1p_span_reserve"].to_numpy(float), res_n
    )
    res_genus_adj = taxonomic_clustering_test(res_resid, resj["genus"], seed=SEED + 9)
    res_genus_adj_unb = taxonomic_clustering_test(res_resid_unb, resj["genus"], seed=SEED + 10)

    overlap = sorted(set(discj["species"].astype(str)) & set(resj["species"].astype(str)))
    descriptives = {
        "discovery_reserve_species_overlap": int(len(overlap)),
        "discovery_species": int(len(discj)),
        "reserve_species": int(len(resj)),
        "discovery_raw_rho_D_spatial": raw_disc_rho,
        "reserve_raw_rho_D_spatial": raw_res_rho,
        "discovery_rho_span_spatial": float(stats.spearmanr(discj["log1p_span_primary"], discj["spatial_observed_rho"]).statistic),
        "reserve_rho_span_spatial": float(stats.spearmanr(resj["log1p_span_reserve"], resj["rho_primary"]).statistic),
        "discovery_rho_n_spatial": float(stats.spearmanr(discj["n_classifiable"], discj["spatial_observed_rho"]).statistic),
        "reserve_rho_n_spatial": float(stats.spearmanr(resj["n_classifiable"], resj["rho_primary"]).statistic),
    }

    if spatial_span_robust:
        spatial_verdict = "SPAN_ROBUST_CONTINUOUS"
    elif disc_cont["partial_rho"] > 0 and res_cont["partial_rho"] > 0:
        spatial_verdict = "POSITIVE_BUT_NOT_TWO_TRANCHE_HOLM_SUPPORTED"
    else:
        spatial_verdict = "SPAN_ADJUSTMENT_WEAKENS_OR_REVERSES_CONTINUOUS_PATTERN"

    result = {
        "analysis": "polymorphism_span_adjusted_robustness_step6",
        "date_jst": "2026-09-10",
        "protocol": PROTOCOL,
        "post_outcome_robustness_diagnostic": True,
        "new_image_acquisition": False,
        "n_permutations": N_PERM,
        "descriptives": descriptives,
        "claim2_continuous_span_adjusted": {
            "discovery": disc_cont,
            "reserve": res_cont,
            "two_tranche_gate_pass": spatial_span_robust,
            "verdict": spatial_verdict,
        },
        "claim2_threshold_span_adjusted": {
            "discovery": disc_thr,
            "reserve": res_thr,
            "two_tranche_gate_pass": threshold_span_robust,
        },
        "step4_reserve_photo_recurrence": {
            "sampled_geographic_span": reserve_span_D,
            "genus_taxonomic_clustering": reserve_genus_raw,
            "holm_supported_attributes": recurrent_attributes,
        },
        "genus_span_adjusted_diagnostic": {
            "discovery_D": disc_genus_adj,
            "discovery_D_unbiased": disc_genus_adj_unb,
            "reserve_D": res_genus_adj,
            "reserve_D_unbiased": res_genus_adj_unb,
        },
        "interpretation_ceiling": (
            "A surviving adjusted effect excludes the measured sampled-span and classifiable-n proxies "
            "as sufficient explanations; it does not establish adaptation, causal range-size effects, "
            "population-genetic differentiation, environmental selection, or a shared global boundary."
        ),
    }

    disc_out = discj[[
        "species", "D", "second_fraction", "primary_second_ge_0_10", "spatial_observed_rho",
        "log1p_span_primary", "n_classifiable", "genus", "family"
    ]].copy()
    res_out = resj[[
        "inat_taxon_id", "species", "D", "second_fraction", "primary_second_ge_0_10", "rho_primary",
        "log1p_span_reserve", "n_classifiable", "genus"
    ]].copy()
    disc_out.to_csv(OUT / "discovery_joined.csv", index=False)
    res_out.to_csv(OUT / "reserve_joined.csv", index=False)
    pd.DataFrame({"species": overlap}).to_csv(OUT / "discovery_reserve_species_overlap.csv", index=False)
    write_json(OUT / "result.json", result)

    lines = [
        "# Polymorphism Step 6 — span-adjusted robustness result",
        "",
        f"**Spatial verdict: `{spatial_verdict}`.**",
        "",
        "## Claim 2: D–spatial organization after sampled-span + n adjustment",
        "",
        f"- discovery: raw rho = **{raw_disc_rho:.6f}**; adjusted partial rho = **{disc_cont['partial_rho']:.6f}**; raw permutation p = **{disc_cont['p_two_sided']:.6g}**; Holm p = **{disc_cont['holm_p_two_tranches']:.6g}**",
        f"- reserve: raw rho = **{raw_res_rho:.6f}**; adjusted partial rho = **{res_cont['partial_rho']:.6f}**; raw permutation p = **{res_cont['p_two_sided']:.6g}**; Holm p = **{res_cont['holm_p_two_tranches']:.6g}**",
        f"- two-tranche continuous gate pass: **{spatial_span_robust}**",
        "",
        "## >=10% second-morph threshold sensitivity",
        "",
        f"- discovery adjusted spatial-rho difference = **{disc_thr['adjusted_raw_spatial_rho_difference']:.6f}**; Holm p = **{disc_thr['holm_p_two_tranches']:.6g}**",
        f"- reserve adjusted spatial-rho difference = **{res_thr['adjusted_raw_spatial_rho_difference']:.6f}**; Holm p = **{res_thr['holm_p_two_tranches']:.6g}**",
        f"- two-tranche threshold gate pass: **{threshold_span_robust}**",
        "",
        "## Step 4 signals in the reserve photo tranche",
        "",
        f"- D vs reserve sampled span: rho = **{reserve_span_D['rho']:.6f}**; raw p = **{reserve_span_D['p_two_sided']:.6g}**; Holm p = **{reserve_span_D['holm_p_two_attribute_recurrence']:.6g}**",
        f"- reserve genus clustering: groups = **{reserve_genus_raw['repeated_groups']}**, species = **{reserve_genus_raw['species_in_repeated_groups']}**, gain = **{reserve_genus_raw.get('clustering_gain', math.nan):.6f}**; raw p = **{reserve_genus_raw['p_lower']:.6g}**; Holm p = **{reserve_genus_raw['holm_p_two_attribute_recurrence']:.6g}**",
        f"- Holm-supported reserve-photo recurrent attributes: **{recurrent_attributes}**",
        "",
        "## Genus clustering after span + n removal",
        "",
        f"- discovery residual-D clustering gain = **{disc_genus_adj.get('clustering_gain', math.nan):.6f}**, p = **{disc_genus_adj['p_lower']:.6g}**; D_unbiased p = **{disc_genus_adj_unb['p_lower']:.6g}**",
        f"- reserve residual-D clustering gain = **{res_genus_adj.get('clustering_gain', math.nan):.6f}**, p = **{res_genus_adj['p_lower']:.6g}**; D_unbiased p = **{res_genus_adj_unb['p_lower']:.6g}**",
        "",
        "## Opportunity descriptives",
        "",
        f"- discovery/reserve species overlap: **{len(overlap)}**",
        f"- rho(span, spatial rho): discovery **{descriptives['discovery_rho_span_spatial']:.6f}**, reserve **{descriptives['reserve_rho_span_spatial']:.6f}**",
        f"- rho(n_classifiable, spatial rho): discovery **{descriptives['discovery_rho_n_spatial']:.6f}**, reserve **{descriptives['reserve_rho_n_spatial']:.6f}**",
        "",
        "This is a post-outcome robustness diagnostic. A surviving result rules out the measured sampled-span and classifiable-n proxies as sufficient explanations; it does not imply causality, adaptation, genetic differentiation, or a shared global boundary.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
