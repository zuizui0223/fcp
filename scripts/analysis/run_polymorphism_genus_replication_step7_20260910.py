#!/usr/bin/env python3
"""Replicate genus-level clustering of flower-colour polymorphism in the species-disjoint reserve."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_genus_replication_step7_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISC_PATH = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
DISC_META_PATH = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "result.json"
RES_PATH = ROOT / "results" / "polymorphism_spatial_span_adjusted_step6_20260910" / "reserve_species_span_adjusted_input.csv"

N_PERM = 20_000
STEP4_SEED = 20260914
STEP7_SEED = 20260917
EXPECTED_DISCOVERY_SPECIES = 369
EXPECTED_RESERVE_SPECIES = 363


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def add_genus(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    derived = out["species"].astype(str).str.strip().str.split().str[0]
    if "genus" in out.columns:
        existing = out["genus"].astype(str).str.strip()
        bad = existing.ne(derived)
        if bad.any():
            raise RuntimeError(f"stored/derived genus disagreement: {out.loc[bad, ['species','genus']].head().to_dict('records')}")
    out["genus"] = derived
    return out


def rank_residual(outcome: np.ndarray, span: np.ndarray, n_classifiable: np.ndarray) -> np.ndarray:
    d = np.asarray(outcome, float)
    s = np.asarray(span, float)
    n = np.asarray(n_classifiable, float)
    if not (np.isfinite(d).all() and np.isfinite(s).all() and np.isfinite(n).all()):
        raise RuntimeError("nonfinite value in opportunity residualization")
    rd = stats.rankdata(d).astype(float)
    rs = stats.rankdata(s).astype(float)
    rn = stats.rankdata(n).astype(float)
    X = np.column_stack([np.ones(len(rd)), rs, rn])
    return rd - X @ np.linalg.lstsq(X, rd, rcond=None)[0]


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
    if ng:
        for _name, idx in groups.items():
            n_pairs = len(idx) * (len(idx) - 1) // 2
            w = 1.0 / (ng * n_pairs)
            for a in range(len(idx) - 1):
                for b in range(a + 1, len(idx)):
                    pair_i.append(remap[int(idx[a])])
                    pair_j.append(remap[int(idx[b])])
                    weights.append(w)
    meta = {
        "repeated_groups": int(ng),
        "species_in_repeated_groups": int(len(repeated_species)),
        "group_sizes": {name: int(len(idx)) for name, idx in groups.items()},
        "gate_at_least_10_repeated_groups": bool(ng >= 10),
        "gate_at_least_30_species_in_repeated_groups": bool(len(repeated_species) >= 30),
        "gate_pass": bool(ng >= 10 and len(repeated_species) >= 30),
    }
    return (
        np.asarray(repeated_species, dtype=int),
        np.asarray(pair_i, dtype=int),
        np.asarray(pair_j, dtype=int),
        np.asarray(weights, dtype=float),
        meta,
    )


def clustering_test(values: np.ndarray, labels: pd.Series, *, seed: int) -> dict[str, Any]:
    idx, pair_i, pair_j, weights, meta = repeated_group_pairs(labels)
    result: dict[str, Any] = dict(meta)
    if not result["gate_pass"]:
        result.update({
            "tested": False,
            "W_observed": None,
            "null_mean_W": None,
            "clustering_gain": None,
            "p_lower": 1.0,
        })
        return result
    z = np.asarray(values, float)[idx]
    if not np.isfinite(z).all():
        raise RuntimeError("nonfinite outcome inside repeated genera")
    if len(weights) != len(pair_i) or not np.isclose(weights.sum(), 1.0):
        raise RuntimeError(f"pair weights invalid: n={len(weights)} sum={weights.sum()}")
    observed = float(np.sum(np.abs(z[pair_i] - z[pair_j]) * weights))
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for start in range(0, N_PERM, 500):
        end = min(N_PERM, start + 500)
        perms = np.stack([rng.permutation(len(z)) for _ in range(end - start)], axis=0)
        vals = z[perms]
        null[start:end] = np.sum(np.abs(vals[:, pair_i] - vals[:, pair_j]) * weights[None, :], axis=1)
    p = float((1 + np.sum(null <= observed + 1e-15)) / (N_PERM + 1))
    null_mean = float(null.mean())
    gain = float(1.0 - observed / null_mean) if null_mean > 0 else math.nan
    result.update({
        "tested": True,
        "pair_count_weighted": int(len(pair_i)),
        "W_observed": observed,
        "null_mean_W": null_mean,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "clustering_gain": gain,
        "p_lower": p,
        "supported_at_0_05": bool(gain > 0 and p < 0.05),
    })
    return result


def spearman_perm_two_sided(x: np.ndarray, y: np.ndarray, *, seed: int) -> dict[str, Any]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) != len(y) or len(x) == 0:
        raise RuntimeError("invalid vectors for cross-tranche Spearman")
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    observed = float(rx @ ry / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = float(rx @ ry[rng.permutation(len(ry))] / denom)
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(x)),
        "rho": observed,
        "p_two_sided": p,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def cross_tranche_genus_means(discovery: pd.DataFrame, reserve: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    dc = discovery.groupby("genus").size()
    rc = reserve.groupby("genus").size()
    shared = sorted(set(dc[dc >= 2].index) & set(rc[rc >= 2].index))
    rows = []
    for genus in shared:
        dg = discovery[discovery["genus"].eq(genus)]
        rg = reserve[reserve["genus"].eq(genus)]
        rows.append({
            "genus": genus,
            "discovery_n": int(len(dg)),
            "reserve_n": int(len(rg)),
            "discovery_mean_D": float(dg["D"].mean()),
            "reserve_mean_D": float(rg["D"].mean()),
        })
    table = pd.DataFrame(rows)
    if len(table) < 10:
        return table, {
            "tested": False,
            "shared_repeated_genera": int(len(table)),
            "gate_at_least_10": False,
            "rho": None,
            "p_two_sided": 1.0,
        }
    test = spearman_perm_two_sided(
        table["discovery_mean_D"].to_numpy(float),
        table["reserve_mean_D"].to_numpy(float),
        seed=STEP7_SEED + 40,
    )
    return table, {"tested": True, "shared_repeated_genera": int(len(table)), "gate_at_least_10": True, **test}


def main() -> None:
    discovery = add_genus(pd.read_csv(DISC_PATH)).sort_values("species").reset_index(drop=True)
    reserve = add_genus(pd.read_csv(RES_PATH)).sort_values("species").reset_index(drop=True)
    step4 = json.loads(DISC_META_PATH.read_text(encoding="utf-8"))

    if len(discovery) != EXPECTED_DISCOVERY_SPECIES or discovery["species"].nunique() != EXPECTED_DISCOVERY_SPECIES:
        raise RuntimeError("discovery species fingerprint mismatch")
    if len(reserve) != EXPECTED_RESERVE_SPECIES or reserve["species"].nunique() != EXPECTED_RESERVE_SPECIES:
        raise RuntimeError("reserve species fingerprint mismatch")
    overlap = sorted(set(discovery["species"].astype(str)) & set(reserve["species"].astype(str)))
    if overlap:
        raise RuntimeError(f"expected species-disjoint tranches, found overlap: {overlap[:10]}")

    # Exact Step-4 reproduction with the original frozen seed.
    disc_raw_step4_seed = clustering_test(discovery["D"].to_numpy(float), discovery["genus"], seed=STEP4_SEED)
    expected = step4["primary"]["genus_taxonomic_clustering"]
    for key in ["repeated_groups", "species_in_repeated_groups"]:
        if int(disc_raw_step4_seed[key]) != int(expected[key]):
            raise RuntimeError(f"Step4 genus {key} reproduction mismatch")
    for key in ["W_observed", "clustering_gain", "p_lower"]:
        if abs(float(disc_raw_step4_seed[key]) - float(expected[key])) > 1e-12:
            raise RuntimeError(f"Step4 genus {key} reproduction mismatch: {disc_raw_step4_seed[key]} vs {expected[key]}")

    disc_raw_step7_seed = clustering_test(discovery["D"].to_numpy(float), discovery["genus"], seed=STEP7_SEED + 1)
    res_raw = clustering_test(reserve["D"].to_numpy(float), reserve["genus"], seed=STEP7_SEED + 2)

    disc_resid = rank_residual(
        discovery["D"].to_numpy(float),
        discovery["log1p_span_primary"].to_numpy(float),
        discovery["n_classifiable"].to_numpy(float),
    )
    res_resid = rank_residual(
        reserve["D"].to_numpy(float),
        reserve["reserve_log1p_span"].to_numpy(float),
        reserve["n_classifiable"].to_numpy(float),
    )
    disc_unb_resid = rank_residual(
        discovery["D"].to_numpy(float) * discovery["n_classifiable"].to_numpy(float) / (discovery["n_classifiable"].to_numpy(float) - 1.0),
        discovery["log1p_span_primary"].to_numpy(float),
        discovery["n_classifiable"].to_numpy(float),
    )
    res_unb_resid = rank_residual(
        reserve["D_unbiased"].to_numpy(float),
        reserve["reserve_log1p_span"].to_numpy(float),
        reserve["n_classifiable"].to_numpy(float),
    )

    disc_adjusted = clustering_test(disc_resid, discovery["genus"], seed=STEP7_SEED + 10)
    res_adjusted = clustering_test(res_resid, reserve["genus"], seed=STEP7_SEED + 11)
    disc_adjusted_unb = clustering_test(disc_unb_resid, discovery["genus"], seed=STEP7_SEED + 12)
    res_adjusted_unb = clustering_test(res_unb_resid, reserve["genus"], seed=STEP7_SEED + 13)

    # Measurement/opportunity descriptors, not primary biological tests.
    disc_n_cluster = clustering_test(discovery["n_classifiable"].to_numpy(float), discovery["genus"], seed=STEP7_SEED + 20)
    res_n_cluster = clustering_test(reserve["n_classifiable"].to_numpy(float), reserve["genus"], seed=STEP7_SEED + 21)
    disc_span_cluster = clustering_test(discovery["log1p_span_primary"].to_numpy(float), discovery["genus"], seed=STEP7_SEED + 22)
    res_span_cluster = clustering_test(reserve["reserve_log1p_span"].to_numpy(float), reserve["genus"], seed=STEP7_SEED + 23)

    shared_table, shared_test = cross_tranche_genus_means(discovery, reserve)

    raw_replication = bool(res_raw.get("supported_at_0_05", False))
    adjusted_replication = bool(raw_replication and res_adjusted.get("supported_at_0_05", False))
    discovery_adjusted_survives = bool(disc_adjusted.get("supported_at_0_05", False))

    result = {
        "analysis": "polymorphism_genus_replication_step7",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_GENUS_REPLICATION_STEP7_PROTOCOL_20260910.md",
        "inferential_role": "post_step4_species_disjoint_replication_and_opportunity_robustness",
        "new_image_acquisition": False,
        "species_disjoint_verified": True,
        "n_permutations_each_test": N_PERM,
        "discovery": {
            "raw_step4_seed_exact_reproduction": disc_raw_step4_seed,
            "raw_step7_seed": disc_raw_step7_seed,
            "adjusted_span_and_n": disc_adjusted,
            "adjusted_span_and_n_D_unbiased": disc_adjusted_unb,
        },
        "reserve": {
            "raw": res_raw,
            "adjusted_span_and_n": res_adjusted,
            "adjusted_span_and_n_D_unbiased": res_adjusted_unb,
        },
        "cross_tranche_genus_mean_concordance": shared_test,
        "descriptive_measurement_controls": {
            "discovery_rho_D_n_classifiable": float(stats.spearmanr(discovery["D"], discovery["n_classifiable"]).statistic),
            "reserve_rho_D_n_classifiable": float(stats.spearmanr(reserve["D"], reserve["n_classifiable"]).statistic),
            "discovery_rho_D_sampled_span": float(stats.spearmanr(discovery["D"], discovery["log1p_span_primary"]).statistic),
            "reserve_rho_D_sampled_span": float(stats.spearmanr(reserve["D"], reserve["reserve_log1p_span"]).statistic),
            "discovery_n_classifiable_genus_clustering": disc_n_cluster,
            "reserve_n_classifiable_genus_clustering": res_n_cluster,
            "discovery_span_genus_clustering": disc_span_cluster,
            "reserve_span_genus_clustering": res_span_cluster,
        },
        "decision": {
            "raw_genus_clustering_replicated_in_species_disjoint_reserve": raw_replication,
            "opportunity_robust_genus_clustering_replicated_in_species_disjoint_reserve": adjusted_replication,
            "discovery_genus_clustering_survives_span_and_n_adjustment": discovery_adjusted_survives,
            "cross_tranche_mean_concordance_supportive": bool(shared_test.get("tested", False) and shared_test.get("rho", 0) is not None and float(shared_test["rho"]) > 0 and float(shared_test["p_two_sided"]) < 0.05),
            "claim_boundary": "Genus-level taxonomic clustering is not formal phylogenetic signal and does not establish genetic determination, adaptation, causal genus effects, or angiosperm-wide prevalence.",
        },
    }

    discovery[["species", "genus", "D", "n_classifiable", "log1p_span_primary"]].assign(opportunity_adjusted_D_rank_residual=disc_resid).to_csv(
        OUT / "discovery_genus_input.csv", index=False
    )
    reserve[["species", "genus", "D", "n_classifiable", "reserve_log1p_span"]].assign(opportunity_adjusted_D_rank_residual=res_resid).to_csv(
        OUT / "reserve_genus_input.csv", index=False
    )
    shared_table.to_csv(OUT / "shared_repeated_genus_means.csv", index=False)
    write_json(OUT / "result.json", result)

    lines = [
        "# Polymorphism genus clustering — Step 7 result",
        "",
        f"**Raw genus clustering replicated in the species-disjoint reserve: `{raw_replication}`.**",
        f"**Opportunity-robust genus clustering replicated: `{adjusted_replication}`.**",
        f"**Discovery genus clustering survives span+n adjustment: `{discovery_adjusted_survives}`.**",
        "",
        "## Discovery",
        "",
        f"- raw Step-4 reproduction: repeated genera **{disc_raw_step4_seed['repeated_groups']}**, species **{disc_raw_step4_seed['species_in_repeated_groups']}**, gain **{disc_raw_step4_seed['clustering_gain']:.6f}**, p_lower **{disc_raw_step4_seed['p_lower']:.6g}**",
        f"- adjusted for sampled span + n: gain **{disc_adjusted['clustering_gain']:.6f}**, p_lower **{disc_adjusted['p_lower']:.6g}**",
        f"- D_unbiased adjusted: gain **{disc_adjusted_unb['clustering_gain']:.6f}**, p_lower **{disc_adjusted_unb['p_lower']:.6g}**",
        "",
        "## Species-disjoint reserve",
        "",
        f"- repeated genera **{res_raw['repeated_groups']}**, species **{res_raw['species_in_repeated_groups']}**",
        f"- raw D: gain **{res_raw['clustering_gain']:.6f}**, p_lower **{res_raw['p_lower']:.6g}**",
        f"- adjusted for sampled span + n: gain **{res_adjusted['clustering_gain']:.6f}**, p_lower **{res_adjusted['p_lower']:.6g}**",
        f"- D_unbiased adjusted: gain **{res_adjusted_unb['clustering_gain']:.6f}**, p_lower **{res_adjusted_unb['p_lower']:.6g}**",
        "",
        "## Cross-tranche repeated-genus means",
        "",
        f"- shared genera with >=2 species in both tranches: **{shared_test['shared_repeated_genera']}**",
        f"- test run: **{shared_test['tested']}**",
    ]
    if shared_test.get("tested", False):
        lines.extend([
            f"- Spearman rho = **{shared_test['rho']:.6f}**, p_two_sided = **{shared_test['p_two_sided']:.6g}**",
        ])
    lines.extend([
        "",
        "## Measurement/opportunity diagnostics",
        "",
        f"- discovery rho(D, n_classifiable) = **{result['descriptive_measurement_controls']['discovery_rho_D_n_classifiable']:.6f}**",
        f"- reserve rho(D, n_classifiable) = **{result['descriptive_measurement_controls']['reserve_rho_D_n_classifiable']:.6f}**",
        f"- discovery rho(D, sampled span) = **{result['descriptive_measurement_controls']['discovery_rho_D_sampled_span']:.6f}**",
        f"- reserve rho(D, sampled span) = **{result['descriptive_measurement_controls']['reserve_rho_D_sampled_span']:.6f}**",
        "",
        "## Claim boundary",
        "",
        "A replicated result is genus-level taxonomic clustering of photo-derived species polymorphism diversity. It is not formal phylogenetic signal, genetic determination, adaptation, or a causal effect of genus membership.",
    ])
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps(result["decision"], indent=2))
    print((OUT / "RESULT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
