#!/usr/bin/env python3
"""Run the prospectively frozen Step 4 species-attribute association family."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910"
OUT.mkdir(parents=True, exist_ok=True)

D_PATH = ROOT / "results" / "polymorphism_directionality_step1_20260909" / "species_metrics.csv"
PANEL_PATH = ROOT / "results" / "polymorphism_species_attributes_step4_preflight_20260910" / "covariate_panel_preoutcome.csv"
PREFLIGHT_JSON = ROOT / "results" / "polymorphism_species_attributes_step4_preflight_20260910" / "result.json"
LIFE_JSON = ROOT / "results" / "polymorphism_species_attributes_step4_lifeform_source_20260910" / "result.json"
LIFE_PANEL = ROOT / "results" / "polymorphism_species_attributes_step4_lifeform_source_20260910" / "life_form_preoutcome.csv"

SEED = 20260910
N_PERM = 20_000
SIX_SLOTS = [
    "sampled_geographic_span",
    "absolute_latitude_centroid",
    "family_taxonomic_clustering",
    "genus_taxonomic_clustering",
    "pollination",
    "life_form",
]


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def rank_center(x: np.ndarray) -> np.ndarray:
    r = stats.rankdata(np.asarray(x, float)).astype(float)
    return r - r.mean()


def spearman_fast(x: np.ndarray, y: np.ndarray) -> float:
    rx = rank_center(x)
    ry = rank_center(y)
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    return float(rx @ ry / denom)


def permutation_spearman_two_sided(
    outcome: np.ndarray,
    predictor: np.ndarray,
    *,
    seed: int,
) -> dict[str, Any]:
    keep = np.isfinite(outcome) & np.isfinite(predictor)
    d = np.asarray(outcome[keep], float)
    x = np.asarray(predictor[keep], float)
    rd = rank_center(d)
    rx = rank_center(x)
    denom = float(np.linalg.norm(rd) * np.linalg.norm(rx))
    observed = float(rd @ rx / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = float(rd[rng.permutation(len(rd))] @ rx / denom)
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(d)),
        "rho": observed,
        "p_two_sided": p,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def partial_spearman(outcome: np.ndarray, predictor: np.ndarray, control: np.ndarray) -> float:
    keep = np.isfinite(outcome) & np.isfinite(predictor) & np.isfinite(control)
    rd = stats.rankdata(outcome[keep]).astype(float)
    rx = stats.rankdata(predictor[keep]).astype(float)
    rc = stats.rankdata(control[keep]).astype(float)
    X = np.column_stack([np.ones(len(rc)), rc])
    ed = rd - X @ np.linalg.lstsq(X, rd, rcond=None)[0]
    ex = rx - X @ np.linalg.lstsq(X, rx, rcond=None)[0]
    return float(np.corrcoef(ed, ex)[0, 1])


def repeated_group_pairs(labels: pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
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
        {**meta, "weights": np.asarray(weights, dtype=float)},
    )


def taxonomic_clustering_test(
    outcome: np.ndarray,
    labels: pd.Series,
    *,
    seed: int,
) -> dict[str, Any]:
    species_idx, pair_i, pair_j, meta = repeated_group_pairs(labels)
    weights = meta.pop("weights")
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

    d = np.asarray(outcome[species_idx], float)
    if not np.isfinite(d).all():
        raise RuntimeError("nonfinite outcome inside repeated taxonomic groups")
    if len(weights) != len(pair_i) or not np.isclose(weights.sum(), 1.0):
        raise RuntimeError(f"taxonomic pair weights invalid: n={len(weights)}, sum={weights.sum()}")
    observed = float(np.sum(np.abs(d[pair_i] - d[pair_j]) * weights))

    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for start in range(0, N_PERM, 500):
        end = min(N_PERM, start + 500)
        perms = np.stack([rng.permutation(len(d)) for _ in range(end - start)], axis=0)
        vals = d[perms]
        null[start:end] = np.sum(
            np.abs(vals[:, pair_i] - vals[:, pair_j]) * weights[None, :],
            axis=1,
        )
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
    })
    return result


def holm_adjust(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    sorted_p = p[order]
    raw_adj = np.array([(m - i) * sorted_p[i] for i in range(m)], float)
    sorted_adj = np.minimum(1.0, np.maximum.accumulate(raw_adj))
    out = np.empty(m, float)
    out[order] = sorted_adj
    return [float(x) for x in out]


def main() -> None:
    # ---- Stage A: validate all pre-outcome covariates and gates BEFORE opening D. ----
    panel = pd.read_csv(PANEL_PATH)
    preflight = json.loads(PREFLIGHT_JSON.read_text(encoding="utf-8"))
    life_meta = json.loads(LIFE_JSON.read_text(encoding="utf-8"))
    life_panel = pd.read_csv(LIFE_PANEL)

    if len(panel) != 369 or panel["species"].nunique() != 369:
        raise RuntimeError("pre-outcome covariate panel fingerprint mismatch")
    if preflight.get("D_association_computed") is not False:
        raise RuntimeError("preflight outcome firewall receipt invalid")
    if life_meta.get("D_association_computed") is not False:
        raise RuntimeError("life-form outcome firewall receipt invalid")
    if bool(preflight["pollination"]["pollination_gate_pass"]):
        raise RuntimeError("protocol expected pollination gate to be closed")
    if bool(life_meta["life_form_gate_pass"]):
        raise RuntimeError("protocol expected life-form gate to be closed")
    if int(panel["log1p_span_primary"].notna().sum()) != 369:
        raise RuntimeError("primary span coverage changed")
    if int(panel["abs_mean_latitude"].notna().sum()) != 369:
        raise RuntimeError("latitude coverage changed")
    if int(panel["family"].notna().sum()) != 345:
        raise RuntimeError("family coverage changed")
    if int(life_panel["life_form"].notna().sum()) != 2:
        raise RuntimeError("life-form coverage changed")

    panel = panel.sort_values("species").reset_index(drop=True)
    family_gate_preview = repeated_group_pairs(panel["family"])[3]
    family_gate_preview.pop("weights")
    genus_gate_preview = repeated_group_pairs(panel["genus"])[3]
    genus_gate_preview.pop("weights")

    # ---- Stage B: now open the frozen D outcome. Covariates are immutable beyond here. ----
    dframe = pd.read_csv(D_PATH, usecols=["species", "n", "D", "second_fraction"])
    if len(dframe) != 369 or dframe["species"].nunique() != 369:
        raise RuntimeError("D frame fingerprint mismatch")
    if abs(float(dframe["D"].max()) - 0.707645) > 1e-5:
        raise RuntimeError(f"D max fingerprint mismatch: {dframe['D'].max()}")
    if abs(float((dframe["second_fraction"] >= 0.10).mean()) - 0.4661) > 5e-4:
        raise RuntimeError("second-morph fingerprint mismatch")

    merged = panel.merge(dframe, on="species", how="inner", validate="one_to_one", suffixes=("", "_d"))
    if len(merged) != 369:
        raise RuntimeError("D/covariate join lost species")
    if not np.array_equal(merged["n_classifiable"].astype(int).to_numpy(), merged["n"].astype(int).to_numpy()):
        raise RuntimeError("n_classifiable mismatch between pre-outcome panel and D frame")

    D = merged["D"].to_numpy(float)
    n = merged["n_classifiable"].to_numpy(float)
    Dunb = D * n / (n - 1.0)
    span = merged["log1p_span_primary"].to_numpy(float)
    span_alt = merged["log1p_span_sensitivity"].to_numpy(float)
    alat = merged["abs_mean_latitude"].to_numpy(float)
    signed_lat = merged["mean_latitude"].to_numpy(float)

    span_primary = permutation_spearman_two_sided(D, span, seed=SEED + 1)
    lat_primary = permutation_spearman_two_sided(D, alat, seed=SEED + 2)
    family_primary = taxonomic_clustering_test(D, merged["family"], seed=SEED + 3)
    genus_primary = taxonomic_clustering_test(D, merged["genus"], seed=SEED + 4)

    # Predeclared sensitivities; no sensitivity changes the six-slot primary family.
    sensitivity = {
        "span_D_unbiased": permutation_spearman_two_sided(Dunb, span, seed=SEED + 11),
        "span_partial_rank_controlling_n_classifiable": partial_spearman(D, span, n),
        "span_alternate_maximum_span_km": permutation_spearman_two_sided(D, span_alt, seed=SEED + 12),
        "latitude_D_unbiased": permutation_spearman_two_sided(Dunb, alat, seed=SEED + 13),
        "latitude_partial_rank_controlling_n_classifiable": partial_spearman(D, alat, n),
        "signed_latitude_descriptive": permutation_spearman_two_sided(D, signed_lat, seed=SEED + 14),
        "family_D_unbiased": taxonomic_clustering_test(Dunb, merged["family"], seed=SEED + 15),
        "genus_D_unbiased": taxonomic_clustering_test(Dunb, merged["genus"], seed=SEED + 16),
    }

    raw_p = {
        "sampled_geographic_span": float(span_primary["p_two_sided"]),
        "absolute_latitude_centroid": float(lat_primary["p_two_sided"]),
        "family_taxonomic_clustering": float(family_primary["p_lower"]),
        "genus_taxonomic_clustering": float(genus_primary["p_lower"]),
        "pollination": 1.0,
        "life_form": 1.0,
    }
    adj = holm_adjust([raw_p[x] for x in SIX_SLOTS])
    holm = {slot: adj[i] for i, slot in enumerate(SIX_SLOTS)}

    primary = {
        "sampled_geographic_span": span_primary,
        "absolute_latitude_centroid": lat_primary,
        "family_taxonomic_clustering": family_primary,
        "genus_taxonomic_clustering": genus_primary,
        "pollination": {
            "tested": False,
            "raw_p": 1.0,
            "reason": "prospective source-backed category gate failed before D was opened",
            "source_backed_species": int(preflight["pollination"]["species_with_admissible_exact_source_backed_pollination"]),
            "mapped_counts": preflight["pollination"]["mapped_counts"],
        },
        "life_form": {
            "tested": False,
            "raw_p": 1.0,
            "reason": "prospective direct GIFT Life_form_1 coverage gate failed before D was opened",
            "direct_species": int(life_meta["fcp_species_with_life_form"]),
            "counts": life_meta["life_form_counts"],
        },
    }

    decisions: dict[str, Any] = {}
    for slot in SIX_SLOTS:
        decisions[slot] = {
            "raw_p": raw_p[slot],
            "holm_p": holm[slot],
            "supported_holm_0_05": bool(holm[slot] < 0.05),
        }

    result = {
        "analysis": "polymorphism_species_attributes_step4_association",
        "protocol": "docs/POLYMORPHISM_SPECIES_ATTRIBUTES_STEP4_ASSOCIATION_PROTOCOL_20260910.md",
        "date_jst": "2026-09-10",
        "new_image_acquisition": False,
        "n_species": 369,
        "permutations_per_active_test": N_PERM,
        "seed": SEED,
        "pre_outcome_gate_preview": {
            "family": family_gate_preview,
            "genus": genus_gate_preview,
            "pollination_gate_pass": False,
            "life_form_gate_pass": False,
        },
        "primary": primary,
        "sensitivities": sensitivity,
        "six_slot_raw_p": raw_p,
        "six_slot_holm_p": holm,
        "decisions": decisions,
        "claim_boundary": (
            "Supported tests are observational species-attribute associations or taxonomic clustering. "
            "They do not establish adaptation, causality, genetic mechanism, true range size, or a shared global colour boundary."
        ),
    }
    write_json(OUT / "result.json", result)

    export_cols = [
        "species", "genus", "family", "n_classifiable", "D", "second_fraction",
        "maximum_span_km_after_observer_cap", "log1p_span_primary",
        "maximum_span_km", "log1p_span_sensitivity", "mean_latitude", "abs_mean_latitude",
        "pollination_guild_raw", "pollination_guild_mapped",
    ]
    merged[export_cols].to_csv(OUT / "species_analysis_table.csv", index=False)

    supported = [slot for slot in SIX_SLOTS if decisions[slot]["supported_holm_0_05"]]
    md = [
        "# Polymorphism species attributes — Step 4 association result",
        "",
        f"**Holm-supported families: `{supported}`.**",
        "",
        "## Primary six-family result",
        "",
        f"- sampled geographic span: rho = **{span_primary['rho']:.6f}**, raw p = **{raw_p['sampled_geographic_span']:.6g}**, Holm p = **{holm['sampled_geographic_span']:.6g}**",
        f"- absolute latitude centroid: rho = **{lat_primary['rho']:.6f}**, raw p = **{raw_p['absolute_latitude_centroid']:.6g}**, Holm p = **{holm['absolute_latitude_centroid']:.6g}**",
        f"- family clustering: groups = **{family_primary['repeated_groups']}**, species = **{family_primary['species_in_repeated_groups']}**, gain = **{family_primary['clustering_gain'] if family_primary['clustering_gain'] is not None else float('nan'):.6f}**, raw p = **{raw_p['family_taxonomic_clustering']:.6g}**, Holm p = **{holm['family_taxonomic_clustering']:.6g}**",
        f"- genus clustering: groups = **{genus_primary['repeated_groups']}**, species = **{genus_primary['species_in_repeated_groups']}**, gain = **{genus_primary['clustering_gain'] if genus_primary['clustering_gain'] is not None else float('nan'):.6f}**, raw p = **{raw_p['genus_taxonomic_clustering']:.6g}**, Holm p = **{holm['genus_taxonomic_clustering']:.6g}**",
        f"- pollination: **not tested; prospective gate failed**, Holm p = **{holm['pollination']:.6g}**",
        f"- life form: **not tested; prospective gate failed**, Holm p = **{holm['life_form']:.6g}**",
        "",
        "## Predeclared sensitivities",
        "",
        f"- span with D_unbiased: rho = **{sensitivity['span_D_unbiased']['rho']:.6f}**, p = **{sensitivity['span_D_unbiased']['p_two_sided']:.6g}**",
        f"- span partial rank controlling n_classifiable: **{sensitivity['span_partial_rank_controlling_n_classifiable']:.6f}**",
        f"- alternate span proxy: rho = **{sensitivity['span_alternate_maximum_span_km']['rho']:.6f}**, p = **{sensitivity['span_alternate_maximum_span_km']['p_two_sided']:.6g}**",
        f"- latitude with D_unbiased: rho = **{sensitivity['latitude_D_unbiased']['rho']:.6f}**, p = **{sensitivity['latitude_D_unbiased']['p_two_sided']:.6g}**",
        f"- latitude partial rank controlling n_classifiable: **{sensitivity['latitude_partial_rank_controlling_n_classifiable']:.6f}**",
        f"- signed latitude descriptive rho = **{sensitivity['signed_latitude_descriptive']['rho']:.6f}**, p = **{sensitivity['signed_latitude_descriptive']['p_two_sided']:.6g}**",
        f"- family D_unbiased gain = **{sensitivity['family_D_unbiased']['clustering_gain'] if sensitivity['family_D_unbiased']['clustering_gain'] is not None else float('nan'):.6f}**, p = **{sensitivity['family_D_unbiased']['p_lower']:.6g}**",
        f"- genus D_unbiased gain = **{sensitivity['genus_D_unbiased']['clustering_gain'] if sensitivity['genus_D_unbiased']['clustering_gain'] is not None else float('nan'):.6f}**, p = **{sensitivity['genus_D_unbiased']['p_lower']:.6g}**",
        "",
        "All claims remain observational. `maximum_span_km_after_observer_cap` is a sampled geographic-opportunity proxy, not true species range size; family/genus results are taxonomic clustering, not phylogenetic signal.",
        "",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()
