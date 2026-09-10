#!/usr/bin/env python3
"""Decompose sampled-span versus post-classification-yield sensitivity after Step 6."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_opportunity_decomposition_step6b_20260910"
OUT.mkdir(parents=True, exist_ok=True)

STEP6_SCRIPT = ROOT / "scripts" / "analysis" / "run_polymorphism_span_adjusted_robustness_step6_20260910.py"
STEP6_JSON = ROOT / "results" / "polymorphism_span_adjusted_robustness_step6_20260910" / "result.json"
DISCOVERY = ROOT / "results" / "polymorphism_span_adjusted_robustness_step6_20260910" / "discovery_joined.csv"
RESERVE = ROOT / "results" / "polymorphism_span_adjusted_robustness_step6_20260910" / "reserve_joined.csv"
PANEL = ROOT / "results" / "polymorphism_species_attributes_step4_preflight_20260910" / "covariate_panel_preoutcome.csv"
RESERVE_GEOM = ROOT / "data" / "frozen" / "rgfca_reserve_geometry_audit_v1.csv"
PROTOCOL = "docs/POLYMORPHISM_OPPORTUNITY_DECOMPOSITION_STEP6B_PROTOCOL_20260910.md"
SEED = 20260910 + 650


def load_step6_module():
    spec = importlib.util.spec_from_file_location("step6_frozen", STEP6_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Step 6 module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def prepare(frame: pd.DataFrame, *, spatial_col: str, span_col: str) -> dict[str, Any]:
    D = frame["D"].to_numpy(float)
    n = frame["n_classifiable"].to_numpy(float)
    Dunb = D * n / (n - 1.0)
    return {
        "D": D,
        "D_unbiased": Dunb,
        "n": n,
        "span": frame[span_col].to_numpy(float),
        "spatial": frame[spatial_col].to_numpy(float),
        "indicator": frame["primary_second_ge_0_10"].astype(int).to_numpy(),
        "genus": frame["genus"],
    }


def main() -> None:
    m = load_step6_module()
    step6 = json.loads(STEP6_JSON.read_text(encoding="utf-8"))
    disc = pd.read_csv(DISCOVERY)
    res = pd.read_csv(RESERVE)
    panel = pd.read_csv(PANEL)
    geom = pd.read_csv(RESERVE_GEOM)

    if len(disc) != 369 or len(res) != 363:
        raise RuntimeError("Step 6 joined-frame fingerprint mismatch")
    if int(panel["latitude_n"].eq(100).sum()) != 369:
        raise RuntimeError("discovery raw-photo opportunity is not fixed at 100")
    reserve_ids = set(res["inat_taxon_id"].astype(int))
    geom_sub = geom.loc[geom["inat_taxon_id"].astype(int).isin(reserve_ids)].copy()
    if len(geom_sub) != 363 or int(geom_sub["raw_photos"].eq(100).sum()) != 363:
        raise RuntimeError("reserve raw-photo opportunity is not fixed at 100")

    d = prepare(disc, spatial_col="spatial_observed_rho", span_col="log1p_span_primary")
    r = prepare(res, spatial_col="rho_primary", span_col="log1p_span_reserve")

    # A. Threat-specific span-only family.
    d_span = m.partial_rank_freedman_lane(d["spatial"], d["D"], [d["span"]], seed=SEED + 1)
    r_span = m.partial_rank_freedman_lane(r["spatial"], r["D"], [r["span"]], seed=SEED + 2)
    span_holm = m.holm_adjust([d_span["p_two_sided"], r_span["p_two_sided"]])
    d_span["holm_p_two_tranches"] = span_holm[0]
    r_span["holm_p_two_tranches"] = span_holm[1]
    span_only_gate = bool(
        d_span["partial_rho"] > 0
        and r_span["partial_rho"] > 0
        and d_span["holm_p_two_tranches"] < 0.05
        and r_span["holm_p_two_tranches"] < 0.05
    )

    # B. Finite-sample unbiased Simpson sensitivity, span only.
    d_span_unb = m.partial_rank_freedman_lane(d["spatial"], d["D_unbiased"], [d["span"]], seed=SEED + 3)
    r_span_unb = m.partial_rank_freedman_lane(r["spatial"], r["D_unbiased"], [r["span"]], seed=SEED + 4)
    span_unb_holm = m.holm_adjust([d_span_unb["p_two_sided"], r_span_unb["p_two_sided"]])
    d_span_unb["holm_p_two_tranches"] = span_unb_holm[0]
    r_span_unb["holm_p_two_tranches"] = span_unb_holm[1]
    span_unbiased_gate = bool(
        d_span_unb["partial_rho"] > 0
        and r_span_unb["partial_rho"] > 0
        and d_span_unb["holm_p_two_tranches"] < 0.05
        and r_span_unb["holm_p_two_tranches"] < 0.05
    )

    # C. >=10% threshold with span only.
    d_thr_span = m.ancova_indicator_freedman_lane(d["spatial"], d["indicator"], [d["span"]], seed=SEED + 5)
    r_thr_span = m.ancova_indicator_freedman_lane(r["spatial"], r["indicator"], [r["span"]], seed=SEED + 6)
    thr_holm = m.holm_adjust([d_thr_span["p_two_sided"], r_thr_span["p_two_sided"]])
    d_thr_span["holm_p_two_tranches"] = thr_holm[0]
    r_thr_span["holm_p_two_tranches"] = thr_holm[1]
    threshold_span_gate = bool(
        d_thr_span["adjusted_raw_spatial_rho_difference"] > 0
        and r_thr_span["adjusted_raw_spatial_rho_difference"] > 0
        and d_thr_span["holm_p_two_tranches"] < 0.05
        and r_thr_span["holm_p_two_tranches"] < 0.05
    )

    # D. n-only and joint diagnostics; these are not new claim gates.
    d_n = m.partial_rank_freedman_lane(d["spatial"], d["D"], [d["n"]], seed=SEED + 7)
    r_n = m.partial_rank_freedman_lane(r["spatial"], r["D"], [r["n"]], seed=SEED + 8)
    d_n_unb = m.partial_rank_freedman_lane(d["spatial"], d["D_unbiased"], [d["n"]], seed=SEED + 9)
    r_n_unb = m.partial_rank_freedman_lane(r["spatial"], r["D_unbiased"], [r["n"]], seed=SEED + 10)
    d_both = m.partial_rank_freedman_lane(d["spatial"], d["D"], [d["span"], d["n"]], seed=SEED + 11)
    r_both = m.partial_rank_freedman_lane(r["spatial"], r["D"], [r["span"], r["n"]], seed=SEED + 12)

    # Fingerprint the joint Step 6 effect sizes, independent of permutation-seed p differences.
    old_d = float(step6["claim2_continuous_span_adjusted"]["discovery"]["partial_rho"])
    old_r = float(step6["claim2_continuous_span_adjusted"]["reserve"]["partial_rho"])
    if abs(d_both["partial_rho"] - old_d) > 1e-12 or abs(r_both["partial_rho"] - old_r) > 1e-12:
        raise RuntimeError("joint Step 6 partial-rho fingerprint mismatch")

    d_D_n = m.permutation_spearman_two_sided(d["D"], d["n"], seed=SEED + 13)
    r_D_n = m.permutation_spearman_two_sided(r["D"], r["n"], seed=SEED + 14)
    d_Dunb_n = m.permutation_spearman_two_sided(d["D_unbiased"], d["n"], seed=SEED + 15)
    r_Dunb_n = m.permutation_spearman_two_sided(r["D_unbiased"], r["n"], seed=SEED + 16)

    # E. Genus clustering decomposition.
    def genus_suite(x: dict[str, Any], base_seed: int) -> dict[str, Any]:
        D = x["D"]
        Du = x["D_unbiased"]
        span = x["span"]
        n = x["n"]
        genus = x["genus"]
        return {
            "raw_D": m.taxonomic_clustering_test(D, genus, seed=base_seed + 1),
            "raw_D_unbiased": m.taxonomic_clustering_test(Du, genus, seed=base_seed + 2),
            "span_only_D": m.taxonomic_clustering_test(
                m.residualize_raw_on_rank_controls(D, span), genus, seed=base_seed + 3
            ),
            "span_only_D_unbiased": m.taxonomic_clustering_test(
                m.residualize_raw_on_rank_controls(Du, span), genus, seed=base_seed + 4
            ),
            "n_only_D": m.taxonomic_clustering_test(
                m.residualize_raw_on_rank_controls(D, n), genus, seed=base_seed + 5
            ),
            "span_plus_n_D": m.taxonomic_clustering_test(
                m.residualize_raw_on_rank_controls(D, span, n), genus, seed=base_seed + 6
            ),
        }

    genus_disc = genus_suite(d, SEED + 100)
    genus_res = genus_suite(r, SEED + 200)

    # Additional descriptives needed to interpret the decomposition.
    descriptives = {
        "raw_photo_opportunity_discovery_fixed_100": True,
        "raw_photo_opportunity_reserve_fixed_100": True,
        "discovery_classifiable_n_min": int(np.min(d["n"])),
        "discovery_classifiable_n_max": int(np.max(d["n"])),
        "reserve_classifiable_n_min": int(np.min(r["n"])),
        "reserve_classifiable_n_max": int(np.max(r["n"])),
        "discovery_rho_D_n": float(stats.spearmanr(d["D"], d["n"]).statistic),
        "reserve_rho_D_n": float(stats.spearmanr(r["D"], r["n"]).statistic),
        "discovery_rho_Dunbiased_n": float(stats.spearmanr(d["D_unbiased"], d["n"]).statistic),
        "reserve_rho_Dunbiased_n": float(stats.spearmanr(r["D_unbiased"], r["n"]).statistic),
    }

    if span_only_gate:
        span_verdict = "SAMPLED_SPAN_INSUFFICIENT_TO_EXPLAIN_CONTINUOUS_PATTERN"
    elif d_span["partial_rho"] > 0 and r_span["partial_rho"] > 0:
        span_verdict = "SPAN_ADJUSTED_EFFECTS_POSITIVE_BUT_NOT_TWO_TRANCHE_HOLM_SUPPORTED"
    else:
        span_verdict = "SAMPLED_SPAN_CAN_ACCOUNT_FOR_OR_REVERSE_AT_LEAST_ONE_TRANCHE"

    result = {
        "analysis": "polymorphism_opportunity_decomposition_step6b",
        "date_jst": "2026-09-10",
        "protocol": PROTOCOL,
        "post_outcome_diagnostic": True,
        "new_data_acquisition": False,
        "span_only_continuous": {
            "discovery": d_span,
            "reserve": r_span,
            "two_tranche_gate_pass": span_only_gate,
            "verdict": span_verdict,
        },
        "span_only_D_unbiased": {
            "discovery": d_span_unb,
            "reserve": r_span_unb,
            "two_tranche_gate_pass": span_unbiased_gate,
        },
        "span_only_threshold_ge_0_10": {
            "discovery": d_thr_span,
            "reserve": r_thr_span,
            "two_tranche_gate_pass": threshold_span_gate,
        },
        "n_only_sensitivity": {
            "raw_D": {"discovery": d_n, "reserve": r_n},
            "D_unbiased": {"discovery": d_n_unb, "reserve": r_n_unb},
            "D_vs_n": {"discovery": d_D_n, "reserve": r_D_n},
            "D_unbiased_vs_n": {"discovery": d_Dunb_n, "reserve": r_Dunb_n},
        },
        "joint_span_plus_n_fingerprint": {"discovery": d_both, "reserve": r_both},
        "genus_clustering_decomposition": {"discovery": genus_disc, "reserve": genus_res},
        "descriptives": descriptives,
        "interpretation_ceiling": (
            "Span-only survival excludes sampled geographic span as a sufficient explanation of the observed "
            "D-spatial relationship. Sensitivity to n_classifiable remains a post-measurement-yield limitation "
            "and does not by itself identify causal measurement bias."
        ),
    }
    write_json(OUT / "result.json", result)

    def gp(suite: dict[str, Any], key: str) -> str:
        x = suite[key]
        gain = x.get("clustering_gain")
        gain_text = "NA" if gain is None else f"{float(gain):.6f}"
        return f"gain={gain_text}, p={float(x['p_lower']):.6g}"

    lines = [
        "# Polymorphism Step 6b — opportunity decomposition result",
        "",
        f"**Span-specific verdict: `{span_verdict}`.**",
        "",
        "## Span-only D–spatial adjustment",
        "",
        f"- discovery raw D: partial rho = **{d_span['partial_rho']:.6f}**, raw p = **{d_span['p_two_sided']:.6g}**, Holm p = **{d_span['holm_p_two_tranches']:.6g}**",
        f"- reserve raw D: partial rho = **{r_span['partial_rho']:.6f}**, raw p = **{r_span['p_two_sided']:.6g}**, Holm p = **{r_span['holm_p_two_tranches']:.6g}**",
        f"- two-tranche span-only gate: **{span_only_gate}**",
        "",
        "## D_unbiased + span-only sensitivity",
        "",
        f"- discovery: partial rho = **{d_span_unb['partial_rho']:.6f}**, Holm p = **{d_span_unb['holm_p_two_tranches']:.6g}**",
        f"- reserve: partial rho = **{r_span_unb['partial_rho']:.6f}**, Holm p = **{r_span_unb['holm_p_two_tranches']:.6g}**",
        f"- two-tranche D_unbiased span-only gate: **{span_unbiased_gate}**",
        "",
        "## >=10% threshold + span-only",
        "",
        f"- discovery adjusted spatial-rho difference = **{d_thr_span['adjusted_raw_spatial_rho_difference']:.6f}**, Holm p = **{d_thr_span['holm_p_two_tranches']:.6g}**",
        f"- reserve adjusted spatial-rho difference = **{r_thr_span['adjusted_raw_spatial_rho_difference']:.6f}**, Holm p = **{r_thr_span['holm_p_two_tranches']:.6g}**",
        f"- two-tranche threshold gate: **{threshold_span_gate}**",
        "",
        "## n_classifiable-only sensitivity",
        "",
        f"- discovery raw D: partial rho = **{d_n['partial_rho']:.6f}**, p = **{d_n['p_two_sided']:.6g}**",
        f"- reserve raw D: partial rho = **{r_n['partial_rho']:.6f}**, p = **{r_n['p_two_sided']:.6g}**",
        f"- discovery D_unbiased: partial rho = **{d_n_unb['partial_rho']:.6f}**, p = **{d_n_unb['p_two_sided']:.6g}**",
        f"- reserve D_unbiased: partial rho = **{r_n_unb['partial_rho']:.6f}**, p = **{r_n_unb['p_two_sided']:.6g}**",
        f"- rho(D, n): discovery **{descriptives['discovery_rho_D_n']:.6f}**, reserve **{descriptives['reserve_rho_D_n']:.6f}**",
        f"- rho(D_unbiased, n): discovery **{descriptives['discovery_rho_Dunbiased_n']:.6f}**, reserve **{descriptives['reserve_rho_Dunbiased_n']:.6f}**",
        "",
        "## Genus clustering decomposition",
        "",
        f"- discovery raw D: {gp(genus_disc, 'raw_D')} ; span-only: {gp(genus_disc, 'span_only_D')} ; n-only: {gp(genus_disc, 'n_only_D')} ; span+n: {gp(genus_disc, 'span_plus_n_D')}",
        f"- discovery D_unbiased raw: {gp(genus_disc, 'raw_D_unbiased')} ; span-only: {gp(genus_disc, 'span_only_D_unbiased')}",
        f"- reserve raw D: {gp(genus_res, 'raw_D')} ; span-only: {gp(genus_res, 'span_only_D')} ; n-only: {gp(genus_res, 'n_only_D')} ; span+n: {gp(genus_res, 'span_plus_n_D')}",
        f"- reserve D_unbiased raw: {gp(genus_res, 'raw_D_unbiased')} ; span-only: {gp(genus_res, 'span_only_D_unbiased')}",
        "",
        "Raw photo opportunity is fixed at 100 photos/species in both tranches. `n_classifiable` is therefore post-classification yield rather than raw sampling effort. Span-only survival does not erase the separate sensitivity to classifiable yield.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
