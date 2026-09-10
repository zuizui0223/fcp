#!/usr/bin/env python3
"""Build deterministic figure-data tables and a synchronized number ledger for FCP paper v0.1."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "paper" / "polymorphism_v0_1"
FIG = OUT / "figure_data"
FIG.mkdir(parents=True, exist_ok=True)

DISC_ATTR = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
DISC_SP = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909" / "species_membership_and_observed_rho.csv"
RES_SP = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "reserve_species_membership_and_observed_metrics.csv"
STEP1 = ROOT / "results" / "polymorphism_directionality_step1_20260909" / "result.json"
STEP4 = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "result.json"
STEP5 = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909" / "result.json"
STEP5B = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "result.json"
STEP8 = ROOT / "results" / "polymorphism_measurement_missingness_step8_analysis_20260910" / "result.json"
STEP8_PRE = ROOT / "results" / "polymorphism_measurement_missingness_step8_preflight_20260910" / "result.json"
STEP9 = ROOT / "results" / "polymorphism_ambiguity_bounds_step9_20260910" / "result.json"
STEP9_DISC_PANEL = ROOT / "results" / "polymorphism_ambiguity_bounds_step9_20260910" / "discovery_ambiguity_endpoint_panel.csv"
STEP9_RES_PANEL = ROOT / "results" / "polymorphism_ambiguity_bounds_step9_20260910" / "reserve_ambiguity_endpoint_panel.csv"
STEP7 = ROOT / "results" / "polymorphism_genus_replication_step7_20260910" / "result.json"
STEP7B = ROOT / "results" / "polymorphism_genus_adjustment_decomposition_step7b_20260910" / "result.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    s1 = load_json(STEP1)
    s4 = load_json(STEP4)
    s5 = load_json(STEP5)
    s5b = load_json(STEP5B)
    s8 = load_json(STEP8)
    s8pre = load_json(STEP8_PRE)
    s9 = load_json(STEP9)
    s7 = load_json(STEP7)
    s7b = load_json(STEP7B)

    disc_attr = pd.read_csv(DISC_ATTR).sort_values("species").reset_index(drop=True)
    disc_sp = pd.read_csv(DISC_SP).sort_values("species").reset_index(drop=True)
    res_sp = pd.read_csv(RES_SP).sort_values("species").reset_index(drop=True)
    dpanel = pd.read_csv(STEP9_DISC_PANEL).sort_values("species").reset_index(drop=True)
    rpanel = pd.read_csv(STEP9_RES_PANEL).sort_values("species").reset_index(drop=True)

    assert len(disc_attr) == 369 and disc_attr["species"].nunique() == 369
    assert len(disc_sp) == 369 and disc_sp["species"].nunique() == 369
    assert len(res_sp) == 363 and res_sp["species"].nunique() == 363
    assert len(dpanel) == 369 and len(rpanel) == 363
    assert s8pre["row_count"] == 50000 and s8pre["unique_species"] == 500
    assert s8pre["every_species_exactly_100_rows"] is True

    # Figure 1: estimand / species-level D distribution.
    fig1 = disc_attr[["species", "D", "second_fraction", "n_classifiable"]].copy()
    fig1["second_ge_0_10"] = fig1["second_fraction"] >= 0.10
    fig1["second_ge_0_20"] = fig1["second_fraction"] >= 0.20
    fig1.to_csv(FIG / "figure1_discovery_species_polymorphism.csv", index=False)

    # Figure 2: discovery and reserve species-level spatial relationships.
    f2d = disc_sp[["species", "D", "second_fraction", "spatial_observed_rho", "primary_second_ge_0_10", "sensitivity_second_ge_0_20"]].copy()
    f2d["tranche"] = "discovery"
    f2d = f2d.rename(columns={"spatial_observed_rho": "spatial_rho"})
    f2r = res_sp[["species", "D", "second_fraction", "rho_primary", "rho_observer_pair_exclusion", "rho_calendar_quarter_stratification", "rho_matched_background_differential", "primary_second_ge_0_10", "sensitivity_second_ge_0_20"]].copy()
    f2r["tranche"] = "reserve"
    f2r = f2r.rename(columns={"rho_primary": "spatial_rho"})
    f2d.to_csv(FIG / "figure2_discovery_spatial.csv", index=False)
    f2r.to_csv(FIG / "figure2_reserve_spatial.csv", index=False)

    effect_rows = [
        {"tranche": "discovery", "metric": "whole_frame_mean_rho", "estimate": s5["upstream"]["whole_frame_observed_mean_rho"], "p": s5["upstream"]["whole_frame_p_upper"], "n": s5["upstream"]["whole_frame_n_species"]},
        {"tranche": "discovery", "metric": "second_ge_0.10_mean_rho", "estimate": s5["primary"]["observed_mean_rho"], "p": s5["primary"]["p_upper"], "n": s5["primary"]["n_species"]},
        {"tranche": "discovery", "metric": "second_lt_0.10_mean_rho", "estimate": s5["complement"]["observed_mean_rho"], "p": s5["complement"]["p_upper"], "n": s5["complement"]["n_species"]},
        {"tranche": "discovery", "metric": "polymorphic_minus_complement", "estimate": s5["dilution_contrast"]["primary_minus_complement_observed_delta"], "p": s5["dilution_contrast"]["p_upper"], "n": 369},
        {"tranche": "discovery", "metric": "rho_D_spatial", "estimate": s5["continuous_D_diagnostic"]["rho_D_species_spatial_rho"], "p": s5["continuous_D_diagnostic"]["p_upper"], "n": 369},
        {"tranche": "reserve", "metric": "second_ge_0.10_primary_mean_rho", "estimate": s5b["primary_subset"]["primary"]["observed_mean_rho"], "p": s5b["primary_subset"]["primary"]["p_upper"], "n": s5b["primary_subset"]["primary"]["n_species"]},
        {"tranche": "reserve", "metric": "second_ge_0.10_observer_exclusion_mean_rho", "estimate": s5b["primary_subset"]["observer_pair_exclusion"]["observed_mean_rho"], "p": s5b["primary_subset"]["observer_pair_exclusion"]["p_upper"], "n": s5b["primary_subset"]["observer_pair_exclusion"]["n_species"]},
        {"tranche": "reserve", "metric": "second_ge_0.10_quarter_mean_rho", "estimate": s5b["primary_subset"]["calendar_quarter_stratification"]["observed_mean_rho"], "p": s5b["primary_subset"]["calendar_quarter_stratification"]["p_upper"], "n": s5b["primary_subset"]["calendar_quarter_stratification"]["n_species"]},
        {"tranche": "reserve", "metric": "second_ge_0.10_matched_background_mean_rho", "estimate": s5b["primary_subset"]["matched_background_differential"]["observed_mean_rho"], "p": s5b["primary_subset"]["matched_background_differential"]["p_upper"], "n": s5b["primary_subset"]["matched_background_differential"]["n_species"]},
        {"tranche": "reserve", "metric": "polymorphic_minus_complement", "estimate": s5b["dilution_contrast"]["primary_minus_complement_observed_delta"], "p": s5b["dilution_contrast"]["p_upper"], "n": 363},
        {"tranche": "reserve", "metric": "rho_D_spatial", "estimate": s5b["continuous_D_diagnostic"]["rho_D_species_primary_spatial_rho"], "p": s5b["continuous_D_diagnostic"]["p_upper"], "n": 363},
    ]
    pd.DataFrame(effect_rows).to_csv(FIG / "figure2_effect_summary.csv", index=False)

    # Figure 3: measurement decomposition and endpoint robustness.
    measurement_counts = s8pre["candidate_status_fields"]["measurement_status"]["value_counts"]
    pd.DataFrame([{"measurement_status": k, "count": v, "fraction": v / s8pre["row_count"]} for k, v in measurement_counts.items()]).to_csv(
        FIG / "figure3_measurement_status_partition.csv", index=False
    )

    robust_rows = [
        {"tranche": "discovery", "response": "primary_spatial", "adjustment": "span_plus_technical_failure", "D_variant": "observed", "partial_rho": s8["discovery_spatial"]["D"]["observed_partial_spearman"], "p": s8["discovery_spatial"]["D"]["p_upper_geometry_preserving_spatial_null"]},
        {"tranche": "reserve", "response": "primary_spatial", "adjustment": "span_plus_technical_failure", "D_variant": "observed", "partial_rho": s8["reserve_spatial_primary"]["D"]["observed_partial_spearman"], "p": s8["reserve_spatial_primary"]["D"]["p_upper_geometry_preserving_spatial_null"]},
        {"tranche": "reserve", "response": "flower_minus_background", "adjustment": "span_plus_technical_failure", "D_variant": "observed", "partial_rho": s8["reserve_spatial_matched_background_differential"]["D"]["observed_partial_spearman"], "p": s8["reserve_spatial_matched_background_differential"]["D"]["p_upper_geometry_preserving_spatial_null"]},
    ]
    for tranche, block in [("discovery", s9["discovery_spatial"]), ("reserve", s9["reserve_spatial_primary"])]:
        for dv in ["D_min4", "D_max4"]:
            robust_rows.append({"tranche": tranche, "response": "primary_spatial", "adjustment": "span_plus_technical_failure", "D_variant": dv, "partial_rho": block[dv]["observed_partial_spearman"], "p": block[dv]["p_upper_geometry_preserving_spatial_null"]})
    for dv in ["D_min4", "D_max4"]:
        block = s9["reserve_spatial_matched_background_differential"][dv]
        robust_rows.append({"tranche": "reserve", "response": "flower_minus_background", "adjustment": "span_plus_technical_failure", "D_variant": dv, "partial_rho": block["observed_partial_spearman"], "p": block["p_upper_geometry_preserving_spatial_null"]})
    pd.DataFrame(robust_rows).to_csv(FIG / "figure3_spatial_robustness_summary.csv", index=False)

    interval_rows = []
    for tranche, desc in [("discovery", s9["discovery_description"]), ("reserve", s9["reserve_description"])]:
        interval_rows.append({"tranche": tranche, **desc})
    pd.DataFrame(interval_rows).to_csv(FIG / "figure3_ambiguity_interval_summary.csv", index=False)
    dpanel[["species", "D", "D_min4", "D_max4", "D_interval_width", "ambiguous_palette_rate", "technical_failure_rate"]].to_csv(FIG / "figure3_discovery_ambiguity_intervals.csv", index=False)
    rpanel[["species", "D", "D_min4", "D_max4", "D_interval_width", "ambiguous_palette_rate", "technical_failure_rate"]].to_csv(FIG / "figure3_reserve_ambiguity_intervals.csv", index=False)

    # Figure 4: secondary genus signal.
    genus_rows = [
        {"tranche": "discovery", "variant": "observed_raw", "gain": s4["primary"]["genus_taxonomic_clustering"]["clustering_gain"], "p": s4["primary"]["genus_taxonomic_clustering"]["p_lower"], "repeated_genera": s4["primary"]["genus_taxonomic_clustering"]["repeated_groups"], "species": s4["primary"]["genus_taxonomic_clustering"]["species_in_repeated_groups"]},
        {"tranche": "reserve", "variant": "observed_raw", "gain": s7["reserve"]["raw"]["clustering_gain"], "p": s7["reserve"]["raw"]["p_lower"], "repeated_genera": s7["reserve"]["raw"]["repeated_groups"], "species": s7["reserve"]["raw"]["species_in_repeated_groups"]},
        {"tranche": "reserve", "variant": "span_only", "gain": s7b["reserve"]["span_only_adjusted_D"]["clustering_gain"], "p": s7b["reserve"]["span_only_adjusted_D"]["p_lower"], "repeated_genera": s7b["reserve"]["span_only_adjusted_D"]["repeated_groups"], "species": s7b["reserve"]["span_only_adjusted_D"]["species_in_repeated_groups"]},
        {"tranche": "reserve", "variant": "span_plus_technical", "gain": s8["reserve_genus"]["span_plus_technical_adjusted_D"]["clustering_gain"], "p": s8["reserve_genus"]["span_plus_technical_adjusted_D"]["p_lower"], "repeated_genera": s8["reserve_genus"]["span_plus_technical_adjusted_D"]["repeated_groups"], "species": s8["reserve_genus"]["span_plus_technical_adjusted_D"]["species_in_repeated_groups"]},
        {"tranche": "reserve", "variant": "D_min4_span_plus_technical", "gain": s9["reserve_genus"]["D_min4"]["clustering_gain"], "p": s9["reserve_genus"]["D_min4"]["p_lower"], "repeated_genera": s9["reserve_genus"]["D_min4"]["repeated_groups"], "species": s9["reserve_genus"]["D_min4"]["species_in_repeated_groups"]},
        {"tranche": "reserve", "variant": "D_max4_span_plus_technical", "gain": s9["reserve_genus"]["D_max4"]["clustering_gain"], "p": s9["reserve_genus"]["D_max4"]["p_lower"], "repeated_genera": s9["reserve_genus"]["D_max4"]["repeated_groups"], "species": s9["reserve_genus"]["D_max4"]["species_in_repeated_groups"]},
    ]
    pd.DataFrame(genus_rows).to_csv(FIG / "figure4_genus_clustering_summary.csv", index=False)

    # Central synchronized number ledger.
    numbers = {
        "paper_version": "v0.1-post-step9",
        "title": "Flower-colour polymorphism is geographically organized within species without a shared global boundary",
        "discovery": {
            "species": 369,
            "D_min": s1["fingerprint"]["D_min"],
            "D_max": s1["fingerprint"]["D_max"],
            "second_ge_0_10_fraction": s1["fingerprint"]["second_ge_0_10"],
            "second_ge_0_20_fraction": s1["fingerprint"]["second_ge_0_20"],
            "rho_D_spatial": s5["continuous_D_diagnostic"]["rho_D_species_spatial_rho"],
            "rho_D_spatial_p": s5["continuous_D_diagnostic"]["p_upper"],
            "technical_adjusted_partial_rho": s8["discovery_spatial"]["D"]["observed_partial_spearman"],
            "technical_adjusted_p": s8["discovery_spatial"]["D"]["p_upper_geometry_preserving_spatial_null"],
        },
        "reserve": {
            "species": 363,
            "rho_D_spatial": s5b["continuous_D_diagnostic"]["rho_D_species_primary_spatial_rho"],
            "rho_D_spatial_p": s5b["continuous_D_diagnostic"]["p_upper"],
            "technical_adjusted_primary_partial_rho": s8["reserve_spatial_primary"]["D"]["observed_partial_spearman"],
            "technical_adjusted_primary_p": s8["reserve_spatial_primary"]["D"]["p_upper_geometry_preserving_spatial_null"],
            "technical_adjusted_background_partial_rho": s8["reserve_spatial_matched_background_differential"]["D"]["observed_partial_spearman"],
            "technical_adjusted_background_p": s8["reserve_spatial_matched_background_differential"]["D"]["p_upper_geometry_preserving_spatial_null"],
            "Dmin_primary_partial_rho": s9["reserve_spatial_primary"]["D_min4"]["observed_partial_spearman"],
            "Dmin_primary_p": s9["reserve_spatial_primary"]["D_min4"]["p_upper_geometry_preserving_spatial_null"],
            "Dmax_primary_partial_rho": s9["reserve_spatial_primary"]["D_max4"]["observed_partial_spearman"],
            "Dmax_primary_p": s9["reserve_spatial_primary"]["D_max4"]["p_upper_geometry_preserving_spatial_null"],
            "Dmin_background_partial_rho": s9["reserve_spatial_matched_background_differential"]["D_min4"]["observed_partial_spearman"],
            "Dmin_background_p": s9["reserve_spatial_matched_background_differential"]["D_min4"]["p_upper_geometry_preserving_spatial_null"],
            "Dmax_background_partial_rho": s9["reserve_spatial_matched_background_differential"]["D_max4"]["observed_partial_spearman"],
            "Dmax_background_p": s9["reserve_spatial_matched_background_differential"]["D_max4"]["p_upper_geometry_preserving_spatial_null"],
        },
        "measurement": {
            "raw_rows": s8pre["row_count"],
            "candidate_species": s8pre["unique_species"],
            "rows_per_species": 100,
            "measurement_status_counts": measurement_counts,
            "discovery_rho_D_technical_failure": s8["discovery_diagnostics"]["rho_D_technical_failure_rate"],
            "reserve_rho_D_technical_failure": s8["reserve_diagnostics"]["rho_D_technical_failure_rate"],
            "discovery_rho_D_ambiguous_palette": s8["discovery_diagnostics"]["rho_D_ambiguous_palette_rate"],
            "reserve_rho_D_ambiguous_palette": s8["reserve_diagnostics"]["rho_D_ambiguous_palette_rate"],
            "discovery_ambiguity_median_width": s9["discovery_description"]["median_interval_width"],
            "reserve_ambiguity_median_width": s9["reserve_description"]["median_interval_width"],
        },
        "genus_secondary": {
            "discovery_gain": s4["primary"]["genus_taxonomic_clustering"]["clustering_gain"],
            "discovery_p": s4["primary"]["genus_taxonomic_clustering"]["p_lower"],
            "discovery_holm_p": s4["decisions"]["genus_taxonomic_clustering"]["holm_p"],
            "reserve_raw_gain": s7["reserve"]["raw"]["clustering_gain"],
            "reserve_raw_p": s7["reserve"]["raw"]["p_lower"],
            "shared_23_rho": s7["cross_tranche_genus_mean_concordance"]["rho"],
            "shared_23_p": s7["cross_tranche_genus_mean_concordance"]["p_two_sided"],
            "reserve_Dmin_gain": s9["reserve_genus"]["D_min4"]["clustering_gain"],
            "reserve_Dmin_p": s9["reserve_genus"]["D_min4"]["p_lower"],
            "reserve_Dmax_gain": s9["reserve_genus"]["D_max4"]["clustering_gain"],
            "reserve_Dmax_p": s9["reserve_genus"]["D_max4"]["p_lower"],
        },
        "rejected_or_nonheadline": {
            "directionality_claim_survives": s1["decision"]["claim2_survives"],
            "directionality_rotation_p": s1["direction_specificity"]["p_rotation_upper"],
            "step4_span_rho": s4["primary"]["sampled_geographic_span"]["rho"],
            "step4_span_holm_p": s4["decisions"]["sampled_geographic_span"]["holm_p"],
            "reserve_span_rho": s8["reserve_diagnostics"]["rho_D_sampled_span"],
            "absolute_latitude_holm_p": s4["decisions"]["absolute_latitude_centroid"]["holm_p"],
            "family_holm_p": s4["decisions"]["family_taxonomic_clustering"]["holm_p"],
        },
    }
    write_json(OUT / "paper_numbers.json", numbers)

    md = [
        "# FCP paper v0.1 synchronized number ledger",
        "",
        "Generated from frozen result receipts. Do not hand-edit numerical values in this file.",
        "",
        "## Discovery",
        "",
        f"- species: **{numbers['discovery']['species']}**",
        f"- D range: **{numbers['discovery']['D_min']:.6f}–{numbers['discovery']['D_max']:.6f}**",
        f"- second morph >=10%: **{100*numbers['discovery']['second_ge_0_10_fraction']:.2f}%**",
        f"- second morph >=20%: **{100*numbers['discovery']['second_ge_0_20_fraction']:.2f}%**",
        f"- rho(D, spatial): **{numbers['discovery']['rho_D_spatial']:.6f}**, p **{numbers['discovery']['rho_D_spatial_p']:.6g}**",
        "",
        "## Reserve flagship robustness",
        "",
        f"- species: **{numbers['reserve']['species']}**",
        f"- raw rho(D, primary spatial): **{numbers['reserve']['rho_D_spatial']:.6f}**, p **{numbers['reserve']['rho_D_spatial_p']:.6g}**",
        f"- span+technical adjusted primary: **{numbers['reserve']['technical_adjusted_primary_partial_rho']:.6f}**, p **{numbers['reserve']['technical_adjusted_primary_p']:.6g}**",
        f"- span+technical adjusted flower-minus-background: **{numbers['reserve']['technical_adjusted_background_partial_rho']:.6f}**, p **{numbers['reserve']['technical_adjusted_background_p']:.6g}**",
        f"- D_min4 primary: **{numbers['reserve']['Dmin_primary_partial_rho']:.6f}**, p **{numbers['reserve']['Dmin_primary_p']:.6g}**",
        f"- D_max4 primary: **{numbers['reserve']['Dmax_primary_partial_rho']:.6f}**, p **{numbers['reserve']['Dmax_primary_p']:.6g}**",
        f"- D_min4 flower-minus-background: **{numbers['reserve']['Dmin_background_partial_rho']:.6f}**, p **{numbers['reserve']['Dmin_background_p']:.6g}**",
        f"- D_max4 flower-minus-background: **{numbers['reserve']['Dmax_background_partial_rho']:.6f}**, p **{numbers['reserve']['Dmax_background_p']:.6g}**",
        "",
        "## Secondary genus result",
        "",
        f"- discovery gain **{numbers['genus_secondary']['discovery_gain']:.6f}**, raw p **{numbers['genus_secondary']['discovery_p']:.6g}**, Holm p **{numbers['genus_secondary']['discovery_holm_p']:.6g}**",
        f"- reserve raw gain **{numbers['genus_secondary']['reserve_raw_gain']:.6f}**, p **{numbers['genus_secondary']['reserve_raw_p']:.6g}**",
        f"- shared 23-genus concordance rho **{numbers['genus_secondary']['shared_23_rho']:.6f}**, p **{numbers['genus_secondary']['shared_23_p']:.6g}**",
        f"- reserve D_min4 gain **{numbers['genus_secondary']['reserve_Dmin_gain']:.6f}**, p **{numbers['genus_secondary']['reserve_Dmin_p']:.6g}**",
        f"- reserve D_max4 gain **{numbers['genus_secondary']['reserve_Dmax_gain']:.6f}**, p **{numbers['genus_secondary']['reserve_Dmax_p']:.6g}**",
        "",
        "## Measurement",
        "",
        f"- raw ledger: **{numbers['measurement']['raw_rows']}** photos across **{numbers['measurement']['candidate_species']}** candidate species, exactly **100 photos/species**",
        f"- discovery rho(D, technical failure): **{numbers['measurement']['discovery_rho_D_technical_failure']:.6f}**",
        f"- reserve rho(D, technical failure): **{numbers['measurement']['reserve_rho_D_technical_failure']:.6f}**",
        f"- discovery rho(D, ambiguous palette): **{numbers['measurement']['discovery_rho_D_ambiguous_palette']:.6f}**",
        f"- reserve rho(D, ambiguous palette): **{numbers['measurement']['reserve_rho_D_ambiguous_palette']:.6f}**",
        f"- ambiguity interval median width: discovery **{numbers['measurement']['discovery_ambiguity_median_width']:.6f}**, reserve **{numbers['measurement']['reserve_ambiguity_median_width']:.6f}**",
    ]
    (OUT / "PAPER_NUMBERS.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    manifest = {
        "paper_version": "v0.1-post-step9",
        "figure_data": sorted(p.name for p in FIG.glob("*.csv")),
        "number_ledger": "paper_numbers.json",
        "human_readable_number_ledger": "PAPER_NUMBERS.md",
        "source_receipts": [str(p.relative_to(ROOT)) for p in [STEP1, STEP4, STEP5, STEP5B, STEP8_PRE, STEP8, STEP9, STEP7, STEP7B]],
    }
    write_json(OUT / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
