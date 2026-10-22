"""Active October 2026 New Phytologist global-4,475-population submission checks.

The historical frozen pre-extension submission is guarded independently in
test_polymorphism_new_phytologist_submission.py. These tests deliberately
validate the CURRENT main manuscript and original-source numerical receipts,
not exact strings or superseded frame denominators in the legacy manuscript.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
SCOPE = ROOT / "study_scope.json"
LEDGER = ROOT / "results" / "fcp_np_all4475_unified_species_ledger_receipt_20261021.json"
W_SOURCE = ROOT / "results" / "fcp_all4475_white_axis_structured_null_actual_20261021.json"
FATE_SOURCE = ROOT / "results" / "fcp_np_all4475_W_geography_climate_fate_20261021.json"
ABIOTIC_SOURCE = ROOT / "results" / "fcp_np_all4475_seven_abiotic_axes_verified_20261021.json"
LEGACY_W_SOURCE = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"


def source(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def manuscript() -> str:
    return MANUSCRIPT.read_text(encoding="utf-8")


def test_global_main_population_and_legacy_are_separate_sources():
    scope = source(SCOPE)
    receipt = source(LEDGER)
    population = scope["current_New_Phytologist_analysis_frame"]
    assert scope["fixed_analysis_frame"]["photographs"] == 149900
    assert scope["fixed_analysis_frame"]["named_species"] == 1499
    assert population["original_source_species"] == 4475
    assert population["original_source_photo_measurements"] == 447500
    assert population["classifiable_four_colour_photos"] == 227107
    assert population["photo_source_sha256"] == PHOTO_SHA
    assert receipt["all_species_preserved"] is True
    assert receipt["n_total_taxa"] == 4475
    assert receipt["n_original_classifiable"] == 227107
    assert receipt["n_result_columns"] >= 100
    assert receipt["n_species_H1_geometry"] == 4475
    assert receipt["n_species_H2_local50_pair"] == 3600
    assert receipt["n_species_H3_BIO5"] == 4112
    assert receipt["n_species_H3_BIO4"] == 4112
    assert receipt["n_species_H4_inherited_morph_verified"] == 0
    assert receipt["n_species_H5_same_photo_individual_molecular_validated"] == 0
    assert receipt["NA_means_not_estimable_not_zero"] is True
    assert receipt["no_genetic_history_adaptation_pigment_cost_proven"] is True


def test_full_4475_white_axis_preserves_frozen_158_prospective_verdict():
    current = source(W_SOURCE)
    values = current["result"]["H1_W"]
    assert current["genuine_photo_release_sha256"] == PHOTO_SHA
    assert values["all_4475"]["n_species"] == 4475
    assert values["all_4475"]["n_null_worlds"] == 199
    assert values["all_4475"]["observed_mean_W"] == 0.47887191104252164
    assert values["all_4475"]["null_median_W"] == 0.4320391234323064
    assert values["all_4475"]["null_upper_p"] == 0.005
    assert current["result"]["no_taxa_excluded_from_full_primary_W"] is True
    assert sum(values["cohort_" + n]["n_species"] for n in
               ("discovery", "validation", "prospective", "wave4", "new_U100")) == 4475
    frozen = source(LEGACY_W_SOURCE)
    assert frozen["decision"]["verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
    assert frozen["thresholds"]["primary_0_10"]["species"] == 158
    assert abs(frozen["thresholds"]["primary_0_10"]["observed_W"] - 0.5172457461053418) < 1e-12
    assert frozen["thresholds"]["primary_0_10"]["structured_null_upper_p"] == 0.001


def test_global_geometry_vs_fate_effect_sizes_and_uncertainty_are_authentic():
    receipt = source(FATE_SOURCE)
    rs = receipt["results"]
    space = rs["W_geometry_vs_geographic_sorting"]
    climate = rs["W_geometry_vs_photo_BIO5_slope"]
    assert receipt["all_species_in_common_parent_roster"] == 4475
    assert receipt["actual_source_species_matrix_sha256"] == source(LEDGER)["merged_species_ledger_sha256"]
    assert space["n_species"] == 3600
    assert climate["n_species"] == 4112
    assert abs(space["partial_rank_r_controlling_photo_opportunity_and_source"] - 0.033745038246833475) < 1e-12
    assert abs(climate["partial_rank_r_controlling_photo_opportunity_and_source"] + 0.009270265653967887) < 1e-12
    for record in (space, climate):
        ci = record["genus_bootstrap_95CI_conditional_fixed_residuals"]
        assert ci[0] < 0 < ci[1]
        assert record["potential_shared_source_photo_measurement_artifact"] is True
        assert record["not_test_of_genetic_IBD_IBE_history_or_selection"] is True
        assert len(record["by_disjoint_species_cohort_raw"]) == 5
    assert receipt["no_genotype_phenotype_fitness_history_causal_estimation"] is True


def test_common_global_abiotic_frame_and_seven_axis_error_control():
    a = source(ABIOTIC_SOURCE)
    assert a["population_frame_4475_species"] is True
    assert a["estimable_environment_species"] == 4112
    assert a["classified_source_photos"] == 227107
    assert a["full_source_447500_photo_measurements_sha256"] == PHOTO_SHA
    assert len(a["axes"]) == 7
    assert all(x["n_species_estimable"] == 4112 for x in a["axes"])
    assert {x["axis"] for x in a["axes"] if x["family7_excludes_zero"]} == {"BIO4","BIO5"}
    assert a["source_4475_species_not_dropped"] is True
    assert a["no_genetic_IBE_or_local_adaptation_claim"] is True


def test_active_new_phytologist_summary_style_and_figure_fidelity():
    text = manuscript()
    assert text.startswith("# Recurrent flower-colour directions")
    headline = text.splitlines()[0][2:]
    assert len(headline) <= 130
    summary = text.split("## Summary\n",1)[1].split("**Key words:**",1)[0]
    bullets = [x for x in summary.splitlines() if x.startswith("- ")]
    assert len(bullets) == 4
    words = re.findall(r"[A-Za-z0-9]+(?:[’'-][A-Za-z0-9]+)*", "\n".join(bullets))
    assert len(words) <= 200
    for token in ("4,475", "447,500", "227,107", "0.4789", "0.4320",
                  "3,600", "4,112", "0.0161", "0.0219", "−0.009"):
        assert token in summary, token
    legends = text.split("## Figure legends\n",1)[1]
    for n in range(1,6):
        assert legends.count(f"**Fig. {n}.") == 1
    scope = source(SCOPE)
    paths = scope["New_Phytologist_main_data_figures"]
    assert len(paths) == 7
    for path in paths:
        svgpath = ROOT / path
        assert svgpath.exists(), path
        # Source-owned validated SVG geometry, not an opaque Markdown figure.
        root = ET.parse(svgpath).getroot()
        assert root.tag.endswith("svg")
        relative = Path(path).relative_to("docs").as_posix()
        assert f"({relative})" in text
    assert "H4 and H5 remain unmeasured" in text
    assert "No native-range restriction or explicit captive/wild filter was imposed" in text
    assert "not a prevalence survey" in text


def test_current_article_hypothesis_falsifiability_keeps_missing_mechanisms_missing():
    s = manuscript().lower()
    assert "why should independently evolving flowering-plant lineages" in s
    for topic in ("phenotypic accessibility", "historical", "environmental",
                  "local maintenance", "petal", "whole-plant"):
        assert topic in s
    for token in ("neutral snp", "not", "fitness", "genetic", "same source"):
        assert token in s
    study = source(SCOPE)
    assert study["current_New_Phytologist_analysis_frame"][
       "historical_genetic_selection_fitness_or_molecular_cost_causes_identified"] is False
    assert "not a natural-population genetic" in s or "not genetic morph coexistence" in s


def test_main_hypotheses_use_complete_4475_frame_not_prospective_158_gate():
    """P1–P4 are tested on the common full-source registry, not a legacy subset."""
    paper = manuscript()
    summary = paper.split("## Summary\n", 1)[1].split("**Key words:**", 1)[0]
    intro = paper.split("## Introduction\n", 1)[1].split("## Materials and Methods\n", 1)[0]
    table1 = paper.split("## Table 1.", 1)[1].split("## Table 2.", 1)[0]
    methods = paper.split("## Materials and Methods\n", 1)[1].split("## Table 1.", 1)[0]
    assert "158-species" not in summary
    assert "158 evaluable species" not in table1
    assert "P1 (recurrent direction)" in intro
    assert "P2 (spatial allocation)" in intro
    assert "P3 (climatic sorting)" in intro
    assert "P4 (coupled outcomes)" in intro
    assert "4,475-species frame" in methods
    assert "post hoc hypothesis test" in methods
    assert "### Test P1:" in methods
    assert "### Test P2:" in methods
    assert "3,600" in table1 and "4,112" in table1
    assert "Historical species-disjoint prospective source checks" in paper
    # This editorial change must not erase the authenticated historical result.
    old = source(LEGACY_W_SOURCE)
    assert old["thresholds"]["primary_0_10"]["species"] == 158
    assert old["decision"]["verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
