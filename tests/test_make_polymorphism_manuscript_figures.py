import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "make_polymorphism_manuscript_figures.py"
SOURCE_DIR = ROOT / "results" / "polymorphism_publication_figure_source_20260918"
SOURCE_MANIFEST = SOURCE_DIR / "manifest.json"
H2_PROSPECTIVE = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
H2_PRIMARY_NULL = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "primary_0_10_structured_null.csv"
H2_STRICT_NULL = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "strict_0_20_structured_null.csv"
SPATIAL = ROOT / "results" / "polymorphism_spatial_organization_clue_20260918" / "result.json"
DISTRIBUTED = ROOT / "results" / "polymorphism_distributed_polymorphism_posthoc_20261007" / "result.json"
DISTRIBUTED_ROBUST = ROOT / "results" / "polymorphism_distributed_polymorphism_robustness_20261007" / "result.json"
IBD_IBE = ROOT / "results" / "polymorphism_phenotypic_IBD_IBE_posthoc_20261007" / "result.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module():
    assert SCRIPT.exists(), "publication figure generator is not implemented"
    spec = importlib.util.spec_from_file_location("polymorphism_figures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reporting_source_freeze_matches_artifact_derived_hashes() -> None:
    manifest = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
    assert manifest["schema"] == "polymorphism_publication_figure_source_v1"
    assert manifest["purpose"].startswith("reporting-only")

    for name, expected in manifest["files"].items():
        path = SOURCE_DIR / name
        assert path.exists()
        assert sha256(path) == expected["sha256"]

    d = pd.read_csv(SOURCE_DIR / "h3a_species_D.csv")
    assert len(d) == 732
    assert d.groupby("cohort").size().to_dict() == {"discovery": 369, "reserve": 363}


def test_prospective_h2_null_arrays_are_complete() -> None:
    result = json.loads(H2_PROSPECTIVE.read_text(encoding="utf-8"))
    primary = pd.read_csv(H2_PRIMARY_NULL)
    strict = pd.read_csv(H2_STRICT_NULL)

    assert len(primary) == 999
    assert len(strict) == 999
    assert list(primary.columns) == ["W"]
    assert list(strict.columns) == ["W"]
    assert result["decision"]["verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
    assert result["thresholds"]["primary_0_10"]["species"] == 158
    assert result["thresholds"]["strict_0_20"]["species"] == 86


def test_spatial_reporting_receipt_is_reporting_only_and_frozen() -> None:
    result = json.loads(SPATIAL.read_text(encoding="utf-8"))
    assert result["schema"] == "polymorphism_spatial_organization_reporting_receipt_v1"
    assert result["new_biological_analysis"] is False
    assert result["reserve"]["span_plus_technical_adjusted_primary"]["partial_rho"] == 0.09928771129095708
    assert result["reserve"]["span_plus_technical_adjusted_flower_minus_background"]["partial_rho"] == 0.1162411363016301


def test_generate_all_publication_figures(tmp_path: Path) -> None:
    module = load_module()
    manifest_path = module.generate_all(ROOT, tmp_path)

    expected_stems = [
        "polymorphism_figure1_measurement_frame",
        "polymorphism_figure2_h1_reproducibility",
        "polymorphism_figure3_h2_target_localization",
        "polymorphism_figure4_prospective_h2",
        "polymorphism_figure5_explanatory_boundaries",
        "polymorphism_figureS9_secondary_mechanism_evidence",
    ]
    for stem in expected_stems:
        png = tmp_path / f"{stem}.png"
        pdf = tmp_path / f"{stem}.pdf"
        assert png.exists() and png.stat().st_size > 10_000
        assert pdf.exists() and pdf.stat().st_size > 5_000

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["schema"] == "polymorphism_manuscript_figure_manifest_v1"
    assert manifest["scientific_claims_changed"] is False
    assert manifest["figures"]["figure4"]["primary"]["species"] == 158
    assert manifest["figures"]["figure4"]["primary"]["p"] == 0.001
    assert manifest["figures"]["figure4"]["strict"]["species"] == 86
    assert manifest["figures"]["figure4"]["strict"]["p"] == 0.001
    f5 = manifest["figures"]["figure5"]
    assert f5["distributed_polymorphism"]["primary_radius_km"] == 50
    assert f5["distributed_polymorphism"]["discovery_depletion"] == 0.020529254583812922
    assert f5["distributed_polymorphism"]["validation_depletion"] == 0.018672971642749295
    assert f5["distributed_polymorphism"]["third_depletion"] == 0.01468491968437185
    assert f5["distributed_polymorphism"]["discovery_p"] == 0.005
    assert f5["distributed_polymorphism"]["validation_p"] == 0.005
    assert f5["distributed_polymorphism"]["third_p"] == 0.005
    assert f5["distributed_polymorphism"]["multiscale_all_three"] is True
    assert f5["robustness"]["different_observer_all_three"] is True
    assert f5["robustness"]["nonwhite_only_all_three"] is True
    assert f5["robustness"]["continuous_nine_colour_all_three"] is True
    assert f5["ibd_ibe"]["IBE_like_all_three"] is True
    assert f5["ibd_ibe"]["IBD_like_all_three"] is True
    assert f5["ibd_ibe"]["IBE_stronger_than_IBD_all_three"] is False
    assert f5["ibd_ibe"]["cohorts"]["discovery"]["IBD_like"] == 0.027604435157351816
    assert f5["ibd_ibe"]["cohorts"]["validation"]["IBE_like"] == 0.00825562350468484
    assert f5["ibd_ibe"]["cohorts"]["third"]["IBE_like"] == 0.0101876235922715

    s9 = manifest["figures"]["supplementary_figure9"]
    assert s9["bio5"]["prospective_holm_p"] == 0.03544867047368517
    assert s9["bio5"]["observer_balanced_median_delta_SD"] == 0.05409791430882366
    assert s9["bio5"]["observer_balanced_wilcoxon_p"] == 0.0749642018520577
    assert s9["bio5"]["observer_paired_median_delta_SD"] == 0.0
    assert s9["bio5"]["observer_paired_wilcoxon_p"] == 0.4849619155258311
    assert s9["bio5"]["observer_paired_conditional_OR"] == 0.7866358392496202
    assert s9["bio5"]["discovery_p"] == 0.7432522901921289
    assert s9["bio5"]["validation_p"] == 0.054066696426422846
    assert s9["molecular_anchor"]["source_doi"] == "10.3389/fpls.2016.00204"
    assert s9["molecular_anchor"]["F3h1_dark_vs_white_fold"] == 49.0
    assert s9["molecular_anchor"]["F3h1_light_vs_white_fold"] == 42.2
    assert s9["molecular_anchor"]["F3h1_shared_significant"] is True
    assert s9["molecular_anchor"]["Myb1a_dark_vs_white_fold"] == 5.1
    assert s9["molecular_anchor"]["expanded_sequence_survey_individuals"] == 38
    assert s9["molecular_anchor"]["causal_variant_identified"] is False
    assert s9["pal_wal"]["silene_PAL_median_percent"] == 15.5
    assert s9["pal_wal"]["silene_WAL_median_percent"] == 0.21
    assert s9["pal_wal"]["cross_PAL_lower_bound_median_percent"] == 5.0
    assert s9["pal_wal"]["cross_WAL_upper_bound_median_percent"] == 0.1
    assert s9["pal_wal"]["phenotype_anchor"] == "silene_HPLC_DAD_MSn_tissue_profiles_crosssystem_source_classified"
    assert s9["claim_boundary"] == "bounded_secondary_empirical_evidence_not_universal_causation"

    f1_layout = manifest["figures"]["figure1"]["layout_contract"]
    assert f1_layout["cohort_topology"] == "inferential_sequence_not_nested_samples"
    assert f1_layout["arrow_direction"] == "top_to_bottom"
    assert f1_layout["stage_necessity"] == "shown_in_each_stage_box"
    assert f1_layout["secondary_followup"] == "dashed_post_h2_annotations"

    f2_layout = manifest["figures"]["figure2"]["layout_contract"]
    assert f2_layout["stress_annotations"] == "offset_no_legend_overlap"

    f3_layout = manifest["figures"]["figure3"]["layout_contract"]
    assert f3_layout["legend"] == "outside_below_axis"

    f4_layout = manifest["figures"]["figure4"]["layout_contract"]
    assert f4_layout["same_universe_nonreplication_note"] == "caption_not_plot_field"

    f5_layout = manifest["figures"]["figure5"]["layout_contract"]
    assert f5_layout["panel_A"] == "three_cohort_local_depletion"
    assert f5_layout["panel_B"] == "three_falsification_tests_across_three_cohorts"
    assert f5_layout["panel_C"] == "paired_IBD_IBE_effects"
