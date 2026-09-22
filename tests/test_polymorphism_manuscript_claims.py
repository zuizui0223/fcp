import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT.md"
LEDGER = ROOT / "docs" / "POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md"
FIGURE_PLAN = ROOT / "docs" / "POLYMORPHISM_FIGURE_PLAN_20260918.md"
README = ROOT / "README.md"
H2 = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
MEASUREMENT = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_measurement_20260917" / "result.json"
SPATIAL = ROOT / "results" / "polymorphism_spatial_organization_clue_20260918" / "result.json"\nVALIDITY = ROOT / "results" / "polymorphism_h2_postconfirmatory_validity_audit_20260922" / "result.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_authoritative_third_cohort_result_is_terminal_and_confirmed() -> None:
    h2 = load_json(H2)
    measurement = load_json(MEASUREMENT)

    assert h2["stage"] == "H2_COMPLETE"
    assert h2["status"] == "untouched_prospective_test_of_previously_frozen_axis"
    assert h2["decision"]["verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
    assert h2["decision"]["primary_confirmed"] is True
    assert h2["decision"]["strict_sensitivity_pass"] is True

    assert measurement["stage"] == "SUPPORT_GATE_COMPLETE"
    assert measurement["support_decision"] == "PASS"
    assert measurement["species"] == 499
    assert measurement["rows"] == 49900
    assert measurement["measurement_evaluable_species"] == 377
    assert measurement["partition_receipts"] == 256
    assert measurement["replacement_species"] == 0
    assert measurement["replacement_rows"] == 0
    assert measurement["persisted_image_pixels"] is False


def test_primary_and_strict_h2_numbers_match_frozen_result() -> None:
    h2 = load_json(H2)
    primary = h2["thresholds"]["primary_0_10"]
    strict = h2["thresholds"]["strict_0_20"]

    assert primary["species"] == 158
    assert primary["observed_W"] == pytest.approx(0.5172457461053418)
    assert primary["structured_null_summary"]["q50"] == pytest.approx(0.45714281500506854)
    assert primary["structured_null_upper_p"] == pytest.approx(0.001)
    assert primary["pass"] is True

    assert strict["species"] == 86
    assert strict["observed_W"] == pytest.approx(0.5329282123135909)
    assert strict["structured_null_summary"]["q50"] == pytest.approx(0.45931966586726264)
    assert strict["structured_null_upper_p"] == pytest.approx(0.001)
    assert strict["pass"] is True


def test_manuscript_reports_authoritative_third_cohort_values() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    for token in (
        "49,900",
        "377",
        "158",
        "0.517",
        "0.457",
        "86",
        "0.533",
        "0.459",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
        "species-disjoint third cohort",
        "same iNaturalist opportunity universe",
        "excess continuous alignment",
        "brightness/exposure",
    ):
        assert token in text


def test_postconfirmatory_h2_validity_audit_does_not_rewrite_frozen_verdict() -> None:
    audit = load_json(VALIDITY)
    assert audit["role"] == "post_confirmatory_validity_diagnostic"
    assert audit["changes_frozen_h2_verdict"] is False
    assert audit["frozen_h2"]["observed_W"] == pytest.approx(0.5172457461053418)
    assert audit["frozen_h2"]["structured_null_median"] == pytest.approx(0.45714281500506854)
    assert audit["white_involvement"]["white_in_primary_or_secondary_coarse_morph_n"] == 137
    assert audit["gate_reapplied_null"]["starting_coarse_gate_species"] == 185
    assert audit["gate_reapplied_null"]["median_W"] == pytest.approx(0.4474950370983093)
    assert audit["background_white_proxy"]["species_with_both_states"] == 461
    assert audit["background_white_proxy"]["wilcoxon_two_sided_p"] == pytest.approx(4.554529382155807e-12)
    assert audit["direct_highlight_control"]["third_cohort_executed"] is False


def test_spatial_organization_receipt_preserves_frozen_positive_clue() -> None:
    result = load_json(SPATIAL)
    assert result["new_biological_analysis"] is False
    assert result["source"]["pr"] == 32
    assert result["source"]["head_sha"] == "f14186590c11ac24c95e1985077908b732132e96"

    assert result["discovery"]["raw_D_spatial"]["rho"] == pytest.approx(0.08921325988911004)
    assert result["discovery"]["raw_D_spatial"]["p"] == pytest.approx(0.034)
    assert result["reserve"]["raw_D_spatial"]["rho"] == pytest.approx(0.10160084472811265)
    assert result["reserve"]["raw_D_spatial"]["p"] == pytest.approx(0.025)

    reserve = result["reserve"]
    assert reserve["span_plus_technical_adjusted_primary"]["partial_rho"] == pytest.approx(0.09928771129095708)
    assert reserve["span_plus_technical_adjusted_primary"]["p_upper_geometry_preserving_spatial_null"] == pytest.approx(0.025)
    assert reserve["span_plus_technical_adjusted_flower_minus_background"]["partial_rho"] == pytest.approx(0.1162411363016301)
    assert reserve["span_plus_technical_adjusted_flower_minus_background"]["p_upper_geometry_preserving_spatial_null"] == pytest.approx(0.01)
    assert reserve["ambiguity_endpoints_primary"]["D_min4"]["p"] == pytest.approx(0.029)
    assert reserve["ambiguity_endpoints_primary"]["D_max4"]["p"] == pytest.approx(0.008)


def test_manuscript_reports_spatial_clue_without_causal_upgrade() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "Greater D is associated with stronger within-species geographic colour organization",
        "partial rho = **0.0992877**, p = **0.025**",
        "partial rho = **0.1162411**, p = **0.010**",
        "structural rather than causal",
        "cannot distinguish among them",
    ):
        assert token in text


def test_claim_ledger_freezes_rgfca_programme_lineage() -> None:
    text = LEDGER.read_text(encoding="utf-8")
    for token in (
        "Repeated Global Flower-Colour Atlas",
        "primary recurrent-field G1 p = **0.070**",
        "species-disjoint commonness p = **0.856**",
        "generality is stronger in **phenotype space than in shared geographic location**",
        "RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md",
        "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md",
    ):
        assert token in text


def test_manuscript_preserves_h1_and_h3_boundaries() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    for token in (
        "Median split Spearman rho was **0.7891**",
        "5th percentile of **0.7652**",
        "rho = **0.7927**",
        "missed its deliberately stricter prespecified rho = 0.80 floor",
        "S1: K = **0.0710190**, p = **0.2716**",
        "S2: K = **0.0601476**, p = **0.4134**",
        "S3: K = **0.0707577**, p = **0.2674**",
        "rho = **-0.0025855**, p = **0.9586021**",
    ):
        assert token in text


def test_claim_boundary_is_explicit_in_manuscript_ledger_and_readme() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")

    for text in (manuscript, ledger, readme):
        assert "not an independent-source replication" in text.lower() or "not described as an independent-source replication" in text.lower()
        assert "global" in text.lower()
        assert "prevalence" in text.lower()

    assert "P500 supplies neither confirmation nor refutation of H2" in manuscript
    assert "P500 supplies **no durable confirmatory evidence for or against H2**" in ledger
    assert "no durable H2 biological verdict" in readme


def test_legacy_post_audit_and_prospective_confirmation_are_separated() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8").lower()
    ledger = LEDGER.read_text(encoding="utf-8").lower()

    for text in (manuscript, ledger):
        assert "after the original" in text or "after the broad" in text
        assert "prospective" in text
        assert "species-disjoint" in text

    assert "those cohorts provide discovery and target-localization evidence" in manuscript
    assert "untouched prospective confirmation" in ledger


def test_figure_plan_uses_authoritative_prospective_denominators() -> None:
    text = FIGURE_PLAN.read_text(encoding="utf-8")

    for token in (
        "499 species × 100 rows",
        "49,900 terminal rows",
        "377 measurement-evaluable",
        "N = **158**",
        "W_obs = **0.5172457461**",
        "N = **86**",
        "W_obs = **0.5329282123**",
    ):
        assert token in text


def test_readme_routes_to_active_polymorphism_paper() -> None:
    text = README.read_text(encoding="utf-8")

    for token in (
        "## Active mainline — global flower-colour polymorphism",
        "docs/POLYMORPHISM_MANUSCRIPT.md",
        "docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md",
        "docs/POLYMORPHISM_FIGURE_PLAN_20260918.md",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
    ):
        assert token in text

    assert "The active research mainline is **RGFCA**" not in text
