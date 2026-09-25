import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT.md"
LEDGER = ROOT / "docs" / "POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md"
FIGURE_PLAN = ROOT / "docs" / "POLYMORPHISM_FIGURE_PLAN_20260918.md"
README = ROOT / "README.md"
D_TRANSPORT = ROOT / "results" / "polymorphism_fresh_D_transport_20260925" / "result.json"
H2 = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
MEASUREMENT = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_measurement_20260917" / "result.json"
SPATIAL = ROOT / "results" / "polymorphism_spatial_organization_clue_20260918" / "result.json"
VALIDITY = ROOT / "results" / "polymorphism_h2_posthoc_validity_diagnostics_20260922" / "result.json"
WHITE_ENV = ROOT / "results" / "polymorphism_white_environment_mechanism_20260925" / "result.json"
BIO5_TRANSPORT = ROOT / "results" / "polymorphism_legacy_white_bio5_replication_20260925" / "result.json"
LINEAGE_MAP = ROOT / "docs" / "POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md"
LINEAGE_AUDIT = ROOT / "results" / "polymorphism_data_lineage_audit_20260925" / "result.json"
WHITE_ENV_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_WHITE_ENVIRONMENT_MECHANISM_PROTOCOL_20260925.md"
BIO5_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_LEGACY_WHITE_BIO5_REPLICATION_PROTOCOL_20260925.md"
H3B_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md"
H3B_RESULT = ROOT / "results" / "polymorphism_h3b_reserve_span_20260912" / "result.json"
WORLDCLIM_RELEASE_RECEIPT = ROOT / "archive" / "fcp_submission_20260925" / "WORLDCLIM_RELEASE_RECEIPT.md"


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
        "excess alignment",
        "exposure/background-context confounding",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
        "species-disjoint third cohort",
        "same iNaturalist opportunity universe",
    ):
        assert token in text



def test_postconfirmatory_validity_receipt_preserves_frozen_verdict_and_caveat() -> None:
    result = load_json(VALIDITY)
    assert result["confirmatory_verdict_changed"] is False
    assert result["frozen_verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED"
    assert result["white_in_coarse_modes"]["species_with_white_as_primary_or_secondary_coarse_morph"] == 137
    assert result["gate_reapplied_structured_null"]["null_replicates"] == 299
    assert result["gate_reapplied_structured_null"]["null_median"] == pytest.approx(0.4474950370983093)
    assert result["gate_reapplied_structured_null"]["plus_one_upper_p"] == pytest.approx(1 / 300)
    assert result["background_white_proxy"]["sample_species"] == 461
    assert result["background_white_proxy"]["wilcoxon_p"] == pytest.approx(4.6e-12)
    assert result["prospective_highlight_control"]["executed_on_third_cohort"] is False
    assert result["disttrait_full_data_audit"]["bitwise_equivalent"] is False
    assert result["disttrait_full_data_audit"]["disttrait_two_mode_axis_W"] == pytest.approx(0.5183899565314756)


def test_bio5_secondary_claim_keeps_positive_third_cohort_and_failed_transport_together() -> None:
    env = load_json(WHITE_ENV)
    transport = load_json(BIO5_TRANSPORT)
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")

    bio5 = next(x for x in env["results"] if x["variable"] == "bio5")
    assert bio5["mechanism_gate_pass"] is True
    assert bio5["wilcoxon_holm_p"] == pytest.approx(0.03544867047368517)
    assert bio5["OR_per_within_species_SD"] == pytest.approx(1.073474538999635)
    assert transport["verdict"] == "LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST"
    assert transport["discovery"]["primary_support"] is False
    assert transport["reserve"]["primary_support"] is False

    for text in (manuscript, ledger):
        assert "281" in text
        assert "0.0354" in text or "0.0354487" in text
        assert "1.073" in text
        assert "LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST" in text

    assert "universal or replicated BIO5" in ledger
    assert "does not support a common cross-cohort BIO5 rule" in manuscript



def test_reader_can_route_each_major_claim_to_provenance() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    lineage = LINEAGE_MAP.read_text(encoding="utf-8")
    audit = load_json(LINEAGE_AUDIT)
    h3b = load_json(H3B_RESULT)

    for path in (
        LINEAGE_MAP,
        LINEAGE_AUDIT,
        WHITE_ENV_PROTOCOL,
        BIO5_PROTOCOL,
        H3B_PROTOCOL,
        H3B_RESULT,
    ):
        assert path.exists()

    for token in (
        "one physical 499-species / 49,900-row measurement cohort",
        "untouched prospective H2 confirmation",
        "post-H2 secondary validity/environmental analyses",
        "5142f7951af0dde5364bb047a566d67e8c479e51",
        "7e538e5c51c05a7cc47b2fcf53eea92634c8a863",
        "10496492307",
        "10292399238",
    ):
        assert token in lineage

    assert audit["status"] == "TRACEABLE_WITH_PERMANENT_GIT_AND_RELEASE_ARCHIVES"
    assert audit["headline_assessment"]["exact_legacy_inputs_recoverable_from_immutable_git"] is True
    assert audit["headline_assessment"]["exact_third_cohort_input_recoverable_from_immutable_git_and_artifact"] is True
    assert audit["headline_assessment"]["artifact_only_intermediates_require_permanent_archive"] is False
    assert audit["headline_assessment"]["worldclim_original_archive_sha256_recorded"] is True
    assert audit["headline_assessment"]["worldclim_checksum_enforced_by_workflows"] is True
    assert audit["headline_assessment"]["worldclim_bytes_mirrored_in_git"] is False
    assert audit["headline_assessment"]["worldclim_bytes_archived_in_release"] is True
    assert audit["headline_assessment"]["all_headline_inputs_permanently_recoverable"] is True
    assert audit["grades"]["highlight_validity"] == "A"
    assert audit["grades"]["H3a"] == "A"
    assert audit["grades"]["H3b"] == "A"
    assert audit["grades"]["secondary_BIO5"] == "A"
    assert audit["grades"]["BIO5_transport"] == "A"

    assert "pre-specified secondary BIO5 association in the prospective H2 cohort" in manuscript
    assert "one physical 499-species / 49,900-row measurement dataset used in two chronologically distinct ways" in manuscript
    assert "Only after H2 was terminalized" in manuscript
    assert h3b["decision"]["verdict"] == "H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED"

    assert WORLDCLIM_RELEASE_RECEIPT.exists()
    release_receipt = WORLDCLIM_RELEASE_RECEIPT.read_text(encoding="utf-8")
    assert "fcp-worldclim-2.1-10m-20260925" in release_receipt
    assert "588322903" in release_receipt
    assert "588322905" in release_receipt


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
        "cross-species generality is established in **phenotype space**",
        "whether geographic realization is shared or species-specific across species",
        "RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md",
        "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md",
    ):
        assert token in text



def test_manuscript_adds_fresh_D_transport_without_importing_v2_counterfactuals() -> None:
    import json

    text = MANUSCRIPT.read_text(encoding="utf-8")
    result = json.loads(D_TRANSPORT.read_text(encoding="utf-8"))

    assert result["overlap_species"] == 136
    assert result["spearman_rho"] == pytest.approx(0.9681064161858759)
    assert result["lin_ccc"] == pytest.approx(0.9720478011581966)
    assert result["calibration_slope"] == pytest.approx(0.96909589220185)
    assert result["absolute_D_change"]["median"] == pytest.approx(0.014761943866009125)

    for token in ("136", "0.968", "0.972", "0.969", "0.0148", "fresh-image"):
        assert token.lower() in text.lower()

    assert "T_white" not in text
    assert "335,994" not in text
    assert "full-pipeline exposure" not in text.lower()


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
