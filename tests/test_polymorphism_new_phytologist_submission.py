import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
CANONICAL = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT.md"
COVER = ROOT / "docs" / "POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md"
RGFCA_INTERPRETATION = ROOT / "docs" / "RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md"
FRAME_PROVENANCE = ROOT / "docs" / "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md"
VALIDITY = ROOT / "results" / "polymorphism_h2_posthoc_validity_diagnostics_20260922" / "result.json"
HIGHLIGHT = ROOT / "results" / "polymorphism_h2_third_cohort_highlight_validity_20260922" / "result.json"
ADJUDICATION = ROOT / "docs" / "POLYMORPHISM_H2_THIRD_COHORT_HIGHLIGHT_DECISION_ADJUDICATION_20260923.md"
D_TRANSPORT = ROOT / "results" / "polymorphism_fresh_D_transport_20260925" / "result.json"
WHITE_ENV = ROOT / "results" / "polymorphism_white_environment_mechanism_20260925" / "result.json"
BIO5_TRANSPORT = ROOT / "results" / "polymorphism_legacy_white_bio5_replication_20260925" / "result.json"
LINEAGE_MAP = ROOT / "docs" / "POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md"
WHITE_ENV_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_WHITE_ENVIRONMENT_MECHANISM_PROTOCOL_20260925.md"
WHITE_ENV_SCRIPT = ROOT / "scripts" / "analysis" / "run_white_environment_mechanism_20260925.py"
BIO5_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_LEGACY_WHITE_BIO5_REPLICATION_PROTOCOL_20260925.md"
BIO5_SCRIPT = ROOT / "scripts" / "analysis" / "run_legacy_white_bio5_replication_20260925.py"
H3B_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md"
H3B_SCRIPT = ROOT / "scripts" / "analysis" / "run_polymorphism_h3b_reserve_span_20260912.R"
H3B_RESULT = ROOT / "results" / "polymorphism_h3b_reserve_span_20260912" / "result.json"


def words(text: str) -> list[str]:
    return re.findall(r"[A-Za-z0-9]+(?:[’'-][A-Za-z0-9]+)*", text)


def test_new_phytologist_front_matter_and_summary_contract() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    title = text.splitlines()[0].removeprefix("# ").strip()
    assert len(title) <= 130

    summary = text.split("## Summary", 1)[1].split("\n\n---", 1)[0]
    bullets = [line for line in summary.splitlines() if line.startswith("- ")]
    assert len(bullets) == 4
    assert len(words(summary)) <= 200

    keyword_line = next(
        line for line in text.splitlines()
        if line.startswith("**Keywords (alphabetical):**")
    )
    keywords = [x.strip() for x in keyword_line.split("**", 2)[-1].split(":", 1)[-1].split(";") if x.strip()]
    assert 5 <= len(keywords) <= 8


def test_new_phytologist_defines_rgfca_and_preserves_the_conceptual_pivot() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "Repeated Global Flower-Colour Atlas (RGFCA)",
        "shared global boundary geography",
        "species-level polymorphism amount, colour-space geometry and within-species spatial organization",
        "From a repeated global atlas to species-level generality",
        "shared-geography estimand did not provide the positive biological spine retained here",
        "whether its spatial realization is common, partially shared or species-specific remains open",
        "Cross-species generality is established most clearly in **phenotype space**",
    ):
        assert token in text


def test_rgfca_interpretation_document_preserves_programme_boundary() -> None:
    assert RGFCA_INTERPRETATION.exists()
    text = RGFCA_INTERPRETATION.read_text(encoding="utf-8")
    for token in (
        "Repeated Global Flower-Colour Atlas",
        "balanced world-map realizations",
        "species-conditioned null",
        "primary recurrent-field G1 concentration: p = 0.070",
        "species-disjoint commonness: p = 0.856",
        "whether those spatial patterns share a common map across species remains unresolved",
        "RGFCA created the global sampling/measurement framework",
        "current paper uses within-species spatial organization as a comparative trait while leaving shared-versus-species-specific mapping open",
    ):
        assert token in text

    assert FRAME_PROVENANCE.exists()
    frame = FRAME_PROVENANCE.read_text(encoding="utf-8")
    for token in (
        "42,111 species",
        "U100 = 4,730 species",
        "third-cohort candidate universe = **3,230 species**",
        "not the denominator for estimating global polymorphism prevalence",
    ):
        assert token in frame


def test_new_phytologist_documents_42111_frame_provenance() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "42,111 unique iNaturalist species",
        "18 × 9 equal-area grid",
        "Twenty metadata-only V2 rounds",
        "3,240 fixed cell-level request attempts",
        "No candidate image pixels or flower-colour outcomes were used",
        "4,730 species",
        "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md",
    ):
        assert token in text


def test_new_phytologist_documents_original_rgfca_acquisition_and_analysis_lineage() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "Acquisition of the original discovery and reserve high-depth cohorts",
        "iNaturalist Research Grade",
        "flowering annotation (term 12, value 13)",
        "positional accuracy no worse than 5 km",
        "at most two retained photographs",
        "deterministic geographic maximin sampling",
        "1,000 species × 100 photographs",
        "did not impose a native-range restriction",
        "Table 1. Data lineage and inferential roles of the high-depth cohorts",
        "D definition/descriptives; H1 diagnostic; legacy H2 target discovery/audit; D–spatial organization; H3b discovery calibration",
        "H1 primary reliability; legacy H2 validation; D–spatial replication and robustness; H3a phylogeny; H3b reserve replication",
    ):
        assert token in text


def test_new_phytologist_required_sections_and_display_items() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for heading in (
        "## Summary",
        "## Introduction",
        "## Materials and Methods",
        "## Table 1.",
        "## Results",
        "## Discussion",
        "## Acknowledgements",
        "## Competing interests",
        "## Author contributions",
        "## Data availability",
        "## Supporting Information",
        "## References",
    ):
        assert heading in text

    assert "- Figures: 5" in text
    assert "- Tables: 1" in text



def test_new_phytologist_imports_only_bounded_fresh_D_transport() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    canonical = CANONICAL.read_text(encoding="utf-8")
    result = json.loads(D_TRANSPORT.read_text(encoding="utf-8"))

    assert result["overlap_species"] == 136
    assert result["spearman_rho"] == pytest.approx(0.9681064161858759)
    assert result["lin_ccc"] == pytest.approx(0.9720478011581966)
    assert result["calibration_slope"] == pytest.approx(0.96909589220185)
    assert result["absolute_D_change"]["median"] == pytest.approx(0.014761943866009125)

    for manuscript in (text, canonical):
        for token in (
            "136",
            "0.968",
            "0.972",
            "0.969",
            "0.0148",
            "fresh-image",
            "not independent-source replication",
        ):
            assert token.lower() in manuscript.lower()

        # FCP v2 remains a separate measurement-validity study. The current
        # paper imports only the D transport receipt.
        assert "T_white" not in manuscript
        assert "335,994" not in manuscript
        assert "full-pipeline exposure" not in manuscript.lower()


def test_new_phytologist_draft_preserves_frozen_h2_claim() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    normalized = text.replace("**", "")
    for token in (
        "49,900",
        "377",
        "158 species",
        "0.517",
        "86 species",
        "0.533",
        "p = 0.001",
        "structured-null median of 0.457",
        "exposure-coupled rather than artifact-cleared",
        "same iNaturalist opportunity universe",
        "not an independent-source replication",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
    ):
        assert token.lower() in normalized.lower()



def test_new_phytologist_preserves_postconfirmatory_validity_boundary() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    result = __import__("json").loads(VALIDITY.read_text(encoding="utf-8"))
    highlight = __import__("json").loads(HIGHLIGHT.read_text(encoding="utf-8"))
    assert result["confirmatory_verdict_changed"] is False
    assert result["inferential_decomposition"]["frozen_observed_W"] == 0.5172457461053418
    assert result["gate_reapplied_structured_null"]["plus_one_upper_p"] == 1 / 300
    assert result["background_white_proxy"]["sample_species"] == 461

    assert highlight["confirmatory_verdict_changed"] is False
    assert highlight["decision"]["state"] == "INDETERMINATE"
    assert highlight["source_identity"]["source_drift_rows"] == 0
    assert highlight["source_identity"]["reacquisition_failed_rows"] == 0
    assert highlight["high_clip"]["rows"] == 2205
    assert highlight["coupling_model"]["odds_ratio"] == __import__("pytest").approx(1.4444932850227639)
    assert highlight["coupling_model"]["or_ci_low"] == __import__("pytest").approx(1.3893320695404507)
    assert highlight["coupling_model"]["or_ci_high"] == __import__("pytest").approx(1.5018445886490097)
    assert highlight["h2_high_clip_sensitivity"]["vector_species"] == 142
    assert highlight["h2_high_clip_sensitivity"]["vector_retention"] == __import__("pytest").approx(142 / 158)
    assert highlight["h2_high_clip_sensitivity"]["observed_W"] == __import__("pytest").approx(0.5034282974026532)
    assert highlight["h2_high_clip_sensitivity"]["structured_null_upper_p"] == __import__("pytest").approx(0.001)
    assert highlight["h2_high_clip_sensitivity"]["support"] is True
    assert ADJUDICATION.exists()

    for token in (
        "increment above a coarse-state-preserving construction baseline",
        "137 (86.7%)",
        "**1.444** (95% CI **1.389–1.502**)",
        "95% CI **1.389–1.502**",
        "142",
        "89.9%",
        "INDETERMINATE",
        "0.968",
        "0.972",
        "exposure-coupled rather than artifact-cleared",
        "not a bitwise numerical reproducer",
    ):
        assert token in text


def test_new_phytologist_reports_bounded_bio5_result_and_failed_transport() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    canonical = CANONICAL.read_text(encoding="utf-8")
    env = json.loads(WHITE_ENV.read_text(encoding="utf-8"))
    transport = json.loads(BIO5_TRANSPORT.read_text(encoding="utf-8"))

    bio5 = next(x for x in env["results"] if x["variable"] == "bio5")
    assert env["eligible_species"] == 281
    assert bio5["median_delta_white_minus_nonwhite_SD"] == pytest.approx(0.0690112924805198)
    assert bio5["wilcoxon_holm_p"] == pytest.approx(0.03544867047368517)
    assert bio5["OR_per_within_species_SD"] == pytest.approx(1.073474538999635)
    assert bio5["p"] == pytest.approx(0.0009188036770296888)
    assert bio5["mechanism_gate_pass"] is True

    assert transport["verdict"] == "LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST"
    assert transport["discovery"]["eligible_species"] == 271
    assert transport["discovery"]["species_level"]["wilcoxon_two_sided_p"] == pytest.approx(0.7432522901921289)
    assert transport["reserve"]["eligible_species"] == 260
    assert transport["reserve"]["species_level"]["wilcoxon_two_sided_p"] == pytest.approx(0.054066696426422846)

    for manuscript in (text, canonical):
        for token in (
            "### Post-confirmatory environmental filter and BIO5 transport test",
            "### A pre-specified secondary BIO5 association in the prospective H2 cohort does not transport",
            "**Holm-adjusted p = 0.0354**",
            "OR = **1.073**",
            "p = **0.000919**",
            "p = **0.743**",
            "p = **0.0541**",
            "LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST",
            "does not support a common cross-cohort BIO5 rule",
        ):
            assert token in manuscript

    assert "Temperature is therefore not supported as a universal cross-species driver" in text
    assert "rather than in one universal BIO5 coefficient" in text



def test_submission_data_lineage_is_reader_traceable() -> None:
    import json

    text = MANUSCRIPT.read_text(encoding="utf-8")
    lineage = LINEAGE_MAP.read_text(encoding="utf-8")

    for path in (
        LINEAGE_MAP,
        WHITE_ENV_PROTOCOL,
        WHITE_ENV_SCRIPT,
        BIO5_PROTOCOL,
        BIO5_SCRIPT,
        H3B_PROTOCOL,
        H3B_SCRIPT,
        H3B_RESULT,
    ):
        assert path.exists(), f"missing reader-traceable provenance file: {path}"

    for token in (
        "Prospective H2 cohort",
        "Secondary environmental follow-up",
        "one physical 499-species measurement dataset used in two chronologically distinct ways",
        "post-H2 secondary validity/environmental analyses",
        "POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md",
    ):
        assert token in text

    for token in (
        "5142f7951af0dde5364bb047a566d67e8c479e51",
        "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
        "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
        "7e538e5c51c05a7cc47b2fcf53eea92634c8a863",
        "10496492307",
        "10709490106",
        "10292218669",
        "10292662493",
        "10292767459",
        "10292399238",
        "Time-limited GitHub Actions artifacts",
        "WorldClim archive bit identity",
    ):
        assert token in lineage

    h3b = json.loads(H3B_RESULT.read_text(encoding="utf-8"))
    assert h3b["decision"]["verdict"] == "H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED"
    assert h3b["source_workflow_run"] == 34677468362
    assert h3b["source_artifact_id"] == 10292399238


def test_new_phytologist_documents_exact_D_spatial_method_and_methodological_scope() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity",
        "rho_i = Spearman(d_geo_ij, d_colour_ij)",
        "999 matched within-species null values",
        "(1 + # {rho_null >= rho_obs}) / 1000",
        "Rank(D) and rank(`rho_i`) are separately residualized",
        "Spearman(d_geo_ij, d_flower_ij - d_background_ij)",
        "It is not the difference between separate flower and background Spearman coefficients",
        "The methodological contribution is architectural rather than a claim to a new standalone statistic",
    ):
        assert token in text

    methods_audit = ROOT / "docs" / "POLYMORPHISM_METHODS_CLASSIFICATION_20260918.md"
    assert methods_audit.exists()
    audit = methods_audit.read_text(encoding="utf-8")
    assert "Standard or widely used components" in audit
    assert "Study-specific design choices" in audit
    assert "The paper should not claim a wholly new statistical method" in audit


def test_new_phytologist_spatial_clue_is_reported_without_causal_upgrade() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "Greater D is associated with stronger within-species geographic colour organization",
        "partial rho = **0.0992877**, p = **0.025**",
        "partial rho = **0.1162411**, p = **0.010**",
        "structural rather than causal",
        "two-layer ecological question",
        "Wessinger & Rausher 2012",
        "Lacey 2026",
        "Narbona et al. 2026",
        "shared pigment-network architecture",
    ):
        assert token in text


def test_new_phytologist_cover_letter_exists_and_preserves_claim_boundary() -> None:
    assert COVER.exists(), "New Phytologist cover letter has not been created"
    text = COVER.read_text(encoding="utf-8")
    for token in (
        "New Phytologist",
        "Within-species flower-colour variation shows achromatic–chromatic alignment beyond coarse colour-state composition",
        "49,900",
        "377",
        "158",
        "0.517",
        "0.533",
        "same iNaturalist opportunity universe",
        "not an independent-source replication",
        "1.444",
        "INDETERMINATE",
    ):
        assert token.lower() in text.lower()


def test_manuscripts_have_no_control_character_math_corruption() -> None:
    for path in (CANONICAL, MANUSCRIPT):
        text = path.read_text(encoding="utf-8")
        bad = [
            ch for ch in text
            if ord(ch) < 32 and ch not in ("\n", "\r")
        ]
        assert not bad, f"{path.name} contains control characters: {[ord(ch) for ch in bad]}"
        assert "W = mean_i (u_i^T q_white)^2" in text
        assert "p = (1 + #(W_null >= W_obs)) / 1000" in text


def test_cover_letter_answers_three_editor_questions_within_50_words() -> None:
    text = COVER.read_text(encoding="utf-8")
    matches = re.findall(
        r"### Question [123][^\n]*\n\n(.+?)(?=\n\n### Question|\n\n## )",
        text,
        flags=re.S,
    )
    assert len(matches) == 3
    for answer in matches:
        assert len(words(answer)) <= 50


def test_new_phytologist_figure_and_supporting_legends_are_present() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for idx in range(1, 6):
        assert f"**Figure {idx}." in text
    for idx in range(1, 10):
        assert f"**Fig. S{idx}." in text
    assert text.index("## References") < text.index("## Supporting Information")


def _section_text(text: str, heading: str, stop_headings: tuple[str, ...]) -> str:
    start_marker = f"## {heading}"
    start = text.index(start_marker) + len(start_marker)
    end = len(text)
    for stop in stop_headings:
        marker = f"## {stop}"
        idx = text.find(marker, start)
        if idx != -1:
            end = min(end, idx)
    return text[start:end]


def _declared_word_count(text: str, label: str) -> int:
    match = re.search(rf"^- {re.escape(label)}: ([0-9,]+) words$", text, flags=re.M)
    assert match, f"missing declared word count for {label}"
    return int(match.group(1).replace(",", ""))


def test_new_phytologist_front_page_counts_are_live() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    intro = len(words(_section_text(text, "Introduction", ("Materials and Methods",))))
    methods = len(words(_section_text(text, "Materials and Methods", ("Table 1.", "Results"))))
    results = len(words(_section_text(text, "Results", ("Discussion",))))
    discussion = len(words(_section_text(text, "Discussion", ("Acknowledgements",))))
    main = intro + methods + results + discussion

    assert _declared_word_count(text, "Introduction") == intro
    assert _declared_word_count(text, "Materials and Methods") == methods
    assert _declared_word_count(text, "Results") == results
    assert _declared_word_count(text, "Discussion") == discussion
    assert _declared_word_count(text, "Main text (Introduction through Discussion)") == main

    assert discussion / main <= 0.30

    figure_match = re.search(r"^- Figures: ([0-9]+)$", text, flags=re.M)
    table_match = re.search(r"^- Tables: ([0-9]+)$", text, flags=re.M)
    assert figure_match and table_match
    display_items = int(figure_match.group(1)) + int(table_match.group(1))
    assert 6 <= display_items <= 8


def test_new_phytologist_keywords_are_actually_alphabetical() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    keyword_line = next(
        line for line in text.splitlines()
        if line.startswith("**Keywords (alphabetical):**")
    )
    keywords = [
        item.strip()
        for item in keyword_line.split("**", 2)[-1].split(":", 1)[-1].split(";")
        if item.strip()
    ]
    assert keywords == sorted(keywords, key=str.casefold)
