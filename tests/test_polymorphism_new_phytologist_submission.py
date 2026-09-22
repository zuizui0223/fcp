import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
CANONICAL = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT.md"
COVER = ROOT / "docs" / "POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md"
RGFCA_INTERPRETATION = ROOT / "docs" / "RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md"
FRAME_PROVENANCE = ROOT / "docs" / "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md"
VALIDITY = ROOT / "results" / "polymorphism_h2_posthoc_validity_diagnostics_20260922" / "result.json"


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
        "species-level polymorphism amount, colour-space geometry and species-specific spatial organization",
        "From a repeated global atlas to species-level generality",
        "shared-geography estimand did not provide the positive biological spine retained here",
        "the geographic realization is allowed to remain species-specific",
        "phenotype space than in geographic space",
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
        "The common rule is more evident in phenotype space than in geographic space",
        "RGFCA created the global sampling/measurement framework",
        "current paper uses that heterogeneity as the biological object of study",
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
        "exposure/background-context confounding remains unresolved",
        "same iNaturalist opportunity universe",
        "not an independent-source replication",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
    ):
        assert token.lower() in normalized.lower()



def test_new_phytologist_preserves_postconfirmatory_validity_boundary() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    result = __import__("json").loads(VALIDITY.read_text(encoding="utf-8"))
    assert result["confirmatory_verdict_changed"] is False
    assert result["inferential_decomposition"]["frozen_observed_W"] == 0.5172457461053418
    assert result["gate_reapplied_structured_null"]["plus_one_upper_p"] == 1 / 300
    assert result["background_white_proxy"]["sample_species"] == 461
    assert result["prospective_highlight_control"]["executed_on_third_cohort"] is False
    for token in (
        "increment above a coarse-state-preserving construction baseline",
        "137 (86.7%)",
        "0.0692 versus 0.0553",
        "cannot detect an artifact that acts upstream",
        "not a bitwise numerical reproducer",
    ):
        assert token in text

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
        "two-stage working model",
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
