import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
CANONICAL = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT.md"
COVER = ROOT / "docs" / "POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md"


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
    for token in (
        "49,900",
        "377",
        "158 species, W = 0.51725, p = 0.001",
        "86 species, W = 0.53293, p = 0.001",
        "same iNaturalist opportunity universe",
        "not an independent-source replication",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
    ):
        assert token.lower() in text.lower()


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
        "A recurrent achromatic–chromatic axis structures within-species flower-colour polymorphism across plant species",
        "49,900",
        "377",
        "158",
        "0.51725",
        "86",
        "0.53293",
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
