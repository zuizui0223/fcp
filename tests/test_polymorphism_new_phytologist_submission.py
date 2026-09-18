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
