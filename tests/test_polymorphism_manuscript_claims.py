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
        "0.5172457461",
        "0.4571428150",
        "86",
        "0.5329282123",
        "0.4593196659",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
        "species-disjoint third cohort",
        "same iNaturalist opportunity universe",
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
        "499 species × 100 rows = 49,900",
        "measurement-evaluable species = **377**",
        "N = **158**",
        "W_obs = **0.5172457461**",
        "N = **86**",
        "W_obs = **0.5329282123**",
    ):
        assert token in text


def test_readme_routes_to_active_polymorphism_paper() -> None:
    text = README.read_text(encoding="utf-8")

    for token in (
        "# FCP — Global flower-colour polymorphism",
        "docs/POLYMORPHISM_MANUSCRIPT.md",
        "docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md",
        "docs/POLYMORPHISM_FIGURE_PLAN_20260918.md",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
    ):
        assert token in text

    assert "The active research mainline is **RGFCA**" not in text
