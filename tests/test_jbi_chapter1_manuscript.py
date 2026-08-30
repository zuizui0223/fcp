import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
MANUSCRIPT = ROOT / "docs" / "JBI_CHAPTER1_MANUSCRIPT.md"
RESULTS = ROOT / "docs" / "JBI_CHAPTER1_RESULTS.md"
PROTOCOL = ROOT / "docs" / "JBI_CHAPTER1_SPATIAL_STATE_DISTRIBUTION_PROTOCOL.md"
FIGURE_PLAN = ROOT / "docs" / "JBI_CHAPTER1_FIGURE_PLAN.md"
STAGE_A = ROOT / "docs" / "supporting" / "jbi_ch1_stage_a_continuous_graph_v1.json"
STAGE_B = ROOT / "docs" / "supporting" / "jbi_ch1_stage_b_shared_transition_concentration_v1.json"
FIGURE_MANIFEST = ROOT / "docs" / "supporting" / "jbi_ch1_figure_manifest_v1.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_frozen_stage_decisions_and_numerical_results() -> None:
    stage_a = load_json(STAGE_A)
    stage_b = load_json(STAGE_B)

    assert stage_a["status"] == "stage_a_evaluation_complete"
    assert stage_a["n_evaluation_records"] == 720
    assert stage_a["primary_k"] == 5
    assert stage_a["sensitivity_k"] == [3, 8]
    assert stage_a["n_permutations"] == 9999
    assert stage_a["primary_rejects_random_labelling_at_0_05"] is True
    assert stage_a["primary_global_result"]["observed"] == pytest.approx(1.3911414429599533)
    assert stage_a["primary_global_result"]["null_mean"] == pytest.approx(1.4294293327808507)
    assert stage_a["primary_global_result"]["standardized_clustering_deficit"] == pytest.approx(2.3113243006332413)
    assert stage_a["primary_global_result"]["p_lower_tail"] == pytest.approx(0.0113)
    assert stage_a["analyses_by_k"]["3"]["global_equal_species_mean_q"]["p_lower_tail"] == pytest.approx(0.0066)
    assert stage_a["analyses_by_k"]["8"]["global_equal_species_mean_q"]["p_lower_tail"] == pytest.approx(0.0065)

    assert stage_b["status"] == "stage_b_evaluation_complete"
    assert stage_b["n_evaluation_records"] == 720
    assert stage_b["primary_permutations"] == 9999
    assert stage_b["geometry_selection_used_colour_values"] is False
    assert stage_b["environment_used"] is False
    assert stage_b["geographic_reference_library_used"] is False
    assert stage_b["primary_rejects_shared_concentration_null_at_0_05"] is False
    assert stage_b["selected_primary_configuration"]["configuration"] == "cap_500km_grid_36x18"
    assert stage_b["primary_result"]["global_concentration"]["observed"] == pytest.approx(0.008231500853652107)
    assert stage_b["primary_result"]["global_concentration"]["null_mean"] == pytest.approx(0.005675688788581615)
    assert stage_b["primary_result"]["global_concentration"]["standardized_concentration_excess"] == pytest.approx(1.438850664549165)
    assert stage_b["primary_result"]["global_concentration"]["p_upper_tail"] == pytest.approx(0.0906)
    assert stage_b["sensitivity_results"]["cap_500km_grid_24x12"]["global_concentration"]["p_upper_tail"] == pytest.approx(0.0445)


def test_manuscript_reports_the_frozen_sample_and_results() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    for token in (
        "1,200 georeferenced photographs",
        "480 calibration and 720 held-out evaluation photographs",
        "All 720 held-out photographs were processed",
        "1.39114",
        "1.42943",
        "2.311",
        "0.0113",
        "0.0066",
        "0.0065",
        "0.0082315",
        "0.0056757",
        "1.4389",
        "0.0906",
        "0.0445",
        "500-km edge cap and 36×18",
    ):
        assert token in text


def test_manuscript_preserves_ordered_gates_and_claim_boundary() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8").lower()

    for token in (
        "stage b was run only if the primary stage-a",
        "the primary shared-concentration null was therefore not rejected",
        "did not show that independent species concentrated their strongest transitions",
        "not a universal boundary or common mechanism",
        "exploratory and cannot rescue the common-boundary hypothesis",
        "not spectrophotometric measurements",
    ):
        assert token in text

    assert (
        "environmental or historical reference layers" in text
        or "environmental or historical overlays" in text
    )

    for forbidden in (
        "the shared-concentration null was rejected",
        "a universal global boundary was confirmed",
        "all six species were individually significant",
        "flower colour was caused by climate",
        "identical selection mechanisms were demonstrated",
    ):
        assert forbidden not in text


def test_manuscript_calls_the_complete_canonical_figure_set() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    plan = FIGURE_PLAN.read_text(encoding="utf-8")

    for label in ("Figure C1", "Figure C2", "Figure C3", "Figure C4", "Figure C-S1", "Figure C-S2"):
        assert label in text
        assert label in plan


def test_results_protocol_and_manuscript_agree_on_the_realized_outcome() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    results = RESULTS.read_text(encoding="utf-8")
    protocol = PROTOCOL.read_text(encoding="utf-8")

    for token in ("p = 0.0113", "p = 0.0906"):
        assert token in results
        assert token in protocol

    assert "without a universal global transition boundary" in manuscript
    assert "does not provide confirmatory evidence" in results
    assert "Stage A primary: supported" in protocol
    assert "Stage B primary: not supported" in protocol


def test_figure_manifest_is_tied_to_the_same_frozen_decisions() -> None:
    manifest = load_json(FIGURE_MANIFEST)

    assert manifest["protocol"] == "jbi-ch1-spatial-figure-manifest-v1"
    assert manifest["status"] == "canonical_figures_generated_from_frozen_stage_a_and_stage_b"
    assert manifest["results"]["stage_a_primary_lower_tail_p"] == pytest.approx(0.0113)
    assert manifest["results"]["stage_a_primary_rejects_random_labelling_at_0_05"] is True
    assert manifest["results"]["stage_b_primary_upper_tail_p"] == pytest.approx(0.0906)
    assert manifest["results"]["stage_b_primary_rejects_shared_concentration_at_0_05"] is False
    assert manifest["environment_used_for_inference"] is False
    assert manifest["geographic_reference_library_used_for_inference"] is False
    assert len(manifest["outputs"]) == 12


def test_readme_separates_the_two_inference_lanes() -> None:
    text = README.read_text(encoding="utf-8")

    for token in (
        "two frozen inferential lanes",
        "docs/JBI_CHAPTER1_MANUSCRIPT.md",
        "docs/JBI_CHAPTER1_SPATIAL_STATUS.md",
        "docs/JBI_CHAPTER1_SPATIAL_STATE_DISTRIBUTION_PROTOCOL.md",
        "docs/JBI_CHAPTER1_RESULTS.md",
        "docs/JBI_CHAPTER1_FIGURE_PLAN.md",
        "docs/jbi_manuscript.md",
        "data/frozen/frozen_34species_five_metric_dataset.csv",
    ):
        assert token in text

    assert "Their samples, response variables, null models and claims are distinct." in text
