import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "data" / "run_jbi_ch1_vision_repeatability.py"
spec = importlib.util.spec_from_file_location("vision_repeatability", SCRIPT)
assert spec is not None and spec.loader is not None
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def record(blind_id, species, pass_index, *, condition="fresh", pattern="approximately_uniform", terms=None):
    return {
        "blind_id": blind_id,
        "species": species,
        "pass_index": pass_index,
        "flower_visibility": "evaluable",
        "visibility_failure_code": "",
        "flower_condition": condition,
        "flower_region": "single_target_clear",
        "within_photo_flower_consistency": "single_flower",
        "segmentation_feasibility": "feasible",
        "target_flower_bbox_pct": [10, 10, 90, 90],
        "primary_petals_visible": True,
        "apparent_petals_colour_terms": terms or ["purple"],
        "colour_pattern": pattern,
        "diagnostic_colour_scope": "petal_background",
        "confidence": 0.9,
        "notes": "test",
    }


def test_repeatability_summary_reports_unanimity_without_thresholding():
    rows = []
    for image_index in range(6):
        for pass_index in range(1, 4):
            rows.append(record(f"id{image_index}", f"Species {image_index}", pass_index))
    summary = mod.summarize(rows)
    assert summary["n_images"] == 6
    assert summary["n_valid_responses"] == 18
    assert summary["field_repeatability"]["flower_visibility"]["unanimous_fraction"] == 1.0
    assert summary["scaleup_decision"].startswith("requires_review")


def test_repeatability_summary_exposes_condition_instability():
    rows = []
    for image_index in range(6):
        for pass_index in range(1, 4):
            condition = "senescent" if image_index == 0 and pass_index == 3 else "fresh"
            rows.append(record(f"id{image_index}", f"Species {image_index}", pass_index, condition=condition))
    summary = mod.summarize(rows)
    metric = summary["field_repeatability"]["flower_condition"]
    assert metric["n_unanimous"] == 5
    assert metric["unanimous_fraction"] == 5 / 6
    first = next(row for row in summary["per_image"] if row["blind_id"] == "id0")
    assert first["flower_condition_unanimous"] is False
    assert first["flower_condition_mode_fraction"] == 2 / 3


def test_colour_term_jaccard_is_reported_separately_from_state_agreement():
    rows = [
        record("id0", "Species 0", 1, terms=["purple", "magenta"]),
        record("id0", "Species 0", 2, terms=["purple"]),
        record("id0", "Species 0", 3, terms=["violet", "purple"]),
    ]
    for image_index in range(1, 6):
        for pass_index in range(1, 4):
            rows.append(record(f"id{image_index}", f"Species {image_index}", pass_index))
    summary = mod.summarize(rows)
    first = next(row for row in summary["per_image"] if row["blind_id"] == "id0")
    assert 0 < first["colour_terms_pairwise_jaccard"] < 1
    assert first["flower_visibility_unanimous"] is True
