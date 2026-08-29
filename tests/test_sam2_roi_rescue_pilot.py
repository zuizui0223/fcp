import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "data" / "run_jbi_ch1_sam2_roi_rescue_pilot.py"
spec = importlib.util.spec_from_file_location("sam2_rescue_pilot", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def synthetic_features():
    rows = []
    species_all = module.EXPECTED_SPECIES + ["Lysimachia arvensis"]
    for species in species_all:
        for ordinal in range(1, 81):
            rows.append({
                "species": species,
                "blind_id": f"{species.replace(' ', '_')}_{ordinal:03d}",
                "photo_id": f"{species.replace(' ', '_')}_photo_{ordinal:03d}",
                "evaluation_row": False,
                "calibration_only": True,
                "final_label": False,
                "feature_status": "ok",
                "selected_bbox": [10.0, 20.0, 90.0, 100.0],
                "downloaded_from": "https://example.invalid/image.jpg",
            })
    return rows


def synthetic_review():
    rescue = {
        "Antirrhinum majus": [7],
        "Dactylorhiza sambucina": [24],
        "Gentiana lutea": [7],
        "Ipomoea purpurea": [70],
        "Raphanus sativus": [10],
        "Lysimachia arvensis": [],
    }
    species = {}
    for name, ordinals in rescue.items():
        species[name] = {
            "rescue_segment": ordinals,
            "invalid": [],
            "ambiguous": [],
            "senescent": [1] if name == "Gentiana lutea" else [],
            "damaged": [],
            "mixed_or_ambiguous": [],
            "not_evaluable": [],
        }
    return {
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "species": species,
    }


def test_selects_one_rescue_and_one_usable_fresh_control_per_rescue_species():
    rows = module.select_pilot_rows(synthetic_features(), synthetic_review())
    assert len(rows) == 10
    assert len({row["blind_id"] for row in rows}) == 10
    assert {row["species"] for row in rows} == set(module.EXPECTED_SPECIES)
    for species in module.EXPECTED_SPECIES:
        group = [row for row in rows if row["species"] == species]
        assert {row["pilot_arm"] for row in group} == {"rescue", "control"}
        rescue = next(row for row in group if row["pilot_arm"] == "rescue")
        control = next(row for row in group if row["pilot_arm"] == "control")
        assert rescue["reviewer1_roi_validity"] == "rescue_segment"
        assert control["reviewer1_roi_validity"] == "usable"
        assert control["reviewer1_condition"] == "fresh"
    gentiana_control = next(
        row for row in rows if row["species"] == "Gentiana lutea" and row["pilot_arm"] == "control"
    )
    assert gentiana_control["review_order_within_species"] == 2


def test_pilot_selection_is_colour_and_metadata_free_and_deterministic():
    features = synthetic_features()
    review = synthetic_review()
    first = module.select_pilot_rows(features, review)
    second = module.select_pilot_rows(features, review)
    assert [(r["pilot_slot"], r["blind_id"]) for r in first] == [
        (r["pilot_slot"], r["blind_id"]) for r in second
    ]
    assert sorted(row["pilot_slot"] for row in first) == list(range(1, 11))
    assert all(row["evaluation_row"] is False for row in first)
    assert all(row["final_label"] is False for row in first)


def test_decision_precedence_matches_frozen_reviewer_manifest_semantics():
    spec = {
        "rescue_segment": [3],
        "invalid": [4],
        "ambiguous": [5],
        "senescent": [3],
        "damaged": [6],
        "mixed_or_ambiguous": [7],
        "not_evaluable": [4, 5],
    }
    assert module.decision_at(spec, 1) == ("usable", "fresh")
    assert module.decision_at(spec, 3) == ("rescue_segment", "senescent")
    assert module.decision_at(spec, 4) == ("invalid", "not_evaluable")
    assert module.decision_at(spec, 5) == ("ambiguous", "not_evaluable")
    assert module.decision_at(spec, 6) == ("usable", "damaged")
    assert module.decision_at(spec, 7) == ("usable", "mixed_or_ambiguous")


def test_mask_metrics_are_descriptive_not_acceptance_thresholds():
    mask = np.zeros((100, 100), dtype=bool)
    mask[20:60, 30:70] = True
    result = module.mask_metrics(mask, [20.0, 10.0, 80.0, 70.0])
    assert result["mask_area_fraction_image"] == 0.16
    assert result["mask_bbox"] == [30, 20, 70, 60]
    assert 0 < result["mask_to_prompt_area_ratio"] < 1
    assert result["mask_pixels_inside_prompt_fraction"] == 1.0
