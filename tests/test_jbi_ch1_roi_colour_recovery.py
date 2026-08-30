from __future__ import annotations

import numpy as np

from scripts.data.build_jbi_ch1_calibration_pseudosplit import build_pseudosplit
from scripts.data.extract_jbi_ch1_roi_colour_features import (
    find_box,
    normalize_box,
    srgb_to_lab,
    trimmed_mean,
)
from scripts.data.finalize_jbi_ch1_continuous_colour_features import (
    apply_species_scalers,
    compute_species_scalers,
)


def test_srgb_to_lab_reference_neutrals() -> None:
    rgb = np.array([[0, 0, 0], [255, 255, 255], [128, 128, 128]], dtype=np.uint8)
    lab = srgb_to_lab(rgb)
    assert lab.shape == (3, 3)
    assert np.allclose(lab[0], [0.0, 0.0, 0.0], atol=1e-5)
    assert np.allclose(lab[1], [100.0, 0.0, 0.0], atol=2e-3)
    assert 53.0 < lab[2, 0] < 54.0
    assert abs(lab[2, 1]) < 2e-3
    assert abs(lab[2, 2]) < 2e-3


def test_nested_florence_box_is_found_by_geometry_only() -> None:
    row = {
        "model_output": {
            "<OPEN_VOCABULARY_DETECTION>": {
                "bboxes": [[10, 10, 30, 30], [5, 5, 80, 70]],
                "bboxes_labels": ["flower", "flower"],
            }
        }
    }
    box, path = find_box(row)
    assert np.array_equal(box, np.array([5.0, 5.0, 80.0, 70.0]))
    assert "bboxes" in path
    assert normalize_box(box, 100, 100, path) == (5, 5, 80, 70)


def test_componentwise_trimmed_mean_is_deterministic() -> None:
    values = np.column_stack(
        (
            np.arange(100, dtype=float),
            np.arange(100, dtype=float) * 2,
            -np.arange(100, dtype=float),
        )
    )
    first = trimmed_mean(values, 0.10)
    second = trimmed_mean(values[::-1], 0.10)
    assert np.array_equal(first, second)
    assert np.allclose(first, [49.5, 99.0, -49.5])


def test_pseudosplit_recovers_all_calibration_and_fixed_padding() -> None:
    rows = []
    for index in range(80):
        rows.append(
            {
                "species": "species a",
                "photo_id": str(index),
                "split": "calibration",
                "split_rank_hash": f"c{index:03d}",
            }
        )
    for index in range(120):
        rows.append(
            {
                "species": "species a",
                "photo_id": str(1000 + index),
                "split": "evaluation",
                "split_rank_hash": f"e{index:03d}",
            }
        )
    pseudo, provenance = build_pseudosplit(rows)
    target = [row for row in pseudo if row["calibration_recovery_role"] == "target_original_calibration"]
    padding = [row for row in pseudo if row["calibration_recovery_role"] == "padding_original_evaluation"]
    unused = [row for row in pseudo if row["calibration_recovery_role"] == "unused_original_evaluation"]
    assert len(target) == 80
    assert len(padding) == 40
    assert len(unused) == 80
    assert all(row["split"] == "evaluation" for row in target + padding)
    assert all(row["split"] == "calibration" for row in unused)
    assert {row["photo_id"] for row in padding} == {str(1000 + i) for i in range(40)}
    assert provenance["species a"]["pseudo_split_counts"] == {
        "calibration": 80,
        "evaluation": 120,
    }


def test_scaler_estimation_has_no_evaluation_argument_or_dependency() -> None:
    calibration_species = ["a"] * 4 + ["b"] * 4
    calibration = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [3.0, 4.0, 5.0],
            [10.0, 11.0, 12.0],
            [11.0, 12.0, 13.0],
            [12.0, 13.0, 14.0],
            [13.0, 14.0, 15.0],
        ]
    )
    scalers = compute_species_scalers(calibration_species, calibration)

    evaluation_species = ["a", "a", "b", "b"]
    evaluation_1 = np.array([[5, 6, 7], [8, 9, 10], [20, 21, 22], [30, 31, 32]], dtype=float)
    evaluation_2 = evaluation_1 * 1_000_000.0 - 12345.0
    z1 = apply_species_scalers(evaluation_species, evaluation_1, scalers)
    z2 = apply_species_scalers(evaluation_species, evaluation_2, scalers)

    # Evaluation values change their own z scores, but the calibration scalers are fixed.
    assert not np.allclose(z1, z2)
    assert scalers == compute_species_scalers(calibration_species, calibration.copy())
    assert scalers["a"]["n_calibration"] == 4
    assert scalers["b"]["n_calibration"] == 4
