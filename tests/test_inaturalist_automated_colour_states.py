import math

import numpy as np

from scripts.data.extract_inaturalist_automated_colour_states import (
    BACKGROUND_REQUIRED_FOR_PHOTO_ADMISSION,
    MODEL_SIZE,
    candidate_weights,
    select_encounters,
    summarize_photo,
    weighted_quantile,
)


def test_select_encounters_is_species_balanced_and_order_invariant():
    rows = [
        {
            "canonical_name": species,
            "encounter_blind_id": f"{species}-{index}",
            "image_files": f"image-{species}-{index}.jpg",
        }
        for species in ("A", "B")
        for index in range(10)
    ]
    first = select_encounters(rows, 3, "seed")
    second = select_encounters(reversed(rows), 3, "seed")
    assert first == second
    assert sum(row["canonical_name"] == "A" for row in first) == 3
    assert sum(row["canonical_name"] == "B" for row in first) == 3


def test_candidate_weight_requires_positive_prompt_to_beat_negative():
    logits = np.zeros((5, 2, 2), dtype=float)
    logits[:3] = 2
    logits[3:] = -2
    prompt, ensemble, background = candidate_weights(logits)
    assert np.all(prompt > 0)
    assert np.all(ensemble > 0)
    assert np.all(background == 0)
    inverse = -logits
    _prompt, inverse_ensemble, inverse_background = candidate_weights(inverse)
    assert np.all(inverse_ensemble == 0)
    assert np.all(inverse_background > 0)


def test_weighted_quantile_is_deterministic():
    values = np.array([3.0, 1.0, 2.0])
    weights = np.array([1.0, 1.0, 2.0])
    assert weighted_quantile(values, weights, 0.5) == 2.0


def test_identical_stable_maps_pass_fixed_gate():
    assert BACKGROUND_REQUIRED_FOR_PHOTO_ADMISSION is False
    rgb = np.zeros((MODEL_SIZE, MODEL_SIZE, 3), dtype=np.uint8)
    rgb[:, :] = (200, 30, 80)
    logits = np.zeros((5, MODEL_SIZE, MODEL_SIZE), dtype=np.float32)
    logits[:3] = 3
    logits[3:] = -3
    summary = summarize_photo(rgb, logits, np.fliplr(rgb), np.flip(logits, axis=2).copy())
    assert summary["automated_colour_state_status"] == "automated_colour_state_admitted"
    assert summary["failure_reasons"] == ""
    assert summary["background_features_available"] is False
    assert math.isclose(summary["flip_delta_e"], 0.0, abs_tol=1e-6)
    assert math.isclose(summary["flip_soft_iou"], 1.0, abs_tol=1e-6)
