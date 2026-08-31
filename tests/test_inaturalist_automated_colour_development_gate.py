import numpy as np

from scripts.data.evaluate_inaturalist_automated_colour_development_gate import (
    TECHNICAL_FOLDS,
    fixed_fold_ids,
    repeatability_permutation,
    ridge_cv_r2,
)


def photo(species, encounter, index, vector):
    return {
        "canonical_name": species,
        "encounter_blind_id": encounter,
        "photo_blind_id": f"{encounter}-{index}",
        "automated_colour_state_status": "automated_colour_state_admitted",
        "flower_L_mean": str(vector[0]),
        "flower_a_mean": str(vector[1]),
        "flower_b_mean": str(vector[2]),
    }


def test_repeatability_permutation_detects_stable_encounters():
    rows = []
    for encounter_index in range(12):
        centre = np.array([encounter_index * 10.0, encounter_index * 3.0, -encounter_index])
        rows.append(photo("Species A", f"E{encounter_index:02d}", 1, centre))
        rows.append(photo("Species A", f"E{encounter_index:02d}", 2, centre + 0.01))
    result = repeatability_permutation(rows, "Species A", permutations=999, seed="test")
    assert result["repeatability_pass"] is True
    assert result["permutation_p_lower"] < 0.05


def test_repeatability_fails_closed_with_too_few_multi_photo_encounters():
    rows = [photo("Species A", f"E{index:02d}", 1, [index, 0, 0]) for index in range(9)]
    result = repeatability_permutation(rows, "Species A", permutations=99, seed="test")
    assert result["repeatability_pass"] is False
    assert result["repeatability_status"].startswith("not_evaluable")


def test_fixed_folds_are_balanced_and_order_invariant():
    ids = [f"E{index:02d}" for index in range(37)]
    first = fixed_fold_ids(ids, "Species A")
    reverse_ids = list(reversed(ids))
    reverse = fixed_fold_ids(reverse_ids, "Species A")
    mapped = dict(zip(reverse_ids, reverse))
    assert np.array_equal(first, np.array([mapped[value] for value in ids]))
    counts = np.bincount(first, minlength=TECHNICAL_FOLDS)
    assert counts.max() - counts.min() <= 1


def test_ridge_cv_r2_separates_technical_prediction_from_noise():
    rng = np.random.default_rng(42)
    x = rng.normal(size=(100, 6))
    folds = np.arange(100) % TECHNICAL_FOLDS
    predictable = np.column_stack((x[:, 0], x[:, 1] * 2, x[:, 2] - x[:, 3]))
    noise = rng.normal(size=(100, 3))
    assert ridge_cv_r2(x, predictable, folds) > 0.95
    assert ridge_cv_r2(x, noise, folds) < 0.20
