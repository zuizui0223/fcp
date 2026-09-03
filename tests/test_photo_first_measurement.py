import numpy as np

from fcp_pipeline.photo_first_measurement import (
    REFERENCE_RGB,
    classify_masked_rgb,
    flower_only_fractions,
    nearest_palette_counts,
)


def pixels(name: str, n: int) -> np.ndarray:
    return np.repeat(np.asarray(REFERENCE_RGB[name], dtype=np.uint8)[None, :], n, axis=0)


def test_reference_anchor_pixels_map_back_to_their_anchor():
    for name in REFERENCE_RGB:
        counts = nearest_palette_counts(pixels(name, 20))
        assert counts[name] == 20
        assert sum(counts.values()) == 20


def test_nuisance_pixels_do_not_create_biological_colour_mass():
    counts = nearest_palette_counts(
        np.vstack([pixels("green", 30), pixels("brown", 20), pixels("black", 10)])
    )
    fractions = flower_only_fractions(counts)
    assert sum(fractions.values()) == 0.0
    result = classify_masked_rgb(
        np.vstack([pixels("green", 50), pixels("brown", 50)]),
        minimum_mask_pixels=50,
    )
    assert result["morph"] == "mixed_uncertain"
    assert result["measurement_status"] == "not_evaluable_no_biological_palette_mass"


def test_four_coarse_groups_are_species_independent():
    cases = {
        "white": pixels("white", 120),
        "yellow_orange": np.vstack([pixels("yellow", 80), pixels("orange", 40)]),
        "red_pink": np.vstack([pixels("red", 70), pixels("pink", 50)]),
        "blue_purple": np.vstack([pixels("blue", 70), pixels("purple", 50)]),
    }
    for expected, rgb in cases.items():
        result = classify_masked_rgb(rgb, minimum_mask_pixels=100)
        assert result["morph"] == expected
        assert result["measurement_status"] == "classified_four_state_morph"


def test_ambiguous_palette_is_structural_missingness():
    rgb = np.vstack([pixels("white", 55), pixels("red", 55)])
    result = classify_masked_rgb(rgb, minimum_mask_pixels=100)
    assert result["morph"] == "mixed_uncertain"
    assert result["measurement_status"] == "not_evaluable_ambiguous_palette_composition"


def test_small_flower_mask_fails_without_reinterpreting_pixels():
    result = classify_masked_rgb(pixels("red", 99), minimum_mask_pixels=100)
    assert result["morph"] == "mixed_uncertain"
    assert result["measurement_status"] == "not_evaluable_insufficient_flower_pixels"
    assert result["mask_pixels"] == 99
