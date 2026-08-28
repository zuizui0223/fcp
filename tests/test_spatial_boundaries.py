import numpy as np
import pytest

from fcp_pipeline.spatial_boundaries import (
    shared_boundary_strength,
    species_conditioned_permutation,
    validate_label_independent_detectability,
)


def test_species_conditioned_permutation_never_crosses_species():
    labels = np.array(["red", "white", "red", "blue", "blue", "yellow"], dtype=object)
    species = np.array(["a", "a", "a", "b", "b", "b"], dtype=object)

    permuted = species_conditioned_permutation(
        labels,
        species,
        rng=np.random.default_rng(1234),
    )

    assert permuted.shape == labels.shape
    for sp in np.unique(species):
        idx = species == sp
        assert sorted(permuted[idx].tolist()) == sorted(labels[idx].tolist())


def test_shared_boundary_strength_uses_label_independent_opportunity_denominator():
    boundary = np.array(
        [
            [True, False, True, False],
            [True, True, False, False],
            [False, True, True, True],
        ]
    )
    detectable = np.array(
        [
            [True, True, False, False],
            [True, False, True, False],
            [True, True, True, False],
        ]
    )

    strength, A = shared_boundary_strength(
        boundary,
        detectable,
        min_detectable_species=2,
    )

    np.testing.assert_array_equal(A, np.array([3, 2, 2, 0]))
    np.testing.assert_allclose(strength[:3], np.array([2 / 3, 1 / 2, 1 / 2]))
    assert np.isnan(strength[3])


def test_boundary_outside_detectability_mask_cannot_contribute():
    boundary = np.array([[True], [True]])
    detectable = np.array([[False], [True]])

    strength, A = shared_boundary_strength(boundary, detectable)

    np.testing.assert_array_equal(A, np.array([1]))
    np.testing.assert_allclose(strength, np.array([1.0]))


def test_detectability_must_not_change_when_labels_are_permuted():
    observed = np.array([[True, False], [True, True]])
    validate_label_independent_detectability(observed, observed.copy())

    changed = np.array([[False, False], [True, True]])
    with pytest.raises(ValueError, match=r"A\(x\)"):
        validate_label_independent_detectability(observed, changed)
