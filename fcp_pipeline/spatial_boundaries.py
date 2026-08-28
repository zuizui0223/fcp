"""Core helpers for species-conditioned global flower-colour boundary tests.

The display layer may hide species identity, but inferential operations must not.
This module therefore keeps species conditioning explicit in the null model and
keeps boundary detectability separate from observed colour labels.
"""

from __future__ import annotations

import numpy as np


def species_conditioned_permutation(
    labels: np.ndarray,
    species: np.ndarray,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute flower-colour labels strictly within species.

    Observation rows and therefore locations remain fixed.  For every species,
    the multiset of labels is preserved exactly.  No label can move between
    species.
    """

    labels = np.asarray(labels)
    species = np.asarray(species)
    if labels.ndim != 1 or species.ndim != 1:
        raise ValueError("labels and species must be one-dimensional")
    if labels.shape[0] != species.shape[0]:
        raise ValueError("labels and species must have equal length")

    out = labels.copy()
    for sp in np.unique(species):
        idx = np.flatnonzero(species == sp)
        out[idx] = rng.permutation(labels[idx])
    return out


def shared_boundary_strength(
    boundary: np.ndarray,
    detectable: np.ndarray,
    *,
    min_detectable_species: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate shared-boundary strength and the opportunity denominator A(x).

    Parameters
    ----------
    boundary
        Boolean matrix with shape (n_species, n_cells). ``True`` means that a
        species-specific colour boundary is detected at a cell.
    detectable
        Boolean matrix of the same shape. ``True`` means that the species had
        enough label-independent spatial support for a boundary to have been
        detectable at that cell.  This matrix must be built without using the
        observed flower-colour labels.
    min_detectable_species
        Cells with fewer detectable species are returned as ``NaN`` rather than
        zero, because they are not sufficiently evaluable.

    Returns
    -------
    strength, A
        ``A`` is the number of detectable species at each cell. ``strength`` is
        the number of detected boundaries among detectable species divided by
        ``A``.  Cells failing the evaluability threshold are ``NaN``.
    """

    boundary = np.asarray(boundary, dtype=bool)
    detectable = np.asarray(detectable, dtype=bool)
    if boundary.ndim != 2 or detectable.ndim != 2:
        raise ValueError("boundary and detectable must be two-dimensional")
    if boundary.shape != detectable.shape:
        raise ValueError("boundary and detectable must have identical shape")
    if min_detectable_species < 1:
        raise ValueError("min_detectable_species must be >= 1")

    # A boundary outside the opportunity mask cannot contribute.
    boundary_evaluable = boundary & detectable
    A = detectable.sum(axis=0).astype(int)
    numerator = boundary_evaluable.sum(axis=0).astype(float)

    strength = np.full(A.shape, np.nan, dtype=float)
    valid = A >= min_detectable_species
    strength[valid] = numerator[valid] / A[valid]
    return strength, A


def validate_label_independent_detectability(
    detectable_observed: np.ndarray,
    detectable_permuted: np.ndarray,
) -> None:
    """Guard that detectability is invariant to colour-label permutations.

    The denominator A(x) is allowed to depend on observation geometry and
    sampling support, but not on the observed colour labels.  A pipeline can
    call this guard after a permutation to ensure that rule is respected.
    """

    observed = np.asarray(detectable_observed, dtype=bool)
    permuted = np.asarray(detectable_permuted, dtype=bool)
    if observed.shape != permuted.shape:
        raise ValueError("detectability matrices must have identical shape")
    if not np.array_equal(observed, permuted):
        raise ValueError(
            "detectability changed after label permutation; A(x) must be "
            "defined independently of flower-colour labels"
        )
