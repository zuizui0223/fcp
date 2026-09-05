"""Prospective community-level flower-colour assembly tests for the global atlas.

This module asks a different question from the within-species G1--G4 atlas:
when plant species are sympatric, are their flower-colour profiles more similar
(convergence) or more different (divergence) than independently matched
allopatric controls?

The matching itself is outcome-blind and is frozen upstream from independent
occurrence, climate and taxonomic metadata.  Each matched set contains one focal
species, one sympatric partner and one or more matched allopatric controls.  The
conditional randomization null treats the sympatric label as exchangeable within
that pre-frozen matched set.  This preserves the focal species and the available
partner colour distances rather than permuting arbitrary species globally.

Sign convention
---------------
    delta = D_colour(sympatric partner) - mean(D_colour(allopatric controls))

    delta < 0  -> colour convergence in sympatry
    delta > 0  -> colour divergence in sympatry

The global statistic is an equal-focal-species mean, so focal species with many
sympatric partners cannot dominate merely because they contribute more matched
sets.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class CommunityAssemblyResult:
    status: str
    mean_focal_delta: float
    median_focal_delta: float
    convergent_focal_fraction: float
    divergent_focal_fraction: float
    null_mean: float
    null_q025: float
    null_q975: float
    p_lower: float
    p_upper: float
    p_two_sided: float
    n_focal_species: int
    n_matched_sets: int
    permutations: int


def _bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment preserving input order."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p_values must be finite and lie in [0, 1]")
    if len(p) == 0:
        return p.copy()
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    m = len(p)
    adjusted = ranked * m / np.arange(1, m + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


def _prepare_sets(
    focal_species: Sequence[int],
    sympatric_colour_distance: Sequence[float],
    allopatric_control_colour_distance: Sequence[Sequence[float]],
    *,
    minimum_controls_per_set: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[np.ndarray]]:
    focal = np.asarray(focal_species)
    sym = np.asarray(sympatric_colour_distance, dtype=float)
    controls = np.asarray(allopatric_control_colour_distance, dtype=float)
    if focal.ndim != 1 or sym.ndim != 1 or len(focal) != len(sym):
        raise ValueError("focal_species and sympatric_colour_distance must be matching 1D arrays")
    if controls.ndim != 2 or controls.shape[0] != len(focal):
        raise ValueError("allopatric_control_colour_distance must be matched_sets x controls")
    minimum_controls_per_set = int(minimum_controls_per_set)
    if minimum_controls_per_set < 1:
        raise ValueError("minimum_controls_per_set must be at least one")
    if np.any(~np.isfinite(sym)):
        raise ValueError("sympatric colour distances must be finite")
    if np.any(sym < 0):
        raise ValueError("colour distances cannot be negative")
    finite_controls: list[np.ndarray] = []
    keep = np.zeros(len(focal), dtype=bool)
    for i in range(len(focal)):
        row = controls[i]
        values = row[np.isfinite(row)]
        if np.any(values < 0):
            raise ValueError("colour distances cannot be negative")
        finite_controls.append(values.astype(float, copy=False))
        if len(values) >= minimum_controls_per_set:
            keep[i] = True
    return focal[keep], sym[keep], controls[keep], [finite_controls[i] for i in np.flatnonzero(keep)]


def _aggregate_equal_focal(
    focal: np.ndarray,
    set_contrast: np.ndarray,
    *,
    minimum_sets_per_focal: int,
) -> tuple[np.ndarray, np.ndarray]:
    minimum_sets_per_focal = int(minimum_sets_per_focal)
    if minimum_sets_per_focal < 1:
        raise ValueError("minimum_sets_per_focal must be at least one")
    focal_ids: list[int] = []
    focal_delta: list[float] = []
    for focal_id in np.unique(focal):
        idx = np.flatnonzero(focal == focal_id)
        values = set_contrast[idx]
        values = values[np.isfinite(values)]
        if len(values) < minimum_sets_per_focal:
            continue
        focal_ids.append(int(focal_id))
        focal_delta.append(float(np.mean(values)))
    return np.asarray(focal_ids, dtype=np.int64), np.asarray(focal_delta, dtype=float)


def matched_sympatry_colour_assembly_test(
    *,
    focal_species: Sequence[int],
    sympatric_colour_distance: Sequence[float],
    allopatric_control_colour_distance: Sequence[Sequence[float]],
    minimum_controls_per_set: int = 2,
    minimum_sets_per_focal: int = 2,
    minimum_focal_species: int = 30,
    permutations: int = 999,
    seed: int = 20260904,
    species_fdr_alpha: float = 0.10,
) -> dict[str, object]:
    """Test and characterize flower-colour assembly in sympatry.

    A matched set is evaluable when it has one finite sympatric colour distance
    and at least ``minimum_controls_per_set`` finite allopatric-control distances.
    The observed set contrast is sympatric minus mean control distance.

    For each null permutation and each matched set, one partner distance is drawn
    uniformly from the complete set consisting of the observed sympatric partner
    plus all matched controls.  That chosen partner is treated as pseudo-sympatric
    and the remainder as pseudo-controls.  Matching membership therefore remains
    fixed while the sympatry label is randomized.

    Global inference averages matched-set contrasts within focal species and then
    gives every evaluable focal species equal weight.  Species-level labels are
    descriptive FDR-controlled characterizations, not separate causal claims.
    """
    minimum_focal_species = int(minimum_focal_species)
    permutations = int(permutations)
    if minimum_focal_species < 2:
        raise ValueError("minimum_focal_species must be at least two")
    if permutations < 19:
        raise ValueError("permutations must be at least 19")
    if not 0 < float(species_fdr_alpha) < 1:
        raise ValueError("species_fdr_alpha must lie in (0, 1)")

    focal, sym, _, control_rows = _prepare_sets(
        focal_species,
        sympatric_colour_distance,
        allopatric_control_colour_distance,
        minimum_controls_per_set=minimum_controls_per_set,
    )
    if len(focal) == 0:
        result = CommunityAssemblyResult(
            status="not_evaluable_matched_set_coverage",
            mean_focal_delta=float("nan"),
            median_focal_delta=float("nan"),
            convergent_focal_fraction=float("nan"),
            divergent_focal_fraction=float("nan"),
            null_mean=float("nan"),
            null_q025=float("nan"),
            null_q975=float("nan"),
            p_lower=float("nan"),
            p_upper=float("nan"),
            p_two_sided=float("nan"),
            n_focal_species=0,
            n_matched_sets=0,
            permutations=permutations,
        )
        return {"global": result, "species": {}, "sign_convention": "negative=convergence; positive=divergence"}

    observed_set = np.asarray(
        [sym[i] - float(np.mean(control_rows[i])) for i in range(len(focal))],
        dtype=float,
    )
    focal_ids, observed_focal = _aggregate_equal_focal(
        focal,
        observed_set,
        minimum_sets_per_focal=minimum_sets_per_focal,
    )
    n_focal = int(len(focal_ids))
    if n_focal < minimum_focal_species:
        result = CommunityAssemblyResult(
            status="not_evaluable_focal_species_coverage",
            mean_focal_delta=float("nan"),
            median_focal_delta=float("nan"),
            convergent_focal_fraction=float("nan"),
            divergent_focal_fraction=float("nan"),
            null_mean=float("nan"),
            null_q025=float("nan"),
            null_q975=float("nan"),
            p_lower=float("nan"),
            p_upper=float("nan"),
            p_two_sided=float("nan"),
            n_focal_species=n_focal,
            n_matched_sets=int(len(focal)),
            permutations=permutations,
        )
        return {"global": result, "species": {}, "sign_convention": "negative=convergence; positive=divergence"}

    focal_position = {int(fid): pos for pos, fid in enumerate(focal_ids)}
    null_focal = np.full((permutations, n_focal), np.nan, dtype=float)
    rng = np.random.default_rng(int(seed))
    for p in range(permutations):
        permuted_set = np.empty(len(focal), dtype=float)
        for i in range(len(focal)):
            pool = np.concatenate(([sym[i]], control_rows[i]))
            chosen_index = int(rng.integers(0, len(pool)))
            chosen = float(pool[chosen_index])
            remaining = np.delete(pool, chosen_index)
            permuted_set[i] = chosen - float(np.mean(remaining))
        perm_ids, perm_values = _aggregate_equal_focal(
            focal,
            permuted_set,
            minimum_sets_per_focal=minimum_sets_per_focal,
        )
        if not np.array_equal(perm_ids, focal_ids):
            raise RuntimeError("null randomization changed the evaluable focal species set")
        null_focal[p] = perm_values

    observed_mean = float(np.mean(observed_focal))
    null_global = np.mean(null_focal, axis=1)
    p_lower = float((1 + np.count_nonzero(null_global <= observed_mean)) / (permutations + 1))
    p_upper = float((1 + np.count_nonzero(null_global >= observed_mean)) / (permutations + 1))
    p_two = float(min(1.0, 2.0 * min(p_lower, p_upper)))
    result = CommunityAssemblyResult(
        status="evaluated",
        mean_focal_delta=observed_mean,
        median_focal_delta=float(np.median(observed_focal)),
        convergent_focal_fraction=float(np.mean(observed_focal < 0)),
        divergent_focal_fraction=float(np.mean(observed_focal > 0)),
        null_mean=float(np.mean(null_global)),
        null_q025=float(np.quantile(null_global, 0.025)),
        null_q975=float(np.quantile(null_global, 0.975)),
        p_lower=p_lower,
        p_upper=p_upper,
        p_two_sided=p_two,
        n_focal_species=n_focal,
        n_matched_sets=int(len(focal)),
        permutations=permutations,
    )

    species_p_two = np.empty(n_focal, dtype=float)
    species_p_lower = np.empty(n_focal, dtype=float)
    species_p_upper = np.empty(n_focal, dtype=float)
    for j in range(n_focal):
        obs = observed_focal[j]
        null_values = null_focal[:, j]
        lower = float((1 + np.count_nonzero(null_values <= obs)) / (permutations + 1))
        upper = float((1 + np.count_nonzero(null_values >= obs)) / (permutations + 1))
        species_p_lower[j] = lower
        species_p_upper[j] = upper
        species_p_two[j] = min(1.0, 2.0 * min(lower, upper))
    species_q = _bh_adjust(species_p_two)

    species_payload: dict[str, dict[str, object]] = {}
    for j, focal_id in enumerate(focal_ids):
        n_sets = int(np.count_nonzero(focal == focal_id))
        delta = float(observed_focal[j])
        if species_q[j] < species_fdr_alpha and delta < 0:
            label = "convergent_in_sympatry"
        elif species_q[j] < species_fdr_alpha and delta > 0:
            label = "divergent_in_sympatry"
        else:
            label = "undetected"
        species_payload[str(int(focal_id))] = {
            "delta_sympatric_minus_allopatric": delta,
            "n_matched_sets": n_sets,
            "p_lower": float(species_p_lower[j]),
            "p_upper": float(species_p_upper[j]),
            "p_two_sided": float(species_p_two[j]),
            "q_bh": float(species_q[j]),
            "label": label,
        }

    return {
        "global": result,
        "species": species_payload,
        "sign_convention": "negative=convergence; positive=divergence",
        "primary_inference": "two-sided deviation from matched allopatric expectation",
        "directional_interpretation": {
            "convergence": bool(observed_mean < 0 and p_lower < 0.05),
            "divergence": bool(observed_mean > 0 and p_upper < 0.05),
        },
        "species_labels_are_descriptive": True,
        "species_fdr_alpha": float(species_fdr_alpha),
    }


__all__ = [
    "CommunityAssemblyResult",
    "matched_sympatry_colour_assembly_test",
]
