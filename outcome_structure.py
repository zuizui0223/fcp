"""Outcome decomposition for floral-colour polymorphism.

The invasion theorems in FCP establish persistence/invasibility at the
metapopulation level. They do NOT, by themselves, imply that both morphs are
locally mixed within the same population.

This module therefore separates:

1. species/metapopulation-level persistence of both morphs;
2. geographic mosaic structure (different patches dominated by different morphs);
3. local coexistence (both morphs present within the same patch);
4. mixed cases containing both mosaic and locally mixed populations.

The classification is descriptive and is intentionally kept separate from the
rare-invasion phase theorem.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple


@dataclass(frozen=True)
class SpatialOutcome:
    label: str
    global_frequency_a: float
    occupied_by_a: int
    occupied_by_b: int
    locally_mixed_patches: int
    patch_count: int


def classify_spatial_outcome(
    patch_frequencies_a: Sequence[float],
    *,
    presence_threshold: float = 1e-3,
    coexistence_threshold: float = 0.05,
) -> SpatialOutcome:
    """Classify realized spatial structure from patch-level morph-A frequencies.

    Parameters
    ----------
    patch_frequencies_a:
        Morph-A frequency in each sampled patch, each in [0, 1].
    presence_threshold:
        Threshold for saying a morph is present somewhere in the metapopulation.
    coexistence_threshold:
        A patch is locally mixed when both morph frequencies are at least this
        value. This threshold should be treated as an analysis/sampling choice,
        not a biological constant.

    Returns
    -------
    SpatialOutcome
        One of: A_monomorphic, B_monomorphic, geographic_mosaic,
        local_coexistence, mixed_spatial_polymorphism.
    """
    if not patch_frequencies_a:
        raise ValueError("at least one patch is required")
    if not (0 <= presence_threshold < 0.5):
        raise ValueError("presence_threshold must be in [0, 0.5)")
    if not (0 < coexistence_threshold <= 0.5):
        raise ValueError("coexistence_threshold must be in (0, 0.5]")
    if coexistence_threshold < presence_threshold:
        raise ValueError("coexistence_threshold must be >= presence_threshold")

    freqs = tuple(float(p) for p in patch_frequencies_a)
    if any(p < 0 or p > 1 for p in freqs):
        raise ValueError("patch frequencies must lie in [0, 1]")

    occupied_a = sum(p > presence_threshold for p in freqs)
    occupied_b = sum((1.0 - p) > presence_threshold for p in freqs)
    mixed = sum(
        coexistence_threshold <= p <= 1.0 - coexistence_threshold
        for p in freqs
    )

    global_a = sum(freqs) / len(freqs)
    a_present = global_a > presence_threshold
    b_present = (1.0 - global_a) > presence_threshold

    if a_present and not b_present:
        label = "A_monomorphic"
    elif b_present and not a_present:
        label = "B_monomorphic"
    elif mixed == len(freqs):
        label = "local_coexistence"
    elif mixed == 0:
        label = "geographic_mosaic"
    else:
        label = "mixed_spatial_polymorphism"

    return SpatialOutcome(
        label=label,
        global_frequency_a=global_a,
        occupied_by_a=occupied_a,
        occupied_by_b=occupied_b,
        locally_mixed_patches=mixed,
        patch_count=len(freqs),
    )


def has_metapopulation_polymorphism(
    patch_frequencies_a: Sequence[float], presence_threshold: float = 1e-3
) -> bool:
    """Return whether both morphs are present at the metapopulation level."""
    outcome = classify_spatial_outcome(
        patch_frequencies_a,
        presence_threshold=presence_threshold,
        coexistence_threshold=max(0.05, presence_threshold),
    )
    return outcome.label not in {"A_monomorphic", "B_monomorphic"}
