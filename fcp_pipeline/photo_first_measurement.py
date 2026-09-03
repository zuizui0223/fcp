"""Pixel-level colour summarization for the prospective photo-first atlas.

This module is deliberately independent of image localization. The measurement
contract supplies a locked flower mask (ROI-v4); this module consumes only the
RGB pixels inside that mask and maps them to the same generic sRGB reference
palette frozen in the completed Chapter-1 Florence pipeline. No species-specific
colour map is permitted.

The nine biological palette anchors are summarized into the four H1 states by
``photo_first_atlas_v2.coarse_morph_from_palette``. Green, brown and black are
nuisance anchors and are excluded from the biological denominator. If no usable
palette mass remains, or if the dominant-group thresholds are not met, the result
is ``mixed_uncertain`` structural measurement missingness.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .photo_first_atlas_v2 import coarse_morph_from_palette


REFERENCE_RGB: dict[str, tuple[int, int, int]] = {
    "white": (245, 245, 245),
    "yellow": (245, 210, 45),
    "orange": (235, 130, 35),
    "red": (195, 40, 45),
    "pink": (235, 115, 165),
    "magenta": (180, 45, 145),
    "purple": (115, 65, 155),
    "blue": (65, 95, 185),
    "bronze": (150, 100, 55),
    "green": (75, 130, 65),
    "brown": (100, 70, 45),
    "black": (25, 25, 25),
}
NUISANCE = frozenset({"green", "brown", "black"})
BIOLOGICAL_PALETTE = tuple(name for name in REFERENCE_RGB if name not in NUISANCE)


def srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    """Convert uint8-like sRGB rows to CIELAB under the frozen D65 transform."""

    x = np.asarray(rgb, dtype=np.float64)
    if x.ndim < 2 or x.shape[-1] != 3:
        raise ValueError("rgb must end with a three-channel dimension")
    if not np.isfinite(x).all() or np.any((x < 0.0) | (x > 255.0)):
        raise ValueError("rgb values must be finite and lie in [0, 255]")
    x = x / 255.0
    x = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
    matrix = np.array(
        [
            [0.4124564, 0.3575761, 0.1804375],
            [0.2126729, 0.7151522, 0.0721750],
            [0.0193339, 0.1191920, 0.9503041],
        ],
        dtype=float,
    )
    xyz = x @ matrix.T
    xyz = xyz / np.array([0.95047, 1.0, 1.08883], dtype=float)
    eps = 216.0 / 24389.0
    kappa = 24389.0 / 27.0
    f = np.where(xyz > eps, np.cbrt(xyz), (kappa * xyz + 16.0) / 116.0)
    L = 116.0 * f[..., 1] - 16.0
    a = 500.0 * (f[..., 0] - f[..., 1])
    b = 200.0 * (f[..., 1] - f[..., 2])
    return np.stack([L, a, b], axis=-1)


REFERENCE_LAB = {
    name: srgb_to_lab(np.asarray(rgb, dtype=float).reshape(1, 3))[0]
    for name, rgb in REFERENCE_RGB.items()
}


def nearest_palette_counts(rgb_pixels: np.ndarray) -> dict[str, int]:
    """Assign each RGB flower-mask pixel to the nearest frozen palette anchor."""

    pixels = np.asarray(rgb_pixels)
    if pixels.ndim != 2 or pixels.shape[1] != 3:
        raise ValueError("rgb_pixels must have shape (n_pixels, 3)")
    if len(pixels) == 0:
        return {name: 0 for name in REFERENCE_RGB}
    lab = srgb_to_lab(pixels.astype(float, copy=False))
    names = tuple(REFERENCE_LAB)
    refs = np.stack([REFERENCE_LAB[name] for name in names])
    distances = ((lab[:, None, :] - refs[None, :, :]) ** 2).sum(axis=2)
    nearest = distances.argmin(axis=1)
    counts = np.bincount(nearest, minlength=len(names))
    return {name: int(counts[index]) for index, name in enumerate(names)}


def flower_only_fractions(counts: Mapping[str, int]) -> dict[str, float]:
    """Normalize the nine non-nuisance anchors after dropping green/brown/black."""

    missing = set(REFERENCE_RGB).difference(counts)
    if missing:
        raise ValueError(f"palette counts missing anchors: {sorted(missing)}")
    usable = {name: int(counts[name]) for name in BIOLOGICAL_PALETTE}
    if any(value < 0 for value in usable.values()):
        raise ValueError("palette counts must be non-negative")
    total = int(sum(usable.values()))
    if total <= 0:
        return {name: 0.0 for name in BIOLOGICAL_PALETTE}
    return {name: value / total for name, value in usable.items()}


def classify_masked_rgb(
    rgb_pixels: np.ndarray,
    *,
    minimum_mask_pixels: int = 100,
    minimum_dominant_fraction: float = 0.50,
    minimum_margin: float = 0.10,
) -> dict[str, object]:
    """Return fixed-palette photo colour state from already-authorized flower pixels."""

    pixels = np.asarray(rgb_pixels)
    if pixels.ndim != 2 or pixels.shape[1] != 3:
        raise ValueError("rgb_pixels must have shape (n_pixels, 3)")
    minimum_mask_pixels = int(minimum_mask_pixels)
    if minimum_mask_pixels < 1:
        raise ValueError("minimum_mask_pixels must be positive")
    if len(pixels) < minimum_mask_pixels:
        return {
            "morph": "mixed_uncertain",
            "measurement_status": "not_evaluable_insufficient_flower_pixels",
            "mask_pixels": int(len(pixels)),
            "palette_counts": {name: 0 for name in REFERENCE_RGB},
            "flower_only_fractions": {name: 0.0 for name in BIOLOGICAL_PALETTE},
        }

    counts = nearest_palette_counts(pixels)
    fractions = flower_only_fractions(counts)
    if sum(fractions.values()) <= 0.0:
        morph = "mixed_uncertain"
        status = "not_evaluable_no_biological_palette_mass"
    else:
        morph = coarse_morph_from_palette(
            fractions,
            minimum_dominant_fraction=minimum_dominant_fraction,
            minimum_margin=minimum_margin,
        )
        status = (
            "classified_four_state_morph"
            if morph != "mixed_uncertain"
            else "not_evaluable_ambiguous_palette_composition"
        )
    return {
        "morph": morph,
        "measurement_status": status,
        "mask_pixels": int(len(pixels)),
        "palette_counts": counts,
        "flower_only_fractions": fractions,
    }


__all__ = [
    "BIOLOGICAL_PALETTE",
    "NUISANCE",
    "REFERENCE_RGB",
    "classify_masked_rgb",
    "flower_only_fractions",
    "nearest_palette_counts",
    "srgb_to_lab",
]
