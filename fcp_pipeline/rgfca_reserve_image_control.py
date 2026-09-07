"""Location-free matched background palette from the unchanged ROI-v4 masks."""
from __future__ import annotations

import numpy as np

from fcp_pipeline.photo_first_measurement import nearest_palette_counts


def background_palette_counts(rgb: np.ndarray, measurement: dict) -> dict[str, int]:
    """Add counts only; never change flower segmentation or admission."""
    flower = np.asarray(measurement["flower_mask"])
    background = np.asarray(measurement["background_mask"])
    if flower.dtype != bool or background.dtype != bool:
        raise ValueError("control masks must be frozen boolean hard masks")
    if flower.shape != rgb.shape[:2] or background.shape != flower.shape:
        raise ValueError("matched control mask shape mismatch")
    if np.any(flower & background):
        raise ValueError("matched annulus includes flower pixels")
    expected = int(measurement["background_effective_pixels"])
    if int(background.sum()) != expected or expected < 100:
        raise ValueError("matched background pixel count changed or is insufficient")
    counts = nearest_palette_counts(rgb[background])
    if sum(counts.values()) != expected:
        raise ValueError("matched background palette lost pixels")
    return counts
