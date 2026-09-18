"""Species-level diversity summaries."""
from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np


def _probabilities(values: Sequence[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) == 0:
        raise ValueError("values must be a non-empty one-dimensional sequence")
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("values must be finite and nonnegative")
    total = float(x.sum())
    if total <= 0:
        raise ValueError("values must have positive total mass")
    return x / total


def gini_simpson(values: Sequence[float]) -> float:
    """Return Gini-Simpson diversity, 1 - sum(p_k^2)."""
    p = _probabilities(values)
    return float(1.0 - np.sum(np.square(p)))


def categorical_diversity(
    labels: Iterable[object],
    *,
    categories: Sequence[object] | None = None,
) -> float:
    """Calculate Gini-Simpson diversity directly from categorical observations."""
    values = list(labels)
    if not values:
        raise ValueError("labels must contain at least one observation")
    if categories is None:
        categories = list(dict.fromkeys(values))
    index = {label: i for i, label in enumerate(categories)}
    counts = np.zeros(len(index), dtype=float)
    for value in values:
        if value not in index:
            raise ValueError(f"label not present in categories: {value!r}")
        counts[index[value]] += 1.0
    return gini_simpson(counts)
