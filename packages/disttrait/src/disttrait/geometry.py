"""Distributional trait geometry."""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TwoModeResult:
    unit_axis: np.ndarray
    cluster_means: np.ndarray
    cluster_sizes: tuple[int, int]
    minor_fraction: float
    labels: np.ndarray


def _normalize_rows(matrix: np.ndarray) -> np.ndarray:
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 2:
        raise ValueError("matrix must have shape (n>=2, p>=2)")
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("matrix must be finite and nonnegative")
    mass = x.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("each row must have positive mass")
    return x / mass[:, None]


def hellinger_rows(matrix: np.ndarray) -> np.ndarray:
    return np.sqrt(_normalize_rows(matrix))


def _initial_farthest_pair(x: np.ndarray) -> tuple[int, int]:
    diff = x[:, None, :] - x[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    np.fill_diagonal(d2, -np.inf)
    flat = int(np.argmax(d2))
    i, j = np.unravel_index(flat, d2.shape)
    return int(min(i, j)), int(max(i, j))


def two_mode_axis(compositions: np.ndarray, *, max_iter: int = 100) -> TwoModeResult:
    """Fit deterministic unlabeled two-means in Hellinger space."""
    p = _normalize_rows(compositions)
    x = np.sqrt(p)
    i, j = _initial_farthest_pair(x)
    centers = np.vstack([x[i], x[j]])
    labels = np.full(len(x), -1, dtype=int)

    for _ in range(int(max_iter)):
        d2 = np.sum((x[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(d2, axis=1).astype(int)
        if len(np.unique(new_labels)) < 2:
            raise ValueError("two-mode clustering collapsed to one cluster")
        new_centers = np.vstack([x[new_labels == k].mean(axis=0) for k in (0, 1)])
        if np.array_equal(new_labels, labels):
            labels = new_labels
            centers = new_centers
            break
        labels = new_labels
        centers = new_centers
    else:
        raise RuntimeError("two-mode clustering did not converge")

    means = np.vstack([p[labels == k].mean(axis=0) for k in (0, 1)])
    sizes = (int(np.sum(labels == 0)), int(np.sum(labels == 1)))
    delta = means[1] - means[0]
    norm = float(np.linalg.norm(delta))
    if norm <= 1e-15:
        raise ValueError("two-mode displacement has zero norm")
    return TwoModeResult(
        unit_axis=delta / norm,
        cluster_means=means,
        cluster_sizes=sizes,
        minor_fraction=float(min(sizes) / len(labels)),
        labels=labels.copy(),
    )


def one_vs_rest_contrast(n_features: int, *, focal_index: int = 0) -> np.ndarray:
    n = int(n_features)
    if n < 2 or not (0 <= int(focal_index) < n):
        raise ValueError("invalid feature count or focal index")
    q = np.full(n, -1.0 / (n - 1), dtype=float)
    q[int(focal_index)] = 1.0
    return q / np.linalg.norm(q)


def alignment_statistic(axes: np.ndarray, contrast: Sequence[float]) -> float:
    a = np.asarray(axes, dtype=float)
    q = np.asarray(contrast, dtype=float)
    if a.ndim != 2 or q.ndim != 1 or a.shape[1] != len(q):
        raise ValueError("axes and contrast dimensions do not match")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(q)):
        raise ValueError("axes and contrast must be finite")
    norms = np.linalg.norm(a, axis=1)
    qnorm = float(np.linalg.norm(q))
    if np.any(norms <= 1e-15) or qnorm <= 1e-15:
        raise ValueError("zero-length axis or contrast")
    return float(np.mean(np.square((a / norms[:, None]) @ (q / qnorm))))


def permute_rows_within_strata(
    values: np.ndarray,
    strata: Sequence[object],
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute complete rows within fixed strata, preserving stratum counts."""
    x = np.asarray(values)
    s = np.asarray(strata, dtype=object)
    if x.ndim < 1 or len(x) != len(s):
        raise ValueError("values and strata must have the same row count")
    out = x.copy()
    for stratum in dict.fromkeys(s.tolist()):
        idx = np.flatnonzero(s == stratum)
        out[idx] = x[rng.permutation(idx)]
    return out
