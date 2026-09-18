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
    separation_ratio: float | None


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


def _fcp_initial_pair(x: np.ndarray) -> tuple[int, int, float]:
    """Frozen FCP initialization: farthest from grand mean, then farthest from it."""
    grand = x.mean(axis=0)
    i0 = int(np.argmax(np.sum((x - grand) ** 2, axis=1)))
    d0 = np.sum((x - x[i0]) ** 2, axis=1)
    i1 = int(np.argmax(d0))
    return i0, i1, float(d0[i1])


def two_mode_axis(
    compositions: np.ndarray,
    *,
    max_iter: int = 200,
    eps: float = 1e-12,
) -> TwoModeResult:
    """Fit the deterministic unlabeled two-means used by the frozen FCP H2 analysis.

    Clustering is performed on Hellinger-transformed normalized rows. The returned
    displacement is oriented from the larger cluster to the smaller cluster, as
    in the FCP implementation. Downstream squared-alignment statistics are sign
    invariant.
    """
    p = _normalize_rows(compositions)
    x = np.sqrt(p)
    i0, i1, initial_separation = _fcp_initial_pair(x)

    if initial_separation <= float(eps):
        labels = np.zeros(len(x), dtype=int)
        labels[len(x) // 2 :] = 1
    else:
        centers = np.vstack([x[i0], x[i1]])
        labels = np.full(len(x), -1, dtype=int)
        for _ in range(int(max_iter)):
            dist = np.sum((x[:, None, :] - centers[None, :, :]) ** 2, axis=2)
            new_labels = np.argmin(dist, axis=1).astype(int)
            if np.all(new_labels == new_labels[0]):
                only = int(new_labels[0])
                other = 1 - only
                far = int(np.argmax(dist[:, only]))
                new_labels[far] = other
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for k in (0, 1):
                members = x[labels == k]
                if len(members) == 0:
                    raise RuntimeError("deterministic two-means produced an empty cluster")
                centers[k] = members.mean(axis=0)
        else:
            raise RuntimeError("deterministic two-means did not converge")

    counts = np.bincount(labels, minlength=2)
    centers_final = np.vstack([x[labels == k].mean(axis=0) for k in (0, 1)])
    between = float(np.linalg.norm(centers_final[0] - centers_final[1]))
    within_ss = float(
        sum(np.sum((x[labels == k] - centers_final[k]) ** 2) for k in (0, 1))
    )
    within_rms = float(np.sqrt(within_ss / len(x))) if within_ss > float(eps) else 0.0
    separation_ratio = None if within_rms <= float(eps) else float(between / within_rms)

    if counts[0] > counts[1]:
        major, minor = 0, 1
    elif counts[1] > counts[0]:
        major, minor = 1, 0
    else:
        c0 = p[labels == 0].mean(axis=0)
        c1 = p[labels == 1].mean(axis=0)
        # Frozen tie orientation: lexicographically smaller composition is major.
        if tuple(c0.tolist()) <= tuple(c1.tolist()):
            major, minor = 0, 1
        else:
            major, minor = 1, 0

    means = np.vstack([p[labels == 0].mean(axis=0), p[labels == 1].mean(axis=0)])
    delta = means[minor] - means[major]
    norm = float(np.linalg.norm(delta))
    if norm <= 1e-15:
        raise ValueError("two-mode displacement has zero norm")
    return TwoModeResult(
        unit_axis=delta / norm,
        cluster_means=means,
        cluster_sizes=(int(counts[0]), int(counts[1])),
        minor_fraction=float(min(counts) / len(labels)),
        labels=labels.copy(),
        separation_ratio=separation_ratio,
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
