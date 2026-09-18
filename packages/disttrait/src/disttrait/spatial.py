"""Species-specific spatial organization and matched null inference."""
from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.stats import rankdata, spearmanr

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class AssociationResult:
    observed: float
    p_upper: float | None
    null: np.ndarray | None


def great_circle_pairwise_km(
    latitude: Sequence[float],
    longitude: Sequence[float],
) -> np.ndarray:
    lat = np.deg2rad(np.asarray(latitude, dtype=float))
    lon = np.deg2rad(np.asarray(longitude, dtype=float))
    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape or len(lat) < 2:
        raise ValueError("latitude and longitude must be equal one-dimensional vectors")
    if np.any(~np.isfinite(lat)) or np.any(~np.isfinite(lon)):
        raise ValueError("coordinates must be finite")
    c = np.cos(lat)
    xyz = np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    dist = np.arccos(dot) * EARTH_RADIUS_KM
    upper = np.triu_indices(len(lat), k=1)
    return dist[upper]


def _normalize_rows(matrix: np.ndarray) -> np.ndarray:
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 2:
        raise ValueError("traits must have shape (n>=2, p>=2)")
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("trait rows must be finite and nonnegative")
    mass = x.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("each trait row must have positive mass")
    return x / mass[:, None]


def jensen_shannon_pairwise(traits: np.ndarray) -> np.ndarray:
    p = _normalize_rows(traits)
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_a = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kl_b = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    jsd = np.clip(0.5 * (kl_a + kl_b), 0.0, 1.0)
    upper = np.triu_indices(len(p), k=1)
    return jsd[upper]


def _rank_pearson(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(np.asarray(x, dtype=float), method="average")
    ry = rankdata(np.asarray(y, dtype=float), method="average")
    xc = rx - rx.mean()
    yc = ry - ry.mean()
    den = float(np.linalg.norm(xc) * np.linalg.norm(yc))
    return float(np.dot(xc, yc) / den) if den > 1e-15 else float("nan")


def spatial_rho(
    latitude: Sequence[float],
    longitude: Sequence[float],
    traits: np.ndarray,
) -> float:
    """Spearman association between pairwise geographic distance and trait JSD."""
    geo = great_circle_pairwise_km(latitude, longitude)
    jsd = jensen_shannon_pairwise(traits)
    if np.ptp(geo) <= 1e-12:
        return float("nan")
    if np.ptp(jsd) <= 1e-15:
        return 0.0
    return _rank_pearson(geo, jsd)


def matched_difference_spatial_rho(
    latitude: Sequence[float],
    longitude: Sequence[float],
    focal_traits: np.ndarray,
    background_traits: np.ndarray,
) -> float:
    """Spearman(geographic distance, focal JSD - matched-background JSD)."""
    geo = great_circle_pairwise_km(latitude, longitude)
    focal = jensen_shannon_pairwise(focal_traits)
    background = jensen_shannon_pairwise(background_traits)
    if focal.shape != background.shape:
        raise ValueError("focal and background pair structures do not match")
    diff = focal - background
    if np.ptp(geo) <= 1e-12:
        return float("nan")
    if np.ptp(diff) <= 1e-15:
        return 0.0
    return _rank_pearson(geo, diff)


def _permutation_seed(seed: int, key: str, permutation_index: int) -> int:
    raw = f"{seed}|{key}|{permutation_index}".encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "little")


def spatial_permutation_null(
    latitude: Sequence[float],
    longitude: Sequence[float],
    traits: np.ndarray,
    *,
    n_permutations: int = 999,
    seed: int = 0,
    key: str = "",
) -> tuple[float, np.ndarray]:
    """Vertex-permutation null preserving coordinates and complete trait rows."""
    p = _normalize_rows(traits)
    geo = great_circle_pairwise_km(latitude, longitude)
    jsd = jensen_shannon_pairwise(p)
    observed = 0.0 if np.ptp(jsd) <= 1e-15 else _rank_pearson(geo, jsd)
    n = len(p)
    upper = np.triu_indices(n, k=1)
    rank_matrix = np.zeros((n, n), dtype=float)
    ranks = rankdata(jsd, method="average")
    rank_matrix[upper] = ranks
    rank_matrix[(upper[1], upper[0])] = ranks
    g_rank = rankdata(geo, method="average")
    gc = g_rank - g_rank.mean()
    gnorm = float(np.linalg.norm(gc))
    cc = ranks - ranks.mean()
    cnorm = float(np.linalg.norm(cc))
    null = np.empty(int(n_permutations), dtype=float)
    if cnorm <= 1e-15:
        null.fill(0.0)
    else:
        u, v = upper
        for idx in range(int(n_permutations)):
            rng = np.random.default_rng(_permutation_seed(seed, key, idx))
            perm = rng.permutation(n)
            vals = rank_matrix[perm[u], perm[v]]
            null[idx] = float(vals @ gc / (cnorm * gnorm))
    return float(observed), null


def distribution_spatial_association(
    diversity: Sequence[float],
    spatial_observed: Sequence[float],
    *,
    spatial_null: np.ndarray | None = None,
) -> AssociationResult:
    """Across-species Spearman association between diversity and spatial rho."""
    d = np.asarray(diversity, dtype=float)
    y = np.asarray(spatial_observed, dtype=float)
    if d.shape != y.shape or d.ndim != 1:
        raise ValueError("diversity and spatial_observed must be equal vectors")
    observed = float(spearmanr(d, y).statistic)
    if spatial_null is None:
        return AssociationResult(observed=observed, p_upper=None, null=None)
    null_y = np.asarray(spatial_null, dtype=float)
    if null_y.ndim != 2 or null_y.shape[0] != len(d):
        raise ValueError("spatial_null must have shape (n_species, n_permutations)")
    vals = np.asarray(
        [spearmanr(d, null_y[:, j]).statistic for j in range(null_y.shape[1])],
        dtype=float,
    )
    p = float((1 + np.sum(vals >= observed)) / (len(vals) + 1))
    return AssociationResult(observed=observed, p_upper=p, null=vals)


def partial_rank_correlation(
    x: Sequence[float],
    y: Sequence[float],
    controls: Sequence[Sequence[float]],
) -> float:
    """Partial Spearman correlation by rank residualization."""
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    if xv.ndim != 1 or yv.shape != xv.shape:
        raise ValueError("x and y must be equal one-dimensional vectors")
    cols = [np.asarray(c, dtype=float) for c in controls]
    if any(c.shape != xv.shape for c in cols):
        raise ValueError("all controls must match x")
    rx = rankdata(xv, method="average")
    ry = rankdata(yv, method="average")
    design = np.column_stack(
        [np.ones(len(xv), dtype=float)] + [rankdata(c, method="average") for c in cols]
    )
    bx, *_ = np.linalg.lstsq(design, rx, rcond=None)
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    ex = rx - design @ bx
    ey = ry - design @ by
    den = float(np.linalg.norm(ex) * np.linalg.norm(ey))
    return float(ex @ ey / den) if den > 1e-15 else float("nan")


def partial_distribution_spatial_association(
    diversity: Sequence[float],
    spatial_observed: Sequence[float],
    *,
    controls: Sequence[Sequence[float]],
    spatial_null: np.ndarray | None = None,
) -> AssociationResult:
    observed = partial_rank_correlation(diversity, spatial_observed, controls)
    if spatial_null is None:
        return AssociationResult(observed=observed, p_upper=None, null=None)
    null_y = np.asarray(spatial_null, dtype=float)
    d = np.asarray(diversity, dtype=float)
    if null_y.ndim != 2 or null_y.shape[0] != len(d):
        raise ValueError("spatial_null must have shape (n_species, n_permutations)")
    vals = np.asarray(
        [
            partial_rank_correlation(d, null_y[:, j], controls)
            for j in range(null_y.shape[1])
        ],
        dtype=float,
    )
    p = float((1 + np.sum(vals >= observed)) / (len(vals) + 1))
    return AssociationResult(observed=observed, p_upper=p, null=vals)
