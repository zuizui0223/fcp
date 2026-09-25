"""Predeclared H2 climate-concordance helpers for the random photo-first atlas.

All climate aggregation is performed independently of flower colour. The primary
predictor is the RMS distance across z-standardized BIO1, BIO4, BIO12 and BIO15.
H2 reuses H1 edge opportunities and the same matched within-species permutation
persistence maps rather than inventing a second spatial null.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

from .photo_first_atlas import adjacent_grid_edges
from .shared_transition_surface import EqualAreaGrid, equal_area_cell_ids


BIO_VARIABLES = ("bio1", "bio4", "bio12", "bio15")
Z_VARIABLES = tuple(f"z_{name}" for name in BIO_VARIABLES)


@dataclass(frozen=True)
class ClimateConcordanceResult:
    statistic: float
    null_statistics: np.ndarray
    p_upper: float
    supported_edges: int
    edge_ids: tuple[str, ...]
    predictor: str
    subset: str


def _deterministic_mode(values: pd.Series) -> str:
    cleaned = values.dropna().astype(str)
    cleaned = cleaned[cleaned.str.len() > 0]
    if len(cleaned) == 0:
        return ""
    counts = cleaned.value_counts()
    maximum = int(counts.max())
    return sorted(counts[counts == maximum].index.astype(str))[0]


def aggregate_climate_to_h1_grid(
    source: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
) -> pd.DataFrame:
    """Aggregate an independently frozen equal-area climate grid to the H1 grid."""

    required = {latitude_col, longitude_col, *BIO_VARIABLES}
    missing = sorted(required.difference(source.columns))
    if missing:
        raise ValueError(f"climate source is missing required columns: {missing}")
    work = source.copy()
    latitude = pd.to_numeric(work[latitude_col], errors="coerce").to_numpy(dtype=float)
    longitude = pd.to_numeric(work[longitude_col], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(latitude).all() or not np.isfinite(longitude).all():
        raise ValueError("climate source coordinates must be finite")
    work["h1_cell_id"] = equal_area_cell_ids(latitude, longitude, grid)
    for variable in BIO_VARIABLES:
        work[variable] = pd.to_numeric(work[variable], errors="coerce")

    climate = (
        work.groupby("h1_cell_id", sort=True)[list(BIO_VARIABLES)]
        .mean()
        .reset_index()
        .rename(columns={"h1_cell_id": "cell_id"})
    )
    counts = work.groupby("h1_cell_id", sort=True).size().rename("source_cell_n")
    climate = climate.merge(counts, left_on="cell_id", right_index=True, how="left")

    for label in ("realm", "biome"):
        if label in work.columns:
            modes = (
                work.groupby("h1_cell_id", sort=True)[label]
                .agg(_deterministic_mode)
                .rename(label)
            )
            climate = climate.merge(
                modes, left_on="cell_id", right_index=True, how="left"
            )
        else:
            climate[label] = ""

    complete = climate[list(BIO_VARIABLES)].notna().all(axis=1)
    if int(complete.sum()) < 2:
        raise ValueError("fewer than two H1 cells have complete macroclimate")
    for variable in BIO_VARIABLES:
        values = climate.loc[complete, variable].to_numpy(dtype=float)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=0))
        if not np.isfinite(sd) or sd <= 0.0:
            raise ValueError(f"cannot standardize constant climate variable: {variable}")
        climate[f"z_{variable}"] = (
            pd.to_numeric(climate[variable], errors="coerce") - mean
        ) / sd
    climate["complete_macroclimate"] = complete
    return climate.sort_values("cell_id").reset_index(drop=True)


def build_edge_climate_contrasts(
    climate_cells: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
) -> pd.DataFrame:
    """Build fixed edge-level climate contrasts on the exact H1 graph."""

    required = {"cell_id", *Z_VARIABLES}
    missing = sorted(required.difference(climate_cells.columns))
    if missing:
        raise ValueError(f"climate cell table is missing required columns: {missing}")
    by_cell = climate_cells.set_index("cell_id", drop=False)
    rows: list[dict[str, object]] = []
    for left, right in adjacent_grid_edges(grid):
        left_i = int(left)
        right_i = int(right)
        left_row = by_cell.loc[left_i] if left_i in by_cell.index else None
        right_row = by_cell.loc[right_i] if right_i in by_cell.index else None
        row: dict[str, object] = {
            "edge_id": f"{left_i}:{right_i}",
            "cell_i": left_i,
            "cell_j": right_i,
            "climate_complete": False,
            "multivariate_climate_distance": np.nan,
        }
        for variable in BIO_VARIABLES:
            row[f"absolute_z_difference_{variable}"] = np.nan
        row["within_biome"] = False
        row["within_realm"] = False
        if left_row is not None and right_row is not None:
            left_z = np.asarray([left_row[name] for name in Z_VARIABLES], dtype=float)
            right_z = np.asarray([right_row[name] for name in Z_VARIABLES], dtype=float)
            if np.isfinite(left_z).all() and np.isfinite(right_z).all():
                delta = np.abs(left_z - right_z)
                row["climate_complete"] = True
                row["multivariate_climate_distance"] = float(
                    np.sqrt(np.mean(delta**2))
                )
                for index, variable in enumerate(BIO_VARIABLES):
                    row[f"absolute_z_difference_{variable}"] = float(delta[index])
            left_biome = str(left_row.get("biome", "") or "")
            right_biome = str(right_row.get("biome", "") or "")
            left_realm = str(left_row.get("realm", "") or "")
            right_realm = str(right_row.get("realm", "") or "")
            row["within_biome"] = bool(left_biome and left_biome == right_biome)
            row["within_realm"] = bool(left_realm and left_realm == right_realm)
        rows.append(row)
    return pd.DataFrame(rows)


def weighted_pearson(
    x: Sequence[float],
    y: Sequence[float],
    weights: Sequence[float],
) -> float:
    """Reliability-weighted Pearson correlation with positive finite weights."""

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    w_arr = np.asarray(weights, dtype=float)
    if x_arr.shape != y_arr.shape or x_arr.shape != w_arr.shape:
        raise ValueError("x, y and weights must have identical shape")
    keep = np.isfinite(x_arr) & np.isfinite(y_arr) & np.isfinite(w_arr) & (w_arr > 0)
    if int(np.count_nonzero(keep)) < 2:
        return float("nan")
    x_arr = x_arr[keep]
    y_arr = y_arr[keep]
    w_arr = w_arr[keep]
    total = float(w_arr.sum())
    x_mean = float(np.sum(w_arr * x_arr) / total)
    y_mean = float(np.sum(w_arr * y_arr) / total)
    x_centered = x_arr - x_mean
    y_centered = y_arr - y_mean
    covariance = float(np.sum(w_arr * x_centered * y_centered) / total)
    x_var = float(np.sum(w_arr * x_centered**2) / total)
    y_var = float(np.sum(w_arr * y_centered**2) / total)
    if x_var <= 0.0 or y_var <= 0.0:
        return float("nan")
    return float(covariance / np.sqrt(x_var * y_var))


def _supported_edge_frame(
    h1_edges: pd.DataFrame,
    climate_edges: pd.DataFrame,
    *,
    predictor: str,
    subset: str,
    minimum_supported_edges: int,
) -> pd.DataFrame:
    required_h1 = {"edge_id", "opportunities", "persistence"}
    if not required_h1.issubset(h1_edges.columns):
        raise ValueError("H1 edge table lacks persistence/opportunity fields")
    if predictor not in climate_edges.columns:
        raise ValueError(f"unknown climate predictor: {predictor}")
    merged = h1_edges.merge(climate_edges, on="edge_id", how="left", validate="one_to_one")
    keep = (
        pd.to_numeric(merged["opportunities"], errors="coerce").fillna(0).gt(0)
        & pd.to_numeric(merged["persistence"], errors="coerce").notna()
        & pd.to_numeric(merged[predictor], errors="coerce").notna()
    )
    if subset == "global":
        pass
    elif subset == "within_biome":
        keep &= merged["within_biome"].fillna(False).astype(bool)
    elif subset == "within_realm":
        keep &= merged["within_realm"].fillna(False).astype(bool)
    else:
        raise ValueError(f"unknown H2 subset: {subset}")
    supported = merged.loc[keep].copy()
    if len(supported) < int(minimum_supported_edges):
        raise ValueError(
            "not_evaluable_h2_supported_edges: "
            f"{len(supported)} < {int(minimum_supported_edges)}"
        )
    return supported


def climate_concordance_test(
    h1_edges: pd.DataFrame,
    null_persistence: np.ndarray,
    null_edge_ids: Sequence[str],
    climate_edges: pd.DataFrame,
    *,
    predictor: str = "multivariate_climate_distance",
    subset: str = "global",
    minimum_supported_edges: int = 20,
) -> ClimateConcordanceResult:
    """Compare observed climate alignment with matched H1 permutation maps."""

    supported = _supported_edge_frame(
        h1_edges,
        climate_edges,
        predictor=predictor,
        subset=subset,
        minimum_supported_edges=minimum_supported_edges,
    )
    edge_ids = tuple(supported["edge_id"].astype(str))
    observed = weighted_pearson(
        supported["persistence"].to_numpy(dtype=float),
        supported[predictor].to_numpy(dtype=float),
        supported["opportunities"].to_numpy(dtype=float),
    )
    if not np.isfinite(observed):
        raise ValueError("not_evaluable_h2_observed_statistic")

    null = np.asarray(null_persistence, dtype=float)
    if null.ndim != 2:
        raise ValueError("null_persistence must have shape permutations x edges")
    null_edge_ids = tuple(str(value) for value in null_edge_ids)
    if null.shape[1] != len(null_edge_ids):
        raise ValueError("null persistence edge dimension does not match edge IDs")
    position = {edge_id: index for index, edge_id in enumerate(null_edge_ids)}
    if len(position) != len(null_edge_ids):
        raise ValueError("null edge IDs must be unique")
    try:
        columns = np.asarray([position[edge_id] for edge_id in edge_ids], dtype=int)
    except KeyError as exc:
        raise ValueError(f"supported edge missing from null matrix: {exc.args[0]}") from exc

    predictor_values = supported[predictor].to_numpy(dtype=float)
    weights = supported["opportunities"].to_numpy(dtype=float)
    null_statistics = np.empty(null.shape[0], dtype=float)
    for index in range(null.shape[0]):
        null_statistics[index] = weighted_pearson(
            null[index, columns], predictor_values, weights
        )
    if not np.isfinite(null_statistics).all():
        raise ValueError("one or more H2 null statistics are not estimable")
    p_upper = float(
        (1 + np.count_nonzero(null_statistics >= observed))
        / (len(null_statistics) + 1)
    )
    return ClimateConcordanceResult(
        statistic=float(observed),
        null_statistics=null_statistics,
        p_upper=p_upper,
        supported_edges=len(supported),
        edge_ids=edge_ids,
        predictor=predictor,
        subset=subset,
    )


def holm_adjust(p_values: Sequence[float]) -> np.ndarray:
    """Holm family-wise adjusted p-values in original order."""

    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any((p < 0.0) | (p > 1.0)):
        raise ValueError("p-values must be a finite one-dimensional vector in [0, 1]")
    m = len(p)
    order = np.argsort(p, kind="stable")
    adjusted_sorted = np.empty(m, dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        value = min(1.0, (m - rank) * float(p[index]))
        running = max(running, value)
        adjusted_sorted[rank] = running
    adjusted = np.empty(m, dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted


def climate_driver_decomposition(
    h1_edges: pd.DataFrame,
    null_persistence: np.ndarray,
    null_edge_ids: Sequence[str],
    climate_edges: pd.DataFrame,
    *,
    minimum_supported_edges: int = 20,
) -> pd.DataFrame:
    """Four predeclared BIO-axis decompositions with Holm correction."""

    rows = []
    for variable in BIO_VARIABLES:
        predictor = f"absolute_z_difference_{variable}"
        result = climate_concordance_test(
            h1_edges,
            null_persistence,
            null_edge_ids,
            climate_edges,
            predictor=predictor,
            subset="global",
            minimum_supported_edges=minimum_supported_edges,
        )
        rows.append(
            {
                "variable": variable,
                "predictor": predictor,
                "weighted_r": result.statistic,
                "p_upper_raw": result.p_upper,
                "supported_edges": result.supported_edges,
            }
        )
    table = pd.DataFrame(rows)
    table["p_upper_holm"] = holm_adjust(table["p_upper_raw"].to_numpy(dtype=float))
    return table


__all__ = [
    "BIO_VARIABLES",
    "ClimateConcordanceResult",
    "aggregate_climate_to_h1_grid",
    "build_edge_climate_contrasts",
    "climate_concordance_test",
    "climate_driver_decomposition",
    "holm_adjust",
    "weighted_pearson",
]
