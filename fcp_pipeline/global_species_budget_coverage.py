"""Equal-area coverage diagnostics for bounded RGFCA species budgets."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CellCoverageSummary:
    taxon_count: int
    occupied_cells: int
    total_taxon_cell_links: int
    cell_species_gini: float
    cell_species_counts: np.ndarray


def gini_nonnegative(values: Iterable[float]) -> float:
    x = np.asarray(list(values), dtype=float)
    if x.ndim != 1 or len(x) == 0:
        return float("nan")
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("gini values must be finite and nonnegative")
    if float(x.sum()) <= 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    return float((2.0 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum())) - (n + 1.0) / n)


def unique_taxon_cell_links(*frames: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for frame in frames:
        missing = {"inat_taxon_id", "cell_id"} - set(frame.columns)
        if missing:
            raise ValueError(f"discovery frame lacks columns: {sorted(missing)}")
        part = frame[["inat_taxon_id", "cell_id"]].dropna().copy()
        part["inat_taxon_id"] = part["inat_taxon_id"].astype(int)
        part["cell_id"] = part["cell_id"].astype(int)
        parts.append(part)
    if not parts:
        return pd.DataFrame(columns=["inat_taxon_id", "cell_id"])
    return pd.concat(parts, ignore_index=True).drop_duplicates(["inat_taxon_id", "cell_id"]).reset_index(drop=True)


def summarize_cell_coverage(
    links: pd.DataFrame,
    taxon_ids: Iterable[int],
    *,
    n_cells: int = 162,
) -> CellCoverageSummary:
    taxa = {int(x) for x in taxon_ids}
    if int(n_cells) < 1:
        raise ValueError("n_cells must be positive")
    if not taxa:
        return CellCoverageSummary(0, 0, 0, 0.0, np.zeros(int(n_cells), dtype=int))
    subset = links.loc[links["inat_taxon_id"].astype(int).isin(taxa), ["inat_taxon_id", "cell_id"]].drop_duplicates()
    if len(subset) and ((subset["cell_id"].astype(int) < 0) | (subset["cell_id"].astype(int) >= int(n_cells))).any():
        raise ValueError("cell_id outside frozen grid")
    counts = np.zeros(int(n_cells), dtype=int)
    if len(subset):
        vc = subset.groupby("cell_id", observed=True)["inat_taxon_id"].nunique()
        for cell, value in vc.items():
            counts[int(cell)] = int(value)
    return CellCoverageSummary(
        taxon_count=len(taxa),
        occupied_cells=int(np.count_nonzero(counts)),
        total_taxon_cell_links=int(len(subset)),
        cell_species_gini=gini_nonnegative(counts),
        cell_species_counts=counts,
    )


def occupied_cell_retention(reference: CellCoverageSummary, subset: CellCoverageSummary) -> float:
    if reference.occupied_cells == 0:
        return float("nan")
    ref = reference.cell_species_counts > 0
    sub = subset.cell_species_counts > 0
    return float(np.count_nonzero(ref & sub) / np.count_nonzero(ref))


def cell_count_correlation(reference: CellCoverageSummary, subset: CellCoverageSummary) -> float:
    x = np.asarray(reference.cell_species_counts, dtype=float)
    y = np.asarray(subset.cell_species_counts, dtype=float)
    keep = (x > 0) | (y > 0)
    if np.count_nonzero(keep) < 3:
        return float("nan")
    x = x[keep]
    y = y[keep]
    if np.ptp(x) <= 0 or np.ptp(y) <= 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


__all__ = [
    "CellCoverageSummary",
    "cell_count_correlation",
    "gini_nonnegative",
    "occupied_cell_retention",
    "summarize_cell_coverage",
    "unique_taxon_cell_links",
]
