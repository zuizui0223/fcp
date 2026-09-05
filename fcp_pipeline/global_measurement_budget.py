"""Outcome-blind species budgeting for global RGFCA flower-colour measurement.

The global metadata/candidate pool can contain thousands of species.  RGFCA is
explicitly a bounded repeated-world-map design, so pixel measurement is capped at a
prospectively fixed number of species before any flower-colour pixels are opened.
"""
from __future__ import annotations

import hashlib
from typing import Iterable

import numpy as np
import pandas as pd


def _rank_key(taxon_id: int, *, seed: int) -> tuple[bytes, int]:
    payload = f"{int(seed)}\x1f{int(taxon_id)}".encode("utf-8")
    return hashlib.sha256(payload).digest(), int(taxon_id)


def select_measurement_taxa(
    taxon_ids: Iterable[int],
    *,
    maximum_species: int = 500,
    seed: int = 20260918,
) -> tuple[int, ...]:
    """Select at most ``maximum_species`` taxa uniformly by deterministic hash rank."""
    maximum = int(maximum_species)
    if maximum < 1:
        raise ValueError("maximum_species must be positive")
    unique = sorted({int(x) for x in taxon_ids})
    if not unique:
        raise ValueError("taxon_ids is empty")
    ranked = sorted(unique, key=lambda x: _rank_key(x, seed=int(seed)))
    return tuple(ranked[: min(maximum, len(ranked))])


def select_measurement_rows(
    candidate: pd.DataFrame,
    *,
    target_photos_per_species: int,
    maximum_species: int = 500,
    seed: int = 20260918,
) -> pd.DataFrame:
    """Return every frozen candidate row for the deterministically selected taxa."""
    required = {"inat_taxon_id", "photo_id"}
    missing = sorted(required - set(candidate.columns))
    if missing:
        raise ValueError(f"candidate frame lacks columns: {missing}")
    frame = candidate.copy()
    frame["inat_taxon_id"] = frame["inat_taxon_id"].astype(int)
    target = int(target_photos_per_species)
    if target < 1:
        raise ValueError("target_photos_per_species must be positive")
    counts = frame.groupby("inat_taxon_id", observed=True).size().astype(int)
    if len(counts) == 0 or not (counts == target).all():
        raise ValueError("candidate frame must contain exactly target rows for every taxon")
    if frame["photo_id"].duplicated().any():
        raise ValueError("candidate photo IDs must be unique")
    selected_taxa = set(
        select_measurement_taxa(counts.index, maximum_species=int(maximum_species), seed=int(seed))
    )
    selected = frame.loc[frame["inat_taxon_id"].isin(selected_taxa)].copy()
    if selected["inat_taxon_id"].nunique() != len(selected_taxa):
        raise RuntimeError("measurement species selection lost a taxon")
    if len(selected) != len(selected_taxa) * target:
        raise RuntimeError("measurement row denominator differs from selected taxa * target")
    sort_cols = [c for c in ("inat_taxon_id", "photo_id") if c in selected.columns]
    return selected.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)


__all__ = ["select_measurement_rows", "select_measurement_taxa"]
