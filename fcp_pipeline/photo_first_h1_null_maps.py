"""Persist matched edge-persistence maps from the cached H1 permutation run.

H2 climate concordance needs the exact edge persistence map for every H1 null
permutation, not only the scalar concentration statistic. This wrapper reuses the
cached H1 plan and stores those maps during the same permutation pass, avoiding a
second 999-permutation execution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

from .photo_first_atlas import PersistenceResult, _require_positive_int
from .photo_first_atlas_v2 import BIOLOGICAL_MORPH_LEVELS, STRUCTURAL_MISSING_LEVELS
from .photo_first_h1_fast import (
    evaluate_label_codes,
    permute_label_codes_within_species,
    prepare_h1_plan,
)
from .shared_transition_surface import EqualAreaGrid


@dataclass(frozen=True)
class CachedH1NullMaps:
    observed: PersistenceResult
    null_concentrations: np.ndarray
    null_persistence: np.ndarray
    edge_ids: tuple[str, ...]
    p_upper: float


def persistence_null_maps_cached(
    photos: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
    target_n: int,
    n_replicates: int,
    species_cap_per_cell: int,
    min_photos_per_cell: int,
    transition_quantile: float = 0.90,
    n_permutations: int = 999,
    sampling_seed: int = 20260903,
    permutation_seed: int = 20260904,
    species_col: str = "species",
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
    morph_col: str = "morph",
    morph_levels: Sequence[str] = BIOLOGICAL_MORPH_LEVELS,
    structural_missing_levels: Sequence[str] = STRUCTURAL_MISSING_LEVELS,
) -> CachedH1NullMaps:
    """Run H1 once and retain the matched edge map for every null permutation."""

    n_permutations = _require_positive_int("n_permutations", n_permutations)
    plan = prepare_h1_plan(
        photos,
        grid=grid,
        target_n=target_n,
        n_replicates=n_replicates,
        species_cap_per_cell=species_cap_per_cell,
        min_photos_per_cell=min_photos_per_cell,
        transition_quantile=transition_quantile,
        random_seed=sampling_seed,
        species_col=species_col,
        latitude_col=latitude_col,
        longitude_col=longitude_col,
        morph_col=morph_col,
        morph_levels=morph_levels,
        structural_missing_levels=structural_missing_levels,
    )
    observed = evaluate_label_codes(plan, plan.observed_label_codes)
    if not np.isfinite(observed.concentration):
        raise ValueError("observed persistence concentration is not estimable")
    edge_ids = tuple(observed.edge_table["edge_id"].astype(str))
    n_edges = len(edge_ids)
    null_concentrations = np.empty(n_permutations, dtype=float)
    null_persistence = np.empty((n_permutations, n_edges), dtype=float)

    seed_sequence = np.random.SeedSequence(int(permutation_seed))
    for position, child_seed in enumerate(seed_sequence.spawn(n_permutations)):
        permuted_codes = permute_label_codes_within_species(
            plan, rng=np.random.default_rng(child_seed)
        )
        result = evaluate_label_codes(plan, permuted_codes)
        if tuple(result.edge_table["edge_id"].astype(str)) != edge_ids:
            raise RuntimeError("H1 null edge order drifted from the observed graph")
        null_concentrations[position] = result.concentration
        null_persistence[position] = result.edge_table["persistence"].to_numpy(dtype=float)

    if not np.isfinite(null_concentrations).all():
        raise ValueError("one or more H1 null concentrations are not estimable")
    # Structural missingness can make unsupported edges NaN in every map. H2
    # selects edges using the observed opportunity denominator before indexing
    # this matrix; supported columns must therefore be finite in all null maps.
    supported = observed.edge_table["opportunities"].to_numpy(dtype=int) > 0
    if not np.isfinite(null_persistence[:, supported]).all():
        raise ValueError("one or more supported H1 null edge persistences are not estimable")
    p_upper = float(
        (1 + np.count_nonzero(null_concentrations >= observed.concentration))
        / (n_permutations + 1)
    )
    return CachedH1NullMaps(
        observed=observed,
        null_concentrations=null_concentrations,
        null_persistence=null_persistence,
        edge_ids=edge_ids,
        p_upper=p_upper,
    )


__all__ = ["CachedH1NullMaps", "persistence_null_maps_cached"]
