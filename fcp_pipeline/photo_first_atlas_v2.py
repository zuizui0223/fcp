"""Quality-safe H1 implementation for the prospective random photo-first atlas.

This module supersedes the pre-data v1 implementation for H1 inference. The key
change is that ``mixed_uncertain`` is treated as structural measurement
missingness, never as a fifth biological flower-colour state. Its geographic
positions are held fixed under the null while only classifiable biological morph
labels are permuted within species.

The photograph remains the sampling unit. Sampling is still colour-blind and is
performed on the complete measured photo table before structural-missing rows are
removed from cell composition summaries.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .photo_first_atlas import (
    PersistenceResult,
    _persistence_concentration,
    _require_positive_int,
    adjacent_grid_edges,
    cell_first_species_capped_sample,
    coarse_morph_from_palette,
    jensen_shannon_divergence,
    prepare_photo_grid,
    species_capped_sampling_capacity,
    validate_photo_table,
)
from .shared_transition_surface import EqualAreaGrid


BIOLOGICAL_MORPH_LEVELS = (
    "white",
    "yellow_orange",
    "red_pink",
    "blue_purple",
)
STRUCTURAL_MISSING_LEVELS = ("mixed_uncertain",)


def _validate_morph_roles(
    photos: pd.DataFrame,
    *,
    morph_col: str,
    morph_levels: Sequence[str],
    structural_missing_levels: Sequence[str],
) -> tuple[tuple[str, ...], frozenset[str]]:
    levels = tuple(str(value) for value in morph_levels)
    missing_levels = frozenset(str(value) for value in structural_missing_levels)
    if len(levels) < 2 or len(set(levels)) != len(levels):
        raise ValueError("morph_levels must contain at least two unique biological states")
    if set(levels).intersection(missing_levels):
        raise ValueError("biological and structural-missing morph levels must be disjoint")
    observed = set(photos[morph_col].astype(str))
    allowed = set(levels).union(missing_levels)
    unknown = sorted(observed.difference(allowed))
    if unknown:
        raise ValueError(f"unrecognized morph labels: {unknown}")
    return levels, missing_levels


def _cell_compositions_quality_safe(
    sampled: pd.DataFrame,
    *,
    morph_levels: tuple[str, ...],
    structural_missing_levels: frozenset[str],
    min_classifiable_photos_per_cell: int,
    morph_col: str,
) -> tuple[dict[int, np.ndarray], dict[int, int], dict[int, int]]:
    """Return biological morph compositions after excluding uncertain measurements."""

    compositions: dict[int, np.ndarray] = {}
    classifiable_n: dict[int, int] = {}
    raw_n: dict[int, int] = {}
    for cell_id, rows in sampled.groupby("cell_id", sort=True):
        raw_n[int(cell_id)] = int(len(rows))
        keep = ~rows[morph_col].astype(str).isin(structural_missing_levels)
        classifiable = rows.loc[keep]
        n = int(len(classifiable))
        classifiable_n[int(cell_id)] = n
        if n < min_classifiable_photos_per_cell:
            continue
        counts = classifiable[morph_col].astype(str).value_counts().reindex(
            morph_levels,
            fill_value=0,
        )
        vector = counts.to_numpy(dtype=float)
        if vector.sum() <= 0:
            continue
        compositions[int(cell_id)] = vector / vector.sum()
    return compositions, classifiable_n, raw_n


def replicate_edge_table(
    sampled: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
    morph_levels: Sequence[str] = BIOLOGICAL_MORPH_LEVELS,
    structural_missing_levels: Sequence[str] = STRUCTURAL_MISSING_LEVELS,
    min_photos_per_cell: int,
    transition_quantile: float,
    rng: np.random.Generator | None = None,
    morph_col: str = "morph",
) -> pd.DataFrame:
    """Build one quality-safe edge table from a fixed colour-blind photo sample."""

    min_photos_per_cell = _require_positive_int("min_photos_per_cell", min_photos_per_cell)
    transition_quantile = float(transition_quantile)
    if not 0.0 < transition_quantile < 1.0:
        raise ValueError("transition_quantile must lie strictly inside (0, 1)")
    if rng is None:
        rng = np.random.default_rng(0)

    levels = tuple(str(value) for value in morph_levels)
    missing = frozenset(str(value) for value in structural_missing_levels)
    compositions, classifiable_n, raw_n = _cell_compositions_quality_safe(
        sampled,
        morph_levels=levels,
        structural_missing_levels=missing,
        min_classifiable_photos_per_cell=min_photos_per_cell,
        morph_col=morph_col,
    )

    rows: list[dict[str, object]] = []
    for left, right in adjacent_grid_edges(grid):
        left_i = int(left)
        right_i = int(right)
        evaluable = left_i in compositions and right_i in compositions
        intensity = (
            jensen_shannon_divergence(compositions[left_i], compositions[right_i])
            if evaluable
            else np.nan
        )
        rows.append(
            {
                "edge_id": f"{left_i}:{right_i}",
                "cell_i": left_i,
                "cell_j": right_i,
                "raw_n_i": int(raw_n.get(left_i, 0)),
                "raw_n_j": int(raw_n.get(right_i, 0)),
                "classifiable_n_i": int(classifiable_n.get(left_i, 0)),
                "classifiable_n_j": int(classifiable_n.get(right_i, 0)),
                "evaluable": bool(evaluable),
                "transition_intensity": float(intensity) if evaluable else np.nan,
                "is_transition": False,
            }
        )
    table = pd.DataFrame(rows)
    evaluable_index = table.index[table["evaluable"]].to_numpy(dtype=int)
    if len(evaluable_index) == 0:
        return table

    intensities = table.loc[evaluable_index, "transition_intensity"].to_numpy(dtype=float)
    tie_break = rng.random(len(evaluable_index))
    order = np.lexsort((tie_break, -intensities))
    n_transition = max(
        1,
        int(np.ceil((1.0 - transition_quantile) * len(evaluable_index))),
    )
    chosen = evaluable_index[order[:n_transition]]
    table.loc[chosen, "is_transition"] = True
    return table


def run_boundary_persistence(
    photos: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
    target_n: int,
    n_replicates: int,
    species_cap_per_cell: int,
    min_photos_per_cell: int,
    transition_quantile: float = 0.90,
    random_seed: int = 20260903,
    species_col: str = "species",
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
    morph_col: str = "morph",
    morph_levels: Sequence[str] = BIOLOGICAL_MORPH_LEVELS,
    structural_missing_levels: Sequence[str] = STRUCTURAL_MISSING_LEVELS,
) -> PersistenceResult:
    """Estimate recurrent transition persistence with uncertain measurements excluded."""

    validate_photo_table(
        photos,
        species_col=species_col,
        latitude_col=latitude_col,
        longitude_col=longitude_col,
        morph_col=morph_col,
    )
    target_n = _require_positive_int("target_n", target_n)
    n_replicates = _require_positive_int("n_replicates", n_replicates)
    species_cap_per_cell = _require_positive_int(
        "species_cap_per_cell", species_cap_per_cell
    )
    min_photos_per_cell = _require_positive_int("min_photos_per_cell", min_photos_per_cell)
    levels, missing = _validate_morph_roles(
        photos,
        morph_col=morph_col,
        morph_levels=morph_levels,
        structural_missing_levels=structural_missing_levels,
    )

    work = prepare_photo_grid(
        photos,
        grid=grid,
        latitude_col=latitude_col,
        longitude_col=longitude_col,
    )
    capacity = species_capped_sampling_capacity(
        work,
        species_cap_per_cell=species_cap_per_cell,
        species_col=species_col,
    )
    if capacity < target_n:
        raise ValueError(
            "not_evaluable_fixed_replicate_size: "
            f"species-capped capacity {capacity} is below target_n {target_n}"
        )

    edge_geometry = adjacent_grid_edges(grid)
    edge_ids = [f"{int(left)}:{int(right)}" for left, right in edge_geometry]
    opportunities = {edge_id: 0 for edge_id in edge_ids}
    transitions = {edge_id: 0 for edge_id in edge_ids}
    sample_sizes: list[int] = []

    seed_sequence = np.random.SeedSequence(int(random_seed))
    for child_seed in seed_sequence.spawn(n_replicates):
        rng = np.random.default_rng(child_seed)
        sampled = cell_first_species_capped_sample(
            work,
            target_n=target_n,
            species_cap_per_cell=species_cap_per_cell,
            rng=rng,
            species_col=species_col,
        )
        if len(sampled) != target_n:
            raise RuntimeError(
                "not_evaluable_fixed_replicate_size: "
                f"replicate sampled {len(sampled)} photos instead of {target_n}"
            )
        sample_sizes.append(int(len(sampled)))
        edge_table = replicate_edge_table(
            sampled,
            grid=grid,
            morph_levels=levels,
            structural_missing_levels=tuple(sorted(missing)),
            min_photos_per_cell=min_photos_per_cell,
            transition_quantile=transition_quantile,
            rng=rng,
            morph_col=morph_col,
        )
        for row in edge_table.itertuples(index=False):
            if bool(row.evaluable):
                opportunities[row.edge_id] += 1
                if bool(row.is_transition):
                    transitions[row.edge_id] += 1

    rows = []
    for left, right in edge_geometry:
        edge_id = f"{int(left)}:{int(right)}"
        opportunity = int(opportunities[edge_id])
        transition_count = int(transitions[edge_id])
        persistence = transition_count / opportunity if opportunity > 0 else np.nan
        rows.append(
            {
                "edge_id": edge_id,
                "cell_i": int(left),
                "cell_j": int(right),
                "opportunities": opportunity,
                "transition_count": transition_count,
                "persistence": float(persistence) if opportunity > 0 else np.nan,
            }
        )
    persistence_table = pd.DataFrame(rows)
    concentration, transition_rate = _persistence_concentration(persistence_table)
    return PersistenceResult(
        edge_table=persistence_table,
        concentration=concentration,
        transition_rate=transition_rate,
        mean_sampled_photos=float(np.mean(sample_sizes)),
        morph_levels=levels,
        n_replicates=n_replicates,
    )


def species_conditioned_morph_permutation(
    photos: pd.DataFrame,
    *,
    rng: np.random.Generator,
    species_col: str = "species",
    morph_col: str = "morph",
    structural_missing_levels: Sequence[str] = STRUCTURAL_MISSING_LEVELS,
) -> pd.DataFrame:
    """Shuffle only classifiable morphs within species; keep uncertainty fixed."""

    out = photos.reset_index(drop=True).copy()
    values = out[morph_col].astype(object).to_numpy(copy=True)
    species = out[species_col].astype(str).to_numpy()
    missing = frozenset(str(value) for value in structural_missing_levels)
    is_missing = np.asarray([str(value) in missing for value in values], dtype=bool)
    for species_name in np.unique(species):
        idx = np.flatnonzero((species == species_name) & ~is_missing)
        if len(idx) > 1:
            values[idx] = values[idx][rng.permutation(len(idx))]
    out[morph_col] = values
    return out


def persistence_null_test(
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
) -> tuple[PersistenceResult, np.ndarray, float]:
    """Species-conditioned null with the structural-missing mask held fixed."""

    n_permutations = _require_positive_int("n_permutations", n_permutations)
    observed = run_boundary_persistence(
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
    if not np.isfinite(observed.concentration):
        raise ValueError("observed persistence concentration is not estimable")

    null = np.empty(n_permutations, dtype=float)
    seed_sequence = np.random.SeedSequence(int(permutation_seed))
    for position, child_seed in enumerate(seed_sequence.spawn(n_permutations)):
        permuted = species_conditioned_morph_permutation(
            photos,
            rng=np.random.default_rng(child_seed),
            species_col=species_col,
            morph_col=morph_col,
            structural_missing_levels=structural_missing_levels,
        )
        result = run_boundary_persistence(
            permuted,
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
            morph_levels=observed.morph_levels,
            structural_missing_levels=structural_missing_levels,
        )
        null[position] = result.concentration
    if not np.isfinite(null).all():
        raise ValueError("one or more null persistence concentrations are not estimable")
    p_upper = float(
        (1 + np.count_nonzero(null >= observed.concentration)) / (n_permutations + 1)
    )
    return observed, null, p_upper


__all__ = [
    "BIOLOGICAL_MORPH_LEVELS",
    "STRUCTURAL_MISSING_LEVELS",
    "PersistenceResult",
    "adjacent_grid_edges",
    "cell_first_species_capped_sample",
    "coarse_morph_from_palette",
    "jensen_shannon_divergence",
    "persistence_null_test",
    "prepare_photo_grid",
    "replicate_edge_table",
    "run_boundary_persistence",
    "species_capped_sampling_capacity",
    "species_conditioned_morph_permutation",
    "validate_photo_table",
]
