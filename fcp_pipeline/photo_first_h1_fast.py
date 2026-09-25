"""Cached execution engine for the prospective random photo-first H1 test.

The frozen H1 design samples photographs without using colour. Under the primary
null, species identity, locations, sampling seeds, and the structural-missing
mask are all fixed; only classifiable morph labels are permuted within species.
This module exploits those invariants without changing the statistic: the 200
replicate samples, edge opportunities, and tie-break random numbers are prepared
once, then reused for observed and null evaluations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

from .photo_first_atlas import (
    PersistenceResult,
    _persistence_concentration,
    _require_positive_int,
    adjacent_grid_edges,
    cell_first_species_capped_sample,
    prepare_photo_grid,
    species_capped_sampling_capacity,
    validate_photo_table,
)
from .photo_first_atlas_v2 import (
    BIOLOGICAL_MORPH_LEVELS,
    STRUCTURAL_MISSING_LEVELS,
    _validate_morph_roles,
)
from .shared_transition_surface import EqualAreaGrid


@dataclass(frozen=True)
class PreparedH1Plan:
    grid: EqualAreaGrid
    target_n: int
    n_replicates: int
    min_photos_per_cell: int
    transition_quantile: float
    morph_levels: tuple[str, ...]
    structural_missing_levels: tuple[str, ...]
    edge_geometry: np.ndarray
    sampled_indices: np.ndarray
    classifiable_sample_indices: np.ndarray
    classifiable_group_base: np.ndarray
    evaluable_edges: np.ndarray
    tie_breaks: np.ndarray
    n_transitions: np.ndarray
    opportunities: np.ndarray
    observed_label_codes: np.ndarray
    species_codes: np.ndarray
    classifiable_global_indices: np.ndarray
    species_target_order: np.ndarray

    @property
    def n_cells(self) -> int:
        return int(self.grid.n_cells)

    @property
    def n_morphs(self) -> int:
        return int(len(self.morph_levels))


def _encode_morphs(
    values: pd.Series,
    *,
    morph_levels: tuple[str, ...],
    structural_missing_levels: frozenset[str],
) -> np.ndarray:
    labels = values.astype(str).to_numpy()
    codes = np.full(len(labels), -1, dtype=np.int16)
    for code, level in enumerate(morph_levels):
        codes[labels == level] = code
    unknown = (codes < 0) & ~np.isin(labels, list(structural_missing_levels))
    if np.any(unknown):
        raise ValueError("morph encoding encountered a non-frozen label")
    return codes


def prepare_h1_plan(
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
) -> PreparedH1Plan:
    """Freeze all colour-independent replicate work for observed and null runs."""

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
    min_photos_per_cell = _require_positive_int(
        "min_photos_per_cell", min_photos_per_cell
    )
    transition_quantile = float(transition_quantile)
    if not 0.0 < transition_quantile < 1.0:
        raise ValueError("transition_quantile must lie strictly inside (0, 1)")

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
    ).reset_index(drop=True)
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

    work = work.copy()
    row_id_col = "__photo_first_fast_row_id"
    if row_id_col in work.columns:
        raise ValueError(f"reserved column already exists: {row_id_col}")
    work[row_id_col] = np.arange(len(work), dtype=np.int64)
    label_codes = _encode_morphs(
        work[morph_col], morph_levels=levels, structural_missing_levels=missing
    )
    species_codes, _ = pd.factorize(work[species_col].astype(str), sort=True)
    species_codes = np.asarray(species_codes, dtype=np.int32)
    classifiable_global = np.flatnonzero(label_codes >= 0).astype(np.int64)
    species_target_order = np.argsort(
        species_codes[classifiable_global], kind="stable"
    ).astype(np.int64)

    edges = adjacent_grid_edges(grid).astype(np.int32)
    n_edges = len(edges)
    sampled_indices = np.empty((n_replicates, target_n), dtype=np.int64)
    evaluable_edges = np.zeros((n_replicates, n_edges), dtype=bool)
    tie_breaks = np.full((n_replicates, n_edges), np.nan, dtype=float)
    n_transitions = np.zeros(n_replicates, dtype=np.int32)
    classifiable_rows: list[np.ndarray] = []
    group_bases: list[np.ndarray] = []
    cell_ids = work["cell_id"].to_numpy(dtype=np.int32)
    left = edges[:, 0]
    right = edges[:, 1]

    seed_sequence = np.random.SeedSequence(int(random_seed))
    for replicate_index, child_seed in enumerate(seed_sequence.spawn(n_replicates)):
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
        indices = sampled[row_id_col].to_numpy(dtype=np.int64)
        sampled_indices[replicate_index] = indices
        keep = label_codes[indices] >= 0
        class_indices = indices[keep]
        class_cells = cell_ids[class_indices]
        classifiable_rows.append(class_indices)
        group_bases.append(
            ((replicate_index * grid.n_cells + class_cells) * len(levels)).astype(
                np.int64
            )
        )
        classifiable_counts = np.bincount(
            class_cells, minlength=grid.n_cells
        )
        cell_ok = classifiable_counts >= min_photos_per_cell
        edge_ok = cell_ok[left] & cell_ok[right]
        evaluable_edges[replicate_index] = edge_ok
        n_evaluable = int(np.count_nonzero(edge_ok))
        if n_evaluable:
            tie_breaks[replicate_index, edge_ok] = rng.random(n_evaluable)
            n_transitions[replicate_index] = max(
                1,
                int(np.ceil((1.0 - transition_quantile) * n_evaluable)),
            )

    classifiable_sample_indices = (
        np.concatenate(classifiable_rows).astype(np.int64, copy=False)
        if classifiable_rows
        else np.empty(0, dtype=np.int64)
    )
    classifiable_group_base = (
        np.concatenate(group_bases).astype(np.int64, copy=False)
        if group_bases
        else np.empty(0, dtype=np.int64)
    )
    opportunities = evaluable_edges.sum(axis=0).astype(np.int32)
    return PreparedH1Plan(
        grid=grid,
        target_n=target_n,
        n_replicates=n_replicates,
        min_photos_per_cell=min_photos_per_cell,
        transition_quantile=transition_quantile,
        morph_levels=levels,
        structural_missing_levels=tuple(sorted(missing)),
        edge_geometry=edges,
        sampled_indices=sampled_indices,
        classifiable_sample_indices=classifiable_sample_indices,
        classifiable_group_base=classifiable_group_base,
        evaluable_edges=evaluable_edges,
        tie_breaks=tie_breaks,
        n_transitions=n_transitions,
        opportunities=opportunities,
        observed_label_codes=label_codes,
        species_codes=species_codes,
        classifiable_global_indices=classifiable_global,
        species_target_order=species_target_order,
    )


def _jsd_rows(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    midpoint = 0.5 * (p + q)

    def kl(a: np.ndarray) -> np.ndarray:
        ratio = np.ones_like(a, dtype=float)
        positive = a > 0.0
        np.divide(a, midpoint, out=ratio, where=positive)
        terms = np.zeros_like(a, dtype=float)
        terms[positive] = a[positive] * np.log2(ratio[positive])
        return terms.sum(axis=-1)

    return np.clip(0.5 * kl(p) + 0.5 * kl(q), 0.0, 1.0)


def evaluate_label_codes(
    plan: PreparedH1Plan,
    label_codes: np.ndarray,
) -> PersistenceResult:
    """Evaluate one observed or permuted label vector on a prepared H1 plan."""

    codes = np.asarray(label_codes, dtype=np.int16)
    if codes.shape != plan.observed_label_codes.shape:
        raise ValueError("label code vector does not match the prepared photo denominator")
    selected_codes = codes[plan.classifiable_sample_indices]
    if np.any(selected_codes < 0) or np.any(selected_codes >= plan.n_morphs):
        raise ValueError("prepared classifiable sample contains invalid morph codes")

    flat_size = plan.n_replicates * plan.n_cells * plan.n_morphs
    counts = np.bincount(
        plan.classifiable_group_base + selected_codes.astype(np.int64),
        minlength=flat_size,
    ).reshape(plan.n_replicates, plan.n_cells, plan.n_morphs)
    totals = counts.sum(axis=2)
    compositions = np.zeros_like(counts, dtype=float)
    np.divide(
        counts,
        totals[:, :, None],
        out=compositions,
        where=totals[:, :, None] > 0,
    )

    left = plan.edge_geometry[:, 0]
    right = plan.edge_geometry[:, 1]
    intensities = _jsd_rows(compositions[:, left, :], compositions[:, right, :])
    transition_counts = np.zeros(len(plan.edge_geometry), dtype=np.int32)
    for replicate_index in range(plan.n_replicates):
        eligible = np.flatnonzero(plan.evaluable_edges[replicate_index])
        n_take = int(plan.n_transitions[replicate_index])
        if not len(eligible) or n_take == 0:
            continue
        order = np.lexsort(
            (
                plan.tie_breaks[replicate_index, eligible],
                -intensities[replicate_index, eligible],
            )
        )
        chosen = eligible[order[:n_take]]
        np.add.at(transition_counts, chosen, 1)

    rows: list[dict[str, object]] = []
    for edge_index, (left_cell, right_cell) in enumerate(plan.edge_geometry):
        opportunity = int(plan.opportunities[edge_index])
        transition_count = int(transition_counts[edge_index])
        persistence = transition_count / opportunity if opportunity else np.nan
        rows.append(
            {
                "edge_id": f"{int(left_cell)}:{int(right_cell)}",
                "cell_i": int(left_cell),
                "cell_j": int(right_cell),
                "opportunities": opportunity,
                "transition_count": transition_count,
                "persistence": float(persistence) if opportunity else np.nan,
            }
        )
    edge_table = pd.DataFrame(rows)
    concentration, transition_rate = _persistence_concentration(edge_table)
    return PersistenceResult(
        edge_table=edge_table,
        concentration=concentration,
        transition_rate=transition_rate,
        mean_sampled_photos=float(plan.target_n),
        morph_levels=plan.morph_levels,
        n_replicates=plan.n_replicates,
    )


def permute_label_codes_within_species(
    plan: PreparedH1Plan,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Uniformly permute classifiable labels within species; keep missing fixed."""

    out = plan.observed_label_codes.copy()
    indices = plan.classifiable_global_indices
    if len(indices) <= 1:
        return out
    species = plan.species_codes[indices]
    random_keys = rng.random(len(indices))
    source_order = np.lexsort((random_keys, species))
    target_order = plan.species_target_order
    out[indices[target_order]] = plan.observed_label_codes[indices[source_order]]
    return out


def persistence_null_test_cached(
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
    """Run the frozen null while caching all colour-independent replicate work."""

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

    null = np.empty(n_permutations, dtype=float)
    seed_sequence = np.random.SeedSequence(int(permutation_seed))
    for position, child_seed in enumerate(seed_sequence.spawn(n_permutations)):
        permuted_codes = permute_label_codes_within_species(
            plan, rng=np.random.default_rng(child_seed)
        )
        result = evaluate_label_codes(plan, permuted_codes)
        null[position] = result.concentration
    if not np.isfinite(null).all():
        raise ValueError("one or more null persistence concentrations are not estimable")
    p_upper = float(
        (1 + np.count_nonzero(null >= observed.concentration))
        / (n_permutations + 1)
    )
    return observed, null, p_upper


__all__ = [
    "PreparedH1Plan",
    "evaluate_label_codes",
    "permute_label_codes_within_species",
    "persistence_null_test_cached",
    "prepare_h1_plan",
]
