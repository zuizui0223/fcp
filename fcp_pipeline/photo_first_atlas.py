"""Photo-first random atlas helpers for recurrent flower-colour boundaries.

This module implements a new prospective analysis line.  It does not reuse the
terminal 60,000-photo experiment as a favourable subset and it does not change
the frozen six-species or 34-species analyses.

The inferential display is species-free, but species identity is retained for
sampling caps and the species-conditioned null.  Boundary persistence is always
computed with an opportunity denominator: an edge contributes to the denominator
only in replicates where both adjacent cells have sufficient sampled photos.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .shared_transition_surface import EqualAreaGrid, equal_area_cell_ids


@dataclass(frozen=True)
class PersistenceResult:
    edge_table: pd.DataFrame
    concentration: float
    transition_rate: float
    mean_sampled_photos: float
    morph_levels: tuple[str, ...]
    n_replicates: int


def _require_positive_int(name: str, value: int) -> int:
    value = int(value)
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def validate_photo_table(
    photos: pd.DataFrame,
    *,
    species_col: str = "species",
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
    morph_col: str = "morph",
) -> None:
    required = {species_col, latitude_col, longitude_col, morph_col}
    missing = sorted(required.difference(photos.columns))
    if missing:
        raise ValueError(f"photo table is missing required columns: {missing}")
    if len(photos) == 0:
        raise ValueError("photo table must contain at least one row")

    lat = pd.to_numeric(photos[latitude_col], errors="coerce").to_numpy(dtype=float)
    lon = pd.to_numeric(photos[longitude_col], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("coordinates must be finite")
    if np.any((lat < -90.0) | (lat > 90.0)):
        raise ValueError("latitude must lie in [-90, 90]")
    if np.any((lon < -180.0) | (lon > 180.0)):
        raise ValueError("longitude must lie in [-180, 180]")

    species = photos[species_col].astype(str).str.strip()
    morph = photos[morph_col].astype(str).str.strip()
    if (species == "").any() or (morph == "").any():
        raise ValueError("species and morph labels must be non-empty")


def coarse_morph_from_palette(
    fractions: Mapping[str, float],
    *,
    minimum_dominant_fraction: float = 0.50,
    minimum_margin: float = 0.10,
) -> str:
    """Collapse the existing nine flower-palette fractions into five coarse states.

    The four chromatic/achromatic groups are fixed independently of species.
    Ambiguous photographs are retained as ``mixed_uncertain`` rather than being
    forced to a colour state.
    """

    minimum_dominant_fraction = float(minimum_dominant_fraction)
    minimum_margin = float(minimum_margin)
    if not 0.0 <= minimum_dominant_fraction <= 1.0:
        raise ValueError("minimum_dominant_fraction must lie in [0, 1]")
    if not 0.0 <= minimum_margin <= 1.0:
        raise ValueError("minimum_margin must lie in [0, 1]")

    expected = {
        "blue", "bronze", "magenta", "orange", "pink",
        "purple", "red", "white", "yellow",
    }
    missing = expected.difference(fractions)
    if missing:
        raise ValueError(f"palette fractions missing: {sorted(missing)}")

    values = {key: float(fractions[key]) for key in expected}
    if any((not np.isfinite(value) or value < 0.0) for value in values.values()):
        raise ValueError("palette fractions must be finite and non-negative")
    total = float(sum(values.values()))
    if total <= 0.0:
        return "mixed_uncertain"
    values = {key: value / total for key, value in values.items()}

    grouped = {
        "white": values["white"],
        "yellow_orange": values["yellow"] + values["orange"] + values["bronze"],
        "red_pink": values["red"] + values["pink"] + values["magenta"],
        "blue_purple": values["blue"] + values["purple"],
    }
    ordered = sorted(grouped.items(), key=lambda item: (-item[1], item[0]))
    winner, top = ordered[0]
    second = ordered[1][1]
    if top < minimum_dominant_fraction or (top - second) < minimum_margin:
        return "mixed_uncertain"
    return winner


def prepare_photo_grid(
    photos: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
    latitude_col: str = "latitude",
    longitude_col: str = "longitude",
) -> pd.DataFrame:
    out = photos.reset_index(drop=True).copy()
    out["cell_id"] = equal_area_cell_ids(
        pd.to_numeric(out[latitude_col], errors="raise").to_numpy(dtype=float),
        pd.to_numeric(out[longitude_col], errors="raise").to_numpy(dtype=float),
        grid,
    )
    return out


def adjacent_grid_edges(grid: EqualAreaGrid) -> np.ndarray:
    """Return unique rook-neighbour edges, wrapping longitude at the dateline."""

    edges: set[tuple[int, int]] = set()
    for row in range(grid.n_sinlat):
        for col in range(grid.n_lon):
            cell = row * grid.n_lon + col
            right = row * grid.n_lon + ((col + 1) % grid.n_lon)
            edges.add(tuple(sorted((cell, right))))
            if row + 1 < grid.n_sinlat:
                above = (row + 1) * grid.n_lon + col
                edges.add(tuple(sorted((cell, above))))
    return np.asarray(sorted(edges), dtype=int)


def cell_first_species_capped_sample(
    photos_with_cells: pd.DataFrame,
    *,
    target_n: int,
    species_cap_per_cell: int,
    rng: np.random.Generator,
    species_col: str = "species",
) -> pd.DataFrame:
    """Sample broadly over cells while limiting any species within each cell.

    First, each cell-by-species pool is randomly reduced to the fixed cap.  Then
    the remaining candidates are drawn round-robin over cells.  Colour/morph is
    never used to decide which row is sampled.
    """

    target_n = _require_positive_int("target_n", target_n)
    species_cap_per_cell = _require_positive_int(
        "species_cap_per_cell", species_cap_per_cell
    )
    if "cell_id" not in photos_with_cells.columns:
        raise ValueError("photos_with_cells must contain cell_id")
    if species_col not in photos_with_cells.columns:
        raise ValueError(f"photos_with_cells must contain {species_col}")

    work = photos_with_cells.reset_index(drop=True)
    pools: dict[int, list[int]] = {}
    for cell_id, cell_rows in work.groupby("cell_id", sort=True):
        retained: list[int] = []
        for _, species_rows in cell_rows.groupby(species_col, sort=True):
            idx = species_rows.index.to_numpy(dtype=int).copy()
            rng.shuffle(idx)
            retained.extend(idx[:species_cap_per_cell].tolist())
        retained_array = np.asarray(retained, dtype=int)
        rng.shuffle(retained_array)
        pools[int(cell_id)] = retained_array.tolist()

    selected: list[int] = []
    cursors = {cell_id: 0 for cell_id in pools}
    active = np.asarray(sorted(pools), dtype=int)
    while len(selected) < target_n and len(active) > 0:
        order = active.copy()
        rng.shuffle(order)
        next_active: list[int] = []
        progressed = False
        for cell_id_raw in order:
            cell_id = int(cell_id_raw)
            cursor = cursors[cell_id]
            pool = pools[cell_id]
            if cursor < len(pool):
                selected.append(pool[cursor])
                cursors[cell_id] = cursor + 1
                progressed = True
                if cursors[cell_id] < len(pool):
                    next_active.append(cell_id)
                if len(selected) >= target_n:
                    break
        if not progressed:
            break
        if len(selected) < target_n:
            # Preserve cells not visited after an early break only when needed.
            active = np.asarray(sorted(set(next_active)), dtype=int)

    return work.iloc[selected].reset_index(drop=True)


def jensen_shannon_divergence(p: Sequence[float], q: Sequence[float]) -> float:
    """Jensen-Shannon divergence in bits, bounded by [0, 1]."""

    p_arr = np.asarray(p, dtype=float)
    q_arr = np.asarray(q, dtype=float)
    if p_arr.ndim != 1 or q_arr.ndim != 1 or p_arr.shape != q_arr.shape:
        raise ValueError("p and q must be equal-length one-dimensional vectors")
    if not np.isfinite(p_arr).all() or not np.isfinite(q_arr).all():
        raise ValueError("p and q must be finite")
    if np.any(p_arr < 0.0) or np.any(q_arr < 0.0):
        raise ValueError("p and q must be non-negative")
    if p_arr.sum() <= 0.0 or q_arr.sum() <= 0.0:
        raise ValueError("p and q must have positive mass")

    p_arr = p_arr / p_arr.sum()
    q_arr = q_arr / q_arr.sum()
    midpoint = 0.5 * (p_arr + q_arr)

    def kl_bits(a: np.ndarray, b: np.ndarray) -> float:
        keep = a > 0.0
        return float(np.sum(a[keep] * np.log2(a[keep] / b[keep])))

    value = 0.5 * kl_bits(p_arr, midpoint) + 0.5 * kl_bits(q_arr, midpoint)
    return float(np.clip(value, 0.0, 1.0))


def _cell_compositions(
    sampled: pd.DataFrame,
    *,
    morph_levels: tuple[str, ...],
    min_photos_per_cell: int,
    morph_col: str,
) -> tuple[dict[int, np.ndarray], dict[int, int]]:
    compositions: dict[int, np.ndarray] = {}
    cell_n: dict[int, int] = {}
    for cell_id, rows in sampled.groupby("cell_id", sort=True):
        n = int(len(rows))
        cell_n[int(cell_id)] = n
        if n < min_photos_per_cell:
            continue
        counts = rows[morph_col].value_counts().reindex(morph_levels, fill_value=0)
        vector = counts.to_numpy(dtype=float)
        compositions[int(cell_id)] = vector / vector.sum()
    return compositions, cell_n


def replicate_edge_table(
    sampled: pd.DataFrame,
    *,
    grid: EqualAreaGrid,
    morph_levels: tuple[str, ...],
    min_photos_per_cell: int,
    transition_quantile: float,
    morph_col: str = "morph",
) -> pd.DataFrame:
    min_photos_per_cell = _require_positive_int(
        "min_photos_per_cell", min_photos_per_cell
    )
    transition_quantile = float(transition_quantile)
    if not 0.0 < transition_quantile < 1.0:
        raise ValueError("transition_quantile must lie strictly inside (0, 1)")

    compositions, _ = _cell_compositions(
        sampled,
        morph_levels=morph_levels,
        min_photos_per_cell=min_photos_per_cell,
        morph_col=morph_col,
    )
    edges = adjacent_grid_edges(grid)
    rows: list[dict[str, object]] = []
    for left, right in edges:
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
                "evaluable": bool(evaluable),
                "transition_intensity": float(intensity) if evaluable else np.nan,
                "is_transition": False,
            }
        )
    table = pd.DataFrame(rows)
    evaluable_values = table.loc[table["evaluable"], "transition_intensity"].to_numpy(
        dtype=float
    )
    if len(evaluable_values) == 0:
        return table
    threshold = float(np.quantile(evaluable_values, transition_quantile))
    mask = table["evaluable"] & (table["transition_intensity"] >= threshold)
    table.loc[mask, "is_transition"] = True
    return table


def _persistence_concentration(edge_table: pd.DataFrame) -> tuple[float, float]:
    supported = edge_table[edge_table["opportunities"] > 0].copy()
    if len(supported) == 0:
        return float("nan"), float("nan")
    total_opportunity = float(supported["opportunities"].sum())
    total_transitions = float(supported["transition_count"].sum())
    transition_rate = total_transitions / total_opportunity
    concentration = float(
        np.sum(
            supported["opportunities"].to_numpy(dtype=float)
            * (supported["persistence"].to_numpy(dtype=float) - transition_rate) ** 2
        )
        / total_opportunity
    )
    return concentration, transition_rate


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
    morph_levels: Sequence[str] | None = None,
) -> PersistenceResult:
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

    if morph_levels is None:
        levels = tuple(sorted(photos[morph_col].astype(str).unique().tolist()))
    else:
        levels = tuple(str(value) for value in morph_levels)
    if len(levels) < 2:
        raise ValueError("at least two morph levels are required")
    observed_levels = set(photos[morph_col].astype(str))
    if not observed_levels.issubset(levels):
        raise ValueError("morph_levels does not cover all observed morph labels")

    work = prepare_photo_grid(
        photos,
        grid=grid,
        latitude_col=latitude_col,
        longitude_col=longitude_col,
    )
    edge_geometry = adjacent_grid_edges(grid)
    edge_ids = [f"{int(left)}:{int(right)}" for left, right in edge_geometry]
    opportunities = {edge_id: 0 for edge_id in edge_ids}
    transitions = {edge_id: 0 for edge_id in edge_ids}
    sample_sizes: list[int] = []

    seed_sequence = np.random.SeedSequence(int(random_seed))
    child_seeds = seed_sequence.spawn(n_replicates)
    for child_seed in child_seeds:
        rng = np.random.default_rng(child_seed)
        sampled = cell_first_species_capped_sample(
            work,
            target_n=target_n,
            species_cap_per_cell=species_cap_per_cell,
            rng=rng,
            species_col=species_col,
        )
        sample_sizes.append(int(len(sampled)))
        edge_table = replicate_edge_table(
            sampled,
            grid=grid,
            morph_levels=levels,
            min_photos_per_cell=min_photos_per_cell,
            transition_quantile=transition_quantile,
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
) -> pd.DataFrame:
    """Shuffle morph labels strictly within species while keeping locations fixed."""

    out = photos.reset_index(drop=True).copy()
    values = out[morph_col].astype(object).to_numpy(copy=True)
    species = out[species_col].astype(str).to_numpy()
    for species_name in np.unique(species):
        idx = np.flatnonzero(species == species_name)
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
    morph_levels: Sequence[str] | None = None,
) -> tuple[PersistenceResult, np.ndarray, float]:
    """Test excess recurrent-boundary concentration against a species-conditioned null."""

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
        )
        null[position] = result.concentration
    if not np.isfinite(null).all():
        raise ValueError("one or more null persistence concentrations are not estimable")
    p_upper = float((1 + np.count_nonzero(null >= observed.concentration)) / (n_permutations + 1))
    return observed, null, p_upper
