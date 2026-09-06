"""Repeated Global Flower-Colour Atlas scheduling and consensus primitives.

This module fixes the computational meaning of the 200 outer realizations before
any RGFCA flower-colour outcome is opened. It is deliberately outcome agnostic:
it schedules species/photo IDs, generates reproducible species-conditioned null
photo assignments, and aggregates realization fields by geographic opportunity.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Sequence

import numpy as np
import pandas as pd

from .global_monte_carlo import balanced_random_schedule


@dataclass(frozen=True)
class RepeatedAtlasSchedule:
    """Frozen balanced species/photo schedule expressed in stable photo IDs."""

    species_labels: tuple[str, ...]
    outer_species: np.ndarray
    outer_photo_ids: np.ndarray
    species_inclusion_counts: np.ndarray
    photo_max_inclusion_imbalance: int
    n_outer: int
    species_per_outer: int
    photos_per_species: int
    species_seed: int
    photo_master_seed: int


@dataclass(frozen=True)
class ConsensusField:
    field: np.ndarray
    aggregate_opportunity: np.ndarray
    evaluable_outer_counts: np.ndarray
    concentration: float
    weighted_mean: float


def stable_seed(master_seed: int, *parts: object) -> int:
    """Derive an input-order-independent NumPy seed from canonical labels."""
    text = ":".join([str(int(master_seed)), *(str(part) for part in parts)])
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    # Stay within signed 63-bit range for broad RNG/backend compatibility.
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def _canonical_photo_frame(
    photo_ids: Sequence[int],
    photo_species: Sequence[object],
) -> pd.DataFrame:
    ids = np.asarray(photo_ids)
    species = np.asarray(photo_species)
    if ids.ndim != 1 or species.ndim != 1 or len(ids) != len(species):
        raise ValueError("photo_ids and photo_species must be equal-length vectors")
    if len(ids) == 0:
        raise ValueError("photo pool must be non-empty")
    try:
        ids = ids.astype(np.int64)
    except (TypeError, ValueError) as exc:
        raise ValueError("photo_ids must be integer-like") from exc
    if len(np.unique(ids)) != len(ids):
        raise ValueError("photo_ids must be globally unique")
    labels = np.asarray([str(value) for value in species], dtype=object)
    if any(label == "" for label in labels):
        raise ValueError("species labels must be non-empty")
    frame = pd.DataFrame({"photo_id": ids, "species": labels})
    return frame.sort_values(["species", "photo_id"], kind="mergesort").reset_index(drop=True)


def build_repeated_atlas_schedule(
    photo_ids: Sequence[int],
    photo_species: Sequence[object],
    *,
    n_outer: int = 200,
    species_per_outer: int = 250,
    photos_per_species: int = 20,
    minimum_pool_photos_per_species: int = 40,
    species_seed: int = 2026090401,
    photo_master_seed: int = 2026090402,
) -> RepeatedAtlasSchedule:
    """Create balanced outer species and within-species photo schedules.

    Every input species must already satisfy the frozen postmeasurement pool gate.
    The function fails closed rather than silently dropping undersampled species.
    Photo schedules are built only for the number of outer realizations in which a
    species actually appears, so inclusion imbalance is audited on used draws.
    """
    n_outer = int(n_outer)
    species_per_outer = int(species_per_outer)
    photos_per_species = int(photos_per_species)
    minimum_pool = int(minimum_pool_photos_per_species)
    if n_outer < 1 or species_per_outer < 1 or photos_per_species < 2:
        raise ValueError("outer/sample sizes must be positive and photos_per_species >=2")
    if minimum_pool < photos_per_species:
        raise ValueError("minimum pool cannot be smaller than photos_per_species")

    frame = _canonical_photo_frame(photo_ids, photo_species)
    counts = frame.groupby("species", sort=True).size()
    if (counts < minimum_pool).any():
        bad = counts[counts < minimum_pool]
        raise ValueError(
            f"{len(bad)} species fail minimum_pool_photos_per_species={minimum_pool}; "
            "postmeasurement filtering must occur before schedule construction"
        )
    species_labels = tuple(sorted(counts.index.astype(str).tolist()))
    if len(species_labels) < species_per_outer:
        raise ValueError("eligible species count is below species_per_outer")

    species_schedule = balanced_random_schedule(
        species_labels,
        n_replicates=n_outer,
        items_per_replicate=species_per_outer,
        seed=int(species_seed),
    )
    outer_species = np.empty((n_outer, species_per_outer), dtype=object)
    for r in range(n_outer):
        outer_species[r] = [species_labels[int(i)] for i in species_schedule.draws[r]]

    outer_photo_ids = np.empty(
        (n_outer, species_per_outer, photos_per_species), dtype=np.int64
    )
    occurrence_cursor = {label: 0 for label in species_labels}
    photo_schedules: dict[str, np.ndarray] = {}
    max_photo_imbalance = 0
    for species_index, label in enumerate(species_labels):
        ids = frame.loc[frame["species"] == label, "photo_id"].astype(np.int64).to_numpy()
        occurrences = int(species_schedule.inclusion_counts[species_index])
        schedule = balanced_random_schedule(
            [str(int(value)) for value in ids],
            n_replicates=occurrences,
            items_per_replicate=photos_per_species,
            seed=stable_seed(int(photo_master_seed), label),
        )
        photo_schedules[label] = np.asarray(
            [[int(schedule.items[int(i)]) for i in row] for row in schedule.draws],
            dtype=np.int64,
        )
        max_photo_imbalance = max(max_photo_imbalance, schedule.max_inclusion_imbalance)

    for r in range(n_outer):
        for slot in range(species_per_outer):
            label = str(outer_species[r, slot])
            occurrence = occurrence_cursor[label]
            outer_photo_ids[r, slot] = photo_schedules[label][occurrence]
            occurrence_cursor[label] = occurrence + 1

    for species_index, label in enumerate(species_labels):
        expected = int(species_schedule.inclusion_counts[species_index])
        if occurrence_cursor[label] != expected:
            raise RuntimeError("internal occurrence accounting failure")
    if max_photo_imbalance > 1:
        raise RuntimeError("used photo inclusion counts differ by more than one within a species")

    # Hard no-duplication audit within every species draw.
    sorted_draws = np.sort(outer_photo_ids, axis=2)
    if np.any(np.diff(sorted_draws, axis=2) == 0):
        raise RuntimeError("a within-species outer draw contains a duplicate photo")

    return RepeatedAtlasSchedule(
        species_labels=species_labels,
        outer_species=outer_species,
        outer_photo_ids=outer_photo_ids,
        species_inclusion_counts=species_schedule.inclusion_counts.copy(),
        photo_max_inclusion_imbalance=int(max_photo_imbalance),
        n_outer=n_outer,
        species_per_outer=species_per_outer,
        photos_per_species=photos_per_species,
        species_seed=int(species_seed),
        photo_master_seed=int(photo_master_seed),
    )


def schedule_audit(schedule: RepeatedAtlasSchedule) -> dict[str, object]:
    counts = np.asarray(schedule.species_inclusion_counts, dtype=np.int64)
    return {
        "eligible_species": len(schedule.species_labels),
        "outer_realizations": schedule.n_outer,
        "species_per_outer": schedule.species_per_outer,
        "photos_per_species": schedule.photos_per_species,
        "total_species_inclusions": int(counts.sum()),
        "minimum_species_inclusions": int(counts.min()),
        "maximum_species_inclusions": int(counts.max()),
        "species_max_inclusion_imbalance": int(counts.max() - counts.min()),
        "photo_max_inclusion_imbalance_within_species": schedule.photo_max_inclusion_imbalance,
        "species_seed": schedule.species_seed,
        "photo_master_seed": schedule.photo_master_seed,
    }


def null_source_photo_ids(
    photo_ids: Sequence[int],
    photo_species: Sequence[object],
    *,
    permutation_index: int,
    master_seed: int = 2026090403,
) -> np.ndarray:
    """Return source-photo IDs assigned to each target photo under one null.

    Output order matches the caller's input order, but the random assignment is
    invariant to that row order because permutations are generated on photo IDs
    sorted within canonical species labels. The same mapping can therefore be
    reused across every outer realization of one null replicate.
    """
    permutation_index = int(permutation_index)
    if permutation_index < 0:
        raise ValueError("permutation_index must be non-negative")
    frame = _canonical_photo_frame(photo_ids, photo_species)
    mapping: dict[int, int] = {}
    for label, group in frame.groupby("species", sort=True):
        ids = group["photo_id"].astype(np.int64).to_numpy()
        rng = np.random.default_rng(stable_seed(int(master_seed), permutation_index, label))
        source = ids[rng.permutation(len(ids))]
        mapping.update({int(target): int(src) for target, src in zip(ids, source)})
    original = np.asarray(photo_ids, dtype=np.int64)
    return np.asarray([mapping[int(value)] for value in original], dtype=np.int64)


def consensus_field(
    fields: Sequence[Sequence[float]],
    opportunities: Sequence[Sequence[float]],
) -> ConsensusField:
    """Aggregate outer fields using the pre-frozen opportunity-weighted rule."""
    values = np.asarray(fields, dtype=float)
    opp = np.asarray(opportunities, dtype=float)
    if values.ndim != 2 or opp.ndim != 2 or values.shape != opp.shape:
        raise ValueError("fields and opportunities must have the same 2D shape")
    if values.shape[0] < 1 or values.shape[1] < 1:
        raise ValueError("fields must contain at least one realization and one cell")
    if np.any(np.isfinite(opp) & (opp < 0)):
        raise ValueError("opportunity cannot be negative")
    evaluable = np.isfinite(values) & np.isfinite(opp) & (opp > 0)
    aggregate_opp = np.where(evaluable, opp, 0.0).sum(axis=0)
    aggregate_num = np.where(evaluable, values * opp, 0.0).sum(axis=0)
    counts = evaluable.sum(axis=0).astype(np.int64)
    out = np.full(values.shape[1], np.nan, dtype=float)
    keep = aggregate_opp > 0
    out[keep] = aggregate_num[keep] / aggregate_opp[keep]
    if not np.any(keep):
        return ConsensusField(
            field=out,
            aggregate_opportunity=aggregate_opp,
            evaluable_outer_counts=counts,
            concentration=float("nan"),
            weighted_mean=float("nan"),
        )
    w = aggregate_opp[keep]
    x = out[keep]
    mean = float(np.average(x, weights=w))
    concentration = float(np.average(np.square(x - mean), weights=w))
    return ConsensusField(
        field=out,
        aggregate_opportunity=aggregate_opp,
        evaluable_outer_counts=counts,
        concentration=concentration,
        weighted_mean=mean,
    )


def running_consensus(
    fields: Sequence[Sequence[float]],
    opportunities: Sequence[Sequence[float]],
    *,
    checkpoints: Sequence[int] = (25, 50, 100, 150, 200),
) -> dict[int, ConsensusField]:
    values = np.asarray(fields, dtype=float)
    opp = np.asarray(opportunities, dtype=float)
    if values.shape != opp.shape or values.ndim != 2:
        raise ValueError("fields/opportunities must be identically shaped matrices")
    output: dict[int, ConsensusField] = {}
    for checkpoint_raw in checkpoints:
        checkpoint = int(checkpoint_raw)
        if checkpoint < 1 or checkpoint > values.shape[0]:
            raise ValueError("running-consensus checkpoint outside realized range")
        output[checkpoint] = consensus_field(values[:checkpoint], opp[:checkpoint])
    return output


def odd_even_consensus(
    fields: Sequence[Sequence[float]],
    opportunities: Sequence[Sequence[float]],
) -> tuple[ConsensusField, ConsensusField]:
    values = np.asarray(fields, dtype=float)
    opp = np.asarray(opportunities, dtype=float)
    if values.shape != opp.shape or values.ndim != 2:
        raise ValueError("fields/opportunities must be identically shaped matrices")
    if values.shape[0] < 2:
        raise ValueError("at least two outer realizations are required")
    # Human-facing realization numbering is 1..R: odd = Python rows 0,2,...
    return consensus_field(values[0::2], opp[0::2]), consensus_field(values[1::2], opp[1::2])


__all__ = [
    "ConsensusField",
    "RepeatedAtlasSchedule",
    "build_repeated_atlas_schedule",
    "consensus_field",
    "null_source_photo_ids",
    "odd_even_consensus",
    "running_consensus",
    "schedule_audit",
    "stable_seed",
]
