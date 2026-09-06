"""Balanced Monte Carlo coverage primitives for the prospective global barrier atlas.

The purpose of this module is computational and design-based: bound each analysis
replicate while spreading inclusion opportunity as evenly as possible across the
full eligible species/photo pool.  It contains no flower-colour outcome logic and
can therefore be validated before fresh biological outcomes are opened.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class BalancedSchedule:
    """A fixed-size repeated sampling schedule with inclusion-count audit."""

    draws: np.ndarray
    items: tuple[str, ...]
    inclusion_counts: np.ndarray
    seed: int

    @property
    def n_replicates(self) -> int:
        return int(self.draws.shape[0])

    @property
    def items_per_replicate(self) -> int:
        return int(self.draws.shape[1])

    @property
    def max_inclusion_imbalance(self) -> int:
        if len(self.inclusion_counts) == 0:
            return 0
        return int(self.inclusion_counts.max() - self.inclusion_counts.min())


def _validate_unique_items(items: Sequence[object]) -> tuple[str, ...]:
    labels = tuple(str(value) for value in items)
    if not labels:
        raise ValueError("items must be non-empty")
    if len(set(labels)) != len(labels):
        raise ValueError("items must be unique")
    return labels


def balanced_random_schedule(
    items: Sequence[object],
    *,
    n_replicates: int,
    items_per_replicate: int,
    seed: int,
) -> BalancedSchedule:
    """Build seeded random draws while keeping long-run inclusion nearly equal.

    Each replicate samples without replacement.  At every replicate, items with
    the lowest current inclusion count receive priority; random keys break ties.
    This gives every eligible item repeated non-zero opportunity while preventing
    a few abundant/easy items from dominating simply through repeated sampling.
    """

    labels = _validate_unique_items(items)
    n = len(labels)
    n_replicates = int(n_replicates)
    items_per_replicate = int(items_per_replicate)
    if n_replicates <= 0:
        raise ValueError("n_replicates must be positive")
    if items_per_replicate <= 0:
        raise ValueError("items_per_replicate must be positive")
    if items_per_replicate > n:
        raise ValueError("items_per_replicate cannot exceed the eligible item count")

    rng = np.random.default_rng(int(seed))
    counts = np.zeros(n, dtype=np.int64)
    draws = np.empty((n_replicates, items_per_replicate), dtype=np.int64)

    for replicate in range(n_replicates):
        random_tie_break = rng.random(n)
        order = np.lexsort((random_tie_break, counts))
        chosen = order[:items_per_replicate]
        draws[replicate] = chosen
        counts[chosen] += 1

    # The greedy lowest-count rule should never allow a gap greater than one.
    if int(counts.max() - counts.min()) > 1:
        raise RuntimeError("balanced schedule inclusion counts differ by more than one")

    return BalancedSchedule(
        draws=draws,
        items=labels,
        inclusion_counts=counts,
        seed=int(seed),
    )


def draw_labels(schedule: BalancedSchedule, replicate: int) -> tuple[str, ...]:
    """Return item labels for one replicate."""

    replicate = int(replicate)
    if replicate < 0 or replicate >= schedule.n_replicates:
        raise IndexError("replicate out of range")
    return tuple(schedule.items[int(index)] for index in schedule.draws[replicate])


def pair_count(n: int) -> int:
    """Number of unordered pairs among ``n`` photographs."""

    n = int(n)
    if n < 0:
        raise ValueError("n must be non-negative")
    return n * (n - 1) // 2


def opportunity_normalized_field(
    numerator: Sequence[float],
    opportunity: Sequence[float],
    *,
    minimum_opportunity: float = 0.0,
) -> np.ndarray:
    """Normalize colour-edge support by geometric opportunity.

    Cells at or below ``minimum_opportunity`` are returned as NaN rather than as
    biological zeros.  Geometry can therefore be held fixed under the null while
    colour contributions are permuted within species.
    """

    num = np.asarray(numerator, dtype=float)
    opp = np.asarray(opportunity, dtype=float)
    if num.shape != opp.shape:
        raise ValueError("numerator and opportunity must have identical shapes")
    if np.any(opp < 0):
        raise ValueError("opportunity cannot be negative")
    threshold = float(minimum_opportunity)
    out = np.full(num.shape, np.nan, dtype=float)
    keep = np.isfinite(num) & np.isfinite(opp) & (opp > threshold)
    out[keep] = num[keep] / opp[keep]
    return out


def schedule_audit(schedule: BalancedSchedule) -> dict[str, object]:
    """Return a JSON-serializable inclusion audit."""

    counts = schedule.inclusion_counts
    return {
        "eligible_items": int(len(schedule.items)),
        "n_replicates": schedule.n_replicates,
        "items_per_replicate": schedule.items_per_replicate,
        "total_item_inclusions": int(counts.sum()),
        "minimum_inclusions_per_item": int(counts.min()),
        "maximum_inclusions_per_item": int(counts.max()),
        "max_inclusion_imbalance": schedule.max_inclusion_imbalance,
        "seed": int(schedule.seed),
    }


__all__ = [
    "BalancedSchedule",
    "balanced_random_schedule",
    "draw_labels",
    "opportunity_normalized_field",
    "pair_count",
    "schedule_audit",
]
