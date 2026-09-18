"""Observer-disjoint reliability for species-level categorical diversity."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .diversity import gini_simpson


@dataclass(frozen=True)
class ReliabilityResult:
    partitions: pd.DataFrame
    summary: dict[str, float]


def _stable_hash_int(seed: int, species: str, observer: str) -> int:
    """FCP-compatible deterministic observer tie-break."""
    text = " | ".join(str(x) for x in (seed, species, observer))
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest(), 16)


def _as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"true", "1", "yes"})
    )


def _ccc(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    cov = float(np.cov(x, y, ddof=1)[0, 1])
    den = vx + vy + (mx - my) ** 2
    return float(2.0 * cov / den) if den > 0 else float("nan")


def _spearman_brown(rho: float) -> float:
    if not np.isfinite(rho) or rho <= -1.0:
        return float("nan")
    return float(2.0 * rho / (1.0 + rho))


def _profile_species(
    frame: pd.DataFrame,
    *,
    species_col: str,
    observer_col: str,
    state_col: str,
    states: Sequence[object],
    classifiable_col: str | None,
    min_full_classifiable: int,
    minimum_distinct_observers: int,
) -> dict[str, list[tuple[str, int, np.ndarray]]]:
    profiles: dict[str, list[tuple[str, int, np.ndarray]]] = {}
    state_index = {state: i for i, state in enumerate(states)}
    work = frame.copy()
    work[species_col] = work[species_col].fillna("").astype(str).str.strip()
    work[observer_col] = work[observer_col].fillna("").astype(str).str.strip()
    work = work.loc[work[species_col].ne("") & work[observer_col].ne("")].copy()
    if classifiable_col is None:
        work["_disttrait_classifiable"] = work[state_col].isin(states)
    else:
        work["_disttrait_classifiable"] = _as_bool(work[classifiable_col]) & work[state_col].isin(states)

    for species, group in work.groupby(species_col, sort=False):
        classifiable = group["_disttrait_classifiable"]
        if int(classifiable.sum()) < int(min_full_classifiable):
            continue
        if group[observer_col].nunique() < int(minimum_distinct_observers):
            continue
        rows: list[tuple[str, int, np.ndarray]] = []
        for observer, og in group.groupby(observer_col, sort=False):
            counts = np.zeros(len(states), dtype=np.int64)
            class_rows = og.loc[og["_disttrait_classifiable"], state_col]
            for value, n in class_rows.value_counts(dropna=False).items():
                if value in state_index:
                    counts[state_index[value]] = int(n)
            rows.append((str(observer), int(len(og)), counts))
        profiles[str(species)] = rows
    return profiles


def _partition_counts(
    profile: list[tuple[str, int, np.ndarray]],
    *,
    species: str,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    items = [
        (observer, n_all, counts, _stable_hash_int(seed, species, observer))
        for observer, n_all, counts in profile
    ]
    items.sort(key=lambda z: (-z[1], z[3]))
    totals = [0, 0]
    state_totals = [
        np.zeros_like(profile[0][2], dtype=np.int64),
        np.zeros_like(profile[0][2], dtype=np.int64),
    ]
    used = [0, 0]
    for _, n_all, counts, h in items:
        if totals[0] < totals[1]:
            side = 0
        elif totals[1] < totals[0]:
            side = 1
        else:
            side = h & 1
        totals[side] += int(n_all)
        state_totals[side] += counts
        used[side] += 1
    if min(used) == 0:
        raise RuntimeError(f"{species}: observer partition produced an empty side")
    return state_totals[0], state_totals[1]


def observer_disjoint_reliability(
    frame: pd.DataFrame,
    *,
    species_col: str,
    observer_col: str,
    state_col: str,
    states: Sequence[object],
    classifiable_col: str | None = None,
    n_partitions: int = 200,
    base_seed: int = 20260913,
    min_full_classifiable: int = 40,
    min_half_classifiable: int = 20,
    minimum_distinct_observers: int = 2,
) -> ReliabilityResult:
    """Repeatedly estimate diversity from completely disjoint observer sets."""
    profiles = _profile_species(
        frame,
        species_col=species_col,
        observer_col=observer_col,
        state_col=state_col,
        states=states,
        classifiable_col=classifiable_col,
        min_full_classifiable=min_full_classifiable,
        minimum_distinct_observers=minimum_distinct_observers,
    )
    rows: list[dict[str, float | int]] = []
    for split_index in range(1, int(n_partitions) + 1):
        seed = int(base_seed) + split_index
        pairs: list[tuple[float, float]] = []
        for species, profile in profiles.items():
            a, b = _partition_counts(profile, species=species, seed=seed)
            if int(a.sum()) < int(min_half_classifiable) or int(b.sum()) < int(min_half_classifiable):
                continue
            pairs.append((gini_simpson(a), gini_simpson(b)))

        if pairs:
            x = np.asarray([p[0] for p in pairs], dtype=float)
            y = np.asarray([p[1] for p in pairs], dtype=float)
            rho = float(spearmanr(x, y).statistic) if len(x) >= 3 else float("nan")
            mae = float(np.mean(np.abs(x - y)))
            bias = float(np.mean(x - y))
            ccc = _ccc(x, y)
        else:
            rho = mae = bias = ccc = float("nan")
        rows.append(
            {
                "split_index": split_index,
                "seed": seed,
                "paired_n": int(len(pairs)),
                "rho": rho,
                "ccc": ccc,
                "spearman_brown": _spearman_brown(rho),
                "mae": mae,
                "bias_a_minus_b": bias,
            }
        )

    partitions = pd.DataFrame(rows)
    defined = partitions["rho"].dropna()
    summary = {
        "eligible_species": float(len(profiles)),
        "partitions": float(len(partitions)),
        "paired_n_median": float(partitions["paired_n"].median()),
        "rho_median": float(defined.median()) if len(defined) else float("nan"),
        "rho_q05": float(defined.quantile(0.05)) if len(defined) else float("nan"),
        "rho_q95": float(defined.quantile(0.95)) if len(defined) else float("nan"),
        "ccc_median": float(partitions["ccc"].median()),
        "spearman_brown_median": float(partitions["spearman_brown"].median()),
        "mae_median": float(partitions["mae"].median()),
    }
    return ReliabilityResult(partitions=partitions, summary=summary)
