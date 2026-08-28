"""Deterministic calibration/holdout freezing for the JBI Chapter 1 photo gate.

The split is intentionally based only on species identity and a stable photograph ID.
No flower-colour, visibility, segmentation, geography, observer, or date field can affect
assignment. This keeps the 480/720 gate independent of downstream phenotype outcomes.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Iterable

import pandas as pd


SPLIT_VERSION = "jbi-ch1-photo-split-v1"
DEFAULT_SALT = "fcp-jbi-ch1-photo-split-v1"

_FORBIDDEN_OUTCOME_COLUMNS = {
    "flower_colour",
    "flower_color",
    "flower_colour_state",
    "flower_color_state",
    "colour_state",
    "color_state",
    "visibility",
    "flower_visibility",
    "segmentation",
    "segmentation_status",
    "evaluable",
    "not_evaluable",
    "unresolved",
}


@dataclass(frozen=True)
class SplitSpec:
    expected_species: int = 6
    photographs_per_species: int = 200
    calibration_per_species: int = 80
    evaluation_per_species: int = 120
    salt: str = DEFAULT_SALT

    def validate(self) -> None:
        if self.expected_species <= 0:
            raise ValueError("expected_species must be positive")
        if self.photographs_per_species <= 0:
            raise ValueError("photographs_per_species must be positive")
        if self.calibration_per_species < 0 or self.evaluation_per_species < 0:
            raise ValueError("split counts must be non-negative")
        if self.calibration_per_species + self.evaluation_per_species != self.photographs_per_species:
            raise ValueError(
                "calibration_per_species + evaluation_per_species must equal photographs_per_species"
            )
        if not self.salt:
            raise ValueError("salt must be non-empty")


def _normalise_column_name(name: object) -> str:
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def find_outcome_columns(columns: Iterable[object]) -> list[str]:
    found: list[str] = []
    for col in columns:
        norm = _normalise_column_name(col)
        if norm in _FORBIDDEN_OUTCOME_COLUMNS:
            found.append(str(col))
    return sorted(found)


def validate_source_manifest(
    frame: pd.DataFrame,
    *,
    species_col: str,
    photo_id_col: str,
    spec: SplitSpec = SplitSpec(),
    reject_outcome_columns: bool = True,
) -> pd.DataFrame:
    """Validate and normalise the 6 x 200 acquisition manifest.

    The returned frame is a copy. Species and photograph identifiers are normalised to
    stripped strings so that split assignment is stable across CSV dtype inference.
    """

    spec.validate()
    missing = [c for c in (species_col, photo_id_col) if c not in frame.columns]
    if missing:
        raise ValueError(f"source manifest missing required columns: {missing}")

    if reject_outcome_columns:
        forbidden = find_outcome_columns(frame.columns)
        if forbidden:
            raise ValueError(
                "source manifest already contains downstream measurement outcome columns: "
                + ", ".join(forbidden)
            )

    out = frame.copy()
    if out[species_col].isna().any() or out[photo_id_col].isna().any():
        raise ValueError("species and photo IDs must be non-missing")

    out[species_col] = out[species_col].astype(str).str.strip()
    out[photo_id_col] = out[photo_id_col].astype(str).str.strip()
    if (out[species_col] == "").any() or (out[photo_id_col] == "").any():
        raise ValueError("species and photo IDs must be non-empty")

    if out[photo_id_col].duplicated().any():
        dup = out.loc[out[photo_id_col].duplicated(keep=False), photo_id_col].iloc[0]
        raise ValueError(f"photo IDs must be globally unique; duplicate={dup}")

    counts = out.groupby(species_col, sort=True)[photo_id_col].size()
    if len(counts) != spec.expected_species:
        raise ValueError(
            f"expected {spec.expected_species} species, found {len(counts)}"
        )
    bad = counts[counts != spec.photographs_per_species]
    if not bad.empty:
        details = ", ".join(f"{idx}={int(n)}" for idx, n in bad.items())
        raise ValueError(
            f"each species must have exactly {spec.photographs_per_species} photos; {details}"
        )

    expected_total = spec.expected_species * spec.photographs_per_species
    if len(out) != expected_total:
        raise ValueError(f"expected {expected_total} rows, found {len(out)}")
    return out


def _assignment_digest(species: str, photo_id: str, salt: str) -> str:
    key = "\x1f".join((salt, species, photo_id)).encode("utf-8")
    return hashlib.sha256(key).hexdigest()


def canonical_id_hash(frame: pd.DataFrame, *, species_col: str, photo_id_col: str) -> str:
    rows = sorted(
        (str(s).strip(), str(p).strip())
        for s, p in frame[[species_col, photo_id_col]].itertuples(index=False, name=None)
    )
    payload = "\n".join(f"{s}\t{p}" for s, p in rows).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def assignment_hash(frame: pd.DataFrame, *, species_col: str, photo_id_col: str) -> str:
    required = {species_col, photo_id_col, "split"}
    if not required.issubset(frame.columns):
        raise ValueError(f"assignment frame missing columns: {sorted(required - set(frame.columns))}")
    rows = sorted(
        (str(s).strip(), str(p).strip(), str(split).strip())
        for s, p, split in frame[[species_col, photo_id_col, "split"]].itertuples(
            index=False, name=None
        )
    )
    payload = "\n".join(f"{s}\t{p}\t{split}" for s, p, split in rows).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def freeze_photo_split(
    frame: pd.DataFrame,
    *,
    species_col: str,
    photo_id_col: str,
    spec: SplitSpec = SplitSpec(),
) -> pd.DataFrame:
    """Assign exactly 80 calibration and 120 held-out photos per species.

    Within each species, rows are ranked by SHA-256(salt, species, photo_id). The
    smallest 80 hashes enter calibration; the remaining 120 are held out. No other
    metadata can change assignment.
    """

    out = validate_source_manifest(
        frame,
        species_col=species_col,
        photo_id_col=photo_id_col,
        spec=spec,
    )
    out["split_rank_hash"] = [
        _assignment_digest(species, photo_id, spec.salt)
        for species, photo_id in out[[species_col, photo_id_col]].itertuples(index=False, name=None)
    ]
    out["split"] = ""

    for species, group in out.groupby(species_col, sort=True):
        ordered = group.sort_values(["split_rank_hash", photo_id_col], kind="mergesort")
        calibration_idx = ordered.index[: spec.calibration_per_species]
        evaluation_idx = ordered.index[spec.calibration_per_species :]
        out.loc[calibration_idx, "split"] = "calibration"
        out.loc[evaluation_idx, "split"] = "evaluation"

    validate_frozen_split(
        out,
        species_col=species_col,
        photo_id_col=photo_id_col,
        spec=spec,
    )
    return out.sort_values([species_col, "split", "split_rank_hash", photo_id_col], kind="mergesort").reset_index(drop=True)


def validate_frozen_split(
    frame: pd.DataFrame,
    *,
    species_col: str,
    photo_id_col: str,
    spec: SplitSpec = SplitSpec(),
) -> None:
    spec.validate()
    required = {species_col, photo_id_col, "split", "split_rank_hash"}
    if not required.issubset(frame.columns):
        raise ValueError(f"frozen split missing columns: {sorted(required - set(frame.columns))}")
    if frame[photo_id_col].astype(str).duplicated().any():
        raise ValueError("frozen split contains duplicate photo IDs")
    if set(frame["split"].astype(str)) != {"calibration", "evaluation"}:
        raise ValueError("split must contain exactly calibration and evaluation")

    counts = frame.groupby([species_col, "split"], sort=True).size().unstack(fill_value=0)
    if len(counts) != spec.expected_species:
        raise ValueError("frozen split has the wrong number of species")
    for species, row in counts.iterrows():
        if int(row.get("calibration", 0)) != spec.calibration_per_species:
            raise ValueError(f"{species}: wrong calibration count")
        if int(row.get("evaluation", 0)) != spec.evaluation_per_species:
            raise ValueError(f"{species}: wrong evaluation count")

    # Recompute every rank hash: a changed salt/species/photo ID must invalidate the freeze.
    expected = [
        _assignment_digest(str(s).strip(), str(p).strip(), spec.salt)
        for s, p in frame[[species_col, photo_id_col]].itertuples(index=False, name=None)
    ]
    observed = frame["split_rank_hash"].astype(str).tolist()
    if observed != expected:
        raise ValueError("split_rank_hash does not match the frozen split specification")
