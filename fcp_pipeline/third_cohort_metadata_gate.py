from __future__ import annotations

import pandas as pd


def validate_selected_manifest(frame: pd.DataFrame, *, expected_n: int = 500) -> None:
    required = {"prospective_rank", "inat_taxon_id", "species"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"selected manifest missing columns: {sorted(missing)}")
    if len(frame) != int(expected_n):
        raise RuntimeError(f"selected manifest row count drift: {len(frame)} != {expected_n}")

    ranks = pd.to_numeric(frame["prospective_rank"], errors="raise").astype(int).tolist()
    if ranks != list(range(1, int(expected_n) + 1)):
        raise RuntimeError("prospective_rank drift")

    taxon = pd.to_numeric(frame["inat_taxon_id"], errors="raise").astype(int)
    species = frame["species"].astype(str)
    if taxon.nunique() != expected_n:
        raise RuntimeError("inat_taxon_id uniqueness drift")
    if species.nunique() != expected_n:
        raise RuntimeError("species uniqueness drift")
    if species.str.strip().eq("").any():
        raise RuntimeError("blank species identity")


def validate_fresh_rows(
    rows: pd.DataFrame,
    exclusion_observation_ids: set[int],
    exclusion_photo_ids: set[int],
) -> None:
    if rows.empty:
        return
    obs = pd.to_numeric(rows["observation_id"], errors="raise").astype(int)
    photo = pd.to_numeric(rows["photo_id"], errors="raise").astype(int)
    if obs.duplicated().any():
        raise RuntimeError("duplicate observation IDs in fresh metadata")
    if photo.duplicated().any():
        raise RuntimeError("duplicate photo IDs in fresh metadata")
    if set(obs) & set(exclusion_observation_ids) or set(photo) & set(exclusion_photo_ids):
        raise RuntimeError("prior experiment ID leaked into fresh metadata")


def seal_exact_denominator(
    audit: pd.DataFrame,
    rows: pd.DataFrame,
    *,
    target_n: int = 100,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "prospective_rank",
        "inat_taxon_id",
        "species",
        "full_fixed_n",
        "retained",
        "request_error",
    }
    missing = required - set(audit.columns)
    if missing:
        raise RuntimeError(f"species audit missing columns: {sorted(missing)}")

    species = audit.loc[audit["full_fixed_n"].astype(bool)].copy()
    species = species.sort_values("prospective_rank", kind="mergesort").reset_index(drop=True)
    if len(species):
        retained = pd.to_numeric(species["retained"], errors="raise").astype(int)
        if not retained.eq(int(target_n)).all():
            raise RuntimeError("full_fixed_n species do not all match exact target")

    ids = set(pd.to_numeric(species["inat_taxon_id"], errors="raise").astype(int))
    row_ids = pd.to_numeric(rows["inat_taxon_id"], errors="raise").astype(int)
    authorized = rows.loc[row_ids.isin(ids)].copy()
    if len(authorized):
        counts = authorized.groupby("inat_taxon_id", sort=False).size()
        if not counts.eq(int(target_n)).all() or len(counts) != len(species):
            raise RuntimeError("authorized row denominator is not exact")
        if {"prospective_rank", "h9_selection_order"}.issubset(authorized.columns):
            authorized = authorized.sort_values(
                ["prospective_rank", "h9_selection_order", "observation_id", "photo_id"],
                kind="mergesort",
                ignore_index=True,
            )
        else:
            authorized = authorized.sort_values(
                ["inat_taxon_id", "observation_id", "photo_id"],
                kind="mergesort",
                ignore_index=True,
            )
    elif len(species):
        raise RuntimeError("authorized denominator lost full-target species")

    return species, authorized


def metadata_gate_verdict(
    request_errors: int,
    request_attempts: int,
    full_species: int,
    *,
    min_full_species: int = 300,
    error_ceiling: float = 0.05,
) -> tuple[str, bool]:
    attempts = int(request_attempts)
    if attempts <= 0:
        raise ValueError("request_attempts must be positive")
    frac = int(request_errors) / attempts
    if frac > float(error_ceiling):
        return "THIRD_COHORT_METADATA_TRANSPORT_NOT_EVALUABLE", False
    if int(full_species) < int(min_full_species):
        return "THIRD_COHORT_METADATA_CAPACITY_NOT_EVALUABLE", False
    return "THIRD_COHORT_METADATA_GATE_PASS", True


__all__ = [
    "metadata_gate_verdict",
    "seal_exact_denominator",
    "validate_fresh_rows",
    "validate_selected_manifest",
]
