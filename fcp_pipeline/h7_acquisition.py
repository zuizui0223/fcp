"""Prospective fresh metadata acquisition for the balanced H7 ITV experiment.

This module never opens image pixels. Species and target cells are already frozen.
Each species x cell target receives exactly one iNaturalist query attempt. A request
error is recorded as zero returned candidates rather than retried or replaced, so a
failed target cannot become a favourable random redraw.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Any, Mapping, Sequence

import pandas as pd

from .random_photo_pool import (
    DEFAULT_ALLOWED_PHOTO_LICENSES,
    ObservationClient,
    inat_query_for_cell,
    parse_candidate_observation,
)
from .shared_transition_surface import EqualAreaGrid


@dataclass(frozen=True)
class H7MetadataFreeze:
    observations: pd.DataFrame
    target_audit: pd.DataFrame
    species_support: pd.DataFrame
    manifest: dict[str, Any]


def _selection_hash(seed: int, taxon_id: int, cell_id: int, observation_id: int) -> str:
    payload = f"{int(seed)}|{int(taxon_id)}|{int(cell_id)}|{int(observation_id)}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_targets(targets: pd.DataFrame, *, expected_species: int, expected_targets: int) -> pd.DataFrame:
    required = {
        "species",
        "inat_taxon_id",
        "target_cell_order",
        "cell_id",
        "cell_center_latitude",
        "cell_center_longitude",
    }
    missing = required - set(targets.columns)
    if missing:
        raise ValueError(f"target table missing columns: {sorted(missing)}")
    out = targets.copy()
    out["inat_taxon_id"] = pd.to_numeric(out["inat_taxon_id"], errors="raise").astype("int64")
    out["cell_id"] = pd.to_numeric(out["cell_id"], errors="raise").astype("int64")
    out["target_cell_order"] = pd.to_numeric(out["target_cell_order"], errors="raise").astype("int64")
    if len(out) != int(expected_targets):
        raise ValueError(f"expected {expected_targets} targets, found {len(out)}")
    if out["inat_taxon_id"].nunique() != int(expected_species):
        raise ValueError(f"expected {expected_species} species")
    if out.duplicated(["inat_taxon_id", "cell_id"]).any():
        raise ValueError("duplicate species x cell target")
    return out.sort_values(["inat_taxon_id", "target_cell_order"], kind="mergesort").reset_index(drop=True)


def _exclusion_sets(exclusions: pd.DataFrame) -> tuple[set[int], set[int]]:
    required = {"observation_id", "photo_id"}
    missing = required - set(exclusions.columns)
    if missing:
        raise ValueError(f"exclusion ledger missing columns: {sorted(missing)}")
    obs = set(pd.to_numeric(exclusions["observation_id"], errors="raise").astype("int64").tolist())
    photos = set(pd.to_numeric(exclusions["photo_id"], errors="raise").astype("int64").tolist())
    return obs, photos


def _observer_capped_selection(
    rows: list[dict[str, object]],
    *,
    observer_cap: int,
    retained_cap: int,
    seed: int,
    taxon_id: int,
    cell_id: int,
) -> list[dict[str, object]]:
    ranked = []
    for row in rows:
        item = dict(row)
        item["selection_hash"] = _selection_hash(seed, taxon_id, cell_id, int(item["observation_id"]))
        ranked.append(item)
    ranked.sort(key=lambda r: (str(r["selection_hash"]), int(r["observation_id"]), int(r["photo_id"])))
    by_observer: dict[str, int] = {}
    selected: list[dict[str, object]] = []
    for row in ranked:
        observer = str(row["observer_id"])
        if by_observer.get(observer, 0) >= int(observer_cap):
            continue
        selected.append(row)
        by_observer[observer] = by_observer.get(observer, 0) + 1
        if len(selected) >= int(retained_cap):
            break
    return selected


def freeze_h7_fresh_metadata(
    *,
    client: ObservationClient,
    targets: pd.DataFrame,
    exclusions: pd.DataFrame,
    grid: EqualAreaGrid,
    per_page: int = 50,
    observer_cap: int = 2,
    retained_cap: int = 30,
    selection_seed: int = 20260909,
    expected_species: int = 52,
    expected_targets: int = 312,
    required_species_for_gate: int = 35,
    required_full_cells_per_species: int = 5,
    maximum_positional_accuracy_m: int = 5000,
    flowering_term_id: int = 12,
    flowering_term_value_id: int = 13,
    allowed_photo_licenses: Sequence[str] = DEFAULT_ALLOWED_PHOTO_LICENSES,
) -> H7MetadataFreeze:
    """Execute exactly one metadata query attempt for every frozen H7 target."""

    targets = _validate_targets(targets, expected_species=expected_species, expected_targets=expected_targets)
    excluded_observations, excluded_photos = _exclusion_sets(exclusions)
    allowed = frozenset(str(x).casefold() for x in allowed_photo_licenses)
    seen_fresh_observations: set[int] = set()
    seen_fresh_photos: set[int] = set()
    kept_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []

    for target in targets.itertuples(index=False):
        taxon_id = int(target.inat_taxon_id)
        cell_id = int(target.cell_id)
        params = inat_query_for_cell(
            grid,
            cell_id,
            per_page=int(per_page),
            taxon_id=taxon_id,
            flowering_term_id=int(flowering_term_id),
            flowering_term_value_id=int(flowering_term_value_id),
            maximum_positional_accuracy_m=int(maximum_positional_accuracy_m),
            allowed_photo_licenses=tuple(allowed_photo_licenses),
        )
        request_error = ""
        try:
            payload = client.observations(params)
            raw_results = payload.get("results") or []
            if not isinstance(raw_results, Sequence) or isinstance(raw_results, (str, bytes)):
                raise RuntimeError("iNaturalist response results is not a sequence")
        except Exception as exc:  # one attempt only: retain failure as zero support
            raw_results = []
            request_error = f"{type(exc).__name__}: {exc}"

        locally_eligible: list[dict[str, object]] = []
        wrong_taxon = 0
        prior_excluded = 0
        fresh_duplicate = 0
        local_rejected = 0
        for observation in raw_results:
            if not isinstance(observation, Mapping):
                local_rejected += 1
                continue
            parsed = parse_candidate_observation(
                observation,
                expected_cell_id=cell_id,
                grid=grid,
                maximum_positional_accuracy_m=float(maximum_positional_accuracy_m),
                allowed_photo_licenses=allowed,
            )
            if parsed is None:
                local_rejected += 1
                continue
            if int(parsed["inat_taxon_id"]) != taxon_id:
                wrong_taxon += 1
                continue
            observation_id = int(parsed["observation_id"])
            photo_id = int(parsed["photo_id"])
            if observation_id in excluded_observations or photo_id in excluded_photos:
                prior_excluded += 1
                continue
            if observation_id in seen_fresh_observations or photo_id in seen_fresh_photos:
                fresh_duplicate += 1
                continue
            locally_eligible.append(parsed)

        selected = _observer_capped_selection(
            locally_eligible,
            observer_cap=int(observer_cap),
            retained_cap=int(retained_cap),
            seed=int(selection_seed),
            taxon_id=taxon_id,
            cell_id=cell_id,
        )
        for row in selected:
            observation_id = int(row["observation_id"])
            photo_id = int(row["photo_id"])
            seen_fresh_observations.add(observation_id)
            seen_fresh_photos.add(photo_id)
            item = dict(row)
            item["species"] = str(target.species)
            item["target_cell_order"] = int(target.target_cell_order)
            item["target_cell_center_latitude"] = float(target.cell_center_latitude)
            item["target_cell_center_longitude"] = float(target.cell_center_longitude)
            kept_rows.append(item)

        audit_rows.append({
            "species": str(target.species),
            "inat_taxon_id": taxon_id,
            "target_cell_order": int(target.target_cell_order),
            "cell_id": cell_id,
            "raw_results": int(len(raw_results)),
            "local_rejected": int(local_rejected),
            "wrong_taxon": int(wrong_taxon),
            "prior_excluded": int(prior_excluded),
            "fresh_duplicate": int(fresh_duplicate),
            "locally_eligible": int(len(locally_eligible)),
            "retained": int(len(selected)),
            "distinct_retained_observers": int(len({str(r["observer_id"]) for r in selected})),
            "request_error": request_error,
        })

    observations = pd.DataFrame(kept_rows)
    if len(observations):
        observations = observations.sort_values(
            ["inat_taxon_id", "target_cell_order", "selection_hash", "observation_id"],
            kind="mergesort",
        ).reset_index(drop=True)
        if observations["observation_id"].duplicated().any() or observations["photo_id"].duplicated().any():
            raise RuntimeError("fresh H7 retained IDs are not globally unique")
        if observations["observation_id"].isin(excluded_observations).any():
            raise RuntimeError("prior observation entered H7")
        if observations["photo_id"].isin(excluded_photos).any():
            raise RuntimeError("prior photo entered H7")

    target_audit = pd.DataFrame(audit_rows).sort_values(
        ["inat_taxon_id", "target_cell_order"], kind="mergesort"
    ).reset_index(drop=True)
    support = (
        target_audit.assign(full_cell=target_audit["retained"] == int(retained_cap))
        .groupby(["species", "inat_taxon_id"], observed=True)
        .agg(
            target_cells=("cell_id", "size"),
            full_cells=("full_cell", "sum"),
            retained_photos=("retained", "sum"),
            request_error_cells=("request_error", lambda s: int((s.astype(str) != "").sum())),
        )
        .reset_index()
    )
    support["passes_premeasurement_species_gate"] = support["full_cells"] >= int(required_full_cells_per_species)
    passing_species = int(support["passes_premeasurement_species_gate"].sum())
    gate_pass = passing_species >= int(required_species_for_gate)

    manifest = {
        "protocol": "random-photo-first-h7-balanced-itv-v1",
        "status": "h7_fresh_metadata_frozen_before_pixels",
        "query_attempts": int(len(target_audit)),
        "expected_query_attempts": int(expected_targets),
        "query_retries": 0,
        "request_error_targets": int((target_audit["request_error"].astype(str) != "").sum()),
        "retained_fresh_photos": int(len(observations)),
        "retained_unique_observation_ids": int(observations["observation_id"].nunique()) if len(observations) else 0,
        "retained_unique_photo_ids": int(observations["photo_id"].nunique()) if len(observations) else 0,
        "species": int(support.shape[0]),
        "target_cells": int(len(target_audit)),
        "premeasurement_gate": {
            "species_requiring_at_least_full_cells": int(required_full_cells_per_species),
            "full_cell_definition_retained_photos": int(retained_cap),
            "passing_species": passing_species,
            "required_species": int(required_species_for_gate),
            "pass": bool(gate_pass),
            "decision": "pixels_may_open" if gate_pass else "not_evaluable_h7_fresh_metadata_support_before_pixels",
        },
        "freshness_firewall": {
            "prior_exclusion_observation_ids": int(len(excluded_observations)),
            "prior_exclusion_photo_ids": int(len(excluded_photos)),
            "prior_observations_retained": 0,
            "prior_photos_retained": 0,
        },
        "candidate_image_pixels_opened": False,
        "colour_used_for_selection": False,
    }
    return H7MetadataFreeze(observations=observations, target_audit=target_audit, species_support=support, manifest=manifest)
