"""Metadata-only scale-out for the repeated FCP atlas cohorts.

This module never opens image pixels.  It converts a complete metadata-passing
species frame into all eight frozen panels and can run a bounded live-API
feasibility audit before the final dated export is acquired.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from .atlas_expansion import draw_disjoint_species_cohorts, validate_expansion_contract
from .image_first_atlas import (
    AtlasMetadataAdapter,
    _balanced_selection,
    _observation_query,
    _prepare_observation,
    _selection_qc,
    _species_query,
    validate_atlas_contract,
)
from .shared_transition_surface import EqualAreaGrid, build_edge_cell_geometry
from .spatial_graph import spherical_knn_edges


@dataclass(frozen=True)
class AtlasScaleoutFreeze:
    panels: tuple[dict[str, Any], ...]
    observations: tuple[dict[str, Any], ...]
    audit: dict[str, Any]


GEOMETRY_AMENDMENT_PROTOCOL = "jbi-atlas-scaleout-geometry-admission-amendment-v1"
GLOBAL_ID_AMENDMENT_PROTOCOL = "jbi-atlas-scaleout-global-id-amendment-v1"


def validate_geometry_admission_amendment(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != GEOMETRY_AMENDMENT_PROTOCOL:
        raise ValueError("unexpected scale-out geometry amendment")
    if (
        amendment.get("scaleout_colour_opened") is not False
        or amendment.get("candidate_image_pixels_opened") is not False
    ):
        raise ValueError("geometry admission amendment was not frozen pre-image")
    if amendment.get("primary_geometry_admission_scale_km") != 100:
        raise ValueError("scale-out primary geometry admission must remain 100 km")
    if amendment.get("recorded_sensitivity_scales_km") != [250, 500]:
        raise ValueError("scale-out geometry sensitivities changed")
    inherited = amendment.get("inherited_geometry_rules", {})
    if (
        inherited.get("knn_k") != 5
        or inherited.get("maximum_edge_km_by_scale") != [100, 250, 500]
        or inherited.get("minimum_edges_per_species_cell") != 2
        or inherited.get("minimum_retained_edges_per_species") != 100
        or inherited.get("minimum_detectable_cells_per_species") != 10
    ):
        raise ValueError("inherited scale-out geometry thresholds changed")


def validate_global_id_amendment(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != GLOBAL_ID_AMENDMENT_PROTOCOL:
        raise ValueError("unexpected scale-out global-ID amendment")
    trigger = amendment.get("trigger", {})
    rules = amendment.get("frozen_reconciliation", {})
    if (
        trigger.get("candidate_image_pixels_opened") is not False
        or trigger.get("continuous_colour_used") is not False
        or rules.get("ownership")
        != "a candidate may claim an observation ID and photo ID only after it passes both the metadata and primary 100-km geometry gates"
        or not str(rules.get("collision_rule", "")).startswith("before selection")
        or not str(rules.get("shortfall_rule", "")).startswith("if fewer than")
    ):
        raise ValueError("scale-out global-ID reconciliation changed")


def exclude_reserved_scaleout_rows(
    rows: Sequence[Mapping[str, Any]],
    reserved_observations: set[str],
    reserved_photos: set[str],
) -> tuple[list[dict[str, Any]], int]:
    """Deterministically remove metadata rows already claimed by earlier species."""

    retained: list[dict[str, Any]] = []
    removed = 0
    for raw in rows:
        observation_id = str(raw.get("observation_id", "")).strip()
        photo_id = str(raw.get("photo_id", "")).strip()
        if not observation_id or not photo_id:
            raise ValueError("scale-out candidates require observation and photo IDs")
        if observation_id in reserved_observations or photo_id in reserved_photos:
            removed += 1
            continue
        retained.append(dict(raw))
    return retained, removed


_OUTCOME_TOKENS = (
    "colour",
    "color",
    "pixel",
    "roi",
    "transition",
    "effect",
    "p_value",
)


def _reject_outcome_fields(row: Mapping[str, Any]) -> None:
    leaked = sorted(
        str(key)
        for key in row
        if any(token in str(key).casefold() for token in _OUTCOME_TOKENS)
    )
    if leaked:
        raise ValueError(f"image/colour outcome fields reached scale-out freeze: {leaked}")


def freeze_scaleout_panels(
    eligible_species: Sequence[Mapping[str, Any]],
    observations_by_taxon: Mapping[str, Sequence[Mapping[str, Any]]],
    expansion_contract: Mapping[str, Any],
    *,
    source_role: str,
) -> AtlasScaleoutFreeze:
    """Freeze 8 x 25 species and exactly 300 metadata rows per species."""

    validate_expansion_contract(expansion_contract)
    panels = draw_disjoint_species_cohorts(eligible_species, expansion_contract)
    target = int(expansion_contract["random_cohort_scaleout"]["observations_per_species"])
    output_rows: list[dict[str, Any]] = []
    seen_observations: set[str] = set()
    seen_photos: set[str] = set()
    for panel in panels:
        taxon_id = str(panel["taxon_id"])
        rows = [dict(row) for row in observations_by_taxon.get(taxon_id, ())]
        if len(rows) != target:
            raise ValueError(
                f"not_evaluable: taxon {taxon_id} requires exactly {target} frozen rows, "
                f"found {len(rows)}"
            )
        for row in rows:
            _reject_outcome_fields(row)
            if str(row.get("inat_taxon_id", "")) != taxon_id:
                raise ValueError(f"taxon {taxon_id} has a mismatched observation row")
            observation_id = str(row.get("observation_id", "")).strip()
            photo_id = str(row.get("photo_id", "")).strip()
            if not observation_id or not photo_id:
                raise ValueError("scale-out rows require observation_id and photo_id")
            if observation_id in seen_observations or photo_id in seen_photos:
                raise ValueError("scale-out observation and photo IDs must be globally unique")
            seen_observations.add(observation_id)
            seen_photos.add(photo_id)
            output_rows.append(
                {
                    "cohort_id": panel["cohort_id"],
                    "cohort_species_index": panel["cohort_species_index"],
                    **row,
                    "candidate_image_pixels_opened": False,
                }
            )

    expected_species = int(expansion_contract["random_cohort_scaleout"]["total_species"])
    expected_observations = int(
        expansion_contract["random_cohort_scaleout"]["total_observations"]
    )
    if len(panels) != expected_species or len(output_rows) != expected_observations:
        raise ValueError("scale-out dimensions do not match the frozen contract")
    audit = {
        "protocol": expansion_contract["protocol"],
        "status": "pass_metadata_only_scaleout_freeze",
        "source_role": source_role,
        "candidate_image_pixels_opened": False,
        "flower_roi_used": False,
        "continuous_colour_used": False,
        "eligible_species": len(eligible_species),
        "frozen_species": len(panels),
        "frozen_observations": len(output_rows),
        "cohorts": int(expansion_contract["random_cohort_scaleout"]["cohorts"]),
        "all_cohorts_required": True,
        "early_stopping": False,
    }
    return AtlasScaleoutFreeze(
        panels=tuple(dict(row) for row in panels),
        observations=tuple(output_rows),
        audit=audit,
    )


def qualify_scaleout_geometry(
    eligible_species: Sequence[Mapping[str, Any]],
    observations_by_taxon: Mapping[str, Sequence[Mapping[str, Any]]],
    atlas_contract: Mapping[str, Any],
    *,
    primary_scale_km: int = 100,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Apply all frozen geometry rules before the stable cohort permutation."""

    validate_atlas_contract(atlas_contract)
    geometry_contract = atlas_contract["geometry_only_scale_selection"]
    criteria = geometry_contract["passing_criteria"]
    candidates = tuple(geometry_contract["candidates"])
    if primary_scale_km not in {int(row["scale_km"]) for row in candidates}:
        raise ValueError("primary scale is absent from the frozen geometry candidates")
    passing: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    for raw in eligible_species:
        taxon_id = str(raw["taxon_id"])
        rows = list(observations_by_taxon.get(taxon_id, ()))
        if not rows:
            raise ValueError(f"metadata-eligible taxon {taxon_id} has no observations")
        latitude = np.asarray([float(row["latitude"]) for row in rows], dtype=float)
        longitude = np.asarray([float(row["longitude"]) for row in rows], dtype=float)
        edges, distance = spherical_knn_edges(
            latitude,
            longitude,
            k=int(geometry_contract["knn_k"]),
        )
        scale_results = []
        for candidate in candidates:
            scale = int(candidate["scale_km"])
            grid = EqualAreaGrid(
                n_lon=int(candidate["n_lon"]),
                n_sinlat=int(candidate["n_sinlat"]),
            )
            try:
                geometry = build_edge_cell_geometry(
                    latitude,
                    longitude,
                    edges,
                    distance,
                    grid=grid,
                    max_edge_km=scale,
                    min_edges_per_cell=int(
                        geometry_contract["minimum_edges_per_species_cell"]
                    ),
                )
                retained_edges = int(len(geometry.retained_edges))
                detectable_cells = int(np.count_nonzero(geometry.detectable))
                detectable_cell_ids = sorted(
                    int(cell)
                    for cell in np.flatnonzero(np.asarray(geometry.detectable, dtype=bool))
                )
            except ValueError:
                retained_edges = 0
                detectable_cells = 0
                detectable_cell_ids = []
            evaluable = (
                retained_edges >= int(criteria["minimum_retained_edges_per_species"])
                and detectable_cells >= int(criteria["minimum_detectable_cells_per_species"])
            )
            scale_results.append(
                {
                    "scale_km": scale,
                    "retained_edges": retained_edges,
                    "detectable_cells": detectable_cells,
                    "detectable_cell_ids": detectable_cell_ids,
                    "geometry_evaluable": evaluable,
                }
            )
        primary_pass = next(
            row["geometry_evaluable"]
            for row in scale_results
            if row["scale_km"] == primary_scale_km
        )
        audit.append(
            {
                "taxon_id": taxon_id,
                "species": str(raw["species"]),
                "status": "geometry_eligible" if primary_pass else "primary_geometry_failed",
                "primary_scale_km": primary_scale_km,
                "scale_results": scale_results,
            }
        )
        if primary_pass:
            passing.append(dict(raw))
    return passing, audit


def live_api_scaleout_feasibility(
    atlas_contract: Mapping[str, Any],
    expansion_contract: Mapping[str, Any],
    geometry_amendment: Mapping[str, Any],
    global_id_amendment: Mapping[str, Any],
    adapter: AtlasMetadataAdapter,
    *,
    candidate_species_pool_size: int = 500,
    maximum_candidates_per_species: int = 1000,
) -> AtlasScaleoutFreeze:
    """Audit whether a bounded live query can supply the frozen 200 x 300 design.

    The output is a feasibility freeze, not the final dated-source cohort.  Every
    candidate in the declared abundance-ranked pool is evaluated before the stable
    SHA-256 cohort draw; the function does not stop after finding 200 passes.
    """

    validate_atlas_contract(atlas_contract)
    validate_expansion_contract(expansion_contract)
    validate_geometry_admission_amendment(geometry_amendment)
    validate_global_id_amendment(global_id_amendment)
    if candidate_species_pool_size < 200:
        raise ValueError("live feasibility needs at least 200 candidate species")
    working = deepcopy(dict(atlas_contract))
    working["metadata_source"]["candidate_species_pool_size"] = int(
        candidate_species_pool_size
    )
    working["admission"]["maximum_candidates_per_species"] = int(
        maximum_candidates_per_species
    )
    working["admission"]["sample_size_tiers_descending"] = [
        int(expansion_contract["random_cohort_scaleout"]["observations_per_species"])
    ]

    excluded_names = {
        str(value).casefold() for value in working["admission"]["excluded_frozen_species"]
    }
    raw_counts = list(adapter.species_counts(_species_query(working)))
    candidates: list[dict[str, Any]] = []
    for record in raw_counts:
        taxon = record.get("taxon") or {}
        try:
            count = int(record["count"])
            taxon_id = int(taxon["id"])
            genus_id = int(taxon["parent_id"])
        except (KeyError, TypeError, ValueError):
            continue
        name = str(taxon.get("name") or "").strip()
        if (
            str(taxon.get("rank") or "").casefold() != "species"
            or taxon.get("is_active") is not True
            or not name
            or name.casefold() in excluded_names
        ):
            continue
        candidates.append(
            {
                "count": count,
                "taxon": {
                    "id": taxon_id,
                    "name": name,
                    "rank": "species",
                    "parent_id": genus_id,
                },
            }
        )
    candidates.sort(
        key=lambda row: (
            -int(row["count"]),
            str(row["taxon"]["name"]),
            int(row["taxon"]["id"]),
        )
    )

    eligible: list[dict[str, Any]] = []
    observations_by_taxon: dict[str, list[dict[str, Any]]] = {}
    species_audit: list[dict[str, Any]] = []
    reserved_observations: set[str] = set()
    reserved_photos: set[str] = set()
    metadata_eligible_species = 0
    for rank, record in enumerate(candidates, start=1):
        taxon = record["taxon"]
        raw_observations = adapter.observations(
            int(taxon["id"]), _observation_query(working)
        )
        prepared_before_reconciliation = [
            row
            for observation in raw_observations
            if (row := _prepare_observation(observation, taxon, working)) is not None
        ]
        prepared, collisions_removed = exclude_reserved_scaleout_rows(
            prepared_before_reconciliation,
            reserved_observations,
            reserved_photos,
        )
        selected = _balanced_selection(prepared, working)
        qc = _selection_qc(selected, len(prepared), working)
        taxon_id = str(taxon["id"])
        audit_row = {
            "metadata_rank": rank,
            "taxon_id": taxon_id,
            "species": str(taxon["name"]),
            "genus": f"inat-genus-{taxon['parent_id']}",
            "flowering_annotated_observation_count": int(record["count"]),
            "pre_reconciliation_candidate_count": len(prepared_before_reconciliation),
            "global_identity_collisions_removed": collisions_removed,
            "post_reconciliation_candidate_count": len(prepared),
            "status": "metadata_eligible" if qc["gate_pass"] else "metadata_gate_failed",
            **qc,
        }
        if not qc["gate_pass"]:
            species_audit.append(audit_row)
            continue
        metadata_eligible_species += 1
        candidate = {
            "taxon_id": taxon_id,
            "species": str(taxon["name"]),
            "genus": f"inat-genus-{taxon['parent_id']}",
        }
        geometry_eligible, geometry_audit = qualify_scaleout_geometry(
            [candidate],
            {taxon_id: selected},
            atlas_contract,
            primary_scale_km=int(expansion_contract["spatial_design"]["primary_scale_km"]),
        )
        geometry = geometry_audit[0]
        audit_row["geometry_status"] = geometry["status"]
        audit_row["geometry_scale_results"] = geometry["scale_results"]
        audit_row["status"] = geometry["status"]
        species_audit.append(audit_row)
        if not geometry_eligible:
            continue
        eligible.append(candidate)
        observations_by_taxon[taxon_id] = selected
        for row in selected:
            observation_id = str(row["observation_id"])
            photo_id = str(row["photo_id"])
            if observation_id in reserved_observations or photo_id in reserved_photos:
                raise ValueError("global-ID reconciliation failed before panel draw")
            reserved_observations.add(observation_id)
            reserved_photos.add(photo_id)
    try:
        frozen = freeze_scaleout_panels(
            eligible,
            observations_by_taxon,
            expansion_contract,
            source_role="live iNaturalist API metadata feasibility only",
        )
    except ValueError as exc:
        if not str(exc).startswith("not_evaluable:"):
            raise
        return AtlasScaleoutFreeze(
            panels=(),
            observations=(),
            audit={
                "protocol": expansion_contract["protocol"],
                "status": "not_evaluable_live_api_scaleout_shortfall",
                "reason": str(exc),
                "source_role": "live iNaturalist API metadata feasibility only",
                "candidate_image_pixels_opened": False,
                "flower_roi_used": False,
                "continuous_colour_used": False,
                "species_count_records_received": len(raw_counts),
                "species_candidates_after_static_filters": len(candidates),
                "eligible_species": metadata_eligible_species,
                "geometry_eligible_species": len(eligible),
                "globally_reserved_observations": len(reserved_observations),
                "globally_reserved_photos": len(reserved_photos),
                "global_id_reconciliation_protocol": global_id_amendment["protocol"],
                "species_results": species_audit,
            },
        )

    audit = {
        **frozen.audit,
        "status": "pass_live_api_scaleout_feasibility",
        "species_count_records_received": len(raw_counts),
        "species_candidates_after_static_filters": len(candidates),
        "species_results": species_audit,
        "metadata_eligible_species": metadata_eligible_species,
        "geometry_eligible_species": len(eligible),
        "globally_reserved_observations": len(reserved_observations),
        "globally_reserved_photos": len(reserved_photos),
        "global_id_reconciliation_protocol": global_id_amendment["protocol"],
        "geometry_admission_preceded_cohort_permutation": True,
        "geometry_admission_protocol": geometry_amendment["protocol"],
        "final_source_still_required": True,
    }
    return AtlasScaleoutFreeze(frozen.panels, frozen.observations, audit)
