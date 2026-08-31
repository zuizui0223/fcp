"""Metadata-only scale-out for the repeated FCP atlas cohorts.

This module never opens image pixels.  It converts a complete metadata-passing
species frame into all eight frozen panels and can run a bounded live-API
feasibility audit before the final dated export is acquired.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

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


@dataclass(frozen=True)
class AtlasScaleoutFreeze:
    panels: tuple[dict[str, Any], ...]
    observations: tuple[dict[str, Any], ...]
    audit: dict[str, Any]


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


def live_api_scaleout_feasibility(
    atlas_contract: Mapping[str, Any],
    expansion_contract: Mapping[str, Any],
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
    for rank, record in enumerate(candidates, start=1):
        taxon = record["taxon"]
        raw_observations = adapter.observations(
            int(taxon["id"]), _observation_query(working)
        )
        prepared = [
            row
            for observation in raw_observations
            if (row := _prepare_observation(observation, taxon, working)) is not None
        ]
        selected = _balanced_selection(prepared, working)
        qc = _selection_qc(selected, len(prepared), working)
        taxon_id = str(taxon["id"])
        species_audit.append(
            {
                "metadata_rank": rank,
                "taxon_id": taxon_id,
                "species": str(taxon["name"]),
                "genus": f"inat-genus-{taxon['parent_id']}",
                "flowering_annotated_observation_count": int(record["count"]),
                "status": "eligible" if qc["gate_pass"] else "metadata_gate_failed",
                **qc,
            }
        )
        if not qc["gate_pass"]:
            continue
        eligible.append(
            {
                "taxon_id": taxon_id,
                "species": str(taxon["name"]),
                "genus": f"inat-genus-{taxon['parent_id']}",
            }
        )
        observations_by_taxon[taxon_id] = selected

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
                "eligible_species": len(eligible),
                "species_results": species_audit,
            },
        )

    audit = {
        **frozen.audit,
        "status": "pass_live_api_scaleout_feasibility",
        "species_count_records_received": len(raw_counts),
        "species_candidates_after_static_filters": len(candidates),
        "species_results": species_audit,
        "final_source_still_required": True,
    }
    return AtlasScaleoutFreeze(frozen.panels, frozen.observations, audit)
