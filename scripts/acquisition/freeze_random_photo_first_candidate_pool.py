#!/usr/bin/env python3
"""Freeze the fresh species-unfixed metadata pool for the random photo-first atlas.

This command opens iNaturalist observation metadata only. It must run before any
candidate image pixels are downloaded or measured. The exact returned IDs are
frozen once; failed licence/geometry rows are not replaced by a second random
query.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from fcp_pipeline.photo_first_atlas import species_capped_sampling_capacity
from fcp_pipeline.random_photo_pool import (
    InaturalistObservationClient,
    freeze_random_photo_candidate_pool,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid


CONTRACT = Path("docs/supporting/random_photo_first_candidate_pool_contract_v1.json")
CANDIDATE_CSV = Path("data/frozen/random_photo_first_candidate_pool_v1.csv")
CELL_AUDIT_CSV = Path("data/frozen/random_photo_first_candidate_pool_cell_audit_v1.csv")
MANIFEST_JSON = Path("docs/supporting/random_photo_first_candidate_pool_manifest_v1.json")


def load_contract() -> dict:
    payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if payload["protocol"] != "random-photo-first-candidate-pool-v1":
        raise RuntimeError("unexpected candidate-pool contract")
    if payload["status"] != "frozen_before_any_fresh_candidate_api_query_or_image_pixel":
        raise RuntimeError("candidate-pool contract is not in the pre-query frozen state")
    if payload["species_policy"]["fixed_or_targeted_species_list"] is not False:
        raise RuntimeError("species targeting is forbidden")
    if payload["legacy_firewall"]["pr21_terminal_60000_photo_records_used"] is not False:
        raise RuntimeError("PR21 terminal records are forbidden")
    return payload


def main() -> None:
    contract = load_contract()
    source = contract["source"]
    sampling = contract["geographic_sampling"]
    grid_spec = sampling["grid"]
    grid = EqualAreaGrid(
        n_lon=int(grid_spec["n_lon"]),
        n_sinlat=int(grid_spec["n_sinlat"]),
    )
    if grid.n_cells != int(grid_spec["cells"]):
        raise RuntimeError("candidate-pool grid cell count drift")

    client = InaturalistObservationClient()
    frozen = freeze_random_photo_candidate_pool(
        client=client,
        grid=grid,
        per_cell_cap=int(sampling["per_cell_random_candidate_cap"]),
        taxon_id=int(source["angiosperm_taxon_id"]),
        flowering_term_id=int(source["flowering_annotation_term_id"]),
        flowering_term_value_id=int(source["flowering_annotation_value_id"]),
        maximum_positional_accuracy_m=int(source["maximum_positional_accuracy_m"]),
        allowed_photo_licenses=tuple(source["allowed_photo_licenses"]),
    )

    h1_target = 10_000
    h1_species_cap = 2
    capped_capacity = species_capped_sampling_capacity(
        frozen.observations,
        species_cap_per_cell=h1_species_cap,
    )
    manifest = dict(frozen.manifest)
    manifest["source_commit"] = os.environ.get("GITHUB_SHA", "local")
    manifest["premmeasurement_h1_gate"] = {
        "fixed_photos_per_replicate": h1_target,
        "species_cap_per_cell_per_replicate": h1_species_cap,
        "candidate_species_capped_capacity": int(capped_capacity),
        "candidate_pool_can_reach_fixed_replicate_size": bool(capped_capacity >= h1_target),
        "pixels_may_open": bool(capped_capacity >= h1_target),
        "reason_if_false": (
            None
            if capped_capacity >= h1_target
            else "not_evaluable_candidate_sampling_capacity_before_pixels"
        ),
    }

    CANDIDATE_CSV.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_JSON.parent.mkdir(parents=True, exist_ok=True)
    frozen.observations.to_csv(CANDIDATE_CSV, index=False, lineterminator="\n")
    frozen.cell_audit.to_csv(CELL_AUDIT_CSV, index=False, lineterminator="\n")
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(manifest, indent=2))
    if capped_capacity < h1_target:
        print(
            "Candidate pool froze successfully but H1 is not evaluable at the fixed "
            "10,000-photo replicate size. Candidate image pixels must remain unopened."
        )
    else:
        print(
            "Candidate pool froze successfully and passes the metadata-only capacity gate. "
            "A separately frozen image-measurement contract may now open candidate pixels."
        )


if __name__ == "__main__":
    main()
