#!/usr/bin/env python3
"""Audit equal-area coverage retained by bounded candidate species before pixels."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_capacity_handoff import resolve_capacity_handoff
from fcp_pipeline.global_species_budget_coverage import (
    cell_count_correlation,
    occupied_cell_retention,
    summarize_cell_coverage,
    unique_taxon_cell_links,
)

ROOT = Path(__file__).resolve().parents[2]
AUTH = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_authorization_v1.json"
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_contract_v1.json"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
CANDIDATE_AUDIT = ROOT / "data/frozen/global_monte_carlo_candidate_species_audit_v1.csv"
OUT = ROOT / "docs/supporting/global_monte_carlo_candidate_species_budget_coverage_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def block(summary) -> dict[str, object]:
    return {
        "taxon_count": int(summary.taxon_count),
        "occupied_equal_area_cells": int(summary.occupied_cells),
        "unique_taxon_cell_links": int(summary.total_taxon_cell_links),
        "cell_species_gini": float(summary.cell_species_gini),
    }


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite coverage audit: {OUT}")
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    handoff = resolve_capacity_handoff(
        ROOT,
        AUTH,
        minimum_species=int(contract["premeasurement_gate"]["minimum_full_target_species"]),
    )
    audit = pd.read_csv(CANDIDATE_AUDIT)
    if len(audit) != handoff.selected_species or audit["inat_taxon_id"].nunique() != handoff.selected_species:
        raise RuntimeError("candidate audit does not equal the bounded queried species set")

    v1 = pd.read_csv(V1, usecols=["inat_taxon_id", "cell_id"])
    v2 = pd.read_csv(V2, usecols=["inat_taxon_id", "cell_id"])
    links = unique_taxon_cell_links(v1, v2)

    capacity_full = pd.read_csv(handoff.selected_path, usecols=["inat_taxon_id"])
    capacity_ids = set(capacity_full["inat_taxon_id"].astype(int))
    queried_ids = set(audit["inat_taxon_id"].astype(int))
    full_target_ids = set(audit.loc[audit["full_target"].astype(bool), "inat_taxon_id"].astype(int))
    if not queried_ids.issubset(capacity_ids) or not full_target_ids.issubset(queried_ids):
        raise RuntimeError("candidate coverage taxon hierarchy drifted")

    cap = summarize_cell_coverage(links, capacity_ids)
    queried = summarize_cell_coverage(links, queried_ids)
    retained = summarize_cell_coverage(links, full_target_ids)
    result = {
        "protocol": "global-monte-carlo-candidate-species-budget-coverage-v1",
        "status": "complete_metadata_only_equal_area_candidate_budget_coverage_audit",
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "capacity_source": handoff.source,
        "selected_raw_photo_target": int(handoff.target),
        "capacity_selected": block(cap),
        "bounded_candidate_queried": block(queried),
        "fresh_full_target_retained": block(retained),
        "bounded_candidate_occupied_cell_retention_vs_capacity": occupied_cell_retention(cap, queried),
        "fresh_full_target_occupied_cell_retention_vs_capacity": occupied_cell_retention(cap, retained),
        "bounded_candidate_cell_species_count_pearson_vs_capacity": cell_count_correlation(cap, queried),
        "fresh_full_target_cell_species_count_pearson_vs_capacity": cell_count_correlation(cap, retained),
        "interpretation": "Descriptive pre-pixel coverage audit only. It measures whether deterministic bounded species stages preserve the equal-area discovery footprint of the full capacity-selected frame; it does not select or rescue flower-colour outcomes.",
        "lineage": {
            "capacity_manifest_sha256": sha256_file(handoff.manifest_path),
            "capacity_selected_species_sha256": sha256_file(handoff.selected_path),
            "candidate_species_audit_sha256": sha256_file(CANDIDATE_AUDIT),
            "v1_discovery_index_sha256": sha256_file(V1),
            "v2_discovery_index_sha256": sha256_file(V2),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
