from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.data.run_jbi_atlas_precolour_environmental_coverage import build_evidence


CONTRACT = json.loads(
    Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json").read_text(
        encoding="utf-8"
    )
)


def inputs() -> tuple[dict, list[dict], dict[int, list[dict]]]:
    panels = []
    species_results = []
    for index in range(200):
        taxon = str(index + 1)
        panels.append(
            {
                "taxon_id": taxon,
                "cohort_id": f"C{index // 25 + 1:02d}",
                "candidate_image_pixels_opened": False,
            }
        )
        species_results.append(
            {
                "taxon_id": taxon,
                "status": "geometry_eligible",
                "geometry_scale_results": [
                    {
                        "scale_km": scale,
                        "detectable_cell_ids": [1, 2, 3, 4],
                    }
                    for scale in (100, 250, 500)
                ],
            }
        )
    feasibility = {
        "status": "pass_live_api_scaleout_feasibility",
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "species_results": species_results,
    }
    boundaries = {
        scale: [
            {
                "cell_id": cell,
                "macroclimate_boundary": 0.1,
                "land_cover_boundary": 0.2,
                "ecoregion_boundary": 0.0,
            }
            for cell in (1, 2, 3, 4)
        ]
        for scale in (100, 250, 500)
    }
    return feasibility, panels, boundaries


def dated_source_pass() -> dict:
    return {
        "status": "pass_dated_source_m2m_scaleout_freeze",
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": 200,
        "selected_photo_assets": 60000,
        "frozen_observations": 60000,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
    }


def test_live_feasibility_cannot_masquerade_as_final_coverage() -> None:
    feasibility, panels, boundaries = inputs()
    result = build_evidence(
        feasibility=feasibility,
        panels=panels,
        boundary_rows_by_scale=boundaries,
        environment_contract=CONTRACT,
        source_stage="live-feasibility",
    )
    assert result["status"] == "pass_live_api_precolour_environmental_coverage_feasibility"
    assert result["coverage_gate_status"] == "pass_precolour_environmental_coverage"
    assert result["final_dated_source_required"] is True
    assert result["image_acquisition_authorized"] is False


def test_final_dated_source_emits_the_protected_join_status() -> None:
    feasibility, panels, boundaries = inputs()
    result = build_evidence(
        feasibility=feasibility,
        panels=panels,
        boundary_rows_by_scale=boundaries,
        environment_contract=CONTRACT,
        source_stage="final-dated-source",
        dated_source_reconciliation=dated_source_pass(),
    )
    assert result["status"] == "pass_precolour_environmental_coverage"
    assert result["final_dated_source_required"] is False


def test_final_coverage_requires_a_separate_complete_dated_source_gate() -> None:
    feasibility, panels, boundaries = inputs()
    with pytest.raises(ValueError, match="requires dated-source"):
        build_evidence(
            feasibility=feasibility,
            panels=panels,
            boundary_rows_by_scale=boundaries,
            environment_contract=CONTRACT,
            source_stage="final-dated-source",
        )

    dated = dated_source_pass()
    dated["frozen_observations"] = 59999
    with pytest.raises(ValueError, match="did not pass unchanged"):
        build_evidence(
            feasibility=feasibility,
            panels=panels,
            boundary_rows_by_scale=boundaries,
            environment_contract=CONTRACT,
            source_stage="final-dated-source",
            dated_source_reconciliation=dated,
        )


def test_live_coverage_rejects_dated_source_outcomes() -> None:
    feasibility, panels, boundaries = inputs()
    with pytest.raises(ValueError, match="must not consume"):
        build_evidence(
            feasibility=feasibility,
            panels=panels,
            boundary_rows_by_scale=boundaries,
            environment_contract=CONTRACT,
            source_stage="live-feasibility",
            dated_source_reconciliation=dated_source_pass(),
        )


def test_panel_or_geometry_drift_fails_closed() -> None:
    feasibility, panels, boundaries = inputs()
    panels.pop()
    with pytest.raises(ValueError, match="exactly 200"):
        build_evidence(
            feasibility=feasibility,
            panels=panels,
            boundary_rows_by_scale=boundaries,
            environment_contract=CONTRACT,
            source_stage="live-feasibility",
        )

    feasibility, panels, boundaries = inputs()
    feasibility["species_results"][0]["status"] = "primary_geometry_failed"
    with pytest.raises(ValueError, match="not geometry eligible"):
        build_evidence(
            feasibility=feasibility,
            panels=panels,
            boundary_rows_by_scale=boundaries,
            environment_contract=CONTRACT,
            source_stage="live-feasibility",
        )
