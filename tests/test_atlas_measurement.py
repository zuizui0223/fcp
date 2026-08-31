from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_measurement import (
    build_measurement_firewall,
    evaluate_scaleout_measurement_gate,
    validate_inference_contract,
)


CONTRACT = Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json")


def test_inference_contract_keeps_stopped_geographic_branch_closed() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_inference_contract(contract)
    contract["ordered_inference_v3"]["branches"][0]["real_colour_test_permitted"] = True
    with pytest.raises(ValueError, match="geographic branch"):
        validate_inference_contract(contract)


def test_coordinate_firewall_strips_location_taxon_date_and_observer() -> None:
    rows = [
        {
            "cohort_id": "C01",
            "species": "Example flower",
            "inat_taxon_id": 9,
            "observation_id": 10,
            "photo_id": 11,
            "photo_url_large": "https://example.test/11/large.jpg",
            "photo_license": "cc-by",
            "attribution": "observer",
            "latitude": 35.0,
            "longitude": 135.0,
            "observed_month": 4,
            "local_solar_quarter": 2,
            "observer_id": 12,
        }
    ]
    split = build_measurement_firewall(rows)
    measurement = split["measurement_manifest"][0]
    assert not {
        "species",
        "inat_taxon_id",
        "latitude",
        "longitude",
        "observed_month",
        "observer_id",
        "attribution",
        "photo_id",
        "photo_url_large",
    }.intersection(measurement)
    assert measurement["image_filename"].startswith("FCPM-")
    assert measurement["image_filename"].endswith(".jpg")
    assert split["sealed_species_key"][0].keys() == {
        "measurement_id",
        "species_blind_id",
        "cohort_id",
    }
    assert split["sealed_coordinate_key"][0]["photo_url_large"].startswith("https://")


def test_coordinate_firewall_rejects_precomputed_pixel_outcomes() -> None:
    row = {
        "cohort_id": "C01",
        "species": "Example flower",
        "inat_taxon_id": 9,
        "observation_id": 10,
        "photo_id": 11,
        "photo_url_large": "https://example.test/11/large.jpg",
        "photo_license": "cc-by",
        "latitude": 35.0,
        "longitude": 135.0,
        "observed_month": 4,
        "local_solar_quarter": 2,
        "observer_id": 12,
        "flower_roi_pixels": 100,
    }
    with pytest.raises(ValueError, match="image outcomes"):
        build_measurement_firewall([row])


def test_measurement_gate_passes_only_complete_eight_cohort_denominator() -> None:
    gate = json.loads(CONTRACT.read_text(encoding="utf-8"))[
        "scaleout_measurement_gate"
    ]
    species_key = []
    results = []
    for cohort in range(1, 9):
        for species in range(25):
            species_id = f"S-{cohort}-{species}"
            for observation in range(300):
                measurement_id = f"M-{cohort}-{species}-{observation}"
                species_key.append(
                    {
                        "measurement_id": measurement_id,
                        "species_blind_id": species_id,
                        "cohort_id": f"C{cohort:02d}",
                    }
                )
                results.append(
                    {
                        "measurement_id": measurement_id,
                        "automated_colour_state_status": (
                            "automated_colour_state_admitted"
                            if observation < 210
                            else "automated_colour_state_not_evaluable"
                        ),
                        "background_features_available": observation < 147,
                    }
                )
    decision = evaluate_scaleout_measurement_gate(results, species_key, gate)
    assert decision["status"] == "pass_scaleout_measurement_completeness"
    assert decision["evaluable_species"] == 200
    assert decision["coordinates_opened"] is False

    # Six failures make C01 fall below the frozen 20-species cohort minimum.
    for species in range(6):
        results[species * 300]["automated_colour_state_status"] = (
            "image_acquisition_failed"
        )
    decision = evaluate_scaleout_measurement_gate(results, species_key, gate)
    assert decision["status"] == "not_evaluable_scaleout_measurement_completeness"
    assert decision["coordinate_join_permitted"] is False
