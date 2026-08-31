from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_measurement import (
    build_measurement_firewall,
    evaluate_scaleout_measurement_gate,
    measurement_shard,
    select_measurement_shard,
    validate_inference_contract,
    validate_measurement_result_rows,
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


def test_worker_shards_are_disjoint_complete_and_reject_leaks() -> None:
    rows = [
        {
            "measurement_id": f"FCPM-{index}",
            "species_blind_id": f"FCPS-{index % 3}",
            "image_filename": f"FCPM-{index}.jpg",
            "photo_license": "cc-by",
        }
        for index in range(50)
    ]
    shards = [
        select_measurement_shard(rows, shard_index=index, shard_count=7)
        for index in range(7)
    ]
    assigned = [row["measurement_id"] for shard in shards for row in shard]
    assert sorted(assigned) == sorted(row["measurement_id"] for row in rows)
    assert len(set(assigned)) == len(rows)
    assert measurement_shard("FCPM-1", 7) == measurement_shard("FCPM-1", 7)

    leaked = [dict(rows[0], latitude=35.0)]
    with pytest.raises(ValueError, match="fields changed or leaked"):
        select_measurement_shard(leaked, shard_index=0, shard_count=1)


def test_measurement_results_are_location_free_terminal_rows() -> None:
    valid = [
        {
            "measurement_id": "FCPM-1",
            "species_blind_id": "FCPS-1",
            "automated_colour_state_status": "automated_colour_state_admitted",
            "background_features_available": True,
            "flower_L_mean": 50.0,
        }
    ]
    assert validate_measurement_result_rows(valid) == valid
    with pytest.raises(ValueError, match="protected fields"):
        validate_measurement_result_rows([dict(valid[0], latitude=35.0)])
    with pytest.raises(ValueError, match="cannot have background"):
        validate_measurement_result_rows(
            [
                {
                    "measurement_id": "FCPM-1",
                    "species_blind_id": "FCPS-1",
                    "automated_colour_state_status": "image_acquisition_failed",
                    "background_features_available": True,
                }
            ]
        )


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
