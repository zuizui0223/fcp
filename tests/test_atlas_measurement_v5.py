from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_acquisition_v5 import validate_v5_acquisition_firewall
from fcp_pipeline.atlas_measurement_v5 import (
    build_v5_measurement_firewall,
    validate_measurement_execution_contract,
    validate_preimage_gates,
)


CONTRACT_PATH = Path("docs/supporting/jbi_atlas_measurement_execution_contract_v5.json")
INFERENCE_PATH = Path("docs/supporting/jbi_image_first_atlas_inference_contract_v5.json")
SOURCE_PROTOCOL = "jbi-atlas-source-role-amendment-v5"
SOURCE_PASS = "pass_frozen_selection_dated_provenance_v5"
SOURCE_ROWS = "scaleout_observation_manifest.csv"


def _contracts():
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    inference = json.loads(INFERENCE_PATH.read_text(encoding="utf-8"))
    return contract, inference


def _passing_gates():
    source = {
        "protocol": SOURCE_PROTOCOL,
        "status": SOURCE_PASS,
        "candidate_image_pixels_opened": False,
        "selected_species": 200,
        "selected_photo_assets": 60000,
        "frozen_observations": 60000,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "current_live_state_used_for_authorization": False,
        "repeat_35gb_stream_used_for_v5_authorization": False,
        "parents": {"scaleout_observation_manifest_sha256": "rows-sha"},
    }
    manifest = {
        "protocol": SOURCE_PROTOCOL,
        "status": SOURCE_PASS,
        "candidate_image_pixels_opened": False,
        "replacement_permitted": False,
        "files": {"source_v5_result.json": "source-result-sha"},
    }
    environment = {
        "status": "pass_precolour_environmental_coverage",
        "coverage_gate_status": "pass_precolour_environmental_coverage",
        "source_stage": "final-source-v5",
        "final_source_v5_required": False,
        "scaleout_colour_opened": False,
        "image_acquisition_authorized": False,
    }
    roi = {
        "status": "pass_roi_v4_locked_test",
        "scaleout_candidate_pixels_opened": False,
        "scaleout_candidate_pixels_permitted": True,
    }
    shared = {
        "status": "pass",
        "scope": {
            "method_gate_only": True,
            "biological_support_claimed": False,
            "candidate_image_pixels_opened": False,
        },
    }
    return source, manifest, environment, roi, shared


def test_v5_measurement_contract_keeps_terminal_inference_on_v5() -> None:
    contract, inference = _contracts()
    validate_measurement_execution_contract(contract, inference)
    assert contract["v5_inference_firewall"]["only_inference_contract"].endswith(
        "jbi_image_first_atlas_inference_contract_v5.json"
    )
    assert contract["pixel_opening_gates"]["dated_source"]["required_protocol"] == SOURCE_PROTOCOL
    assert contract["pixel_opening_gates"]["dated_source"]["required_status"] == SOURCE_PASS
    assert contract["pixel_opening_gates"]["environmental_coverage"]["required_source_stage"] == "final-source-v5"
    assert contract["v5_inference_firewall"][
        "superseded_v3_ordered_inference_must_not_authorize_or_classify_terminal_results"
    ] is True


def test_all_preimage_gates_are_jointly_required() -> None:
    contract, _ = _contracts()
    source, manifest, environment, roi, shared = _passing_gates()
    validate_preimage_gates(
        dated_reconciliation=source,
        dated_manifest=manifest,
        observation_manifest_name=SOURCE_ROWS,
        observation_manifest_sha256="rows-sha",
        environmental_coverage=environment,
        roi_locked_result=roi,
        shared_qualification_result=shared,
        contract=contract,
    )

    broken = dict(source)
    broken["status"] = "not_evaluable_source_role_v5"
    with pytest.raises(RuntimeError, match="source role v5"):
        validate_preimage_gates(
            dated_reconciliation=broken,
            dated_manifest=manifest,
            observation_manifest_name=SOURCE_ROWS,
            observation_manifest_sha256="rows-sha",
            environmental_coverage=environment,
            roi_locked_result=roi,
            shared_qualification_result=shared,
            contract=contract,
        )

    wrong_protocol = dict(source)
    wrong_protocol["protocol"] = "jbi-atlas-dated-source-uuid-bucket-amendment-v4"
    with pytest.raises(RuntimeError, match="source role v5"):
        validate_preimage_gates(
            dated_reconciliation=wrong_protocol,
            dated_manifest=manifest,
            observation_manifest_name=SOURCE_ROWS,
            observation_manifest_sha256="rows-sha",
            environmental_coverage=environment,
            roi_locked_result=roi,
            shared_qualification_result=shared,
            contract=contract,
        )

    wrong_row_hash = json.loads(json.dumps(source))
    wrong_row_hash["parents"]["scaleout_observation_manifest_sha256"] = "wrong"
    with pytest.raises(RuntimeError, match="source role v5"):
        validate_preimage_gates(
            dated_reconciliation=wrong_row_hash,
            dated_manifest=manifest,
            observation_manifest_name=SOURCE_ROWS,
            observation_manifest_sha256="rows-sha",
            environmental_coverage=environment,
            roi_locked_result=roi,
            shared_qualification_result=shared,
            contract=contract,
        )

    broken_environment = dict(environment)
    broken_environment["source_stage"] = "final-dated-source"
    with pytest.raises(RuntimeError, match="environmental coverage"):
        validate_preimage_gates(
            dated_reconciliation=source,
            dated_manifest=manifest,
            observation_manifest_name=SOURCE_ROWS,
            observation_manifest_sha256="rows-sha",
            environmental_coverage=broken_environment,
            roi_locked_result=roi,
            shared_qualification_result=shared,
            contract=contract,
        )

    broken_shared = json.loads(json.dumps(shared))
    broken_shared["scope"]["biological_support_claimed"] = True
    with pytest.raises(RuntimeError, match="shared-transition"):
        validate_preimage_gates(
            dated_reconciliation=source,
            dated_manifest=manifest,
            observation_manifest_name=SOURCE_ROWS,
            observation_manifest_sha256="rows-sha",
            environmental_coverage=environment,
            roi_locked_result=roi,
            shared_qualification_result=broken_shared,
            contract=contract,
        )


def test_v5_worker_packet_is_location_blind_and_uses_unchanged_salt() -> None:
    contract, _ = _contracts()
    rows = [
        {
            "cohort_id": "C01",
            "species": "Plantus alpha",
            "inat_taxon_id": "11",
            "observation_id": "101",
            "photo_id": "201",
            "photo_url_large": "https://example.invalid/201.jpg",
            "photo_license": "cc-by",
            "latitude": "1.0",
            "longitude": "2.0",
            "observed_month": "6",
            "local_solar_quarter": "2",
            "observer_id": "301",
        },
        {
            "cohort_id": "C02",
            "species": "Plantus beta",
            "inat_taxon_id": "12",
            "observation_id": "102",
            "photo_id": "202",
            "photo_url_large": "https://example.invalid/202.jpg",
            "photo_license": "cc-by-sa",
            "latitude": "3.0",
            "longitude": "4.0",
            "observed_month": "7",
            "local_solar_quarter": "3",
            "observer_id": "302",
        },
    ]
    split = build_v5_measurement_firewall(rows, contract)
    allowed = set(contract["location_blind_measurement"]["measurement_worker_allowed_fields"])
    assert len(split["measurement_manifest"]) == 2
    assert all(set(row) == allowed for row in split["measurement_manifest"])
    assert all("latitude" not in row and "longitude" not in row for row in split["measurement_manifest"])
    assert all("photo_url_large" not in row and "species" not in row for row in split["measurement_manifest"])
    assert len(split["sealed_coordinate_key"]) == 2
    assert split["measurement_manifest"][0]["measurement_id"].startswith("FCPM-")
    assert split["measurement_manifest"][0]["species_blind_id"].startswith("FCPS-")


def test_v5_acquisition_rejects_legacy_firewall_and_requires_all_gate_hashes() -> None:
    contract, inference = _contracts()
    firewall = {
        "status": "pass_scaleout_measurement_firewall_v5",
        "protocol": contract["protocol"],
        "inference_version": inference["version"],
        "frozen_measurements": 60000,
        "candidate_image_pixels_opened": False,
        "terminal_scaleout_colour_measured": False,
        "coordinate_key_opened_by_measurement_worker": False,
        "superseded_v3_ordered_inference_used": False,
        "sealed_keys": {"acquisition_coordinate_key.csv": "key-sha"},
        "preimage_gate_sha256": {
            "source_v5_manifest": "a",
            "source_v5_result": "b",
            "environmental_coverage": "c",
            "roi_locked_result": "d",
            "shared_transition_qualification": "e",
        },
    }
    validate_v5_acquisition_firewall(
        firewall,
        key_name="acquisition_coordinate_key.csv",
        key_sha256="key-sha",
        contract=contract,
        inference_v5=inference,
    )

    legacy = dict(firewall)
    legacy["status"] = "pass_scaleout_measurement_firewall"
    with pytest.raises(RuntimeError, match="does not authorize"):
        validate_v5_acquisition_firewall(
            legacy,
            key_name="acquisition_coordinate_key.csv",
            key_sha256="key-sha",
            contract=contract,
            inference_v5=inference,
        )

    incomplete = json.loads(json.dumps(firewall))
    del incomplete["preimage_gate_sha256"]["environmental_coverage"]
    with pytest.raises(RuntimeError, match="gate hashes"):
        validate_v5_acquisition_firewall(
            incomplete,
            key_name="acquisition_coordinate_key.csv",
            key_sha256="key-sha",
            contract=contract,
            inference_v5=inference,
        )
