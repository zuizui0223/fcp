from __future__ import annotations

import pytest

from scripts.data.acquire_jbi_atlas_blinded_images import (
    validate_firewall_for_acquisition,
)
from scripts.data.build_jbi_atlas_measurement_firewall import (
    validate_preimage_firewall_gates,
)


def test_firewall_requires_final_dated_source_and_environmental_coverage() -> None:
    dated = {
        "status": "pass_dated_source_m2m_scaleout_freeze",
        "candidate_image_pixels_opened": False,
        "files": {"dated_source_observation_manifest.csv": "abc"},
    }
    coverage = {
        "status": "pass_precolour_environmental_coverage",
        "coverage_gate_status": "pass_precolour_environmental_coverage",
        "source_stage": "final-dated-source",
        "final_dated_source_required": False,
        "scaleout_colour_opened": False,
        "image_acquisition_authorized": False,
    }
    validate_preimage_firewall_gates(
        dated,
        coverage,
        observation_manifest_name="dated_source_observation_manifest.csv",
        observation_manifest_sha256="abc",
    )
    coverage["source_stage"] = "live-feasibility"
    with pytest.raises(RuntimeError, match="environmental coverage"):
        validate_preimage_firewall_gates(
            dated,
            coverage,
            observation_manifest_name="dated_source_observation_manifest.csv",
            observation_manifest_sha256="abc",
        )


def test_acquisition_requires_the_exact_sealed_key_and_upstream_gate_hashes() -> None:
    firewall = {
        "status": "pass_scaleout_measurement_firewall",
        "candidate_image_pixels_opened": False,
        "coordinate_key_opened_by_measurement_worker": False,
        "sealed_keys": {"acquisition_coordinate_key.csv": "keyhash"},
        "dated_source_gate_sha256": "datedhash",
        "environmental_coverage_gate_sha256": "coveragehash",
    }
    validate_firewall_for_acquisition(
        firewall,
        key_name="acquisition_coordinate_key.csv",
        key_sha256="keyhash",
    )
    firewall["sealed_keys"]["acquisition_coordinate_key.csv"] = "other"
    with pytest.raises(RuntimeError, match="does not authorize"):
        validate_firewall_for_acquisition(
            firewall,
            key_name="acquisition_coordinate_key.csv",
            key_sha256="keyhash",
        )
