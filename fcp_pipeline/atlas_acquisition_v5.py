"""Authorization guard for terminal atlas v5 blinded image acquisition."""

from __future__ import annotations

from typing import Any, Mapping

from .atlas_measurement_v5 import validate_measurement_execution_contract


FIREWALL_STATUS = "pass_scaleout_measurement_firewall_v5"


def validate_v5_acquisition_firewall(
    firewall: Mapping[str, Any],
    *,
    key_name: str,
    key_sha256: str,
    contract: Mapping[str, Any],
    inference_v5: Mapping[str, Any],
) -> None:
    """Require the exact v5 firewall and sealed acquisition key before pixels open."""

    validate_measurement_execution_contract(contract, inference_v5)
    if (
        firewall.get("status") != FIREWALL_STATUS
        or firewall.get("protocol") != contract.get("protocol")
        or firewall.get("inference_version") != inference_v5.get("version")
        or firewall.get("candidate_image_pixels_opened") is not False
        or firewall.get("terminal_scaleout_colour_measured") is not False
        or firewall.get("coordinate_key_opened_by_measurement_worker") is not False
        or firewall.get("superseded_v3_ordered_inference_used") is not False
        or firewall.get("sealed_keys", {}).get(key_name) != key_sha256
    ):
        raise RuntimeError("v5 measurement firewall does not authorize this acquisition key")

    gate_hashes = firewall.get("preimage_gate_sha256", {})
    required = {
        "source_v5_manifest",
        "source_v5_result",
        "environmental_coverage",
        "roi_locked_result",
        "shared_transition_qualification",
    }
    if set(gate_hashes) != required or any(not str(gate_hashes[key]) for key in required):
        raise RuntimeError("v5 acquisition firewall lacks exact preimage gate hashes")

    if firewall.get("frozen_measurements") != 60000:
        raise RuntimeError("v5 acquisition firewall denominator changed")


__all__ = ["FIREWALL_STATUS", "validate_v5_acquisition_firewall"]
