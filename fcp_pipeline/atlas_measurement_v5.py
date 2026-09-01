"""Pre-image authorization and location-blind measurement rules for atlas v5.

This module deliberately reuses only low-level blinding/sharding helpers from the
legacy v3 measurement implementation.  It never validates or imports the
superseded v3 ordered inference tree; terminal biological decisions are governed
exclusively by ``atlas_inference_cascade`` and the v5 contract.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from .atlas_inference_cascade import BRANCHES, validate_contract as validate_v5_inference
from .atlas_measurement import build_measurement_firewall, measurement_shard


PROTOCOL = "jbi-atlas-measurement-execution-v5"


class MeasurementV5ContractError(ValueError):
    pass


def validate_measurement_execution_contract(
    contract: Mapping[str, Any], inference_v5: Mapping[str, Any]
) -> None:
    """Fail closed if the v5 pre-image execution contract drifts."""

    validate_v5_inference(inference_v5)
    if contract.get("protocol") != PROTOCOL:
        raise MeasurementV5ContractError("unexpected measurement execution protocol")
    if contract.get("status") != (
        "prospectively_frozen_before_any_terminal_scaleout_candidate_pixel"
    ):
        raise MeasurementV5ContractError("measurement execution was not frozen pre-pixel")
    firewall = contract.get("outcome_firewall", {})
    if any(value is not False for value in firewall.values()):
        raise MeasurementV5ContractError("v5 measurement contract contains opened outcomes")

    v5_firewall = contract.get("v5_inference_firewall", {})
    if (
        v5_firewall.get("only_inference_contract")
        != "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"
        or v5_firewall.get(
            "superseded_v3_ordered_inference_must_not_authorize_or_classify_terminal_results"
        )
        is not True
        or tuple(v5_firewall.get("confirmatory_sequence", ())) != BRANCHES
        or v5_firewall.get("not_evaluable_never_advances_confirmatory") is not True
        or v5_firewall.get("same_frozen_colour_field_must_be_reused_downstream")
        is not True
    ):
        raise MeasurementV5ContractError("v5 inference firewall changed")

    measurement = contract.get("location_blind_measurement", {})
    expected_measurement = {
        "blinding_salt": "fcp-atlas-v3-coordinate-firewall",
        "salt_inherited_unchanged_from_v3": True,
        "minimum_admitted_fraction_per_species": 0.7,
        "minimum_admitted_observations_per_species": 210,
        "minimum_background_control_fraction_among_admitted": 0.7,
        "minimum_evaluable_species_per_cohort": 20,
        "minimum_evaluable_species_total": 160,
        "all_eight_cohorts_must_finish": True,
    }
    for key, value in expected_measurement.items():
        if measurement.get(key) != value:
            raise MeasurementV5ContractError(f"v5 measurement rule changed: {key}")

    execution = contract.get("technical_execution", {})
    expected_execution = {
        "acquisition_shard_count": 16,
        "measurement_shard_count": 16,
        "maximum_concurrent_workers": 8,
        "acquisition_retries_per_image": 4,
        "acquisition_timeout_seconds": 60.0,
        "acquisition_pause_seconds_after_each_terminal_image": 0.05,
        "measurement_torch_threads_per_worker": 2,
        "all_shards_required": True,
        "early_stopping": False,
    }
    for key, value in expected_execution.items():
        if execution.get(key) != value:
            raise MeasurementV5ContractError(f"v5 technical execution changed: {key}")

    gates = contract.get("pixel_opening_gates", {})
    if gates.get("all_required") is not True:
        raise MeasurementV5ContractError("all v5 pixel-opening gates must be required")
    if gates.get("dated_source", {}).get("required_status") != (
        "pass_dated_source_m2m_scaleout_freeze"
    ):
        raise MeasurementV5ContractError("dated-source pass label changed")
    if gates.get("environmental_coverage", {}).get("required_status") != (
        "pass_precolour_environmental_coverage"
    ):
        raise MeasurementV5ContractError("environmental coverage pass label changed")
    if gates.get("roi_v4_locked", {}).get("required_status") != (
        "pass_roi_v4_locked_test"
    ):
        raise MeasurementV5ContractError("ROI v4 pass label changed")
    if gates.get("shared_transition_method", {}).get("required_status") != "pass":
        raise MeasurementV5ContractError("shared-transition method gate changed")


def validate_preimage_gates(
    *,
    dated_reconciliation: Mapping[str, Any],
    dated_manifest: Mapping[str, Any],
    observation_manifest_name: str,
    observation_manifest_sha256: str,
    environmental_coverage: Mapping[str, Any],
    roi_locked_result: Mapping[str, Any],
    shared_qualification_result: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> None:
    """Require every independent pre-image gate before building the sealed key."""

    gates = contract["pixel_opening_gates"]
    dated_gate = gates["dated_source"]
    if (
        dated_reconciliation.get("status") != dated_gate["required_status"]
        or dated_reconciliation.get("candidate_image_pixels_opened") is not False
        or dated_reconciliation.get("selected_species") != dated_gate["required_species"]
        or dated_reconciliation.get("selected_photo_assets")
        != dated_gate["required_photo_assets"]
        or dated_reconciliation.get("frozen_observations")
        != dated_gate["required_frozen_observations"]
        or dated_reconciliation.get("replacement_permitted") is not False
        or dated_reconciliation.get("image_acquisition_authorized") is not False
    ):
        raise RuntimeError("dated-source v5 gate did not pass unchanged")
    if (
        dated_manifest.get("status") != dated_gate["required_status"]
        or dated_manifest.get("candidate_image_pixels_opened") is not False
        or dated_manifest.get("files", {}).get(observation_manifest_name)
        != observation_manifest_sha256
    ):
        raise RuntimeError("dated-source manifest does not authorize exact frozen rows")

    environment_gate = gates["environmental_coverage"]
    if (
        environmental_coverage.get("status") != environment_gate["required_status"]
        or environmental_coverage.get("coverage_gate_status")
        != environment_gate["required_status"]
        or environmental_coverage.get("source_stage")
        != environment_gate["required_source_stage"]
        or environmental_coverage.get("final_dated_source_required") is not False
        or environmental_coverage.get("scaleout_colour_opened") is not False
        or environmental_coverage.get("image_acquisition_authorized") is not False
    ):
        raise RuntimeError("final pre-colour environmental coverage did not pass")

    roi_gate = gates["roi_v4_locked"]
    if (
        roi_locked_result.get("status") != roi_gate["required_status"]
        or roi_locked_result.get("scaleout_candidate_pixels_opened") is not False
        or roi_locked_result.get("scaleout_candidate_pixels_permitted") is not True
    ):
        raise RuntimeError("ROI v4 locked gate did not authorize scaleout")

    shared_gate = gates["shared_transition_method"]
    shared_scope = shared_qualification_result.get("scope", {})
    if (
        shared_qualification_result.get("status") != shared_gate["required_status"]
        or shared_scope.get("method_gate_only") is not True
        or shared_scope.get("biological_support_claimed") is not False
        or shared_scope.get("candidate_image_pixels_opened") is not False
    ):
        raise RuntimeError("shared-transition preimage method qualification did not pass")


def build_v5_measurement_firewall(
    rows: Sequence[Mapping[str, Any]], contract: Mapping[str, Any]
) -> dict[str, list[dict[str, Any]]]:
    """Use the unchanged pre-pixel blinding transformation under a v5 contract."""

    salt = str(contract["location_blind_measurement"]["blinding_salt"])
    return build_measurement_firewall(rows, salt=salt)


__all__ = [
    "build_v5_measurement_firewall",
    "measurement_shard",
    "validate_measurement_execution_contract",
    "validate_preimage_gates",
]
