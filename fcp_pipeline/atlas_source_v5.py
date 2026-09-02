"""Immutable selection + dated provenance gate for terminal atlas source v5."""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Sequence
from urllib.parse import urlparse

PROTOCOL = "jbi-atlas-source-role-amendment-v5"
PASS_LABEL = "pass_frozen_selection_dated_provenance_v5"


class SourceV5Error(ValueError):
    pass


def validate_source_v5_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != PROTOCOL:
        raise SourceV5Error("unexpected source v5 protocol")
    if contract.get("status") != "prospectively_frozen_after_v4_stop_and_before_any_terminal_scaleout_candidate_pixel":
        raise SourceV5Error("source v5 was not frozen before terminal pixels")
    trigger = contract.get("trigger", {})
    if (
        trigger.get("candidate_image_pixels_opened") is not False
        or trigger.get("terminal_scaleout_colour_measured") is not False
        or trigger.get("replacement_or_resampling_permitted") is not False
        or trigger.get("v3_stop", {}).get("git_blob_sha") != "0742186c425cb506d16b1cda4c274534cd8d5a8d"
        or trigger.get("v4_stop", {}).get("git_blob_sha") != "3f8feb8fa79728606d0a69b4cdbf8719f732eabd"
    ):
        raise SourceV5Error("source v5 trigger or preserved STOP identity changed")
    selection = contract.get("immutable_scientific_selection", {})
    if (
        selection.get("artifact_id") != 9769165047
        or selection.get("required_species") != 200
        or selection.get("required_cohorts") != 8
        or selection.get("required_observations") != 60000
        or selection.get("required_unique_observation_ids") != 60000
        or selection.get("required_unique_photo_ids") != 60000
        or selection.get("required_unique_source_urls") != 60000
        or selection.get("selection_used_candidate_pixels") is not False
        or selection.get("selection_used_continuous_colour") is not False
        or selection.get("selection_must_not_be_rerun_from_current_live_state") is not True
    ):
        raise SourceV5Error("immutable 200-species scientific selection changed")
    asset = contract.get("frozen_asset_reference_audit", {})
    if (
        asset.get("required_scheme") != "https"
        or asset.get("required_host") != "inaturalist-open-data.s3.amazonaws.com"
        or asset.get("required_path_prefix_rule") != "/photos/{photo_id}/large."
        or asset.get("known_frozen_url_extension_missing") != 14
        or asset.get("required_bad_host_or_photo_id_path_count") != 0
    ):
        raise SourceV5Error("frozen asset-reference audit changed")
    snapshot = contract.get("dated_snapshot_provenance", {})
    if (
        snapshot.get("snapshot_date") != "2026-08-27"
        or snapshot.get("content_length_bytes") != 35093052336
        or snapshot.get("sha256") != "c98202c07796b275fe41fc1518fc394ac09caf2dede370a4ee64ce6d68b0c50d"
        or snapshot.get("required_identity_passed") is not True
        or snapshot.get("required_bytes_read") != 35093052336
        or snapshot.get("required_computed_sha256") != "c98202c07796b275fe41fc1518fc394ac09caf2dede370a4ee64ce6d68b0c50d"
        or snapshot.get("repeat_35gb_stream_required_for_v5_authorization") is not False
    ):
        raise SourceV5Error("dated snapshot provenance changed")
    current = contract.get("current_live_api_role", {})
    if (
        current.get("authorization_role") != "none"
        or current.get("allowed_role") != "diagnostic only"
        or any(current.get(key) is not False for key in (
            "may_redefine_frozen_species",
            "may_redefine_frozen_photo",
            "may_redefine_frozen_license_record",
            "may_remove_a_frozen_row",
            "may_trigger_replacement_or_resampling",
        ))
    ):
        raise SourceV5Error("current-live API regained an authorization role")
    acquisition = contract.get("terminal_acquisition_role", {})
    if (
        acquisition.get("image_access_occurs_only_after_all_independent_preimage_gates_pass") is not True
        or acquisition.get("download_failure_is_terminal_measurement_failure") is not True
        or acquisition.get("download_failure_remains_in_60000_denominator") is not True
        or acquisition.get("download_failure_does_not_trigger_replacement") is not True
        or acquisition.get("acquisition_success_or_failure_cannot_change_species_or_cohort_membership") is not True
    ):
        raise SourceV5Error("terminal acquisition role changed")
    passed = contract.get("source_v5_pass_rule", {})
    if (
        passed.get("pass_label") != PASS_LABEL
        or passed.get("current_live_state_not_used_for_authorization") is not True
        or passed.get("candidate_image_pixels_opened") is not False
        or passed.get("replacement_or_resampling_permitted") is not False
        or passed.get("image_acquisition_authorized_by_this_gate_alone") is not False
    ):
        raise SourceV5Error("source v5 pass rule changed")


def audit_terminal_source_rows(
    rows: Sequence[Mapping[str, Any]],
    panels: Sequence[Mapping[str, Any]],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    validate_source_v5_contract(contract)
    selection = contract["immutable_scientific_selection"]
    if len(rows) != selection["required_observations"]:
        raise SourceV5Error("terminal source row denominator changed")
    if len(panels) != selection["required_species"]:
        raise SourceV5Error("terminal source species panel changed")
    cohorts = {str(row.get("cohort_id") or "") for row in panels}
    if cohorts != {f"C{index:02d}" for index in range(1, 9)}:
        raise SourceV5Error("terminal source cohort set changed")
    if any(str(row.get("candidate_image_pixels_opened", "")).casefold() != "false" for row in rows):
        raise SourceV5Error("terminal source rows report opened pixels")

    observation_ids = [str(row.get("observation_id") or "") for row in rows]
    photo_ids = [str(row.get("photo_id") or "") for row in rows]
    urls = [str(row.get("photo_url_large") or "") for row in rows]
    if (
        "" in observation_ids
        or "" in photo_ids
        or "" in urls
        or len(set(observation_ids)) != selection["required_unique_observation_ids"]
        or len(set(photo_ids)) != selection["required_unique_photo_ids"]
        or len(set(urls)) != selection["required_unique_source_urls"]
    ):
        raise SourceV5Error("terminal source identifier uniqueness changed")

    asset = contract["frozen_asset_reference_audit"]
    allowed = {str(value).casefold() for value in asset["allowed_frozen_photo_licenses"]}
    license_counts: Counter[str] = Counter()
    bad_source_rows: list[str] = []
    extension_missing = 0
    for row in rows:
        photo_id = str(row["photo_id"])
        parsed = urlparse(str(row["photo_url_large"]))
        if (
            parsed.scheme != asset["required_scheme"]
            or parsed.netloc != asset["required_host"]
            or not parsed.path.startswith(f"/photos/{photo_id}/large.")
        ):
            bad_source_rows.append(photo_id)
        basename = parsed.path.rsplit("/", 1)[-1]
        if basename == "large.":
            extension_missing += 1
        license_code = str(row.get("photo_license") or "").casefold()
        if license_code not in allowed:
            raise SourceV5Error(f"terminal frozen licence is outside the allowed set: {license_code!r}")
        license_counts[license_code] += 1
    if bad_source_rows:
        raise SourceV5Error(f"terminal source host/path identity changed for {len(bad_source_rows)} rows")
    if dict(sorted(license_counts.items())) != dict(sorted(asset["expected_license_counts"].items())):
        raise SourceV5Error("terminal frozen licence-count ledger changed")
    if extension_missing != asset["known_frozen_url_extension_missing"]:
        raise SourceV5Error("known frozen missing-extension count changed")

    return {
        "rows": len(rows),
        "species": len(panels),
        "cohorts": len(cohorts),
        "unique_observation_ids": len(set(observation_ids)),
        "unique_photo_ids": len(set(photo_ids)),
        "unique_source_urls": len(set(urls)),
        "official_bucket_host_path_rows": len(rows),
        "bad_source_host_or_photo_id_path_rows": 0,
        "frozen_license_counts": dict(sorted(license_counts.items())),
        "known_frozen_url_extension_missing": extension_missing,
    }


def validate_snapshot_identity_evidence(evidence: Mapping[str, Any], contract: Mapping[str, Any]) -> None:
    validate_source_v5_contract(contract)
    snapshot = contract["dated_snapshot_provenance"]
    identity = evidence.get("snapshot_identity", {})
    if (
        identity.get("identity_passed") is not True
        or identity.get("bytes_read") != snapshot["required_bytes_read"]
        or identity.get("expected_bytes") != snapshot["required_bytes_read"]
        or identity.get("computed_sha256") != snapshot["required_computed_sha256"]
        or identity.get("expected_sha256") != snapshot["required_computed_sha256"]
        or identity.get("archive_persisted") is not False
    ):
        raise SourceV5Error("exact dated snapshot identity evidence does not pass")


def build_source_v5_result(
    *,
    source_audit: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    validate_source_v5_contract(contract)
    required = contract["immutable_scientific_selection"]
    if (
        source_audit.get("rows") != required["required_observations"]
        or source_audit.get("species") != required["required_species"]
        or source_audit.get("cohorts") != required["required_cohorts"]
        or source_audit.get("bad_source_host_or_photo_id_path_rows") != 0
    ):
        raise SourceV5Error("source audit cannot satisfy the frozen pass rule")
    return {
        "protocol": PROTOCOL,
        "status": PASS_LABEL,
        "selected_species": 200,
        "selected_photo_assets": 60000,
        "frozen_observations": 60000,
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "current_live_state_used_for_authorization": False,
        "repeat_35gb_stream_used_for_v5_authorization": False,
        "source_audit": dict(source_audit),
        "claim_ceiling": contract["claim_ceiling"],
    }


__all__ = [
    "PASS_LABEL",
    "PROTOCOL",
    "SourceV5Error",
    "audit_terminal_source_rows",
    "build_source_v5_result",
    "validate_snapshot_identity_evidence",
    "validate_source_v5_contract",
]
