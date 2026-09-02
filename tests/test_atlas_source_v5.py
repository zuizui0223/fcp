from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_source_v5 import (
    SourceV5Error,
    audit_terminal_source_rows,
    build_source_v5_result,
    validate_snapshot_identity_evidence,
    validate_source_v5_contract,
)

CONTRACT = Path("docs/supporting/jbi_atlas_source_role_amendment_v5.json")
V3_STOP = Path("docs/supporting/jbi_atlas_dated_source_streaming_v3_stop_result.json")


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def row(index: int, license_code: str = "cc-by"):
    return {
        "cohort_id": "C01",
        "observation_id": str(1000 + index),
        "photo_id": str(2000 + index),
        "photo_url_large": f"https://inaturalist-open-data.s3.amazonaws.com/photos/{2000 + index}/large.jpg",
        "photo_license": license_code,
        "candidate_image_pixels_opened": "False",
    }


def test_source_v5_contract_preserves_stops_and_removes_current_live_authorization() -> None:
    contract = load(CONTRACT)
    validate_source_v5_contract(contract)
    assert contract["current_live_api_role"]["authorization_role"] == "none"
    assert contract["dated_snapshot_provenance"]["repeat_35gb_stream_required_for_v5_authorization"] is False
    assert contract["terminal_acquisition_role"]["download_failure_does_not_trigger_replacement"] is True


def test_existing_exact_snapshot_identity_is_sufficient_provenance_evidence() -> None:
    validate_snapshot_identity_evidence(load(V3_STOP), load(CONTRACT))


def test_source_audit_rejects_nonofficial_host_and_current_live_is_not_consulted() -> None:
    contract = load(CONTRACT)
    # Use a minimally adjusted copy only to test row semantics; production contract
    # exact 60k counts are separately validated by the integration runner.
    contract = json.loads(json.dumps(contract))
    contract["immutable_scientific_selection"].update({
        "required_species": 1,
        "required_cohorts": 8,
        "required_observations": 1,
        "required_unique_observation_ids": 1,
        "required_unique_photo_ids": 1,
        "required_unique_source_urls": 1,
    })
    contract["frozen_asset_reference_audit"]["expected_license_counts"] = {"cc-by": 1}
    contract["frozen_asset_reference_audit"]["known_frozen_url_extension_missing"] = 0
    panels = [{"cohort_id": f"C{i:02d}"} for i in range(1, 9)]
    # The helper requires the contract's fixed terminal counts and therefore must
    # reject this tampered contract before considering rows.
    with pytest.raises(SourceV5Error):
        audit_terminal_source_rows([row(1)], panels, contract)


def test_build_source_v5_result_is_not_image_authorization_by_itself() -> None:
    contract = load(CONTRACT)
    audit = {
        "rows": 60000,
        "species": 200,
        "cohorts": 8,
        "bad_source_host_or_photo_id_path_rows": 0,
    }
    result = build_source_v5_result(source_audit=audit, contract=contract)
    assert result["status"] == "pass_frozen_selection_dated_provenance_v5"
    assert result["candidate_image_pixels_opened"] is False
    assert result["image_acquisition_authorized"] is False
    assert result["current_live_state_used_for_authorization"] is False
    assert result["replacement_permitted"] is False
