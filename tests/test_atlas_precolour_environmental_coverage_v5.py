from __future__ import annotations

import pytest

from scripts.data.run_jbi_atlas_precolour_environmental_coverage_v5 import (
    SOURCE_PASS,
    SOURCE_PROTOCOL,
    validate_source_v5_result,
)


def passing_source_v5():
    return {
        "protocol": SOURCE_PROTOCOL,
        "status": SOURCE_PASS,
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": 200,
        "selected_photo_assets": 60000,
        "frozen_observations": 60000,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "current_live_state_used_for_authorization": False,
        "repeat_35gb_stream_used_for_v5_authorization": False,
        "dated_snapshot_identity": {
            "identity_passed": True,
            "reused_existing_full_stream_proof": True,
        },
    }


def test_final_environmental_coverage_requires_source_v5() -> None:
    validate_source_v5_result(passing_source_v5())

    legacy = passing_source_v5()
    legacy["protocol"] = "jbi-atlas-dated-source-uuid-bucket-amendment-v4"
    legacy["status"] = "pass_dated_source_uuid_bucket_scaleout_freeze"
    with pytest.raises(ValueError, match="source role v5"):
        validate_source_v5_result(legacy)

    opened = passing_source_v5()
    opened["candidate_image_pixels_opened"] = True
    with pytest.raises(ValueError, match="source role v5"):
        validate_source_v5_result(opened)

    current_live = passing_source_v5()
    current_live["current_live_state_used_for_authorization"] = True
    with pytest.raises(ValueError, match="source role v5"):
        validate_source_v5_result(current_live)

    no_snapshot_proof = passing_source_v5()
    no_snapshot_proof["dated_snapshot_identity"]["identity_passed"] = False
    with pytest.raises(ValueError, match="source role v5"):
        validate_source_v5_result(no_snapshot_proof)
