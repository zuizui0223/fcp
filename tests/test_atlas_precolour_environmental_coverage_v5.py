from __future__ import annotations

import pytest

from scripts.data.run_jbi_atlas_precolour_environmental_coverage_v5 import (
    DATED_PASS,
    DATED_PROTOCOL,
    validate_v4_dated_source,
)


def passing_v4():
    return {
        "protocol": DATED_PROTOCOL,
        "status": DATED_PASS,
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": 200,
        "selected_photo_assets": 60000,
        "frozen_observations": 60000,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
    }


def test_final_environmental_coverage_requires_uuid_source_v4() -> None:
    validate_v4_dated_source(passing_v4())

    legacy = passing_v4()
    legacy["protocol"] = "jbi-atlas-dated-source-streaming-execution-amendment-v3"
    legacy["status"] = "pass_dated_source_m2m_scaleout_freeze"
    with pytest.raises(ValueError, match="UUID dated-source v4"):
        validate_v4_dated_source(legacy)

    opened = passing_v4()
    opened["candidate_image_pixels_opened"] = True
    with pytest.raises(ValueError, match="UUID dated-source v4"):
        validate_v4_dated_source(opened)

    replaced = passing_v4()
    replaced["replacement_permitted"] = True
    with pytest.raises(ValueError, match="UUID dated-source v4"):
        validate_v4_dated_source(replaced)
