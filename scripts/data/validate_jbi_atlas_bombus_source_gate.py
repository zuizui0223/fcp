#!/usr/bin/env python3
"""Validate the fail-closed pre-colour Bombus source decision."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


RESULT = Path("docs/supporting/jbi_atlas_bombus_source_gate_result_v1.json")
INFERENCE = Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json")
EXPECTED_INFERENCE_SHA256 = "598b34bd0c996f2744e7dd9444d10673bdbe39092ba7a3625ce1f9fc003c604d"


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    if canonical_sha256(INFERENCE) != EXPECTED_INFERENCE_SHA256:
        raise RuntimeError("parent inference contract changed")
    if (
        result.get("protocol") != "jbi-atlas-bombus-source-gate-v1"
        or result.get("status")
        != "pollinator_biogeographic_concordance_not_evaluable_precolour_source_gate"
        or result.get("gbif_taxon", {}).get("usage_key") != 1340278
        or int(result.get("public_occurrence_count_at_check", 0)) <= int(
            result.get("public_search_api_hard_record_ceiling", 0)
        )
        or result.get("complete_download_requires_registered_user_authentication") is not True
        or result.get("download_requested") is not False
        or result.get("download_doi_available") is not False
        or result.get("truncated_search_subset_used") is not False
        or result.get("candidate_image_pixels_opened") is not False
        or result.get("scaleout_colour_opened") is not False
        or result.get("environment_colour_join_performed") is not False
    ):
        raise RuntimeError("Bombus source gate changed or admitted a truncated source")
    print(
        json.dumps(
            {
                "status": result["status"],
                "gbif_usage_key": result["gbif_taxon"]["usage_key"],
                "occurrence_count_at_check": result["public_occurrence_count_at_check"],
                "truncated_search_subset_used": False,
                "candidate_image_pixels_opened": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
