from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.atlas_scaleout import freeze_scaleout_panels


CONTRACT = Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json")


def contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def eligible_rows() -> list[dict[str, str]]:
    return [
        {
            "taxon_id": str(index),
            "species": f"Species {index}",
            "genus": f"Genus{index}",
        }
        for index in range(205)
    ]


def observation_rows() -> dict[str, list[dict[str, object]]]:
    return {
        str(taxon): [
            {
                "species": f"Species {taxon}",
                "inat_taxon_id": taxon,
                "observation_id": taxon * 1000 + index,
                "photo_id": str(taxon * 1000 + index),
                "latitude": float(index % 90),
                "longitude": float((index * 7) % 180),
            }
            for index in range(300)
        ]
        for taxon in range(205)
    }


def test_scaleout_freezes_all_panels_without_opening_pixels() -> None:
    frozen = freeze_scaleout_panels(
        eligible_rows(),
        observation_rows(),
        contract(),
        source_role="unit-test metadata",
    )
    assert len(frozen.panels) == 200
    assert len(frozen.observations) == 60_000
    assert {row["cohort_id"] for row in frozen.panels} == {
        f"C{index:02d}" for index in range(1, 9)
    }
    assert not any(row["candidate_image_pixels_opened"] for row in frozen.observations)
    assert frozen.audit["status"] == "pass_metadata_only_scaleout_freeze"


def test_scaleout_rejects_short_species_and_outcome_leakage() -> None:
    observations = observation_rows()
    observations["0"] = observations["0"][:-1]
    with pytest.raises(ValueError, match="requires exactly 300"):
        freeze_scaleout_panels(
            eligible_rows(), observations, contract(), source_role="unit-test metadata"
        )

    observations = observation_rows()
    observations["0"][0]["colour_L"] = 50.0
    with pytest.raises(ValueError, match="outcome fields"):
        freeze_scaleout_panels(
            eligible_rows(), observations, contract(), source_role="unit-test metadata"
        )
