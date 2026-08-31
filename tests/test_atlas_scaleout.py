from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import fcp_pipeline.atlas_scaleout as scaleout_module
from fcp_pipeline.atlas_scaleout import (
    AtlasScaleoutFreeze,
    exclude_reserved_scaleout_rows,
    freeze_scaleout_panels,
    live_api_scaleout_feasibility,
    qualify_scaleout_geometry,
)


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


def test_scaleout_geometry_precedes_cohort_draw() -> None:
    atlas = json.loads(
        Path("docs/supporting/jbi_image_first_atlas_contract_v1.json").read_text(
            encoding="utf-8"
        )
    )
    clustered = []
    for cluster in range(15):
        for within in range(20):
            clustered.append(
                {
                    "latitude": 20.0 + cluster * 1.2 + (within % 5) * 0.005,
                    "longitude": -120.0 + (within // 5) * 0.005,
                }
            )
    passing, audit = qualify_scaleout_geometry(
        [{"taxon_id": "1", "species": "Passing species", "genus": "Passing"}],
        {"1": clustered},
        atlas,
    )
    assert len(passing) == 1
    assert audit[0]["status"] == "geometry_eligible"
    assert {row["scale_km"] for row in audit[0]["scale_results"]} == {100, 250, 500}
    assert all(
        row["detectable_cells"] == len(row["detectable_cell_ids"])
        for row in audit[0]["scale_results"]
    )

    failing, failed_audit = qualify_scaleout_geometry(
        [{"taxon_id": "2", "species": "Failed species", "genus": "Failed"}],
        {"2": [{"latitude": 0.0, "longitude": 0.0} for _ in range(300)]},
        atlas,
    )
    assert failing == []
    assert failed_audit[0]["status"] == "primary_geometry_failed"

    observations = observation_rows()
    observations["0"][0]["colour_L"] = 50.0
    with pytest.raises(ValueError, match="outcome fields"):
        freeze_scaleout_panels(
            eligible_rows(), observations, contract(), source_role="unit-test metadata"
        )


def test_global_identity_reconciliation_is_deterministic() -> None:
    rows = [
        {"observation_id": "1", "photo_id": "10"},
        {"observation_id": "2", "photo_id": "20"},
        {"observation_id": "3", "photo_id": "10"},
        {"observation_id": "1", "photo_id": "40"},
    ]
    retained, removed = exclude_reserved_scaleout_rows(rows, {"1"}, {"10"})
    assert retained == [{"observation_id": "2", "photo_id": "20"}]
    assert removed == 3


def test_live_scaleout_geometry_receives_the_unmodified_atlas_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    atlas = json.loads(
        Path("docs/supporting/jbi_image_first_atlas_contract_v1.json").read_text(
            encoding="utf-8"
        )
    )
    expansion = contract()
    geometry_amendment = json.loads(
        Path(
            "docs/supporting/jbi_atlas_scaleout_geometry_admission_amendment_v1.json"
        ).read_text(encoding="utf-8")
    )
    global_id_amendment = json.loads(
        Path(
            "docs/supporting/jbi_atlas_scaleout_global_id_amendment_v1.json"
        ).read_text(encoding="utf-8")
    )

    class Adapter:
        def species_counts(self, _query: dict) -> list[dict]:
            return [
                {
                    "count": 1000 - index,
                    "taxon": {
                        "id": index + 1,
                        "name": f"Species {index + 1}",
                        "rank": "species",
                        "parent_id": index + 10_000,
                        "is_active": True,
                    },
                }
                for index in range(200)
            ]

        def observations(self, taxon_id: int, _query: dict) -> list[dict]:
            return [{"inat_taxon_id": taxon_id}]

    monkeypatch.setattr(
        scaleout_module,
        "_prepare_observation",
        lambda observation, taxon, _working: {
            "inat_taxon_id": str(taxon["id"]),
            "seed_taxon_id": int(observation["inat_taxon_id"]),
            "observation_id": f"seed-{taxon['id']}",
            "photo_id": f"seed-{taxon['id']}",
        },
    )

    def balanced(prepared: list[dict], _working: dict) -> list[dict]:
        taxon_id = int(prepared[0]["seed_taxon_id"])
        return [
            {
                "inat_taxon_id": str(taxon_id),
                "observation_id": str(taxon_id * 1000 + index),
                "photo_id": str(taxon_id * 1000 + index),
                "latitude": float(index),
                "longitude": float(index),
            }
            for index in range(300)
        ]

    monkeypatch.setattr(scaleout_module, "_balanced_selection", balanced)
    monkeypatch.setattr(
        scaleout_module,
        "_selection_qc",
        lambda selected, _candidate_count, _working: {
            "gate_pass": len(selected) == 300
        },
    )

    def qualify(candidates: list[dict], _rows: dict, received: dict, **_kwargs):
        assert received["admission"]["sample_size_tiers_descending"] == [500, 400, 300]
        candidate = dict(candidates[0])
        return [candidate], [
            {
                "taxon_id": candidate["taxon_id"],
                "status": "geometry_eligible",
                "scale_results": [],
            }
        ]

    monkeypatch.setattr(scaleout_module, "qualify_scaleout_geometry", qualify)
    monkeypatch.setattr(
        scaleout_module,
        "freeze_scaleout_panels",
        lambda eligible, _rows, _contract, **_kwargs: AtlasScaleoutFreeze(
            (), (), {"status": "pass_metadata_only_scaleout_freeze", "eligible_species": len(eligible)}
        ),
    )
    result = live_api_scaleout_feasibility(
        atlas,
        expansion,
        geometry_amendment,
        global_id_amendment,
        Adapter(),
        candidate_species_pool_size=200,
        maximum_candidates_per_species=1000,
    )
    assert result.audit["status"] == "pass_live_api_scaleout_feasibility"
    assert result.audit["geometry_eligible_species"] == 200
