from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from fcp_pipeline.image_first_atlas import (
    freeze_atlas_geometry,
    freeze_atlas_metadata,
    species_free_display_rows,
    validate_atlas_contract,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = json.loads(
    (ROOT / "docs" / "supporting" / "jbi_image_first_atlas_contract_v1.json").read_text(
        encoding="utf-8"
    )
)


class FakeMetadataAdapter:
    def __init__(self, species: int = 55):
        self.species = species

    def species_counts(self, query):
        assert query["term_id"] == 12
        assert query["term_value_id"] == 13
        assert query["acc_below"] == 5000
        return [
            {
                "count": 100_000 - index,
                "taxon": {
                    "id": 10_000 + index,
                    "name": f"Atlas species {index:02d}",
                    "rank": "species",
                    "parent_id": 20_000 + index,
                    "is_active": True,
                },
            }
            for index in range(self.species)
        ]

    def observations(self, taxon_id, query):
        assert query["limit"] == 1000
        index = taxon_id - 10_000
        rows = []
        for observation_index in range(600):
            cell = observation_index % 120
            latitude = -45.0 + (cell % 30) * 0.5 + (observation_index // 120) * 0.01
            longitude = -150.0 + (cell // 30) * 30.0 + (cell % 5) * 0.5
            month = observation_index % 12 + 1
            observation_id = taxon_id * 10_000 + observation_index
            rows.append(
                {
                    "id": observation_id,
                    "positional_accuracy": 1000,
                    "obscured": False,
                    "geoprivacy": None,
                    "geojson": {"coordinates": [longitude, latitude]},
                    "observed_on": f"2020-{month:02d}-{observation_index % 27 + 1:02d}",
                    "photos": [
                        {
                            "id": observation_id * 10,
                            "license_code": "cc-by",
                            "url": f"https://example.invalid/{observation_id}/medium.jpg",
                            "attribution": "synthetic test",
                        }
                    ],
                    "user": {
                        "id": f"observer-{index}-{observation_index % 100}",
                        "login": f"observer-{observation_index % 100}",
                    },
                }
            )
        return rows


def test_contract_freezes_image_first_mainline_before_pixels():
    validate_atlas_contract(CONTRACT)
    assert CONTRACT["scientific_roles"]["active_mainline"] == "image-first global flower-colour atlas"
    assert CONTRACT["outcome_firewall"]["candidate_image_pixels_opened"] is False
    assert CONTRACT["outcome_firewall"]["literature_classification_used_for_admission"] is False
    assert CONTRACT["parent_freeze"]["stage_a_primary_p_lower"] == pytest.approx(0.0113)
    assert CONTRACT["parent_freeze"]["stage_b_primary_p_upper"] == pytest.approx(0.0906)


def test_readme_promotes_atlas_and_demotes_literature_classification_to_support():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert text.startswith("# FCP — image-first global flower-colour atlas")
    assert "species-conditioned transition boundary" in text
    assert "species-free" in text
    assert "34-species classifications are no longer the atlas cohort selector" in text


def test_metadata_freeze_admits_50_species_at_300_to_500_without_colour():
    freeze = freeze_atlas_metadata(CONTRACT, FakeMetadataAdapter())
    assert freeze.audit["status"] == "pass_50_species_metadata_only"
    assert len(freeze.cohort) == 50
    assert len(freeze.observations) == 25_000
    assert {row["selected_photographs"] for row in freeze.cohort} == {500}
    assert all(float(row["positional_accuracy_m"]) <= 5000 for row in freeze.observations)
    assert all(row["literature_classification_used_for_admission"] is False for row in freeze.cohort)
    assert freeze.audit["continuous_colour_used"] is False


def test_metadata_shortfall_stops_without_relaxing_rules():
    freeze = freeze_atlas_metadata(CONTRACT, FakeMetadataAdapter(species=40))
    assert freeze.audit["status"] == "not_evaluable_insufficient_species"
    assert len(freeze.cohort) == 40
    assert freeze.audit["next_gate"].startswith("STOP")


def test_colour_dependent_admission_contract_fails_closed():
    bad = copy.deepcopy(CONTRACT)
    bad["outcome_firewall"]["continuous_colour_used"] = True
    with pytest.raises(ValueError, match="outcome firewall"):
        validate_atlas_contract(bad)


class ColourAccessTrap(dict):
    def __getitem__(self, key):
        if key in {"colour", "values", "roi"}:
            raise AssertionError("geometry selection accessed an image-derived field")
        return super().__getitem__(key)


def test_geometry_selects_finest_passing_scale_without_colour_access():
    contract = copy.deepcopy(CONTRACT)
    criteria = contract["geometry_only_scale_selection"]["passing_criteria"]
    criteria.update(
        {
            "minimum_retained_edges_per_species": 20,
            "minimum_detectable_cells_per_species": 1,
            "minimum_evaluable_species": 50,
            "minimum_cells_A_ge_3": 1,
            "minimum_species_with_shared_opportunity": 50,
        }
    )
    rows = []
    for species_index in range(50):
        for row_index in range(60):
            rows.append(
                ColourAccessTrap(
                    {
                        "species": f"Species {species_index:02d}",
                        "latitude": 35.0 + (row_index % 10) * 0.01,
                        "longitude": 135.0 + (row_index // 10) * 0.01,
                        "values": object(),
                    }
                )
            )
    result = freeze_atlas_geometry(rows, contract)
    assert result["status"] == "geometry_scale_frozen"
    assert result["selected_primary_scale_km"] == 100
    assert result["selection_used_continuous_colour"] is False


def test_species_free_display_strips_taxonomy_but_keeps_roi_and_colour():
    displayed = species_free_display_rows(
        [
            {
                "display_id": "atlas-1",
                "photo_id": "123",
                "species": "Hidden species",
                "inat_taxon_id": 42,
                "latitude": 35.0,
                "longitude": 135.0,
                "roi_thumbnail_path": "roi/123.webp",
                "colour_L": 50.0,
                "colour_a": 20.0,
                "colour_b": -10.0,
                "colour_hex": "#aa3377",
                "photo_bar_order": 1,
            }
        ]
    )
    assert displayed[0]["roi_thumbnail_path"] == "roi/123.webp"
    assert "species" not in displayed[0]
    assert "inat_taxon_id" not in displayed[0]
