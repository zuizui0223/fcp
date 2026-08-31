from __future__ import annotations

from scripts.data.export_inaturalist_locked_analysis_table import (
    ADMITTED,
    EXPECTED_ROWS_PER_SPECIES,
    export_rows,
    observer_group,
)


def row(species: str, index: int, status: str = ADMITTED) -> dict[str, str]:
    return {
        "canonical_name": species,
        "encounter_blind_id": f"{species}-{index}",
        "encounter_status": status,
        "latitude": "35.0",
        "longitude": "135.0",
        "observer_id": f"person-{index % 4}",
    }


def test_observer_groups_are_stable_scoped_and_non_identifying() -> None:
    first = observer_group("Species one", "person-1", 2)
    assert first == observer_group("Species one", "person-1", 2)
    assert first != observer_group("Species two", "person-1", 2)
    assert "person-1" not in first
    assert first.startswith("observer_0002_")


def test_export_preserves_admitted_grouping_and_closes_non_admitted_rows() -> None:
    species = ("Species one", "Species two", "Species three")
    rows = [
        row(name, index, "not_evaluable" if index == 0 else ADMITTED)
        for name in species
        for index in range(EXPECTED_ROWS_PER_SPECIES)
    ]
    exported = export_rows(rows)
    assert len(exported) == 3 * EXPECTED_ROWS_PER_SPECIES
    rejected = [item for item in exported if item["encounter_status"] != ADMITTED]
    assert all(not item[field] for item in rejected for field in ("latitude", "longitude", "observer_id"))
    admitted = [item for item in exported if item["encounter_status"] == ADMITTED]
    assert all(item["observer_id"].startswith("observer_") for item in admitted)
    assert len({item["observer_id"] for item in admitted if item["canonical_name"] == species[0]}) == 4
