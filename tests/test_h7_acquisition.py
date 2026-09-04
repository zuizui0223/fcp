from __future__ import annotations

import pandas as pd

from fcp_pipeline.h7_acquisition import _observer_capped_selection, freeze_h7_fresh_metadata
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


def _obs(obs_id: int, photo_id: int, taxon_id: int, species: str, lat: float, lon: float, observer: int):
    return {
        "id": obs_id,
        "quality_grade": "research",
        "positional_accuracy": 10.0,
        "obscured": False,
        "geoprivacy": None,
        "geojson": {"coordinates": [lon, lat]},
        "taxon": {"id": taxon_id, "rank": "species", "name": species},
        "photos": [
            {
                "id": photo_id,
                "url": f"https://static.inaturalist.org/photos/{photo_id}/square.jpg",
                "license_code": "cc-by",
                "attribution": "test",
            }
        ],
        "user": {"id": observer, "login": f"u{observer}"},
        "observed_on": "2026-01-01",
    }


class SequentialClient:
    def __init__(self, payloads):
        self.payloads = list(payloads)
        self.calls = []

    def observations(self, params):
        self.calls.append(dict(params))
        payload = self.payloads[len(self.calls) - 1]
        if isinstance(payload, Exception):
            raise payload
        return payload


def _targets():
    grid = EqualAreaGrid(4, 2)
    ids, lats, lons = equal_area_cell_centers(grid)
    centers = {int(i): (float(la), float(lo)) for i, la, lo in zip(ids, lats, lons)}
    rows = []
    for species, taxon, cells in [("Alpha beta", 101, [1, 2]), ("Gamma delta", 202, [5, 6])]:
        for order, cell in enumerate(cells, 1):
            lat, lon = centers[cell]
            rows.append(
                {
                    "species": species,
                    "inat_taxon_id": taxon,
                    "target_cell_order": order,
                    "cell_id": cell,
                    "cell_center_latitude": lat,
                    "cell_center_longitude": lon,
                }
            )
    return grid, pd.DataFrame(rows), centers


def test_observer_cap_selection_is_deterministic_and_capped():
    rows = [
        {"observation_id": 1, "photo_id": 11, "observer_id": "x"},
        {"observation_id": 2, "photo_id": 12, "observer_id": "x"},
        {"observation_id": 3, "photo_id": 13, "observer_id": "y"},
    ]
    a = _observer_capped_selection(rows, observer_cap=1, retained_cap=2, seed=7, taxon_id=10, cell_id=3)
    b = _observer_capped_selection(list(reversed(rows)), observer_cap=1, retained_cap=2, seed=7, taxon_id=10, cell_id=3)
    assert [x["observation_id"] for x in a] == [x["observation_id"] for x in b]
    assert len(a) == 2
    assert len({x["observer_id"] for x in a}) == 2


def test_fresh_freeze_excludes_prior_ids_and_never_retries_request_errors():
    grid, targets, centers = _targets()
    payloads = []
    next_id = 1000
    for species, taxon, cells in [("Alpha beta", 101, [1, 2]), ("Gamma delta", 202, [5, 6])]:
        for cell in cells:
            lat, lon = centers[cell]
            results = [
                _obs(next_id, next_id + 10000, taxon, species, lat, lon, observer=next_id),
                _obs(next_id + 1, next_id + 10001, taxon, species, lat, lon, observer=next_id + 1),
                _obs(next_id + 2, next_id + 10002, taxon, species, lat, lon, observer=next_id + 2),
            ]
            payloads.append({"results": results})
            next_id += 10
    # Make the final target a hard request error. It must be recorded, not retried.
    payloads[-1] = RuntimeError("network failure")
    client = SequentialClient(payloads)
    # Exclude one returned observation and a different returned photo from the first target.
    exclusions = pd.DataFrame(
        {
            "observation_id": [1000, 999999],
            "photo_id": [999998, 11001],
        }
    )

    result = freeze_h7_fresh_metadata(
        client=client,
        targets=targets,
        exclusions=exclusions,
        grid=grid,
        per_page=3,
        observer_cap=1,
        retained_cap=1,
        selection_seed=77,
        expected_species=2,
        expected_targets=4,
        required_species_for_gate=1,
        required_full_cells_per_species=2,
    )

    assert len(client.calls) == 4
    assert result.manifest["query_attempts"] == 4
    assert result.manifest["query_retries"] == 0
    assert result.manifest["request_error_targets"] == 1
    assert result.target_audit.iloc[-1]["retained"] == 0
    assert "network failure" in result.target_audit.iloc[-1]["request_error"]
    assert not result.observations["observation_id"].isin({1000, 999999}).any()
    assert not result.observations["photo_id"].isin({999998, 11001}).any()
    assert result.observations["observation_id"].is_unique
    assert result.observations["photo_id"].is_unique
    assert result.manifest["premeasurement_gate"]["pass"] is True


def test_wrong_taxon_is_rejected_locally():
    grid, targets, centers = _targets()
    one = targets.iloc[[0]].copy()
    lat, lon = centers[int(one.iloc[0]["cell_id"])]
    client = SequentialClient(
        [{"results": [_obs(1, 11, 999, "Wrong taxon", lat, lon, observer=1)]}]
    )
    result = freeze_h7_fresh_metadata(
        client=client,
        targets=one,
        exclusions=pd.DataFrame({"observation_id": [], "photo_id": []}),
        grid=grid,
        per_page=1,
        observer_cap=1,
        retained_cap=1,
        expected_species=1,
        expected_targets=1,
        required_species_for_gate=1,
        required_full_cells_per_species=1,
    )
    assert result.target_audit.iloc[0]["wrong_taxon"] == 1
    assert result.target_audit.iloc[0]["retained"] == 0
    assert result.manifest["premeasurement_gate"]["pass"] is False
