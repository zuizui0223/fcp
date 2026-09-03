import pandas as pd

from fcp_pipeline.random_photo_pool import (
    equal_area_cell_bounds,
    freeze_random_photo_candidate_pool,
    inat_query_for_cell,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


class FakeClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def observations(self, params):
        self.calls.append(dict(params))
        return self.responses.pop(0)


def observation(obs_id, photo_id, species, lat, lon, *, license_code="cc-by"):
    return {
        "id": obs_id,
        "quality_grade": "research",
        "positional_accuracy": 100,
        "obscured": False,
        "geoprivacy": None,
        "observed_on": "2026-06-01",
        "geojson": {"coordinates": [lon, lat]},
        "taxon": {"id": obs_id + 1000, "rank": "species", "name": species},
        "user": {"id": obs_id + 2000, "login": f"observer{obs_id}"},
        "photos": [
            {
                "id": photo_id,
                "url": f"https://example.org/photos/{photo_id}/medium.jpg",
                "license_code": license_code,
                "attribution": "CC photo",
            }
        ],
    }


def test_equal_area_bounds_cover_full_world_without_latitude_overlap_gap():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    first = equal_area_cell_bounds(grid, 0)
    last = equal_area_cell_bounds(grid, grid.n_cells - 1)
    assert first["swlat"] == -90.0
    assert first["swlng"] == -180.0
    assert last["nelat"] == 90.0
    assert last["nelng"] == 180.0
    assert equal_area_cell_bounds(grid, 0)["nelat"] == equal_area_cell_bounds(grid, 4)["swlat"]


def test_cell_query_is_species_unfixed_flowering_photo_query():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    query = inat_query_for_cell(grid, 3)
    assert query["taxon_id"] == 47125
    assert query["term_id"] == 12
    assert query["term_value_id"] == 13
    assert query["rank"] == "species"
    assert query["order_by"] == "random"
    assert query["per_page"] == 200
    assert "species_id" not in query


def test_freeze_queries_each_cell_once_and_never_uses_colour():
    grid = EqualAreaGrid(n_lon=2, n_sinlat=2)
    cell_id, latitudes, longitudes = equal_area_cell_centers(grid)
    responses = []
    for cid, lat, lon in zip(cell_id, latitudes, longitudes, strict=True):
        responses.append(
            {
                "total_results": 1,
                "results": [
                    observation(
                        100 + int(cid),
                        500 + int(cid),
                        f"Genus species{int(cid)}",
                        float(lat),
                        float(lon),
                    )
                ],
            }
        )
    client = FakeClient(responses)
    frozen = freeze_random_photo_candidate_pool(
        client=client,
        grid=grid,
        per_cell_cap=20,
    )
    assert len(client.calls) == grid.n_cells
    assert len(frozen.observations) == grid.n_cells
    assert frozen.observations["cell_id"].nunique() == grid.n_cells
    assert frozen.manifest["counts"]["species"] == grid.n_cells
    assert frozen.manifest["outcome_firewall"]["candidate_image_pixels_opened"] is False
    assert frozen.manifest["outcome_firewall"]["morph_used_for_selection"] is False
    for call in client.calls:
        assert "morph" not in call
        assert "colour" not in call
        assert "color" not in call


def test_freeze_rejects_disallowed_license_and_wrong_rank_without_replacement():
    grid = EqualAreaGrid(n_lon=2, n_sinlat=2)
    cell_id, latitudes, longitudes = equal_area_cell_centers(grid)
    responses = []
    for cid, lat, lon in zip(cell_id, latitudes, longitudes, strict=True):
        valid = observation(
            100 + int(cid), 500 + int(cid), f"Genus species{int(cid)}", float(lat), float(lon)
        )
        if int(cid) == 0:
            valid["photos"][0]["license_code"] = None
        if int(cid) == 1:
            valid["taxon"]["rank"] = "genus"
        responses.append({"total_results": 1, "results": [valid]})
    client = FakeClient(responses)
    frozen = freeze_random_photo_candidate_pool(client=client, grid=grid, per_cell_cap=10)
    assert len(frozen.observations) == 2
    audit = frozen.cell_audit.set_index("cell_id")
    assert int(audit.loc[0, "rejected"]) == 1
    assert int(audit.loc[1, "rejected"]) == 1
    assert int(audit.loc[0, "accepted"]) == 0
    assert int(audit.loc[1, "accepted"]) == 0


def test_manifest_hash_is_order_stable_for_frozen_rows():
    grid = EqualAreaGrid(n_lon=2, n_sinlat=2)
    _, latitudes, longitudes = equal_area_cell_centers(grid)
    def make_client():
        return FakeClient([
            {
                "total_results": 1,
                "results": [observation(100 + i, 500 + i, f"Genus species{i}", float(latitudes[i]), float(longitudes[i]))],
            }
            for i in range(grid.n_cells)
        ])
    first = freeze_random_photo_candidate_pool(client=make_client(), grid=grid, per_cell_cap=10)
    second = freeze_random_photo_candidate_pool(client=make_client(), grid=grid, per_cell_cap=10)
    assert first.manifest["candidate_table_sha256"] == second.manifest["candidate_table_sha256"]
    pd.testing.assert_frame_equal(first.observations, second.observations)
