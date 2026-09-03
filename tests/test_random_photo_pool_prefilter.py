from fcp_pipeline.random_photo_pool import inat_query_for_cell
from fcp_pipeline.shared_transition_surface import EqualAreaGrid


def test_one_shot_query_prefilters_public_coordinates_and_allowed_photo_licenses():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    query = inat_query_for_cell(grid, 3)
    assert query["obscuration"] == "none"
    assert set(str(query["photo_license"]).split(",")) == {
        "cc0",
        "cc-by",
        "cc-by-sa",
        "cc-by-nc",
        "cc-by-nc-sa",
    }
    assert query["order_by"] == "random"
    assert query["per_page"] == 200


def test_custom_license_allowlist_is_normalized_before_random_page():
    grid = EqualAreaGrid(n_lon=4, n_sinlat=2)
    query = inat_query_for_cell(
        grid,
        0,
        allowed_photo_licenses=("CC-BY-SA", "cc0", "cc-by-sa"),
    )
    assert query["photo_license"] == "cc-by-sa,cc0"
