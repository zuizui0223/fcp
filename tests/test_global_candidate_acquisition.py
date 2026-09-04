from __future__ import annotations

from fcp_pipeline.global_candidate_acquisition import (
    deterministic_candidate_pages,
    stable_candidate_query,
)


LICENSES = ("cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa")


def test_candidate_pages_are_deterministic_unique_and_bounded():
    a = deterministic_candidate_pages(12345, 12500, seed=20260917)
    b = deterministic_candidate_pages(12345, 12500, seed=20260917)
    assert a == b
    assert len(a) == 3
    assert len(set(a)) == 3
    assert all(1 <= page <= 50 for page in a)


def test_candidate_pages_use_all_available_when_fewer_than_three():
    pages = deterministic_candidate_pages(12345, 350, per_page=200, seed=20260917)
    assert set(pages) == {1, 2}
    assert len(pages) == 2
    assert deterministic_candidate_pages(12345, 0) == ()


def test_candidate_page_choice_varies_with_taxon_but_not_outcome_data():
    left = deterministic_candidate_pages(1001, 20000, seed=20260917)
    right = deterministic_candidate_pages(1002, 20000, seed=20260917)
    assert left != right
    assert len(left) == len(right) == 3


def test_stable_candidate_query_never_uses_random_ordering():
    query = stable_candidate_query(
        12345,
        page=7,
        allowed_photo_licenses=LICENSES,
    )
    assert query["order_by"] == "id"
    assert query["order"] == "desc"
    assert query["page"] == 7
    assert query["per_page"] == 200
    assert query["taxon_id"] == 12345
    assert "random" not in {str(value).casefold() for value in query.values()}
