"""Outcome-blind candidate acquisition primitives for the global Monte Carlo atlas.

These helpers deliberately avoid iNaturalist random ordering.  A stable taxon query is
paired with a deterministic hash ordering of available pages, so the actual candidate
image pool is a prospectively reproducible draw rather than a cache-sensitive random
page.  No colour outcome enters page choice.
"""
from __future__ import annotations

import hashlib
import math
from typing import Sequence


DEFAULT_PAGE_SEED = 20260917


def deterministic_candidate_pages(
    taxon_id: int,
    total_results: int,
    *,
    per_page: int = 200,
    maximum_api_page: int = 50,
    pages_per_species: int = 3,
    seed: int = DEFAULT_PAGE_SEED,
) -> tuple[int, ...]:
    """Return a fixed hash-ranked set of candidate pages for one species.

    The set depends only on taxon identity, the metadata-only result count and frozen
    constants.  It never depends on colour, observer composition, geographic spread or
    whether an earlier page would already have supplied the requested target.
    """
    taxon_id = int(taxon_id)
    total_results = int(total_results)
    per_page = int(per_page)
    maximum_api_page = int(maximum_api_page)
    pages_per_species = int(pages_per_species)
    if taxon_id < 1:
        raise ValueError("taxon_id must be positive")
    if total_results < 0:
        raise ValueError("total_results cannot be negative")
    if per_page < 1 or per_page > 200:
        raise ValueError("per_page must lie in 1..200")
    if maximum_api_page < 1:
        raise ValueError("maximum_api_page must be positive")
    if pages_per_species < 1:
        raise ValueError("pages_per_species must be positive")
    if total_results == 0:
        return ()
    available = min(maximum_api_page, int(math.ceil(total_results / per_page)))
    pages = list(range(1, available + 1))
    pages.sort(
        key=lambda page: hashlib.sha256(
            f"{int(seed)}|{taxon_id}|{page}".encode("utf-8")
        ).hexdigest()
    )
    return tuple(pages[: min(pages_per_species, available)])


def stable_candidate_query(
    taxon_id: int,
    *,
    page: int,
    per_page: int = 200,
    maximum_positional_accuracy_m: int = 5000,
    flowering_term_id: int = 12,
    flowering_term_value_id: int = 13,
    allowed_photo_licenses: Sequence[str],
) -> dict[str, object]:
    """Build the cache-resistant metadata query used for candidate acquisition."""
    page = int(page)
    if page < 1:
        raise ValueError("page must be positive")
    per_page = int(per_page)
    if per_page < 1 or per_page > 200:
        raise ValueError("per_page must lie in 1..200")
    allowed = sorted({str(value).casefold() for value in allowed_photo_licenses})
    if not allowed:
        raise ValueError("allowed_photo_licenses cannot be empty")
    return {
        "taxon_id": int(taxon_id),
        "quality_grade": "research",
        "photos": "true",
        "geo": "true",
        "rank": "species",
        "term_id": int(flowering_term_id),
        "term_value_id": int(flowering_term_value_id),
        "acc_below": int(maximum_positional_accuracy_m),
        "obscuration": "none",
        "photo_license": ",".join(allowed),
        "order_by": "id",
        "order": "desc",
        "per_page": per_page,
        "page": page,
    }


__all__ = [
    "DEFAULT_PAGE_SEED",
    "deterministic_candidate_pages",
    "stable_candidate_query",
]
