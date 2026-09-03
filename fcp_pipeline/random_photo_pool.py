"""Metadata-only freeze for the prospective random photo-first atlas.

The candidate pool is sampled by equal-area geographic cell, never by a fixed
species list and never using image pixels or colour.  The remote query uses one
randomly ordered iNaturalist page per cell; the exact returned observation/photo
IDs are then frozen as the experiment's candidate manifest.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import time
from typing import Any, Mapping, Protocol, Sequence
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from .shared_transition_surface import EqualAreaGrid, equal_area_cell_ids


INAT_API = "https://api.inaturalist.org/v1/observations"
DEFAULT_USER_AGENT = "fcp-random-photo-first-atlas/1.0 (github.com/zuizui0223/fcp)"


@dataclass(frozen=True)
class PoolFreeze:
    observations: pd.DataFrame
    cell_audit: pd.DataFrame
    manifest: dict[str, Any]


class ObservationClient(Protocol):
    def observations(self, params: Mapping[str, object]) -> Mapping[str, Any]:
        """Return one iNaturalist-compatible observation response."""


class InaturalistObservationClient:
    """Small unauthenticated iNaturalist v1 client with API-friendly pacing."""

    def __init__(
        self,
        *,
        api_url: str = INAT_API,
        user_agent: str = DEFAULT_USER_AGENT,
        request_interval_seconds: float = 1.05,
        timeout_seconds: float = 45.0,
        max_retries: int = 4,
    ) -> None:
        self.api_url = str(api_url)
        self.user_agent = str(user_agent)
        self.request_interval_seconds = float(request_interval_seconds)
        self.timeout_seconds = float(timeout_seconds)
        self.max_retries = int(max_retries)
        self._last_request_started = 0.0

    def observations(self, params: Mapping[str, object]) -> Mapping[str, Any]:
        query = urlencode({key: value for key, value in params.items() if value is not None})
        url = f"{self.api_url}?{query}"
        for attempt in range(self.max_retries + 1):
            elapsed = time.monotonic() - self._last_request_started
            if elapsed < self.request_interval_seconds:
                time.sleep(self.request_interval_seconds - elapsed)
            self._last_request_started = time.monotonic()
            request = Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": self.user_agent,
                },
            )
            try:
                with urlopen(request, timeout=self.timeout_seconds) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                if not isinstance(payload, Mapping):
                    raise RuntimeError("iNaturalist returned a non-object JSON payload")
                return payload
            except Exception:
                if attempt >= self.max_retries:
                    raise
                time.sleep(min(30.0, 2.0 ** attempt))
        raise RuntimeError("unreachable retry state")


def equal_area_cell_bounds(grid: EqualAreaGrid, cell_id: int) -> dict[str, float]:
    cell_id = int(cell_id)
    if cell_id < 0 or cell_id >= grid.n_cells:
        raise ValueError("cell_id is outside the grid")
    row = cell_id // grid.n_lon
    col = cell_id % grid.n_lon
    swlng = -180.0 + col * (360.0 / grid.n_lon)
    nelng = -180.0 + (col + 1) * (360.0 / grid.n_lon)
    sin_south = -1.0 + row * (2.0 / grid.n_sinlat)
    sin_north = -1.0 + (row + 1) * (2.0 / grid.n_sinlat)
    swlat = math.degrees(math.asin(float(np.clip(sin_south, -1.0, 1.0))))
    nelat = math.degrees(math.asin(float(np.clip(sin_north, -1.0, 1.0))))
    return {
        "swlat": float(swlat),
        "swlng": float(swlng),
        "nelat": float(nelat),
        "nelng": float(nelng),
    }


def inat_query_for_cell(
    grid: EqualAreaGrid,
    cell_id: int,
    *,
    per_page: int = 200,
    taxon_id: int = 47125,
    flowering_term_id: int = 12,
    flowering_term_value_id: int = 13,
    maximum_positional_accuracy_m: int = 5000,
) -> dict[str, object]:
    per_page = int(per_page)
    if per_page < 1 or per_page > 200:
        raise ValueError("per_page must lie in 1..200")
    params: dict[str, object] = {
        "taxon_id": int(taxon_id),
        "quality_grade": "research",
        "photos": "true",
        "geo": "true",
        "rank": "species",
        "term_id": int(flowering_term_id),
        "term_value_id": int(flowering_term_value_id),
        "acc_below": int(maximum_positional_accuracy_m),
        "order_by": "random",
        "per_page": per_page,
        "page": 1,
    }
    params.update(equal_area_cell_bounds(grid, cell_id))
    return params


def _coordinates(observation: Mapping[str, Any]) -> tuple[float | None, float | None]:
    geojson = observation.get("geojson") or {}
    coordinates = geojson.get("coordinates") or []
    if isinstance(coordinates, Sequence) and not isinstance(coordinates, (str, bytes)) and len(coordinates) >= 2:
        try:
            return float(coordinates[1]), float(coordinates[0])
        except (TypeError, ValueError):
            return None, None
    return None, None


def _eligible_photo(
    observation: Mapping[str, Any],
    *,
    allowed_photo_licenses: frozenset[str],
) -> Mapping[str, Any] | None:
    photos = observation.get("photos") or []
    eligible: list[Mapping[str, Any]] = []
    for photo in photos:
        if not isinstance(photo, Mapping):
            continue
        license_code = str(photo.get("license_code") or "").casefold()
        if license_code not in allowed_photo_licenses:
            continue
        if photo.get("id") in (None, "") or not photo.get("url"):
            continue
        eligible.append(photo)
    if not eligible:
        return None
    return min(eligible, key=lambda photo: int(photo["id"]))


def _large_photo_url(url: str) -> str:
    for size in ("square", "small", "medium"):
        token = f"/{size}."
        if token in url:
            return url.replace(token, "/large.")
    return url


def parse_candidate_observation(
    observation: Mapping[str, Any],
    *,
    expected_cell_id: int,
    grid: EqualAreaGrid,
    maximum_positional_accuracy_m: float,
    allowed_photo_licenses: frozenset[str],
) -> dict[str, object] | None:
    try:
        observation_id = int(observation["id"])
        positional_accuracy = float(observation["positional_accuracy"])
    except (KeyError, TypeError, ValueError):
        return None
    if not math.isfinite(positional_accuracy) or not (0.0 <= positional_accuracy <= maximum_positional_accuracy_m):
        return None
    if bool(observation.get("obscured")):
        return None
    if observation.get("geoprivacy") not in (None, "", "open"):
        return None
    if str(observation.get("quality_grade") or "") != "research":
        return None

    taxon = observation.get("taxon") or {}
    if not isinstance(taxon, Mapping) or str(taxon.get("rank") or "") != "species":
        return None
    species = str(taxon.get("name") or "").strip()
    if not species or " " not in species:
        return None

    latitude, longitude = _coordinates(observation)
    if latitude is None or longitude is None:
        return None
    if not (-90.0 <= latitude <= 90.0 and -180.0 <= longitude <= 180.0):
        return None
    actual_cell_id = int(
        equal_area_cell_ids(
            np.asarray([latitude], dtype=float),
            np.asarray([longitude], dtype=float),
            grid,
        )[0]
    )
    # Make cell membership half-open and deterministic even when the API bbox
    # includes a boundary point in more than one request.
    if actual_cell_id != int(expected_cell_id):
        return None

    photo = _eligible_photo(
        observation,
        allowed_photo_licenses=allowed_photo_licenses,
    )
    if photo is None:
        return None
    user = observation.get("user") or {}
    observer_id = user.get("id") if isinstance(user, Mapping) else None
    if observer_id in (None, ""):
        return None

    return {
        "cell_id": actual_cell_id,
        "observation_id": observation_id,
        "photo_id": int(photo["id"]),
        "photo_url_large": _large_photo_url(str(photo["url"])),
        "photo_license": str(photo.get("license_code") or ""),
        "attribution": str(photo.get("attribution") or ""),
        "species": species,
        "inat_taxon_id": int(taxon.get("id") or 0),
        "latitude": float(latitude),
        "longitude": float(longitude),
        "positional_accuracy_m": float(positional_accuracy),
        "observed_on": str(observation.get("observed_on") or ""),
        "observer_id": str(observer_id),
        "observer": str(user.get("login") or "") if isinstance(user, Mapping) else "",
    }


def _sha256_rows(rows: pd.DataFrame) -> str:
    canonical = rows.sort_values(["cell_id", "observation_id", "photo_id"]).to_csv(index=False, lineterminator="\n")
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def freeze_random_photo_candidate_pool(
    *,
    client: ObservationClient,
    grid: EqualAreaGrid,
    per_cell_cap: int = 200,
    taxon_id: int = 47125,
    flowering_term_id: int = 12,
    flowering_term_value_id: int = 13,
    maximum_positional_accuracy_m: int = 5000,
    allowed_photo_licenses: Sequence[str] = (
        "cc0",
        "cc-by",
        "cc-by-sa",
        "cc-by-nc",
        "cc-by-nc-sa",
    ),
) -> PoolFreeze:
    """Query each equal-area cell once and freeze the exact random metadata pool."""

    per_cell_cap = int(per_cell_cap)
    if per_cell_cap < 1 or per_cell_cap > 200:
        raise ValueError("per_cell_cap must lie in 1..200")
    allowed = frozenset(str(value).casefold() for value in allowed_photo_licenses)
    if not allowed:
        raise ValueError("allowed_photo_licenses cannot be empty")

    rows: list[dict[str, object]] = []
    audit: list[dict[str, object]] = []
    seen_observation_ids: set[int] = set()
    seen_photo_ids: set[int] = set()

    for cell_id in range(grid.n_cells):
        params = inat_query_for_cell(
            grid,
            cell_id,
            per_page=per_cell_cap,
            taxon_id=taxon_id,
            flowering_term_id=flowering_term_id,
            flowering_term_value_id=flowering_term_value_id,
            maximum_positional_accuracy_m=maximum_positional_accuracy_m,
        )
        payload = client.observations(params)
        raw_results = payload.get("results") or []
        if not isinstance(raw_results, Sequence):
            raise RuntimeError(f"cell {cell_id} response has non-list results")
        accepted = 0
        duplicate = 0
        rejected = 0
        for observation in raw_results:
            if not isinstance(observation, Mapping):
                rejected += 1
                continue
            parsed = parse_candidate_observation(
                observation,
                expected_cell_id=cell_id,
                grid=grid,
                maximum_positional_accuracy_m=float(maximum_positional_accuracy_m),
                allowed_photo_licenses=allowed,
            )
            if parsed is None:
                rejected += 1
                continue
            observation_id = int(parsed["observation_id"])
            photo_id = int(parsed["photo_id"])
            if observation_id in seen_observation_ids or photo_id in seen_photo_ids:
                duplicate += 1
                continue
            seen_observation_ids.add(observation_id)
            seen_photo_ids.add(photo_id)
            rows.append(parsed)
            accepted += 1
        audit.append(
            {
                "cell_id": cell_id,
                "api_total_results": int(payload.get("total_results") or len(raw_results)),
                "api_returned": int(len(raw_results)),
                "accepted": accepted,
                "rejected": rejected,
                "duplicate": duplicate,
                **equal_area_cell_bounds(grid, cell_id),
            }
        )

    observations = pd.DataFrame(rows)
    if len(observations) == 0:
        raise RuntimeError("candidate pool contains zero eligible observations")
    observations = observations.sort_values(
        ["cell_id", "observation_id", "photo_id"]
    ).reset_index(drop=True)
    cell_audit = pd.DataFrame(audit).sort_values("cell_id").reset_index(drop=True)

    species_counts = observations["species"].value_counts()
    cell_counts = observations["cell_id"].value_counts()
    manifest = {
        "protocol": "random-photo-first-candidate-pool-v1",
        "status": "metadata_pool_frozen_before_candidate_image_pixels",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "grid": {
            "n_lon": grid.n_lon,
            "n_sinlat": grid.n_sinlat,
            "n_cells": grid.n_cells,
            "cell_area_km2": grid.cell_area_km2,
        },
        "query": {
            "taxon_id": int(taxon_id),
            "flowering_term_id": int(flowering_term_id),
            "flowering_term_value_id": int(flowering_term_value_id),
            "quality_grade": "research",
            "rank": "species",
            "photos": True,
            "geo": True,
            "order_by": "random",
            "per_cell_cap": per_cell_cap,
            "maximum_positional_accuracy_m": int(maximum_positional_accuracy_m),
            "allowed_photo_licenses": sorted(allowed),
            "api_requests": grid.n_cells,
            "pages_per_cell": 1,
        },
        "outcome_firewall": {
            "candidate_image_pixels_opened": False,
            "morph_used_for_selection": False,
            "continuous_colour_used_for_selection": False,
            "species_list_fixed_or_targeted": False,
            "legacy_pr21_terminal_records_used": False,
        },
        "counts": {
            "observations": int(len(observations)),
            "unique_observation_ids": int(observations["observation_id"].nunique()),
            "unique_photo_ids": int(observations["photo_id"].nunique()),
            "species": int(observations["species"].nunique()),
            "occupied_cells": int(observations["cell_id"].nunique()),
            "maximum_photos_per_species": int(species_counts.max()),
            "maximum_photos_per_cell": int(cell_counts.max()),
        },
        "candidate_table_sha256": _sha256_rows(observations),
    }
    return PoolFreeze(
        observations=observations,
        cell_audit=cell_audit,
        manifest=manifest,
    )
