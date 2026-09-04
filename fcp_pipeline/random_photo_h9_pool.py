"""Fresh metadata acquisition for prospective H9 individual-distance replication."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .random_photo_pool import (
    DEFAULT_ALLOWED_PHOTO_LICENSES,
    InaturalistObservationClient,
)


@dataclass(frozen=True)
class H9Freeze:
    observations: pd.DataFrame
    species_audit: pd.DataFrame
    manifest: dict[str, Any]


def h9_query_for_species(
    taxon_id: int,
    *,
    per_page: int = 200,
    maximum_positional_accuracy_m: int = 5000,
    flowering_term_id: int = 12,
    flowering_term_value_id: int = 13,
    allowed_photo_licenses: Sequence[str] = DEFAULT_ALLOWED_PHOTO_LICENSES,
) -> dict[str, object]:
    per_page = int(per_page)
    if per_page < 1 or per_page > 200:
        raise ValueError("per_page must lie in 1..200")
    allowed = sorted({str(x).casefold() for x in allowed_photo_licenses})
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
        "order_by": "random",
        "per_page": per_page,
        "page": 1,
    }


def _coords(observation: Mapping[str, Any]) -> tuple[float | None, float | None]:
    geojson = observation.get("geojson") or {}
    coords = geojson.get("coordinates") or []
    if isinstance(coords, Sequence) and not isinstance(coords, (str, bytes)) and len(coords) >= 2:
        try:
            return float(coords[1]), float(coords[0])
        except (TypeError, ValueError):
            pass
    return None, None


def _photo(observation: Mapping[str, Any], allowed: frozenset[str]) -> Mapping[str, Any] | None:
    eligible = []
    for p in observation.get("photos") or []:
        if not isinstance(p, Mapping):
            continue
        if str(p.get("license_code") or "").casefold() not in allowed:
            continue
        if p.get("id") in (None, "") or not p.get("url"):
            continue
        eligible.append(p)
    if not eligible:
        return None
    return min(eligible, key=lambda x: int(x["id"]))


def _large_url(url: str) -> str:
    for size in ("square", "small", "medium"):
        token = f"/{size}."
        if token in url:
            return url.replace(token, "/large.")
    return url


def parse_h9_observation(
    observation: Mapping[str, Any],
    *,
    expected_taxon_id: int,
    maximum_positional_accuracy_m: float,
    allowed_photo_licenses: frozenset[str],
) -> dict[str, object] | None:
    try:
        obs_id = int(observation["id"])
        accuracy = float(observation["positional_accuracy"])
    except (KeyError, TypeError, ValueError):
        return None
    if not math.isfinite(accuracy) or not (0.0 <= accuracy <= maximum_positional_accuracy_m):
        return None
    if bool(observation.get("obscured")):
        return None
    if observation.get("geoprivacy") not in (None, "", "open"):
        return None
    if str(observation.get("quality_grade") or "") != "research":
        return None
    taxon = observation.get("taxon") or {}
    if not isinstance(taxon, Mapping):
        return None
    if str(taxon.get("rank") or "") != "species" or int(taxon.get("id") or -1) != int(expected_taxon_id):
        return None
    species = str(taxon.get("name") or "").strip()
    if not species:
        return None
    lat, lon = _coords(observation)
    if lat is None or lon is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
        return None
    p = _photo(observation, allowed_photo_licenses)
    if p is None:
        return None
    user = observation.get("user") or {}
    if not isinstance(user, Mapping) or user.get("id") in (None, ""):
        return None
    return {
        "species": species,
        "inat_taxon_id": int(expected_taxon_id),
        "observation_id": obs_id,
        "photo_id": int(p["id"]),
        "photo_url_large": _large_url(str(p["url"])),
        "photo_license": str(p.get("license_code") or ""),
        "attribution": str(p.get("attribution") or ""),
        "latitude": float(lat),
        "longitude": float(lon),
        "positional_accuracy_m": float(accuracy),
        "observed_on": str(observation.get("observed_on") or ""),
        "observer_id": str(user["id"]),
        "observer": str(user.get("login") or ""),
    }


def _stable_hash(row: Mapping[str, object]) -> str:
    return hashlib.sha256(f"{int(row['observation_id'])}:{int(row['photo_id'])}".encode()).hexdigest()


def observer_cap(rows: pd.DataFrame, cap: int) -> pd.DataFrame:
    if cap < 1:
        raise ValueError("observer cap must be positive")
    if rows.empty:
        return rows.copy()
    x = rows.copy()
    x["row_hash"] = [
        _stable_hash({"observation_id": o, "photo_id": p})
        for o, p in zip(x["observation_id"], x["photo_id"])
    ]
    x = x.sort_values(["observer_id", "row_hash"], kind="mergesort")
    x = x.groupby("observer_id", sort=False, observed=True).head(int(cap)).copy()
    return x.sort_values("row_hash", kind="mergesort").reset_index(drop=True)


def _distance_matrix(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    phi = np.deg2rad(np.asarray(lat, dtype=float))
    lam = np.deg2rad(np.asarray(lon, dtype=float))
    c = np.cos(phi)
    xyz = np.column_stack([c * np.cos(lam), c * np.sin(lam), np.sin(phi)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * 6371.0088


def geographic_maximin(rows: pd.DataFrame, n: int) -> pd.DataFrame:
    n = int(n)
    if n < 1:
        raise ValueError("n must be positive")
    if len(rows) <= n:
        return rows.copy().reset_index(drop=True)
    x = rows.copy().reset_index(drop=True)
    if "row_hash" not in x.columns:
        x["row_hash"] = [
            _stable_hash({"observation_id": o, "photo_id": p})
            for o, p in zip(x["observation_id"], x["photo_id"])
        ]
    d = _distance_matrix(x["latitude"].to_numpy(), x["longitude"].to_numpy())
    mean_d = d.mean(axis=1)
    best = np.flatnonzero(np.isclose(mean_d, mean_d.max(), rtol=0, atol=1e-12))
    first = min(best.tolist(), key=lambda i: str(x.loc[i, "row_hash"]))
    selected = [first]
    remaining = set(range(len(x))) - {first}
    while len(selected) < n:
        indices = sorted(remaining)
        min_d = d[np.ix_(indices, selected)].min(axis=1)
        max_val = float(min_d.max())
        candidates = [indices[j] for j, v in enumerate(min_d) if math.isclose(float(v), max_val, rel_tol=0, abs_tol=1e-12)]
        chosen = min(candidates, key=lambda i: str(x.loc[i, "row_hash"]))
        selected.append(chosen)
        remaining.remove(chosen)
    out = x.loc[selected].copy().reset_index(drop=True)
    out["h9_selection_order"] = np.arange(1, len(out) + 1, dtype=int)
    return out


def freeze_h9_metadata(
    *,
    client: object,
    species_frame: pd.DataFrame,
    exclusion_observation_ids: set[int],
    exclusion_photo_ids: set[int],
    per_page: int = 200,
    observer_cap_n: int = 2,
    fixed_raw_photos: int = 60,
    maximum_positional_accuracy_m: int = 5000,
    allowed_photo_licenses: Sequence[str] = DEFAULT_ALLOWED_PHOTO_LICENSES,
) -> H9Freeze:
    allowed = frozenset(str(x).casefold() for x in allowed_photo_licenses)
    retained_frames = []
    audits = []
    seen_obs: set[int] = set()
    seen_photo: set[int] = set()
    request_errors = 0
    for row in species_frame.sort_values("inat_taxon_id", kind="mergesort").itertuples(index=False):
        taxon_id = int(row.inat_taxon_id)
        species = str(row.species)
        params = h9_query_for_species(
            taxon_id,
            per_page=per_page,
            maximum_positional_accuracy_m=maximum_positional_accuracy_m,
            allowed_photo_licenses=allowed,
        )
        try:
            payload = client.observations(params)
            results = payload.get("results") or []
            error = ""
        except Exception as exc:
            results = []
            error = f"{type(exc).__name__}:{str(exc)[:200]}"
            request_errors += 1
        parsed_rows = []
        rejected = wrong_taxon = prior = duplicate = 0
        for obs in results:
            if not isinstance(obs, Mapping):
                rejected += 1
                continue
            taxon = obs.get("taxon") or {}
            if isinstance(taxon, Mapping) and int(taxon.get("id") or -1) != taxon_id:
                wrong_taxon += 1
                continue
            parsed = parse_h9_observation(
                obs,
                expected_taxon_id=taxon_id,
                maximum_positional_accuracy_m=float(maximum_positional_accuracy_m),
                allowed_photo_licenses=allowed,
            )
            if parsed is None:
                rejected += 1
                continue
            oid, pid = int(parsed["observation_id"]), int(parsed["photo_id"])
            if oid in exclusion_observation_ids or pid in exclusion_photo_ids:
                prior += 1
                continue
            if oid in seen_obs or pid in seen_photo:
                duplicate += 1
                continue
            parsed_rows.append(parsed)
        eligible = pd.DataFrame(parsed_rows)
        capped = observer_cap(eligible, observer_cap_n) if len(eligible) else eligible
        selected = geographic_maximin(capped, fixed_raw_photos) if len(capped) else capped
        if len(selected):
            selected["query_species"] = species
            retained_frames.append(selected)
            seen_obs.update(selected["observation_id"].astype(int).tolist())
            seen_photo.update(selected["photo_id"].astype(int).tolist())
        audits.append({
            "species": species,
            "inat_taxon_id": taxon_id,
            "raw_results": int(len(results)),
            "locally_eligible_before_exclusion": int(len(parsed_rows) + prior),
            "local_rejected": int(rejected),
            "wrong_taxon": int(wrong_taxon),
            "prior_excluded": int(prior),
            "fresh_duplicate": int(duplicate),
            "after_observer_cap": int(len(capped)),
            "retained": int(len(selected)),
            "full_fixed_n": bool(len(selected) == int(fixed_raw_photos)),
            "request_error": error,
        })
    observations = pd.concat(retained_frames, ignore_index=True) if retained_frames else pd.DataFrame()
    audit = pd.DataFrame(audits)
    full_species = int(audit["full_fixed_n"].sum()) if len(audit) else 0
    manifest = {
        "query_attempts": int(len(species_frame)),
        "request_errors": int(request_errors),
        "retained_fresh_photos": int(len(observations)),
        "selected_species": int(len(species_frame)),
        "full_fixed_n_species": full_species,
    }
    return H9Freeze(observations=observations, species_audit=audit, manifest=manifest)


__all__ = [
    "H9Freeze",
    "freeze_h9_metadata",
    "geographic_maximin",
    "h9_query_for_species",
    "observer_cap",
    "parse_h9_observation",
]
