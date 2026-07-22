#!/usr/bin/env python3
"""Build reproducible geography covariates from GBIF occurrence search results.

This script intentionally avoids inferring island status or biological populations from
country names or occurrence clusters. It derives quantities directly supported by
coordinate-bearing GBIF occurrences: occurrence count, robust geographic ranges,
nearest-neighbour spacing, spatial extent, occupied grid cells, and distance-threshold
connected components. Results are checkpointed after every species so interrupted CI
runs retain completed work.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

API = "https://api.gbif.org/v1/occurrence/search"
EARTH_RADIUS_KM = 6371.0088
FIELDS = [
    "canonical_name",
    "family",
    "gbif_occurrences",
    "gbif_coordinate_records_sampled",
    "absolute_latitude",
    "median_latitude",
    "median_longitude",
    "latitudinal_range",
    "longitudinal_range",
    "latitude_q05",
    "latitude_q95",
    "longitude_q05",
    "longitude_q95",
    "occupied_grid_cells_1deg",
    "median_nearest_neighbor_km",
    "spatial_extent_q95_km",
    "components_50km",
    "largest_component_fraction_50km",
    "components_100km",
    "largest_component_fraction_100km",
    "gbif_query_status",
]


def read_rows(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def request_json(url: str, timeout: int, retries: int, backoff: float) -> dict:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "fcp-gbif-geography/1.1 (research workflow)",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            last = exc
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            if retry_after:
                try:
                    wait = max(float(retry_after), backoff * (2**attempt))
                except ValueError:
                    wait = backoff * (2**attempt)
            else:
                wait = backoff * (2**attempt)
        except Exception as exc:  # network failures are recorded per species
            last = exc
            wait = backoff * (2**attempt)
        if attempt + 1 < retries:
            time.sleep(wait)
    raise RuntimeError(str(last))


def pairwise_haversine_km(latitudes: np.ndarray, longitudes: np.ndarray) -> np.ndarray:
    """Return a dense pairwise great-circle distance matrix in kilometres."""
    lat = np.radians(latitudes)
    lon = np.radians(longitudes)
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def connected_component_summary(distances: np.ndarray, threshold_km: float) -> tuple[int, float]:
    """Summarise occurrence connectivity without labelling components as populations."""
    n = distances.shape[0]
    if n == 0:
        return 0, float("nan")
    parent = list(range(n))
    size = [1] * n

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri == rj:
            return
        if size[ri] < size[rj]:
            ri, rj = rj, ri
        parent[rj] = ri
        size[ri] += size[rj]

    rows, cols = np.where(np.triu(distances <= threshold_km, k=1))
    for i, j in zip(rows.tolist(), cols.tolist()):
        union(i, j)
    component_sizes: dict[int, int] = {}
    for i in range(n):
        root = find(i)
        component_sizes[root] = component_sizes.get(root, 0) + 1
    return len(component_sizes), max(component_sizes.values()) / n


def summarize(name: str, family: str, payload: dict) -> dict[str, object]:
    total = int(payload.get("count") or 0)
    coordinates: list[tuple[float, float]] = []
    for item in payload.get("results") or []:
        if not isinstance(item, dict):
            continue
        try:
            lat = float(item.get("decimalLatitude"))
            lon = float(item.get("decimalLongitude"))
        except (TypeError, ValueError):
            continue
        if math.isfinite(lat) and math.isfinite(lon) and -90 <= lat <= 90 and -180 <= lon <= 180:
            coordinates.append((lat, lon))

    if coordinates:
        arr = np.asarray(coordinates, dtype=float)
        latitudes = arr[:, 0]
        longitudes = arr[:, 1]
        lat_q05, median_lat, lat_q95 = np.quantile(latitudes, [0.05, 0.5, 0.95])
        lon_q05, median_lon, lon_q95 = np.quantile(longitudes, [0.05, 0.5, 0.95])
        absolute_latitude = abs(float(median_lat))
        latitudinal_range = float(lat_q95 - lat_q05)
        longitudinal_range = float(lon_q95 - lon_q05)
        occupied_grid_cells = len({(math.floor(lat), math.floor(lon)) for lat, lon in coordinates})
        distances = pairwise_haversine_km(latitudes, longitudes)
        if len(coordinates) >= 2:
            nearest = distances.copy()
            np.fill_diagonal(nearest, np.inf)
            median_nearest_neighbor = float(np.median(np.min(nearest, axis=1)))
        else:
            median_nearest_neighbor = 0.0
        centre_distances = pairwise_haversine_km(
            np.asarray([median_lat, *latitudes]),
            np.asarray([median_lon, *longitudes]),
        )[0, 1:]
        spatial_extent_q95 = float(np.quantile(centre_distances, 0.95))
        components_50, largest_50 = connected_component_summary(distances, 50.0)
        components_100, largest_100 = connected_component_summary(distances, 100.0)
        status = "complete"
    else:
        lat_q05 = median_lat = lat_q95 = absolute_latitude = latitudinal_range = ""
        lon_q05 = median_lon = lon_q95 = longitudinal_range = ""
        occupied_grid_cells = median_nearest_neighbor = spatial_extent_q95 = ""
        components_50 = largest_50 = components_100 = largest_100 = ""
        status = "no_coordinate_records"

    return {
        "canonical_name": name,
        "family": family,
        "gbif_occurrences": total,
        "gbif_coordinate_records_sampled": len(coordinates),
        "absolute_latitude": absolute_latitude,
        "median_latitude": median_lat,
        "median_longitude": median_lon,
        "latitudinal_range": latitudinal_range,
        "longitudinal_range": longitudinal_range,
        "latitude_q05": lat_q05,
        "latitude_q95": lat_q95,
        "longitude_q05": lon_q05,
        "longitude_q95": lon_q95,
        "occupied_grid_cells_1deg": occupied_grid_cells,
        "median_nearest_neighbor_km": median_nearest_neighbor,
        "spatial_extent_q95_km": spatial_extent_q95,
        "components_50km": components_50,
        "largest_component_fraction_50km": largest_50,
        "components_100km": components_100,
        "largest_component_fraction_100km": largest_100,
        "gbif_query_status": status,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", default="data/global_flower_colour_species_ranked_resolved.csv")
    parser.add_argument("--out", default="data/species_gbif_geography_covariates.csv")
    parser.add_argument("--qc-out", default="data/species_gbif_geography_covariates_qc.json")
    parser.add_argument("--limit", type=int, default=300)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--backoff", type=float, default=2.0)
    parser.add_argument("--request-delay", type=float, default=0.25)
    parser.add_argument("--max-species", type=int, default=0)
    args = parser.parse_args()

    species = read_rows(args.species)
    if args.max_species > 0:
        species = species[: args.max_species]
    out_path = Path(args.out)
    completed: dict[str, dict[str, object]] = {}
    if out_path.exists():
        completed = {row["canonical_name"]: row for row in read_rows(str(out_path))}

    errors: list[dict[str, str]] = []
    ordered: list[dict[str, object]] = []
    for index, row in enumerate(species, start=1):
        name = row.get("canonical_name", "").strip()
        family = row.get("family", "").strip()
        if not name:
            continue
        cached = completed.get(name)
        cache_is_current = cached and all(field in cached for field in FIELDS)
        if cache_is_current and cached.get("gbif_query_status") in {"complete", "no_coordinate_records"}:
            result = cached
        else:
            params = {
                "scientificName": name,
                "hasCoordinate": "true",
                "occurrenceStatus": "present",
                "limit": min(max(args.limit, 1), 300),
            }
            url = API + "?" + urllib.parse.urlencode(params)
            try:
                payload = request_json(url, args.timeout, args.retries, args.backoff)
                result = summarize(name, family, payload)
            except RuntimeError as exc:
                result = {field: "" for field in FIELDS}
                result.update(
                    {
                        "canonical_name": name,
                        "family": family,
                        "gbif_coordinate_records_sampled": 0,
                        "gbif_query_status": "request_failed",
                    }
                )
                errors.append({"canonical_name": name, "error": str(exc)[:500]})
            if index < len(species):
                time.sleep(max(0.0, args.request_delay))
        completed[name] = result
        ordered = [completed[r["canonical_name"]] for r in species[:index] if r.get("canonical_name") in completed]
        write_rows(out_path, ordered)
        qc = {
            "species_requested": len(species),
            "species_completed": index,
            "complete_coordinate_species": sum(r.get("gbif_query_status") == "complete" for r in ordered),
            "no_coordinate_species": sum(r.get("gbif_query_status") == "no_coordinate_records" for r in ordered),
            "failed_species": sum(r.get("gbif_query_status") == "request_failed" for r in ordered),
            "errors": errors,
            "record_limit_per_species": min(max(args.limit, 1), 300),
            "complete": index == len(species),
            "method": (
                "GBIF occurrence search; robust coordinate ranges, 1-degree occupied cells, median nearest-neighbour "
                "distance, 95th-percentile distance from the coordinate median, and 50/100-km threshold graph "
                "components from up to 300 coordinate records per species"
            ),
            "interpretation_guardrail": (
                "Distance-threshold components describe sampled occurrence connectivity and are not inferred biological populations."
            ),
        }
        Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(qc, ensure_ascii=False))


if __name__ == "__main__":
    main()
