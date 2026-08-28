#!/usr/bin/env python3
"""Compute robust geographic extent covariates from quality-filtered GBIF coordinates.

These metrics are intentionally separate from occupied-climate breadth. They are used to
ask whether any association between climatic breadth and floral-colour organization
persists after accounting for the simple fact that geographically widespread species
usually occupy more climatic conditions.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

EARTH_RADIUS_KM = 6371.0088


def read_csv(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def circular_longitude_span(lons: list[float]) -> float:
    if len(lons) <= 1:
        return 0.0
    vals = sorted((lon + 360.0) % 360.0 for lon in lons)
    gaps = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
    gaps.append((vals[0] + 360.0) - vals[-1])
    return max(0.0, 360.0 - max(gaps))


def spherical_centroid(lats: list[float], lons: list[float]) -> tuple[float, float]:
    xs = ys = zs = 0.0
    for lat, lon in zip(lats, lons):
        phi = math.radians(lat)
        lam = math.radians(lon)
        cp = math.cos(phi)
        xs += cp * math.cos(lam)
        ys += cp * math.sin(lam)
        zs += math.sin(phi)
    n = float(len(lats))
    xs /= n
    ys /= n
    zs /= n
    lon = math.degrees(math.atan2(ys, xs))
    hyp = math.hypot(xs, ys)
    lat = math.degrees(math.atan2(zs, hyp))
    return lat, lon


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dphi = p2 - p1
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    a = min(1.0, max(0.0, a))
    return 2.0 * EARTH_RADIUS_KM * math.asin(math.sqrt(a))


def quantile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    x = sorted(values)
    if len(x) == 1:
        return x[0]
    pos = (len(x) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return x[lo]
    f = pos - lo
    return x[lo] * (1 - f) + x[hi] * f


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--occurrences", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--qc-out", required=True)
    args = p.parse_args()

    rows = read_csv(args.occurrences)
    by_species: dict[str, list[tuple[float, float]]] = defaultdict(list)
    family: dict[str, str] = {}
    for row in rows:
        name = str(row.get("canonical_name") or "").strip()
        if not name:
            continue
        try:
            lat = float(row["decimalLatitude"])
            lon = float(row["decimalLongitude"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (math.isfinite(lat) and math.isfinite(lon)):
            continue
        by_species[name].append((lat, lon))
        family[name] = str(row.get("family") or "").strip()

    output = []
    for name in sorted(by_species):
        coords = by_species[name]
        lats = [x[0] for x in coords]
        lons = [x[1] for x in coords]
        c_lat, c_lon = spherical_centroid(lats, lons)
        radii = [haversine_km(c_lat, c_lon, lat, lon) for lat, lon in coords]
        lat_span = max(lats) - min(lats)
        lon_span = circular_longitude_span(lons)
        mean_lat = sum(lats) / len(lats)
        lat_km = lat_span * 111.32
        lon_km = lon_span * 111.32 * max(0.0, math.cos(math.radians(mean_lat)))
        bbox_area = lat_km * lon_km
        rms = math.sqrt(sum(r * r for r in radii) / len(radii))
        output.append({
            "canonical_name": name,
            "family": family.get(name, ""),
            "n_geographic_records": len(coords),
            "centroid_latitude": c_lat,
            "centroid_longitude": c_lon,
            "latitudinal_span_deg": lat_span,
            "minimal_longitude_span_deg": lon_span,
            "geographic_radius_median_km": quantile(radii, 0.5),
            "geographic_radius_95_km": quantile(radii, 0.95),
            "geographic_dispersion_rms_km": rms,
            "geographic_bbox_area_km2_approx": bbox_area,
        })

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fields = list(output[0].keys()) if output else ["canonical_name"]
    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)

    qc = {
        "status": "complete",
        "occurrence_rows": len(rows),
        "species_with_geographic_metrics": len(output),
        "primary_geographic_covariate": "log1p(geographic_radius_95_km)",
        "secondary_geographic_covariates": [
            "log1p(geographic_dispersion_rms_km)",
            "log1p(geographic_bbox_area_km2_approx)",
        ],
        "antimeridian_handling": "minimal circular longitude span",
        "semantic_guard": "Geographic extent is a confounder/sensitivity covariate, not a climatic niche metric.",
    }
    Path(args.qc_out).write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))

    if not output:
        raise SystemExit("No geographic extent metrics produced")


if __name__ == "__main__":
    main()
