#!/usr/bin/env python3
"""Build a paginated, quality-filtered GBIF occurrence sample for reviewer sensitivity analyses.

The script resolves each input binomial against the GBIF backbone, pages through the
Occurrence Search API, excludes records flagged with geospatial issues, and applies
explicit coordinate and basis-of-record filters. It writes the same core columns used
by compute_climatic_niche_metrics.py.
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
from collections import Counter
from pathlib import Path
from typing import Any

MATCH_API = "https://api.gbif.org/v1/species/match"
OCCURRENCE_API = "https://api.gbif.org/v1/occurrence/search"
ACCEPTED_BASIS = {
    "PRESERVED_SPECIMEN",
    "MATERIAL_SAMPLE",
    "HUMAN_OBSERVATION",
    "MACHINE_OBSERVATION",
    "OBSERVATION",
}
OUTPUT_FIELDS = [
    "canonical_name",
    "family",
    "role",
    "focal_species",
    "match_level",
    "gbif_key",
    "gbif_taxon_key",
    "gbif_accepted_name",
    "decimalLatitude",
    "decimalLongitude",
    "coordinateUncertaintyInMeters",
    "year",
    "basisOfRecord",
    "datasetKey",
]


def read_rows(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def request_json(url: str, *, timeout: int, retries: int) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "fcp-jbi-reviewer-sensitivity/1.0 (GitHub: zuizui0223/fcp)",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"Unexpected response type: {type(payload)!r}")
            return payload
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"Request failed after {retries} attempts: {url}") from last


def resolve_taxon(name: str, *, timeout: int, retries: int) -> dict[str, Any]:
    params = urllib.parse.urlencode(
        {"name": name, "kingdom": "Plantae", "strict": "true"}
    )
    payload = request_json(f"{MATCH_API}?{params}", timeout=timeout, retries=retries)
    key = payload.get("acceptedUsageKey") or payload.get("speciesKey") or payload.get("usageKey")
    match_type = str(payload.get("matchType") or "")
    rank = str(payload.get("rank") or "")
    scientific_name = str(
        payload.get("scientificName") or payload.get("canonicalName") or ""
    )
    accepted_name = str(
        payload.get("acceptedScientificName")
        or payload.get("scientificName")
        or payload.get("canonicalName")
        or ""
    )
    status = str(payload.get("status") or "")
    if key is None or match_type not in {"EXACT", "CONFIDENCE"}:
        raise RuntimeError(
            f"No strict GBIF backbone match for {name!r}: "
            f"matchType={match_type!r}, rank={rank!r}, name={scientific_name!r}"
        )
    return {
        "taxon_key": int(key),
        "match_type": match_type,
        "rank": rank,
        "scientific_name": scientific_name,
        "accepted_name": accepted_name,
        "status": status,
    }


def finite_coordinate(record: dict[str, Any]) -> tuple[float, float] | None:
    try:
        lat = float(record.get("decimalLatitude"))
        lon = float(record.get("decimalLongitude"))
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(lat) and math.isfinite(lon)):
        return None
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        return None
    if abs(lat) < 1e-12 and abs(lon) < 1e-12:
        return None
    return lat, lon


def uncertainty_ok(record: dict[str, Any], maximum: float) -> bool:
    value = record.get("coordinateUncertaintyInMeters")
    if value in (None, ""):
        return True
    try:
        uncertainty = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(uncertainty) and 0 <= uncertainty <= maximum


def record_key(record: dict[str, Any], lat: float, lon: float, decimals: int) -> tuple[Any, ...]:
    return (
        round(lat, decimals),
        round(lon, decimals),
    )


def occurrence_pages(
    taxon_key: int,
    *,
    max_records: int,
    page_size: int,
    timeout: int,
    retries: int,
    delay: float,
) -> tuple[list[dict[str, Any]], int, int]:
    records: list[dict[str, Any]] = []
    offset = 0
    requests = 0
    reported_count = 0
    while len(records) < max_records and offset < 100_000:
        limit = min(page_size, max_records - len(records), 300)
        params = [
            ("taxonKey", str(taxon_key)),
            ("hasCoordinate", "true"),
            ("hasGeospatialIssue", "false"),
            ("occurrenceStatus", "present"),
            ("limit", str(limit)),
            ("offset", str(offset)),
        ]
        url = f"{OCCURRENCE_API}?{urllib.parse.urlencode(params)}"
        payload = request_json(url, timeout=timeout, retries=retries)
        requests += 1
        batch = payload.get("results", [])
        if not isinstance(batch, list):
            raise RuntimeError(f"Occurrence results were not a list for taxonKey={taxon_key}")
        reported_count = int(payload.get("count") or reported_count or 0)
        records.extend(x for x in batch if isinstance(x, dict))
        offset += len(batch)
        if not batch or bool(payload.get("endOfRecords")):
            break
        time.sleep(delay)
    return records[:max_records], requests, reported_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True, help="CSV containing canonical_name and family")
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--taxon-audit-out", required=True)
    parser.add_argument("--max-records", type=int, default=3000)
    parser.add_argument("--page-size", type=int, default=300)
    parser.add_argument("--dedup-decimals", type=int, default=3)
    parser.add_argument("--max-coordinate-uncertainty-m", type=float, default=20_000)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--delay", type=float, default=0.10)
    args = parser.parse_args()

    rows = read_rows(args.species)
    species: dict[str, dict[str, str]] = {}
    for row in rows:
        name = str(row.get("canonical_name") or "").strip()
        if not name:
            continue
        species[name] = {
            "canonical_name": name,
            "family": str(row.get("family") or row.get("family_class") or "").strip(),
        }
    if not species:
        raise ValueError("No species found in input")

    output: list[dict[str, Any]] = []
    taxon_audit: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    exclusion_counts: Counter[str] = Counter()
    request_count = 0
    total_reported_occurrences = 0

    for index, item in enumerate(sorted(species.values(), key=lambda x: x["canonical_name"])):
        name = item["canonical_name"]
        try:
            match = resolve_taxon(name, timeout=args.timeout, retries=args.retries)
            raw, requests, reported = occurrence_pages(
                match["taxon_key"],
                max_records=args.max_records,
                page_size=args.page_size,
                timeout=args.timeout,
                retries=args.retries,
                delay=args.delay,
            )
            request_count += requests + 1
            total_reported_occurrences += reported
        except Exception as exc:  # noqa: BLE001
            errors.append({"canonical_name": name, "error": f"{type(exc).__name__}: {exc}"[:500]})
            continue

        seen: set[tuple[Any, ...]] = set()
        retained = 0
        basis_counts: Counter[str] = Counter()
        for record in raw:
            basis = str(record.get("basisOfRecord") or "")
            if basis not in ACCEPTED_BASIS:
                exclusion_counts["basis_of_record"] += 1
                continue
            coordinate = finite_coordinate(record)
            if coordinate is None:
                exclusion_counts["invalid_coordinate"] += 1
                continue
            lat, lon = coordinate
            if not uncertainty_ok(record, args.max_coordinate_uncertainty_m):
                exclusion_counts["coordinate_uncertainty"] += 1
                continue
            key = record_key(record, lat, lon, args.dedup_decimals)
            if key in seen:
                exclusion_counts["spatial_duplicate"] += 1
                continue
            seen.add(key)
            retained += 1
            basis_counts[basis] += 1
            output.append(
                {
                    **item,
                    "role": "focal",
                    "focal_species": name,
                    "match_level": "self",
                    "gbif_key": record.get("key", ""),
                    "gbif_taxon_key": match["taxon_key"],
                    "gbif_accepted_name": match["accepted_name"],
                    "decimalLatitude": lat,
                    "decimalLongitude": lon,
                    "coordinateUncertaintyInMeters": record.get(
                        "coordinateUncertaintyInMeters", ""
                    ),
                    "year": record.get("year", ""),
                    "basisOfRecord": basis,
                    "datasetKey": record.get("datasetKey", ""),
                }
            )

        taxon_audit.append(
            {
                "canonical_name": name,
                "family": item["family"],
                **match,
                "gbif_reported_count": reported,
                "records_retrieved_before_filters": len(raw),
                "records_retained": retained,
                "basis_counts": json.dumps(dict(sorted(basis_counts.items())), sort_keys=True),
            }
        )
        if index + 1 < len(species):
            time.sleep(args.delay)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(output)

    audit_path = Path(args.taxon_audit_out)
    with audit_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "canonical_name",
            "family",
            "taxon_key",
            "match_type",
            "rank",
            "scientific_name",
            "accepted_name",
            "status",
            "gbif_reported_count",
            "records_retrieved_before_filters",
            "records_retained",
            "basis_counts",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(taxon_audit)

    per_species = Counter(str(row["canonical_name"]) for row in output)
    qc = {
        "species_requested": len(species),
        "species_resolved": len(taxon_audit),
        "species_with_retained_records": len(per_species),
        "species_ge20_retained_records": sum(value >= 20 for value in per_species.values()),
        "retained_coordinate_rows": len(output),
        "gbif_requests_including_name_resolution": request_count,
        "sum_gbif_reported_occurrence_counts": total_reported_occurrences,
        "max_records_retrieved_per_species": args.max_records,
        "page_size": min(args.page_size, 300),
        "hasGeospatialIssue_query_value": False,
        "accepted_basis_of_record": sorted(ACCEPTED_BASIS),
        "maximum_coordinate_uncertainty_m": args.max_coordinate_uncertainty_m,
        "coordinate_deduplication": (
            f"one record per rounded latitude/longitude cell at {args.dedup_decimals} decimals"
        ),
        "exclusions": dict(sorted(exclusion_counts.items())),
        "failed_species": errors,
        "method_guard": (
            "Paginated GBIF occurrence-search sensitivity sample; not a DOI-bearing GBIF "
            "download and not a complete census when a taxon exceeds the configured record cap."
        ),
    }
    Path(args.qc_out).write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
