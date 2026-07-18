#!/usr/bin/env python3
"""Build reproducible geography covariates from GBIF occurrence search results.

This script intentionally avoids inferring island status from country names. It derives only
quantities directly supported by coordinate-bearing GBIF occurrences: occurrence count,
median absolute latitude, and robust latitudinal range. Results are checkpointed after every
species so interrupted CI runs retain completed work.
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
FIELDS = [
    "canonical_name",
    "family",
    "gbif_occurrences",
    "gbif_coordinate_records_sampled",
    "absolute_latitude",
    "median_latitude",
    "latitudinal_range",
    "latitude_q05",
    "latitude_q95",
    "gbif_query_status",
]


def read_rows(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
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
                    "User-Agent": "fcp-gbif-geography/1.0 (research workflow)",
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


def summarize(name: str, family: str, payload: dict) -> dict[str, object]:
    total = int(payload.get("count") or 0)
    latitudes: list[float] = []
    for item in payload.get("results") or []:
        if not isinstance(item, dict):
            continue
        value = item.get("decimalLatitude")
        try:
            lat = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(lat) and -90 <= lat <= 90:
            latitudes.append(lat)
    if latitudes:
        arr = np.asarray(latitudes, dtype=float)
        q05, median, q95 = np.quantile(arr, [0.05, 0.5, 0.95])
        absolute_latitude = abs(float(median))
        latitudinal_range = float(q95 - q05)
        status = "complete"
    else:
        q05 = median = q95 = absolute_latitude = latitudinal_range = ""
        status = "no_coordinate_records"
    return {
        "canonical_name": name,
        "family": family,
        "gbif_occurrences": total,
        "gbif_coordinate_records_sampled": len(latitudes),
        "absolute_latitude": absolute_latitude,
        "median_latitude": median,
        "latitudinal_range": latitudinal_range,
        "latitude_q05": q05,
        "latitude_q95": q95,
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
        if name in completed and completed[name].get("gbif_query_status") in {"complete", "no_coordinate_records"}:
            result = completed[name]
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
                result = {
                    "canonical_name": name,
                    "family": family,
                    "gbif_occurrences": "",
                    "gbif_coordinate_records_sampled": 0,
                    "absolute_latitude": "",
                    "median_latitude": "",
                    "latitudinal_range": "",
                    "latitude_q05": "",
                    "latitude_q95": "",
                    "gbif_query_status": "request_failed",
                }
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
            "method": "GBIF occurrence search; median latitude and 5th-95th percentile range from up to 300 coordinate records per species",
        }
        Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(qc, ensure_ascii=False))


if __name__ == "__main__":
    main()
