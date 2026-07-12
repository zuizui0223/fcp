#!/usr/bin/env python3
"""Build a maximal angiosperm species census from the GBIF taxonomic backbone.

This script is deliberately census-first: it does not preselect species based on
known flower-colour polymorphism. Every accepted species returned under the
resolved angiosperm root is written with unknown empirical fields that can be
filled later without changing the sampling universe.

Network access is required when the script is run.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional

API = "https://api.gbif.org/v1"
FIELDS = [
    "gbif_key", "scientific_name", "canonical_name", "family", "genus",
    "order", "class", "phylum", "kingdom", "taxonomic_status",
    "pollination_mode", "display_opportunity", "outcome_class",
    "evidence_status", "sampling_effort", "notes",
]


def get_json(path: str, params: Optional[Dict[str, object]] = None,
             retries: int = 5, timeout: int = 60) -> Dict[str, object]:
    query = "?" + urllib.parse.urlencode(params or {}, doseq=True) if params else ""
    url = API + path + query
    last_error: Optional[Exception] = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "fcp-angiosperm-census/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.load(response)
        except Exception as exc:  # noqa: BLE001 - retry transient HTTP/network failures
            last_error = exc
            time.sleep(min(2 ** attempt, 20))
    raise RuntimeError(f"Failed after {retries} attempts: {url}") from last_error


def resolve_angiosperm_root() -> Dict[str, object]:
    """Resolve an angiosperm backbone taxon, trying common names in order."""
    candidates = ["Magnoliophyta", "Angiospermae"]
    for name in candidates:
        data = get_json("/species/match", {"name": name, "verbose": "true"})
        if data.get("usageKey") and str(data.get("matchType", "")).upper() != "NONE":
            return data
    raise RuntimeError("Could not resolve an angiosperm root in the GBIF backbone")


def iter_species(root_key: int, page_size: int = 1000,
                 max_records: Optional[int] = None) -> Iterator[Dict[str, object]]:
    """Page through accepted species below the resolved root."""
    offset = 0
    yielded = 0
    while True:
        data = get_json(
            "/species/search",
            {
                "highertaxon_key": root_key,
                "rank": "SPECIES",
                "status": "ACCEPTED",
                "limit": page_size,
                "offset": offset,
            },
        )
        results = data.get("results", [])
        if not isinstance(results, list) or not results:
            break
        for row in results:
            if not isinstance(row, dict):
                continue
            yield row
            yielded += 1
            if max_records is not None and yielded >= max_records:
                return
        if bool(data.get("endOfRecords", False)):
            break
        offset += len(results)


def normalize(row: Dict[str, object]) -> Dict[str, object]:
    return {
        "gbif_key": row.get("key", ""),
        "scientific_name": row.get("scientificName", ""),
        "canonical_name": row.get("canonicalName", ""),
        "family": row.get("family", ""),
        "genus": row.get("genus", ""),
        "order": row.get("order", ""),
        "class": row.get("class", ""),
        "phylum": row.get("phylum", ""),
        "kingdom": row.get("kingdom", ""),
        "taxonomic_status": row.get("taxonomicStatus", row.get("status", "")),
        "pollination_mode": "unknown",
        "display_opportunity": "unknown",
        "outcome_class": "unknown",
        "evidence_status": "unreviewed",
        "sampling_effort": "unknown",
        "notes": "",
    }


def write_census(rows: Iterable[Dict[str, object]], out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    seen = set()
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for raw in rows:
            row = normalize(raw)
            key = row["gbif_key"] or row["scientific_name"]
            if key in seen:
                continue
            seen.add(key)
            writer.writerow(row)
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="data/angiosperm_census.csv")
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--max-records", type=int, default=None,
                        help="Optional smoke-test cap; omit for the full census")
    args = parser.parse_args()

    root = resolve_angiosperm_root()
    root_key = int(root["usageKey"])
    count = write_census(
        iter_species(root_key, page_size=args.page_size, max_records=args.max_records),
        Path(args.out),
    )
    print(json.dumps({"root": root, "records_written": count, "out": args.out}, ensure_ascii=False))


if __name__ == "__main__":
    main()
