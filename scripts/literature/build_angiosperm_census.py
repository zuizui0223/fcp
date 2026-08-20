#!/usr/bin/env python3
"""Build a maximal angiosperm species census from the GBIF backbone.

GBIF Species Search becomes unreliable at very deep offsets.  This builder avoids
that by traversing the hierarchy only until FAMILY rank, then paging accepted
species separately within each family.  No angiosperm family is expected to
approach the global 100,000-offset failure seen for the whole clade.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from collections import deque
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence

API = "https://api.gbif.org/v1"
FIELDS = [
    "gbif_key", "scientific_name", "canonical_name", "family", "genus",
    "order", "class", "phylum", "kingdom", "taxonomic_status",
    "pollination_mode", "display_opportunity", "outcome_class",
    "evidence_status", "sampling_effort", "notes",
]
ROOT_NAMES = (
    "Angiospermae",
    "Magnoliophyta",
    "Magnoliopsida",
    "Liliopsida",
    "Magnoliidae",
)


def get_json(path: str, params: Optional[Dict[str, object]] = None,
             retries: int = 6, timeout: int = 90):
    query = "?" + urllib.parse.urlencode(params or {}, doseq=True) if params else ""
    url = API + path + query
    last_error: Optional[Exception] = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "fcp-angiosperm-census/4.0 (GitHub: zuizui0223/fcp)"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.load(response)
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"Failed after {retries} attempts: {url}") from last_error


def candidate_keys(match: Dict[str, object]) -> List[int]:
    keys: List[int] = []
    for field in ("acceptedUsageKey", "usageKey"):
        value = match.get(field)
        if value is None:
            continue
        try:
            key = int(value)
        except (TypeError, ValueError):
            continue
        if key not in keys:
            keys.append(key)
    return keys


def search_page(root_key: int, limit: int, offset: int) -> Dict[str, object]:
    payload = get_json(
        "/species/search",
        {
            "highertaxonKey": root_key,
            "rank": "SPECIES",
            "status": "ACCEPTED",
            "limit": limit,
            "offset": offset,
        },
    )
    if not isinstance(payload, dict):
        raise RuntimeError(f"Unexpected Species Search response type: {type(payload)!r}")
    return payload


def probe_root(root_key: int) -> int:
    data = search_page(root_key, limit=1, offset=0)
    results = data.get("results", [])
    if not isinstance(results, list) or not results:
        return 0
    try:
        return int(data.get("count", len(results)))
    except (TypeError, ValueError):
        return len(results)


def resolve_angiosperm_roots() -> List[Dict[str, object]]:
    roots: List[Dict[str, object]] = []
    seen_keys = set()
    diagnostics: List[Dict[str, object]] = []
    for name in ROOT_NAMES:
        match = get_json("/species/match", {"name": name, "verbose": "true"})
        if not isinstance(match, dict):
            diagnostics.append({"name": name, "match": "invalid_response"})
            continue
        if str(match.get("matchType", "")).upper() == "NONE":
            diagnostics.append({"name": name, "match": "NONE"})
            continue
        for key in candidate_keys(match):
            if key in seen_keys:
                continue
            count = probe_root(key)
            diagnostics.append({"name": name, "key": key, "descendant_count": count})
            if count > 0:
                seen_keys.add(key)
                roots.append({"query_name": name, "key": key, "count": count, "match": match})
    if not roots:
        raise RuntimeError(
            "Could not resolve any productive angiosperm root in GBIF; diagnostics="
            + json.dumps(diagnostics, ensure_ascii=False)
        )
    return roots


def iter_children(parent_key: int, page_size: int = 1000) -> Iterator[Dict[str, object]]:
    offset = 0
    while True:
        payload = get_json(
            f"/species/{parent_key}/children",
            {"limit": page_size, "offset": offset},
        )
        if isinstance(payload, dict):
            results = payload.get("results", [])
            end = bool(payload.get("endOfRecords", False))
        elif isinstance(payload, list):
            results = payload
            end = len(results) < page_size
        else:
            raise RuntimeError(
                f"Unexpected children response for parent {parent_key}: {type(payload)!r}"
            )
        if not isinstance(results, list) or not results:
            break
        for row in results:
            if isinstance(row, dict):
                yield row
        if end or len(results) < page_size:
            break
        offset += len(results)


def collect_family_keys(root_keys: Sequence[int], page_size: int = 1000) -> List[int]:
    """Traverse only high taxonomy and return unique angiosperm family keys."""
    queue = deque(int(key) for key in root_keys)
    visited = set()
    families = set()

    while queue:
        parent_key = queue.popleft()
        if parent_key in visited:
            continue
        visited.add(parent_key)

        for row in iter_children(parent_key, page_size=page_size):
            raw_key = row.get("key", row.get("usageKey"))
            try:
                key = int(raw_key)
            except (TypeError, ValueError):
                continue
            rank = str(row.get("rank", "")).upper()
            if rank == "FAMILY":
                families.add(key)
            elif rank not in {"GENUS", "SPECIES", "SUBSPECIES", "VARIETY", "FORM"}:
                queue.append(key)

    if not families:
        raise RuntimeError("No angiosperm family keys were discovered from resolved roots")
    return sorted(families)


def iter_species_in_family(family_key: int, page_size: int = 1000) -> Iterator[Dict[str, object]]:
    offset = 0
    while True:
        data = search_page(family_key, limit=page_size, offset=offset)
        results = data.get("results", [])
        if not isinstance(results, list) or not results:
            break
        for row in results:
            if isinstance(row, dict):
                yield row
        if bool(data.get("endOfRecords", False)) or len(results) < page_size:
            break
        offset += len(results)
        if offset >= 100000:
            raise RuntimeError(
                f"Family {family_key} exceeded safe GBIF paging depth at offset {offset}"
            )


def iter_species(root_keys: Sequence[int], page_size: int = 1000,
                 max_records: Optional[int] = None) -> Iterator[Dict[str, object]]:
    family_keys = collect_family_keys(root_keys, page_size=page_size)
    seen_species = set()
    yielded = 0

    print(json.dumps({"families_discovered": len(family_keys)}, ensure_ascii=False), flush=True)

    for index, family_key in enumerate(family_keys, start=1):
        for row in iter_species_in_family(family_key, page_size=page_size):
            raw_key = row.get("key", row.get("usageKey"))
            try:
                key = int(raw_key)
            except (TypeError, ValueError):
                continue
            if key in seen_species:
                continue
            seen_species.add(key)
            yield row
            yielded += 1
            if max_records is not None and yielded >= max_records:
                return
        if index % 25 == 0:
            print(
                json.dumps(
                    {"families_processed": index, "families_total": len(family_keys), "species_yielded": yielded},
                    ensure_ascii=False,
                ),
                flush=True,
            )


def normalize(row: Dict[str, object]) -> Dict[str, object]:
    return {
        "gbif_key": row.get("key", row.get("usageKey", "")),
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
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for raw in rows:
            writer.writerow(normalize(raw))
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="data/angiosperm_census.csv")
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--max-records", type=int, default=None,
                        help="Optional smoke-test cap; omit for the full census")
    parser.add_argument("--min-records", type=int, default=10000,
                        help="Fail if fewer rows are written (capped by --max-records)")
    args = parser.parse_args()

    roots = resolve_angiosperm_roots()
    root_keys = [int(root["key"]) for root in roots]
    count = write_census(
        iter_species(root_keys, page_size=args.page_size, max_records=args.max_records),
        Path(args.out),
    )
    required = min(args.min_records, args.max_records) if args.max_records is not None else args.min_records
    if count < required:
        Path(args.out).unlink(missing_ok=True)
        raise RuntimeError(
            f"Angiosperm census too small: wrote {count}, required at least {required}; "
            f"roots={[(r['query_name'], r['key'], r['count']) for r in roots]}"
        )

    print(json.dumps({
        "retrieval_mode": "family_partitioned_species_search",
        "roots": [{"query_name": r["query_name"], "key": r["key"], "count": r["count"]} for r in roots],
        "records_written": count,
        "out": args.out,
    }, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
