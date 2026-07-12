#!/usr/bin/env python3
"""Build a maximal angiosperm species census from the GBIF backbone.

The builder deliberately keeps all empirical fields unknown.  Taxonomy retrieval is
separate from later flower-colour screening.  Network access is required.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

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
             retries: int = 5, timeout: int = 60) -> Dict[str, object]:
    query = "?" + urllib.parse.urlencode(params or {}, doseq=True) if params else ""
    url = API + path + query
    last_error: Optional[Exception] = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "fcp-angiosperm-census/2.0 (GitHub: zuizui0223/fcp)"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
                if not isinstance(payload, dict):
                    raise RuntimeError(f"Unexpected non-object response from {url}")
                return payload
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(min(2 ** attempt, 20))
    raise RuntimeError(f"Failed after {retries} attempts: {url}") from last_error


def candidate_keys(match: Dict[str, object]) -> List[int]:
    """Return accepted and usage keys, accepted first, without duplicates."""
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
    """Query descendants using the GBIF Species Search parameter spelling."""
    return get_json(
        "/species/search",
        {
            "highertaxon_key": root_key,
            "rank": "SPECIES",
            "status": "ACCEPTED",
            "limit": limit,
            "offset": offset,
        },
    )


def probe_root(root_key: int) -> int:
    """Return the reported descendant count, or zero when the root is unusable."""
    data = search_page(root_key, limit=1, offset=0)
    results = data.get("results", [])
    if not isinstance(results, list) or not results:
        return 0
    try:
        return int(data.get("count", len(results)))
    except (TypeError, ValueError):
        return len(results)


def resolve_angiosperm_roots() -> List[Dict[str, object]]:
    """Resolve one or more GBIF taxa whose union covers flowering plants.

    GBIF may resolve names such as Magnoliophyta to a synonym.  We therefore prefer
    acceptedUsageKey, probe every key for actual accepted-species descendants, and
    retain all productive roots.  Overlap is removed later by GBIF species key.
    """
    roots: List[Dict[str, object]] = []
    seen_keys = set()
    diagnostics: List[Dict[str, object]] = []
    for name in ROOT_NAMES:
        match = get_json("/species/match", {"name": name, "verbose": "true"})
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


def iter_species(root_keys: Sequence[int], page_size: int = 1000,
                 max_records: Optional[int] = None) -> Iterator[Dict[str, object]]:
    """Page through accepted species below all roots, deduplicating by GBIF key."""
    yielded = 0
    seen = set()
    for root_key in root_keys:
        offset = 0
        while True:
            data = search_page(root_key, limit=page_size, offset=offset)
            results = data.get("results", [])
            if not isinstance(results, list) or not results:
                break
            for row in results:
                if not isinstance(row, dict):
                    continue
                key = row.get("key") or row.get("usageKey") or row.get("scientificName")
                if key in seen:
                    continue
                seen.add(key)
                yield row
                yielded += 1
                if max_records is not None and yielded >= max_records:
                    return
            if bool(data.get("endOfRecords", False)):
                break
            offset += len(results)


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
    required = args.min_records
    if args.max_records is not None:
        required = min(required, args.max_records)
    if count < required:
        Path(args.out).unlink(missing_ok=True)
        raise RuntimeError(
            f"Angiosperm census too small: wrote {count}, required at least {required}; "
            f"roots={[(r['query_name'], r['key'], r['count']) for r in roots]}"
        )
    print(json.dumps({
        "roots": [{"query_name": r["query_name"], "key": r["key"], "count": r["count"]} for r in roots],
        "records_written": count,
        "out": args.out,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
