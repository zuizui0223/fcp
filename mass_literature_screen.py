#!/usr/bin/env python3
"""Fast first-pass discovery of flower-colour polymorphism literature.

The previous implementation issued 12,000 sequential HTTP requests for a
1,000-species chunk (6 queries x OpenAlex + Crossref).  A single slow request
could consume 45 s x 4 retries, while failures were silently discarded.  This
version uses one combined query per species, bounded concurrency, short retries,
and explicit progress/failure accounting.  Crossref is intentionally reserved
for a later candidate-validation pass rather than duplicated across the entire
angiosperm census.
"""
from __future__ import annotations

import argparse
import csv
import json
import threading
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

OPENALEX = "https://api.openalex.org/works"

HIGH_TERMS = (
    "flower color polymorphism", "flower colour polymorphism",
    "floral color polymorphism", "floral colour polymorphism",
    "color morph", "colour morph", "petal color", "petal colour",
)
MEDIUM_TERMS = (
    "flower color variation", "flower colour variation",
    "floral color variation", "floral colour variation",
    "geographic variation", "frequency-dependent", "pollinator preference",
)
FIELDS = [
    "canonical_name", "family", "query", "source", "title", "year", "doi",
    "landing_url", "candidate_score", "evidence_status", "notes",
]


@dataclass(frozen=True)
class Hit:
    title: str
    year: str
    doi: str
    landing_url: str
    source: str
    query: str


class RequestThrottle:
    """Thread-safe minimum spacing between outbound requests."""

    def __init__(self, requests_per_second: float) -> None:
        self.interval = 0.0 if requests_per_second <= 0 else 1.0 / requests_per_second
        self.lock = threading.Lock()
        self.next_at = 0.0

    def wait(self) -> None:
        if self.interval <= 0:
            return
        with self.lock:
            now = time.monotonic()
            delay = max(0.0, self.next_at - now)
            self.next_at = max(now, self.next_at) + self.interval
        if delay:
            time.sleep(delay)


def get_json(
    url: str,
    *,
    throttle: RequestThrottle,
    retries: int,
    timeout: int,
) -> Dict[str, object]:
    last: Exception | None = None
    for attempt in range(retries):
        throttle.wait()
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "fcp-mass-screen/2.0 (GitHub: zuizui0223/fcp)",
                    "Accept": "application/json",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"Unexpected response type: {type(payload)!r}")
            return payload
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt + 1 < retries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"request failed after {retries} attempts: {url}") from last


def combined_query(name: str) -> str:
    # OpenAlex `search` is relevance ranked. One species-specific query replaces
    # six nearly duplicate phrase searches in the census-scale discovery pass.
    return f'"{name}" flower floral color colour polymorphism morph variation petal'


def search_openalex(
    query: str,
    *,
    per_page: int,
    throttle: RequestThrottle,
    retries: int,
    timeout: int,
) -> List[Hit]:
    params = urllib.parse.urlencode({"search": query, "per-page": per_page})
    data = get_json(
        f"{OPENALEX}?{params}",
        throttle=throttle,
        retries=retries,
        timeout=timeout,
    )
    hits: List[Hit] = []
    results = data.get("results", [])
    if not isinstance(results, list):
        return hits
    for item in results:
        if not isinstance(item, dict):
            continue
        title = str(item.get("title") or "")
        year = str(item.get("publication_year") or "")
        doi = str(item.get("doi") or "").replace("https://doi.org/", "")
        url = str(item.get("id") or "")
        primary = item.get("primary_location") or {}
        if isinstance(primary, dict):
            url = str(primary.get("landing_page_url") or url)
        hits.append(Hit(title, year, doi, url, "openalex", query))
    return hits


def score_hit(species: str, title: str) -> int:
    text = title.lower()
    score = 0
    if species.lower() in text:
        score += 5
    score += 5 * sum(term in text for term in HIGH_TERMS)
    score += 2 * sum(term in text for term in MEDIUM_TERMS)
    if "cultivar" in text or "horticultur" in text:
        score -= 3
    if "ontogenetic" in text or "color change" in text or "colour change" in text:
        score -= 2
    return score


def iter_census(path: Path, start: int, limit: int | None) -> Iterable[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        for idx, row in enumerate(rows):
            if idx < start:
                continue
            if limit is not None and idx >= start + limit:
                break
            yield row


def dedupe_hits(hits: Sequence[Hit]) -> List[Hit]:
    seen = set()
    out: List[Hit] = []
    for hit in hits:
        key = hit.doi.lower().strip() if hit.doi else hit.title.lower().strip()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(hit)
    return out


def screen_species(
    row: Dict[str, str],
    *,
    per_page: int,
    throttle: RequestThrottle,
    retries: int,
    timeout: int,
) -> Tuple[Dict[str, str], str, List[Hit], str | None]:
    name = (row.get("canonical_name") or row.get("scientific_name") or "").strip()
    if not name:
        return row, "", [], "missing_name"
    query = combined_query(name)
    try:
        hits = search_openalex(
            query,
            per_page=per_page,
            throttle=throttle,
            retries=retries,
            timeout=timeout,
        )
        return row, query, hits, None
    except Exception as exc:  # noqa: BLE001
        return row, query, [], f"{type(exc).__name__}: {exc}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--out", default="data/mass_screen_hits.csv")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None)
    # Kept for workflow compatibility; census-scale first pass is deliberately one query.
    parser.add_argument("--queries-per-species", type=int, default=1)
    parser.add_argument("--sleep", type=float, default=0.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--requests-per-second", type=float, default=8.0)
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--per-page", type=int, default=10)
    parser.add_argument("--max-failure-rate", type=float, default=0.20)
    args = parser.parse_args()

    rows = list(iter_census(Path(args.census), args.start, args.limit))
    if not rows:
        raise RuntimeError("No census rows selected; refusing to write an empty successful batch")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    throttle = RequestThrottle(args.requests_per_second)
    processed = emitted = failures = 0
    started = time.monotonic()

    with out_path.open("w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=FIELDS)
        writer.writeheader()

        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futures = [
                pool.submit(
                    screen_species,
                    row,
                    per_page=args.per_page,
                    throttle=throttle,
                    retries=args.retries,
                    timeout=args.timeout,
                )
                for row in rows
            ]
            for future in as_completed(futures):
                row, query, hits, error = future.result()
                processed += 1
                if error:
                    failures += 1
                else:
                    name = (row.get("canonical_name") or row.get("scientific_name") or "").strip()
                    family = (row.get("family") or "").strip()
                    for hit in dedupe_hits(hits):
                        score = score_hit(name, hit.title)
                        if score <= 0:
                            continue
                        writer.writerow({
                            "canonical_name": name,
                            "family": family,
                            "query": query,
                            "source": hit.source,
                            "title": hit.title,
                            "year": hit.year,
                            "doi": hit.doi,
                            "landing_url": hit.landing_url,
                            "candidate_score": score,
                            "evidence_status": "machine_candidate",
                            "notes": "First-pass OpenAlex candidate; requires primary-source review.",
                        })
                        emitted += 1
                if processed % 25 == 0 or processed == len(rows):
                    out_handle.flush()
                    print(json.dumps({
                        "processed_species": processed,
                        "total_species": len(rows),
                        "emitted_hits": emitted,
                        "request_failures": failures,
                        "elapsed_seconds": round(time.monotonic() - started, 1),
                    }), flush=True)

    failure_rate = failures / processed if processed else 1.0
    summary = {
        "processed_species": processed,
        "emitted_hits": emitted,
        "request_failures": failures,
        "failure_rate": round(failure_rate, 4),
        "elapsed_seconds": round(time.monotonic() - started, 1),
        "out": str(out_path),
        "mode": "one_openalex_query_per_species_concurrent",
    }
    print(json.dumps(summary), flush=True)
    if failure_rate > args.max_failure_rate:
        raise RuntimeError(
            f"OpenAlex failure rate {failure_rate:.1%} exceeded allowed "
            f"{args.max_failure_rate:.1%}; batch is not trustworthy"
        )


if __name__ == "__main__":
    main()
