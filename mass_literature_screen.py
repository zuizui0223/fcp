#!/usr/bin/env python3
"""Fast first-pass discovery of flower-colour polymorphism literature.

OpenAlex currently allows only about 100 search calls/day without an API key and
about 1,000/day with a free key. This script therefore refuses large unauthenticated
batches before making requests, includes an API key when available, and reports
quota/rate-limit failures explicitly instead of silently consuming the whole batch.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

OPENALEX = "https://api.openalex.org/works"
UNAUTHENTICATED_SAFE_LIMIT = 90

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


class OpenAlexQuotaError(RuntimeError):
    """OpenAlex rejected the request because of rate or daily-budget limits."""


@dataclass(frozen=True)
class Hit:
    title: str
    year: str
    doi: str
    landing_url: str
    source: str
    query: str


class RequestThrottle:
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
                    "User-Agent": "fcp-mass-screen/3.0 (GitHub: zuizui0223/fcp)",
                    "Accept": "application/json",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"Unexpected response type: {type(payload)!r}")
            return payload
        except urllib.error.HTTPError as exc:
            body = ""
            try:
                body = exc.read().decode("utf-8", errors="replace")[:500]
            except Exception:
                pass
            retry_after = exc.headers.get("Retry-After", "") if exc.headers else ""
            if exc.code in {401, 402, 403, 429}:
                raise OpenAlexQuotaError(
                    f"OpenAlex HTTP {exc.code}; retry_after={retry_after!r}; body={body!r}"
                ) from exc
            last = RuntimeError(f"OpenAlex HTTP {exc.code}; body={body!r}")
        except Exception as exc:  # noqa: BLE001
            last = exc
        if attempt + 1 < retries:
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"request failed after {retries} attempts: {url}") from last


def combined_query(name: str) -> str:
    return f'"{name}" flower floral color colour polymorphism morph variation petal'


def search_openalex(
    query: str,
    *,
    api_key: str,
    per_page: int,
    throttle: RequestThrottle,
    retries: int,
    timeout: int,
) -> List[Hit]:
    params: Dict[str, object] = {"search": query, "per-page": per_page}
    if api_key:
        params["api_key"] = api_key
    data = get_json(
        f"{OPENALEX}?{urllib.parse.urlencode(params)}",
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
    score = 5 if species.lower() in text else 0
    score += 5 * sum(term in text for term in HIGH_TERMS)
    score += 2 * sum(term in text for term in MEDIUM_TERMS)
    if "cultivar" in text or "horticultur" in text:
        score -= 3
    if "ontogenetic" in text or "color change" in text or "colour change" in text:
        score -= 2
    return score


def iter_census(path: Path, start: int, limit: int | None) -> Iterable[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        for idx, row in enumerate(csv.DictReader(handle)):
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
        if key and key not in seen:
            seen.add(key)
            out.append(hit)
    return out


def screen_species(
    row: Dict[str, str],
    *,
    api_key: str,
    per_page: int,
    throttle: RequestThrottle,
    retries: int,
    timeout: int,
) -> Tuple[Dict[str, str], str, List[Hit], str | None, bool]:
    name = (row.get("canonical_name") or row.get("scientific_name") or "").strip()
    if not name:
        return row, "", [], "missing_name", False
    query = combined_query(name)
    try:
        hits = search_openalex(
            query,
            api_key=api_key,
            per_page=per_page,
            throttle=throttle,
            retries=retries,
            timeout=timeout,
        )
        return row, query, hits, None, False
    except OpenAlexQuotaError as exc:
        return row, query, [], str(exc), True
    except Exception as exc:  # noqa: BLE001
        return row, query, [], f"{type(exc).__name__}: {exc}", False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--out", default="data/mass_screen_hits.csv")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--requests-per-second", type=float, default=4.0)
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--per-page", type=int, default=10)
    parser.add_argument("--max-failure-rate", type=float, default=0.20)
    parser.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    args = parser.parse_args()

    rows = list(iter_census(Path(args.census), args.start, args.limit))
    if not rows:
        raise RuntimeError("No census rows selected")

    api_key = os.environ.get(args.api_key_env, "").strip()
    if not api_key and len(rows) > UNAUTHENTICATED_SAFE_LIMIT:
        raise RuntimeError(
            f"OpenAlex API key is required for {len(rows)} searches. Without a key, "
            f"this script permits at most {UNAUTHENTICATED_SAFE_LIMIT} searches. "
            f"Add repository secret {args.api_key_env} or lower --limit."
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    throttle = RequestThrottle(args.requests_per_second)
    processed = emitted = failures = 0
    started = time.monotonic()
    first_error = ""

    with out_path.open("w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=FIELDS)
        writer.writeheader()
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futures = [
                pool.submit(
                    screen_species,
                    row,
                    api_key=api_key,
                    per_page=args.per_page,
                    throttle=throttle,
                    retries=args.retries,
                    timeout=args.timeout,
                )
                for row in rows
            ]
            quota_error = ""
            for future in as_completed(futures):
                row, query, hits, error, is_quota = future.result()
                processed += 1
                if error:
                    failures += 1
                    if not first_error:
                        first_error = error
                    if is_quota:
                        quota_error = error
                        for pending in futures:
                            pending.cancel()
                        break
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
                        "authenticated": bool(api_key),
                        "first_error": first_error[:300],
                        "elapsed_seconds": round(time.monotonic() - started, 1),
                    }), flush=True)
            if quota_error:
                raise OpenAlexQuotaError(
                    "OpenAlex quota/rate limit reached; stopped immediately. " + quota_error
                )

    failure_rate = failures / processed if processed else 1.0
    summary = {
        "processed_species": processed,
        "emitted_hits": emitted,
        "request_failures": failures,
        "failure_rate": round(failure_rate, 4),
        "authenticated": bool(api_key),
        "first_error": first_error[:500],
        "elapsed_seconds": round(time.monotonic() - started, 1),
        "out": str(out_path),
    }
    print(json.dumps(summary), flush=True)
    if failure_rate > args.max_failure_rate:
        raise RuntimeError(
            f"OpenAlex failure rate {failure_rate:.1%} exceeded allowed "
            f"{args.max_failure_rate:.1%}; first_error={first_error!r}"
        )


if __name__ == "__main__":
    main()
