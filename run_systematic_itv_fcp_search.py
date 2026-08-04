#!/usr/bin/env python3
"""Run the systematic search with public-API compatibility fallbacks.

OpenAlex keys are used when available. Anonymous access remains possible, but
requests are paced and HTTP 429 responses honour `Retry-After` with exponential
backoff. The wrapper also removes unsupported Crossref select fields and uses the
largest supported OpenAlex page size to reduce request pressure without changing
search scope.
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

import systematic_itv_fcp_search as search

SENTINEL = "__OPENALEX_ANONYMOUS__"


def compatible_url(url: str) -> str:
    parts = urllib.parse.urlsplit(url)
    query: list[tuple[str, str]] = []
    for key, value in urllib.parse.parse_qsl(parts.query, keep_blank_values=True):
        if key == "api_key" and value == SENTINEL:
            continue
        if parts.netloc == "api.crossref.org" and key == "select":
            selected = [field for field in value.split(",") if field and field != "language"]
            value = ",".join(selected)
        if parts.netloc == "api.openalex.org" and key == "per_page" and value == "100":
            value = "200"
        query.append((key, value))
    return urllib.parse.urlunsplit(
        (parts.scheme, parts.netloc, parts.path, urllib.parse.urlencode(query), parts.fragment)
    )


def request_json_rate_limited(
    url: str,
    *,
    timeout: int,
    retries: int,
    headers: dict[str, str] | None = None,
    delay: float = 0.0,
):
    cleaned = compatible_url(url)
    parts = urllib.parse.urlsplit(cleaned)
    merged_headers = {"User-Agent": search.USER_AGENT, "Accept": "application/json"}
    if headers:
        merged_headers.update(headers)

    attempts = max(retries, 8)
    last_error: Exception | None = None
    for attempt in range(attempts):
        host_delay = 0.35 if parts.netloc == "api.openalex.org" else 0.15
        time.sleep(max(delay, host_delay))
        try:
            request = urllib.request.Request(cleaned, headers=merged_headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"Expected JSON object from {cleaned}")
            return payload
        except urllib.error.HTTPError as exc:
            last_error = exc
            body = exc.read().decode("utf-8", errors="replace")[:500]
            if exc.code in {400, 401, 403, 404}:
                raise RuntimeError(f"HTTP {exc.code} for {cleaned}: {body}") from exc
            if exc.code == 429:
                retry_after = exc.headers.get("Retry-After", "")
                try:
                    wait_seconds = float(retry_after)
                except (TypeError, ValueError):
                    wait_seconds = min(180.0, 20.0 * (2 ** attempt))
                time.sleep(max(10.0, wait_seconds))
                continue
            if 500 <= exc.code < 600:
                time.sleep(min(120.0, 3.0 * (2 ** attempt)))
                continue
            raise RuntimeError(f"HTTP {exc.code} for {cleaned}: {body}") from exc
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(min(120.0, 2.0 * (2 ** attempt)))

    raise RuntimeError(f"Request failed after {attempts} attempts: {cleaned}") from last_error


if not os.environ.get("OPENALEX_API_KEY", "").strip():
    os.environ["OPENALEX_API_KEY"] = SENTINEL

search.request_json = request_json_rate_limited
search.main()
