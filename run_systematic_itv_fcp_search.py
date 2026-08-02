#!/usr/bin/env python3
"""Run the systematic search with public-API compatibility fallbacks.

OpenAlex keys are used when available. When the repository secret is absent,
the wrapper removes only the sentinel key parameter so the same search can use
anonymous public access. It also removes unsupported fields from Crossref's
`select` parameter without changing the query or pagination logic.
"""
from __future__ import annotations

import os
import urllib.parse

import systematic_itv_fcp_search as search

SENTINEL = "__OPENALEX_ANONYMOUS__"
_original_request_json = search.request_json


def request_json_with_api_compatibility(url: str, **kwargs):
    parts = urllib.parse.urlsplit(url)
    query: list[tuple[str, str]] = []
    for key, value in urllib.parse.parse_qsl(parts.query, keep_blank_values=True):
        if key == "api_key" and value == SENTINEL:
            continue
        if parts.netloc == "api.crossref.org" and key == "select":
            selected = [field for field in value.split(",") if field and field != "language"]
            value = ",".join(selected)
        query.append((key, value))
    cleaned = urllib.parse.urlunsplit(
        (parts.scheme, parts.netloc, parts.path, urllib.parse.urlencode(query), parts.fragment)
    )
    return _original_request_json(cleaned, **kwargs)


if not os.environ.get("OPENALEX_API_KEY", "").strip():
    os.environ["OPENALEX_API_KEY"] = SENTINEL

search.request_json = request_json_with_api_compatibility
search.main()
