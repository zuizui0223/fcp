#!/usr/bin/env python3
"""Run the systematic search with optional anonymous OpenAlex access.

OpenAlex keys are used when available. When the repository secret is absent,
this wrapper removes only the sentinel empty-key parameter and preserves all
other query parameters, allowing the same reproducible search at public limits.
"""
from __future__ import annotations

import os
import urllib.parse

import systematic_itv_fcp_search as search

SENTINEL = "__OPENALEX_ANONYMOUS__"
_original_request_json = search.request_json


def request_json_without_sentinel(url: str, **kwargs):
    parts = urllib.parse.urlsplit(url)
    query = [
        (key, value)
        for key, value in urllib.parse.parse_qsl(parts.query, keep_blank_values=True)
        if not (key == "api_key" and value == SENTINEL)
    ]
    cleaned = urllib.parse.urlunsplit(
        (parts.scheme, parts.netloc, parts.path, urllib.parse.urlencode(query), parts.fragment)
    )
    return _original_request_json(cleaned, **kwargs)


if not os.environ.get("OPENALEX_API_KEY", "").strip():
    os.environ["OPENALEX_API_KEY"] = SENTINEL
    search.request_json = request_json_without_sentinel

search.main()
