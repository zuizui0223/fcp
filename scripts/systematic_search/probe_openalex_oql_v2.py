#!/usr/bin/env python3
"""Probe OpenAlex OQL v2 title/abstract query counts without downloading result sets.

This is a search-completeness diagnostic only. It issues one small request per query,
records OpenAlex's reported hit count, and compares the four legacy OpenAlex shards that
were truncated under broad title+abstract+full-text search.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

USER_AGENT = "fcp-jbi-search-v2-probe/1.0 (https://github.com/zuizui0223/fcp)"
LEGACY_TRUNCATED_OPENALEX = {
    "floral_reflectance": 411733,
    "japanese_chinese": 11831,
    "orchid_labellum": 26105,
    "pigment_variation": 50201,
}


def request_json(url: str, timeout: int, retries: int) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError("OpenAlex response was not a JSON object")
            return payload
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(url) from last


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="literature/itv_fcp_search_config_v2.json")
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    queries = config.get("queries") or []
    if len(queries) != 15:
        raise SystemExit(f"Expected 15 v2 query blocks, found {len(queries)}")

    api_key = os.environ.get(args.api_key_env, "").strip()
    rows: list[dict[str, Any]] = []
    for item in queries:
        query_id = str(item["id"])
        oql = str(item["oql"])
        params = {
            "oql": oql,
            "per-page": 1,
        }
        if api_key:
            params["api_key"] = api_key
        url = "https://api.openalex.org/?" + urllib.parse.urlencode(params)
        payload = request_json(url, args.timeout, args.retries)
        meta = payload.get("meta") if isinstance(payload.get("meta"), dict) else {}
        try:
            count = int(meta.get("count"))
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"No valid OpenAlex count for {query_id}: {meta}") from exc
        legacy = LEGACY_TRUNCATED_OPENALEX.get(query_id)
        rows.append({
            "query_id": query_id,
            "v2_title_abstract_count": count,
            "legacy_openalex_broad_count": legacy if legacy is not None else "",
            "legacy_was_truncated_openalex": int(query_id in LEGACY_TRUNCATED_OPENALEX),
            "v2_to_legacy_ratio": (count / legacy) if legacy else "",
            "oql": oql,
        })
        print({"query_id": query_id, "v2_title_abstract_count": count}, flush=True)
        time.sleep(0.05)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    truncated_comparison = [row for row in rows if row["legacy_was_truncated_openalex"]]
    summary = {
        "status": "complete",
        "protocol_version": config.get("protocol_version"),
        "query_blocks": len(rows),
        "legacy_truncated_openalex_blocks": len(truncated_comparison),
        "v2_counts_for_legacy_truncated_openalex": {
            row["query_id"]: int(row["v2_title_abstract_count"])
            for row in truncated_comparison
        },
        "all_v2_counts": {row["query_id"]: int(row["v2_title_abstract_count"]) for row in rows},
        "crossref_role": config.get("crossref_role"),
        "semantic_guard": "Reported counts measure retrieval scope only and cannot be used as biological evidence.",
    }
    Path(args.summary_out).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
