#!/usr/bin/env python3
"""Retrieve and deduplicate all OpenAlex OQL v2 title/abstract search blocks.

This is a search-surface repair experiment for the JBI upstream re-audit. It does not
perform biological inclusion or spatial classification. The script records query-level
counts, deduplicates by normalized DOI then title+year, and audits recovery of the exact
34 historical classification sources.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

USER_AGENT = "fcp-jbi-search-v2/1.0 (https://github.com/zuizui0223/fcp)"
SPACE_RE = re.compile(r"\s+")
TAG_RE = re.compile(r"<[^>]+>")


def clean_text(value: Any) -> str:
    text = html.unescape(str(value or ""))
    text = TAG_RE.sub(" ", text)
    return SPACE_RE.sub(" ", text).strip()


def normalize_doi(value: Any) -> str:
    doi = clean_text(value).lower()
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi)
    doi = re.sub(r"^doi:\s*", "", doi)
    return doi.strip().rstrip(".")


def normalize_title(value: Any) -> str:
    title = unicodedata.normalize("NFKD", clean_text(value)).lower()
    title = "".join(ch for ch in title if not unicodedata.combining(ch))
    title = re.sub(r"[^a-z0-9\u3040-\u30ff\u3400-\u9fff]+", " ", title)
    return SPACE_RE.sub(" ", title).strip()


def normalize_openalex(value: Any) -> str:
    text = clean_text(value).rstrip("/")
    if "openalex.org/" in text.lower():
        return text.rsplit("/", 1)[-1].upper()
    if re.fullmatch(r"W\d+", text, re.I):
        return text.upper()
    return ""


def stable_key(row: dict[str, Any]) -> str:
    doi = normalize_doi(row.get("doi"))
    if doi:
        return "doi:" + doi
    return f"title:{normalize_title(row.get('title'))}|year:{row.get('year') or ''}"


def reconstruct_abstract(inverted: Any) -> str:
    if not isinstance(inverted, dict):
        return ""
    positions: list[tuple[int, str]] = []
    for word, indexes in inverted.items():
        if not isinstance(indexes, list):
            continue
        for index in indexes:
            if isinstance(index, int):
                positions.append((index, str(word)))
    positions.sort()
    return clean_text(" ".join(word for _, word in positions))


def request_json(url: str, timeout: int, retries: int) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
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


def retrieve_query(query_id: str, oql: str, api_key: str, timeout: int, retries: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    cursor = "*"
    rows: list[dict[str, Any]] = []
    reported_count: int | None = None
    pages = 0
    while True:
        params = {
            "oql": oql,
            "per-page": 100,
            "cursor": cursor,
            "select": "id,doi,title,display_name,publication_year,publication_date,type,language,cited_by_count,abstract_inverted_index,primary_location",
        }
        if api_key:
            params["api_key"] = api_key
        url = "https://api.openalex.org/?" + urllib.parse.urlencode(params)
        payload = request_json(url, timeout, retries)
        pages += 1
        meta = payload.get("meta") if isinstance(payload.get("meta"), dict) else {}
        if reported_count is None:
            try:
                reported_count = int(meta.get("count"))
            except (TypeError, ValueError):
                reported_count = None
        results = payload.get("results")
        if not isinstance(results, list) or not results:
            break
        for item in results:
            if not isinstance(item, dict):
                continue
            primary = item.get("primary_location") if isinstance(item.get("primary_location"), dict) else {}
            source = primary.get("source") if isinstance(primary.get("source"), dict) else {}
            rows.append({
                "query_id": query_id,
                "record_id": clean_text(item.get("id")),
                "doi": normalize_doi(item.get("doi")),
                "title": clean_text(item.get("title") or item.get("display_name")),
                "abstract": reconstruct_abstract(item.get("abstract_inverted_index")),
                "year": item.get("publication_year") or "",
                "publication_date": clean_text(item.get("publication_date")),
                "work_type": clean_text(item.get("type")),
                "language": clean_text(item.get("language")),
                "cited_by_count": item.get("cited_by_count") or 0,
                "journal": clean_text(source.get("display_name")),
                "url": clean_text(primary.get("landing_page_url") or item.get("id")),
            })
        next_cursor = meta.get("next_cursor")
        if not next_cursor or str(next_cursor) == cursor:
            break
        cursor = str(next_cursor)
        if pages % 20 == 0:
            print({"query_id": query_id, "pages": pages, "retrieved": len(rows), "reported": reported_count}, flush=True)
        time.sleep(0.02)

    return rows, {
        "query_id": query_id,
        "reported_count": reported_count,
        "retrieved_count": len(rows),
        "pages": pages,
        "truncated": bool(reported_count is not None and len(rows) < reported_count),
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        if fields is None:
            fields = ["empty"]
        with path.open("w", newline="", encoding="utf-8") as handle:
            csv.DictWriter(handle, fieldnames=fields).writeheader()
        return
    fields = fields or list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def historical_source_key(source_id: str) -> tuple[str, str]:
    doi = normalize_doi(source_id)
    if doi.startswith("10."):
        return "doi", doi
    oa = normalize_openalex(source_id)
    if oa:
        return "openalex", oa
    return "raw", clean_text(source_id).lower()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="literature/itv_fcp_search_config_v2.json")
    parser.add_argument("--historical-manifest", default="docs/supporting/frozen_classification_manifest.csv")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    queries = config.get("queries") or []
    if len(queries) != 15:
        raise SystemExit(f"Expected 15 v2 query blocks, found {len(queries)}")
    api_key = os.environ.get(args.api_key_env, "").strip()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw: list[dict[str, Any]] = []
    logs: list[dict[str, Any]] = []
    for item in queries:
        query_id = str(item["id"])
        rows, log = retrieve_query(query_id, str(item["oql"]), api_key, args.timeout, args.retries)
        raw.extend(rows)
        logs.append(log)
        print(log, flush=True)

    dedup: dict[str, dict[str, Any]] = {}
    query_membership: dict[str, set[str]] = {}
    for row in raw:
        key = stable_key(row)
        query_membership.setdefault(key, set()).add(str(row["query_id"]))
        previous = dedup.get(key)
        if previous is None or len(str(row.get("abstract") or "")) > len(str(previous.get("abstract") or "")):
            dedup[key] = dict(row)
    dedup_rows: list[dict[str, Any]] = []
    for key, row in sorted(dedup.items()):
        out = dict(row)
        out["dedup_key"] = key
        out["query_ids"] = ";".join(sorted(query_membership[key]))
        dedup_rows.append(out)

    with Path(args.historical_manifest).open(newline="", encoding="utf-8") as handle:
        historical = list(csv.DictReader(handle))
    if len(historical) != 34:
        raise SystemExit("Historical manifest must contain 34 rows")

    dois = {normalize_doi(row.get("doi")) for row in dedup_rows if normalize_doi(row.get("doi"))}
    oa_ids = {normalize_openalex(row.get("record_id")) for row in dedup_rows if normalize_openalex(row.get("record_id"))}
    recovery: list[dict[str, Any]] = []
    for row in historical:
        kind, value = historical_source_key(str(row.get("source_id") or ""))
        recovered = (kind == "doi" and value in dois) or (kind == "openalex" and value in oa_ids)
        recovery.append({
            "canonical_name": row.get("canonical_name", ""),
            "historical_source_id": row.get("source_id", ""),
            "source_kind": kind,
            "recovered_in_openalex_v2_direct_queries": int(recovered),
        })

    raw_path = outdir / "openalex_oql_v2_raw.csv"
    dedup_path = outdir / "openalex_oql_v2_deduplicated.csv"
    log_path = outdir / "openalex_oql_v2_query_log.csv"
    recovery_path = outdir / "openalex_oql_v2_historical34_recovery.csv"
    write_csv(raw_path, raw)
    write_csv(dedup_path, dedup_rows)
    write_csv(log_path, logs)
    write_csv(recovery_path, recovery)

    truncated = [row for row in logs if row["truncated"]]
    summary = {
        "status": "complete",
        "protocol_version": config.get("protocol_version"),
        "query_blocks": len(logs),
        "raw_query_memberships": len(raw),
        "deduplicated_works": len(dedup_rows),
        "duplicates_removed": len(raw) - len(dedup_rows),
        "truncated_v2_query_blocks": len(truncated),
        "truncated_query_ids": [row["query_id"] for row in truncated],
        "historical_34_exact_sources_recovered_by_direct_v2_queries": sum(row["recovered_in_openalex_v2_direct_queries"] for row in recovery),
        "historical_34_exact_sources_not_recovered_by_direct_v2_queries": [
            row["canonical_name"] for row in recovery if not row["recovered_in_openalex_v2_direct_queries"]
        ],
        "crossref_role": config.get("crossref_role"),
        "semantic_guard": "V2 retrieval measures search coverage only; biological inclusion and spatial coding require blinded source review.",
    }
    (outdir / "openalex_oql_v2_retrieval_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

    if truncated:
        raise SystemExit(f"V2 direct query retrieval unexpectedly truncated: {truncated}")


if __name__ == "__main__":
    main()
