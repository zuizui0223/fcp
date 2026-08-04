#!/usr/bin/env python3
"""Reproducible multi-database search for natural floral-colour ITV and strict FCP.

The retrieval stage is deliberately broad. Deterministic screening creates an
ordered human-review queue but never promotes a record to biological inclusion.
Final eligibility and spatial-scale coding require source-level human review.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

USER_AGENT = "fcp-systematic-itv-fcp-search/1.0 (https://github.com/zuizui0223/fcp)"
BINOMIAL_RE = re.compile(r"\b([A-Z][a-z-]{2,})\s+([a-z][a-z-]{2,})\b")
TAG_RE = re.compile(r"<[^>]+>")
SPACE_RE = re.compile(r"\s+")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


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


def stable_key(record: dict[str, Any]) -> str:
    doi = normalize_doi(record.get("doi"))
    if doi:
        return f"doi:{doi}"
    title = normalize_title(record.get("title"))
    year = str(record.get("year") or "")
    return f"title:{title}|year:{year}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def request_json(
    url: str,
    *,
    timeout: int,
    retries: int,
    headers: dict[str, str] | None = None,
    delay: float = 0.0,
) -> dict[str, Any]:
    merged_headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if headers:
        merged_headers.update(headers)
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            if delay:
                time.sleep(delay)
            req = urllib.request.Request(url, headers=merged_headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"Expected JSON object from {url}")
            return payload
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in {400, 401, 403, 404}:
                body = exc.read().decode("utf-8", errors="replace")[:500]
                raise RuntimeError(f"HTTP {exc.code} for {url}: {body}") from exc
        except Exception as exc:  # noqa: BLE001
            last_error = exc
        if attempt + 1 < retries:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Request failed after {retries} attempts: {url}") from last_error


def reconstruct_openalex_abstract(inverted: Any) -> str:
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


def record_base(database: str, query_id: str) -> dict[str, Any]:
    return {
        "database": database,
        "query_ids": query_id,
        "record_id": "",
        "doi": "",
        "title": "",
        "abstract": "",
        "year": "",
        "publication_date": "",
        "work_type": "",
        "language": "",
        "journal": "",
        "authors": "",
        "url": "",
        "is_retracted": False,
        "cited_by_count": "",
        "retrieved_at_utc": utc_now(),
    }


def openalex_search(
    query_id: str,
    query: str,
    *,
    api_key: str,
    max_records: int,
    timeout: int,
    retries: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    cursor = "*"
    total_count: int | None = None
    pages = 0
    while len(records) < max_records:
        params = {
            "api_key": api_key,
            "search": query,
            "filter": "is_retracted:false",
            "per_page": 100,
            "cursor": cursor,
            "select": "id,doi,title,display_name,publication_year,publication_date,type,language,cited_by_count,is_retracted,abstract_inverted_index,authorships,primary_location,referenced_works",
        }
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
        payload = request_json(url, timeout=timeout, retries=retries, delay=0.05)
        pages += 1
        meta = payload.get("meta") if isinstance(payload.get("meta"), dict) else {}
        if total_count is None:
            try:
                total_count = int(meta.get("count"))
            except (TypeError, ValueError):
                total_count = None
        results = payload.get("results")
        if not isinstance(results, list) or not results:
            break
        for item in results:
            if not isinstance(item, dict):
                continue
            row = record_base("openalex", query_id)
            authorships = item.get("authorships") if isinstance(item.get("authorships"), list) else []
            authors: list[str] = []
            for authorship in authorships:
                if not isinstance(authorship, dict):
                    continue
                author = authorship.get("author") if isinstance(authorship.get("author"), dict) else {}
                name = clean_text(author.get("display_name"))
                if name:
                    authors.append(name)
            primary = item.get("primary_location") if isinstance(item.get("primary_location"), dict) else {}
            source = primary.get("source") if isinstance(primary.get("source"), dict) else {}
            row.update(
                {
                    "record_id": clean_text(item.get("id")),
                    "doi": normalize_doi(item.get("doi")),
                    "title": clean_text(item.get("title") or item.get("display_name")),
                    "abstract": reconstruct_openalex_abstract(item.get("abstract_inverted_index")),
                    "year": item.get("publication_year") or "",
                    "publication_date": clean_text(item.get("publication_date")),
                    "work_type": clean_text(item.get("type")),
                    "language": clean_text(item.get("language")),
                    "journal": clean_text(source.get("display_name")),
                    "authors": "; ".join(authors),
                    "url": clean_text(primary.get("landing_page_url") or item.get("id")),
                    "is_retracted": bool(item.get("is_retracted")),
                    "cited_by_count": item.get("cited_by_count") or 0,
                    "referenced_works": ";".join(str(x) for x in (item.get("referenced_works") or [])),
                }
            )
            records.append(row)
            if len(records) >= max_records:
                break
        next_cursor = meta.get("next_cursor")
        if not next_cursor:
            break
        cursor = str(next_cursor)
    return records, {
        "database": "openalex",
        "query_id": query_id,
        "query": query,
        "reported_count": total_count,
        "retrieved_count": len(records),
        "pages": pages,
        "truncated": total_count is not None and len(records) < total_count,
    }


def crossref_search(
    query_id: str,
    query: str,
    *,
    mailto: str,
    max_records: int,
    timeout: int,
    retries: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    cursor = "*"
    total_count: int | None = None
    pages = 0
    while len(records) < max_records:
        params: dict[str, Any] = {
            "query.bibliographic": query,
            "rows": 1000,
            "cursor": cursor,
            "select": "DOI,title,abstract,published,created,type,language,container-title,author,URL,is-referenced-by-count",
        }
        if mailto:
            params["mailto"] = mailto
        url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
        payload = request_json(url, timeout=timeout, retries=retries, delay=0.1)
        pages += 1
        message = payload.get("message") if isinstance(payload.get("message"), dict) else {}
        if total_count is None:
            try:
                total_count = int(message.get("total-results"))
            except (TypeError, ValueError):
                total_count = None
        items = message.get("items")
        if not isinstance(items, list) or not items:
            break
        for item in items:
            if not isinstance(item, dict):
                continue
            row = record_base("crossref", query_id)
            titles = item.get("title") if isinstance(item.get("title"), list) else []
            journals = item.get("container-title") if isinstance(item.get("container-title"), list) else []
            authors: list[str] = []
            for author in item.get("author") or []:
                if isinstance(author, dict):
                    name = clean_text(f"{author.get('given', '')} {author.get('family', '')}")
                    if name:
                        authors.append(name)
            date_parts: list[Any] = []
            for date_field in ("published-print", "published-online", "published", "created"):
                value = item.get(date_field)
                if isinstance(value, dict):
                    parts = value.get("date-parts")
                    if isinstance(parts, list) and parts and isinstance(parts[0], list):
                        date_parts = parts[0]
                        break
            year = date_parts[0] if date_parts else ""
            publication_date = "-".join(str(x).zfill(2) for x in date_parts) if date_parts else ""
            row.update(
                {
                    "record_id": normalize_doi(item.get("DOI")) or clean_text(item.get("URL")),
                    "doi": normalize_doi(item.get("DOI")),
                    "title": clean_text(titles[0] if titles else ""),
                    "abstract": clean_text(item.get("abstract")),
                    "year": year,
                    "publication_date": publication_date,
                    "work_type": clean_text(item.get("type")),
                    "language": clean_text(item.get("language")),
                    "journal": clean_text(journals[0] if journals else ""),
                    "authors": "; ".join(authors),
                    "url": clean_text(item.get("URL")),
                    "cited_by_count": item.get("is-referenced-by-count") or 0,
                }
            )
            records.append(row)
            if len(records) >= max_records:
                break
        if len(items) < 1000:
            break
        next_cursor = message.get("next-cursor")
        if not next_cursor or str(next_cursor) == cursor:
            break
        cursor = str(next_cursor)
    return records, {
        "database": "crossref",
        "query_id": query_id,
        "query": query,
        "reported_count": total_count,
        "retrieved_count": len(records),
        "pages": pages,
        "truncated": total_count is not None and len(records) < total_count,
    }


def europe_pmc_search(
    query_id: str,
    query: str,
    *,
    email: str,
    max_records: int,
    timeout: int,
    retries: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    cursor = "*"
    total_count: int | None = None
    pages = 0
    while len(records) < max_records:
        params: dict[str, Any] = {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": 1000,
            "cursorMark": cursor,
        }
        if email:
            params["email"] = email
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(params)
        payload = request_json(url, timeout=timeout, retries=retries, delay=0.05)
        pages += 1
        try:
            total_count = int(payload.get("hitCount")) if total_count is None else total_count
        except (TypeError, ValueError):
            total_count = total_count
        result_list = payload.get("resultList") if isinstance(payload.get("resultList"), dict) else {}
        results = result_list.get("result")
        if not isinstance(results, list) or not results:
            break
        for item in results:
            if not isinstance(item, dict):
                continue
            row = record_base("europe_pmc", query_id)
            author_list = item.get("authorList") if isinstance(item.get("authorList"), dict) else {}
            authors = []
            for author in author_list.get("author") or []:
                if isinstance(author, dict):
                    name = clean_text(author.get("fullName") or f"{author.get('firstName', '')} {author.get('lastName', '')}")
                    if name:
                        authors.append(name)
            doi = normalize_doi(item.get("doi"))
            source = clean_text(item.get("source"))
            ext_id = clean_text(item.get("id") or item.get("extId"))
            journal_info = item.get("journalInfo") if isinstance(item.get("journalInfo"), dict) else {}
            publication_date = item.get("firstPublicationDate") or journal_info.get("printPublicationDate") or ""
            row.update(
                {
                    "record_id": f"{source}:{ext_id}" if source or ext_id else doi,
                    "doi": doi,
                    "title": clean_text(item.get("title")),
                    "abstract": clean_text(item.get("abstractText")),
                    "year": item.get("pubYear") or "",
                    "publication_date": clean_text(publication_date),
                    "work_type": clean_text(item.get("pubType")),
                    "language": clean_text(item.get("language")),
                    "journal": clean_text(item.get("journalTitle")),
                    "authors": "; ".join(authors),
                    "url": f"https://europepmc.org/article/{source}/{ext_id}" if source and ext_id else "",
                    "cited_by_count": item.get("citedByCount") or 0,
                }
            )
            records.append(row)
            if len(records) >= max_records:
                break
        next_cursor = payload.get("nextCursorMark")
        if not next_cursor or str(next_cursor) == cursor:
            break
        cursor = str(next_cursor)
    return records, {
        "database": "europe_pmc",
        "query_id": query_id,
        "query": query,
        "reported_count": total_count,
        "retrieved_count": len(records),
        "pages": pages,
        "truncated": total_count is not None and len(records) < total_count,
    }


def openalex_get_by_doi(doi: str, *, api_key: str, timeout: int, retries: int) -> dict[str, Any] | None:
    params = {
        "api_key": api_key,
        "filter": f"doi:{normalize_doi(doi)}",
        "per_page": 1,
    }
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
    try:
        payload = request_json(url, timeout=timeout, retries=retries)
    except RuntimeError:
        return None
    results = payload.get("results")
    return results[0] if isinstance(results, list) and results and isinstance(results[0], dict) else None


def openalex_get_many_by_ids(ids: list[str], *, api_key: str, timeout: int, retries: int) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for start in range(0, len(ids), 50):
        batch = ids[start:start + 50]
        if not batch:
            continue
        short_ids = [x.rsplit("/", 1)[-1] for x in batch]
        params = {
            "api_key": api_key,
            "filter": "openalex:" + "|".join(short_ids),
            "per_page": 100,
        }
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
        payload = request_json(url, timeout=timeout, retries=retries, delay=0.05)
        results = payload.get("results")
        if isinstance(results, list):
            records.extend(x for x in results if isinstance(x, dict))
    return records


def openalex_citing(seed_id: str, *, api_key: str, max_records: int, timeout: int, retries: int) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    cursor = "*"
    short_id = seed_id.rsplit("/", 1)[-1]
    while len(records) < max_records:
        params = {
            "api_key": api_key,
            "filter": f"cites:{short_id},is_retracted:false",
            "per_page": 100,
            "cursor": cursor,
        }
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
        payload = request_json(url, timeout=timeout, retries=retries, delay=0.05)
        results = payload.get("results")
        if not isinstance(results, list) or not results:
            break
        records.extend(x for x in results if isinstance(x, dict))
        meta = payload.get("meta") if isinstance(payload.get("meta"), dict) else {}
        next_cursor = meta.get("next_cursor")
        if not next_cursor:
            break
        cursor = str(next_cursor)
    return records[:max_records]


def openalex_item_to_record(item: dict[str, Any], query_id: str) -> dict[str, Any]:
    row = record_base("openalex_citation", query_id)
    authors = []
    for authorship in item.get("authorships") or []:
        if isinstance(authorship, dict) and isinstance(authorship.get("author"), dict):
            name = clean_text(authorship["author"].get("display_name"))
            if name:
                authors.append(name)
    primary = item.get("primary_location") if isinstance(item.get("primary_location"), dict) else {}
    source = primary.get("source") if isinstance(primary.get("source"), dict) else {}
    row.update(
        {
            "record_id": clean_text(item.get("id")),
            "doi": normalize_doi(item.get("doi")),
            "title": clean_text(item.get("title") or item.get("display_name")),
            "abstract": reconstruct_openalex_abstract(item.get("abstract_inverted_index")),
            "year": item.get("publication_year") or "",
            "publication_date": clean_text(item.get("publication_date")),
            "work_type": clean_text(item.get("type")),
            "language": clean_text(item.get("language")),
            "journal": clean_text(source.get("display_name")),
            "authors": "; ".join(authors),
            "url": clean_text(primary.get("landing_page_url") or item.get("id")),
            "is_retracted": bool(item.get("is_retracted")),
            "cited_by_count": item.get("cited_by_count") or 0,
            "referenced_works": ";".join(str(x) for x in (item.get("referenced_works") or [])),
        }
    )
    return row


def citation_chase(
    seeds: list[dict[str, str]],
    *,
    api_key: str,
    max_records_per_seed: int,
    timeout: int,
    retries: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    logs: list[dict[str, Any]] = []
    for seed in seeds:
        doi = normalize_doi(seed.get("doi"))
        seed_item = openalex_get_by_doi(doi, api_key=api_key, timeout=timeout, retries=retries)
        if not seed_item:
            logs.append({"database": "openalex_citation", "query_id": f"seed:{doi}", "query": doi, "reported_count": None, "retrieved_count": 0, "pages": 0, "truncated": False, "status": "seed_not_resolved"})
            continue
        seed_id = clean_text(seed_item.get("id"))
        backward_ids = [str(x) for x in (seed_item.get("referenced_works") or [])]
        backward_items = openalex_get_many_by_ids(backward_ids, api_key=api_key, timeout=timeout, retries=retries)
        forward_items = openalex_citing(seed_id, api_key=api_key, max_records=max_records_per_seed, timeout=timeout, retries=retries)
        records.append(openalex_item_to_record(seed_item, f"seed:{doi}"))
        records.extend(openalex_item_to_record(x, f"backward:{doi}") for x in backward_items)
        records.extend(openalex_item_to_record(x, f"forward:{doi}") for x in forward_items)
        logs.append(
            {
                "database": "openalex_citation",
                "query_id": f"seed:{doi}",
                "query": seed.get("reason", ""),
                "reported_count": len(backward_ids) + len(forward_items) + 1,
                "retrieved_count": len(backward_items) + len(forward_items) + 1,
                "pages": None,
                "truncated": len(forward_items) >= max_records_per_seed,
                "status": "complete",
            }
        )
    return records, logs


def load_census(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    names: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = clean_text(row.get("canonical_name"))
            if len(name.split()) == 2:
                names.add(name)
    return names


def candidate_species(text: str, census: set[str]) -> list[str]:
    found = {" ".join(match) for match in BINOMIAL_RE.findall(text)}
    if census:
        found = {name for name in found if name in census}
    return sorted(found)


def any_term(text: str, terms: Iterable[str]) -> tuple[bool, list[str]]:
    lowered = text.lower()
    matched = sorted({term for term in terms if term.lower() in lowered})
    return bool(matched), matched


def screen_record(record: dict[str, Any], terms: dict[str, list[str]], census: set[str]) -> dict[str, Any]:
    text = clean_text(f"{record.get('title', '')} {record.get('abstract', '')}")
    signals: dict[str, bool] = {}
    matches: dict[str, list[str]] = {}
    for label, vocab in terms.items():
        signals[label], matches[label] = any_term(text, vocab)

    species = candidate_species(text, census)
    likely_hard_exclusion = (
        (signals.get("cultivated") or signals.get("induced") or signals.get("ontogenetic"))
        and not signals.get("natural")
        and not signals.get("within")
        and not signals.get("among")
    )
    core_relevance = signals.get("display_colour", False) and signals.get("variation", False)

    if core_relevance and signals.get("natural") and not likely_hard_exclusion:
        priority = "P1_high_natural_itv"
    elif core_relevance and (signals.get("within") or signals.get("among")) and not likely_hard_exclusion:
        priority = "P1_high_population_itv"
    elif core_relevance and not likely_hard_exclusion:
        priority = "P2_possible_itv"
    elif core_relevance and likely_hard_exclusion:
        priority = "P3_likely_exclusion_review"
    elif signals.get("display_colour"):
        priority = "P4_flower_colour_context_only"
    else:
        priority = "P5_low_relevance"

    if signals.get("within") and signals.get("among"):
        provisional_scale = "mixed_signal"
    elif signals.get("within"):
        provisional_scale = "within_signal"
    elif signals.get("among"):
        provisional_scale = "among_signal"
    else:
        provisional_scale = "unresolved"

    if signals.get("discrete"):
        provisional_form = "discrete_signal"
    elif signals.get("continuous"):
        provisional_form = "continuous_signal"
    else:
        provisional_form = "unspecified"

    return {
        **record,
        "dedup_key": stable_key(record),
        "candidate_species": ";".join(species),
        "candidate_species_count": len(species),
        "screen_priority": priority,
        "provisional_scale_signal": provisional_scale,
        "provisional_variation_form": provisional_form,
        "display_colour_signal": int(signals.get("display_colour", False)),
        "variation_signal": int(signals.get("variation", False)),
        "within_signal": int(signals.get("within", False)),
        "among_signal": int(signals.get("among", False)),
        "natural_signal": int(signals.get("natural", False)),
        "cultivated_signal": int(signals.get("cultivated", False)),
        "induced_signal": int(signals.get("induced", False)),
        "ontogenetic_signal": int(signals.get("ontogenetic", False)),
        "non_display_floral_signal": int(signals.get("non_display_floral", False)),
        "matched_terms": json.dumps(matches, ensure_ascii=False, sort_keys=True),
        "review_status": "unreviewed",
        "reviewer_1_decision": "",
        "reviewer_1_reason": "",
        "reviewer_2_decision": "",
        "reviewer_2_reason": "",
        "adjudicated_decision": "",
        "natural_status": "",
        "variation_form": "",
        "spatial_scale": "",
        "display_organ": "",
        "evidence_quote": "",
        "evidence_page_or_section": "",
    }


def deduplicate(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        title = clean_text(record.get("title"))
        if not title:
            continue
        grouped[stable_key(record)].append(record)

    output: list[dict[str, Any]] = []
    for key, items in grouped.items():
        items.sort(
            key=lambda x: (
                bool(clean_text(x.get("abstract"))),
                len(clean_text(x.get("abstract"))),
                int(x.get("cited_by_count") or 0),
            ),
            reverse=True,
        )
        best = dict(items[0])
        best["database"] = ";".join(sorted({clean_text(x.get("database")) for x in items if clean_text(x.get("database"))}))
        best["query_ids"] = ";".join(sorted({q for x in items for q in clean_text(x.get("query_ids")).split(";") if q}))
        best["record_id"] = ";".join(sorted({clean_text(x.get("record_id")) for x in items if clean_text(x.get("record_id"))}))
        best["duplicate_record_count"] = len(items)
        best["dedup_key"] = key
        output.append(best)
    output.sort(key=lambda x: (-int(x.get("cited_by_count") or 0), str(x.get("year") or ""), normalize_title(x.get("title"))))
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--openalex-api-key-env", default="OPENALEX_API_KEY")
    parser.add_argument("--contact-email-env", default="LITERATURE_SEARCH_CONTACT_EMAIL")
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--max-records-per-query", type=int)
    parser.add_argument("--max-citations-per-seed", type=int, default=5000)
    parser.add_argument("--databases", default="openalex,crossref,europe_pmc")
    parser.add_argument("--skip-citation-chasing", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    openalex_key = os.environ.get(args.openalex_api_key_env, "").strip()
    contact_email = os.environ.get(args.contact_email_env, "").strip()
    databases = {x.strip() for x in args.databases.split(",") if x.strip()}
    max_records = args.max_records_per_query or int(config.get("max_records_per_query", 10000))
    if "openalex" in databases and not openalex_key:
        raise RuntimeError(f"{args.openalex_api_key_env} is required for OpenAlex")

    started = utc_now()
    raw_records: list[dict[str, Any]] = []
    search_logs: list[dict[str, Any]] = []
    for query in config.get("queries", []):
        query_id = str(query["id"])
        if "openalex" in databases:
            records, log = openalex_search(
                query_id,
                str(query["openalex"]),
                api_key=openalex_key,
                max_records=max_records,
                timeout=args.timeout,
                retries=args.retries,
            )
            raw_records.extend(records)
            search_logs.append(log)
            print(json.dumps(log, ensure_ascii=False), flush=True)
        if "crossref" in databases:
            records, log = crossref_search(
                query_id,
                str(query["crossref"]),
                mailto=contact_email,
                max_records=max_records,
                timeout=args.timeout,
                retries=args.retries,
            )
            raw_records.extend(records)
            search_logs.append(log)
            print(json.dumps(log, ensure_ascii=False), flush=True)
        if "europe_pmc" in databases:
            records, log = europe_pmc_search(
                query_id,
                str(query["europe_pmc"]),
                email=contact_email,
                max_records=max_records,
                timeout=args.timeout,
                retries=args.retries,
            )
            raw_records.extend(records)
            search_logs.append(log)
            print(json.dumps(log, ensure_ascii=False), flush=True)

    if not args.skip_citation_chasing:
        citation_records, citation_logs = citation_chase(
            list(config.get("seed_reviews", [])),
            api_key=openalex_key,
            max_records_per_seed=args.max_citations_per_seed,
            timeout=args.timeout,
            retries=args.retries,
        )
        raw_records.extend(citation_records)
        search_logs.extend(citation_logs)

    raw_path = outdir / "itv_fcp_search_raw_records.csv"
    log_path = outdir / "itv_fcp_search_log.csv"
    write_csv(raw_path, raw_records)
    write_csv(log_path, search_logs)

    deduped = deduplicate(raw_records)
    census = load_census(Path(args.census))
    terms = config.get("screening_terms")
    if not isinstance(terms, dict):
        raise ValueError("screening_terms missing from config")
    screened = [screen_record(row, terms, census) for row in deduped]
    priority_order = {
        "P1_high_natural_itv": 0,
        "P1_high_population_itv": 1,
        "P2_possible_itv": 2,
        "P3_likely_exclusion_review": 3,
        "P4_flower_colour_context_only": 4,
        "P5_low_relevance": 5,
    }
    screened.sort(
        key=lambda x: (
            priority_order.get(str(x.get("screen_priority")), 99),
            -int(x.get("candidate_species_count") or 0),
            -int(x.get("cited_by_count") or 0),
            normalize_title(x.get("title")),
        )
    )

    dedup_path = outdir / "itv_fcp_search_deduplicated_records.csv"
    queue_path = outdir / "itv_fcp_human_screening_queue.csv"
    write_csv(dedup_path, deduped)
    write_csv(queue_path, screened)

    priority_counts: dict[str, int] = defaultdict(int)
    for row in screened:
        priority_counts[str(row.get("screen_priority"))] += 1
    truncated = [log for log in search_logs if log.get("truncated")]
    manifest = {
        "status": "complete",
        "protocol_version": config.get("protocol_version"),
        "started_at_utc": started,
        "completed_at_utc": utc_now(),
        "databases": sorted(databases),
        "query_count": len(config.get("queries", [])),
        "seed_review_count": len(config.get("seed_reviews", [])),
        "raw_records": len(raw_records),
        "deduplicated_records": len(deduped),
        "duplicate_records_removed": len(raw_records) - len(deduped),
        "screen_priority_counts": dict(sorted(priority_counts.items())),
        "records_with_candidate_binomials": sum(int(x.get("candidate_species_count") or 0) > 0 for x in screened),
        "queries_truncated": len(truncated),
        "truncated_query_ids": [f"{x.get('database')}:{x.get('query_id')}" for x in truncated],
        "human_review_required": True,
        "semantic_guard": (
            "Automated retrieval and screening identify candidate publications only. "
            "Natural ITV, strict FCP, spatial scale and species eligibility require source-level human review."
        ),
        "files": {
            raw_path.name: {"sha256": sha256(raw_path), "rows": len(raw_records)},
            log_path.name: {"sha256": sha256(log_path), "rows": len(search_logs)},
            dedup_path.name: {"sha256": sha256(dedup_path), "rows": len(deduped)},
            queue_path.name: {"sha256": sha256(queue_path), "rows": len(screened)},
            config_path.name: {"sha256": sha256(config_path)},
        },
    }
    (outdir / "itv_fcp_search_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
