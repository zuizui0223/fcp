#!/usr/bin/env python3
"""Discover flower-colour polymorphism literature first, then map works to the census.

This avoids one API request per angiosperm species. A small set of broad OpenAlex
queries retrieves a literature corpus, and exact canonical binomials found in titles
or abstracts are reconciled against the GBIF-derived angiosperm census.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

OPENALEX = "https://api.openalex.org/works"
QUERIES = (
    '"flower color polymorphism"',
    '"flower colour polymorphism"',
    '"floral color polymorphism"',
    '"floral colour polymorphism"',
    '"flower color variation" pollinator',
    '"flower colour variation" pollinator',
    '"floral color morph"',
    '"floral colour morph"',
)
WORK_FIELDS = [
    "openalex_id", "query", "title", "year", "doi", "landing_url",
    "matched_species", "matched_families", "candidate_score",
]
SPECIES_FIELDS = [
    "rank", "canonical_name", "family", "n_works", "max_score", "total_score",
    "best_title", "best_doi", "best_openalex_id", "review_status",
]
TOKEN_RE = re.compile(r"\b[A-Z][a-z-]+\s+[a-z][a-z-]+\b")


def reconstruct_abstract(inverted: object) -> str:
    if not isinstance(inverted, dict):
        return ""
    positions: List[Tuple[int, str]] = []
    for word, indexes in inverted.items():
        if not isinstance(indexes, list):
            continue
        for index in indexes:
            if isinstance(index, int):
                positions.append((index, str(word)))
    positions.sort()
    return " ".join(word for _, word in positions)


def load_census(path: Path) -> Tuple[Dict[str, str], Dict[str, set[str]]]:
    name_to_family: Dict[str, str] = {}
    genus_to_names: Dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = (row.get("canonical_name") or "").strip()
            if len(name.split()) != 2:
                continue
            genus = name.split()[0]
            name_to_family[name] = (row.get("family") or "").strip()
            genus_to_names[genus].add(name)
    if len(name_to_family) < 10000:
        raise RuntimeError(f"Census mapping unexpectedly small: {len(name_to_family)} binomials")
    return name_to_family, genus_to_names


def request_json(url: str, *, timeout: int, retries: int) -> Dict[str, object]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "fcp-global-literature-discovery/1.0 (GitHub: zuizui0223/fcp)",
                    "Accept": "application/json",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"Unexpected OpenAlex response: {type(payload)!r}")
            return payload
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code in {401, 402, 403, 429}:
                body = exc.read().decode("utf-8", errors="replace")[:500]
                raise RuntimeError(f"OpenAlex HTTP {exc.code}: {body}") from exc
        except Exception as exc:  # noqa: BLE001
            last = exc
        if attempt + 1 < retries:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"OpenAlex request failed after {retries} attempts: {url}") from last


def find_species(text: str, name_to_family: Dict[str, str], genus_to_names: Dict[str, set[str]]) -> List[str]:
    matches: set[str] = set()
    for candidate in TOKEN_RE.findall(text):
        genus = candidate.split()[0]
        if candidate in genus_to_names.get(genus, ()) and candidate in name_to_family:
            matches.add(candidate)
    return sorted(matches)


def work_score(title: str, query: str, matched_species: List[str]) -> int:
    text = f"{title} {query}".lower()
    score = 4 * len(matched_species)
    if "polymorphism" in text:
        score += 6
    if "morph" in text:
        score += 3
    if "variation" in text:
        score += 2
    if "pollinator" in text or "pollination" in text:
        score += 2
    if "cultivar" in text or "horticultur" in text:
        score -= 4
    return score


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--works-out", default="data/global_flower_colour_works.csv")
    parser.add_argument("--species-out", default="data/global_flower_colour_species_ranked.csv")
    parser.add_argument("--pages-per-query", type=int, default=2)
    parser.add_argument("--per-page", type=int, default=200)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    args = parser.parse_args()

    name_to_family, genus_to_names = load_census(Path(args.census))
    api_key = os.environ.get(args.api_key_env, "").strip()
    seen_works: set[str] = set()
    works: List[Dict[str, object]] = []
    request_count = 0

    for query in QUERIES:
        cursor = "*"
        for page in range(args.pages_per_query):
            params: Dict[str, object] = {
                "search": query,
                "per-page": min(max(args.per_page, 1), 200),
                "cursor": cursor,
            }
            if api_key:
                params["api_key"] = api_key
            url = OPENALEX + "?" + urllib.parse.urlencode(params)
            payload = request_json(url, timeout=args.timeout, retries=args.retries)
            request_count += 1
            results = payload.get("results", [])
            if not isinstance(results, list):
                raise RuntimeError(f"OpenAlex results was not a list for query {query!r}")
            for item in results:
                if not isinstance(item, dict):
                    continue
                openalex_id = str(item.get("id") or "")
                if not openalex_id or openalex_id in seen_works:
                    continue
                seen_works.add(openalex_id)
                title = str(item.get("title") or "")
                abstract = reconstruct_abstract(item.get("abstract_inverted_index"))
                matched = find_species(f"{title} {abstract}", name_to_family, genus_to_names)
                if not matched:
                    continue
                primary = item.get("primary_location") or {}
                landing = openalex_id
                if isinstance(primary, dict):
                    landing = str(primary.get("landing_page_url") or landing)
                doi = str(item.get("doi") or "").replace("https://doi.org/", "")
                score = work_score(title, query, matched)
                works.append({
                    "openalex_id": openalex_id,
                    "query": query,
                    "title": title,
                    "year": item.get("publication_year") or "",
                    "doi": doi,
                    "landing_url": landing,
                    "matched_species": ";".join(matched),
                    "matched_families": ";".join(sorted({name_to_family[n] for n in matched if name_to_family[n]})),
                    "candidate_score": score,
                })
            meta = payload.get("meta", {})
            next_cursor = meta.get("next_cursor") if isinstance(meta, dict) else None
            print(json.dumps({
                "query": query,
                "page": page + 1,
                "returned": len(results),
                "mapped_works_so_far": len(works),
                "requests_used": request_count,
            }), flush=True)
            if not next_cursor or not results:
                break
            cursor = str(next_cursor)

    works.sort(key=lambda row: (-int(row["candidate_score"]), str(row["title"])))
    works_path = Path(args.works_out)
    works_path.parent.mkdir(parents=True, exist_ok=True)
    with works_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=WORK_FIELDS)
        writer.writeheader()
        writer.writerows(works)

    by_species: Dict[str, Dict[str, object]] = {}
    for work in works:
        for name in str(work["matched_species"]).split(";"):
            if not name:
                continue
            entry = by_species.setdefault(name, {
                "canonical_name": name,
                "family": name_to_family.get(name, ""),
                "n_works": 0,
                "max_score": 0,
                "total_score": 0,
                "best_title": "",
                "best_doi": "",
                "best_openalex_id": "",
                "review_status": "unreviewed",
            })
            score = int(work["candidate_score"])
            entry["n_works"] = int(entry["n_works"]) + 1
            entry["total_score"] = int(entry["total_score"]) + score
            if score > int(entry["max_score"]):
                entry["max_score"] = score
                entry["best_title"] = work["title"]
                entry["best_doi"] = work["doi"]
                entry["best_openalex_id"] = work["openalex_id"]

    ranked = sorted(
        by_species.values(),
        key=lambda row: (-int(row["max_score"]), -int(row["total_score"]), -int(row["n_works"]), str(row["canonical_name"])),
    )
    for rank, row in enumerate(ranked, 1):
        row["rank"] = rank
    species_path = Path(args.species_out)
    with species_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SPECIES_FIELDS)
        writer.writeheader()
        writer.writerows(ranked)

    summary = {
        "requests_used": request_count,
        "unique_works_examined": len(seen_works),
        "mapped_works": len(works),
        "candidate_species": len(ranked),
        "works_out": str(works_path),
        "species_out": str(species_path),
    }
    print(json.dumps(summary), flush=True)
    if not works or not ranked:
        raise RuntimeError("Global literature discovery produced no mapped works/species; inspect query or name mapping")


if __name__ == "__main__":
    main()
