#!/usr/bin/env python3
"""Mass-screen angiosperm species for flower-colour polymorphism literature.

This is a discovery pipeline, not an outcome classifier. It searches OpenAlex and
Crossref for every species in a census and writes a ranked candidate table.
No search hit is promoted to polymorphic without manual/primary-source review.

Example:
    python mass_literature_screen.py \
      --census data/angiosperm_census.csv \
      --out data/mass_screen_hits.csv \
      --start 0 --limit 5000

Network access is required.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

OPENALEX = "https://api.openalex.org/works"
CROSSREF = "https://api.crossref.org/works"

QUERY_TEMPLATES = (
    '"{name}" "flower color polymorphism"',
    '"{name}" "flower colour polymorphism"',
    '"{name}" "floral color polymorphism"',
    '"{name}" "floral colour polymorphism"',
    '"{name}" "color morph" flower',
    '"{name}" "colour morph" flower',
    '"{name}" "petal color variant"',
    '"{name}" "petal colour variant"',
    '"{name}" "geographic flower color variation"',
    '"{name}" "geographic flower colour variation"',
)

HIGH_TERMS = (
    "flower color polymorphism", "flower colour polymorphism",
    "floral color polymorphism", "floral colour polymorphism",
    "color morph", "colour morph", "petal color variant", "petal colour variant",
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


def _get_json(url: str, retries: int = 4, timeout: int = 45) -> Dict[str, object]:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "fcp-mass-screen/1.0 (research synthesis)"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.load(response)
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(min(2 ** attempt, 10))
    raise RuntimeError(f"Failed after {retries} attempts: {url}") from last


def search_openalex(query: str, per_page: int = 5) -> List[Hit]:
    params = urllib.parse.urlencode({"search": query, "per-page": per_page})
    data = _get_json(f"{OPENALEX}?{params}")
    hits: List[Hit] = []
    for item in data.get("results", []):
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


def search_crossref(query: str, rows: int = 5) -> List[Hit]:
    params = urllib.parse.urlencode({"query.bibliographic": query, "rows": rows})
    data = _get_json(f"{CROSSREF}?{params}")
    message = data.get("message", {})
    items = message.get("items", []) if isinstance(message, dict) else []
    hits: List[Hit] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        titles = item.get("title") or []
        title = str(titles[0]) if isinstance(titles, list) and titles else ""
        year = ""
        issued = item.get("issued") or {}
        if isinstance(issued, dict):
            parts = issued.get("date-parts") or []
            if parts and isinstance(parts[0], list) and parts[0]:
                year = str(parts[0][0])
        doi = str(item.get("DOI") or "")
        url = str(item.get("URL") or "")
        hits.append(Hit(title, year, doi, url, "crossref", query))
    return hits


def score_hit(species: str, title: str, query: str) -> int:
    text = f"{title} {query}".lower()
    score = 0
    if species.lower() in text:
        score += 4
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--out", default="data/mass_screen_hits.csv")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--queries-per-species", type=int, default=len(QUERY_TEMPLATES))
    parser.add_argument("--sleep", type=float, default=0.15)
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if out_path.exists() else "w"
    write_header = mode == "w"

    processed = 0
    emitted = 0
    with out_path.open(mode, newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()

        for row in iter_census(Path(args.census), args.start, args.limit):
            name = (row.get("canonical_name") or row.get("scientific_name") or "").strip()
            if not name:
                continue
            family = (row.get("family") or "").strip()
            all_hits: List[Hit] = []
            for template in QUERY_TEMPLATES[: max(1, args.queries_per_species)]:
                query = template.format(name=name)
                try:
                    all_hits.extend(search_openalex(query))
                except Exception:
                    pass
                try:
                    all_hits.extend(search_crossref(query))
                except Exception:
                    pass
                time.sleep(max(0.0, args.sleep))

            for hit in dedupe_hits(all_hits):
                score = score_hit(name, hit.title, hit.query)
                if score <= 0:
                    continue
                writer.writerow({
                    "canonical_name": name,
                    "family": family,
                    "query": hit.query,
                    "source": hit.source,
                    "title": hit.title,
                    "year": hit.year,
                    "doi": hit.doi,
                    "landing_url": hit.landing_url,
                    "candidate_score": score,
                    "evidence_status": "machine_candidate",
                    "notes": "Requires manual primary-source review before outcome coding.",
                })
                emitted += 1
            processed += 1
            if processed % 100 == 0:
                print(json.dumps({"processed_species": processed, "emitted_hits": emitted}))

    print(json.dumps({"processed_species": processed, "emitted_hits": emitted, "out": str(out_path)}))


if __name__ == "__main__":
    main()
