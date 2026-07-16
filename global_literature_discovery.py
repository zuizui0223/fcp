#!/usr/bin/env python3
"""Discover flower-colour polymorphism literature and map only context-supported species.

Broad OpenAlex queries retrieve a small literature corpus. Candidate binomials are
matched to the GBIF census, but a species is retained only when its mention is in the
title or near flower-colour / morph / polymorphism language in the abstract. This
prevents long species lists in reviews, floras, pollinator surveys, and genomics papers
from being promoted wholesale.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

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
    "title_matched_species", "context_matched_species", "raw_species_mentions",
    "match_evidence",
]
SPECIES_FIELDS = [
    "rank", "canonical_name", "family", "n_works", "max_score", "total_score",
    "n_title_matches", "n_context_matches", "best_title", "best_doi",
    "best_openalex_id", "best_match_evidence", "review_status",
]
TOKEN_RE = re.compile(r"\b[A-Z][a-z-]+\s+[a-z][a-z-]+\b")
TAG_RE = re.compile(r"<[^>]+>")
POSITIVE_RE = re.compile(
    r"\b(?:flower|floral|petal|corolla|colour|color|pigment|anthocyanin|"
    r"polymorph(?:ism|ic)?|morphs?|variation|variant|white[- ]flowered|"
    r"pink[- ]flowered|blue[- ]flowered|yellow[- ]flowered)\b",
    re.IGNORECASE,
)
STRONG_RE = re.compile(
    r"\b(?:flower|floral|petal|corolla)\b.{0,80}\b(?:colour|color|morph|polymorph|variation|variant)\b|"
    r"\b(?:colour|color|morph|polymorph|variation|variant)\b.{0,80}\b(?:flower|floral|petal|corolla)\b",
    re.IGNORECASE,
)
NEGATIVE_TITLE_RE = re.compile(
    r"\b(?:genome|genomic|transcriptome|phylogenom|RAPD|checklist|flora of|"
    r"pollinator enhancement|yield|invasive species|genetic diversity|"
    r"metabarcoding|community assembl|species richness)\b",
    re.IGNORECASE,
)


def clean_text(value: object) -> str:
    return re.sub(r"\s+", " ", TAG_RE.sub(" ", html.unescape(str(value or "")))).strip()


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
    return clean_text(" ".join(word for _, word in positions))


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
                    "User-Agent": "fcp-global-literature-discovery/2.0 (GitHub: zuizui0223/fcp)",
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


def census_species_mentions(text: str, name_to_family: Dict[str, str], genus_to_names: Dict[str, set[str]]) -> List[str]:
    matches: set[str] = set()
    for candidate in TOKEN_RE.findall(text):
        genus = candidate.split()[0]
        if candidate in genus_to_names.get(genus, ()) and candidate in name_to_family:
            matches.add(candidate)
    return sorted(matches)


def evidence_for_species(name: str, title: str, abstract: str, raw_count: int) -> Tuple[int, str, str] | None:
    """Return score, basis, and evidence snippet when the mention is context-supported."""
    title_lower = title.lower()
    name_lower = name.lower()
    title_has_name = name_lower in title_lower
    title_relevant = bool(STRONG_RE.search(title) or (POSITIVE_RE.search(title) and "flower" in title_lower))
    negative_title = bool(NEGATIVE_TITLE_RE.search(title))

    best_context = ""
    context_strength = 0
    for match in re.finditer(re.escape(name), abstract, flags=re.IGNORECASE):
        lo = max(0, match.start() - 220)
        hi = min(len(abstract), match.end() + 220)
        context = abstract[lo:hi]
        if STRONG_RE.search(context):
            strength = 2
        elif POSITIVE_RE.search(context):
            strength = 1
        else:
            strength = 0
        if strength > context_strength:
            context_strength = strength
            best_context = context

    if raw_count > 5 and not title_has_name and context_strength < 2:
        return None
    if not title_has_name and context_strength == 0:
        return None
    if negative_title and not title_has_name and context_strength < 2:
        return None

    score = 0
    basis: List[str] = []
    evidence = ""
    if title_has_name:
        score += 12
        basis.append("title_name")
        evidence = title
    if title_relevant:
        score += 8
        basis.append("title_flower_colour_context")
    if context_strength == 2:
        score += 10
        basis.append("strong_abstract_context")
        evidence = best_context or evidence
    elif context_strength == 1:
        score += 5
        basis.append("weak_abstract_context")
        evidence = best_context or evidence
    if "polymorph" in f"{title} {best_context}".lower():
        score += 6
        basis.append("polymorphism_term")
    if "morph" in f"{title} {best_context}".lower():
        score += 3
        basis.append("morph_term")
    if negative_title:
        score -= 8
        basis.append("negative_title_penalty")
    if raw_count > 5:
        score -= min(10, raw_count - 5)
        basis.append("species_list_penalty")
    if score < 8:
        return None
    return score, "+".join(basis), clean_text(evidence)[:500]


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
    rejected_mentions = 0

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
                title = clean_text(item.get("title"))
                abstract = reconstruct_abstract(item.get("abstract_inverted_index"))
                raw_mentions = census_species_mentions(f"{title} {abstract}", name_to_family, genus_to_names)
                retained: List[str] = []
                title_matches: List[str] = []
                context_matches: List[str] = []
                evidence_parts: List[str] = []
                species_scores: List[int] = []
                for name in raw_mentions:
                    evidence = evidence_for_species(name, title, abstract, len(raw_mentions))
                    if evidence is None:
                        rejected_mentions += 1
                        continue
                    score, basis, snippet = evidence
                    retained.append(name)
                    species_scores.append(score)
                    if "title_name" in basis:
                        title_matches.append(name)
                    if "abstract_context" in basis:
                        context_matches.append(name)
                    evidence_parts.append(f"{name}|{score}|{basis}|{snippet}")
                if not retained:
                    continue
                primary = item.get("primary_location") or {}
                landing = openalex_id
                if isinstance(primary, dict):
                    landing = str(primary.get("landing_page_url") or landing)
                doi = str(item.get("doi") or "").replace("https://doi.org/", "")
                works.append({
                    "openalex_id": openalex_id,
                    "query": query,
                    "title": title,
                    "year": item.get("publication_year") or "",
                    "doi": doi,
                    "landing_url": landing,
                    "matched_species": ";".join(retained),
                    "matched_families": ";".join(sorted({name_to_family[n] for n in retained if name_to_family[n]})),
                    "candidate_score": max(species_scores),
                    "title_matched_species": ";".join(title_matches),
                    "context_matched_species": ";".join(context_matches),
                    "raw_species_mentions": len(raw_mentions),
                    "match_evidence": " || ".join(evidence_parts),
                })
            meta = payload.get("meta", {})
            next_cursor = meta.get("next_cursor") if isinstance(meta, dict) else None
            print(json.dumps({
                "query": query,
                "page": page + 1,
                "returned": len(results),
                "retained_works_so_far": len(works),
                "rejected_mentions_so_far": rejected_mentions,
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
        title_set = set(str(work["title_matched_species"]).split(";"))
        context_set = set(str(work["context_matched_species"]).split(";"))
        evidence_by_name: Dict[str, str] = {}
        score_by_name: Dict[str, int] = {}
        for part in str(work["match_evidence"]).split(" || "):
            bits = part.split("|", 3)
            if len(bits) == 4:
                evidence_by_name[bits[0]] = bits[3]
                try:
                    score_by_name[bits[0]] = int(bits[1])
                except ValueError:
                    score_by_name[bits[0]] = int(work["candidate_score"])
        for name in str(work["matched_species"]).split(";"):
            if not name:
                continue
            score = score_by_name.get(name, int(work["candidate_score"]))
            entry = by_species.setdefault(name, {
                "canonical_name": name,
                "family": name_to_family.get(name, ""),
                "n_works": 0,
                "max_score": 0,
                "total_score": 0,
                "n_title_matches": 0,
                "n_context_matches": 0,
                "best_title": "",
                "best_doi": "",
                "best_openalex_id": "",
                "best_match_evidence": "",
                "review_status": "unreviewed",
            })
            entry["n_works"] = int(entry["n_works"]) + 1
            entry["total_score"] = int(entry["total_score"]) + score
            entry["n_title_matches"] = int(entry["n_title_matches"]) + int(name in title_set)
            entry["n_context_matches"] = int(entry["n_context_matches"]) + int(name in context_set)
            if score > int(entry["max_score"]):
                entry["max_score"] = score
                entry["best_title"] = work["title"]
                entry["best_doi"] = work["doi"]
                entry["best_openalex_id"] = work["openalex_id"]
                entry["best_match_evidence"] = evidence_by_name.get(name, "")

    ranked = sorted(
        by_species.values(),
        key=lambda row: (
            -int(row["n_title_matches"]), -int(row["max_score"]),
            -int(row["total_score"]), -int(row["n_works"]), str(row["canonical_name"]),
        ),
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
        "retained_works": len(works),
        "candidate_species": len(ranked),
        "rejected_species_mentions": rejected_mentions,
        "works_out": str(works_path),
        "species_out": str(species_path),
        "mode": "context_supported_species_mapping_v2",
    }
    print(json.dumps(summary), flush=True)
    if not works or not ranked:
        raise RuntimeError("Context filtering produced no mapped works/species; inspect thresholds")


if __name__ == "__main__":
    main()
