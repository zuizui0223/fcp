#!/usr/bin/env python3
"""Targeted OpenAlex follow-up for strong deferred flower-colour candidates.

The follow-up is deliberately best-effort: external API failures are recorded in
QC and do not abort the full discovery workflow. Candidate evidence remains
unverified until primary review.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List

OPENALEX = "https://api.openalex.org/works"
TAG_RE = re.compile(r"<[^>]+>")
RELEVANT_RE = re.compile(
    r"\b(?:flower|floral|petal|corolla)\b.{0,120}\b(?:colou?r|morph|polymorph|dimorph|variation)\b|"
    r"\b(?:colou?r|morph|polymorph|dimorph|variation)\b.{0,120}\b(?:flower|floral|petal|corolla)\b",
    re.I,
)
NATURAL_RE = re.compile(
    r"\b(?:natural population|wild population|polymorphic population|within[- ]population|"
    r"morph frequenc|geographic variation|spatial variation|balancing selection|"
    r"frequency[- ]dependent selection|pollinator[- ]mediated selection)\b",
    re.I,
)
NEGATIVE_RE = re.compile(
    r"\b(?:cultivar|breeding|ornamental|qtl|transgenic|gene editing|mutagenesis|"
    r"crop improvement|market value|horticultur|floricultur)\b",
    re.I,
)
WORK_FIELDS = [
    "canonical_name", "family", "openalex_id", "title", "year", "doi",
    "landing_url", "score", "evidence_basis", "evidence_snippet",
]


def clean(value: object) -> str:
    return re.sub(r"\s+", " ", TAG_RE.sub(" ", html.unescape(str(value or "")))).strip()


def abstract_text(inverted: object) -> str:
    if not isinstance(inverted, dict):
        return ""
    positions = [
        (idx, str(word))
        for word, indexes in inverted.items()
        if isinstance(indexes, list)
        for idx in indexes
        if isinstance(idx, int)
    ]
    positions.sort()
    return clean(" ".join(word for _, word in positions))


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def request_json(url: str, timeout: int, retries: int) -> Dict[str, object]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "fcp-targeted-followup/3.0 (GitHub: zuizui0223/fcp)",
                    "Accept": "application/json",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError("OpenAlex response was not an object")
            return payload
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(str(last) if last else "OpenAlex request failed")


def choose_targets(rows: List[Dict[str, str]], limit: int) -> List[Dict[str, str]]:
    eligible = []
    for row in rows:
        cls = row.get("evidence_class", "")
        title_matches = int(row.get("n_title_matches") or 0)
        context_matches = int(row.get("n_context_matches") or 0)
        natural_signals = int(row.get("natural_signal_count") or 0)
        if cls == "possible_polymorphism":
            tier = 0
        elif cls == "insufficient_direct_evidence" and (title_matches > 0 or natural_signals > 0):
            tier = 1
        else:
            continue
        eligible.append((
            tier, -title_matches, -natural_signals, -context_matches,
            -int(row.get("max_score") or 0), row.get("canonical_name", ""), row,
        ))
    eligible.sort(key=lambda item: item[:-1])
    return [item[-1] for item in eligible[:limit]]


def score_work(name: str, title: str, abstract: str) -> tuple[int, str, str] | None:
    text = f"{title} {abstract}"
    if name.lower() not in text.lower() or not RELEVANT_RE.search(text):
        return None
    score = 10
    basis = ["exact_species"]
    snippet = title
    if name.lower() in title.lower():
        score += 10
        basis.append("title_species")
    if RELEVANT_RE.search(title):
        score += 8
        basis.append("title_colour_context")
    if NATURAL_RE.search(text):
        score += 10
        basis.append("natural_population_signal")
    if "polymorph" in text.lower() or "dimorph" in text.lower():
        score += 6
        basis.append("polymorphism_term")
    if NEGATIVE_RE.search(text):
        score -= 10
        basis.append("artificial_penalty")
    match = re.search(re.escape(name), abstract, re.I)
    if match:
        snippet = abstract[max(0, match.start() - 240):min(len(abstract), match.end() + 240)]
    if score < 14:
        return None
    return score, "+".join(basis), clean(snippet)[:500]


def grouped_targets(targets: List[Dict[str, str]], size: int) -> List[List[Dict[str, str]]]:
    return [targets[i:i + size] for i in range(0, len(targets), size)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deferred", default="data/global_flower_colour_deferred.csv")
    parser.add_argument("--species", default="data/global_flower_colour_species_ranked.csv")
    parser.add_argument("--works-out", default="data/targeted_followup_works.csv")
    parser.add_argument("--qc-out", default="data/targeted_followup_qc.json")
    parser.add_argument("--max-targets", type=int, default=80)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--per-page", type=int, default=100)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    args = parser.parse_args()

    deferred = read_rows(Path(args.deferred))
    species_rows = read_rows(Path(args.species))
    targets = choose_targets(deferred, args.max_targets)
    api_key = os.environ.get(args.api_key_env, "").strip()
    by_name = {row["canonical_name"]: row for row in species_rows}
    target_by_name = {row["canonical_name"]: row for row in targets}

    works: List[Dict[str, object]] = []
    targets_with_evidence: set[str] = set()
    seen_pairs: set[tuple[str, str]] = set()
    errors: List[Dict[str, object]] = []
    request_count = 0
    batches = grouped_targets(targets, max(1, args.batch_size))

    for batch_index, batch in enumerate(batches, 1):
        genera = sorted({row["canonical_name"].split()[0] for row in batch})
        query = " OR ".join(f'"{genus}"' for genus in genera)
        query = f'({query}) "flower color"'
        params: Dict[str, object] = {
            "search": query,
            "per-page": min(max(args.per_page, 1), 200),
        }
        if api_key:
            params["api_key"] = api_key
        request_count += 1
        try:
            payload = request_json(
                OPENALEX + "?" + urllib.parse.urlencode(params),
                args.timeout,
                args.retries,
            )
        except RuntimeError as exc:
            errors.append({"batch": batch_index, "genera": genera, "error": str(exc)[:500]})
            continue

        results = payload.get("results", [])
        if not isinstance(results, list):
            errors.append({"batch": batch_index, "genera": genera, "error": "results_not_list"})
            continue

        best_by_name: Dict[str, Dict[str, object]] = {}
        names = [row["canonical_name"] for row in batch]
        for item in results:
            if not isinstance(item, dict):
                continue
            openalex_id = str(item.get("id") or "")
            if not openalex_id:
                continue
            title = clean(item.get("title"))
            abstract = abstract_text(item.get("abstract_inverted_index"))
            text_lower = f"{title} {abstract}".lower()
            for name in names:
                if name.lower() not in text_lower:
                    continue
                evidence = score_work(name, title, abstract)
                if evidence is None or (name, openalex_id) in seen_pairs:
                    continue
                seen_pairs.add((name, openalex_id))
                score, basis, snippet = evidence
                primary = item.get("primary_location") or {}
                landing = str(primary.get("landing_page_url") or openalex_id) if isinstance(primary, dict) else openalex_id
                row = {
                    "canonical_name": name,
                    "family": target_by_name[name].get("family", ""),
                    "openalex_id": openalex_id,
                    "title": title,
                    "year": item.get("publication_year") or "",
                    "doi": str(item.get("doi") or "").replace("https://doi.org/", ""),
                    "landing_url": landing,
                    "score": score,
                    "evidence_basis": basis,
                    "evidence_snippet": snippet,
                }
                works.append(row)
                targets_with_evidence.add(name)
                if name not in best_by_name or score > int(best_by_name[name]["score"]):
                    best_by_name[name] = row

        for name, best in best_by_name.items():
            row = by_name[name]
            row["n_works"] = str(int(row.get("n_works") or 0) + 1)
            row["n_title_matches"] = str(int(row.get("n_title_matches") or 0) + int("title_species" in str(best["evidence_basis"])))
            row["n_context_matches"] = str(int(row.get("n_context_matches") or 0) + 1)
            row["total_score"] = str(int(row.get("total_score") or 0) + int(best["score"]))
            if int(best["score"]) > int(row.get("max_score") or 0):
                row["max_score"] = str(best["score"])
                row["best_title"] = str(best["title"])
                row["best_doi"] = str(best["doi"])
                row["best_openalex_id"] = str(best["openalex_id"])
                row["best_match_evidence"] = str(best["evidence_snippet"])

    species_rows.sort(key=lambda row: (
        -int(row.get("n_title_matches") or 0), -int(row.get("max_score") or 0),
        -int(row.get("total_score") or 0), -int(row.get("n_works") or 0),
        row.get("canonical_name", ""),
    ))
    for rank, row in enumerate(species_rows, 1):
        row["rank"] = str(rank)

    with Path(args.species).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(species_rows[0].keys()))
        writer.writeheader()
        writer.writerows(species_rows)

    with Path(args.works_out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=WORK_FIELDS)
        writer.writeheader()
        writer.writerows(sorted(works, key=lambda row: (-int(row["score"]), str(row["canonical_name"]))))

    qc = {
        "eligible_targets": len(targets),
        "target_batches": len(batches),
        "requests_used": request_count,
        "failed_batches": len(errors),
        "errors": errors,
        "targets_with_new_evidence": len(targets_with_evidence),
        "retained_followup_works": len(works),
        "api_key_present": bool(api_key),
        "mode": "targeted_deferred_followup_resilient_v3",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc), flush=True)


if __name__ == "__main__":
    main()
