#!/usr/bin/env python3
"""Targeted literature enrichment for ambiguous spatial-scale evidence.

This script does not discover new candidate species. It enriches the confirmed review set
with additional exact-species OpenAlex records that may distinguish within-population
polymorphism from among-population colour differentiation.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

WITHIN = re.compile(r"\b(within[- ]population|same population|coexist|co-occurr|morph frequenc|frequency[- ]dependent|polymorphic population|multiple (?:colou?r|flower) morphs|colour morphs|color morphs)\b", re.I)
GEOGRAPHIC = re.compile(r"\b(geographic|spatial variation|among populations|between populations|population differentiation|cline|hybrid zone|range edge|local adaptation|regional variation)\b", re.I)
DIRECT = re.compile(r"(?:flower|floral|petal|corolla|bract).{0,100}(?:colou?r).{0,100}(?:polymorph|morph|variation|dimorph)|(?:polymorph|dimorph|morph).{0,100}(?:flower|floral|petal|corolla|bract).{0,100}(?:colou?r)", re.I)
NEGATIVE = re.compile(r"\b(cultivar|breeding|ornamental|transgenic|mutagenesis|horticultur|floricultur)\b", re.I)
TERMS = '("flower color" OR "flower colour" OR "floral polymorphism" OR "color morph" OR "colour morph" OR "morph frequency" OR "geographic variation" OR "population differentiation" OR "local adaptation")'
FIELDS = ["canonical_name","family","openalex_id","title","year","doi","landing_url","query_mode","within_signal","among_signal","direct_colour_signal","artificial_signal","score","evidence_snippet"]


def clean(v: object) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html.unescape(str(v or "")))).strip()


def abstract(inv: object) -> str:
    if not isinstance(inv, dict):
        return ""
    pairs = [(i, str(word)) for word, idxs in inv.items() if isinstance(idxs, list) for i in idxs if isinstance(i, int)]
    return clean(" ".join(word for _, word in sorted(pairs)))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def request_json(url: str, timeout: int, retries: int, backoff: float) -> dict:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"fcp-spatial-scale-enrichment/1.0","Accept":"application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.load(response)
        except Exception as exc:
            last = exc
            if attempt + 1 < retries:
                time.sleep(backoff * (2 ** attempt))
    raise RuntimeError(str(last))


def baseline_scale(row: dict[str, str]) -> str:
    text = " ".join([row.get("best_title", ""), row.get("best_match_evidence", ""), row.get("review_reason", "")])
    w = bool(WITHIN.search(text)); g = bool(GEOGRAPHIC.search(text))
    if w and not g: return "within_population"
    if g and not w: return "among_population"
    if w and g: return "mixed"
    return "unclear"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review", required=True)
    ap.add_argument("--out", default="analysis_outputs/evidence_scale/enriched_openalex_works.csv")
    ap.add_argument("--species-summary", default="analysis_outputs/evidence_scale/enriched_species_summary.csv")
    ap.add_argument("--qc-out", default="analysis_outputs/evidence_scale/enriched_openalex_qc.json")
    ap.add_argument("--per-page", type=int, default=50)
    ap.add_argument("--timeout", type=int, default=40)
    ap.add_argument("--retries", type=int, default=4)
    ap.add_argument("--backoff", type=float, default=2.0)
    ap.add_argument("--request-delay", type=float, default=0.25)
    args = ap.parse_args()

    review = read_rows(Path(args.review))
    targets = [r for r in review if baseline_scale(r) in {"mixed", "unclear"}]
    out_rows: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for index, row in enumerate(targets):
        name = row.get("canonical_name", "").strip()
        if index:
            time.sleep(max(0.0, args.request_delay))
        params = {"search": f'"{name}" {TERMS}', "per-page": min(max(args.per_page, 1), 200)}
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
        try:
            payload = request_json(url, args.timeout, args.retries, args.backoff)
        except RuntimeError as exc:
            errors.append({"canonical_name": name, "error": str(exc)[:500]})
            continue
        for item in payload.get("results", []) if isinstance(payload, dict) else []:
            if not isinstance(item, dict):
                continue
            oid = str(item.get("id") or "")
            title = clean(item.get("title"))
            abs_text = abstract(item.get("abstract_inverted_index"))
            text = f"{title} {abs_text}"
            if name.lower() not in text.lower() or (name, oid) in seen:
                continue
            direct = bool(DIRECT.search(text)); within = bool(WITHIN.search(text)); among = bool(GEOGRAPHIC.search(text)); artificial = bool(NEGATIVE.search(text))
            if not direct and not within and not among:
                continue
            score = (12 if name.lower() in title.lower() else 5) + (10 if direct else 0) + (8 if within else 0) + (8 if among else 0) - (10 if artificial else 0)
            if score < 10:
                continue
            seen.add((name, oid))
            primary = item.get("primary_location") or {}
            landing = str(primary.get("landing_page_url") or oid) if isinstance(primary, dict) else oid
            doi = str(item.get("doi") or "").replace("https://doi.org/", "")
            snippet = abs_text[:1000] if abs_text else title
            out_rows.append({
                "canonical_name": name,
                "family": row.get("family", ""),
                "openalex_id": oid,
                "title": title,
                "year": item.get("publication_year") or "",
                "doi": doi,
                "landing_url": landing,
                "query_mode": "exact_species_spatial_scale_enrichment",
                "within_signal": int(within),
                "among_signal": int(among),
                "direct_colour_signal": int(direct),
                "artificial_signal": int(artificial),
                "score": score,
                "evidence_snippet": clean(snippet),
            })

    out_path = Path(args.out); out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS); writer.writeheader(); writer.writerows(sorted(out_rows, key=lambda r: (str(r["canonical_name"]), -int(r["score"]))))

    by_species: dict[str, dict[str, object]] = {}
    for target in targets:
        name = target.get("canonical_name", "")
        works = [r for r in out_rows if r["canonical_name"] == name]
        w = any(int(r["within_signal"]) for r in works); g = any(int(r["among_signal"]) for r in works)
        enriched = "mixed" if w and g else "within_population" if w else "among_population" if g else baseline_scale(target)
        by_species[name] = {
            "canonical_name": name,
            "family": target.get("family", ""),
            "baseline_scale": baseline_scale(target),
            "enriched_scale": enriched,
            "n_enrichment_works": len(works),
            "n_within_signal_works": sum(int(r["within_signal"]) for r in works),
            "n_among_signal_works": sum(int(r["among_signal"]) for r in works),
            "best_enrichment_score": max([int(r["score"]) for r in works], default=0),
            "manual_review_required": enriched in {"mixed", "unclear"},
        }
    summary_fields = list(next(iter(by_species.values())).keys()) if by_species else []
    with Path(args.species_summary).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields); writer.writeheader(); writer.writerows(by_species.values())

    qc = {
        "confirmed_species": len(review),
        "ambiguous_targets": len(targets),
        "targets_with_retained_works": len({str(r["canonical_name"]) for r in out_rows}),
        "retained_works": len(out_rows),
        "baseline_counts": dict(__import__("collections").Counter(baseline_scale(r) for r in review)),
        "enriched_target_counts": dict(__import__("collections").Counter(str(r["enriched_scale"]) for r in by_species.values())),
        "failed_requests": len(errors),
        "errors": errors,
        "semantic_guard": "enrichment applies only to already-confirmed species and does not create new FCP cases",
        "complete": len(errors) == 0,
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False))


if __name__ == "__main__":
    main()
