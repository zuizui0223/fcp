#!/usr/bin/env python3
"""Audit historical benchmark recovery after one-generation OpenAlex citation chasing.

Starts from the seven prespecified review seeds used by the systematic-map protocol,
collects each seed, its OpenAlex references, and all direct citers, then tests exact
DOI/OpenAlex recovery of the historical 34 classification sources. This is a search
completeness diagnostic only; it does not classify biological eligibility or spatial state.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

USER_AGENT = "fcp-jbi-v2-citation-audit/1.0 (https://github.com/zuizui0223/fcp)"
SEEDS = [
    "10.1111/plb.12575",
    "10.1016/j.tree.2021.01.011",
    "10.1086/705589",
    "10.3389/fpls.2021.617851",
    "10.1073/pnas.97.13.7016",
    "10.1046/j.1469-8137.2001.00159.x",
    "10.7818/ECOS.2014.23-3.06",
]


def clean(v: Any) -> str:
    return " ".join(str(v or "").split()).strip()


def norm_doi(v: Any) -> str:
    s = clean(v).lower()
    s = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", s)
    s = re.sub(r"^doi:\s*", "", s)
    return s.rstrip(".")


def norm_oa(v: Any) -> str:
    s = clean(v).rstrip("/")
    if "openalex.org/" in s.lower():
        return s.rsplit("/", 1)[-1].upper()
    if re.fullmatch(r"W\d+", s, re.I):
        return s.upper()
    return ""


def request(url: str, timeout: int, retries: int) -> dict[str, Any]:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = json.load(r)
            if not isinstance(data, dict):
                raise RuntimeError("non-object response")
            return data
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(url) from last


def works_by_ids(ids: list[str], api_key: str, timeout: int, retries: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for i in range(0, len(ids), 50):
        batch = [norm_oa(x) for x in ids[i:i+50] if norm_oa(x)]
        if not batch:
            continue
        params = {"filter": "openalex:" + "|".join(batch), "per-page": 100, "select": "id,doi,title"}
        if api_key: params["api_key"] = api_key
        data = request("https://api.openalex.org/works?" + urllib.parse.urlencode(params), timeout, retries)
        results = data.get("results")
        if isinstance(results, list): out.extend(x for x in results if isinstance(x, dict))
    return out


def resolve_seed(doi: str, api_key: str, timeout: int, retries: int) -> dict[str, Any] | None:
    params = {"filter": "doi:" + doi, "per-page": 1, "select": "id,doi,title,referenced_works"}
    if api_key: params["api_key"] = api_key
    data = request("https://api.openalex.org/works?" + urllib.parse.urlencode(params), timeout, retries)
    results = data.get("results")
    return results[0] if isinstance(results, list) and results and isinstance(results[0], dict) else None


def citers(seed_id: str, api_key: str, timeout: int, retries: int) -> list[dict[str, Any]]:
    short = norm_oa(seed_id)
    cursor = "*"
    out: list[dict[str, Any]] = []
    while True:
        params = {"filter": f"cites:{short},is_retracted:false", "per-page": 100, "cursor": cursor, "select": "id,doi,title"}
        if api_key: params["api_key"] = api_key
        data = request("https://api.openalex.org/works?" + urllib.parse.urlencode(params), timeout, retries)
        results = data.get("results")
        if not isinstance(results, list) or not results: break
        out.extend(x for x in results if isinstance(x, dict))
        meta = data.get("meta") if isinstance(data.get("meta"), dict) else {}
        nxt = meta.get("next_cursor")
        if not nxt or str(nxt) == cursor: break
        cursor = str(nxt)
        time.sleep(0.03)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--historical-manifest", default="docs/supporting/frozen_classification_manifest.csv")
    p.add_argument("--direct-recovery", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    p.add_argument("--timeout", type=int, default=45)
    p.add_argument("--retries", type=int, default=3)
    a = p.parse_args()
    api_key = os.environ.get(a.api_key_env, "").strip()

    with Path(a.historical_manifest).open(newline="", encoding="utf-8") as h:
        hist = list(csv.DictReader(h))
    with Path(a.direct_recovery).open(newline="", encoding="utf-8") as h:
        direct = list(csv.DictReader(h))
    direct_by_name = {r["canonical_name"]: int(r["recovered_in_openalex_v2_direct_queries"]) for r in direct}

    neighborhood: dict[str, dict[str, Any]] = {}
    seed_logs = []
    for doi in SEEDS:
        seed = resolve_seed(doi, api_key, a.timeout, a.retries)
        if not seed:
            seed_logs.append({"seed_doi": doi, "resolved": False, "backward": 0, "forward": 0})
            continue
        refs = [str(x) for x in (seed.get("referenced_works") or [])]
        back = works_by_ids(refs, api_key, a.timeout, a.retries)
        forward = citers(str(seed.get("id") or ""), api_key, a.timeout, a.retries)
        all_items = [seed] + back + forward
        for item in all_items:
            oid = norm_oa(item.get("id"))
            doi2 = norm_doi(item.get("doi"))
            key = oid or ("DOI:" + doi2 if doi2 else "")
            if key: neighborhood[key] = item
        seed_logs.append({"seed_doi": doi, "resolved": True, "backward": len(back), "forward": len(forward)})
        print(seed_logs[-1], flush=True)

    neighborhood_dois = {norm_doi(x.get("doi")) for x in neighborhood.values() if norm_doi(x.get("doi"))}
    neighborhood_oa = {norm_oa(x.get("id")) for x in neighborhood.values() if norm_oa(x.get("id"))}
    rows = []
    for r in hist:
        name = clean(r.get("canonical_name"))
        sid = clean(r.get("source_id"))
        doi = norm_doi(sid)
        oa = norm_oa(sid)
        cited = (doi.startswith("10.") and doi in neighborhood_dois) or (bool(oa) and oa in neighborhood_oa)
        direct_hit = direct_by_name.get(name, 0)
        rows.append({
            "canonical_name": name,
            "historical_source_id": sid,
            "direct_v2_recovered": direct_hit,
            "citation_chase_recovered": int(cited),
            "direct_or_citation_recovered": int(bool(direct_hit or cited)),
        })

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    total_direct = sum(r["direct_v2_recovered"] for r in rows)
    total_cited = sum(r["citation_chase_recovered"] for r in rows)
    total_union = sum(r["direct_or_citation_recovered"] for r in rows)
    summary = {
        "status": "complete",
        "seed_reviews": len(SEEDS),
        "resolved_seeds": sum(bool(x["resolved"]) for x in seed_logs),
        "seed_logs": seed_logs,
        "unique_citation_neighborhood_works": len(neighborhood),
        "historical_direct_recovered": total_direct,
        "historical_citation_recovered": total_cited,
        "historical_direct_or_citation_recovered": total_union,
        "remaining_historical_benchmark_misses": [r["canonical_name"] for r in rows if not r["direct_or_citation_recovered"]],
        "semantic_guard": "Benchmark recovery diagnoses search recall only and does not validate biological eligibility or spatial coding."
    }
    Path(a.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == "__main__": main()
