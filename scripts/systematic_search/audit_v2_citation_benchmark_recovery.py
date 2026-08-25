#!/usr/bin/env python3
"""Audit historical benchmark recovery through one-generation seed citation chasing.

Instead of enumerating every citer of every seed, resolve only the seven prespecified
seed reviews and the 34 historical benchmark sources. A historical source is reachable
by one-generation chasing if it is a seed, is referenced by a seed (backward), or cites
a seed (forward). This is exactly the benchmark-recovery question and avoids huge citer
enumeration. No biological classification is performed.
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

USER_AGENT = "fcp-jbi-v2-citation-audit/2.0 (https://github.com/zuizui0223/fcp)"
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


def resolve_by_doi(doi: str, api_key: str, timeout: int, retries: int) -> dict[str, Any] | None:
    params = {"filter": "doi:" + doi, "per-page": 1, "select": "id,doi,title,referenced_works"}
    if api_key:
        params["api_key"] = api_key
    data = request("https://api.openalex.org/works?" + urllib.parse.urlencode(params), timeout, retries)
    results = data.get("results")
    return results[0] if isinstance(results, list) and results and isinstance(results[0], dict) else None


def resolve_by_oa(oa: str, api_key: str, timeout: int, retries: int) -> dict[str, Any] | None:
    params = {"filter": "openalex:" + oa, "per-page": 1, "select": "id,doi,title,referenced_works"}
    if api_key:
        params["api_key"] = api_key
    data = request("https://api.openalex.org/works?" + urllib.parse.urlencode(params), timeout, retries)
    results = data.get("results")
    return results[0] if isinstance(results, list) and results and isinstance(results[0], dict) else None


def resolve_source(source_id: str, api_key: str, timeout: int, retries: int) -> dict[str, Any] | None:
    doi = norm_doi(source_id)
    if doi.startswith("10."):
        return resolve_by_doi(doi, api_key, timeout, retries)
    oa = norm_oa(source_id)
    if oa:
        return resolve_by_oa(oa, api_key, timeout, retries)
    return None


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
    if len(hist) != 34:
        raise SystemExit("Historical manifest must have 34 rows")
    direct_by_name = {r["canonical_name"]: int(r["recovered_in_openalex_v2_direct_queries"]) for r in direct}

    seed_works: dict[str, dict[str, Any]] = {}
    seed_logs = []
    for doi in SEEDS:
        work = resolve_by_doi(doi, api_key, a.timeout, a.retries)
        if work:
            seed_works[norm_oa(work.get("id"))] = work
            seed_logs.append({"seed_doi": doi, "resolved": True, "openalex_id": norm_oa(work.get("id")), "references": len(work.get("referenced_works") or [])})
        else:
            seed_logs.append({"seed_doi": doi, "resolved": False, "openalex_id": "", "references": 0})
        print(seed_logs[-1], flush=True)

    seed_ids = set(seed_works)
    seed_reference_ids = {norm_oa(ref) for work in seed_works.values() for ref in (work.get("referenced_works") or []) if norm_oa(ref)}

    rows = []
    unresolved_targets = []
    for r in hist:
        name = clean(r.get("canonical_name"))
        sid = clean(r.get("source_id"))
        target = resolve_source(sid, api_key, a.timeout, a.retries)
        direct_hit = direct_by_name.get(name, 0)
        if target:
            target_id = norm_oa(target.get("id"))
            target_refs = {norm_oa(x) for x in (target.get("referenced_works") or []) if norm_oa(x)}
            is_seed = target_id in seed_ids
            backward = target_id in seed_reference_ids
            forward = bool(target_refs.intersection(seed_ids))
            citation_hit = bool(is_seed or backward or forward)
            relation = "seed" if is_seed else "backward_from_seed" if backward else "forward_cites_seed" if forward else "none"
        else:
            target_id = ""
            citation_hit = False
            relation = "target_unresolved"
            unresolved_targets.append(name)
        rows.append({
            "canonical_name": name,
            "historical_source_id": sid,
            "target_openalex_id": target_id,
            "direct_v2_recovered": direct_hit,
            "citation_chase_recovered": int(citation_hit),
            "citation_relation": relation,
            "direct_or_citation_recovered": int(bool(direct_hit or citation_hit)),
        })
        time.sleep(0.02)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    total_direct = sum(r["direct_v2_recovered"] for r in rows)
    total_cited = sum(r["citation_chase_recovered"] for r in rows)
    total_union = sum(r["direct_or_citation_recovered"] for r in rows)
    relation_counts: dict[str, int] = {}
    for row in rows:
        relation_counts[row["citation_relation"]] = relation_counts.get(row["citation_relation"], 0) + 1
    summary = {
        "status": "complete",
        "seed_reviews": len(SEEDS),
        "resolved_seeds": sum(bool(x["resolved"]) for x in seed_logs),
        "seed_logs": seed_logs,
        "historical_targets_resolved_in_openalex": 34 - len(unresolved_targets),
        "unresolved_historical_targets": unresolved_targets,
        "citation_relation_counts": relation_counts,
        "historical_direct_recovered": total_direct,
        "historical_citation_recovered": total_cited,
        "historical_direct_or_citation_recovered": total_union,
        "remaining_historical_benchmark_misses": [r["canonical_name"] for r in rows if not r["direct_or_citation_recovered"]],
        "semantic_guard": "Benchmark recovery diagnoses search recall only and does not validate biological eligibility or spatial coding."
    }
    Path(a.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
