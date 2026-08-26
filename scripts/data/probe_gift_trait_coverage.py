#!/usr/bin/env python3
"""Audit GIFT trait coverage for the focal-consistent FCP core.

This script calls the public GIFT 3.2 API directly and never interprets missing
trait records as biological absence. It resolves focal names through GIFT's
name-matching table, identifies a small hypothesis-driven trait set from the
trait metadata table, retrieves aggregated species-level values, and writes
coverage/QC tables for downstream C/S analyses.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

API = "https://gift.uni-goettingen.de/api/extended/index3.2.php"
TARGET_PATTERNS = {
    "self_fertilization": [r"self[- _]?fertili[sz]ation"],
    "lifecycle": [r"^life ?cycle$", r"lifecycle"],
    "dispersal_syndrome": [r"dispersal.*syndrome", r"dispersal syndrome"],
    "flowering_start": [r"flowering.*start", r"start.*flower"],
    "flowering_end": [r"flowering.*end", r"end.*flower"],
}


def get_json(params: dict[str, Any], retries: int = 4, timeout: int = 60) -> Any:
    url = API + "?" + urlencode(params)
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = Request(url, headers={"Accept": "application/json", "User-Agent": "fcp-jbi-gift-coverage/1.0"})
            with urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"GIFT request failed after {retries} attempts: {url}") from last


def choose_traits(meta: pd.DataFrame) -> dict[str, dict[str, str]]:
    chosen: dict[str, dict[str, str]] = {}
    for key, patterns in TARGET_PATTERNS.items():
        hits = []
        for row in meta.itertuples(index=False):
            label = str(getattr(row, "Trait2", "") or "").strip()
            text = label.lower()
            if any(re.search(p, text, flags=re.I) for p in patterns):
                hits.append(row)
        if not hits:
            raise RuntimeError(f"No GIFT trait metadata match for {key}")
        hits.sort(key=lambda r: float(getattr(r, "count", 0) or 0), reverse=True)
        row = hits[0]
        chosen[key] = {
            "trait_ID": str(getattr(row, "Lvl3")),
            "Trait2": str(getattr(row, "Trait2")),
            "type": str(getattr(row, "type")),
            "count": str(getattr(row, "count")),
        }
    return chosen


def resolve_name(name: str) -> dict[str, Any]:
    parts = name.split()
    if len(parts) < 2:
        return {"canonical_name": name, "status": "invalid_binomial"}
    genus, epithet = parts[0], parts[1]
    rows = get_json({"query": "names_matched_unique", "genus": genus, "epithet": epithet})
    if not rows:
        return {"canonical_name": name, "status": "not_found"}
    q = pd.DataFrame(rows)
    if "work_species" in q.columns:
        exact = q.loc[q.work_species.astype(str).str.casefold().eq(name.casefold())]
        if len(exact):
            q = exact
    for col in ["overallscore", "accepted", "resolved", "matched"]:
        if col in q.columns:
            q[col] = pd.to_numeric(q[col], errors="coerce")
    sortcols = [c for c in ["overallscore", "accepted", "resolved", "matched"] if c in q.columns]
    if sortcols:
        q = q.sort_values(sortcols, ascending=False, na_position="last")
    r = q.iloc[0]
    work_id = pd.to_numeric(pd.Series([r.get("work_ID")]), errors="coerce").iloc[0]
    return {
        "canonical_name": name,
        "status": "resolved" if pd.notna(work_id) else "unresolved",
        "work_ID": int(work_id) if pd.notna(work_id) else None,
        "gift_work_species": str(r.get("work_species") or ""),
        "overallscore": r.get("overallscore"),
        "accepted": r.get("accepted"),
        "resolved": r.get("resolved"),
        "synonym": r.get("synonym"),
    }


def fetch_trait(trait_id: str, expected_count: int) -> pd.DataFrame:
    pages = []
    start = 0
    max_pages = max(1, int(math.ceil(max(expected_count, 1) / 10000)) + 2)
    for _ in range(max_pages):
        rows = get_json({"query": "traits", "traitid": trait_id, "biasref": 0, "biasderiv": 0, "startat": start})
        if not rows:
            break
        q = pd.DataFrame(rows)
        if q.empty:
            break
        pages.append(q)
        if len(q) < 10000:
            break
        start += 10000
    return pd.concat(pages, ignore_index=True) if pages else pd.DataFrame()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--core-species", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    core = pd.read_csv(a.core_species)
    if len(core) != 74:
        raise SystemExit(f"Expected 74 v6 core species, got {len(core)}")

    meta = pd.DataFrame(get_json({"query": "traits_meta"}))
    meta["count"] = pd.to_numeric(meta.get("count"), errors="coerce")
    chosen = choose_traits(meta)
    (out / "gift_target_trait_metadata.json").write_text(json.dumps(chosen, indent=2) + "\n")

    resolved = pd.DataFrame([resolve_name(str(x)) for x in core.canonical_name])
    resolved.to_csv(out / "gift_core_name_resolution.csv", index=False)
    work_to_name = {
        int(r.work_ID): r.canonical_name
        for r in resolved.itertuples(index=False)
        if getattr(r, "status") == "resolved" and pd.notna(getattr(r, "work_ID", None))
    }

    joined = core[["canonical_name", "family", "organization_state"]].copy()
    raw_rows = []
    coverage_rows = []
    for key, info in chosen.items():
        q = fetch_trait(info["trait_ID"], int(float(info["count"] or 0)))
        if q.empty:
            joined[key] = pd.NA
            coverage_rows.append({"trait": key, "trait_ID": info["trait_ID"], "matched_core_species": 0, "matched_informative_species": 0})
            continue
        q["work_ID"] = pd.to_numeric(q["work_ID"], errors="coerce")
        q = q.loc[q.work_ID.isin(work_to_name)].copy()
        q["canonical_name"] = q.work_ID.astype(int).map(work_to_name)
        if "agreement" in q.columns:
            ag = pd.to_numeric(q.agreement, errors="coerce")
            q = q.loc[ag.isna() | ag.ge(0.66)].copy()
        q["target_trait"] = key
        raw_rows.append(q)
        qq = q.sort_values(["canonical_name"]).drop_duplicates("canonical_name", keep="first")
        mapping = dict(zip(qq.canonical_name.astype(str), qq.trait_value))
        joined[key] = joined.canonical_name.map(mapping)
        informative = joined.organization_state.isin(["local_coexistence_only", "spatial_segregation_only", "coexistence_and_segregation"])
        coverage_rows.append({
            "trait": key,
            "trait_ID": info["trait_ID"],
            "gift_label": info["Trait2"],
            "gift_type": info["type"],
            "matched_core_species": int(joined[key].notna().sum()),
            "matched_informative_species": int(joined.loc[informative, key].notna().sum()),
            "C_only_n": int(joined.loc[joined.organization_state.eq("local_coexistence_only"), key].notna().sum()),
            "S_only_n": int(joined.loc[joined.organization_state.eq("spatial_segregation_only"), key].notna().sum()),
            "mixed_n": int(joined.loc[joined.organization_state.eq("coexistence_and_segregation"), key].notna().sum()),
            "unique_nonmissing_values": int(joined[key].dropna().astype(str).nunique()),
        })

    joined.to_csv(out / "gift_core_trait_coverage.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(out / "gift_core_trait_coverage_summary.csv", index=False)
    if raw_rows:
        pd.concat(raw_rows, ignore_index=True).to_csv(out / "gift_core_trait_records.csv", index=False)
    qc = {
        "status": "complete",
        "gift_version": "3.2",
        "api": API,
        "core_species": int(len(core)),
        "resolved_names": int(resolved.status.eq("resolved").sum()),
        "target_traits": chosen,
        "missing_semantics": "no GIFT value is treated as missing documentation, never biological absence",
        "categorical_agreement_threshold": 0.66,
        "bias_ref": False,
        "bias_deriv": False,
    }
    (out / "gift_core_trait_coverage_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    print(pd.DataFrame(coverage_rows).to_csv(index=False))


if __name__ == "__main__":
    main()
