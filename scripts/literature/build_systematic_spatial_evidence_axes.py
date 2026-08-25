#!/usr/bin/env python3
"""Build a taxon-validated, mixed-preserving evidence-axis review queue.

This script is intentionally upstream of ecological analysis. Automated record flags are
used for navigation only; final natural eligibility and spatial evidence axes remain human
review fields until adjudication.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

TRUTHY = {"1", "true", "yes"}
BINOMIAL = re.compile(r"^[A-Z][a-z-]{2,} [a-z][a-z-]{2,}$")
DISCOVERY_PRIORITIES = {
    "P1_high_natural_itv",
    "P1_high_population_itv",
    "P2_possible_itv",
    "P3_likely_exclusion_review",
}
PRIORITY_ORDER = {
    "P1_high_natural_itv": 0,
    "P1_high_population_itv": 1,
    "P2_possible_itv": 2,
    "P3_likely_exclusion_review": 3,
    "P4_flower_colour_context_only": 4,
    "P5_low_relevance": 5,
}


def truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in TRUTHY


def clean(value: Any) -> str:
    return " ".join(str(value or "").split())


def integer(value: Any) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def contextual_candidate_names(row: dict[str, str]) -> list[str]:
    """Return full-binomial candidates explicitly occurring in title or abstract.

    Candidate strings remain provisional until GBIF validation. The whole abstract is
    searched; the earlier exploratory code only inspected its first 1,200 characters.
    """
    text = f"{row.get('title', '')} {row.get('abstract', '')}".lower()
    names: list[str] = []
    seen: set[str] = set()
    for token in str(row.get("candidate_species") or "").split(";"):
        name = clean(token)
        if not BINOMIAL.fullmatch(name) or name in seen:
            continue
        if name.lower() in text:
            seen.add(name)
            names.append(name)
    return names


def state_from_counts(within_records: int, among_records: int) -> str:
    if within_records > 0 and among_records > 0:
        return "mixed_evidence"
    if within_records > 0:
        return "within_evidence_only"
    if among_records > 0:
        return "among_evidence_only"
    return "unresolved"


def get_json(url: str, timeout: int, retries: int) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "fcp-upstream-spatial-reaudit/1.0",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError("GBIF response was not a JSON object")
            return payload
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(url) from last


def resolve_gbif_name(name: str, timeout: int, retries: int) -> dict[str, Any]:
    params = urllib.parse.urlencode({"name": name, "kingdom": "Plantae", "strict": "true"})
    match = get_json("https://api.gbif.org/v1/species/match?" + params, timeout, retries)
    match_type = str(match.get("matchType") or "")
    rank = str(match.get("rank") or "")
    kingdom = str(match.get("kingdom") or "")
    confidence = integer(match.get("confidence"))
    accepted = (
        match_type == "EXACT"
        or (match_type == "CONFIDENCE" and confidence >= 95)
    ) and rank == "SPECIES" and kingdom == "Plantae"

    if not accepted:
        return {
            "input_name": name,
            "accepted": False,
            "match_type": match_type,
            "rank": rank,
            "kingdom": kingdom,
            "confidence": confidence,
            "accepted_name": "",
            "family": str(match.get("family") or ""),
            "usage_key": match.get("acceptedUsageKey") or match.get("usageKey") or "",
            "reason": "not_strict_plant_species_match",
        }

    usage_key = match.get("acceptedUsageKey") or match.get("usageKey")
    accepted_payload = match
    if match.get("acceptedUsageKey") and match.get("acceptedUsageKey") != match.get("usageKey"):
        try:
            accepted_payload = get_json(
                f"https://api.gbif.org/v1/species/{int(match['acceptedUsageKey'])}",
                timeout,
                retries,
            )
        except Exception:  # keep the strict match if accepted-name fetch fails
            accepted_payload = match

    accepted_name = clean(
        accepted_payload.get("canonicalName")
        or accepted_payload.get("scientificName")
        or match.get("canonicalName")
        or name
    )
    accepted_name = " ".join(accepted_name.split()[:2])
    family = clean(accepted_payload.get("family") or match.get("family"))
    accepted = bool(family) and bool(BINOMIAL.fullmatch(accepted_name))
    return {
        "input_name": name,
        "accepted": accepted,
        "match_type": match_type,
        "rank": rank,
        "kingdom": kingdom,
        "confidence": confidence,
        "accepted_name": accepted_name if accepted else "",
        "family": family,
        "usage_key": usage_key or "",
        "reason": "accepted_species" if accepted else "accepted_name_or_family_unresolved",
    }


def best_record(records: Iterable[dict[str, str]], axis: str | None = None) -> dict[str, str] | None:
    rows = list(records)
    if axis == "within":
        rows = [row for row in rows if truthy(row.get("within_signal"))]
    elif axis == "among":
        rows = [row for row in rows if truthy(row.get("among_signal"))]
    if not rows:
        return None
    return max(
        rows,
        key=lambda row: (
            -PRIORITY_ORDER.get(str(row.get("screen_priority")), 99),
            truthy(row.get("natural_signal")),
            integer(row.get("cited_by_count")),
            len(str(row.get("abstract") or "")),
        ),
    )


def source_id(row: dict[str, str] | None) -> str:
    if not row:
        return ""
    return clean(row.get("doi") or row.get("record_id") or row.get("url"))


def excerpt(row: dict[str, str] | None, limit: int = 1600) -> str:
    if not row:
        return ""
    return clean(row.get("abstract"))[:limit]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--priorities", default=",".join(sorted(DISCOVERY_PRIORITIES)))
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--delay", type=float, default=0.03)
    parser.add_argument("--gbif-cache", default="")
    parser.add_argument("--frozen-manifest", default="docs/supporting/frozen_classification_manifest.csv")
    args = parser.parse_args()

    selected_priorities = {x.strip() for x in args.priorities.split(",") if x.strip()}
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    diag: Counter[str] = Counter()
    record_links: list[dict[str, Any]] = []
    names_to_records: dict[str, list[dict[str, str]]] = defaultdict(list)

    with Path(args.queue).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            diag["records_total"] += 1
            priority = str(row.get("screen_priority") or "")
            diag[f"priority_{priority}"] += 1
            if priority not in selected_priorities:
                continue
            diag["records_in_discovery_priorities"] += 1
            names = contextual_candidate_names(row)
            if not names:
                diag["records_without_contextual_candidate"] += 1
                continue
            diag["records_with_contextual_candidate"] += 1
            for name in names:
                names_to_records[name].append(row)
                record_links.append({
                    "input_name": name,
                    "record_id": row.get("record_id", ""),
                    "doi": row.get("doi", ""),
                    "title": row.get("title", ""),
                    "screen_priority": priority,
                    "within_signal": int(truthy(row.get("within_signal"))),
                    "among_signal": int(truthy(row.get("among_signal"))),
                    "natural_signal": int(truthy(row.get("natural_signal"))),
                    "cultivated_signal": int(truthy(row.get("cultivated_signal"))),
                    "induced_signal": int(truthy(row.get("induced_signal"))),
                    "ontogenetic_signal": int(truthy(row.get("ontogenetic_signal"))),
                    "non_display_floral_signal": int(truthy(row.get("non_display_floral_signal"))),
                    "provisional_variation_form": row.get("provisional_variation_form", ""),
                    "evidence_excerpt": excerpt(row),
                })

    diag["unique_contextual_input_names"] = len(names_to_records)

    cache_path = Path(args.gbif_cache) if args.gbif_cache else outdir / "gbif_resolution_cache.json"
    cache: dict[str, dict[str, Any]] = {}
    if cache_path.exists():
        try:
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                cache = payload
        except Exception:
            cache = {}

    resolution_rows: list[dict[str, Any]] = []
    for index, name in enumerate(sorted(names_to_records), start=1):
        if name in cache:
            resolved = cache[name]
            diag["gbif_cache_hits"] += 1
        else:
            try:
                resolved = resolve_gbif_name(name, args.timeout, args.retries)
            except Exception as exc:  # noqa: BLE001
                resolved = {
                    "input_name": name,
                    "accepted": False,
                    "match_type": "",
                    "rank": "",
                    "kingdom": "",
                    "confidence": 0,
                    "accepted_name": "",
                    "family": "",
                    "usage_key": "",
                    "reason": f"resolution_error:{type(exc).__name__}",
                }
            cache[name] = resolved
            if args.delay:
                time.sleep(args.delay)
        resolution_rows.append(resolved)
        if index % 100 == 0:
            cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    accepted_by_input = {
        row["input_name"]: row
        for row in resolution_rows
        if row.get("accepted") and row.get("accepted_name")
    }
    diag["gbif_accepted_input_names"] = len(accepted_by_input)
    diag["gbif_rejected_input_names"] = len(resolution_rows) - len(accepted_by_input)

    species_records: dict[str, list[dict[str, str]]] = defaultdict(list)
    species_inputs: dict[str, set[str]] = defaultdict(set)
    species_family: dict[str, str] = {}
    for input_name, records in names_to_records.items():
        resolved = accepted_by_input.get(input_name)
        if not resolved:
            continue
        accepted_name = str(resolved["accepted_name"])
        species_inputs[accepted_name].add(input_name)
        species_family[accepted_name] = str(resolved.get("family") or "")
        species_records[accepted_name].extend(records)

    species_rows: list[dict[str, Any]] = []
    for accepted_name in sorted(species_records):
        unique: dict[str, dict[str, str]] = {}
        for row in species_records[accepted_name]:
            key = clean(row.get("doi") or row.get("record_id") or row.get("title"))
            if not key:
                continue
            previous = unique.get(key)
            if previous is None or PRIORITY_ORDER.get(str(row.get("screen_priority")), 99) < PRIORITY_ORDER.get(str(previous.get("screen_priority")), 99):
                unique[key] = row
        records = list(unique.values())
        within_count = sum(truthy(row.get("within_signal")) for row in records)
        among_count = sum(truthy(row.get("among_signal")) for row in records)
        natural_count = sum(truthy(row.get("natural_signal")) for row in records)
        cultivated_count = sum(truthy(row.get("cultivated_signal")) for row in records)
        induced_count = sum(truthy(row.get("induced_signal")) for row in records)
        ontogenetic_count = sum(truthy(row.get("ontogenetic_signal")) for row in records)
        non_display_count = sum(truthy(row.get("non_display_floral_signal")) for row in records)
        discrete_count = sum(str(row.get("provisional_variation_form")) == "discrete_signal" for row in records)
        best_any = best_record(records)
        best_within = best_record(records, "within")
        best_among = best_record(records, "among")
        state = state_from_counts(within_count, among_count)
        species_rows.append({
            "canonical_name": accepted_name,
            "family": species_family.get(accepted_name, ""),
            "input_names": ";".join(sorted(species_inputs[accepted_name])),
            "n_supporting_records": len(records),
            "within_signal_records": within_count,
            "among_signal_records": among_count,
            "natural_signal_records": natural_count,
            "cultivated_signal_records": cultivated_count,
            "induced_signal_records": induced_count,
            "ontogenetic_signal_records": ontogenetic_count,
            "non_display_floral_signal_records": non_display_count,
            "discrete_signal_records": discrete_count,
            "automated_evidence_state": state,
            "automated_eligibility_conflict": int(
                (natural_count > 0 and (cultivated_count + induced_count + ontogenetic_count > 0))
                or (natural_count == 0)
            ),
            "best_source_id": source_id(best_any),
            "best_title": best_any.get("title", "") if best_any else "",
            "best_evidence_excerpt": excerpt(best_any),
            "best_within_source_id": source_id(best_within),
            "best_within_title": best_within.get("title", "") if best_within else "",
            "best_within_evidence_excerpt": excerpt(best_within),
            "best_among_source_id": source_id(best_among),
            "best_among_title": best_among.get("title", "") if best_among else "",
            "best_among_evidence_excerpt": excerpt(best_among),
            "review_status": "unreviewed",
            "reviewer_1_natural_eligibility": "",
            "reviewer_1_local_coexistence": "",
            "reviewer_1_geographic_structure": "",
            "reviewer_1_variation_form": "",
            "reviewer_1_notes": "",
            "reviewer_2_natural_eligibility": "",
            "reviewer_2_local_coexistence": "",
            "reviewer_2_geographic_structure": "",
            "reviewer_2_variation_form": "",
            "reviewer_2_notes": "",
            "adjudicated_natural_eligibility": "",
            "adjudicated_local_coexistence": "",
            "adjudicated_geographic_structure": "",
            "adjudicated_variation_form": "",
            "adjudicated_state": "",
            "adjudication_notes": "",
        })

    state_counts = Counter(row["automated_evidence_state"] for row in species_rows)
    diag["gbif_accepted_species"] = len(species_rows)
    for state, count in state_counts.items():
        diag[f"state_{state}"] = count

    frozen_rows: list[dict[str, str]] = []
    frozen_path = Path(args.frozen_manifest)
    if frozen_path.exists():
        with frozen_path.open(newline="", encoding="utf-8") as handle:
            frozen_rows = list(csv.DictReader(handle))
    species_by_name = {row["canonical_name"]: row for row in species_rows}
    recovery_rows: list[dict[str, Any]] = []
    for row in frozen_rows:
        name = row.get("canonical_name", "")
        new = species_by_name.get(name)
        recovery_rows.append({
            "canonical_name": name,
            "historical_spatial_scale": row.get("spatial_scale", ""),
            "recovered_in_systematic_candidate_queue": int(new is not None),
            "automated_evidence_state": new.get("automated_evidence_state", "") if new else "",
            "within_signal_records": new.get("within_signal_records", "") if new else "",
            "among_signal_records": new.get("among_signal_records", "") if new else "",
            "note": "diagnostic_only_not_an_adjudication",
        })
    diag["historical_34_recovered"] = sum(row["recovered_in_systematic_candidate_queue"] for row in recovery_rows)

    write_csv(outdir / "systematic_species_review_queue.csv", species_rows)
    write_csv(outdir / "systematic_record_evidence_links.csv", record_links)
    write_csv(outdir / "gbif_candidate_resolution_audit.csv", resolution_rows)
    write_csv(outdir / "historical_34_recovery_diagnostic.csv", recovery_rows)

    qc = {
        "status": "complete",
        "source_queue": str(args.queue),
        "selected_priorities": sorted(selected_priorities),
        "diagnostics": dict(sorted(diag.items())),
        "automated_state_counts": dict(sorted(state_counts.items())),
        "semantic_guard": (
            "Automated evidence states are navigation aids only. Taxon validation precedes state aggregation; "
            "mixed evidence is retained; natural eligibility and both spatial axes require blinded source-level review."
        ),
    }
    (outdir / "upstream_reaudit_qc.json").write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False, indent=2))

    if not species_rows:
        raise SystemExit("No GBIF-validated candidate species were produced")
    if state_counts.get("mixed_evidence", 0) == 0:
        raise SystemExit("Mixed evidence disappeared entirely; inspect extraction/resolution before proceeding")


if __name__ == "__main__":
    main()
