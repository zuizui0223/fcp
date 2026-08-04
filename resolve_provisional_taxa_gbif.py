#!/usr/bin/env python3
"""Resolve provisional candidate names with GBIF and reaggregate evidence by accepted species.

The important invariant is that spatial-scale labels are never selected before synonym
resolution.  All provisional rows mapping to the same accepted GBIF species are combined;
accepted species supported by both within- and among-population evidence are written to
the audit as mixed and excluded from the binary ecological analysis.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any


def get_json(url: str, timeout: int = 45, retries: int = 4) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "fcp-exploratory-taxon-resolver/2.0",
                },
            )
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError("GBIF response was not a JSON object")
            return payload
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(url) from last


def integer(value: Any) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def stable_source_id(row: dict[str, str]) -> str:
    return str(row.get("best_doi") or row.get("best_openalex_id") or row.get("best_title") or "").strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args()

    input_rows = list(csv.DictReader(Path(args.input).open(encoding="utf-8")))
    resolved_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    audit: list[dict[str, Any]] = []

    for row in input_rows:
        input_name = row["canonical_name"].strip()
        params = urllib.parse.urlencode({"name": input_name, "kingdom": "Plantae", "strict": "true"})
        try:
            match = get_json(
                "https://api.gbif.org/v1/species/match?" + params,
                args.timeout,
                args.retries,
            )
        except Exception as exc:  # noqa: BLE001
            audit.append({
                "input_name": input_name,
                "accepted": False,
                "reason": type(exc).__name__,
                "post_resolution_status": "resolution_error",
            })
            continue

        match_type = str(match.get("matchType") or "")
        rank = str(match.get("rank") or "")
        kingdom = str(match.get("kingdom") or "")
        family = str(match.get("family") or "")
        accepted_name = " ".join(
            str(match.get("canonicalName") or match.get("scientificName") or input_name).split()[:2]
        )
        accepted = (
            match_type in {"EXACT", "CONFIDENCE"}
            and rank == "SPECIES"
            and kingdom == "Plantae"
            and bool(family)
        )
        base_audit = {
            "input_name": input_name,
            "accepted": accepted,
            "match_type": match_type,
            "rank": rank,
            "kingdom": kingdom,
            "family": family,
            "accepted_name": accepted_name,
            "usage_key": match.get("acceptedUsageKey") or match.get("usageKey") or "",
            "input_spatial_scale": row.get("spatial_scale", ""),
            "source_id": stable_source_id(row),
        }
        if not accepted:
            base_audit["post_resolution_status"] = "rejected_taxon"
            audit.append(base_audit)
            continue

        resolved = dict(row)
        resolved["original_candidate_name"] = input_name
        resolved["canonical_name"] = accepted_name
        resolved["family"] = family
        resolved["gbif_taxon_validation"] = "strict_species_match"
        resolved_groups[accepted_name].append(resolved)
        base_audit["post_resolution_status"] = "pending_species_reaggregation"
        audit.append(base_audit)
        time.sleep(0.05)

    kept: list[dict[str, Any]] = []
    mixed_names: set[str] = set()
    for accepted_name, rows in sorted(resolved_groups.items()):
        scales = {row.get("spatial_scale", "") for row in rows if row.get("spatial_scale")}
        if scales != {"within_population"} and scales != {"among_population"}:
            mixed_names.add(accepted_name)
            continue

        # Deduplicate supporting sources before summing evidence effort.
        unique_sources: dict[str, dict[str, str]] = {}
        for index, row in enumerate(rows):
            source_key = stable_source_id(row) or f"row:{index}:{row.get('original_candidate_name', '')}"
            previous = unique_sources.get(source_key)
            if previous is None or integer(row.get("n_supporting_p1_records")) > integer(previous.get("n_supporting_p1_records")):
                unique_sources[source_key] = row
        evidence_rows = list(unique_sources.values())
        best = max(
            evidence_rows,
            key=lambda row: (
                integer(row.get("n_supporting_p1_records")),
                len(str(row.get("best_match_evidence") or "")),
            ),
        )
        output = dict(best)
        output["canonical_name"] = accepted_name
        output["family"] = rows[0]["family"]
        output["spatial_scale"] = next(iter(scales))
        output["classification_source"] = "exploratory_automated"
        output["automated_classification_source"] = "systematic_P1_contextual_binomial_post_GBIF_reaggregated"
        output["n_supporting_p1_records"] = sum(integer(row.get("n_supporting_p1_records")) for row in evidence_rows)
        output["n_distinct_supporting_sources"] = len(evidence_rows)
        output["n_input_names_resolved"] = len({row.get("original_candidate_name", "") for row in rows})
        output["resolved_input_names"] = ";".join(sorted({row.get("original_candidate_name", "") for row in rows}))
        output["post_gbif_scale_conflict"] = False
        kept.append(output)

    for row in audit:
        if row.get("accepted_name") in mixed_names:
            row["post_resolution_status"] = "mixed_after_synonym_resolution_excluded"
        elif row.get("post_resolution_status") == "pending_species_reaggregation":
            row["post_resolution_status"] = "accepted_unambiguous_after_reaggregation"

    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if kept:
        fields = sorted({key for row in kept for key in row})
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(kept)

    audit_path = Path(args.audit)
    with audit_path.open("w", newline="", encoding="utf-8") as handle:
        fields = sorted({key for row in audit for key in row})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(audit)

    summary = {
        "requested_input_names": len(input_rows),
        "accepted_species_before_scale_reaggregation": len(resolved_groups),
        "mixed_species_excluded_after_synonym_resolution": len(mixed_names),
        "strict_species_accepted": len(kept),
        "within": sum(row["spatial_scale"] == "within_population" for row in kept),
        "among": sum(row["spatial_scale"] == "among_population" for row in kept),
        "families": len({row["family"] for row in kept}),
        "classification_source_preserved": "exploratory_automated",
    }
    print(json.dumps(summary, indent=2))
    if len(kept) < 20 or len({row["spatial_scale"] for row in kept}) < 2:
        raise SystemExit(json.dumps(summary))


if __name__ == "__main__":
    main()
