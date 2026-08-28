#!/usr/bin/env python3
"""Build a bounded GBIF resolution cache for the upstream spatial re-audit.

The systematic-search artifact contains permissive two-word candidate strings. This
helper therefore separates *source attribution* from taxon validation before any
biological classification is attempted:

1. collect every contextual two-word candidate seen by the mixed-preserving builder;
2. designate a primary source candidate only when the name occurs in the title or in
   the first 2,000 abstract characters (the source-attribution window);
3. query GBIF concurrently only for those primary candidates;
4. write explicit non-primary rejection rows for the remaining candidate strings so
   downstream code never silently falls back to thousands of sequential GBIF calls.

No natural-status or spatial-state decision is made here.
"""
from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from build_systematic_spatial_evidence_axes import (
    DISCOVERY_PRIORITIES,
    contextual_candidate_names,
    resolve_gbif_name,
)


def clean(value: Any) -> str:
    return " ".join(str(value or "").split())


def primary_source_candidates(row: dict[str, str], abstract_chars: int = 2000) -> list[str]:
    """Return candidates plausibly attributable to the focal publication.

    Title occurrence is preferred. If none of the candidates is in the title, retain
    names occurring near the beginning of the abstract. This is only a record-to-taxon
    attribution rule; full-source review later decides biological relevance and spatial
    evidence.
    """
    candidates = contextual_candidate_names(row)
    if not candidates:
        return []
    title = clean(row.get("title")).lower()
    title_hits = [name for name in candidates if name.lower() in title]
    if title_hits:
        return title_hits[:5]
    abstract = clean(row.get("abstract"))[:abstract_chars].lower()
    return [name for name in candidates if name.lower() in abstract][:5]


def collect_names(
    queue: Path,
    priorities: set[str],
) -> tuple[set[str], set[str], dict[str, int]]:
    all_contextual: set[str] = set()
    primary: set[str] = set()
    diagnostics = {
        "records_total": 0,
        "records_selected": 0,
        "records_with_contextual_candidate": 0,
        "records_with_primary_source_candidate": 0,
    }
    with queue.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            diagnostics["records_total"] += 1
            if str(row.get("screen_priority") or "") not in priorities:
                continue
            diagnostics["records_selected"] += 1
            contextual = contextual_candidate_names(row)
            if contextual:
                diagnostics["records_with_contextual_candidate"] += 1
                all_contextual.update(contextual)
            focal = primary_source_candidates(row)
            if focal:
                diagnostics["records_with_primary_source_candidate"] += 1
                primary.update(focal)
    diagnostics["unique_contextual_input_names"] = len(all_contextual)
    diagnostics["unique_primary_source_candidates"] = len(primary)
    diagnostics["contextual_names_deferred_as_nonprimary"] = len(all_contextual - primary)
    return all_contextual, primary, diagnostics


def resolve_one(name: str, timeout: int, retries: int) -> tuple[str, dict[str, Any]]:
    try:
        resolved = resolve_gbif_name(name, timeout, retries)
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
    return name, resolved


def nonprimary_row(name: str) -> dict[str, Any]:
    return {
        "input_name": name,
        "accepted": False,
        "match_type": "",
        "rank": "",
        "kingdom": "",
        "confidence": 0,
        "accepted_name": "",
        "family": "",
        "usage_key": "",
        "reason": "deferred_nonprimary_source_candidate",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--priorities", default=",".join(sorted(DISCOVERY_PRIORITIES)))
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--timeout", type=int, default=25)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    priorities = {x.strip() for x in args.priorities.split(",") if x.strip()}
    all_names, primary_names, diagnostics = collect_names(Path(args.queue), priorities)
    if not primary_names:
        raise SystemExit("No primary source candidate names found")

    # Populate the complete downstream cache up front. Names outside the focal
    # source-attribution window remain visible as deferred, not silently discarded.
    cache: dict[str, dict[str, Any]] = {
        name: nonprimary_row(name) for name in sorted(all_names - primary_names)
    }

    workers = max(1, min(args.workers, 48))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(resolve_one, name, args.timeout, args.retries): name
            for name in sorted(primary_names)
        }
        for index, future in enumerate(as_completed(futures), start=1):
            name, resolved = future.result()
            cache[name] = resolved
            if index % 250 == 0:
                print({"gbif_resolved": index, "gbif_total": len(primary_names)}, flush=True)

    primary_rows = [cache[name] for name in primary_names]
    accepted = [row for row in primary_rows if bool(row.get("accepted"))]
    errors = [row for row in primary_rows if str(row.get("reason", "")).startswith("resolution_error:")]
    qc = {
        "status": "complete",
        **diagnostics,
        "workers": workers,
        "gbif_queries_sent": len(primary_names),
        "gbif_accepted_input_names": len(accepted),
        "gbif_rejected_primary_input_names": len(primary_names) - len(accepted),
        "gbif_resolution_errors": len(errors),
        "accepted_species_after_synonym_collapse": len({str(row.get("accepted_name")) for row in accepted if row.get("accepted_name")}),
        "semantic_guard": "Source attribution and taxon validation only; no biological eligibility or spatial classification is inferred here.",
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2), flush=True)

    if len(errors) > max(10, int(0.03 * len(primary_names))):
        raise SystemExit(f"Too many GBIF resolution errors: {len(errors)}/{len(primary_names)}")
    if len(accepted) < 40:
        raise SystemExit(f"Unexpectedly few accepted plant species: {len(accepted)}")


if __name__ == "__main__":
    main()
