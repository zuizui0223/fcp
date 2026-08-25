#!/usr/bin/env python3
"""Build a concurrent GBIF resolution cache for upstream spatial re-audit.

This helper only validates candidate taxon names. It does not classify biological
eligibility or spatial organization. Candidate names are extracted from the frozen
systematic-screening queue using the exact contextual-name function used by the
mixed-preserving review builder, then resolved concurrently against GBIF.
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


def collect_names(queue: Path, priorities: set[str]) -> tuple[set[str], dict[str, int]]:
    names: set[str] = set()
    diagnostics = {
        "records_total": 0,
        "records_selected": 0,
        "records_with_contextual_candidate": 0,
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
                names.update(contextual)
    diagnostics["unique_contextual_input_names"] = len(names)
    return names, diagnostics


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
    names, diagnostics = collect_names(Path(args.queue), priorities)
    if not names:
        raise SystemExit("No contextual candidate names found")

    cache: dict[str, dict[str, Any]] = {}
    workers = max(1, min(args.workers, 48))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(resolve_one, name, args.timeout, args.retries): name
            for name in sorted(names)
        }
        for index, future in enumerate(as_completed(futures), start=1):
            name, resolved = future.result()
            cache[name] = resolved
            if index % 1000 == 0:
                print({"resolved": index, "total": len(names)}, flush=True)

    accepted = [row for row in cache.values() if bool(row.get("accepted"))]
    errors = [row for row in cache.values() if str(row.get("reason", "")).startswith("resolution_error:")]
    qc = {
        "status": "complete",
        **diagnostics,
        "workers": workers,
        "gbif_accepted_input_names": len(accepted),
        "gbif_rejected_input_names": len(cache) - len(accepted),
        "gbif_resolution_errors": len(errors),
        "accepted_species_after_synonym_collapse": len({str(row.get("accepted_name")) for row in accepted if row.get("accepted_name")}),
        "semantic_guard": "Taxon validation only; no biological eligibility or spatial classification is inferred here.",
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2), flush=True)

    # A few network failures are tolerable because the re-audit is a review queue,
    # but a broad outage should fail loudly rather than silently shrink the sample.
    if len(errors) > max(25, int(0.02 * len(cache))):
        raise SystemExit(f"Too many GBIF resolution errors: {len(errors)}/{len(cache)}")
    if len(accepted) < 50:
        raise SystemExit(f"Unexpectedly few accepted plant species: {len(accepted)}")


if __name__ == "__main__":
    main()
