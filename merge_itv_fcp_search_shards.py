#!/usr/bin/env python3
"""Merge independently retrieved ITV/FCP query and citation shards."""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import systematic_itv_fcp_search as search


def read_csv(path: Path) -> list[dict[str, Any]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--census", required=True)
    parser.add_argument("--shards-root", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    config_path = Path(args.config)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    shards_root = Path(args.shards_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    statuses: list[dict[str, Any]] = []
    for path in sorted(shards_root.rglob("shard_status.json")):
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:  # noqa: BLE001
            row = {"status": "invalid", "reason": f"{type(exc).__name__}: {exc}"}
        row["status_file"] = str(path.relative_to(shards_root))
        statuses.append(row)

    raw_records: list[dict[str, Any]] = []
    search_logs: list[dict[str, Any]] = []
    for path in sorted(shards_root.rglob("itv_fcp_search_raw_records.csv")):
        raw_records.extend(read_csv(path))
    for path in sorted(shards_root.rglob("itv_fcp_search_log.csv")):
        search_logs.extend(read_csv(path))

    raw_path = outdir / "itv_fcp_search_raw_records.csv"
    log_path = outdir / "itv_fcp_search_log.csv"
    search.write_csv(raw_path, raw_records)
    search.write_csv(log_path, search_logs)

    deduped = search.deduplicate(raw_records)
    census = search.load_census(Path(args.census))
    terms = config.get("screening_terms")
    if not isinstance(terms, dict):
        raise ValueError("screening_terms missing from config")
    screened = [search.screen_record(row, terms, census) for row in deduped]
    priority_order = {
        "P1_high_natural_itv": 0,
        "P1_high_population_itv": 1,
        "P2_possible_itv": 2,
        "P3_likely_exclusion_review": 3,
        "P4_flower_colour_context_only": 4,
        "P5_low_relevance": 5,
    }
    screened.sort(
        key=lambda x: (
            priority_order.get(str(x.get("screen_priority")), 99),
            -int(x.get("candidate_species_count") or 0),
            -int(x.get("cited_by_count") or 0),
            search.normalize_title(x.get("title")),
        )
    )

    dedup_path = outdir / "itv_fcp_search_deduplicated_records.csv"
    queue_path = outdir / "itv_fcp_human_screening_queue.csv"
    status_path = outdir / "itv_fcp_search_shard_status.csv"
    search.write_csv(dedup_path, deduped)
    search.write_csv(queue_path, screened)
    search.write_csv(status_path, statuses)

    expected_query_keys = {
        (str(query["id"]), database)
        for query in config.get("queries", [])
        for database in ("openalex", "crossref", "europe_pmc")
    }
    found_query_keys = {
        (str(row.get("query_id")), str(row.get("database")))
        for row in statuses
        if row.get("shard_type") == "query_database"
    }
    expected_seed_indexes = set(range(len(config.get("seed_reviews", []))))
    found_seed_indexes = {
        int(row["seed_index"])
        for row in statuses
        if row.get("shard_type") == "citation_seed" and str(row.get("seed_index", "")).isdigit()
    }
    failed_statuses = [row for row in statuses if row.get("status") != "success"]
    missing_queries = sorted(expected_query_keys - found_query_keys)
    missing_seeds = sorted(expected_seed_indexes - found_seed_indexes)

    priority_counts: dict[str, int] = defaultdict(int)
    for row in screened:
        priority_counts[str(row.get("screen_priority"))] += 1

    truncated = []
    for log in search_logs:
        value = str(log.get("truncated", "")).strip().lower()
        if value in {"true", "1", "yes"}:
            truncated.append(log)

    complete = not failed_statuses and not missing_queries and not missing_seeds
    manifest = {
        "status": "complete" if complete else "partial",
        "protocol_version": config.get("protocol_version"),
        "completed_at_utc": search.utc_now(),
        "databases": ["crossref", "europe_pmc", "openalex"],
        "query_count": len(config.get("queries", [])),
        "seed_review_count": len(config.get("seed_reviews", [])),
        "expected_query_database_shards": len(expected_query_keys),
        "expected_citation_shards": len(expected_seed_indexes),
        "found_shards": len(statuses),
        "successful_shards": sum(row.get("status") == "success" for row in statuses),
        "failed_or_incomplete_shards": len(failed_statuses),
        "missing_query_database_shards": [f"{query}:{database}" for query, database in missing_queries],
        "missing_citation_seed_indexes": missing_seeds,
        "raw_records": len(raw_records),
        "deduplicated_records": len(deduped),
        "duplicate_records_removed": len(raw_records) - len(deduped),
        "screen_priority_counts": dict(sorted(priority_counts.items())),
        "records_with_candidate_binomials": sum(int(row.get("candidate_species_count") or 0) > 0 for row in screened),
        "queries_truncated": len(truncated),
        "truncated_query_ids": sorted({f"{row.get('database')}:{row.get('query_id')}" for row in truncated}),
        "human_review_required": True,
        "semantic_guard": (
            "Automated retrieval and screening identify candidate publications only. "
            "Natural ITV, strict FCP, spatial scale and species eligibility require duplicate source-level human review and adjudication."
        ),
        "files": {
            raw_path.name: {"sha256": search.sha256(raw_path), "rows": len(raw_records)},
            log_path.name: {"sha256": search.sha256(log_path), "rows": len(search_logs)},
            dedup_path.name: {"sha256": search.sha256(dedup_path), "rows": len(deduped)},
            queue_path.name: {"sha256": search.sha256(queue_path), "rows": len(screened)},
            status_path.name: {"sha256": search.sha256(status_path), "rows": len(statuses)},
            config_path.name: {"sha256": search.sha256(config_path)},
        },
    }
    (outdir / "itv_fcp_search_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
