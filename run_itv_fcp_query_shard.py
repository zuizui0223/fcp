#!/usr/bin/env python3
"""Run one query/database shard of the ITV/FCP systematic search.

The shard writes a status file before retrieval so GitHub Actions can preserve
provenance even when an external service fails or the job is cancelled.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--census", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--query-id", required=True)
    parser.add_argument("--database", required=True, choices=["openalex", "crossref", "europe_pmc"])
    parser.add_argument("--max-records", type=int, required=True)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    status_path = outdir / "shard_status.json"
    status = {
        "status": "running",
        "shard_type": "query_database",
        "query_id": args.query_id,
        "database": args.database,
        "max_records": args.max_records,
        "started_at_utc": now(),
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    selected = [query for query in config.get("queries", []) if str(query.get("id")) == args.query_id]
    if len(selected) != 1:
        status.update({"status": "failed", "reason": f"expected one query, found {len(selected)}", "completed_at_utc": now()})
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise SystemExit(2)

    shard_config = dict(config)
    shard_config["queries"] = selected
    shard_config["seed_reviews"] = []
    shard_config_path = outdir / "shard_config.json"
    shard_config_path.write_text(json.dumps(shard_config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    command = [
        sys.executable,
        "run_systematic_itv_fcp_search.py",
        "--config", str(shard_config_path),
        "--census", args.census,
        "--outdir", str(outdir),
        "--max-records-per-query", str(args.max_records),
        "--max-citations-per-seed", "1",
        "--databases", args.database,
        "--timeout", str(args.timeout),
        "--retries", str(args.retries),
        "--skip-citation-chasing",
    ]
    completed = subprocess.run(command, env=os.environ.copy(), check=False)
    status.update({
        "status": "success" if completed.returncode == 0 else "failed",
        "return_code": completed.returncode,
        "completed_at_utc": now(),
    })
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
