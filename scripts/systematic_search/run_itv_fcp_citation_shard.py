#!/usr/bin/env python3
"""Run citation chasing for one prespecified seed review."""
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
    parser.add_argument("--seed-index", type=int, required=True)
    parser.add_argument("--max-citations", type=int, default=2000)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    status_path = outdir / "shard_status.json"
    status = {
        "status": "running",
        "shard_type": "citation_seed",
        "seed_index": args.seed_index,
        "max_citations": args.max_citations,
        "started_at_utc": now(),
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    seeds = list(config.get("seed_reviews", []))
    if args.seed_index < 0 or args.seed_index >= len(seeds):
        status.update({"status": "failed", "reason": "seed index out of range", "completed_at_utc": now()})
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise SystemExit(2)

    seed = seeds[args.seed_index]
    status["seed_doi"] = seed.get("doi", "")
    shard_config = dict(config)
    shard_config["queries"] = []
    shard_config["seed_reviews"] = [seed]
    shard_config_path = outdir / "shard_config.json"
    shard_config_path.write_text(json.dumps(shard_config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

    command = [
        sys.executable,
        "run_systematic_itv_fcp_search.py",
        "--config", str(shard_config_path),
        "--census", args.census,
        "--outdir", str(outdir),
        "--max-records-per-query", "1",
        "--max-citations-per-seed", str(args.max_citations),
        "--databases", "openalex",
        "--timeout", str(args.timeout),
        "--retries", str(args.retries),
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
