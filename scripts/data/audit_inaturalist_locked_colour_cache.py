#!/usr/bin/env python3
"""Audit the 717-photo locked automated-colour cache."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

if __package__:
    from . import audit_inaturalist_automated_colour_cache as shared
else:
    import audit_inaturalist_automated_colour_cache as shared  # type: ignore[no-redef]


EXPECTED_LOCKED_PHOTOS = 717


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--locked-packet", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    shared.EXPECTED_TOTAL_PHOTOS = EXPECTED_LOCKED_PHOTOS
    report = shared.audit_cache(
        args.locked_packet.resolve(), args.cache_dir.resolve(), args.model_dir.resolve()
    )
    report["partition"] = "locked_60_development_passing_species_only"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] == "invalid_cache":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
