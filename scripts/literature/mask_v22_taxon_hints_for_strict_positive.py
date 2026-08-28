#!/usr/bin/env python3
"""Mask automatic taxon hints outside strict spatial-positive records.

This is a performance-only preprocessing step for the strict single-pass builder.
All 12,064 records are retained. The coordinator key is copied unchanged except that
`detected_binomial_strings` is blanked when the record does not satisfy an explicit
local or geographic positive-evidence pattern. Historical source->taxon rescue is
handled separately by the downstream builder.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from build_v22_single_pass_strict_evidence import (
    DISPLAY_RE,
    DISCRETE_RE,
    GEO_PATTERNS,
    LOCAL_PATTERNS,
    clean,
    match_snippet,
)


def read(path: Path):
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--blind", required=True)
    p.add_argument("--key", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    blind = read(Path(args.blind))
    key = read(Path(args.key))
    bmap = {r["record_review_id"]: r for r in blind}
    if len(blind) != 12064 or len(key) != 12064 or set(bmap) != {r["record_review_id"] for r in key}:
        raise SystemExit("blind/key boundary mismatch")

    kept = 0
    masked = 0
    for row in key:
        b = bmap[row["record_review_id"]]
        text = clean(f"{b.get('title','')} {b.get('abstract','')}")
        display = bool(DISPLAY_RE.search(text))
        discrete = bool(DISCRETE_RE.search(text))
        local = bool(match_snippet(LOCAL_PATTERNS, text)) if display and discrete else False
        geo = bool(match_snippet(GEO_PATTERNS, text)) if display else False
        if local or geo:
            kept += 1
        else:
            row["detected_binomial_strings"] = ""
            masked += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(key[0].keys()))
        w.writeheader()
        w.writerows(key)
    print({"records": len(key), "strict_positive_taxon_hints_kept": kept, "other_taxon_hints_masked": masked})


if __name__ == "__main__":
    main()
