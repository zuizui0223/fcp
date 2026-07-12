#!/usr/bin/env python3
"""Merge parallel mass-screen CSV chunks and produce a species ranking.

The workflow writes each chunk independently, then this script combines them in a
single aggregation job. Rows are deduplicated conservatively by the complete row
content so repeated API hits do not inflate scores.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


def read_rows(paths: Sequence[Path]) -> Tuple[List[str], List[Dict[str, str]]]:
    fieldnames: List[str] = []
    rows: List[Dict[str, str]] = []
    seen = set()
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames and not fieldnames:
                fieldnames = list(reader.fieldnames)
            for row in reader:
                key = tuple((name, row.get(name, "")) for name in (reader.fieldnames or []))
                if key in seen:
                    continue
                seen.add(key)
                rows.append(dict(row))
    return fieldnames, rows


def write_rows(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="artifacts")
    parser.add_argument("--glob", default="**/mass_screen_hits_batch_*.csv")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    paths = sorted(Path(args.input_dir).glob(args.glob))
    if not paths:
        raise SystemExit(f"No chunk files matched {args.input_dir}/{args.glob}")
    fieldnames, rows = read_rows(paths)
    if not fieldnames:
        raise SystemExit("Chunk files had no CSV header")
    write_rows(Path(args.out), fieldnames, rows)
    print(f"merged {len(paths)} chunks and {len(rows)} unique rows into {args.out}")


if __name__ == "__main__":
    main()
