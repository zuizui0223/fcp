#!/usr/bin/env python3
"""Build a deterministic blinded calibration wave from the canonical v2.2 record screen."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def clean(value):
    return " ".join(str(value or "").split()).strip()


def rank(seed: str, record_id: str) -> str:
    return hashlib.sha256(f"{seed}|{record_id}".encode()).hexdigest()


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blind", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--blind-out", required=True)
    parser.add_argument("--key-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--seed", default="fcp-jbi-v22-calibration-wave0-v1")
    args = parser.parse_args()

    blind = read_csv(args.blind)
    key = read_csv(args.key)
    bmap = {row["record_review_id"]: row for row in blind}
    kmap = {row["record_review_id"]: row for row in key}
    if set(bmap) != set(kmap):
        raise SystemExit("Blind sheet and coordinator key are not aligned")

    selected: list[str] = []
    strata: dict[str, str] = {}

    def take(name: str, pool: list[str], n: int) -> None:
        available = [rid for rid in pool if rid not in strata]
        available = sorted(available, key=lambda rid: rank(args.seed + name, rid))
        chosen = available[:n]
        if len(chosen) < n:
            raise SystemExit(f"insufficient {name}: {len(chosen)}/{n}")
        for rid in chosen:
            strata[rid] = name
        selected.extend(chosen)

    all_ids = list(bmap)
    historical = [
        rid for rid in all_ids
        if str(kmap[rid].get("historical_34_exact_source_match", "0")) == "1"
    ]
    if len(historical) != 34:
        raise SystemExit(f"historical benchmark stratum: {len(historical)}/34")
    for rid in sorted(historical, key=lambda rid: rank(args.seed + "hist", rid)):
        strata[rid] = "historical_benchmark_blinded"
        selected.append(rid)

    no_binomial = [
        rid for rid in all_ids
        if not clean(kmap[rid].get("detected_binomial_strings"))
    ]
    take("no_detected_binomial", no_binomial, 100)

    with_binomial = [
        rid for rid in all_ids
        if clean(kmap[rid].get("detected_binomial_strings"))
    ]
    take("detected_binomial", with_binomial, 100)

    non_english = [
        rid for rid in all_ids
        if clean(bmap[rid].get("language")).lower() not in {"", "en"}
    ]
    take("non_english", non_english, 100)

    missing_abstract = [
        rid for rid in all_ids
        if not clean(bmap[rid].get("abstract"))
    ]
    take("missing_abstract", missing_abstract, 50)

    ordered = sorted(selected, key=lambda rid: rank(args.seed + "review-order", rid))
    blind_out = []
    key_out = []
    for order, rid in enumerate(ordered, 1):
        row = dict(bmap[rid])
        row["calibration_order"] = order
        row["calibration_wave"] = "W0"
        blind_out.append(row)

        hidden = dict(kmap[rid])
        hidden["calibration_stratum"] = strata[rid]
        hidden["calibration_order"] = order
        key_out.append(hidden)

    blind_fields = ["calibration_wave", "calibration_order"] + list(blind[0].keys())
    key_fields = list(key[0].keys()) + ["calibration_stratum", "calibration_order"]
    write_csv(args.blind_out, blind_out, blind_fields)
    write_csv(args.key_out, key_out, key_fields)

    summary = {
        "status": "complete",
        "calibration_records": len(blind_out),
        "strata": dict(Counter(strata.values())),
        "historical_status_visible_to_reviewer": False,
        "query_membership_visible_to_reviewer": False,
        "automated_taxon_hint_visible_to_reviewer": False,
        "purpose": (
            "Calibrate duplicate blind record-screening instructions and inter-rater agreement "
            "before full 12,064-record screening."
        ),
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))

    if len(blind_out) != 384:
        raise SystemExit(f"Expected 384 calibration records, found {len(blind_out)}")


if __name__ == "__main__":
    main()
