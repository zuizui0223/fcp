#!/usr/bin/env python3
"""Finalize duplicate v2.2 record screening after explicit disagreement adjudication.

Reviewer consensus is copied automatically only when both independent values agree.
Disagreements require an explicit adjudicated value. This script does not infer taxa or
spatial states.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

CORE = [
    "record_relevance",
    "natural_intraspecific_variation",
    "floral_display_colour",
    "full_text_required",
    "focal_taxon_text",
]


def clean(value):
    return " ".join(str(value or "").split()).strip()


def norm(value):
    return clean(value).lower()


def read(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{f: r.get(f, "") for f in fields} for r in rows])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reviewed", required=True)
    p.add_argument("--adjudication", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--retained-out", required=True)
    p.add_argument("--fulltext-out", required=True)
    p.add_argument("--excluded-out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    reviewed = read(args.reviewed)
    adjudication = read(args.adjudication)
    amap = {r.get("record_review_id", ""): r for r in adjudication if r.get("record_review_id")}
    if len(amap) != len(adjudication):
        raise SystemExit("Duplicate or missing record_review_id in adjudication sheet")

    finalized = []
    unresolved = []
    for row in reviewed:
        rid = row.get("record_review_id", "")
        out = dict(row)
        adj = amap.get(rid, {})
        unresolved_fields = []
        for field in CORE:
            a = clean(row.get(f"reviewer_1_{field}"))
            b = clean(row.get(f"reviewer_2_{field}"))
            if a and b and norm(a) == norm(b):
                value = a
            elif a and b:
                value = clean(adj.get(f"adjudicated_{field}"))
                if not value:
                    unresolved_fields.append(field)
            else:
                value = ""
                unresolved_fields.append(field)
            out[f"adjudicated_{field}"] = value

        # Exclusion reason is required only for an adjudicated exclusion.
        r1_reason = clean(row.get("reviewer_1_exclusion_reason"))
        r2_reason = clean(row.get("reviewer_2_exclusion_reason"))
        if norm(out.get("adjudicated_record_relevance")) == "exclude":
            if r1_reason and r2_reason and norm(r1_reason) == norm(r2_reason):
                reason = r1_reason
            else:
                reason = clean(adj.get("adjudicated_exclusion_reason"))
            if not reason:
                unresolved_fields.append("exclusion_reason")
            out["adjudicated_exclusion_reason"] = reason
        else:
            out["adjudicated_exclusion_reason"] = ""

        out["adjudication_notes"] = clean(adj.get("adjudication_notes"))
        out["review_status"] = "adjudicated" if not unresolved_fields else "adjudication_incomplete"
        out["unresolved_adjudication_fields"] = ";".join(sorted(set(unresolved_fields)))
        finalized.append(out)
        if unresolved_fields:
            unresolved.append(out)

    base_fields = list(finalized[0].keys()) if finalized else []
    write(args.out, finalized, base_fields)

    retained = [r for r in finalized if norm(r.get("adjudicated_record_relevance")) == "include"]
    fulltext = [
        r for r in finalized
        if norm(r.get("adjudicated_full_text_required")) == "yes"
        or norm(r.get("adjudicated_record_relevance")) == "uncertain"
        or norm(r.get("adjudicated_natural_intraspecific_variation")) == "uncertain"
        or norm(r.get("adjudicated_floral_display_colour")) == "uncertain"
    ]
    excluded = [r for r in finalized if norm(r.get("adjudicated_record_relevance")) == "exclude"]
    write(args.retained_out, retained, base_fields)
    write(args.fulltext_out, fulltext, base_fields)
    write(args.excluded_out, excluded, base_fields)

    counts = {
        "include": len(retained),
        "exclude": len(excluded),
        "uncertain": sum(norm(r.get("adjudicated_record_relevance")) == "uncertain" for r in finalized),
        "blank_or_unresolved_relevance": sum(not norm(r.get("adjudicated_record_relevance")) for r in finalized),
    }
    summary = {
        "status": "complete" if not unresolved else "not_ready",
        "records": len(finalized),
        "adjudication_incomplete_records": len(unresolved),
        "record_relevance_counts": counts,
        "retained_records": len(retained),
        "full_text_queue_records": len(fulltext),
        "excluded_records": len(excluded),
        "ready_for_taxon_stage": not unresolved,
        "semantic_guard": (
            "Record screening only. Retention does not establish natural eligibility, accepted taxonomy, "
            "local coexistence, geographic structure, or mixed state."
        ),
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
