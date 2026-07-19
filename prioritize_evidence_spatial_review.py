#!/usr/bin/env python3
"""Prioritize manual review of ambiguous flower-colour evidence scale.

This script never converts unclear evidence into a biological class. It only ranks rows
for manual review using evidence richness, conflicting signals, review priority, and the
potential to change the within- versus among-population comparison.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PRIORITY_WEIGHT = {"P0": 4, "P1": 3, "P2": 2, "P3": 1}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screening", required=True)
    ap.add_argument("--outdir", default="analysis_outputs/evidence_scale_review")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(args.screening)
    q = d[d["manual_review_required"].astype(str).str.lower().isin(["true", "1"])].copy()

    for col in ["best_title", "best_match_evidence"]:
        q[col] = q[col].fillna("").astype(str)
    q["evidence_chars"] = q["best_match_evidence"].str.len()
    q["priority_weight"] = q["review_priority"].map(PRIORITY_WEIGHT).fillna(0)
    q["conflict_weight"] = (q["evidence_spatial_scale"] == "mixed").astype(int) * 4
    q["maintenance_weight"] = pd.to_numeric(q["maintenance_mechanism_signal"], errors="coerce").fillna(0) * 2
    q["doi_weight"] = q["best_doi"].fillna("").astype(str).str.len().gt(0).astype(int)
    q["evidence_weight"] = (q["evidence_chars"].clip(upper=800) / 200).round(3)
    q["manual_review_score"] = (
        q["priority_weight"] + q["conflict_weight"] + q["maintenance_weight"]
        + q["doi_weight"] + q["evidence_weight"]
    )
    q["review_question"] = q["evidence_spatial_scale"].map({
        "mixed": "Does the primary study document morph coexistence within populations, geographic differentiation among populations, or both?",
        "unclear": "Does the primary study explicitly show multiple colour morphs coexisting within at least one natural population?",
    })
    q["manual_label"] = ""
    q["manual_notes"] = ""
    q = q.sort_values(["manual_review_score", "review_priority", "canonical_name"], ascending=[False, True, True])
    q.insert(0, "manual_review_rank", range(1, len(q) + 1))

    cols = [
        "manual_review_rank", "manual_review_score", "review_priority", "canonical_name", "family",
        "evidence_spatial_scale", "within_population_signal", "geographic_differentiation_signal",
        "maintenance_mechanism_signal", "review_question", "best_title", "best_doi",
        "best_openalex_id", "best_match_evidence", "manual_label", "manual_notes",
    ]
    q[cols].to_csv(out / "evidence_spatial_scale_prioritized_review.csv", index=False)
    q.head(25)[cols].to_csv(out / "evidence_spatial_scale_top25_review.csv", index=False)

    summary = {
        "manual_review_rows": int(len(q)),
        "mixed_rows": int((q["evidence_spatial_scale"] == "mixed").sum()),
        "unclear_rows": int((q["evidence_spatial_scale"] == "unclear").sum()),
        "top25_with_doi": int(q.head(25)["best_doi"].fillna("").astype(str).str.len().gt(0).sum()),
        "top25_maintenance_signal": int(pd.to_numeric(q.head(25)["maintenance_mechanism_signal"], errors="coerce").fillna(0).sum()),
        "semantic_guard": "ranking only; no ambiguous row is reclassified automatically",
        "recommended_labels": ["within_population", "among_population", "mixed", "not_natural_polymorphism", "insufficient_evidence"],
    }
    (out / "evidence_spatial_scale_review_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
