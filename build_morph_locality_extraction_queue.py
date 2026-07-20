#!/usr/bin/env python3
"""Build a reproducible manual extraction queue for morph-labelled locality evidence.

The queue combines each species' retained baseline evidence with eligible enrichment works.
It never infers morph localities from generic GBIF records.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

OUTPUT_COLUMNS = [
    "priority_rank", "canonical_name", "family", "spatial_scale",
    "classification_source", "evidence_origin", "openalex_id", "title", "year",
    "doi", "landing_url", "evidence_score", "evidence_snippet",
    "extraction_status", "country_or_region", "locality_name", "latitude",
    "longitude", "coordinate_source", "population_label", "colour_state",
    "colour_state_verbatim", "sample_size", "locality_evidence_quote",
    "page_or_table", "extractor_notes",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True)
    parser.add_argument("--works", required=True)
    parser.add_argument("--review", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    classification = pd.read_csv(args.classification)
    works = pd.read_csv(args.works)
    review = pd.read_csv(args.review)

    required_class = {"canonical_name", "family", "enriched_scale", "classification_source"}
    required_works = {
        "canonical_name", "openalex_id", "title", "year", "doi", "landing_url",
        "score", "evidence_snippet", "direct_colour_signal", "artificial_signal",
    }
    required_review = {
        "canonical_name", "family", "best_title", "best_doi", "best_openalex_id",
        "best_match_evidence", "max_score",
    }
    missing = sorted(
        (required_class - set(classification.columns))
        | (required_works - set(works.columns))
        | (required_review - set(review.columns))
    )
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    classified = classification[
        classification["enriched_scale"].isin(["within_population", "among_population"])
    ][["canonical_name", "family", "enriched_scale", "classification_source"]].drop_duplicates("canonical_name")

    eligible = works.copy()
    for column in ("score", "direct_colour_signal", "artificial_signal"):
        eligible[column] = pd.to_numeric(eligible[column], errors="coerce").fillna(0)
    eligible = eligible[
        (eligible["direct_colour_signal"] == 1) & (eligible["artificial_signal"] == 0)
    ].copy()
    eligible = eligible.rename(columns={"score": "evidence_score"})
    eligible["evidence_origin"] = "openalex_enrichment"

    baseline = review[list(required_review)].drop_duplicates("canonical_name").copy()
    baseline = baseline.rename(columns={
        "best_title": "title",
        "best_doi": "doi",
        "best_openalex_id": "openalex_id",
        "best_match_evidence": "evidence_snippet",
        "max_score": "evidence_score",
    })
    baseline["year"] = ""
    baseline["landing_url"] = baseline["openalex_id"].fillna("")
    baseline["evidence_origin"] = "baseline_best_evidence"

    evidence_columns = [
        "canonical_name", "openalex_id", "title", "year", "doi", "landing_url",
        "evidence_score", "evidence_snippet", "evidence_origin",
    ]
    evidence = pd.concat(
        [baseline[evidence_columns], eligible[evidence_columns]],
        ignore_index=True,
    )
    evidence["openalex_id"] = evidence["openalex_id"].fillna("").astype(str)
    evidence["doi"] = evidence["doi"].fillna("").astype(str)
    evidence["title"] = evidence["title"].fillna("").astype(str)
    evidence["evidence_snippet"] = evidence["evidence_snippet"].fillna("").astype(str)
    evidence["evidence_score"] = pd.to_numeric(evidence["evidence_score"], errors="coerce").fillna(0)
    evidence["dedupe_key"] = evidence["openalex_id"].where(
        evidence["openalex_id"].str.len() > 0,
        evidence["doi"].where(evidence["doi"].str.len() > 0, evidence["title"].str.lower()),
    )
    evidence = evidence.sort_values(
        ["canonical_name", "evidence_score", "evidence_origin"],
        ascending=[True, False, True],
        kind="stable",
    ).drop_duplicates(["canonical_name", "dedupe_key"])

    queue = classified.merge(evidence.drop(columns="dedupe_key"), on="canonical_name", how="left")
    queue["spatial_scale"] = queue["enriched_scale"]
    queue["scale_priority"] = queue["spatial_scale"].map(
        {"among_population": 0, "within_population": 1}
    ).fillna(2)
    queue["source_priority"] = queue["classification_source"].map(
        {"baseline_unambiguous": 0, "high_confidence_enrichment": 1}
    ).fillna(2)
    queue["origin_priority"] = queue["evidence_origin"].map(
        {"baseline_best_evidence": 0, "openalex_enrichment": 1}
    ).fillna(2)
    queue = queue.sort_values(
        ["scale_priority", "source_priority", "origin_priority", "evidence_score", "canonical_name"],
        ascending=[True, True, True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    queue["priority_rank"] = range(1, len(queue) + 1)

    defaults = {
        "extraction_status": "not_started", "country_or_region": "", "locality_name": "",
        "latitude": "", "longitude": "", "coordinate_source": "",
        "population_label": "", "colour_state": "", "colour_state_verbatim": "",
        "sample_size": "", "locality_evidence_quote": "", "page_or_table": "",
        "extractor_notes": "",
    }
    for column, value in defaults.items():
        queue[column] = value
    for column in OUTPUT_COLUMNS:
        if column not in queue:
            queue[column] = ""

    out = queue[OUTPUT_COLUMNS]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    rows_without = int(
        (
            out["openalex_id"].fillna("").astype(str).str.len().eq(0)
            & out["doi"].fillna("").astype(str).str.len().eq(0)
            & out["title"].fillna("").astype(str).str.len().eq(0)
        ).sum()
    )
    manifest = {
        "classified_species": int(classified["canonical_name"].nunique()),
        "queue_rows": int(len(out)),
        "species_with_queue_rows": int(out["canonical_name"].nunique()),
        "among_population_species": int((classified["enriched_scale"] == "among_population").sum()),
        "within_population_species": int((classified["enriched_scale"] == "within_population").sum()),
        "among_population_queue_rows": int((out["spatial_scale"] == "among_population").sum()),
        "baseline_evidence_rows": int((out["evidence_origin"] == "baseline_best_evidence").sum()),
        "enrichment_evidence_rows": int((out["evidence_origin"] == "openalex_enrichment").sum()),
        "rows_without_retained_work": rows_without,
        "required_manual_fields": [
            "locality_name", "latitude", "longitude", "population_label",
            "colour_state", "locality_evidence_quote",
        ],
        "semantic_guard": (
            "The queue requests explicit morph-labelled locality evidence from literature; generic "
            "GBIF occurrences must not be assigned colour states by proximity or clustering."
        ),
    }
    Path(args.manifest).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
