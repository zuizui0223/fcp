#!/usr/bin/env python3
"""Build a reproducible manual extraction queue for morph-labelled locality evidence.

This step does not infer morph localities from generic GBIF records. It ranks literature records
for manual extraction of explicit colour-state, population, and locality information needed to
test environmental sorting directly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

OUTPUT_COLUMNS = [
    "priority_rank",
    "canonical_name",
    "family",
    "spatial_scale",
    "classification_source",
    "openalex_id",
    "title",
    "year",
    "doi",
    "landing_url",
    "evidence_score",
    "evidence_snippet",
    "extraction_status",
    "country_or_region",
    "locality_name",
    "latitude",
    "longitude",
    "coordinate_source",
    "population_label",
    "colour_state",
    "colour_state_verbatim",
    "sample_size",
    "locality_evidence_quote",
    "page_or_table",
    "extractor_notes",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True)
    parser.add_argument("--works", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    classification = pd.read_csv(args.classification)
    works = pd.read_csv(args.works)

    required_class = {"canonical_name", "family", "enriched_scale", "classification_source"}
    required_works = {
        "canonical_name",
        "openalex_id",
        "title",
        "year",
        "doi",
        "landing_url",
        "score",
        "evidence_snippet",
        "direct_colour_signal",
        "artificial_signal",
    }
    missing = sorted((required_class - set(classification.columns)) | (required_works - set(works.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    classified = classification[
        classification["enriched_scale"].isin(["within_population", "among_population"])
    ][["canonical_name", "family", "enriched_scale", "classification_source"]].drop_duplicates("canonical_name")

    eligible = works.copy()
    eligible["score"] = pd.to_numeric(eligible["score"], errors="coerce").fillna(0)
    eligible["direct_colour_signal"] = pd.to_numeric(
        eligible["direct_colour_signal"], errors="coerce"
    ).fillna(0)
    eligible["artificial_signal"] = pd.to_numeric(
        eligible["artificial_signal"], errors="coerce"
    ).fillna(0)
    eligible = eligible[
        (eligible["direct_colour_signal"] == 1)
        & (eligible["artificial_signal"] == 0)
    ].copy()

    queue = classified.merge(eligible, on="canonical_name", how="left", suffixes=("_class", "_work"))
    queue["family"] = queue.get("family_class", queue.get("family", "")).fillna(
        queue.get("family_work", "")
    )
    queue["spatial_scale"] = queue["enriched_scale"]
    queue["evidence_score"] = pd.to_numeric(queue["score"], errors="coerce").fillna(0)

    # Among-population records are primary because they are most informative for testing whether
    # documented colour states occupy distinct environments. Within-population records remain as
    # negative/contrast evidence and for checking co-occurrence claims.
    queue["scale_priority"] = queue["spatial_scale"].map(
        {"among_population": 0, "within_population": 1}
    ).fillna(2)
    queue["source_priority"] = queue["classification_source"].map(
        {"baseline_unambiguous": 0, "high_confidence_enrichment": 1}
    ).fillna(2)
    queue = queue.sort_values(
        ["scale_priority", "source_priority", "evidence_score", "canonical_name"],
        ascending=[True, True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    queue["priority_rank"] = range(1, len(queue) + 1)

    defaults = {
        "extraction_status": "not_started",
        "country_or_region": "",
        "locality_name": "",
        "latitude": "",
        "longitude": "",
        "coordinate_source": "",
        "population_label": "",
        "colour_state": "",
        "colour_state_verbatim": "",
        "sample_size": "",
        "locality_evidence_quote": "",
        "page_or_table": "",
        "extractor_notes": "",
    }
    for column, value in defaults.items():
        queue[column] = value

    for column in ("openalex_id", "title", "year", "doi", "landing_url", "evidence_snippet"):
        if column not in queue:
            queue[column] = ""

    out = queue[OUTPUT_COLUMNS]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    manifest = {
        "classified_species": int(classified["canonical_name"].nunique()),
        "queue_rows": int(len(out)),
        "species_with_queue_rows": int(out["canonical_name"].nunique()),
        "among_population_species": int(
            classified.loc[classified["enriched_scale"] == "among_population", "canonical_name"].nunique()
        ),
        "within_population_species": int(
            classified.loc[classified["enriched_scale"] == "within_population", "canonical_name"].nunique()
        ),
        "among_population_queue_rows": int((out["spatial_scale"] == "among_population").sum()),
        "rows_without_retained_work": int(out["openalex_id"].isna().sum()),
        "required_manual_fields": [
            "locality_name",
            "latitude",
            "longitude",
            "population_label",
            "colour_state",
            "locality_evidence_quote",
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
