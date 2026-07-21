#!/usr/bin/env python3
"""Build a source-level audit table for high-confidence enrichment classifications.

The output exposes every retained work that contributed a within- or among-population
signal, so classifications can be checked against titles and evidence snippets rather
than inferred from species-level aggregates alone.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True)
    parser.add_argument("--works", required=True)
    parser.add_argument("--scale-dataset")
    parser.add_argument("--out", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    classification = pd.read_csv(args.classification)
    works = pd.read_csv(args.works)

    required_class = {
        "canonical_name", "family", "enriched_scale", "classification_source"
    }
    required_works = {
        "canonical_name", "openalex_id", "title", "year", "doi", "landing_url",
        "score", "evidence_snippet", "within_signal", "among_signal",
        "direct_colour_signal", "artificial_signal",
    }
    missing = sorted(
        (required_class - set(classification.columns))
        | (required_works - set(works.columns))
    )
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    targets = classification.loc[
        classification["classification_source"].eq("high_confidence_enrichment"),
        ["canonical_name", "family", "enriched_scale", "classification_source"],
    ].drop_duplicates("canonical_name")

    if args.scale_dataset:
        scale = pd.read_csv(args.scale_dataset)
        required_scale = {"canonical_name", "n_climate_cells"}
        missing_scale = sorted(required_scale - set(scale.columns))
        if missing_scale:
            raise ValueError(f"Missing scale columns: {missing_scale}")
        targets = targets.merge(
            scale[["canonical_name", "n_climate_cells"]].drop_duplicates("canonical_name"),
            on="canonical_name",
            how="left",
        )
    else:
        targets["n_climate_cells"] = pd.NA

    for column in (
        "score", "within_signal", "among_signal", "direct_colour_signal", "artificial_signal"
    ):
        works[column] = pd.to_numeric(works[column], errors="coerce").fillna(0)

    retained = works.loc[
        works["canonical_name"].isin(targets["canonical_name"])
        & works["direct_colour_signal"].eq(1)
        & works["artificial_signal"].eq(0)
        & works["score"].ge(20)
        & (works["within_signal"].eq(1) | works["among_signal"].eq(1))
    ].copy()

    audit = targets.merge(retained, on="canonical_name", how="left", suffixes=("", "_work"))
    audit["claimed_signal"] = audit["enriched_scale"].map({
        "within_population": "within_signal",
        "among_population": "among_signal",
    })
    audit["signal_matches_classification"] = (
        (audit["enriched_scale"].eq("within_population") & audit["within_signal"].eq(1))
        | (audit["enriched_scale"].eq("among_population") & audit["among_signal"].eq(1))
    )
    audit["audit_status"] = audit["signal_matches_classification"].map(
        {True: "matching_source_present", False: "manual_source_review_required"}
    )
    audit["scale_priority"] = audit["enriched_scale"].map(
        {"among_population": 0, "within_population": 1}
    ).fillna(2)
    audit = audit.sort_values(
        ["scale_priority", "audit_status", "n_climate_cells", "canonical_name", "score"],
        ascending=[True, False, False, True, False],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)
    audit["audit_priority_rank"] = range(1, len(audit) + 1)

    columns = [
        "audit_priority_rank", "canonical_name", "family", "enriched_scale",
        "classification_source", "n_climate_cells", "claimed_signal",
        "signal_matches_classification", "audit_status", "openalex_id", "title", "year", "doi",
        "landing_url", "score", "within_signal", "among_signal",
        "evidence_snippet",
    ]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    audit[columns].to_csv(out_path, index=False)

    matched_species = audit.loc[
        audit["signal_matches_classification"].fillna(False), "canonical_name"
    ].nunique()
    target_species = int(targets["canonical_name"].nunique())
    manifest = {
        "target_species": target_species,
        "audit_rows": int(len(audit)),
        "species_with_matching_source": int(matched_species),
        "species_requiring_manual_source_review": int(target_species - matched_species),
        "among_population_targets": int(targets["enriched_scale"].eq("among_population").sum()),
        "within_population_targets": int(targets["enriched_scale"].eq("within_population").sum()),
        "climate_priority_available": bool(args.scale_dataset),
        "semantic_guard": (
            "This table exposes retained source-level evidence for manual verification; "
            "a regex signal is not treated as a substitute for reading the source, and "
            "missing matching sources are retained as explicit manual-review targets."
        ),
    }
    Path(args.manifest).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
