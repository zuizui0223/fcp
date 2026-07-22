#!/usr/bin/env python3
"""Integrate baseline and high-confidence literature evidence for spatial scale.

Existing unambiguous classifications are preserved. High-confidence enrichment is used
only to resolve baseline mixed/unclear cases. Ambiguous or conflicting evidence remains
available for manual review rather than being forced into a binary class.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

WITHIN = re.compile(
    r"\b(within[- ]population|same population|coexist|co-occurr|morph frequenc|"
    r"frequency[- ]dependent|polymorphic population|multiple (?:colou?r|flower) morphs|"
    r"colour morphs|color morphs)\b",
    re.I,
)
GEOGRAPHIC = re.compile(
    r"\b(geographic|spatial variation|among populations|between populations|"
    r"population differentiation|cline|hybrid zone|range edge|local adaptation|"
    r"translocation|regional variation)\b",
    re.I,
)


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if sd and np.isfinite(sd) else x * 0


def baseline_scale(row: pd.Series) -> str:
    text = " ".join(
        str(row.get(column, "") or "")
        for column in ("best_title", "best_match_evidence", "review_reason")
    )
    within = bool(WITHIN.search(text))
    among = bool(GEOGRAPHIC.search(text))
    if within and not among:
        return "within_population"
    if among and not within:
        return "among_population"
    if within and among:
        return "mixed"
    return "unclear"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review", required=True)
    parser.add_argument("--works", required=True)
    parser.add_argument("--geography", required=True)
    parser.add_argument("--outdir", default="analysis_outputs/evidence_scale_enriched")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    review = pd.read_csv(args.review)
    works = pd.read_csv(args.works)
    geography = pd.read_csv(args.geography)

    required_review = {"canonical_name", "family"}
    required_works = {
        "canonical_name",
        "openalex_id",
        "direct_colour_signal",
        "artificial_signal",
        "score",
        "within_signal",
        "among_signal",
    }
    required_geo = {
        "canonical_name",
        "latitudinal_range",
        "absolute_latitude",
        "gbif_occurrences",
    }
    for label, frame, required in (
        ("review", review, required_review),
        ("works", works, required_works),
        ("geography", geography, required_geo),
    ):
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"{label} input missing required columns: {missing}")

    review = review.drop_duplicates("canonical_name").copy()
    review["baseline_scale"] = review.apply(baseline_scale, axis=1)

    for column in ("direct_colour_signal", "artificial_signal", "score", "within_signal", "among_signal"):
        works[column] = pd.to_numeric(works[column], errors="coerce").fillna(0)

    eligible = works[
        (works["direct_colour_signal"] == 1)
        & (works["artificial_signal"] == 0)
        & (works["score"] >= 20)
    ].copy()

    if eligible.empty:
        aggregated = pd.DataFrame(
            columns=["canonical_name", "n_eligible", "enrichment_within", "enrichment_among"]
        )
    else:
        aggregated = (
            eligible.groupby("canonical_name", as_index=False)
            .agg(
                n_eligible=("openalex_id", "nunique"),
                enrichment_within=("within_signal", "max"),
                enrichment_among=("among_signal", "max"),
            )
        )

    classified = review.merge(aggregated, on="canonical_name", how="left")
    for column in ("n_eligible", "enrichment_within", "enrichment_among"):
        classified[column] = classified[column].fillna(0).astype(int)

    def integrate(row: pd.Series) -> tuple[str, str]:
        baseline = row["baseline_scale"]
        within = int(row["enrichment_within"]) == 1
        among = int(row["enrichment_among"]) == 1

        if baseline in {"within_population", "among_population"}:
            conflicting = (
                baseline == "within_population" and among
            ) or (
                baseline == "among_population" and within
            )
            if conflicting:
                return "mixed", "baseline_enrichment_conflict"
            return baseline, "baseline_unambiguous"

        if within and not among:
            return "within_population", "high_confidence_enrichment"
        if among and not within:
            return "among_population", "high_confidence_enrichment"
        if within and among:
            return "mixed", "high_confidence_enrichment_mixed"
        return baseline, "unresolved"

    integrated = classified.apply(integrate, axis=1, result_type="expand")
    classified["enriched_scale"] = integrated[0]
    classified["classification_source"] = integrated[1]
    classified["manual_review_required"] = classified["enriched_scale"].isin(["mixed", "unclear"])

    classification_path = outdir / "high_confidence_enriched_classification.csv"
    manual_path = outdir / "high_confidence_manual_review.csv"
    classified.to_csv(classification_path, index=False)
    classified.loc[classified["manual_review_required"]].to_csv(manual_path, index=False)

    merged = classified.merge(
        geography,
        on="canonical_name",
        how="left",
        suffixes=("_review", "_geo"),
    )
    model_data = merged[
        merged["enriched_scale"].isin(["within_population", "among_population"])
    ].dropna(subset=["latitudinal_range", "absolute_latitude", "gbif_occurrences"]).copy()

    counts = model_data["enriched_scale"].value_counts().to_dict()
    result: dict[str, object] = {
        "status": "not_run",
        "reason": "fewer than 15 rows in either unambiguous class",
        "n_rows": int(len(model_data)),
        "class_counts": counts,
    }

    if counts.get("within_population", 0) >= 15 and counts.get("among_population", 0) >= 15:
        model_data["is_within"] = (
            model_data["enriched_scale"] == "within_population"
        ).astype(int)
        model_data["range_z"] = zscore(
            np.log1p(pd.to_numeric(model_data["latitudinal_range"], errors="coerce").clip(lower=0))
        )
        model_data["latitude_z"] = zscore(model_data["absolute_latitude"])
        model_data["gbif_z"] = zscore(
            np.log1p(pd.to_numeric(model_data["gbif_occurrences"], errors="coerce").clip(lower=0))
        )
        family_column = "family_review" if "family_review" in model_data.columns else "family"
        try:
            model = smf.glm(
                "is_within ~ range_z + latitude_z + gbif_z",
                model_data,
                family=sm.families.Binomial(),
            ).fit(
                cov_type="cluster",
                cov_kwds={"groups": model_data[family_column].fillna("unknown")},
            )
            ci = model.conf_int()
            table = pd.DataFrame(
                {
                    "term": model.params.index,
                    "estimate": model.params.values,
                    "std_error": model.bse.values,
                    "odds_ratio": np.exp(model.params.values),
                    "ci_low": np.exp(ci[0].values),
                    "ci_high": np.exp(ci[1].values),
                    "p_value": model.pvalues.values,
                }
            )
            table.to_csv(outdir / "high_confidence_within_vs_among_model.csv", index=False)
            indexed = table.set_index("term")
            result = {
                "status": "complete",
                "n_rows": int(len(model_data)),
                "class_counts": counts,
                "formula": "is_within ~ range_z + latitude_z + gbif_z",
                "cluster_column": family_column,
                "range_odds_ratio": float(indexed.loc["range_z", "odds_ratio"]),
                "range_p_value": float(indexed.loc["range_z", "p_value"]),
            }
        except Exception as exc:
            result = {
                "status": "not_run",
                "reason": f"model fitting failed: {type(exc).__name__}: {exc}",
                "n_rows": int(len(model_data)),
                "class_counts": counts,
            }

    manifest = {
        "n_confirmed": int(len(review)),
        "eligible_natural_direct_works": int(len(eligible)),
        "species_with_eligible_works": int(eligible["canonical_name"].nunique()),
        "baseline_counts": review["baseline_scale"].value_counts().to_dict(),
        "scale_counts": classified["enriched_scale"].value_counts().to_dict(),
        "classification_source_counts": classified["classification_source"].value_counts().to_dict(),
        "manual_review_rows": int(classified["manual_review_required"].sum()),
        "model": result,
        "semantic_guard": (
            "existing unambiguous evidence is preserved; enrichment requires direct natural "
            "flower-colour evidence with score >=20; ambiguous evidence is not forced"
        ),
    }
    (outdir / "high_confidence_enriched_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
