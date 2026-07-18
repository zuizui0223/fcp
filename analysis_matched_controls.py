#!/usr/bin/env python3
"""Compare confirmed FCP species with taxonomically matched screening controls.

Controls are not asserted to be monomorphic. The estimand is whether confirmed FCP
species differ geographically from same-genus/family species outside the candidate set.
Unknown controls therefore bias contrasts toward the null rather than becoming absences.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if sd and np.isfinite(sd) else x * 0


def fit_one(df: pd.DataFrame, label: str) -> dict:
    x = df.copy()
    x["range_z"] = zscore(np.log1p(x["latitudinal_range"].clip(lower=0)))
    x["latitude_z"] = zscore(x["absolute_latitude"])
    x["gbif_z"] = zscore(np.log1p(x["gbif_occurrences"].clip(lower=0)))
    exog = x[["range_z", "latitude_z", "gbif_z"]]
    model = ConditionalLogit(x["is_focal"], exog, groups=x["focal_species"]).fit(disp=False, maxiter=300)
    ci = model.conf_int()
    table = pd.DataFrame({
        "specification": label,
        "term": model.params.index,
        "estimate": model.params.values,
        "std_error": model.bse.values,
        "odds_ratio": np.exp(model.params.values),
        "ci_low": np.exp(ci[0].values),
        "ci_high": np.exp(ci[1].values),
        "p_value": model.pvalues.values,
        "n_rows": len(x),
        "n_strata": x["focal_species"].nunique(),
    })
    return {"table": table, "model": model, "n_rows": len(x), "n_strata": x["focal_species"].nunique()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review", required=True)
    ap.add_argument("--candidate-species", required=True)
    ap.add_argument("--controls", default="data/control_candidate_pool.csv")
    ap.add_argument("--geography", required=True)
    ap.add_argument("--outdir", default="analysis_outputs/matched_controls")
    args = ap.parse_args()

    review = pd.read_csv(args.review)
    candidates = set(pd.read_csv(args.candidate_species)["canonical_name"].dropna().astype(str))
    controls = pd.read_csv(args.controls)
    geo = pd.read_csv(args.geography)
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    focals = review[["canonical_name", "family"]].drop_duplicates().rename(columns={"canonical_name": "focal_species", "family": "focal_family"})
    controls = controls[controls["focal_species"].isin(set(focals["focal_species"]))].copy()
    controls["excluded_current_candidate"] = controls["control_candidate"].isin(candidates)
    screened = controls[~controls["excluded_current_candidate"]].copy()

    focal_rows = focals.assign(species=lambda d: d["focal_species"], is_focal=1, match_level="focal")
    control_rows = screened.rename(columns={"control_candidate": "species"}).assign(is_focal=0)
    cols = ["focal_species", "focal_family", "species", "is_focal", "match_level"]
    dat = pd.concat([focal_rows[cols], control_rows[cols]], ignore_index=True)
    dat = dat.merge(geo, left_on="species", right_on="canonical_name", how="left")
    for c in ["absolute_latitude", "latitudinal_range", "gbif_occurrences", "gbif_coordinate_records_sampled"]:
        dat[c] = pd.to_numeric(dat.get(c), errors="coerce")
    dat.to_csv(out / "matched_control_dataset.csv", index=False)

    complete = dat.dropna(subset=["absolute_latitude", "latitudinal_range", "gbif_occurrences"]).copy()
    valid = complete.groupby("focal_species")["is_focal"].agg(["sum", "count"])
    keep = set(valid[(valid["sum"] == 1) & (valid["count"] >= 2)].index)
    complete = complete[complete["focal_species"].isin(keep)].copy()

    specs = []
    specs.append(("all_taxonomic_controls", complete))
    genus_strata = set(screened.loc[screened["match_level"] == "same_genus", "focal_species"])
    genus = complete[(complete["match_level"].isin(["focal", "same_genus"])) & complete["focal_species"].isin(genus_strata)].copy()
    gv = genus.groupby("focal_species")["is_focal"].agg(["sum", "count"])
    genus = genus[genus["focal_species"].isin(set(gv[(gv["sum"] == 1) & (gv["count"] >= 2)].index))]
    specs.append(("same_genus_only", genus))
    if "gbif_coordinate_records_sampled" in complete:
        high = complete[(complete["is_focal"] == 1) | (complete["gbif_coordinate_records_sampled"] >= 50)].copy()
        hv = high.groupby("focal_species")["is_focal"].agg(["sum", "count"])
        high = high[high["focal_species"].isin(set(hv[(hv["sum"] == 1) & (hv["count"] >= 2)].index))]
        specs.append(("controls_min_50_coordinates", high))

    tables = []
    summaries = []
    for label, frame in specs:
        if frame["focal_species"].nunique() < 20:
            summaries.append({"specification": label, "status": "not_run", "n_strata": int(frame["focal_species"].nunique())})
            continue
        result = fit_one(frame, label)
        tables.append(result["table"])
        term = result["table"].set_index("term").loc["range_z"]
        summaries.append({
            "specification": label, "status": "complete", "n_rows": result["n_rows"], "n_strata": result["n_strata"],
            "range_odds_ratio": float(term["odds_ratio"]), "range_ci_low": float(term["ci_low"]),
            "range_ci_high": float(term["ci_high"]), "range_p_value": float(term["p_value"]),
        })
    if tables:
        pd.concat(tables, ignore_index=True).to_csv(out / "matched_control_conditional_logit.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "matched_control_summary.csv", index=False)

    qc = {
        "focal_species_requested": int(len(focals)),
        "control_rows_original": int(len(controls)),
        "control_rows_excluded_as_current_candidates": int(controls["excluded_current_candidate"].sum()),
        "control_rows_screening_comparators": int(len(screened)),
        "unique_screening_controls": int(screened["control_candidate"].nunique()),
        "complete_strata": int(complete["focal_species"].nunique()),
        "semantic_guard": "screening_controls_are_not_asserted_monomorphic",
        "estimand": "confirmed_FCP_species_vs_taxonomically_matched_species_outside_candidate_set",
        "specifications": summaries,
    }
    (out / "matched_control_manifest.json").write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(qc, ensure_ascii=False))


if __name__ == "__main__":
    main()
