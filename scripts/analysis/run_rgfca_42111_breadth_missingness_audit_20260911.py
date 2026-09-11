#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "docs/RGFCA_42111_BREADTH_POSTMEASUREMENT_MISSINGNESS_AUDIT_PROTOCOL_20260911.md"
CENSUS = ROOT / "results/rgfca_42111_breadth_depth_step8b_20260911/species_capacity_census.csv.gz"
BIOLOGICAL = {"white", "yellow_orange", "red_pink", "blue_purple"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--measured", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def capacity_bin(n: int) -> str:
    if n <= 0:
        return "0"
    if n == 1:
        return "1"
    if n <= 4:
        return "2_4"
    if n <= 9:
        return "5_9"
    if n <= 19:
        return "10_19"
    return "20_plus"


def breadth_block(rank: int) -> str:
    lo = ((int(rank) - 1) // 5000) * 5000 + 1
    hi = min(42111, lo + 4999)
    return f"{lo}_{hi}"


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt((p * (1 - p) / n) + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def source_class(row: pd.Series) -> str:
    v1 = bool(row["source_v1"])
    v2 = bool(row["source_v2"])
    if v1 and v2:
        return "both"
    if v1:
        return "v1_only"
    if v2:
        return "v2_only"
    raise RuntimeError("species lacks V1/V2 source membership")


def summarize_family(x: pd.DataFrame, family: str, level_col: str) -> tuple[pd.DataFrame, dict[str, object]]:
    rows: list[dict[str, object]] = []
    for level, g in x.groupby(level_col, sort=True, observed=True):
        n = len(g)
        class_n = int(g["classifiable"].sum())
        acq_n = int(g["acquisition_failure"].sum())
        roi_n = int(g["roi_nonadmission"].sum())
        mixed_n = int(g["mixed_uncertain"].sum())
        lo, hi = wilson(class_n, n)
        rows.append({
            "family": family,
            "level": str(level),
            "denominator": int(n),
            "image_acquisition_failure_n": acq_n,
            "image_acquisition_failure_fraction": float(acq_n / n),
            "roi_nonadmission_n": roi_n,
            "roi_nonadmission_fraction": float(roi_n / n),
            "mixed_uncertain_n": mixed_n,
            "mixed_uncertain_fraction": float(mixed_n / n),
            "classifiable_n": class_n,
            "classifiable_fraction": float(class_n / n),
            "classifiable_wilson_low": lo,
            "classifiable_wilson_high": hi,
        })
    out = pd.DataFrame(rows)
    values = out["classifiable_fraction"].to_numpy(float)
    spread = float(values.max() - values.min()) if len(values) else float("nan")
    summary = {
        "family": family,
        "levels": int(len(out)),
        "max_absolute_classifiable_fraction_difference": spread,
        "minimum_level_classifiable_fraction": float(values.min()) if len(values) else None,
        "maximum_level_classifiable_fraction": float(values.max()) if len(values) else None,
    }
    return out, summary


def make_design(x: pd.DataFrame, categorical: list[str]) -> tuple[np.ndarray, list[str]]:
    cols = [np.ones(len(x), dtype=float)]
    names = ["intercept"]
    for c in categorical:
        levels = sorted(x[c].astype(str).unique())
        for level in levels[1:]:
            cols.append((x[c].astype(str).to_numpy() == level).astype(float))
            names.append(f"{c}={level}")
    return np.column_stack(cols), names


def fit_logit_irls(X: np.ndarray, y: np.ndarray, ridge: float = 1e-6, max_iter: int = 100) -> np.ndarray:
    beta = np.zeros(X.shape[1], dtype=float)
    penalty = np.eye(X.shape[1]) * ridge
    penalty[0, 0] = 0.0
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -30.0, 30.0)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p * (1.0 - p), 1e-6, None)
        score = X.T @ (y - p) - penalty @ beta
        h = (X.T * w) @ X + penalty
        step = np.linalg.solve(h, score)
        beta_new = beta + step
        if np.max(np.abs(step)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new
    return beta


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not PROTOCOL.exists():
        raise RuntimeError("frozen missingness protocol is missing")
    x = pd.read_csv(args.measured)
    c = pd.read_csv(CENSUS, usecols=[
        "inat_taxon_id", "after_observer_cap", "source_v1", "source_v2"
    ])
    if len(x) != 42111 or x["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("measured table is not exact 42,111 species")
    if len(c) != 42111 or c["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("metadata census is not exact 42,111 species")
    x = x.merge(c, on="inat_taxon_id", how="left", validate="one_to_one")
    if x[["after_observer_cap", "source_v1", "source_v2"]].isna().any().any():
        raise RuntimeError("metadata merge left missing strata")

    x["classifiable"] = (
        x["morph"].astype(str).isin(BIOLOGICAL)
        & x["measurement_status"].astype(str).eq("classified_four_state_morph")
    )
    x["acquisition_failure"] = x["measurement_status"].astype(str).isin([
        "image_acquisition_failed", "anchor_resolution_failed"
    ])
    x["roi_nonadmission"] = (
        ~x["acquisition_failure"]
        & ~x["classifiable"]
        & x["roi_status"].astype(str).ne("")
    )
    x["mixed_uncertain"] = x["morph"].astype(str).eq("mixed_uncertain")
    x["discovery_source"] = x.apply(source_class, axis=1)
    x["capacity_bin"] = pd.to_numeric(x["after_observer_cap"], errors="raise").astype(int).map(capacity_bin)
    x["opened_group"] = np.where(x["opened_existing_species"].astype(str).str.lower().isin(["true", "1"]), "previously_opened_1000", "not_previously_opened")
    x["breadth_rank_block"] = pd.to_numeric(x["breadth_rank"], errors="raise").astype(int).map(breadth_block)
    x["genus_label"] = x["species"].astype(str).str.split().str[0]
    genus_n = x.groupby("genus_label", observed=True)["inat_taxon_id"].transform("size")
    x["genus_size_class"] = np.where(genus_n.eq(1), "singleton_genus_label", "non_singleton_genus_label")

    families = [
        ("discovery_source", "discovery_source"),
        ("capacity_bin", "capacity_bin"),
        ("previously_opened", "opened_group"),
        ("breadth_rank_block", "breadth_rank_block"),
        ("genus_size_class", "genus_size_class"),
    ]
    tables = []
    family_summaries: dict[str, object] = {}
    for family, col in families:
        t, s = summarize_family(x, family, col)
        tables.append(t)
        family_summaries[family] = s
    strata = pd.concat(tables, ignore_index=True)

    # Prespecified metadata-only classifiability model. No colour category enters X.
    model_cols = ["discovery_source", "capacity_bin", "opened_group", "breadth_rank_block", "genus_size_class"]
    X, names = make_design(x, model_cols)
    y = x["classifiable"].astype(int).to_numpy(float)
    beta = fit_logit_irls(X, y)
    p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30.0, 30.0)))
    p = np.clip(p, 0.01, 0.99)
    x["classifiability_probability_metadata_only"] = p
    x["ipw_classifiable"] = np.where(x["classifiable"], 1.0 / p, 0.0)

    # Sensitivity composition: model-based IPW among the same four admitted states.
    ipw_num = {m: float(x.loc[x["morph"].astype(str).eq(m), "ipw_classifiable"].sum()) for m in sorted(BIOLOGICAL)}
    ipw_den = float(sum(ipw_num.values()))
    ipw_comp = {m: (v / ipw_den if ipw_den > 0 else None) for m, v in ipw_num.items()}

    # Complete-stratum standardization. Standardization cells use only frozen metadata.
    cell_cols = model_cols
    cell = x.groupby(cell_cols, observed=True, dropna=False).agg(
        denominator=("inat_taxon_id", "size"),
        classifiable=("classifiable", "sum"),
    ).reset_index()
    for m in sorted(BIOLOGICAL):
        z = x.assign(_is=(x["morph"].astype(str) == m) & x["classifiable"]).groupby(cell_cols, observed=True, dropna=False)["_is"].sum().reset_index(name=f"n_{m}")
        cell = cell.merge(z, on=cell_cols, how="left", validate="one_to_one")
    complete = cell.loc[cell["classifiable"] > 0].copy()
    standardized_num = {m: float(((complete[f"n_{m}"] / complete["classifiable"]) * complete["denominator"]).sum()) for m in sorted(BIOLOGICAL)}
    standardized_den = float(complete["denominator"].sum())
    standardized_comp = {m: (v / standardized_den if standardized_den > 0 else None) for m, v in standardized_num.items()}

    strata.to_csv(args.output_dir / "missingness_by_frozen_stratum.csv", index=False, lineterminator="\n")
    cell.to_csv(args.output_dir / "standardization_cells.csv", index=False, lineterminator="\n")
    pd.DataFrame({"term": names, "coefficient": beta}).to_csv(args.output_dir / "classifiability_logit_coefficients.csv", index=False, lineterminator="\n")

    full_class = int(x["classifiable"].sum())
    result = {
        "analysis": "rgfca_42111_breadth_postmeasurement_missingness_audit",
        "status": "complete_prespecified_missingness_audit",
        "species_denominator": 42111,
        "classifiable_species": full_class,
        "classifiable_fraction": float(full_class / 42111),
        "family_selection_diagnostics": family_summaries,
        "metadata_only_classifiability_model_terms": names,
        "ipw_four_state_composition_sensitivity": ipw_comp,
        "complete_stratum_standardized_four_state_composition_sensitivity": standardized_comp,
        "complete_stratum_denominator": int(standardized_den),
        "primary_estimator_replaced": False,
        "flower_colour_used_as_classifiability_predictor": False,
        "anchor_replacement_performed": False,
        "claim_boundary": "Missingness diagnostic only. One photo per species does not identify modal colour, polymorphism prevalence, D, or C*/S*."
    }
    (args.output_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.output_dir / "RESULT.md").write_text(
        "# RGFCA 42,111 breadth — postmeasurement missingness audit\n\n"
        f"- species denominator: **42,111**\n"
        f"- classifiable: **{full_class:,} ({full_class/42111:.3%})**\n"
        "- primary estimator replaced: **false**\n"
        "- colour used in classifiability model: **false**\n"
        "- failed anchors replaced: **false**\n\n"
        "The audit diagnoses selection into the classifiable subset using only frozen metadata strata; it does not change the primary species-equal observed-state composition.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
