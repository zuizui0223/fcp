#!/usr/bin/env python3
"""Test whether the polymorphism-spatial gradient survives sampled-span opportunity controls."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_spatial_span_adjusted_step6_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISC_ATTR = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
DISC_SPATIAL = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909" / "species_membership_and_observed_rho.csv"
DISC_JSON = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909" / "result.json"
RES_SPATIAL = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "reserve_species_membership_and_observed_metrics.csv"
RES_JSON = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "result.json"
RES_RAW = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"

SEED = 20260910
N_PERM = 20_000
N_BINS = 5
EARTH_KM = 6371.0088


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(stats.spearmanr(np.asarray(x, float), np.asarray(y, float)).statistic)


def rank_residual(values: np.ndarray, controls: np.ndarray) -> np.ndarray:
    y = stats.rankdata(np.asarray(values, float)).astype(float)
    cols = [np.ones(len(y), float)]
    for j in range(controls.shape[1]):
        cols.append(stats.rankdata(np.asarray(controls[:, j], float)).astype(float))
    X = np.column_stack(cols)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    return y - X @ beta


def partial_spearman_permutation(
    d: np.ndarray,
    y: np.ndarray,
    span: np.ndarray,
    n_classifiable: np.ndarray,
    *,
    seed: int,
) -> dict[str, Any]:
    arrays = [np.asarray(v, float) for v in [d, y, span, n_classifiable]]
    keep = np.logical_and.reduce([np.isfinite(v) for v in arrays])
    d0, y0, s0, n0 = [v[keep] for v in arrays]
    controls = np.column_stack([s0, n0])
    ed = rank_residual(d0, controls)
    ey = rank_residual(y0, controls)
    denom = float(np.linalg.norm(ed) * np.linalg.norm(ey))
    if denom <= 0:
        raise RuntimeError("zero residual norm in partial Spearman")
    observed = float(ed @ ey / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        ep = ed[rng.permutation(len(ed))]
        null[i] = float(ep @ ey / (np.linalg.norm(ep) * np.linalg.norm(ey)))
    p_two = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    p_upper = float((1 + np.sum(null >= observed - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(d0)),
        "partial_spearman": observed,
        "p_two_sided": p_two,
        "p_upper_directional": p_upper,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def equal_count_bins(span: np.ndarray, species: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    order = np.lexsort((np.asarray(species, str), np.asarray(span, float)))
    bins = np.empty(len(order), int)
    for rank, idx in enumerate(order):
        bins[idx] = min(n_bins - 1, (rank * n_bins) // len(order))
    counts = np.bincount(bins, minlength=n_bins)
    if np.any(counts < 2):
        raise RuntimeError(f"span bins unexpectedly sparse: {counts.tolist()}")
    return bins


def conditioned_binary_contrast(
    response: np.ndarray,
    membership: np.ndarray,
    span: np.ndarray,
    species: np.ndarray,
    *,
    seed: int,
) -> dict[str, Any]:
    y = np.asarray(response, float)
    m = np.asarray(membership, bool)
    s = np.asarray(span, float)
    names = np.asarray(species, str)
    keep = np.isfinite(y) & np.isfinite(s)
    y, m, s, names = y[keep], m[keep], s[keep], names[keep]
    bins = equal_count_bins(s, names)
    if m.sum() == 0 or (~m).sum() == 0:
        raise RuntimeError("binary contrast has empty side")
    observed = float(y[m].mean() - y[~m].mean())
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    bin_indices = [np.flatnonzero(bins == b) for b in range(N_BINS)]
    for i in range(N_PERM):
        mp = m.copy()
        for idx in bin_indices:
            mp[idx] = m[idx][rng.permutation(len(idx))]
        null[i] = float(y[mp].mean() - y[~mp].mean())
    p_upper = float((1 + np.sum(null >= observed - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(y)),
        "n_polymorphic": int(m.sum()),
        "n_complement": int((~m).sum()),
        "observed_mean_difference": observed,
        "p_upper": p_upper,
        "span_bin_counts": {str(i): int(np.sum(bins == i)) for i in range(N_BINS)},
        "span_bin_polymorphic_counts": {str(i): int(np.sum(m[bins == i])) for i in range(N_BINS)},
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def max_great_circle_km(lat: np.ndarray, lon: np.ndarray) -> float:
    lat = np.radians(np.asarray(lat, float))
    lon = np.radians(np.asarray(lon, float))
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2.0) ** 2
    a = np.clip(a, 0.0, 1.0)
    dist = 2.0 * EARTH_KM * np.arcsin(np.sqrt(a))
    return float(np.nanmax(dist))


def reserve_spans(eligible_species: set[str]) -> pd.DataFrame:
    raw = pd.read_csv(RES_RAW, usecols=["species", "inat_taxon_id", "latitude", "longitude"])
    raw = raw[raw["species"].astype(str).isin(eligible_species)].copy()
    raw["latitude"] = pd.to_numeric(raw["latitude"], errors="coerce")
    raw["longitude"] = pd.to_numeric(raw["longitude"], errors="coerce")
    rows: list[dict[str, Any]] = []
    for (taxon, species), g in raw.groupby(["inat_taxon_id", "species"], sort=True):
        if len(g) != 100:
            raise RuntimeError(f"reserve geometry row count for {species}: {len(g)} != 100")
        if not np.isfinite(g[["latitude", "longitude"]].to_numpy(float)).all():
            raise RuntimeError(f"nonfinite reserve geometry for {species}")
        span = max_great_circle_km(g["latitude"].to_numpy(float), g["longitude"].to_numpy(float))
        rows.append({
            "inat_taxon_id": int(taxon),
            "species": str(species),
            "reserve_photo_rows": int(len(g)),
            "reserve_maximum_span_km_all_fixed_photos": span,
            "reserve_log1p_span": float(np.log1p(span)),
        })
    out = pd.DataFrame(rows)
    if len(out) != len(eligible_species) or out["species"].nunique() != len(eligible_species):
        raise RuntimeError(f"reserve span species mismatch: {len(out)} vs {len(eligible_species)}")
    return out


def describe_frame(df: pd.DataFrame, *, dcol: str, ycol: str, scol: str, ncol: str) -> dict[str, Any]:
    return {
        "n_species": int(len(df)),
        "rho_D_span": spearman(df[dcol].to_numpy(float), df[scol].to_numpy(float)),
        "rho_spatial_span": spearman(df[ycol].to_numpy(float), df[scol].to_numpy(float)),
        "rho_D_spatial": spearman(df[dcol].to_numpy(float), df[ycol].to_numpy(float)),
        "span_km_min": float(np.expm1(df[scol].min())),
        "span_km_median": float(np.expm1(df[scol].median())),
        "span_km_max": float(np.expm1(df[scol].max())),
        "n_classifiable_min": int(df[ncol].min()),
        "n_classifiable_max": int(df[ncol].max()),
    }


def run_frame(
    df: pd.DataFrame,
    *,
    dcol: str,
    second_col: str,
    ycol: str,
    scol: str,
    ncol: str,
    seed_base: int,
) -> dict[str, Any]:
    D = df[dcol].to_numpy(float)
    n = df[ncol].to_numpy(float)
    Dunb = D * n / (n - 1.0)
    y = df[ycol].to_numpy(float)
    span = df[scol].to_numpy(float)
    membership = df[second_col].to_numpy(float) >= 0.10
    species = df["species"].astype(str).to_numpy()
    return {
        "descriptive": describe_frame(df, dcol=dcol, ycol=ycol, scol=scol, ncol=ncol),
        "continuous_partial": partial_spearman_permutation(D, y, span, n, seed=seed_base + 1),
        "continuous_partial_D_unbiased": partial_spearman_permutation(Dunb, y, span, n, seed=seed_base + 2),
        "span_conditioned_second_ge_0_10_contrast": conditioned_binary_contrast(
            y, membership, span, species, seed=seed_base + 3
        ),
    }


def main() -> None:
    disc_attr = pd.read_csv(DISC_ATTR)
    disc_sp = pd.read_csv(DISC_SPATIAL)
    disc_meta = json.loads(DISC_JSON.read_text(encoding="utf-8"))
    res_sp = pd.read_csv(RES_SPATIAL)
    res_meta = json.loads(RES_JSON.read_text(encoding="utf-8"))

    if len(disc_attr) != 369 or len(disc_sp) != 369:
        raise RuntimeError("discovery species count drift")
    if len(res_sp) != 363:
        raise RuntimeError("reserve species count drift")

    disc = disc_sp.merge(
        disc_attr[["species", "log1p_span_primary", "n_classifiable"]],
        on="species",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_attr"),
    )
    if len(disc) != 369:
        raise RuntimeError("discovery span/spatial join lost species")
    if not np.array_equal(disc["n_classifiable"].astype(int), disc["n_classifiable_attr"].astype(int)):
        raise RuntimeError("discovery n_classifiable mismatch")
    disc = disc.drop(columns=["n_classifiable_attr"]).sort_values("species").reset_index(drop=True)

    raw_disc_rho = spearman(disc["D"].to_numpy(float), disc["spatial_observed_rho"].to_numpy(float))
    expected_disc = float(disc_meta["continuous_D_diagnostic"]["rho_D_species_spatial_rho"])
    if abs(raw_disc_rho - expected_disc) > 1e-12:
        raise RuntimeError(f"discovery raw D-spatial anchor mismatch: {raw_disc_rho} vs {expected_disc}")

    rspan = reserve_spans(set(res_sp["species"].astype(str)))
    reserve = res_sp.merge(rspan, on=["inat_taxon_id", "species"], how="inner", validate="one_to_one")
    if len(reserve) != 363:
        raise RuntimeError("reserve span/spatial join lost species")
    reserve = reserve.sort_values("species").reset_index(drop=True)

    raw_res_rho = spearman(reserve["D"].to_numpy(float), reserve["rho_primary"].to_numpy(float))
    expected_res = float(res_meta["continuous_D_diagnostic"]["rho_D_species_primary_spatial_rho"])
    if abs(raw_res_rho - expected_res) > 1e-12:
        raise RuntimeError(f"reserve raw D-spatial anchor mismatch: {raw_res_rho} vs {expected_res}")

    discovery = run_frame(
        disc,
        dcol="D",
        second_col="second_fraction",
        ycol="spatial_observed_rho",
        scol="log1p_span_primary",
        ncol="n_classifiable",
        seed_base=SEED + 100,
    )
    reserve_primary = run_frame(
        reserve,
        dcol="D",
        second_col="second_fraction",
        ycol="rho_primary",
        scol="reserve_log1p_span",
        ncol="n_classifiable",
        seed_base=SEED + 200,
    )
    reserve_background = run_frame(
        reserve,
        dcol="D",
        second_col="second_fraction",
        ycol="rho_matched_background_differential",
        scol="reserve_log1p_span",
        ncol="n_classifiable",
        seed_base=SEED + 300,
    )

    disc_partial = float(discovery["continuous_partial"]["partial_spearman"])
    res_partial = float(reserve_primary["continuous_partial"]["partial_spearman"])
    bg_partial = float(reserve_background["continuous_partial"]["partial_spearman"])
    span_alt_not_sufficient = bool(disc_partial > 0 and res_partial > 0 and bg_partial >= 0)
    strong_adjusted_replication = bool(
        disc_partial > 0
        and discovery["continuous_partial"]["p_two_sided"] < 0.05
        and res_partial > 0
        and reserve_primary["continuous_partial"]["p_two_sided"] < 0.05
    )
    flower_specific_adjusted_support = bool(
        bg_partial > 0 and reserve_background["continuous_partial"]["p_two_sided"] < 0.05
    )

    result = {
        "analysis": "polymorphism_spatial_span_adjusted_step6",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_SPATIAL_SPAN_ADJUSTED_STEP6_PROTOCOL_20260910.md",
        "inferential_role": "post_step4_robustness_diagnostic",
        "new_image_acquisition": False,
        "new_spatial_permutations_generated": False,
        "new_species_level_permutations": N_PERM,
        "discovery": discovery,
        "reserve_primary": reserve_primary,
        "reserve_matched_background_differential": reserve_background,
        "decision": {
            "sampled_span_alternative_not_sufficient_directionally": span_alt_not_sufficient,
            "strong_adjusted_discovery_plus_reserve_replication_p_lt_0_05": strong_adjusted_replication,
            "flower_specific_adjusted_support_p_lt_0_05": flower_specific_adjusted_support,
            "claim_boundary": "Adjustment can show that sampled span and classifiable-photo count do not fully account for the observed association; it does not identify causal range-size effects, adaptation, selection, population-genetic differentiation, or a shared boundary.",
        },
    }

    disc_out = disc.copy()
    disc_out["D_unbiased"] = disc_out["D"] * disc_out["n_classifiable"] / (disc_out["n_classifiable"] - 1.0)
    reserve_out = reserve.copy()
    reserve_out["D_unbiased"] = reserve_out["D"] * reserve_out["n_classifiable"] / (reserve_out["n_classifiable"] - 1.0)
    disc_out.to_csv(OUT / "discovery_species_span_adjusted_input.csv", index=False)
    reserve_out.to_csv(OUT / "reserve_species_span_adjusted_input.csv", index=False)
    write_json(OUT / "result.json", result)

    lines = [
        "# Polymorphism spatial organization — Step 6 span-adjusted result",
        "",
        f"**Sampled-span alternative not sufficient directionally: `{span_alt_not_sufficient}`.**",
        f"**Strong adjusted discovery + reserve replication (two-sided p<0.05 in both): `{strong_adjusted_replication}`.**",
        f"**Reserve flower-minus-background adjusted support (two-sided p<0.05): `{flower_specific_adjusted_support}`.**",
        "",
        "## Continuous partial Spearman controlling sampled span + n_classifiable",
        "",
        f"- discovery: partial rho = **{disc_partial:.6f}**, two-sided p = **{discovery['continuous_partial']['p_two_sided']:.6g}**",
        f"- reserve primary: partial rho = **{res_partial:.6f}**, two-sided p = **{reserve_primary['continuous_partial']['p_two_sided']:.6g}**",
        f"- reserve matched-background differential: partial rho = **{bg_partial:.6f}**, two-sided p = **{reserve_background['continuous_partial']['p_two_sided']:.6g}**",
        "",
        "## Span-conditioned second-morph >=10% contrasts",
        "",
        f"- discovery primary delta = **{discovery['span_conditioned_second_ge_0_10_contrast']['observed_mean_difference']:.6f}**, p_upper = **{discovery['span_conditioned_second_ge_0_10_contrast']['p_upper']:.6g}**",
        f"- reserve primary delta = **{reserve_primary['span_conditioned_second_ge_0_10_contrast']['observed_mean_difference']:.6f}**, p_upper = **{reserve_primary['span_conditioned_second_ge_0_10_contrast']['p_upper']:.6g}**",
        f"- reserve matched-background delta = **{reserve_background['span_conditioned_second_ge_0_10_contrast']['observed_mean_difference']:.6f}**, p_upper = **{reserve_background['span_conditioned_second_ge_0_10_contrast']['p_upper']:.6g}**",
        "",
        "## Raw opportunity diagnostics",
        "",
        f"- discovery rho(D, span) = **{discovery['descriptive']['rho_D_span']:.6f}**; rho(spatial, span) = **{discovery['descriptive']['rho_spatial_span']:.6f}**; raw rho(D, spatial) = **{discovery['descriptive']['rho_D_spatial']:.6f}**",
        f"- reserve rho(D, span) = **{reserve_primary['descriptive']['rho_D_span']:.6f}**; rho(spatial, span) = **{reserve_primary['descriptive']['rho_spatial_span']:.6f}**; raw rho(D, spatial) = **{reserve_primary['descriptive']['rho_D_spatial']:.6f}**",
        "",
        "## Claim boundary",
        "",
        "This is a post-Step-4 robustness diagnostic. Survival supports only that sampled geographic span and classifiable-photo count do not fully account for the polymorphism-spatial association in these fixed photo-derived frames. It does not establish adaptation, selection, causal range-size effects, population-genetic differentiation, or a shared global boundary.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "decision": result["decision"],
        "discovery_partial": discovery["continuous_partial"],
        "reserve_primary_partial": reserve_primary["continuous_partial"],
        "reserve_background_partial": reserve_background["continuous_partial"],
    }, indent=2))


if __name__ == "__main__":
    main()
