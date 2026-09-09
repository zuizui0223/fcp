#!/usr/bin/env python3
"""Aggregate frozen reserve spatial null arrays over fixed polymorphism subsets."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
OUT = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909"
OUT.mkdir(parents=True, exist_ok=True)
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
METRICS = ["primary", "observer_pair_exclusion", "calendar_quarter_stratification", "matched_background_differential"]
OBS_COLS = {
    "primary": "rho_primary",
    "observer_pair_exclusion": "rho_observer_pair_exclusion",
    "calendar_quarter_stratification": "rho_calendar_quarter_stratification",
    "matched_background_differential": "rho_matched_background_differential",
}
EXPECTED_SPECIES = 363
EXPECTED_NULL = 999
PRIMARY_THRESHOLD = 0.10
SENS_THRESHOLD = 0.20
EXPECTED_WHOLE = {
    "primary": (0.025482606069841617, 0.001),
    "observer_pair_exclusion": (0.02521802297988336, 0.001),
    "calendar_quarter_stratification": (0.025482606069841617, 0.001),
    "matched_background_differential": (0.004477251893133873, 0.087),
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-root", type=Path, required=True)
    return p.parse_args()


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def reconstruct_polymorphism() -> pd.DataFrame:
    df = pd.read_csv(RESERVE, usecols=["inat_taxon_id", "species", "morph", "global_classifiable"])
    keep = as_bool(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    rows = []
    for (taxon, sp), g in df.loc[keep].groupby(["inat_taxon_id", "species"], sort=True):
        if len(g) < 40:
            continue
        cnt = Counter(g["morph"].astype(str))
        p = np.array([cnt.get(m, 0) / len(g) for m in MORPHS], float)
        ps = np.sort(p)[::-1]
        rows.append({
            "inat_taxon_id": int(taxon),
            "species": str(sp),
            "n_classifiable": int(len(g)),
            "D": float(1.0 - np.sum(p * p)),
            "second_fraction": float(ps[1]),
        })
    out = pd.DataFrame(rows).sort_values("inat_taxon_id").reset_index(drop=True)
    if len(out) != EXPECTED_SPECIES:
        raise RuntimeError(f"reserve polymorphism frame mismatch: {len(out)} != {EXPECTED_SPECIES}")
    dmax = float(out["D"].max())
    s10 = float((out["second_fraction"] >= 0.10).mean())
    s20 = float((out["second_fraction"] >= 0.20).mean())
    if abs(dmax - 0.680272) > 1e-5:
        raise RuntimeError(f"reserve D max fingerprint mismatch: {dmax}")
    if abs(s10 - 0.4160) > 0.001 or abs(s20 - 0.2479) > 0.001:
        raise RuntimeError(f"reserve second-morph fingerprint mismatch: {s10}, {s20}")
    return out


def load_artifacts(root: Path) -> tuple[pd.DataFrame, np.ndarray, dict]:
    shard_dirs = sorted([p for p in root.glob("reserve-inference-shard-*") if p.is_dir()])
    if len(shard_dirs) != 20:
        raise RuntimeError(f"expected 20 reserve shard dirs, found {len(shard_dirs)}")
    species_frames = []
    null_by_taxon: dict[int, np.ndarray] = {}
    for d in shard_dirs:
        species_path = d / "species.csv"
        null_path = d / "null.npz"
        receipt_path = d / "receipt.json"
        if not species_path.exists() or not null_path.exists() or not receipt_path.exists():
            raise RuntimeError(f"incomplete reserve shard artifact: {d}")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("metrics") != METRICS or int(receipt.get("permutations", -1)) != EXPECTED_NULL:
            raise RuntimeError(f"reserve shard receipt drift: {d}")
        sp = pd.read_csv(species_path)
        species_frames.append(sp)
        z = np.load(null_path)
        tids = z["taxon_ids"].astype(int)
        arr = z["null"].astype(float)
        if arr.shape != (len(tids), len(METRICS), EXPECTED_NULL):
            raise RuntimeError(f"null array shape mismatch {d}: {arr.shape}")
        for i, tid in enumerate(tids):
            if int(tid) in null_by_taxon:
                raise RuntimeError(f"duplicate taxon in null arrays: {tid}")
            null_by_taxon[int(tid)] = arr[i]
    species = pd.concat(species_frames, ignore_index=True)
    if species["inat_taxon_id"].duplicated().any() or len(species) != EXPECTED_SPECIES:
        raise RuntimeError(f"reserve observed species duplication/count mismatch: {len(species)}")
    species = species.sort_values("inat_taxon_id").reset_index(drop=True)
    tids = species["inat_taxon_id"].astype(int).to_numpy()
    if set(tids) != set(null_by_taxon) or len(null_by_taxon) != EXPECTED_SPECIES:
        raise RuntimeError("reserve observed/null taxon identity mismatch")
    null = np.stack([null_by_taxon[int(t)] for t in tids], axis=0)

    full_dir = root / "rgfca-reserve-replication-full-result-v1"
    full_result_path = full_dir / "result.json"
    full_species_path = full_dir / "species.csv"
    full_null_path = full_dir / "global_null.csv"
    if not (full_result_path.exists() and full_species_path.exists() and full_null_path.exists()):
        raise RuntimeError("missing reserve full-result artifact")
    full_result = json.loads(full_result_path.read_text(encoding="utf-8"))
    full_species = pd.read_csv(full_species_path).sort_values("inat_taxon_id").reset_index(drop=True)
    if len(full_species) != EXPECTED_SPECIES:
        raise RuntimeError("reserve full species count mismatch")
    for metric in METRICS:
        err = float(np.max(np.abs(species[OBS_COLS[metric]].to_numpy(float) - full_species[OBS_COLS[metric]].to_numpy(float))))
        if err > 1e-14:
            raise RuntimeError(f"reserve shard/full observed mismatch for {metric}: {err}")

    # Whole-frame reproduction gate for every metric.
    saved_global_null = pd.read_csv(full_null_path).sort_values("permutation_index")
    for j, metric in enumerate(METRICS):
        observed = float(species[OBS_COLS[metric]].mean())
        null_mean = null[:, j, :].mean(axis=0)
        p = float((1 + np.sum(null_mean >= observed)) / (EXPECTED_NULL + 1))
        exp_mean, exp_p = EXPECTED_WHOLE[metric]
        if abs(observed - exp_mean) > 1e-14 or p != exp_p:
            raise RuntimeError(f"whole reserve result mismatch for {metric}: {observed}, {p}")
        err = float(np.max(np.abs(saved_global_null[metric].to_numpy(float) - null_mean)))
        if err > 1e-14:
            raise RuntimeError(f"reserve full global-null mismatch for {metric}: {err}")
        saved = full_result["metrics"][metric]
        if abs(float(saved["mean_rho"]) - observed) > 1e-14 or float(saved["p_upper"]) != p:
            raise RuntimeError(f"reserve result.json mismatch for {metric}")
    return species, null, full_result


def subset_metric(species: pd.DataFrame, null: np.ndarray, members: np.ndarray, metric_index: int, label: str) -> tuple[dict, pd.DataFrame]:
    metric = METRICS[metric_index]
    obs_values = species.loc[members, OBS_COLS[metric]].to_numpy(float)
    observed = float(obs_values.mean())
    null_values = null[members, metric_index, :].mean(axis=0)
    p = float((1 + np.sum(null_values >= observed)) / (EXPECTED_NULL + 1))
    res = {
        "subset": label,
        "metric": metric,
        "n_species": int(np.sum(members)),
        "observed_mean_rho": observed,
        "observed_median_species_rho": float(np.median(obs_values)),
        "observed_positive_species_fraction": float(np.mean(obs_values > 0)),
        "null_mean": float(np.mean(null_values)),
        "null_sd": float(np.std(null_values, ddof=1)),
        "null_q025": float(np.quantile(null_values, 0.025)),
        "null_q975": float(np.quantile(null_values, 0.975)),
        "p_upper": p,
        "positive_supported_at_0_05": bool(observed > 0 and p < 0.05),
    }
    df = pd.DataFrame({"permutation_index": np.arange(EXPECTED_NULL), "mean_rho": null_values, "subset": label, "metric": metric})
    return res, df


def rho(x: np.ndarray, y: np.ndarray) -> float:
    return float(stats.spearmanr(x, y).statistic)


def main() -> None:
    args = parse_args()
    poly = reconstruct_polymorphism()
    species, null, _ = load_artifacts(args.artifact_root)
    joined = species.merge(poly, on=["inat_taxon_id", "species"], how="inner", validate="one_to_one")
    if len(joined) != EXPECTED_SPECIES:
        raise RuntimeError(f"reserve spatial/polymorphism join mismatch: {len(joined)}")
    # Reorder null arrays to the joined/sorted species order. Both are taxon-sorted by construction.
    if not np.array_equal(joined["inat_taxon_id"].to_numpy(int), species["inat_taxon_id"].to_numpy(int)):
        raise RuntimeError("reserve join changed taxon order")

    primary_mask = joined["second_fraction"].to_numpy(float) >= PRIMARY_THRESHOLD
    sens_mask = joined["second_fraction"].to_numpy(float) >= SENS_THRESHOLD
    complement_mask = ~primary_mask

    primary_results = {}
    sens_results = {}
    null_frames = []
    for j, metric in enumerate(METRICS):
        r, ndf = subset_metric(joined, null, primary_mask, j, "second_ge_0.10")
        primary_results[metric] = r
        null_frames.append(ndf)
        rs, ndfs = subset_metric(joined, null, sens_mask, j, "second_ge_0.20")
        sens_results[metric] = rs
        null_frames.append(ndfs)

    comp_primary, comp_null_df = subset_metric(joined, null, complement_mask, 0, "second_lt_0.10")
    null_frames.append(comp_null_df)
    pnull = null[primary_mask, 0, :].mean(axis=0)
    cnull = null[complement_mask, 0, :].mean(axis=0)
    delta_obs = float(primary_results["primary"]["observed_mean_rho"] - comp_primary["observed_mean_rho"])
    delta_null = pnull - cnull
    p_delta = float((1 + np.sum(delta_null >= delta_obs)) / (EXPECTED_NULL + 1))
    dilution = {
        "primary_minus_complement_observed_delta": delta_obs,
        "null_mean": float(delta_null.mean()),
        "null_q025": float(np.quantile(delta_null, 0.025)),
        "null_q975": float(np.quantile(delta_null, 0.975)),
        "p_upper": p_delta,
        "replicated_at_0_05": bool(delta_obs > 0 and p_delta < 0.05),
    }

    d = joined["D"].to_numpy(float)
    obs_primary_species = joined[OBS_COLS["primary"]].to_numpy(float)
    rho_obs = rho(d, obs_primary_species)
    rho_null = np.array([rho(d, null[:, 0, k]) for k in range(EXPECTED_NULL)], float)
    p_grad = float((1 + np.sum(rho_null >= rho_obs)) / (EXPECTED_NULL + 1))
    gradient = {
        "rho_D_species_primary_spatial_rho": rho_obs,
        "null_mean": float(rho_null.mean()),
        "null_q025": float(np.quantile(rho_null, 0.025)),
        "null_q975": float(np.quantile(rho_null, 0.975)),
        "p_upper": p_grad,
        "positive_supported_at_0_05": bool(rho_obs > 0 and p_grad < 0.05),
    }

    robust = all(primary_results[m]["positive_supported_at_0_05"] for m in ["primary", "observer_pair_exclusion", "calendar_quarter_stratification"])
    bg_support = bool(primary_results["matched_background_differential"]["positive_supported_at_0_05"])
    dilution_pass = bool(dilution["replicated_at_0_05"])
    if robust:
        verdict = "RESERVE_POLYMORPHIC_SUBSET_REPLICATED_ROBUSTLY"
        if dilution_pass:
            verdict += "_WITH_DILUTION_REPLICATION"
        if bg_support:
            verdict += "_WITH_MATCHED_BACKGROUND_SUPPORT"
    else:
        verdict = "RESERVE_POLYMORPHIC_SUBSET_NOT_ROBUSTLY_REPLICATED"

    joined["primary_second_ge_0_10"] = primary_mask
    joined["sensitivity_second_ge_0_20"] = sens_mask
    joined.to_csv(OUT / "reserve_species_membership_and_observed_metrics.csv", index=False)
    pd.concat(null_frames, ignore_index=True).to_csv(OUT / "reserve_subset_nulls.csv", index=False)
    pd.DataFrame({"permutation_index": np.arange(EXPECTED_NULL), "primary_minus_complement_delta": delta_null}).to_csv(OUT / "reserve_dilution_delta_null.csv", index=False)
    pd.DataFrame({"permutation_index": np.arange(EXPECTED_NULL), "rho_D_primary_spatial_rho": rho_null}).to_csv(OUT / "reserve_D_gradient_null.csv", index=False)

    result = {
        "analysis": "polymorphism_spatial_reserve_step5b",
        "date_jst": "2026-09-09",
        "protocol": "docs/POLYMORPHISM_SPATIAL_RESERVE_STEP5B_PROTOCOL_20260909.md",
        "upstream_run_id": 34178957447,
        "new_permutations_generated": False,
        "reserve_polymorphism_fingerprint": {
            "n_species": int(len(joined)),
            "D_min": float(joined["D"].min()),
            "D_max": float(joined["D"].max()),
            "fraction_second_ge_0_10": float(primary_mask.mean()),
            "fraction_second_ge_0_20": float(sens_mask.mean()),
        },
        "primary_subset": primary_results,
        "sensitivity_subset": sens_results,
        "primary_complement": comp_primary,
        "dilution_contrast": dilution,
        "continuous_D_diagnostic": gradient,
        "decision": {
            "robust_primary_observer_quarter_replication": bool(robust),
            "dilution_contrast_replicated": dilution_pass,
            "matched_background_subset_supported": bg_support,
            "continuous_gradient_supported_as_secondary": bool(gradient["positive_supported_at_0_05"]),
            "verdict": verdict,
            "common_boundary_claim_allowed": False,
            "adaptation_claim_allowed": False,
        },
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    p10 = primary_results
    p20 = sens_results
    lines = [
        "# Polymorphism spatial reserve Step 5b — result",
        "",
        f"**Frozen verdict: `{verdict}`.**",
        "",
        "All reserve subset statistics reuse the completed run 34178957447 species outputs and 999 null arrays; no spatial permutation was regenerated.",
        "",
        "## Primary reserve subset: second morph >=10%",
        "",
        f"- species: **{p10['primary']['n_species']}**",
        f"- primary mean rho **{p10['primary']['observed_mean_rho']:.9f}**, p **{p10['primary']['p_upper']:.6g}**",
        f"- observer-pair exclusion mean rho **{p10['observer_pair_exclusion']['observed_mean_rho']:.9f}**, p **{p10['observer_pair_exclusion']['p_upper']:.6g}**",
        f"- quarter-stratified mean rho **{p10['calendar_quarter_stratification']['observed_mean_rho']:.9f}**, p **{p10['calendar_quarter_stratification']['p_upper']:.6g}**",
        f"- matched-background differential mean rho **{p10['matched_background_differential']['observed_mean_rho']:.9f}**, p **{p10['matched_background_differential']['p_upper']:.6g}**",
        "",
        "## Sensitivity reserve subset: second morph >=20%",
        "",
        f"- species: **{p20['primary']['n_species']}**",
        f"- primary mean rho **{p20['primary']['observed_mean_rho']:.9f}**, p **{p20['primary']['p_upper']:.6g}**",
        f"- observer-pair exclusion p **{p20['observer_pair_exclusion']['p_upper']:.6g}**; quarter p **{p20['calendar_quarter_stratification']['p_upper']:.6g}**; matched-background p **{p20['matched_background_differential']['p_upper']:.6g}**",
        "",
        "## Replication of the dilution/gradient pattern",
        "",
        f"- <10% complement species: **{comp_primary['n_species']}**, primary mean rho **{comp_primary['observed_mean_rho']:.9f}**",
        f"- >=10% minus <10% delta: **{delta_obs:.9f}**, p_delta **{p_delta:.6g}**",
        f"- Spearman rho(D, species primary spatial rho): **{rho_obs:.6f}**, randomization p **{p_grad:.6g}**",
        "",
        "## Claim boundary",
        "",
        f"- primary + observer + quarter robust replication: **{robust}**",
        f"- dilution contrast replicated: **{dilution_pass}**",
        f"- matched-background subset support: **{bg_support}**",
        f"- continuous D gradient supported as secondary: **{gradient['positive_supported_at_0_05']}**",
        "- common-boundary / adaptation / causal environmental language: **not allowed**",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
