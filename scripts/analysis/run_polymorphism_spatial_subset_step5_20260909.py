#!/usr/bin/env python3
"""Aggregate frozen within-species spatial permutations over polymorphism subsets."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
OUT = ROOT / "results" / "polymorphism_spatial_subset_step5_20260909"
OUT.mkdir(parents=True, exist_ok=True)
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
EXPECTED_SPECIES = 369
EXPECTED_NULL = 999
PRIMARY_THRESHOLD = 0.10
SENS_THRESHOLD = 0.20
EXPECTED_WHOLE_MEAN = 0.02702130374565584


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-root", type=Path, required=True)
    return p.parse_args()


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def reconstruct_polymorphism() -> pd.DataFrame:
    df = pd.read_csv(DISCOVERY, usecols=["species", "morph", "global_classifiable"])
    keep = as_bool(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    rows = []
    for sp, g in df.loc[keep].groupby("species", sort=True):
        if len(g) < 40:
            continue
        cnt = Counter(g["morph"].astype(str))
        p = np.array([cnt.get(m, 0) / len(g) for m in MORPHS], float)
        ps = np.sort(p)[::-1]
        rows.append({
            "species": str(sp),
            "n_classifiable": int(len(g)),
            "D": float(1.0 - np.sum(p * p)),
            "second_fraction": float(ps[1]),
        })
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    if len(out) != EXPECTED_SPECIES:
        raise RuntimeError(f"discovery polymorphism frame mismatch: {len(out)} != {EXPECTED_SPECIES}")
    dmax = float(out["D"].max())
    s10 = float((out["second_fraction"] >= 0.10).mean())
    s20 = float((out["second_fraction"] >= 0.20).mean())
    if abs(dmax - 0.707645) > 1e-5:
        raise RuntimeError(f"D max fingerprint mismatch: {dmax}")
    if abs(s10 - 0.4661) > 0.001 or abs(s20 - 0.2656) > 0.001:
        raise RuntimeError(f"second-morph fingerprint mismatch: {s10}, {s20}")
    return out


def find_artifact_files(root: Path) -> tuple[list[Path], Path | None, Path | None]:
    shard = sorted(root.glob("global-rgfca-within-species-omnibus-shard-*/species_rho_permutations.csv"))
    full_species = next(iter(root.glob("global-rgfca-within-species-spatial-omnibus-full-v1/species_results.csv")), None)
    full_null = next(iter(root.glob("global-rgfca-within-species-spatial-omnibus-full-v1/global_null.csv")), None)
    if len(shard) != 20:
        raise RuntimeError(f"expected 20 shard permutation files, found {len(shard)}")
    return shard, full_species, full_null


def load_frozen_permutations(root: Path) -> pd.DataFrame:
    shard_files, full_species, full_null = find_artifact_files(root)
    pieces = [pd.read_csv(p, usecols=["species", "permutation_index", "rho"]) for p in shard_files]
    long = pd.concat(pieces, ignore_index=True)
    long["species"] = long["species"].astype(str)
    long["permutation_index"] = long["permutation_index"].astype(int)
    long["rho"] = long["rho"].astype(float)
    if long.duplicated(["species", "permutation_index"]).any():
        raise RuntimeError("duplicate species/permutation rows in frozen shards")
    species = sorted(long["species"].unique())
    if len(species) != EXPECTED_SPECIES:
        raise RuntimeError(f"frozen spatial species mismatch: {len(species)}")
    expected_indices = set([-1] + list(range(EXPECTED_NULL)))
    for sp, g in long.groupby("species", sort=False):
        idx = set(g["permutation_index"].tolist())
        if idx != expected_indices:
            raise RuntimeError(f"incomplete permutation index set for {sp}: {len(idx)}")

    # Reproduce the already published whole-frame result before subsetting.
    observed = long.loc[long["permutation_index"].eq(-1), "rho"].mean()
    if abs(float(observed) - EXPECTED_WHOLE_MEAN) > 1e-14:
        raise RuntimeError(f"whole-frame observed mean mismatch: {observed}")
    null = long.loc[long["permutation_index"].ge(0)].groupby("permutation_index")["rho"].mean().sort_index()
    if len(null) != EXPECTED_NULL:
        raise RuntimeError(f"whole-frame null count mismatch: {len(null)}")
    whole_p = float((1 + np.sum(null.to_numpy() >= observed)) / (EXPECTED_NULL + 1))
    if whole_p != 0.001:
        raise RuntimeError(f"whole-frame p mismatch: {whole_p}")

    if full_species is not None:
        saved = pd.read_csv(full_species, usecols=["species", "observed_rho"]).sort_values("species")
        obs = long.loc[long["permutation_index"].eq(-1), ["species", "rho"]].sort_values("species")
        merged = saved.merge(obs, on="species", validate="one_to_one")
        err = float(np.max(np.abs(merged["observed_rho"].to_numpy() - merged["rho"].to_numpy())))
        if err > 1e-14:
            raise RuntimeError(f"full artifact vs shard observed mismatch: {err}")
    if full_null is not None:
        saved_null = pd.read_csv(full_null).sort_values("permutation_index")
        err = float(np.max(np.abs(saved_null["mean_rho"].to_numpy() - null.to_numpy())))
        if err > 1e-14:
            raise RuntimeError(f"full artifact vs shard null mismatch: {err}")
    return long


def subset_summary(long: pd.DataFrame, members: set[str], label: str) -> tuple[dict, pd.DataFrame]:
    sub = long[long["species"].isin(members)]
    n = sub["species"].nunique()
    if n != len(members):
        raise RuntimeError(f"{label} membership mismatch {n} != {len(members)}")
    obs_rows = sub[sub["permutation_index"].eq(-1)]
    observed = float(obs_rows["rho"].mean())
    observed_median = float(obs_rows["rho"].median())
    observed_positive = float((obs_rows["rho"] > 0).mean())
    null = (
        sub[sub["permutation_index"].ge(0)]
        .groupby("permutation_index", sort=True)["rho"].mean()
        .rename("mean_rho").reset_index()
    )
    p = float((1 + np.sum(null["mean_rho"].to_numpy() >= observed)) / (EXPECTED_NULL + 1))
    result = {
        "label": label,
        "n_species": int(n),
        "observed_mean_rho": observed,
        "observed_median_species_rho": observed_median,
        "observed_positive_species_fraction": observed_positive,
        "null_mean": float(null["mean_rho"].mean()),
        "null_sd": float(null["mean_rho"].std(ddof=1)),
        "null_q025": float(null["mean_rho"].quantile(0.025)),
        "null_q975": float(null["mean_rho"].quantile(0.975)),
        "p_upper": p,
        "supported_at_0_05": bool(observed > 0 and p < 0.05),
    }
    null["subset"] = label
    return result, null


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    return float(stats.spearmanr(x, y).statistic)


def continuous_diagnostic(long: pd.DataFrame, poly: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    wide = long.pivot(index="species", columns="permutation_index", values="rho").sort_index()
    ptab = poly.set_index("species").loc[wide.index]
    d = ptab["D"].to_numpy(float)
    observed = spearman_rho(d, wide[-1].to_numpy(float))
    vals = []
    for idx in range(EXPECTED_NULL):
        vals.append(spearman_rho(d, wide[idx].to_numpy(float)))
    vals = np.asarray(vals, float)
    p = float((1 + np.sum(vals >= observed)) / (EXPECTED_NULL + 1))
    res = {
        "rho_D_species_spatial_rho": observed,
        "null_mean": float(vals.mean()),
        "null_q025": float(np.quantile(vals, 0.025)),
        "null_q975": float(np.quantile(vals, 0.975)),
        "p_upper": p,
        "positive_gradient_supported_at_0_05": bool(observed > 0 and p < 0.05),
    }
    return res, pd.DataFrame({"permutation_index": np.arange(EXPECTED_NULL), "rho_D_spatial_rho": vals})


def main() -> None:
    args = parse_args()
    poly = reconstruct_polymorphism()
    long = load_frozen_permutations(args.artifact_root)
    if set(poly["species"]) != set(long["species"]):
        missing_poly = sorted(set(long["species"]) - set(poly["species"]))
        missing_spatial = sorted(set(poly["species"]) - set(long["species"]))
        raise RuntimeError(f"species identity mismatch: poly_missing={missing_poly[:10]}, spatial_missing={missing_spatial[:10]}")

    primary_members = set(poly.loc[poly["second_fraction"].ge(PRIMARY_THRESHOLD), "species"])
    sens_members = set(poly.loc[poly["second_fraction"].ge(SENS_THRESHOLD), "species"])
    complement = set(poly["species"]) - primary_members

    primary, null_primary = subset_summary(long, primary_members, "second_ge_0.10")
    sensitivity, null_sens = subset_summary(long, sens_members, "second_ge_0.20")
    comp, null_comp = subset_summary(long, complement, "second_lt_0.10")

    pn = null_primary.set_index("permutation_index")["mean_rho"]
    cn = null_comp.set_index("permutation_index")["mean_rho"]
    delta_obs = float(primary["observed_mean_rho"] - comp["observed_mean_rho"])
    delta_null = (pn - cn).sort_index()
    p_delta = float((1 + np.sum(delta_null.to_numpy() >= delta_obs)) / (EXPECTED_NULL + 1))
    dilution = {
        "primary_minus_complement_observed_delta": delta_obs,
        "null_mean": float(delta_null.mean()),
        "null_q025": float(delta_null.quantile(0.025)),
        "null_q975": float(delta_null.quantile(0.975)),
        "p_upper": p_delta,
        "dilution_contrast_supported_at_0_05": bool(delta_obs > 0 and p_delta < 0.05),
    }

    continuous, continuous_null = continuous_diagnostic(long, poly)

    primary_pass = bool(primary["supported_at_0_05"])
    if primary_pass and dilution["dilution_contrast_supported_at_0_05"]:
        verdict = "POLYMORPHIC_SUBSET_SPATIALLY_STRUCTURED_WITH_DILUTION_CONTRAST"
    elif primary_pass:
        verdict = "POLYMORPHIC_SUBSET_SPATIALLY_STRUCTURED_NO_DILUTION_CONTRAST"
    else:
        verdict = "POLYMORPHIC_SUBSET_NOT_SUPPORTED"

    membership = poly.copy()
    membership["primary_second_ge_0_10"] = membership["second_fraction"].ge(PRIMARY_THRESHOLD)
    membership["sensitivity_second_ge_0_20"] = membership["second_fraction"].ge(SENS_THRESHOLD)
    observed_species = long.loc[long["permutation_index"].eq(-1), ["species", "rho"]].rename(columns={"rho": "spatial_observed_rho"})
    membership = membership.merge(observed_species, on="species", validate="one_to_one")

    null_all = pd.concat([null_primary, null_sens, null_comp], ignore_index=True)
    delta_df = delta_null.rename("primary_minus_complement_delta").reset_index()

    result = {
        "analysis": "polymorphism_spatial_subset_step5",
        "date_jst": "2026-09-09",
        "inferential_role": "postoutcome_exploratory_subset_of_existing_exact_randomization",
        "protocol": "docs/POLYMORPHISM_SPATIAL_SUBSET_STEP5_PROTOCOL_20260909.md",
        "upstream": {
            "github_run_id": 34088925008,
            "spatial_protocol": "global-rgfca-within-species-spatial-omnibus-v1",
            "whole_frame_n_species": 369,
            "whole_frame_observed_mean_rho": EXPECTED_WHOLE_MEAN,
            "whole_frame_p_upper": 0.001,
            "permutations_reused": 999,
            "new_spatial_permutations_generated": False,
        },
        "polymorphism_fingerprint": {
            "n_species": int(len(poly)),
            "D_min": float(poly["D"].min()),
            "D_max": float(poly["D"].max()),
            "fraction_second_ge_0_10": float((poly["second_fraction"] >= 0.10).mean()),
            "fraction_second_ge_0_20": float((poly["second_fraction"] >= 0.20).mean()),
        },
        "primary": primary,
        "sensitivity": sensitivity,
        "complement": comp,
        "dilution_contrast": dilution,
        "continuous_D_diagnostic": continuous,
        "decision": {
            "verdict": verdict,
            "claim3_subset_strengthened": primary_pass,
            "dilution_language_allowed": bool(dilution["dilution_contrast_supported_at_0_05"]),
            "continuous_gradient_language_allowed_as_exploratory": bool(continuous["positive_gradient_supported_at_0_05"]),
            "common_boundary_claim_allowed": False,
            "adaptation_claim_allowed": False,
        },
    }

    membership.to_csv(OUT / "species_membership_and_observed_rho.csv", index=False)
    null_all.to_csv(OUT / "subset_global_nulls.csv", index=False)
    delta_df.to_csv(OUT / "primary_vs_complement_delta_null.csv", index=False)
    continuous_null.to_csv(OUT / "D_spatial_rho_null.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Polymorphism spatial subset Step 5 — result",
        "",
        f"**Frozen verdict: `{verdict}`.**",
        "",
        "This is a post-outcome exploratory subset analysis of the already completed 369-species exact randomization; all 999 upstream per-species permutations were reused without regeneration.",
        "",
        "## Primary: second morph >=10%",
        "",
        f"- species: **{primary['n_species']}**",
        f"- observed equal-species mean rho: **{primary['observed_mean_rho']:.9f}**",
        f"- null 2.5–97.5%: **[{primary['null_q025']:.9f}, {primary['null_q975']:.9f}]**",
        f"- upper-tail p: **{primary['p_upper']:.6g}**",
        "",
        "## Sensitivity: second morph >=20%",
        "",
        f"- species: **{sensitivity['n_species']}**",
        f"- observed equal-species mean rho: **{sensitivity['observed_mean_rho']:.9f}**",
        f"- null 2.5–97.5%: **[{sensitivity['null_q025']:.9f}, {sensitivity['null_q975']:.9f}]**",
        f"- upper-tail p: **{sensitivity['p_upper']:.6g}**",
        "",
        "## Primary subset versus complement",
        "",
        f"- complement species: **{comp['n_species']}**",
        f"- complement observed mean rho: **{comp['observed_mean_rho']:.9f}**",
        f"- observed delta: **{delta_obs:.9f}**",
        f"- delta null 2.5–97.5%: **[{dilution['null_q025']:.9f}, {dilution['null_q975']:.9f}]**",
        f"- one-sided p_delta: **{p_delta:.6g}**",
        "",
        "## Continuous D diagnostic",
        "",
        f"- Spearman rho(D, species spatial rho): **{continuous['rho_D_species_spatial_rho']:.6f}**",
        f"- randomization p: **{continuous['p_upper']:.6g}**",
        "",
        "## Claim boundary",
        "",
        f"- Claim 3 sharpened to the >=10% polymorphic subset: **{primary_pass}**",
        f"- explicit dilution interpretation supported: **{dilution['dilution_contrast_supported_at_0_05']}**",
        f"- positive D–spatial-structure gradient supported as exploratory: **{continuous['positive_gradient_supported_at_0_05']}**",
        "- common boundary / adaptation / causal environmental claim: **not allowed**",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
