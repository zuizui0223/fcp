"""Exploratory W/D/geographic-depletion associations with frozen species-level environment contrasts.

This does not estimate new climate effects from photos and cannot establish any
mechanism. All 18 possible association tests are reported with familywise Holm
adjustment rather than selecting nominally significant results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ENV_NAMES = ("bio5", "bio14", "srad_mean")
TRAITS = ("D", "w_alignment_sq", "depletion_50km")
MODES = ("signed", "absolute")
FROZEN_ENV_SOURCE_SHA256 = "631d450fb31f224dc8d50fd86706394e9d8b741d51dfdbceefa1afacb753a81f"
ENV_SOURCE_ARTIFACT_ID = 10860659652


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_frozen_environment(source: pd.DataFrame, reference: dict) -> None:
    needed = {"species", "variable", "delta_white_minus_nonwhite_SD"}
    if not needed.issubset(source):
        raise ValueError("source missing frozen species/variable/white delta columns")
    if source.duplicated(["species", "variable"]).any() or source["species"].isna().any():
        raise ValueError("duplicate species×environment or missing source species")
    if set(source["variable"].dropna()) != set(ENV_NAMES):
        raise ValueError("incorrect environmental variable family")
    primary = {r["variable"]: r for r in reference["results"]}
    for var in ENV_NAMES:
        row = source.loc[source["variable"].eq(var)]
        if len(row) != 281 or len(row) != int(primary[var]["n_species_delta"]):
            raise ValueError(f"{var}: frozen species-count mismatch")
        a = pd.to_numeric(row["delta_white_minus_nonwhite_SD"], errors="raise").to_numpy(float)
        if not np.isfinite(a).all():
            raise ValueError(f"{var}: nonfinite white environmental delta")
        for target, observed in (("mean_delta_white_minus_nonwhite_SD", float(a.mean())),
                                 ("median_delta_white_minus_nonwhite_SD", float(np.median(a)))):
            if not np.isclose(observed, primary[var][target], rtol=0, atol=1e-12):
                raise ValueError(f"{var}: source mismatch with frozen 2026-09-25 environmental result")


def holm_adjust(p: np.ndarray) -> np.ndarray:
    """Holm over one fixed 18-test family; missing tests count as p=1."""
    p = np.asarray(p, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = len(p)
    order = np.argsort(p, kind="stable")
    sorted_adjusted = np.minimum(1.0, np.maximum.accumulate((n - np.arange(n)) * p[order]))
    adjusted = np.empty_like(sorted_adjusted)
    adjusted[order] = sorted_adjusted
    return adjusted


def compare_species(registry: pd.DataFrame, env: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = {"inat_taxon_id", "species", "cohort", "D", "w_alignment_sq", "depletion_50km"}
    if not columns.issubset(registry):
        raise ValueError("U100 registry missing required D/W/spatial columns")
    if registry["species"].duplicated().any() or registry["inat_taxon_id"].duplicated().any():
        raise ValueError("U100 registry contains duplicated species")
    if len(registry) != 4730:
        raise ValueError("U100 registry size must be exactly 4730")
    env_wide = env.pivot(index="species", columns="variable", values="delta_white_minus_nonwhite_SD")
    if not env_wide.index.isin(registry["species"]).all():
        raise ValueError("frozen environmental species absent from U100 frame")
    merged = registry.merge(env_wide.reset_index(), on="species", how="left", validate="one_to_one")
    for var in ENV_NAMES:
        named = merged.loc[merged[var].notna(), "cohort"].unique()
        if set(named) != {"prospective"}:
            raise ValueError(f"{var}: environmental contrast leaks to non-prospective cohort")
    records = []
    for env_name in ENV_NAMES:
        for trait in TRAITS:
            for mode in MODES:
                a = merged.loc[merged["cohort"].eq("prospective"), [trait, env_name]].copy()
                a[trait] = pd.to_numeric(a[trait], errors="coerce")
                a[env_name] = pd.to_numeric(a[env_name], errors="coerce")
                a = a.replace([np.inf, -np.inf], np.nan).dropna()
                x, y = a[trait].to_numpy(float), a[env_name].to_numpy(float)
                if mode == "absolute":
                    y = np.abs(y)
                estimable = len(a) >= 30 and np.ptp(x) > 0 and np.ptp(y) > 0
                rho, p = spearmanr(x, y) if estimable else (np.nan, np.nan)
                records.append({"environment": env_name, "trait": trait,
                                "contrast": mode, "n_species": len(a),
                                "rho": float(rho), "p_nominal_exploratory": float(p),
                                "status": "exploratory" if estimable else "not_estimable"})
    result = pd.DataFrame(records)
    assert len(result) == len(ENV_NAMES) * len(TRAITS) * len(MODES) == 18
    result["holm_18_p"] = holm_adjust(result["p_nominal_exploratory"].to_numpy(float))
    result["holm_18_supported"] = result["holm_18_p"].lt(.05) & result["status"].eq("exploratory")
    return result, merged.loc[merged["cohort"].eq("prospective"),
                              ["inat_taxon_id", "species", "D", "w_alignment_sq", "depletion_50km", *ENV_NAMES]].copy()


def verified_run(source_file: Path, reference_file: Path, registry_file: Path, outdir: Path) -> dict:
    if file_sha256(source_file) != FROZEN_ENV_SOURCE_SHA256:
        raise ValueError("original frozen environmental species-delta SHA256 mismatch")
    env = pd.read_csv(source_file)
    validate_frozen_environment(env, json.loads(reference_file.read_text(encoding="utf-8")))
    comparisons, joined = compare_species(pd.read_csv(registry_file, low_memory=False), env)
    outdir.mkdir(parents=True, exist_ok=True)
    comparison_path = outdir / "frozen_species_environment_18_comparisons.csv"
    joined_path = outdir / "prospective_frozen_environment_join.csv"
    comparisons.to_csv(comparison_path, index=False, float_format="%.17g")
    joined.to_csv(joined_path, index=False, float_format="%.17g")
    receipt = {"status": "posthoc_18_associations_not_environmental_mechanism_confirmation",
               "source_artifact_id": ENV_SOURCE_ARTIFACT_ID,
               "source_sha256": FROZEN_ENV_SOURCE_SHA256,
               "reference_sha256": file_sha256(reference_file),
               "registry_sha256": file_sha256(registry_file),
               "num_environment_variables": 3, "num_traits": 3, "num_modes": 2,
               "prespecified_as_confirmatory": False,
               "num_comparisons": len(comparisons),
               "min_holm_18_p": float(comparisons["holm_18_p"].min()),
               "holm_18_supported_count": int(comparisons["holm_18_supported"].sum()),
               "source_environment_species": 281,
               "original_observer_sensitivity_still_applies": True,
               "original_H2_or_spatial_decisions_changed": False,
               "results_sha256": file_sha256(comparison_path),
               "join_sha256": file_sha256(joined_path)}
    (outdir / "frozen_environment_screen_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--registry", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args(argv)
    print(json.dumps(verified_run(a.source, a.reference, a.registry, a.outdir)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
