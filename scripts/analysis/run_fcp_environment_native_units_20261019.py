#!/usr/bin/env python3
"""Exact common-unit re-expression of the frozen joint FCP environment fits.

Section A of FCP_ENVIRONMENT_TRANSPORT_PROTOCOL_20261019.md. This does not
refit the photograph models, trim responses, revise traits, or overwrite the
standardized analysis. Equal-species fitted-slope means are the primary target.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.special import expit

from scripts.analysis.run_fcp_447500_environment_model_20261018 import (
    PREDICTORS, load_panel, reml_meta,
)
from scripts.analysis.run_fcp_gift_modifiers_20261018 import (
    CONTRASTS, FAMILY_SIZE, summarize_contrast,
)


PHYSICAL_INCREMENTS = {"BIO5": 1., "BIO14": 10., "SRADmean": 1000.}
UNIT_LABELS = {"BIO5": "1 degree C", "BIO14": "10 mm",
               "SRADmean": "1000 kJ m^-2 day^-1"}
REPETITIONS = 5000
MEAN_BOOTSTRAP_SEED = 20261019
TRAIT_BASE_SEED = 20261018  # Identical draws to the frozen SD-scale contrasts.
LINEAR_META_KEYS = (
    "reml_mean", "reml_mean_se_mhksj", "reml_mean_ci_low", "reml_mean_ci_high",
    "tau", "working_new_species_prediction_low", "working_new_species_prediction_high",
    "equal_species_unshrunk_mean", "equal_species_unshrunk_se_independence",
    "equal_species_unshrunk_ci_low", "equal_species_unshrunk_ci_high",
    "equal_species_eb_mean_conditional",
)
SQUARED_META_KEYS = ("tau2", "tau2_profile_ci_low", "tau2_profile_ci_high")


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def verify_file(path: Path, expected: str, provenance: str) -> dict:
    path = Path(path).resolve()
    actual = sha(path)
    require(actual == expected, f"SHA256 mismatch: {path}")
    return {"path": str(path), "receipt": provenance, "expected_sha256": expected,
            "actual_sha256": actual, "bytes": path.stat().st_size, "verified": True}


def authenticate_inputs(root: Path, inputs: Path) -> tuple[list[dict], dict]:
    """Verify every file/code/cache digest named by the three source receipts."""
    environment = root / "results/fcp_real_447500_environment_model_20261018"
    gift = root / "results/fcp_gift_modifiers_20261018"
    env_receipt_path = environment / "execution_receipt.json"
    gift_receipt_path = gift / "profile.json"
    modifier_receipt_path = gift / "modifier_profile.json"
    env = json.loads(env_receipt_path.read_text())
    extraction = json.loads(gift_receipt_path.read_text())
    modifier = json.loads(modifier_receipt_path.read_text())
    records = []

    def check(path, digest, provenance):
        records.append(verify_file(path, digest, str(provenance)))

    for folder, receipt, receipt_path in [
        (environment, env, env_receipt_path), (gift, extraction, gift_receipt_path),
        (gift, modifier, modifier_receipt_path),
    ]:
        for name, digest in receipt["output_sha256"].items():
            check(folder / name, digest, receipt_path)
        for name, digest in receipt.get("source_sha256", {}).items():
            check(root / name, digest, receipt_path)
    for relative, digest, receipt_path in [
        ("scripts/analysis/run_fcp_447500_environment_model_20261018.py", env["script_sha256"], env_receipt_path),
        ("scripts/analysis/summarize_fcp_447500_within_species_environment_20261018.py", env["imported_sampling_script_sha256"], env_receipt_path),
        ("scripts/data/fetch_fcp_gift_modifiers_20261018.py", extraction["code_sha256"], gift_receipt_path),
        ("scripts/analysis/run_fcp_gift_modifiers_20261018.py", modifier["code_sha256"], modifier_receipt_path),
    ]:
        check(root / relative, digest, receipt_path)
    photos = inputs / "full447500/measured/fcp_all4730_original_method_measured_photos.csv.gz"
    check(photos, env["photo_sha256"], env_receipt_path)
    worldclim = inputs / "worldclim_frozen"
    for name, digest in env["worldclim_zip_sha256"].items():
        check(worldclim / name, digest, env_receipt_path)
    for raster in env["rasters"]:
        subdir = "bio" if "_bio_" in raster["file"] else "srad"
        check(worldclim / subdir / raster["file"], raster["sha256"], env_receipt_path)
    for page in extraction["public_api_pages"]:
        check(inputs / "gift_3_2_20261018" / page["cache_file"], page["sha256"], gift_receipt_path)
    require(modifier["bootstrap_repetitions"] == REPETITIONS, "Frozen bootstrap count changed")
    require(modifier["contrasts"] == {k: list(v) for k, v in CONTRASTS.items()}, "Frozen trait levels changed")
    require(modifier["family_size"] == FAMILY_SIZE == 9, "Frozen contrast family changed")
    return records, {
        "environment": environment, "gift": gift, "photos": photos, "worldclim": worldclim,
        "source_receipt_sha256": {str(p): sha(p) for p in [env_receipt_path, gift_receipt_path, modifier_receipt_path]},
    }


def convert_primary_slopes(primary: pd.DataFrame) -> pd.DataFrame:
    """Apply b_j*d_j/s_j and V_jk*d_j*d_k/(s_j*s_k), retaining SD originals."""
    q = primary.copy()
    require(not q.empty, "No primary slopes")
    require(q.model.eq("primary").all() and q.status.eq("ok").all(), "Only successful primary slopes are admitted")
    require(not q.duplicated(["inat_taxon_id", "variable"]).any(), "Duplicate species/axis")
    require(set(q.variable) == set(PREDICTORS), "Wrong predictor set")
    groups = q.groupby("inat_taxon_id")
    require(groups.size().eq(3).all(), "Each species must have all three joint slopes")
    require(groups.species.nunique().eq(1).all(), "Species identity varies between axes")
    s = q.primary_within_species_sd_native.to_numpy(float)
    require((np.isfinite(s) & (s > 0)).all(), "Environmental SD must be finite and positive")
    require(np.isfinite(q[["slope", "variance_approx", "se_approx"]]).all().all(), "Nonfinite source coefficient")
    require(q.variance_approx.gt(0).all(), "Nonpositive source sampling variance")
    np.testing.assert_allclose(q.se_approx**2, q.variance_approx, rtol=2e-10, atol=1e-12)
    for tid, g in groups:
        g = g.set_index("variable").loc[list(PREDICTORS)]
        covariance = g[["cov_with_" + v for v in PREDICTORS]].to_numpy(float)
        require(np.isfinite(covariance).all(), f"Nonfinite covariance: {tid}")
        require(np.allclose(covariance, covariance.T, rtol=2e-10, atol=1e-12), f"Asymmetric covariance: {tid}")
        require((np.linalg.eigvalsh(covariance) > 0).all(), f"Nonpositive covariance: {tid}")
        require(np.allclose(np.diag(covariance), g.variance_approx, rtol=2e-10, atol=1e-12), f"Covariance diagonal differs: {tid}")
    q["physical_increment"] = q.variable.map(PHYSICAL_INCREMENTS)
    q["physical_increment_label"] = q.variable.map(UNIT_LABELS)
    q["conversion_factor"] = q.physical_increment / q.primary_within_species_sd_native
    for column in ["slope", "se_approx", "variance_approx", "ci_low_wald", "ci_high_wald"]:
        q[column + "_standardized"] = q[column]
        power = 2 if column == "variance_approx" else 1
        q[column] = q[column] * q.conversion_factor**power
    factors = q.pivot(index="inat_taxon_id", columns="variable", values="conversion_factor")
    for variable in PREDICTORS:
        name = "cov_with_" + variable
        q["cov_standardized_with_" + variable] = q[name]
        q[name] = q[name] * q.conversion_factor * q.inat_taxon_id.map(factors[variable])
    np.testing.assert_allclose(q.slope / q.physical_increment, q.slope_per_native_unit, rtol=2e-10, atol=1e-12)
    q["native_sd_divided_by_increment"] = q.primary_within_species_sd_native / q.physical_increment
    q["observed_span_native"] = q.primary_predictor_max - q.primary_predictor_min
    q["observed_span_divided_by_increment"] = q.observed_span_native / q.physical_increment
    q["increment_exceeds_observed_span"] = q.observed_span_divided_by_increment.lt(1)
    q["native_sd_below_0p1_increment"] = q.native_sd_divided_by_increment.lt(.1)
    threshold = 32 * np.finfo(float).eps * np.maximum(1., np.maximum(abs(q.primary_predictor_min), abs(q.primary_predictor_max)))
    q["near_zero_sd_numerical_threshold"] = threshold
    q["near_zero_sd_numerical"] = q.primary_within_species_sd_native.le(threshold)
    require(~q.near_zero_sd_numerical.any(), "Source primary numerical SD gate no longer holds")
    q["axis_order"] = q.variable.map({v: j for j, v in enumerate(PREDICTORS)})
    return q.sort_values(["inat_taxon_id", "axis_order"]).drop(columns="axis_order").reset_index(drop=True)


def reml_in_physical_units(y, variance, numerical_scale=None):
    """Use unchanged source REML after a common numerical rescaling, then undo it.

    A common (not species-specific) scale does not change the physical estimand
    or relative weights. It avoids absolute optimizer tolerances at tiny tau2.
    """
    y, variance = np.asarray(y, float), np.asarray(variance, float)
    scale = float(np.sqrt(np.median(variance))) if numerical_scale is None else float(numerical_scale)
    require(np.isfinite(scale) and scale > 0, "Invalid common REML numerical scale")
    summary, means, variances, probabilities = reml_meta(y / scale, variance / scale**2)
    for key in LINEAR_META_KEYS:
        summary[key] *= scale
    for key in SQUARED_META_KEYS:
        summary[key] *= scale**2
    summary["reml_common_numerical_scale"] = scale
    return summary, means * scale, variances * scale**2, probabilities


def summarize_native(native):
    summaries, checks = [], []
    cohorts = ["pooled", *sorted(native.source_cohort.unique())]
    for cohort in cohorts:
        q = native if cohort == "pooled" else native.loc[native.source_cohort.eq(cohort)]
        for variable in PREDICTORS:
            g = q.loc[q.variable.eq(variable)]
            result, *_ = reml_in_physical_units(g.slope, g.variance_approx)
            base = {"model": "primary", "group": cohort, "variable": variable,
                    "physical_increment": PHYSICAL_INCREMENTS[variable], "physical_increment_label": UNIT_LABELS[variable]}
            summaries.append({**base, **result})
            common_scale = result["reml_common_numerical_scale"]
            for tested_scale in [1., common_scale / 10, common_scale * 10]:
                other, *_ = reml_in_physical_units(g.slope, g.variance_approx, tested_scale)
                for name in ["reml_mean", "reml_mean_se_mhksj", *SQUARED_META_KEYS]:
                    power = 2 if name in SQUARED_META_KEYS else 1
                    difference = abs(other[name] - result[name])
                    tolerance = 1e-4 * abs(result[name]) + 1e-8 * common_scale**power
                    passed = difference <= tolerance
                    checks.append({**base, "statistic": name, "reference_numerical_scale": common_scale,
                                   "tested_numerical_scale": tested_scale, "reference": result[name],
                                   "test": other[name], "absolute_difference": difference,
                                   "absolute_tolerance": tolerance, "passed": passed})
                    require(passed, f"REML numerical-scale check failed: {cohort}/{variable}/{name}")
    return pd.DataFrame(summaries), pd.DataFrame(checks)


def means_from_multiplicities(multiplicities, species_counts, slope_sums):
    denominator = np.asarray(multiplicities) @ np.asarray(species_counts)
    require((denominator > 0).all(), "Empty genus bootstrap mean")
    return (np.asarray(multiplicities) @ np.asarray(slope_sums)) / denominator[:, None]


def mean_bootstrap(wide, repetitions=REPETITIONS, seed=MEAN_BOOTSTRAP_SEED):
    q = wide.assign(genus=wide.species.str.split().str[0])
    aggregate = q.groupby("genus", sort=True)
    counts = aggregate.size().to_numpy()
    sums = aggregate[list(PREDICTORS)].sum().to_numpy()
    n_genera = len(counts)
    rng = np.random.default_rng(seed)
    multiplicities = rng.multinomial(n_genera, np.full(n_genera, 1 / n_genera), size=repetitions)
    values = means_from_multiplicities(multiplicities, counts, sums)
    rows = []
    for j, variable in enumerate(PREDICTORS):
        low, high = np.quantile(values[:, j], [.025, .975])
        rows.append({"variable": variable, "n_species": len(q), "n_genera": n_genera,
                     "repetitions": repetitions, "seed": seed, "equal_species_mean": float(q[variable].mean()),
                     "genus_cluster_bootstrap_sd": float(values[:, j].std(ddof=1)),
                     "genus_cluster_percentile_low": float(low), "genus_cluster_percentile_high": float(high),
                     "shared_draws_across_axes": True})
    return pd.DataFrame(rows), values


def observed_prediction_and_covariance_audit(old, native, panel):
    """Check all authentic photo predictors without refitting any photo model."""
    old_groups = {int(tid): g.set_index("variable").loc[list(PREDICTORS)]
                  for tid, g in old.groupby("inat_taxon_id")}
    new_groups = {int(tid): g.set_index("variable").loc[list(PREDICTORS)]
                  for tid, g in native.groupby("inat_taxon_id")}
    increments = np.array([PHYSICAL_INCREMENTS[v] for v in PREDICTORS])
    rows = []
    for tid, photos in panel.groupby("inat_taxon_id", sort=True):
        a, b = old_groups[int(tid)], new_groups[int(tid)]
        require(len(photos) == int(a.n_primary.iloc[0]), "Primary photo count differs")
        require(photos.species.eq(a.species.iloc[0]).all(), "Photo species identity differs")
        raw = photos[list(PREDICTORS)].to_numpy(float)
        sd = a.primary_within_species_sd_native.to_numpy()
        np.testing.assert_allclose(raw.std(axis=0), sd, rtol=2e-10, atol=1e-10)
        np.testing.assert_allclose(raw.min(axis=0), a.primary_predictor_min, rtol=2e-10, atol=1e-10)
        np.testing.assert_allclose(raw.max(axis=0), a.primary_predictor_max, rtol=2e-10, atol=1e-10)
        centered = raw - raw.mean(axis=0)
        old_x, new_x = centered / sd, centered / increments
        old_eta = a.intercept.iloc[0] + old_x @ a.slope.to_numpy()
        new_eta = b.intercept.iloc[0] + new_x @ b.slope.to_numpy()
        old_cov = a[["cov_with_" + v for v in PREDICTORS]].to_numpy()
        new_cov = b[["cov_with_" + v for v in PREDICTORS]].to_numpy()
        factor = increments / sd
        old_corr = old_cov / np.sqrt(np.outer(np.diag(old_cov), np.diag(old_cov)))
        new_corr = new_cov / np.sqrt(np.outer(np.diag(new_cov), np.diag(new_cov)))
        old_score_variance = np.einsum("ij,jk,ik->i", old_x, old_cov, old_x)
        new_score_variance = np.einsum("ij,jk,ik->i", new_x, new_cov, new_x)
        np.testing.assert_allclose(old_eta, new_eta, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(old_score_variance, new_score_variance, rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(old_corr, new_corr, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(np.diag(new_cov), b.variance_approx, rtol=1e-10, atol=1e-10)
        require(np.linalg.eigvalsh(new_cov).min() > 0, "Converted covariance lost positive definiteness")
        rows.append({"inat_taxon_id": int(tid), "species": a.species.iloc[0], "n_photos": len(photos),
                     "max_abs_linear_predictor_difference": float(np.max(abs(old_eta - new_eta))),
                     "max_abs_fitted_probability_difference": float(np.max(abs(expit(old_eta) - expit(new_eta)))),
                     "max_abs_slope_score_variance_difference": float(np.max(abs(old_score_variance - new_score_variance))),
                     "max_relative_slope_score_variance_difference": float(np.max(abs(old_score_variance - new_score_variance) / np.maximum(1., abs(old_score_variance)))),
                     "max_abs_covariance_roundtrip_difference": float(np.max(abs(new_cov / np.outer(factor, factor) - old_cov))),
                     "max_abs_correlation_difference": float(np.max(abs(new_corr - old_corr))),
                     "max_abs_wald_z_difference": float(np.max(abs(a.slope / a.se_approx - b.slope / b.se_approx))),
                     "old_covariance_min_eigenvalue": float(np.linalg.eigvalsh(old_cov).min()),
                     "converted_covariance_min_eigenvalue": float(np.linalg.eigvalsh(new_cov).min()),
                     "passed": True})
    require(len(rows) == 1574 and sum(r["n_photos"] for r in rows) == 86114, "Prediction audit panel differs")
    return pd.DataFrame(rows)


def range_diagnostics(native):
    summaries, contributions, flagged = [], [], []
    for variable in PREDICTORS:
        q = native.loc[native.variable.eq(variable)].copy()
        d = PHYSICAL_INCREMENTS[variable]
        s = q.primary_within_species_sd_native.to_numpy()
        y, v = q.slope.to_numpy(), q.variance_approx.to_numpy()
        span = q.observed_span_native.to_numpy()
        summary = {"variable": variable, "physical_increment": d, "physical_increment_label": UNIT_LABELS[variable],
                   "n_species": len(q), "n_numerical_near_zero_sd": int(q.near_zero_sd_numerical.sum()),
                   "n_native_sd_below_0p1_increment": int(np.sum(s < .1 * d)),
                   "n_increment_exceeds_observed_span": int(np.sum(span < d)),
                   "n_species_trimmed": 0}
        for name, values in [("sd_native", s), ("span_native", span), ("slope", y), ("se_approx", np.sqrt(v))]:
            for label, probability in [("min", 0), ("p01", .01), ("p05", .05), ("median", .5), ("p95", .95), ("p99", .99), ("max", 1)]:
                summary[name + "_" + label] = float(np.quantile(values, probability))
        summaries.append(summary)
        q["lowest_1pct_sd"] = s <= np.quantile(s, .01)
        q["lowest_5pct_sd"] = s <= np.quantile(s, .05)
        for label, mask in [("lowest_1pct_sd", q.lowest_1pct_sd.to_numpy()),
                            ("lowest_5pct_sd", q.lowest_5pct_sd.to_numpy()),
                            ("sd_below_0p1_increment", s < .1 * d)]:
            contributions.append({"variable": variable, "tail_definition": label, "n_tail_species": int(mask.sum()),
                                  "n_full_panel_species": len(q), "n_trimmed": 0,
                                  "contribution_to_full_equal_species_mean": float(y[mask].sum() / len(q)),
                                  "share_absolute_slope_sum": float(abs(y[mask]).sum() / abs(y).sum()),
                                  "share_independence_sampling_variance_sum": float(v[mask].sum() / v.sum()),
                                  "full_equal_species_mean": float(y.mean())})
        flagged.append(q.loc[q.lowest_5pct_sd | q.native_sd_below_0p1_increment].sort_values("primary_within_species_sd_native"))
    return pd.DataFrame(summaries), pd.DataFrame(contributions), pd.concat(flagged, ignore_index=True)


def trait_comparisons(old_wide, new_wide, frozen_contrasts):
    standardized, converted = [], []
    for i, trait in enumerate(CONTRASTS):
        standardized.extend(summarize_contrast(old_wide, trait, REPETITIONS, TRAIT_BASE_SEED + i))
        converted.extend(summarize_contrast(new_wide, trait, REPETITIONS, TRAIT_BASE_SEED + i))
    old, new = pd.DataFrame(standardized), pd.DataFrame(converted)
    keys = ["trait", "variable"]
    actual = old.set_index(keys).sort_index()
    frozen = frozen_contrasts.set_index(keys).sort_index()
    require(actual.index.equals(frozen.index), "Frozen trait comparisons differ")
    for column in frozen.columns:
        if pd.api.types.is_numeric_dtype(frozen[column]):
            np.testing.assert_allclose(actual[column], frozen[column], rtol=1e-11, atol=1e-12, equal_nan=True)
        else:
            require(actual[column].fillna("").equals(frozen[column].fillna("")), "Frozen trait field differs: " + column)
    coverage = ["n_first", "n_reference", "genera_first", "genera_reference", "n_excluded_missing_or_other_level", "genus_clusters_resampled", "invalid_bootstrap_draws"]
    for column in coverage:
        require(old[column].equals(new[column]), "Native conversion changed GIFT coverage: " + column)
    new["physical_increment"] = new.variable.map(PHYSICAL_INCREMENTS)
    new["physical_increment_label"] = new.variable.map(UNIT_LABELS)
    new["nominal_interval_excludes_zero"] = (new.percentile_ci_low > 0) | (new.percentile_ci_high < 0)
    new["family9_interval_excludes_zero"] = (new.family9_percentile_ci_low > 0) | (new.family9_percentile_ci_high < 0)
    comparison = new.merge(old[keys + ["difference_first_minus_reference", "percentile_ci_low", "percentile_ci_high", "family9_percentile_ci_low", "family9_percentile_ci_high"]].rename(columns=lambda c: c if c in keys else "sd_scale_" + c), on=keys, validate="one_to_one")
    return new, comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs-root", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("results/fcp_environment_native_units_20261019"))
    args = parser.parse_args()
    started = datetime.now(timezone.utc).isoformat()
    root = Path(__file__).resolve().parents[2]
    inputs = args.inputs_root.resolve() if args.inputs_root else root.parent / "fcp_inputs"
    outdir = (root / args.outdir).resolve()
    require(outdir == root / "results/fcp_environment_native_units_20261019", "Use only the assigned new results directory")
    outdir.mkdir(parents=True, exist_ok=True)
    protocol = root / "docs/FCP_ENVIRONMENT_TRANSPORT_PROTOCOL_20261019.md"
    protocol_text = protocol.read_text()
    (outdir / "protocol_at_execution.md").write_text(protocol_text)
    authentication, paths = authenticate_inputs(root, inputs)
    print(f"AUTHENTICATED {len(authentication)} frozen file/code hash references", flush=True)
    source = pd.read_csv(paths["environment"] / "species_slopes.csv")
    primary = source.loc[source.model.eq("primary")].copy()
    require(len(primary) == 1574 * 3 and primary.inat_taxon_id.nunique() == 1574, "Primary slope panel changed")
    require(primary.groupby("inat_taxon_id").n_primary.first().sum() == 86114, "Primary photo panel changed")
    native = convert_primary_slopes(primary)
    old_wide = pd.read_csv(paths["gift"] / "species_slopes_with_traits.csv")
    require(len(old_wide) == 1574 and old_wide.inat_taxon_id.is_unique, "Frozen GIFT panel changed")
    original_pivot = primary.pivot(index="inat_taxon_id", columns="variable", values="slope")
    ordered_old = old_wide.set_index("inat_taxon_id").loc[original_pivot.index]
    np.testing.assert_allclose(ordered_old[list(PREDICTORS)], original_pivot[list(PREDICTORS)], rtol=1e-12, atol=1e-12)
    primary_names = primary.groupby("inat_taxon_id").species.first()
    require(ordered_old.species.equals(primary_names), "Exact GIFT species identities changed")
    traits = pd.read_csv(paths["gift"] / "species_trait_join.csv").set_index("inat_taxon_id")
    require(traits.index.is_unique and len(traits) == 4475, "GIFT roster changed")
    for trait in CONTRASTS:
        require(ordered_old[trait].fillna("").equals(traits.loc[ordered_old.index, trait].fillna("")), "Frozen GIFT category differs: " + trait)
    new_wide = old_wide.copy()
    converted_pivot = native.pivot(index="inat_taxon_id", columns="variable", values="slope")
    for variable in PREDICTORS:
        new_wide[variable] = new_wide.inat_taxon_id.map(converted_pivot[variable])
    # Reload exactly the authenticated panel to check the fitted predictor for
    # every observed photograph. No optimization or raw photograph model refit.
    panel, sampled_rasters, missing = load_panel(paths["photos"], paths["worldclim"])
    prediction_audit = observed_prediction_and_covariance_audit(primary, native, panel)
    print("EXACT_PREDICTION_AUDIT 1574 species / 86114 photos passed", flush=True)
    summary, scale_checks = summarize_native(native)
    genus_summary, genus_draws = mean_bootstrap(new_wide)
    old_genus_summary, old_genus_draws = mean_bootstrap(old_wide)
    frozen_contrasts = pd.read_csv(paths["gift"] / "slope_modifier_contrasts.csv")
    contrasts, contrast_comparison = trait_comparisons(old_wide, new_wide, frozen_contrasts)
    ranges, contributions, tail_species = range_diagnostics(native)
    old_meta = pd.read_csv(paths["environment"] / "meta_summary.csv")
    old_meta = old_meta.loc[old_meta.model.eq("primary") & old_meta.group.eq("pooled")]
    primary_summary = summary.loc[summary.group.eq("pooled")].merge(genus_summary.drop(columns="n_species"), on="variable", validate="one_to_one")
    selected_meta = ["reml_mean", "reml_mean_ci_low", "reml_mean_ci_high", "tau2", "tau2_profile_ci_low", "tau2_profile_ci_high", "equal_species_unshrunk_mean"]
    scale_comparison = primary_summary.merge(old_meta[["variable", *selected_meta]].rename(columns=lambda c: c if c == "variable" else "sd_scale_" + c), on="variable", validate="one_to_one")
    scale_comparison = scale_comparison.merge(old_genus_summary[["variable", "genus_cluster_percentile_low", "genus_cluster_percentile_high"]].rename(columns=lambda c: c if c == "variable" else "sd_scale_" + c), on="variable", validate="one_to_one")
    frames = {
        "input_verification.csv": pd.DataFrame(authentication),
        "species_slopes.csv": native, "species_slopes_with_traits.csv": new_wide,
        "prediction_covariance_audit.csv": prediction_audit,
        "meta_summary.csv": summary, "reml_numerical_scale_checks.csv": scale_checks,
        "genus_cluster_bootstrap.csv": genus_summary,
        "primary_summary.csv": primary_summary, "overall_scale_comparison.csv": scale_comparison,
        "slope_modifier_contrasts.csv": contrasts, "gift_scale_comparison.csv": contrast_comparison,
        "environment_range_diagnostics.csv": ranges, "range_tail_contributions.csv": contributions,
        "range_tail_species.csv": tail_species,
    }
    for name, frame in frames.items():
        frame.to_csv(outdir / name, index=False)
    np.savez_compressed(outdir / "equal_species_bootstrap_draws.npz", variables=np.array(PREDICTORS),
                        native=genus_draws, standardized=old_genus_draws)
    profile = {
        "schema": "fcp_environment_common_physical_units_v1", "protocol_section": "A",
        "analysis_type": "post_hoc_internal_reexpression_of_existing_joint_Firth_fits",
        "n_species": 1574, "n_photos": 86114, "n_joint_slope_rows": len(native),
        "n_species_trimmed": 0, "physical_increments": PHYSICAL_INCREMENTS, "unit_labels": UNIT_LABELS,
        "response": "white=1 versus other three classified colour groups=0",
        "transform": "b_j*d_j/s_j; V_jk*d_j*d_k/(s_j*s_k); intercept and fitted probabilities unchanged",
        "primary_estimand": "equal-species arithmetic mean of unshrunk log-odds slopes per common physical increment",
        "secondary_reml_estimand": "marginal normal-normal REML hypermean with inverse total-variance weights in common physical units",
        "reml_method": "unchanged frozen reml_meta, with a common median-SE numerical normalization undone in all dimensional outputs; modified HKSJ and restricted-likelihood tau2 profiles",
        "reml_numerical_checks": {"n_checks": len(scale_checks), "all_passed": bool(scale_checks.passed.all())},
        "equal_mean_independence_interval": "finite-panel approximate sampling variance sum(V_s)/K^2, with a normal interval; assumes independent photos/species",
        "overall_bootstrap": {"repetitions": REPETITIONS, "seed": MEAN_BOOTSTRAP_SEED,
                              "shared_axes_and_unit_scales": True, "resampling": "genus clusters, retaining all species and equal species weight within each draw"},
        "gift_bootstrap": {"repetitions_per_trait": REPETITIONS, "seeds": {trait: TRAIT_BASE_SEED + i for i, trait in enumerate(CONTRASTS)},
                           "shared_trait_groups_axes_and_unit_scales": True, "contrasts": CONTRASTS,
                           "same_coverage_and_categories_as_frozen_analysis": True, "standardized_analysis_exactly_reproduced": True,
                           "nominal_interval": "2.5th to 97.5th bootstrap percentiles",
                           "family9_interval": "0.05/(2*9) to 1-0.05/(2*9) bootstrap quantiles; approximate Bonferroni family9 intervals",
                           "bootstrap_expected_draws_per_family9_tail": REPETITIONS * .05 / (2 * FAMILY_SIZE)},
        "prediction_and_covariance_checks": {"all_species_passed": bool(prediction_audit.passed.all()),
             "n_observed_photos": int(prediction_audit.n_photos.sum()),
             **{c: float(prediction_audit[c].max()) for c in prediction_audit.columns if c.startswith("max_")},
             "score_variance_scope": "slope-only centered linear-predictor variance; original files do not contain intercept-slope covariances"},
        "range_diagnostics": "all slopes retained; numerical near-zero gate rechecked; SD<0.1 specified increment and lowest 1%/5% SD groups are diagnostic only",
        "range_summary": ranges.to_dict(orient="records"),
        "interpretation": [
            "Each individual slope sign, Wald statistic, covariance correlation and observed fitted probability is invariant to re-expression.",
            "Mean(b_s/SD_s) differs from mean(b_s); standardized and common-unit analyses average distinct quantities.",
            "REML inverse-total-variance weights also change when species-specific scales change; differing summaries are not refutations of one another.",
            "A physical increment larger than a species observed climate span is a slope unit, not validation of predictions over that increment.",
            "Very small sampled environmental ranges can amplify both unshrunk physical-unit estimates and their uncertainty.",
        ],
        "limitations": [
            "The fixed 1574-species panel was selected using full-data response counts and seven-variable climate completeness; this analysis is conditional on that selected panel.",
            "Genus bootstrap intervals describe empirical stability of fitted slopes, without photo refitting, within-fit variance simulation, or full phylogenetic covariance.",
            "The nine GIFT intervals have only about fourteen expected bootstrap draws per extreme family9 tail and therefore appreciable Monte Carlo approximation.",
            "REML is an approximate two-stage marginal normal model using approximate Firth Fisher sampling covariances; cross-axis between-species covariance is not modeled.",
            "Tau=0 is a fitted heterogeneity boundary with a nonzero profile upper bound; it does not establish species equivalence or absence of heterogeneity.",
            "GIFT missingness remains unknown; self-fertilization presence is not realized selfing rate. No causal adaptation, molecular mechanism, pigment cost or fitness effect is identified.",
        ],
    }
    (outdir / "profile.json").write_text(json.dumps(profile, indent=2, allow_nan=False) + "\n")
    artifact_names = [*frames, "equal_species_bootstrap_draws.npz", "protocol_at_execution.md", "profile.json"]
    receipt = {
        "started_utc": started, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv, "python": sys.version, "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "git_head_at_execution": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "script_sha256": sha(Path(__file__)),
        "test_sha256": sha(root / "tests/test_fcp_environment_native_units_20261019.py"),
        "protocol_source_sha256_at_execution": hashlib.sha256(protocol_text.encode()).hexdigest(),
        "source_receipt_sha256": paths["source_receipt_sha256"],
        "input_hash_references_verified": len(authentication),
        "output_sha256": {name: sha(outdir / name) for name in artifact_names},
    }
    (outdir / "execution_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(primary_summary[["variable", "equal_species_unshrunk_mean", "genus_cluster_percentile_low", "genus_cluster_percentile_high", "reml_mean", "reml_mean_ci_low", "reml_mean_ci_high", "tau2", "tau2_profile_ci_high"]].to_string(index=False), flush=True)
    print(contrasts[["trait", "variable", "difference_first_minus_reference", "percentile_ci_low", "percentile_ci_high", "family9_percentile_ci_low", "family9_percentile_ci_high"]].to_string(index=False), flush=True)
    print("FCP_NATIVE_UNITS_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
