#!/usr/bin/env python3
"""Exact post hoc white/nonwhite decomposition on the original 50-km pair graph.

Coarse photographic states do not identify pigments, genes or wild morph
frequencies. The original four-state result is an authentication target, not a
new confirmation. Only the pooled two component tests receive new p-values.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_447500_robust_matched_null_20261018 import (
    COHORTS, COUNTS, CLASSIFIED, LOCAL_COUNTS, SEEDS, INPUTS, canonical_photos, sha,
)
from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    MORPHS, MORPH_CODE, stable_seed, pairwise_geo_km, as_bool, PERMUTATIONS,
    PRIMARY_RADIUS, MIN_LOCAL_PAIRS,
)

NAMES = ("total", "white_nonwhite", "nonwhite_nonwhite")
BASE_NULL_SHA = "e0623b881379fec3533558b302deeccc0bc4c69b132f068ca7ab9f2e9951e80b"
BASE_RESULT_SHA = "dad77a486f9477f8cad0e2aaf53859d3ec9fb2880b6f968b9854ef61cc1b765e"
TOL = 1e-12


def demand(ok, message):
    if not ok:
        raise ValueError("FCP_WHITE_DECOMPOSITION_INVALID: " + message)


def component_worlds(labels, permutations, u, v, local):
    """Return all-pair global, observed local/delta, and null local/delta.

    Components always share the full original pair denominator, including
    equal-colour pairs. No white vertex is removed or eligibility reselected.
    """
    labels = np.asarray(labels, np.int8)
    permutations = np.asarray(permutations, np.int8)
    demand(permutations.ndim == 2 and permutations.shape[1] == len(labels),
           "label/null shape differs")
    demand(len(u) == len(v) == len(local) and np.any(local), "invalid local pair graph")
    a, b = labels[u], labels[v]
    total = a != b
    white = (a == 0) != (b == 0)
    chromatic = total & (a != 0) & (b != 0)
    demand(np.array_equal(total.astype(int), white.astype(int) + chromatic.astype(int)),
           "observed pair partition failed")
    global_values = np.array([x.mean() for x in (total, white, chromatic)], float)
    local_values = np.array([x[local].mean() for x in (total, white, chromatic)], float)
    aa, bb = permutations[:, u[local]], permutations[:, v[local]]
    ztotal = aa != bb
    zwhite = (aa == 0) != (bb == 0)
    zchromatic = ztotal & (aa != 0) & (bb != 0)
    demand(np.array_equal(ztotal.astype(np.int8), zwhite.astype(np.int8) + zchromatic.astype(np.int8)),
           "null pair partition failed")
    null_local = np.column_stack([x.mean(axis=1) for x in (ztotal, zwhite, zchromatic)])
    delta = global_values - local_values
    null_delta = global_values[None, :] - null_local
    return global_values, local_values, delta, null_local, null_delta


def check_reconstruction(delta, null_delta, old_observed, old_null, context):
    demand(np.allclose(delta[0], delta[1] + delta[2], atol=TOL, rtol=0),
           context + ": observed additive identity failed")
    demand(np.allclose(null_delta[:, 0], null_delta[:, 1] + null_delta[:, 2], atol=TOL, rtol=0),
           context + ": per-world additive identity failed")
    demand(np.allclose(delta[0], old_observed, atol=TOL, rtol=0),
           context + ": original observed total differs")
    demand(np.allclose(null_delta[:, 0], old_null, atol=TOL, rtol=0),
           context + ": original 199 vertex worlds differ")


def composition_table(photo):
    identity = photo.groupby("inat_taxon_id", sort=True).agg(
        species=("species", "first"), source_cohort=("fcp_source_cohort", "first"),
        measured_photos=("photo_id", "size"))
    valid = photo.loc[as_bool(photo.global_classifiable) & photo.morph.isin(MORPHS)]
    count = pd.crosstab(valid.inat_taxon_id, valid.morph).reindex(
        index=identity.index, columns=MORPHS, fill_value=0).fillna(0).astype(int)
    out = identity.copy()
    out["n_classifiable"] = count.sum(axis=1)
    out["n_states"] = (count > 0).sum(axis=1)
    out["n_white"] = count.white
    out["n_nonwhite"] = out.n_classifiable - out.n_white
    out["white_nonwhite_pairs"] = out.n_white * out.n_nonwhite
    # Integer counts of unordered unequal pairs avoid finite-sample ambiguity.
    out["all_unequal_pairs"] = ((out.n_classifiable ** 2 - (count ** 2).sum(axis=1)) // 2).astype(int)
    out["nonwhite_unequal_pairs"] = out.all_unequal_pairs - out.white_nonwhite_pairs
    out["fraction_unequal_pairs_involving_white"] = np.where(
        out.all_unequal_pairs > 0, out.white_nonwhite_pairs / out.all_unequal_pairs, np.nan)
    return out.reset_index()


def describe_composition(table):
    out = {}
    for name, subset in (("all_measured_taxa", table),
                         ("ge40_classifiable", table.loc[table.n_classifiable >= 40])):
        poly = subset.n_states >= 2
        white_poly = poly & (subset.n_white > 0) & (subset.n_nonwhite > 0)
        supported = (subset.n_white >= 5) & (subset.n_nonwhite >= 5)
        any_unequal = subset.all_unequal_pairs > 0
        fractions = (subset.loc[any_unequal, "white_nonwhite_pairs"] /
                     subset.loc[any_unequal, "all_unequal_pairs"])
        n, npoly = len(subset), int(poly.sum())
        den = int(subset.all_unequal_pairs.sum())
        out[name] = {
            "n_species": n, "n_classifiable_photos": int(subset.n_classifiable.sum()),
            "observed_multistate_species": npoly,
            "white_in_observed_polymorphism_species": int(white_poly.sum()),
            "white_in_observed_polymorphism_fraction": float(white_poly.sum() / n) if n else None,
            "white_fraction_among_observed_multistate_species": float(white_poly.sum() / npoly) if npoly else None,
            "white_ge5_nonwhite_ge5_species": int(supported.sum()),
            "white_ge5_nonwhite_ge5_fraction_of_all_species_in_panel": float(supported.sum() / n) if n else None,
            "all_unequal_pairs": den,
            "white_nonwhite_pairs": int(subset.white_nonwhite_pairs.sum()),
            "pooled_photo_pair_weighted_white_fraction_of_unequal_pairs":
                float(subset.white_nonwhite_pairs.sum() / den) if den else None,
            "species_equal_mean_white_fraction_of_unequal_pairs_among_multistate_species":
                float(fractions.mean()) if len(fractions) else None,
            "species_median_white_fraction_of_unequal_pairs_among_multistate_species":
                float(fractions.median()) if len(fractions) else None,
            "fraction_average_denominator_multistate_species": int(any_unequal.sum()),
        }
    return out


def component_summary(global_values, local_values, null_local, *, inferential=True):
    global_values, local_values, null_local = map(np.asarray, (global_values, local_values, null_local))
    demand(global_values.shape == local_values.shape and global_values.shape[1] == 3,
           "summary component shape differs")
    demand(null_local.shape == (len(global_values), PERMUTATIONS, 3), "summary original-null shape differs")
    null_delta = global_values[:, None, :] - null_local
    delta = global_values - local_values
    mean_null_delta = null_delta.mean(axis=0)
    mean_null_local = null_local.mean(axis=0)
    summary = {}
    for k, name in enumerate(NAMES):
        p = float((1 + np.count_nonzero(mean_null_delta[:, k] >= delta[:, k].mean())) / 200)
        summary[name] = {
            "n_species": len(global_values),
            "observed_mean_global_mismatch": float(global_values[:, k].mean()),
            "observed_mean_local_mismatch": float(local_values[:, k].mean()),
            "observed_mean_depletion": float(delta[:, k].mean()),
            "observed_positive_depletion_species": int((delta[:, k] > 0).sum()),
            "null_mean_global_mismatch": float(global_values[:, k].mean()),
            "null_mean_local_mismatch": float(mean_null_local[:, k].mean()),
            "null_mean_local_q025": float(np.quantile(mean_null_local[:, k], .025)),
            "null_mean_local_q975": float(np.quantile(mean_null_local[:, k], .975)),
            "null_mean_depletion": float(mean_null_delta[:, k].mean()),
            "null_mean_depletion_q025": float(np.quantile(mean_null_delta[:, k], .025)),
            "null_mean_depletion_q975": float(np.quantile(mean_null_delta[:, k], .975)),
            "p_upper_posthoc": p if inferential and k > 0 else None,
            "holm_two_component_p_posthoc": None,
        }
    if inferential:
        order = sorted(NAMES[1:], key=lambda name: summary[name]["p_upper_posthoc"])
        floor = 0.
        for rank, name in enumerate(order):
            floor = max(floor, min(1., (2 - rank) * summary[name]["p_upper_posthoc"]))
            summary[name]["holm_two_component_p_posthoc"] = floor
    return summary


def authenticate_inputs(root, original_dir):
    paths = {}
    for name, (relative, digest) in INPUTS.items():
        p = root / relative
        demand(p.is_file() and sha(p) == digest, "source SHA mismatch: " + relative)
        paths[name] = p
    baseline = original_dir / "matched_vertex_nulls.npz"
    baseline_result = original_dir / "result.json"
    demand(sha(baseline) == BASE_NULL_SHA and sha(baseline_result) == BASE_RESULT_SHA,
           "original authenticated 1761 x 199 null/result SHA mismatch")
    return paths, baseline, baseline_result


def calculate(photo, spatial, original_npz):
    demand(PERMUTATIONS == 199 and PRIMARY_RADIUS == 50 and MIN_LOCAL_PAIRS == 30,
           "original method constants changed")
    photo = canonical_photos(photo)
    demand(len(photo) == 447500 and photo.inat_taxon_id.nunique() == 4475,
           "447500/4475 denominator differs")
    demand((photo.groupby("inat_taxon_id").size() == 100).all(), "not 100 original photos per species")
    demand(photo.groupby("fcp_source_cohort").inat_taxon_id.nunique().to_dict() == COUNTS,
           "source cohort counts differ")
    valid = photo.loc[as_bool(photo.global_classifiable) & photo.morph.isin(MORPHS)].copy()
    demand(valid.groupby("fcp_source_cohort").size().to_dict() == CLASSIFIED,
           "227107 classified rows differ")
    demand(valid[["latitude", "longitude"]].notna().all().all(), "classified coordinates missing")
    profile = composition_table(photo)
    demand(len(profile.loc[profile.n_classifiable >= 40]) == 3297, "ge40 species count differs")
    ids = original_npz["taxon_ids"]
    demand(len(ids) == 1761 and len(np.unique(ids)) == 1761, "original null taxon identity differs")
    eligible = spatial.loc[spatial.depletion_50km.notna()].set_index("inat_taxon_id", verify_integrity=True)
    demand(set(ids) == set(eligible.index), "saved spatial/null taxa differ")
    demand(eligible.groupby("source_cohort").size().to_dict() == LOCAL_COUNTS,
           "local cohort denominator differs")
    groups = valid.loc[valid.inat_taxon_id.isin(ids)].groupby("inat_taxon_id", sort=True)
    global_rows, local_rows, null_rows, records = [], [], [], []
    max_original_error = max_additivity_error = 0.
    for index, taxon in enumerate(ids):
        g = groups.get_group(int(taxon))
        source = str(g.fcp_source_cohort.iloc[0])
        old = eligible.loc[taxon]
        demand(source == str(old.source_cohort) == COHORTS[int(original_npz["cohort_codes"][index])],
               f"taxon {taxon}: original source cohort differs")
        demand(len(g) == int(old.n_classifiable) and len(g) >= 40 and g.photo_id.is_monotonic_increasing,
               f"taxon {taxon}: original classifiability/order differs")
        labels = np.array([MORPH_CODE[str(x)] for x in g.morph], np.int8)
        rng = np.random.default_rng(stable_seed(SEEDS[source], int(taxon)))
        permutations = np.stack([rng.permutation(labels) for _ in range(199)])
        u, v = np.triu_indices(len(labels), 1)
        distance = pairwise_geo_km(g.latitude.to_numpy(float), g.longitude.to_numpy(float))[u, v]
        local = distance <= 50
        demand(int(local.sum()) == int(old.local_pairs_50km) and local.sum() >= 30,
               f"taxon {taxon}: original local graph differs")
        glob, loc, delta, null_local, null_delta = component_worlds(labels, permutations, u, v, local)
        check_reconstruction(delta, null_delta, original_npz["observed"][index],
                             original_npz["nulls"][index], f"taxon {taxon}")
        demand(abs(glob[0] - old.D_pair) <= TOL and abs(loc[0] - old.D_local_50km) <= TOL,
               f"taxon {taxon}: saved original global/local mismatch differs")
        max_original_error = max(max_original_error, float(np.max(np.abs(null_delta[:, 0] - original_npz["nulls"][index]))))
        max_additivity_error = max(max_additivity_error, float(np.max(np.abs(null_delta[:, 0] - null_delta[:, 1] - null_delta[:, 2]))))
        global_rows.append(glob); local_rows.append(loc); null_rows.append(null_local)
        row = {"inat_taxon_id": int(taxon), "species": str(g.species.iloc[0]),
               "source_cohort": source, "n_classifiable": len(g), "local_pairs_50km": int(local.sum())}
        for k, name in enumerate(NAMES):
            row.update({f"global_{name}": float(glob[k]), f"local_{name}": float(loc[k]),
                        f"depletion_{name}": float(delta[k]),
                        f"null_mean_depletion_{name}": float(null_delta[:, k].mean())})
        records.append(row)
        if (index + 1) % 100 == 0 or index + 1 == len(ids):
            print(f"FCP_WHITE_DECOMPOSITION_PROGRESS {index+1}/1761", flush=True)
    glob, loc, null_loc = map(np.stack, (global_rows, local_rows, null_rows))
    species = pd.DataFrame(records)
    pooled = component_summary(glob, loc, null_loc)
    by_cohort = {}
    for cohort in COHORTS:
        ii = species.source_cohort.eq(cohort).to_numpy()
        by_cohort[cohort] = component_summary(glob[ii], loc[ii], null_loc[ii], inferential=False)
    total, white, chrom = [pooled[name]["observed_mean_depletion"] for name in NAMES]
    result = {
        "schema": "fcp_real_447500_exact_white_axis_spatial_decomposition_v1",
        "date": "2026-10-18", "status": "complete_source_and_original_null_authenticated",
        "method": {"radius_km": 50, "radius_post_hoc": True, "original_vertex_worlds": 199,
                   "same_1761_species_same_all_pair_denominators": True,
                   "white_vertices_removed": False, "equal_species_weighting": True,
                   "seed_mapping": SEEDS,
                   "new_test_family": "two pooled one-sided positive-depletion component tests; Holm2",
                   "cohort_summaries": "descriptive matched-null summaries; no new cohort p-values",
                   "total_result": "original result authenticated, not a new confirmation",
                   "unequal_pair_partition": "I(a!=b) = I(exactly_one_white) + I(a!=b and neither_white)"},
        "composition_definitions": {
            "observed_multistate": "at least two coarse photo states among retained classifiable photos",
            "white_in_observed_polymorphism": "at least one white and at least one nonwhite photo",
            "support_sensitivity": "at least five white and at least five nonwhite photos",
            "main_panel": "3297 species with >=40 four-state-classifiable photos; all4475 also shown",
            "pair_fraction": "unordered white-nonwhite pairs / all unordered unequal-state pairs",
            "zero_unequal_pairs": "monomorphic species excluded from pair-fraction averages only",
            "no_population_prevalence": True},
        "composition": describe_composition(profile),
        "pooled_spatial": pooled, "by_cohort_descriptive": by_cohort,
        "pooled_effect_decomposition": {
            "total": total, "white_nonwhite": white, "nonwhite_nonwhite": chrom,
            "white_nonwhite_over_positive_pooled_total": white / total if total > 0 else None,
            "nonwhite_nonwhite_over_positive_pooled_total": chrom / total if total > 0 else None,
            "ratio_definition": "ratio of species-equal pooled mean effects, not mean of species ratios",
            "ratio_caveat": "components may be negative or cancel; this is descriptive and not a bounded causal variance fraction",
            "species_ratios_not_computed": True},
        "exact_reconstruction": {"species": 1761, "worlds_per_species": 199,
                                 "original_total_null_values_checked": 1761 * 199,
                                 "max_absolute_original_total_null_error": max_original_error,
                                 "max_absolute_additivity_error": max_additivity_error,
                                 "absolute_tolerance": TOL},
        "nonclaims": ["Post hoc decomposition of existing photographic structure, not a new prospective test",
                      "Coarse white/nonwhite state is not pigment chemistry, anthocyanin loss or a genotype",
                      "No fitness, pigment production cost, selection or local adaptation measured",
                      "Photo state proportions are not natural plant/population morph frequencies",
                      "This pairwise decomposition is distinct from sign-invariant continuous H2 W",
                      "Same source and exposure/observer biases remain"]}
    return result, species, profile, glob, loc, null_loc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extracted-root", type=Path, required=True)
    parser.add_argument("--original-null-dir", type=Path,
                        default=Path("results/fcp_real_447500_robust_matched_null_20261018"))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    paths, baseline, baseline_result = authenticate_inputs(args.extracted_root, args.original_null_dir)
    print("FCP_WHITE_DECOMPOSITION_SOURCE_HASHES_AUTHENTICATED", flush=True)
    columns = ["inat_taxon_id", "species", "photo_id", "fcp_source_cohort", "morph",
               "global_classifiable", "latitude", "longitude"]
    photo = pd.read_csv(paths["photos"], usecols=columns, low_memory=False)
    original = np.load(baseline, allow_pickle=False)
    result, species, profile, glob, loc, null_loc = calculate(photo, pd.read_csv(paths["species"]), original)
    result["source"] = {"original_real_measurement_release": "fcp-all4730-original-method-actual-colour-results-20261017",
                        "input_files": {key: {"path": str(paths[key].resolve()), "sha256": val[1]}
                                        for key, val in INPUTS.items()},
                        "original_total_null_npz_sha256": BASE_NULL_SHA,
                        "original_total_robust_result_sha256": BASE_RESULT_SHA}
    args.outdir.mkdir(parents=True, exist_ok=True)
    species.to_csv(args.outdir / "species_spatial_decomposition.csv", index=False)
    profile.to_csv(args.outdir / "species_photo_composition.csv", index=False)
    summary = [{"group": "pooled", "component": name, **vals} for name, vals in result["pooled_spatial"].items()]
    summary += [{"group": cohort, "component": name, **vals}
                for cohort, comps in result["by_cohort_descriptive"].items() for name, vals in comps.items()]
    pd.DataFrame(summary).to_csv(args.outdir / "component_summary.csv", index=False)
    np.savez_compressed(args.outdir / "component_vertex_nulls.npz", taxon_ids=species.inat_taxon_id.to_numpy(np.int64),
                        global_components=glob, local_components=loc, null_local_components=null_loc)
    output_names = ("species_spatial_decomposition.csv", "species_photo_composition.csv", "component_summary.csv",
                    "component_vertex_nulls.npz")
    result["output_sha256"] = {name: sha(args.outdir / name) for name in output_names}
    result_path = args.outdir / "result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    repo = Path(__file__).resolve().parents[2]
    source_names = [Path(__file__).resolve(),
                    repo / "scripts/analysis/run_fcp_447500_robust_matched_null_20261018.py",
                    repo / "scripts/analysis/run_polymorphism_distributed_polymorphism_posthoc_20261007.py"]
    receipt = {"schema": "fcp_white_axis_decomposition_execution_receipt_v1", "actual_inference_replay": True,
               "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
               "python": platform.python_version(), "packages": {"numpy": np.__version__, "pandas": pd.__version__},
               "command": sys.argv, "source_code_sha256": {str(p.relative_to(repo)): sha(p) for p in source_names},
               "unit_test_source_sha256": sha(repo / "tests/test_run_fcp_447500_white_axis_decomposition_20261018.py"),
               "result_json_sha256": sha(result_path),
               "original_npz_sha256": sha(baseline), "original_result_sha256": sha(baseline_result)}
    (args.outdir / "execution_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print("FCP_WHITE_AXIS_DECOMPOSITION_COMPLETE " + json.dumps(result["pooled_effect_decomposition"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
