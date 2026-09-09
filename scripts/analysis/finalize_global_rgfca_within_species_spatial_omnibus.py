#!/usr/bin/env python3
"""Finalize exact RGFCA-only within-species spatial omnibus species shards."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_global_rgfca_within_species_spatial_omnibus as run
from fcp_pipeline.global_g3 import species_distance_colour_rho
from fcp_pipeline.global_rgfca_engine import COLOUR_COLUMNS


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def direct_spot_checks(measured_path: Path, table: pd.DataFrame) -> list[dict[str, object]]:
    pool = run.load_rgfca_pool(measured_path)
    species_labels = sorted(pool.species.unique())
    selected_species = [species_labels[0], species_labels[len(species_labels)//2], species_labels[-1]]
    checks: list[dict[str, object]] = []
    for species in selected_species:
        q = pool.loc[pool.species == species]
        lat = q.latitude.to_numpy(float)
        lon = q.longitude.to_numpy(float)
        colours = q[list(COLOUR_COLUMNS)].to_numpy(float)
        observed_direct = species_distance_colour_rho(lat, lon, colours)
        observed_saved = float(table[(table.species == species) & (table.permutation_index == -1)].rho.iloc[0])
        error = abs(observed_direct - observed_saved)
        checks.append({
            "species": species,
            "permutation_index": -1,
            "saved_rho": observed_saved,
            "direct_scipy_rho": observed_direct,
            "absolute_error": error,
        })
        if error > 2e-12:
            raise RuntimeError(f"observed direct Spearman mismatch for {species}: {error}")
        for permutation_index in (0, 499, 998):
            perm = np.random.default_rng(run.seed_for(species, permutation_index)).permutation(len(q))
            direct = species_distance_colour_rho(lat, lon, colours[perm])
            saved = float(table[(table.species == species) & (table.permutation_index == permutation_index)].rho.iloc[0])
            error = abs(direct - saved)
            checks.append({
                "species": species,
                "permutation_index": permutation_index,
                "saved_rho": saved,
                "direct_scipy_rho": direct,
                "absolute_error": error,
            })
            if error > 2e-12:
                raise RuntimeError(f"permuted direct Spearman mismatch for {species}/{permutation_index}: {error}")
    return checks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=Path, required=True)
    ap.add_argument("--measured", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh finalizer output")

    result_files = sorted(args.shards.rglob("species_rho_permutations.csv"))
    audit_files = sorted(args.shards.rglob("species_audit.csv"))
    meta_files = sorted(args.shards.rglob("meta.json"))
    if not (len(result_files) == len(audit_files) == len(meta_files) == 20):
        raise RuntimeError(f"expected 20 shard triplets, got {len(result_files)}/{len(audit_files)}/{len(meta_files)}")
    metas = [json.loads(p.read_text()) for p in meta_files]
    if {int(m["shard_index"]) for m in metas} != set(range(20)):
        raise RuntimeError("species shard index census mismatch")
    if any(int(m["n_shards"]) != 20 for m in metas):
        raise RuntimeError("species shard count drift")
    if any(m["specification_commit"] != run.SPECIFICATION_COMMIT for m in metas):
        raise RuntimeError("specification drift")
    if any(m["measured_table_sha256"] != run.EXPECTED_MEASURED_SHA256 for m in metas):
        raise RuntimeError("measured-table hash drift across shards")
    if len({m["runner_sha256"] for m in metas}) != 1:
        raise RuntimeError("runner differs between species shards")
    if any(m.get("six_species_used") is not False or m.get("thirty_four_species_used") is not False for m in metas):
        raise RuntimeError("non-RGFCA evidence firewall failure")
    if any(m.get("environmental_variables_used") is not False for m in metas):
        raise RuntimeError("environmental-variable firewall failure")
    if any(m.get("parent_G1_reclassified") is not False for m in metas):
        raise RuntimeError("G1 reclassification guard failure")

    table = pd.concat([pd.read_csv(p) for p in result_files], ignore_index=True)
    audits = pd.concat([pd.read_csv(p) for p in audit_files], ignore_index=True)
    if not np.isfinite(table.rho.to_numpy(float)).all():
        raise RuntimeError("nonfinite species rho")
    if table.duplicated(["species", "permutation_index"]).any():
        raise RuntimeError("duplicate species/permutation row")
    species = sorted(table.species.unique())
    if len(species) != run.EXPECTED_SPECIES or len(audits) != run.EXPECTED_SPECIES or audits.species.nunique() != run.EXPECTED_SPECIES:
        raise RuntimeError("369-species census failure")
    if sum(int(m["species_count"]) for m in metas) != run.EXPECTED_SPECIES:
        raise RuntimeError("shard species total mismatch")
    expected_index = set([-1] + list(range(run.N_PERM)))
    for species_label, q in table.groupby("species", sort=False):
        if len(q) != run.N_PERM + 1 or set(q.permutation_index.astype(int)) != expected_index:
            raise RuntimeError(f"permutation census failure for {species_label}")

    observed = table[table.permutation_index == -1].set_index("species").loc[species].rho.astype(float)
    null = table[table.permutation_index >= 0].copy()
    null_global = null.groupby("permutation_index", sort=True).rho.mean()
    if len(null_global) != run.N_PERM or list(null_global.index.astype(int)) != list(range(run.N_PERM)):
        raise RuntimeError("global null permutation index census failure")
    observed_mean = float(observed.mean())
    observed_median = float(observed.median())
    observed_positive_fraction = float((observed > 0).mean())
    null_mean_values = null_global.to_numpy(float)
    primary_p = float((1 + np.count_nonzero(null_mean_values >= observed_mean)) / (run.N_PERM + 1))

    null_secondary = []
    pivot = null.pivot(index="permutation_index", columns="species", values="rho").loc[:, species]
    for permutation_index, row in pivot.iterrows():
        values = row.to_numpy(float)
        null_secondary.append({
            "permutation_index": int(permutation_index),
            "mean_rho": float(np.mean(values)),
            "median_rho": float(np.median(values)),
            "positive_fraction": float(np.mean(values > 0)),
        })
    null_secondary_frame = pd.DataFrame(null_secondary)

    species_rows = []
    for species_label in species:
        obs = float(observed.loc[species_label])
        values = null.loc[null.species == species_label].sort_values("permutation_index").rho.to_numpy(float)
        if len(values) != run.N_PERM:
            raise RuntimeError("species null vector length mismatch")
        p = float((1 + np.count_nonzero(values >= obs)) / (run.N_PERM + 1))
        audit = audits.loc[audits.species == species_label].iloc[0]
        species_rows.append({
            "species": species_label,
            "photos": int(audit.photos),
            "pairs": int(audit.pairs),
            "observed_rho": obs,
            "null_mean_rho": float(np.mean(values)),
            "null_q025_rho": float(np.quantile(values, 0.025)),
            "null_q975_rho": float(np.quantile(values, 0.975)),
            "p_upper_exploratory": p,
        })
    species_frame = pd.DataFrame(species_rows)
    species_frame["p_bh_exploratory"] = bh_adjust(species_frame.p_upper_exploratory.to_numpy(float))
    species_frame["bh_detectable_0_05"] = species_frame.p_bh_exploratory < 0.05
    bh_count = int(species_frame.bh_detectable_0_05.sum())

    spot_checks = direct_spot_checks(args.measured, table)
    max_spot_error = float(max(x["absolute_error"] for x in spot_checks))

    args.output.mkdir(parents=True)
    species_frame.to_csv(args.output / "species_results.csv", index=False, lineterminator="\n")
    null_secondary_frame.to_csv(args.output / "global_null.csv", index=False, lineterminator="\n")
    audits.sort_values("species", kind="stable").to_csv(args.output / "species_audit.csv", index=False, lineterminator="\n")
    pd.DataFrame(spot_checks).to_csv(args.output / "independent_spot_checks.csv", index=False, lineterminator="\n")

    summary = {
        "protocol": run.PROTOCOL,
        "status": "complete_pending_external_ci_confirmation",
        "specification_commit": run.SPECIFICATION_COMMIT,
        "inferential_role": "postoutcome_exploratory_rgfca_only_exact_randomization",
        "eligible_species": run.EXPECTED_SPECIES,
        "classifiable_rows": run.EXPECTED_ROWS,
        "permutations": run.N_PERM,
        "primary": {
            "statistic": "equal_species_mean_distance_colour_spearman_rho",
            "observed_mean_rho": observed_mean,
            "null_mean": float(np.mean(null_mean_values)),
            "null_sd": float(np.std(null_mean_values, ddof=1)),
            "null_q025": float(np.quantile(null_mean_values, 0.025)),
            "null_q975": float(np.quantile(null_mean_values, 0.975)),
            "p_upper": primary_p,
            "supported_at_0_05": bool(primary_p < 0.05),
        },
        "secondary_descriptive": {
            "observed_median_species_rho": observed_median,
            "observed_positive_species_fraction": observed_positive_fraction,
            "null_positive_fraction_mean": float(null_secondary_frame.positive_fraction.mean()),
            "null_positive_fraction_q025": float(null_secondary_frame.positive_fraction.quantile(0.025)),
            "null_positive_fraction_q975": float(null_secondary_frame.positive_fraction.quantile(0.975)),
            "null_median_rho_mean": float(null_secondary_frame.median_rho.mean()),
            "bh_detectable_species_count_0_05": bh_count,
            "bh_detectable_species_fraction_0_05": bh_count / run.EXPECTED_SPECIES,
            "bh_label": "exploratory detectable-species count only; not biological prevalence",
        },
        "verification": {
            "twenty_nonoverlapping_species_shards": True,
            "all_369_species_present_once_per_observed_and_null_index": True,
            "all_999_global_null_statistics_reconstructed": True,
            "primary_p_recomputed_from_global_null": True,
            "bh_adjustment_recomputed": True,
            "direct_scipy_spot_checks": len(spot_checks),
            "maximum_spot_check_absolute_error": max_spot_error,
            "measured_table_sha256": digest(args.measured),
            "runner_sha256": next(iter({m["runner_sha256"] for m in metas})),
            "species_results_sha256": digest(args.output / "species_results.csv"),
            "global_null_sha256": digest(args.output / "global_null.csv"),
        },
        "six_species_used": False,
        "thirty_four_species_used": False,
        "environmental_variables_used": False,
        "parent_G1_reclassified": False,
        "shared_boundary_claim_allowed": False,
        "biological_prevalence_claim_allowed": False,
        "claim_ceiling": "Postoutcome exploratory RGFCA-only exact randomization evidence for a species-equal geographic-distance versus flower-colour-dissimilarity association. The BH species count is detectable-species count, not true prevalence; no shared-boundary, environment, causal, six-species, or 34-species claim is allowed."
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
