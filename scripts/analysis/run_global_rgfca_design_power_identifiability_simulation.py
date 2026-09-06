#!/usr/bin/env python3
"""Post-outcome methodological simulations for the frozen RGFCA design.

This script never changes or recomputes an observed biological support decision.
It uses frozen photo geometry and outer schedules to quantify detectability,
classifiability attrition, partial sharing, outer-field convergence, and Holm cost.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from fcp_pipeline.global_barrier_field import equal_area_grid_centers
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule, stable_seed
from fcp_pipeline.global_rgfca_engine import canonical_colour_pool, prepare_sparse_outer_geometry

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_design_power_identifiability_simulation_contract_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
OBSERVED_OUTER = ROOT / "data/derived/global_rgfca_observed_outer_fields_v1.npz"

BLOCK_RESULTS = {
    "terrain_structure": ROOT / "docs/supporting/global_rgfca_earthenv_terrain_block_result_v1.json",
    "thermal_regime": ROOT / "docs/supporting/global_rgfca_thermal_regime_block_result_v1.json",
    "water_balance": ROOT / "docs/supporting/global_rgfca_water_balance_block_result_v1.json",
    "atmospheric_energy_dryness": ROOT / "docs/supporting/global_rgfca_atmospheric_energy_dryness_block_result_v1.json",
    "edaphic_regime": ROOT / "docs/supporting/global_rgfca_edaphic_regime_block_result_v1.json",
}
EARTH_RADIUS_KM = 6371.0088


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--s1-only-frame", choices=["actual_classifiable", "synthetic_complete_measured"], default=None)
    return p.parse_args()


def _bool_flag(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.casefold().isin({"true", "1", "yes"})


def _dummy_colour_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame[["photo_id", "species", "latitude", "longitude"]].copy()
    out["colour_white"] = 1.0
    out["colour_yellow_orange"] = 0.0
    out["colour_red_pink"] = 0.0
    out["colour_blue_purple"] = 0.0
    return out


def load_geometry_frame(frame_id: str) -> pd.DataFrame:
    cols = ["photo_id", "species", "latitude", "longitude", "global_classifiable"]
    raw = pd.read_csv(MEASURED, usecols=cols)
    raw["species"] = raw["species"].astype(str)
    raw["latitude"] = pd.to_numeric(raw["latitude"], errors="coerce")
    raw["longitude"] = pd.to_numeric(raw["longitude"], errors="coerce")
    raw = raw[np.isfinite(raw["latitude"]) & np.isfinite(raw["longitude"])].copy()
    if frame_id == "actual_classifiable":
        raw = raw[_bool_flag(raw["global_classifiable"])].copy()
    elif frame_id != "synthetic_complete_measured":
        raise ValueError(frame_id)
    counts = raw.groupby("species", observed=True).size()
    eligible = set(counts[counts >= 40].index.astype(str))
    raw = raw[raw["species"].isin(eligible)].copy()
    if raw["photo_id"].duplicated().any():
        raise RuntimeError(f"{frame_id}: duplicate photo_id")
    if raw["species"].nunique() < 250:
        raise RuntimeError(f"{frame_id}: fewer than 250 eligible species")
    return _dummy_colour_frame(raw)


def _xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat_r = np.deg2rad(lat)
    lon_r = np.deg2rad(lon)
    c = np.cos(lat_r)
    return np.column_stack([c * np.cos(lon_r), c * np.sin(lon_r), np.sin(lat_r)])


def _unit_vectors(rng: np.random.Generator, n: int) -> np.ndarray:
    x = rng.normal(size=(int(n), 3))
    length = np.linalg.norm(x, axis=1)
    while np.any(length == 0):
        bad = length == 0
        x[bad] = rng.normal(size=(int(np.sum(bad)), 3))
        length = np.linalg.norm(x, axis=1)
    return x / length[:, None]


def _species_index_by_row(pool) -> np.ndarray:
    out = np.empty(len(pool.photo_ids), dtype=np.int32)
    for i, (start, stop) in enumerate(pool.species_slices):
        out[start:stop] = i
    return out


def prepare_geometries(frame_id: str) -> tuple[object, object, list[object], np.ndarray, np.ndarray, np.ndarray]:
    frame = load_geometry_frame(frame_id)
    pool = canonical_colour_pool(frame)
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=200,
        species_per_outer=250,
        photos_per_species=20,
        minimum_pool_photos_per_species=40,
        species_seed=2026090401,
        photo_master_seed=2026090402,
    )
    grid = equal_area_grid_centers(36, 18)
    geometries = []
    aggregate_opp = np.zeros(grid.n_cells, dtype=np.float64)
    for outer in range(200):
        geom = prepare_sparse_outer_geometry(
            pool,
            schedule.outer_photo_ids[outer],
            schedule.outer_species[outer],
            grid=grid,
            k=3,
            kernel_km=500.0,
            cutoff_multiplier=3.0,
            minimum_distinct_species=5,
        )
        geometries.append(geom)
        aggregate_opp[geom.evaluable] += geom.opportunity[geom.evaluable]
    xyz = _xyz(pool.latitude, pool.longitude)
    row_species = _species_index_by_row(pool)
    return pool, schedule, geometries, aggregate_opp, xyz, row_species


def simulate_latent_batch(
    *,
    frame_id: str,
    xyz: np.ndarray,
    row_species: np.ndarray,
    n_species: int,
    shared_fraction: float,
    amplitude: float,
    replicate_indices: np.ndarray,
    master_seed: int,
) -> np.ndarray:
    y = np.empty((len(replicate_indices), len(xyz)), dtype=np.float32)
    n_shared = int(round(float(shared_fraction) * int(n_species)))
    n_shared = min(max(n_shared, 0), int(n_species))
    for r, replicate in enumerate(replicate_indices):
        rng = np.random.default_rng(
            stable_seed(master_seed, frame_id, f"shared={shared_fraction:.6f}", f"amp={amplitude:.6f}", int(replicate))
        )
        if amplitude == 0.0:
            y[r] = rng.normal(size=len(xyz)).astype(np.float32)
            continue
        normals = _unit_vectors(rng, n_species)
        common = _unit_vectors(rng, 1)[0]
        if n_shared:
            shared_idx = rng.choice(n_species, size=n_shared, replace=False)
            normals[shared_idx] = common
        response_sign = rng.choice(np.array([-1.0, 1.0]), size=n_species)
        dot = np.einsum("ij,ij->i", xyz, normals[row_species])
        signed_km = np.arcsin(np.clip(dot, -1.0, 1.0)) * EARTH_RADIUS_KM
        signal = float(amplitude) * response_sign[row_species] * np.tanh(signed_km / 500.0)
        y[r] = (signal + rng.normal(size=len(xyz))).astype(np.float32)
    return y


def _rank_edge_matrix(raw: np.ndarray, slices: tuple[tuple[int, int], ...]) -> np.ndarray:
    out = np.empty(raw.shape, dtype=np.float64)
    for start, stop in slices:
        n = stop - start
        ranks = rankdata(raw[:, start:stop], method="average", axis=1)
        out[:, start:stop] = (ranks - 0.5) / float(n)
    return out


def concentration_batch(latent: np.ndarray, geometries: list[object], aggregate_opp: np.ndarray) -> np.ndarray:
    n_rep = latent.shape[0]
    n_cells = len(aggregate_opp)
    aggregate_num = np.zeros((n_rep, n_cells), dtype=np.float64)
    for geom in geometries:
        left = geom.edge_nodes[:, 0]
        right = geom.edge_nodes[:, 1]
        raw = np.abs(latent[:, left] - latent[:, right])
        ranks = _rank_edge_matrix(raw, geom.edge_species_slices)
        numerator = np.asarray(geom.weighted_kernel.T.dot(ranks.T).T)
        keep = geom.evaluable
        aggregate_num[:, keep] += numerator[:, keep]
    keep = aggregate_opp > 0
    w = aggregate_opp[keep]
    field = aggregate_num[:, keep] / w[None, :]
    mean = np.sum(field * w[None, :], axis=1) / np.sum(w)
    return np.sum(np.square(field - mean[:, None]) * w[None, :], axis=1) / np.sum(w)


def run_s1(frame_id: str, contract: dict[str, object]) -> tuple[pd.DataFrame, dict[str, object]]:
    pool, schedule, geometries, aggregate_opp, xyz, row_species = prepare_geometries(frame_id)
    s1 = contract["S1_spatial_design_simulation"]
    base_seed = int(s1["scenario_grid"]["seed"])
    null_n = int(s1["calibration"]["null_replicates"])
    batch_size = 100
    null_parts = []
    for start in range(0, null_n, batch_size):
        idx = np.arange(start, min(start + batch_size, null_n), dtype=int)
        latent = simulate_latent_batch(
            frame_id=frame_id,
            xyz=xyz,
            row_species=row_species,
            n_species=len(pool.species_labels),
            shared_fraction=0.0,
            amplitude=0.0,
            replicate_indices=idx,
            master_seed=base_seed,
        )
        null_parts.append(concentration_batch(latent, geometries, aggregate_opp))
    null = np.concatenate(null_parts)
    calibration = null[: null_n // 2]
    evaluation = null[null_n // 2 :]
    threshold = float(np.quantile(calibration, 0.95, method="higher"))
    fpr = float(np.mean(evaluation >= threshold))
    null_mean = float(np.mean(calibration))

    rows = []
    fractions = [float(x) for x in s1["scenario_grid"]["shared_species_fraction"]]
    amplitudes = [float(x) for x in s1["scenario_grid"]["signal_amplitude"]]
    reps = int(s1["scenario_grid"]["simulation_replicates_per_nonredundant_scenario"])
    for amplitude in amplitudes:
        scenario_fractions = [0.0] if amplitude == 0.0 else fractions
        for shared in scenario_fractions:
            idx = np.arange(reps, dtype=int)
            latent = simulate_latent_batch(
                frame_id=frame_id,
                xyz=xyz,
                row_species=row_species,
                n_species=len(pool.species_labels),
                shared_fraction=shared,
                amplitude=amplitude,
                replicate_indices=idx,
                master_seed=stable_seed(base_seed, "alternative"),
            )
            conc = concentration_batch(latent, geometries, aggregate_opp)
            ratio = conc / null_mean
            rows.append(
                {
                    "frame": frame_id,
                    "eligible_species": len(pool.species_labels),
                    "pool_rows": len(pool.photo_ids),
                    "shared_species_fraction": shared,
                    "signal_amplitude": amplitude,
                    "simulation_replicates": reps,
                    "synthetic_null_threshold_95": threshold,
                    "synthetic_null_calibration_mean": null_mean,
                    "power_ge_threshold": float(np.mean(conc >= threshold)),
                    "median_concentration_ratio_to_null_mean": float(np.median(ratio)),
                    "q10_concentration_ratio_to_null_mean": float(np.quantile(ratio, 0.10)),
                    "q90_concentration_ratio_to_null_mean": float(np.quantile(ratio, 0.90)),
                }
            )
    summary = {
        "frame": frame_id,
        "eligible_species": len(pool.species_labels),
        "pool_rows": len(pool.photo_ids),
        "schedule": {
            "outer_realizations": schedule.n_outer,
            "species_per_outer": schedule.species_per_outer,
            "photos_per_species": schedule.photos_per_species,
            "minimum_species_inclusions": int(np.min(schedule.species_inclusion_counts)),
            "maximum_species_inclusions": int(np.max(schedule.species_inclusion_counts)),
            "photo_max_inclusion_imbalance": int(schedule.photo_max_inclusion_imbalance),
        },
        "null": {
            "replicates": null_n,
            "calibration_replicates": len(calibration),
            "evaluation_replicates": len(evaluation),
            "threshold_95_higher": threshold,
            "calibration_mean": null_mean,
            "evaluation_type1_fraction": fpr,
        },
    }
    return pd.DataFrame(rows), summary


def _rowwise_corr(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.full(a.shape[0], np.nan, dtype=float)
    for i in range(a.shape[0]):
        keep = np.isfinite(a[i]) & np.isfinite(b[i])
        if np.sum(keep) < 3:
            continue
        x = a[i, keep]
        y = b[i, keep]
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        out[i] = float(np.corrcoef(x, y)[0, 1])
    return out


def _bootstrap_consensus(counts: np.ndarray, num: np.ndarray, opp: np.ndarray) -> np.ndarray:
    agg_num = counts @ num
    agg_opp = counts @ opp
    out = np.full(agg_num.shape, np.nan, dtype=float)
    keep = agg_opp > 0
    out[keep] = agg_num[keep] / agg_opp[keep]
    return out


def run_s2(contract: dict[str, object]) -> pd.DataFrame:
    data = np.load(OBSERVED_OUTER)
    fields = np.asarray(data["observed_outer_fields"], dtype=float)
    opp = np.asarray(data["observed_outer_opportunities"], dtype=float)
    if fields.shape != opp.shape or fields.shape[0] != 200:
        raise RuntimeError("unexpected observed outer-field shape")
    evaluable = np.isfinite(fields) & np.isfinite(opp) & (opp > 0)
    num = np.where(evaluable, fields * opp, 0.0)
    eff_opp = np.where(evaluable, opp, 0.0)
    full_num = np.sum(num, axis=0)
    full_opp = np.sum(eff_opp, axis=0)
    full = np.full(fields.shape[1], np.nan, dtype=float)
    keep = full_opp > 0
    full[keep] = full_num[keep] / full_opp[keep]

    s2 = contract["S2_outer_stability_simulation"]
    reps = int(s2["bootstrap_replicates"])
    rng = np.random.default_rng(int(s2["seed"]))
    p = np.full(200, 1.0 / 200.0)
    rows = []
    for n in [int(x) for x in s2["outer_counts"]]:
        c1 = rng.multinomial(n, p, size=reps)
        c2 = rng.multinomial(n, p, size=reps)
        f1 = _bootstrap_consensus(c1, num, eff_opp)
        f2 = _bootstrap_consensus(c2, num, eff_opp)
        pair_r = _rowwise_corr(f1, f2)
        full_stack = np.broadcast_to(full, f1.shape)
        to_full_r = _rowwise_corr(f1, full_stack)
        pair_r = pair_r[np.isfinite(pair_r)]
        to_full_r = to_full_r[np.isfinite(to_full_r)]
        rows.append(
            {
                "outer_count": n,
                "bootstrap_replicates": reps,
                "independent_pair_median_r": float(np.median(pair_r)),
                "independent_pair_q10_r": float(np.quantile(pair_r, 0.10)),
                "independent_pair_q90_r": float(np.quantile(pair_r, 0.90)),
                "independent_pair_prob_r_ge_0_9": float(np.mean(pair_r >= 0.9)),
                "independent_pair_prob_r_ge_0_95": float(np.mean(pair_r >= 0.95)),
                "to_frozen_200_median_r": float(np.median(to_full_r)),
                "to_frozen_200_q10_r": float(np.quantile(to_full_r, 0.10)),
                "to_frozen_200_q90_r": float(np.quantile(to_full_r, 0.90)),
                "to_frozen_200_prob_r_ge_0_9": float(np.mean(to_full_r >= 0.9)),
                "to_frozen_200_prob_r_ge_0_95": float(np.mean(to_full_r >= 0.95)),
            }
        )
    return pd.DataFrame(rows)


def holm_adjust_matrix(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = p.shape[1]
    order = np.argsort(p, axis=1)
    sorted_p = np.take_along_axis(p, order, axis=1)
    factors = np.arange(m, 0, -1, dtype=float)[None, :]
    adjusted_sorted = np.maximum.accumulate(sorted_p * factors, axis=1)
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)
    out = np.empty_like(adjusted_sorted)
    rows = np.arange(p.shape[0])[:, None]
    out[rows, order] = adjusted_sorted
    return out


def run_s3_main_blocks(contract: dict[str, object]) -> tuple[pd.DataFrame, dict[str, object]]:
    names = list(BLOCK_RESULTS)
    null_mean = []
    sigma = []
    for name in names:
        x = json.loads(BLOCK_RESULTS[name].read_text(encoding="utf-8"))
        q025 = float(x["null_q025"])
        q975 = float(x["null_q975"])
        null_mean.append(float(x.get("null_mean", 0.0)))
        sigma.append((q975 - q025) / (2.0 * norm.ppf(0.975)))
    null_mean = np.asarray(null_mean)
    sigma = np.asarray(sigma)
    s3_family = next(x for x in contract["S3_multiplicity_detectability"]["families"] if x["id"] == "five_environmental_parent_blocks")
    effects = [float(x) for x in s3_family["effect_grid"]]
    reps = int(contract["S3_multiplicity_detectability"]["replicates"])
    base_seed = int(contract["S3_multiplicity_detectability"]["seed"])
    rows = []
    for corr in [float(x) for x in contract["S3_multiplicity_detectability"]["correlation_scenarios"]]:
        corr_matrix = np.full((5, 5), corr, dtype=float)
        np.fill_diagonal(corr_matrix, 1.0)
        cov = corr_matrix * sigma[:, None] * sigma[None, :]
        for target_index, target in enumerate(names):
            rng = np.random.default_rng(stable_seed(base_seed, "five-block", corr, target))
            noise = rng.multivariate_normal(np.zeros(5), cov, size=reps)
            for effect in effects:
                mean = null_mean.copy()
                mean[target_index] += effect
                stat = mean[None, :] + noise
                z = (stat - null_mean[None, :]) / sigma[None, :]
                p_upper = norm.sf(z)
                p_holm = holm_adjust_matrix(p_upper)
                target_reject = p_holm[:, target_index] < 0.05
                one_test = p_upper[:, target_index] < 0.05
                any_reject = np.any(p_holm < 0.05, axis=1)
                rows.append(
                    {
                        "family": "five_environmental_parent_blocks",
                        "correlation": corr,
                        "target_block": target,
                        "true_target_effect": effect,
                        "replicates": reps,
                        "target_holm_power": float(np.mean(target_reject)),
                        "family_any_discovery_probability": float(np.mean(any_reject)),
                        "target_one_test_power": float(np.mean(one_test)),
                        "holm_power_cost": float(np.mean(one_test) - np.mean(target_reject)),
                    }
                )
    meta = {
        "blocks": names,
        "null_mean": {name: float(value) for name, value in zip(names, null_mean)},
        "gaussian_sigma_from_reported_null_95_interval": {name: float(value) for name, value in zip(names, sigma)},
        "approximation": "Gaussian null approximation fixed from each raw block's reported 2.5/97.5% null quantiles; no observed block is removed from the five-test family.",
    }
    return pd.DataFrame(rows), meta


def summarize_minimum_shared(power: pd.DataFrame) -> list[dict[str, object]]:
    out = []
    for (frame, amp), group in power[power["signal_amplitude"] > 0].groupby(["frame", "signal_amplitude"]):
        g = group.sort_values("shared_species_fraction")
        hit = g[g["power_ge_threshold"] >= 0.8]
        out.append(
            {
                "frame": str(frame),
                "signal_amplitude": float(amp),
                "minimum_shared_fraction_for_power_ge_0_8": None if hit.empty else float(hit.iloc[0]["shared_species_fraction"]),
            }
        )
    return out


def main() -> int:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_before_any_new_design_power_simulation_result":
        raise RuntimeError("simulation contract is not in frozen pre-result state")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frames = [args.s1_only_frame] if args.s1_only_frame else ["actual_classifiable", "synthetic_complete_measured"]
    s1_frames = []
    s1_tables = []
    for frame_id in frames:
        table, summary = run_s1(frame_id, contract)
        s1_tables.append(table)
        s1_frames.append(summary)
    s1 = pd.concat(s1_tables, ignore_index=True)
    s1.to_csv(args.output_dir / "global_rgfca_design_power_s1_v1.csv", index=False)

    s2 = run_s2(contract)
    s2.to_csv(args.output_dir / "global_rgfca_outer_stability_s2_v1.csv", index=False)

    s3, s3_meta = run_s3_main_blocks(contract)
    s3.to_csv(args.output_dir / "global_rgfca_environmental_holm_power_s3_v1.csv", index=False)

    payload = {
        "protocol": contract["protocol"],
        "status": "complete_postoutcome_methodological_simulation_no_observed_decision_change",
        "role": contract["role"],
        "S1": {
            "frames": s1_frames,
            "minimum_shared_fraction": summarize_minimum_shared(s1),
            "scenario_rows": int(len(s1)),
        },
        "S2": {
            "rows": s2.to_dict(orient="records"),
            "frozen_observed_gate_reference": {"r_100_vs_200": 0.8809461021778187, "threshold": 0.9},
        },
        "S3": {
            "implemented_family": "five_environmental_parent_blocks",
            "meta": s3_meta,
            "rows": int(len(s3)),
            "ten_interaction_family_status": "deferred_until_exact_interaction_null_distribution_is_recovered_from_immutable_artifact; no approximation from observed pair p-values is used",
        },
        "decision_guard": {
            "G1_primary_reclassified": false,
            "environmental_parent_blocks_reclassified": false,
            "interaction_family_reclassified": false,
            "variable_level_decomposition_opened": false,
        },
    }
    (args.output_dir / "global_rgfca_design_power_identifiability_simulation_result_v1.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
