#!/usr/bin/env python3
"""Finalize exact matched-background recovery and the frozen paired RGFCA diagnostic."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from fcp_pipeline.photo_first_measurement import REFERENCE_RGB

EARTH_RADIUS_KM = 6371.0088
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "docs/supporting/global_rgfca_within_species_spatial_background_control_contract_v2.json"
DEFAULT_AMENDMENT = ROOT / "docs/supporting/global_rgfca_background_control_postoutcome_status_amendment_v2a.json"
DEFAULT_MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.casefold().isin({"true", "1"})


def pairwise_geo_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(np.asarray(lat, float))
    lonr = np.deg2rad(np.asarray(lon, float))
    c = np.cos(latr)
    xyz = np.column_stack([c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def pairwise_jsd_matrix(prob: np.ndarray) -> np.ndarray:
    p = np.asarray(prob, float)
    if p.ndim != 2 or np.any(~np.isfinite(p)) or np.any(p < 0):
        raise ValueError("palette probabilities must be finite nonnegative rows")
    mass = p.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("palette probabilities require positive row mass")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    return np.clip(0.5 * (ka + kb), 0.0, 1.0)


def direct_rho(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if np.ptp(y) <= 1e-15:
        return 0.0
    if np.ptp(x) <= 1e-12:
        return float("nan")
    value = float(spearmanr(x, y).statistic)
    return value if np.isfinite(value) else float("nan")


def seed_for(master: int, species: str, permutation: int) -> int:
    payload = f"{int(master)}\x1f{species}\x1f{int(permutation)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def observed_and_null_differential(
    geo: np.ndarray,
    differential: np.ndarray,
    *,
    species: str,
    master_seed: int,
    permutations: int,
    batch_size: int = 64,
) -> tuple[float, np.ndarray]:
    n = geo.shape[0]
    if geo.shape != (n, n) or differential.shape != (n, n):
        raise ValueError("pairwise matrices must be square and equal")
    ui, uj = np.triu_indices(n, k=1)
    gv = geo[ui, uj]
    dv = differential[ui, uj]
    observed = direct_rho(gv, dv)
    if not np.isfinite(observed):
        raise RuntimeError(f"non-evaluable geographic distances for species {species}")
    if np.ptp(dv) <= 1e-15:
        return 0.0, np.zeros(int(permutations), dtype=float)

    gr = rankdata(gv, method="average").astype(float)
    dr = rankdata(dv, method="average").astype(float)
    gc = gr - gr.mean()
    dc = dr - dr.mean()
    denom = float(np.sqrt(np.dot(gc, gc) * np.dot(dc, dc)))
    if denom <= 0:
        return 0.0, np.zeros(int(permutations), dtype=float)
    rank_matrix = np.zeros((n, n), dtype=float)
    rank_matrix[ui, uj] = dr
    rank_matrix[uj, ui] = dr
    null = np.empty(int(permutations), dtype=float)
    for start in range(0, int(permutations), int(batch_size)):
        stop = min(int(permutations), start + int(batch_size))
        perms = np.stack([
            np.random.default_rng(seed_for(master_seed, species, p)).permutation(n)
            for p in range(start, stop)
        ])
        values = rank_matrix[perms[:, ui], perms[:, uj]]
        null[start:stop] = (values @ gc) / denom
    return observed, null


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recovery-dir", type=Path, required=True)
    ap.add_argument("--measured", type=Path, default=DEFAULT_MEASURED)
    ap.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    ap.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    amendment = json.loads(args.amendment.read_text(encoding="utf-8"))
    permutations = int(contract["joint_exact_randomization_null"]["permutations"])
    master_seed = int(contract["joint_exact_randomization_null"]["master_seed"])
    if amendment.get("status") != "frozen_after_flower_only_omnibus_outcome_but_before_any_background_colour_recovery":
        raise RuntimeError("background timing amendment drifted")
    if sha256_file(args.measured) != contract["matched_frame"]["measured_table_sha256"]:
        raise RuntimeError("measured table differs from frozen matched control")

    files = sorted(args.recovery_dir.rglob("recovery_s*_p*.csv"))
    if len(files) != 128:
        raise RuntimeError(f"expected 128 background recovery partitions, found {len(files)}")
    recovery = pd.concat([pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in files], ignore_index=True)
    if len(recovery) != 21424 or recovery["measurement_id"].nunique() != 21424:
        raise RuntimeError("recovery partitions do not cover exactly 21,424 unique matched photos")

    exact = recovery["recovery_status"].eq("exact_matched_background_recovered")
    status_counts = recovery["recovery_status"].value_counts().sort_index().to_dict()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    recovery.to_csv(out / "background_recovery_rows_v2.csv", index=False, lineterminator="\n")
    if not bool(exact.all()):
        summary = {
            "protocol": contract["protocol"],
            "status": "not_evaluable_incomplete_exact_background_recovery",
            "inferential_role": "postoutcome_falsification_diagnostic",
            "matched_rows_required": 21424,
            "exact_rows": int(exact.sum()),
            "failed_rows": int((~exact).sum()),
            "status_counts": {str(k): int(v) for k, v in status_counts.items()},
            "primary_statistic_computed": False,
            "replacement_photos_used": False,
            "denominator_adapted_after_recovery": False,
            "claim": "No flower-background diagnostic may be interpreted because exact same-image/ROI reproduction was incomplete.",
            "lineage": {
                "contract_sha256": sha256_file(args.contract),
                "amendment_sha256": sha256_file(args.amendment),
                "measured_sha256": sha256_file(args.measured),
                "recovery_rows_sha256": sha256_file(out / "background_recovery_rows_v2.csv"),
            }
        }
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 0

    measured = pd.read_csv(args.measured, dtype={"measurement_id": str}).fillna("")
    cf = measured.loc[bool_series(measured["global_classifiable"])].copy()
    counts = cf.groupby("species", observed=True).size()
    frame = cf.loc[cf["species"].isin(counts[counts >= 40].index)].copy()
    if len(frame) != 21424 or frame["species"].nunique() != 369:
        raise RuntimeError("matched statistical frame drifted")
    merged = frame.merge(recovery, on="measurement_id", how="inner", validate="one_to_one")
    if len(merged) != 21424:
        raise RuntimeError("matched flower/background join is incomplete")

    names = tuple(REFERENCE_RGB)
    flower_cols = [f"palette_count_{name}" for name in names]
    background_cols = [f"background_palette_count_{name}" for name in names]
    flower_counts = merged[flower_cols].apply(pd.to_numeric, errors="raise").to_numpy(float)
    background_counts = merged[background_cols].apply(pd.to_numeric, errors="raise").to_numpy(float)
    if np.any(flower_counts.sum(axis=1) <= 0) or np.any(background_counts.sum(axis=1) <= 0):
        raise RuntimeError("matched twelve-anchor vectors contain zero mass")
    flower_prob = flower_counts / flower_counts.sum(axis=1, keepdims=True)
    background_prob = background_counts / background_counts.sum(axis=1, keepdims=True)

    species_rows = []
    all_null = []
    spot_rows = []
    sorted_species = sorted(merged["species"].astype(str).unique())
    spot_species = set(sorted_species[i] for i in np.linspace(0, len(sorted_species)-1, 12, dtype=int))
    for species, group in merged.groupby("species", sort=True, observed=True):
        idx = group.index.to_numpy()
        # merged retains original index only implicitly; map by measurement IDs to row positions.
        positions = merged.index[merged["species"].eq(species)].to_numpy()
        lat = pd.to_numeric(merged.loc[positions, "latitude"], errors="raise").to_numpy(float)
        lon = pd.to_numeric(merged.loc[positions, "longitude"], errors="raise").to_numpy(float)
        fp = flower_prob[positions]
        bp = background_prob[positions]
        geo = pairwise_geo_km(lat, lon)
        fj = pairwise_jsd_matrix(fp)
        bj = pairwise_jsd_matrix(bp)
        diff = fj - bj
        ui, uj = np.triu_indices(len(positions), k=1)
        flower_rho = direct_rho(geo[ui, uj], fj[ui, uj])
        background_rho = direct_rho(geo[ui, uj], bj[ui, uj])
        observed, null = observed_and_null_differential(
            geo, diff, species=str(species), master_seed=master_seed, permutations=permutations
        )
        if not (np.isfinite(flower_rho) and np.isfinite(background_rho) and np.isfinite(observed) and np.isfinite(null).all()):
            raise RuntimeError(f"non-finite matched statistic for species {species}")
        species_rows.append({
            "species": str(species),
            "photos": int(len(positions)),
            "flower12_rho": float(flower_rho),
            "background12_rho": float(background_rho),
            "differential_rho": float(observed),
        })
        all_null.append(null)
        if str(species) in spot_species:
            for p in (0, permutations - 1):
                perm = np.random.default_rng(seed_for(master_seed, str(species), p)).permutation(len(positions))
                direct = direct_rho(geo[ui, uj], diff[np.ix_(perm, perm)][ui, uj])
                spot_rows.append({
                    "species": str(species),
                    "permutation": int(p),
                    "fast_rho": float(null[p]),
                    "direct_scipy_rho": float(direct),
                    "absolute_error": abs(float(null[p]) - float(direct)),
                })

    species_frame = pd.DataFrame(species_rows).sort_values("species", kind="mergesort").reset_index(drop=True)
    if len(species_frame) != 369:
        raise RuntimeError("species result denominator is not 369")
    null_matrix = np.stack(all_null, axis=0)
    if null_matrix.shape != (369, permutations):
        raise RuntimeError("species-by-permutation null matrix shape drifted")
    global_null = null_matrix.mean(axis=0)
    observed_global = float(species_frame["differential_rho"].mean())
    p_upper = (1 + int(np.count_nonzero(global_null >= observed_global - 1e-15))) / (permutations + 1)
    spot = pd.DataFrame(spot_rows)
    max_spot_error = float(spot["absolute_error"].max()) if len(spot) else float("nan")
    if not np.isfinite(max_spot_error) or max_spot_error > 2e-12:
        raise RuntimeError(f"independent SciPy differential spot-check failed: {max_spot_error}")

    species_frame.to_csv(out / "species_results_v2.csv", index=False, lineterminator="\n")
    pd.DataFrame({"permutation": np.arange(permutations), "equal_species_mean_differential_rho": global_null}).to_csv(
        out / "global_null_v2.csv", index=False, lineterminator="\n"
    )
    spot.to_csv(out / "independent_spot_checks_v2.csv", index=False, lineterminator="\n")

    flower_mean = float(species_frame["flower12_rho"].mean())
    background_mean = float(species_frame["background12_rho"].mean())
    summary = {
        "protocol": contract["protocol"],
        "status": "complete_postoutcome_background_falsification_diagnostic",
        "inferential_role": "postoutcome_falsification_diagnostic_not_confirmatory",
        "known_flower_only_four_group_result": amendment["known_flower_only_result_before_background_control"],
        "matched_rows": 21424,
        "matched_species": 369,
        "background_recovery_exact_rows": 21424,
        "primary": {
            "statistic": "equal_species_mean_distance_vs_flower12_minus_background12_jsd_spearman_rho",
            "observed_mean_differential_rho": observed_global,
            "null_mean": float(np.mean(global_null)),
            "null_sd": float(np.std(global_null, ddof=1)),
            "null_q025": float(np.quantile(global_null, 0.025)),
            "null_q975": float(np.quantile(global_null, 0.975)),
            "p_upper": float(p_upper),
            "supported_at_0_05": bool(p_upper < 0.05),
        },
        "mandatory_negative_control": {
            "equal_species_mean_flower12_rho": flower_mean,
            "equal_species_mean_background12_rho": background_mean,
            "flower_minus_background_mean_rho": flower_mean - background_mean,
            "background_comparable_or_larger_than_flower": bool(background_mean >= flower_mean),
        },
        "species_descriptive": {
            "median_differential_rho": float(species_frame["differential_rho"].median()),
            "positive_differential_fraction": float((species_frame["differential_rho"] > 0).mean()),
        },
        "interpretation": (
            "Flower-specific spatial differentiation exceeds the frozen matched local-background diagnostic, but this remains postoutcome exploratory and is not sharedness or causal evidence."
            if p_upper < 0.05
            else "The matched background-adjusted diagnostic is unsupported; the flower-only p=.001 result must not be promoted as biological evidence of flower-specific spatial organization."
        ),
        "sharedness_main_line_unchanged": True,
        "g1_reclassified": False,
        "verification": {
            "permutations": permutations,
            "joint_flower_background_vertex_permutation": True,
            "direct_scipy_spot_checks": int(len(spot)),
            "maximum_spot_check_absolute_error": max_spot_error,
        },
        "lineage": {
            "contract_sha256": sha256_file(args.contract),
            "amendment_sha256": sha256_file(args.amendment),
            "measured_sha256": sha256_file(args.measured),
            "recovery_rows_sha256": sha256_file(out / "background_recovery_rows_v2.csv"),
            "species_results_sha256": sha256_file(out / "species_results_v2.csv"),
            "global_null_sha256": sha256_file(out / "global_null_v2.csv"),
        }
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
