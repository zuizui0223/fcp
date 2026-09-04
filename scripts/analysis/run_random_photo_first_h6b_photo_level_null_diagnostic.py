#!/usr/bin/env python3
"""Diagnose the only nominally positive H6 sparse-frame sensitivity.

The completed H6 cell-label null swaps species-cell mean colour compositions. H6b
uses a stricter photo-level null: individual classifiable photo colour vectors are
permuted within species while cell assignments and exact per-cell photo counts stay
fixed, then species-cell mean compositions are rebuilt. This diagnostic can never
rescue any primary result.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h6b_photo_level_null_diagnostic_contract_v1.json"
MEASURED = ROOT / "data/derived/random_photo_first_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_measurement_result_v1.json"
CAPACITY = ROOT / "data/derived/random_photo_first_h4_within_species_capacity_v1.csv"
H6_RESULT = ROOT / "docs/supporting/random_photo_first_h6_species_specific_spatial_structure_result_v1.json"

OUT_RESULT = ROOT / "docs/supporting/random_photo_first_h6b_photo_level_null_diagnostic_result_v1.json"
OUT_SPECIES = ROOT / "data/derived/random_photo_first_h6b_target_species_v1.csv"
OUT_NULL = ROOT / "data/derived/random_photo_first_h6b_photo_level_null_v1.csv"

BIOLOGICAL_MORPHS = {"white", "yellow_orange", "red_pink", "blue_purple"}
PALETTE = [
    "flower_fraction_white", "flower_fraction_yellow", "flower_fraction_orange",
    "flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta",
    "flower_fraction_purple", "flower_fraction_blue", "flower_fraction_bronze",
]
COLOUR_COLS = ["soft_white", "soft_yellow_orange", "soft_red_pink", "soft_blue_purple"]
N_LON = 18
N_SINLAT = 9


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    dx = x - float(np.mean(x)); dy = y - float(np.mean(y))
    nx = float(np.sqrt(np.sum(dx * dx))); ny = float(np.sqrt(np.sum(dy * dy)))
    if nx <= 0.0 or ny <= 0.0: return float("nan")
    return float(np.clip(np.sum(dx * dy) / (nx * ny), -1.0, 1.0))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rankdata(np.asarray(x, dtype=float), method="average"), rankdata(np.asarray(y, dtype=float), method="average"))


def fixed_rank_unit(x: np.ndarray) -> np.ndarray:
    r = rankdata(np.asarray(x, dtype=float), method="average")
    r -= float(np.mean(r)); n = float(np.sqrt(np.sum(r * r)))
    if n <= 0.0: raise RuntimeError("constant geographic rank vector")
    return r / n


def spearman_against_fixed_unit(fixed: np.ndarray, y: np.ndarray) -> float:
    r = rankdata(np.asarray(y, dtype=float), method="average")
    r -= float(np.mean(r)); n = float(np.sqrt(np.sum(r * r)))
    if n <= 0.0: return float("nan")
    return float(np.clip(np.sum(fixed * (r / n)), -1.0, 1.0))


def stable_species_permutation(n: int, *, seed: int, permutation: int, species: str) -> np.ndarray:
    token = f"{seed}|{permutation}|{species}".encode("utf-8")
    child = int.from_bytes(hashlib.sha256(token).digest()[:8], "big", signed=False)
    return np.random.default_rng(child).permutation(n)


def cell_centroid(cell_id: int) -> tuple[float, float]:
    row = int(cell_id) // N_LON; col = int(cell_id) % N_LON
    lon = -180.0 + (col + 0.5) * (360.0 / N_LON)
    sin_lo = -1.0 + row * (2.0 / N_SINLAT); sin_hi = -1.0 + (row + 1) * (2.0 / N_SINLAT)
    lat = float(np.degrees(np.arcsin(np.clip(0.5 * (sin_lo + sin_hi), -1.0, 1.0))))
    return lat, lon


def great_circle_km(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    lat1 = np.radians(np.asarray(lat1, dtype=float)); lat2 = np.radians(np.asarray(lat2, dtype=float))
    lon1 = np.radians(np.asarray(lon1, dtype=float)); lon2 = np.radians(np.asarray(lon2, dtype=float))
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def jsd_pairs(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float); q = np.asarray(q, dtype=float); m = 0.5 * (p + q)
    t1 = np.zeros_like(p); t2 = np.zeros_like(q)
    k = p > 0; t1[k] = p[k] * np.log2(p[k] / m[k])
    k = q > 0; t2[k] = q[k] * np.log2(q[k] / m[k])
    return np.clip(0.5 * (t1.sum(axis=1) + t2.sum(axis=1)), 0.0, 1.0)


def add_soft_colour(photos: pd.DataFrame) -> pd.DataFrame:
    out = photos.copy()
    raw = out[PALETTE].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).all(): raise RuntimeError("non-finite classifiable palette fraction")
    soft = np.column_stack([
        raw[:, 0], raw[:, 1] + raw[:, 2] + raw[:, 8],
        raw[:, 3] + raw[:, 4] + raw[:, 5], raw[:, 6] + raw[:, 7],
    ])
    mass = soft.sum(axis=1)
    if np.any(mass <= 0.0): raise RuntimeError("zero biological colour mass")
    soft /= mass[:, None]
    for i, c in enumerate(COLOUR_COLS): out[c] = soft[:, i]
    return out


def prepare_species_photo_null(species: str, group: pd.DataFrame, *, permutations: int, seed: int) -> tuple[dict[str, Any], np.ndarray]:
    group = group.sort_values(["cell_id", "measurement_id"], kind="mergesort").reset_index(drop=True)
    cell_ids, counts = np.unique(group["cell_id"].to_numpy(dtype=int), return_counts=True)
    if len(cell_ids) < 3: raise RuntimeError(f"target species lost required cells: {species}")
    colour = group[COLOUR_COLS].to_numpy(dtype=float)
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    means = np.add.reduceat(colour, starts, axis=0) / counts[:, None]
    means /= means.sum(axis=1)[:, None]
    ii, jj = np.triu_indices(len(cell_ids), k=1)
    centers = np.array([cell_centroid(int(c)) for c in cell_ids], dtype=float)
    geo = great_circle_km(centers[ii, 0], centers[ii, 1], centers[jj, 0], centers[jj, 1])
    fixed = fixed_rank_unit(geo)
    observed_colour_distance = jsd_pairs(means[ii], means[jj])
    observed_rho = spearman(geo, observed_colour_distance)
    if not np.isfinite(observed_rho): raise RuntimeError(f"target species non-evaluable: {species}")

    null = np.empty(permutations, dtype=float)
    for p in range(permutations):
        order = stable_species_permutation(len(colour), seed=seed, permutation=p, species=species)
        permuted = colour[order]
        pmeans = np.add.reduceat(permuted, starts, axis=0) / counts[:, None]
        pmeans /= pmeans.sum(axis=1)[:, None]
        colour_distance = jsd_pairs(pmeans[ii], pmeans[jj])
        rho = spearman_against_fixed_unit(fixed, colour_distance)
        if not np.isfinite(rho): raise RuntimeError(f"photo-null non-evaluable for {species} permutation {p}")
        null[p] = rho

    record = {
        "species": species,
        "classifiable_photos": int(len(group)),
        "occupied_h1_cells": int(len(cell_ids)),
        "cell_pairs": int(len(ii)),
        "observed_geographic_colour_spearman_rho": float(observed_rho),
        "min_photos_per_cell": int(np.min(counts)),
        "max_photos_per_cell": int(np.max(counts)),
        "mean_photos_per_cell": float(np.mean(counts)),
    }
    return record, null


def main() -> int:
    contract = load_json(CONTRACT); measurement = load_json(MEASUREMENT_RESULT); h6 = load_json(H6_RESULT)
    if contract.get("status") != "diagnostic_frozen_after_h6_before_any_photo_level_null_retest": raise RuntimeError("H6b contract not frozen")
    if h6.get("decision") != "no_support_general_species_specific_spatial_structuring": raise RuntimeError("H6 primary decision drifted")
    target_h6 = next(x for x in h6["sensitivities"] if x["name"] == "n5_cells3_morphs2")
    if int(target_h6["evaluable_species"]) != 188 or float(target_h6["p_upper"]) != 0.001: raise RuntimeError("H6 target sensitivity drifted")
    if measurement.get("measurement_table_sha256") != sha256(MEASURED): raise RuntimeError("measurement SHA drifted")
    if int(measurement.get("classified_rows", -1)) != 10103: raise RuntimeError("classified denominator drifted")

    permutations = int(contract["stricter_null"]["permutations"]); seed = int(contract["stricter_null"]["seed"])
    if permutations != 999 or seed != 20260908: raise RuntimeError("H6b null settings drifted")

    measured = pd.read_csv(MEASURED)
    required = {"measurement_id", "species", "cell_id", "morph", *PALETTE}
    if not required.issubset(measured.columns): raise RuntimeError("measurement table lacks H6b fields")
    photos = measured.loc[measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS)].copy()
    photos["species"] = photos["species"].astype(str)
    photos = add_soft_colour(photos)

    capacity = pd.read_csv(CAPACITY); capacity["species"] = capacity["species"].astype(str)
    target_species = sorted(capacity.loc[(capacity["classifiable_photos"] >= 5) & (capacity["occupied_h1_cells"] >= 3) & (capacity["morph_levels"] >= 2), "species"].astype(str))
    reliable_species = set(capacity.loc[(capacity["classifiable_photos"] >= 10) & (capacity["occupied_h1_cells"] >= 5) & (capacity["morph_levels"] >= 2), "species"].astype(str))
    if len(target_species) != 188 or len(reliable_species) != 74: raise RuntimeError("H6b target/reliable frame capacity drifted")

    records: list[dict[str, Any]] = []; nulls: dict[str, np.ndarray] = {}
    for species in target_species:
        rec, null = prepare_species_photo_null(species, photos.loc[photos["species"] == species], permutations=permutations, seed=seed)
        records.append(rec); nulls[species] = null
    table = pd.DataFrame(records).sort_values("species", kind="mergesort").reset_index(drop=True)

    observed = float(table["observed_geographic_colour_spearman_rho"].mean())
    if not np.isclose(observed, float(target_h6["observed_mean_spearman_rho"]), atol=5e-4, rtol=0.0):
        raise RuntimeError("H6b observed statistic does not replay completed H6 target sensitivity")
    matrix = np.vstack([nulls[s] for s in table["species"].astype(str)])
    equal_null = matrix.mean(axis=0)
    p_equal = float((1 + np.count_nonzero(equal_null >= observed)) / (permutations + 1))

    pair_weights = table["cell_pairs"].to_numpy(dtype=float); pair_weights /= pair_weights.sum()
    observed_pair_weighted = float(np.sum(pair_weights * table["observed_geographic_colour_spearman_rho"].to_numpy(dtype=float)))
    pair_null = pair_weights @ matrix
    p_pair = float((1 + np.count_nonzero(pair_null >= observed_pair_weighted)) / (permutations + 1))

    reliable = table.loc[table["species"].astype(str).isin(reliable_species)].copy().reset_index(drop=True)
    if len(reliable) != 74: raise RuntimeError("H6b reliable frame lost species")
    reliable_obs = float(reliable["observed_geographic_colour_spearman_rho"].mean())
    if not np.isclose(reliable_obs, float(h6["primary"]["observed_mean_spearman_rho"]), atol=5e-4, rtol=0.0):
        raise RuntimeError("H6b reliable-frame observed statistic does not replay H6 primary")
    reliable_matrix = np.vstack([nulls[s] for s in reliable["species"].astype(str)])
    reliable_null = reliable_matrix.mean(axis=0)
    p_reliable = float((1 + np.count_nonzero(reliable_null >= reliable_obs)) / (permutations + 1))

    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    pd.DataFrame({
        "permutation": np.arange(permutations, dtype=int),
        "equal_species_mean_rho": equal_null,
        "pair_count_weighted_mean_rho": pair_null,
        "reliable_74_equal_species_mean_rho": reliable_null,
    }).to_csv(OUT_NULL, index=False, lineterminator="\n")

    target_survives = bool(observed > 0 and p_equal < float(contract["stricter_null"]["alpha"]))
    result = {
        "protocol": contract["protocol"],
        "status": "complete_h6b_diagnostic_evaluable",
        "introduced_after_h6_outcomes": True,
        "claim_role": contract["claim_role"],
        "target_sparse_frame": {
            "evaluable_species": 188,
            "observed_mean_spearman_rho": observed,
            "completed_h6_cell_label_null_p_upper": float(target_h6["p_upper"]),
            "photo_level_null_p_upper": p_equal,
            "photo_level_null_mean": float(np.mean(equal_null)),
            "photo_level_null_q025": float(np.quantile(equal_null, 0.025)),
            "photo_level_null_q975": float(np.quantile(equal_null, 0.975)),
            "survives_photo_level_null": target_survives,
        },
        "diagnostic_sensitivities": [
            {
                "name": "pair_count_weighted_global_mean",
                "observed_mean_spearman_rho": observed_pair_weighted,
                "photo_level_null_p_upper": p_pair,
                "photo_level_null_mean": float(np.mean(pair_null)),
                "nominal_positive": bool(observed_pair_weighted > 0 and p_pair < 0.05),
                "can_rescue": False,
            },
            {
                "name": "reliable_n10_cells5_morphs2_frame",
                "evaluable_species": 74,
                "observed_mean_spearman_rho": reliable_obs,
                "photo_level_null_p_upper": p_reliable,
                "photo_level_null_mean": float(np.mean(reliable_null)),
                "nominal_positive": bool(reliable_obs > 0 and p_reliable < 0.05),
                "can_rescue": False,
            },
        ],
        "diagnostic_decision": "sparse_equal_species_signal_survives_photo_count_preserving_null_but_is_not_information_weight_robust" if target_survives and p_pair >= 0.05 else ("sparse_signal_not_robust_to_photo_level_null" if not target_survives else "sparse_signal_survives_multiple_diagnostics"),
        "h6_primary_decision_unchanged": h6["decision"],
        "lineage": {"contract_sha256": sha256(CONTRACT), "measurement_table_sha256": sha256(MEASURED), "h6_result_sha256": sha256(H6_RESULT), "permutations": permutations, "seed": seed},
        "files": {"target_species": str(OUT_SPECIES.relative_to(ROOT)), "photo_level_null": str(OUT_NULL.relative_to(ROOT))},
        "claim_ceiling": contract["claim_ceiling"],
    }
    OUT_RESULT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
