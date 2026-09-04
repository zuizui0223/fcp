#!/usr/bin/env python3
"""Run exploratory H6: species-specific spatial structuring of flower colour.

H6 asks whether colour composition tends to diverge with geographic separation
within species even though the preregistered H1 found no shared global boundary.
It is explicitly exploratory and cannot rescue H1/H2/H4a/H5.
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
CONTRACT = ROOT / "docs/supporting/random_photo_first_h6_species_specific_spatial_structure_contract_v1.json"
MEASURED = ROOT / "data/derived/random_photo_first_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_measurement_result_v1.json"
CAPACITY = ROOT / "data/derived/random_photo_first_h4_within_species_capacity_v1.csv"
CAPACITY_RESULT = ROOT / "docs/supporting/random_photo_first_h4_within_species_capacity_v1.json"
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"
H2_RESULT = ROOT / "docs/supporting/random_photo_first_h2_result_v1.json"
H4A_RESULT = ROOT / "docs/supporting/random_photo_first_h4a_result_v1.json"
H5_RESULT = ROOT / "docs/supporting/random_photo_first_h5_scale_moderation_result_v1.json"

OUT_RESULT = ROOT / "docs/supporting/random_photo_first_h6_species_specific_spatial_structure_result_v1.json"
OUT_SPECIES = ROOT / "data/derived/random_photo_first_h6_species_primary_v1.csv"
OUT_NULL = ROOT / "data/derived/random_photo_first_h6_null_primary_v1.csv"
OUT_SENS = ROOT / "data/derived/random_photo_first_h6_sensitivities_v1.csv"

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


def jsd_bits(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float); q = np.asarray(q, dtype=float)
    p = p / p.sum(); q = q / q.sum(); m = 0.5 * (p + q)
    out = 0.0
    kp = p > 0
    if kp.any(): out += 0.5 * float(np.sum(p[kp] * np.log2(p[kp] / m[kp])))
    kq = q > 0
    if kq.any(): out += 0.5 * float(np.sum(q[kq] * np.log2(q[kq] / m[kq])))
    return float(np.clip(out, 0.0, 1.0))


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
    r = r - float(np.mean(r)); n = float(np.sqrt(np.sum(r * r)))
    if n <= 0.0: raise RuntimeError("constant geographic rank vector")
    return r / n


def spearman_against_fixed_unit(fixed: np.ndarray, y: np.ndarray) -> float:
    r = rankdata(np.asarray(y, dtype=float), method="average")
    r = r - float(np.mean(r)); n = float(np.sqrt(np.sum(r * r)))
    if n <= 0.0: return float("nan")
    return float(np.clip(np.sum(fixed * (r / n)), -1.0, 1.0))


def stable_species_permutation(n: int, *, seed: int, permutation: int, species: str) -> np.ndarray:
    token = f"{seed}|{permutation}|{species}".encode("utf-8")
    child = int.from_bytes(hashlib.sha256(token).digest()[:8], "big", signed=False)
    return np.random.default_rng(child).permutation(n)


def cell_centroid(cell_id: int) -> tuple[float, float]:
    row = int(cell_id) // N_LON; col = int(cell_id) % N_LON
    lon = -180.0 + (col + 0.5) * (360.0 / N_LON)
    sin_lo = -1.0 + row * (2.0 / N_SINLAT)
    sin_hi = -1.0 + (row + 1) * (2.0 / N_SINLAT)
    lat = float(np.degrees(np.arcsin(np.clip(0.5 * (sin_lo + sin_hi), -1.0, 1.0))))
    return lat, lon


def great_circle_km(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    lat1 = np.radians(np.asarray(lat1, dtype=float)); lat2 = np.radians(np.asarray(lat2, dtype=float))
    lon1 = np.radians(np.asarray(lon1, dtype=float)); lon2 = np.radians(np.asarray(lon2, dtype=float))
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def build_species_cells(measured: pd.DataFrame) -> pd.DataFrame:
    work = measured.loc[measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS)].copy()
    raw = work[PALETTE].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).all(): raise RuntimeError("non-finite classifiable palette fractions")
    four = np.column_stack([
        raw[:, 0], raw[:, 1] + raw[:, 2] + raw[:, 8],
        raw[:, 3] + raw[:, 4] + raw[:, 5], raw[:, 6] + raw[:, 7],
    ])
    mass = four.sum(axis=1)
    if np.any(mass <= 0.0): raise RuntimeError("zero biological colour mass")
    four = four / mass[:, None]
    for i, name in enumerate(COLOUR_COLS): work[name] = four[:, i]
    agg = work.groupby(["species", "cell_id"], sort=True, observed=True).agg(
        classifiable_photos=("measurement_id", "size"), **{c: (c, "mean") for c in COLOUR_COLS}
    ).reset_index()
    agg.loc[:, COLOUR_COLS] = agg[COLOUR_COLS].to_numpy(dtype=float) / agg[COLOUR_COLS].sum(axis=1).to_numpy(dtype=float)[:, None]
    return agg


def colour_matrix(colour: np.ndarray) -> np.ndarray:
    n = len(colour); out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = jsd_bits(colour[i], colour[j]); out[i, j] = d; out[j, i] = d
    return out


def prepare_species(species: str, group: pd.DataFrame, caprow: pd.Series, *, permutations: int, seed: int) -> tuple[dict[str, Any], np.ndarray] | None:
    group = group.sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    n = len(group)
    if n < 3: return None
    colour = group[COLOUR_COLS].to_numpy(dtype=float)
    ii, jj = np.triu_indices(n, k=1)
    cmat = colour_matrix(colour); colour_dist = cmat[ii, jj]
    centers = np.array([cell_centroid(int(c)) for c in group["cell_id"].to_numpy()], dtype=float)
    geo = great_circle_km(centers[ii, 0], centers[ii, 1], centers[jj, 0], centers[jj, 1])
    rho = spearman(geo, colour_dist)
    if not np.isfinite(rho): return None
    fixed = fixed_rank_unit(geo)
    null = np.empty(permutations, dtype=float)
    for p in range(permutations):
        order = stable_species_permutation(n, seed=seed, permutation=p, species=species)
        prho = spearman_against_fixed_unit(fixed, cmat[order[ii], order[jj]])
        if not np.isfinite(prho): raise RuntimeError(f"{species} non-evaluable under H6 permutation {p}")
        null[p] = prho
    rec = {
        "species": species,
        "pre_frame_classifiable_photos": int(caprow["classifiable_photos"]),
        "pre_frame_occupied_h1_cells": int(caprow["occupied_h1_cells"]),
        "pre_frame_morph_levels": int(caprow["morph_levels"]),
        "spatial_cells": int(n),
        "classifiable_photos_in_spatial_cells": int(group["classifiable_photos"].sum()),
        "cell_pairs": int(len(ii)),
        "geographic_colour_spearman_rho": float(rho),
        "mean_pairwise_h1_cell_centroid_great_circle_km": float(np.mean(geo)),
        "mean_colour_jsd_bits": float(np.mean(colour_dist)),
    }
    return rec, null


def frame_mask(df: pd.DataFrame, min_photos: int, min_cells: int, min_morphs: int) -> np.ndarray:
    return (
        (df["pre_frame_classifiable_photos"].to_numpy(dtype=int) >= min_photos)
        & (df["pre_frame_occupied_h1_cells"].to_numpy(dtype=int) >= min_cells)
        & (df["pre_frame_morph_levels"].to_numpy(dtype=int) >= min_morphs)
        & (df["spatial_cells"].to_numpy(dtype=int) >= min_cells)
    )


def analyze_frame(frame: pd.DataFrame, null_by_species: dict[str, np.ndarray], *, permutations: int) -> tuple[dict[str, Any], np.ndarray]:
    obs = float(frame["geographic_colour_spearman_rho"].mean())
    matrix = np.vstack([null_by_species[s] for s in frame["species"].astype(str)])
    null = matrix.mean(axis=0)
    p_upper = float((1 + np.count_nonzero(null >= obs)) / (permutations + 1))
    return {
        "evaluable_species": int(len(frame)),
        "observed_mean_spearman_rho": obs,
        "observed_median_spearman_rho": float(frame["geographic_colour_spearman_rho"].median()),
        "positive_species_fraction": float((frame["geographic_colour_spearman_rho"] > 0).mean()),
        "p_upper": p_upper,
        "null_mean": float(np.mean(null)),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }, null


def main() -> int:
    contract = load_json(CONTRACT); measurement = load_json(MEASUREMENT_RESULT)
    capacity_result = load_json(CAPACITY_RESULT); h1 = load_json(H1_RESULT); h2 = load_json(H2_RESULT)
    h4a = load_json(H4A_RESULT); h5 = load_json(H5_RESULT)
    if contract.get("status") != "exploratory_frozen_after_h5_before_any_h6_spatial_colour_test": raise RuntimeError("H6 contract not frozen")
    if measurement.get("measurement_table_sha256") != sha256(MEASURED): raise RuntimeError("measurement SHA drifted")
    if int(measurement.get("classified_rows", -1)) != 10103: raise RuntimeError("classified denominator drifted")
    if capacity_result.get("status") != "complete_exploratory_capacity_audit_no_h4_model_fit": raise RuntimeError("capacity lineage drifted")
    if h1.get("decision") != "no_support_excess_recurrent_boundary_concentration": raise RuntimeError("H1 decision drifted")
    if h2.get("hierarchical_decision") != "diagnostic_only_h1_not_supported_no_climate_mechanism_claim": raise RuntimeError("H2 decision drifted")
    if h4a.get("decision") != "no_support_within_species_climate_colour_divergence_do_not_open_h4b": raise RuntimeError("H4a decision drifted")
    if h5.get("decision") != "no_support_scale_dependent_weakening": raise RuntimeError("H5 decision drifted")

    permutations = int(contract["null"]["permutations"]); seed = int(contract["null"]["seed"])
    if permutations != 999 or seed != 20260907: raise RuntimeError("H6 null settings drifted")
    measured = pd.read_csv(MEASURED)
    required = {"measurement_id", "species", "cell_id", "morph", *PALETTE}
    missing = sorted(required.difference(measured.columns))
    if missing: raise RuntimeError(f"measurement table lacks H6 fields: {missing}")
    if int(measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS).sum()) != 10103: raise RuntimeError("classifiable replay failed")

    capacity = pd.read_csv(CAPACITY); capacity["species"] = capacity["species"].astype(str)
    capidx = capacity.set_index("species", drop=False)
    species_cells = build_species_cells(measured)
    union = set(capacity.loc[(capacity["classifiable_photos"] >= 5) & (capacity["occupied_h1_cells"] >= 3) & (capacity["morph_levels"] >= 2), "species"].astype(str))
    candidate = species_cells.loc[species_cells["species"].astype(str).isin(union)].copy()

    records: list[dict[str, Any]] = []; null_by_species: dict[str, np.ndarray] = {}
    for species, group in candidate.groupby("species", sort=True, observed=True):
        name = str(species)
        prepared = prepare_species(name, group, capidx.loc[name], permutations=permutations, seed=seed)
        if prepared is None: continue
        rec, null = prepared; records.append(rec); null_by_species[name] = null
    all_species = pd.DataFrame(records).sort_values("species", kind="mergesort").reset_index(drop=True)

    specs = {
        "primary_n10_cells5_morphs2": (10, 5, 2),
        "n5_cells3_morphs2": (5, 3, 2),
        "n10_cells3_morphs2": (10, 3, 2),
        "n20_cells3_morphs2": (20, 3, 2),
    }
    frames = {name: all_species.loc[frame_mask(all_species, *vals)].copy().reset_index(drop=True) for name, vals in specs.items()}
    expected = {"primary_n10_cells5_morphs2": 74, "n5_cells3_morphs2": 188, "n10_cells3_morphs2": 80, "n20_cells3_morphs2": 30}
    for name, n in expected.items():
        if len(frames[name]) != n: raise RuntimeError(f"H6 capacity-frame drift for {name}: {len(frames[name])} != {n}")

    primary, primary_null = analyze_frame(frames["primary_n10_cells5_morphs2"], null_by_species, permutations=permutations)
    primary.update({
        "name": "primary_n10_cells5_morphs2", "primary": True,
        "minimum_evaluable_species_required": int(contract["primary_statistic"]["minimum_evaluable_species"]),
        "permutations": permutations, "seed": seed, "alpha": float(contract["null"]["alpha"]),
    })
    support = primary["evaluable_species"] >= primary["minimum_evaluable_species_required"] and primary["observed_mean_spearman_rho"] > 0 and primary["p_upper"] < primary["alpha"]
    primary["support"] = bool(support)

    sens_rows: list[dict[str, Any]] = []
    for name in ["n5_cells3_morphs2", "n10_cells3_morphs2", "n20_cells3_morphs2"]:
        row, _ = analyze_frame(frames[name], null_by_species, permutations=permutations)
        row.update({"name": name, "primary": False, "nominal_positive": bool(row["observed_mean_spearman_rho"] > 0 and row["p_upper"] < 0.05), "can_rescue_primary": False})
        sens_rows.append(row)

    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    frames["primary_n10_cells5_morphs2"].to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    pd.DataFrame({"permutation": np.arange(permutations, dtype=int), "mean_species_spearman_rho": primary_null}).to_csv(OUT_NULL, index=False, lineterminator="\n")
    pd.DataFrame(sens_rows).to_csv(OUT_SENS, index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"], "status": "complete_exploratory_h6_evaluable",
        "introduced_after_h1_h2_h4a_h5_outcomes": True, "claim_role": contract["claim_role"],
        "primary": primary, "sensitivities": sens_rows,
        "decision": "support_species_specific_spatial_structuring_exploratory_only" if support else "no_support_general_species_specific_spatial_structuring",
        "shared_global_boundary_claim_opened": False, "climate_mechanism_claim_opened": False,
        "frozen_primary_decisions_unchanged": {"h1": h1["decision"], "h2": h2["hierarchical_decision"], "h4a": h4a["decision"], "h5": h5["decision"]},
        "lineage": {"contract_sha256": sha256(CONTRACT), "measurement_table_sha256": sha256(MEASURED), "capacity_table_sha256": sha256(CAPACITY), "permutations": permutations, "seed": seed},
        "files": {"primary_species": str(OUT_SPECIES.relative_to(ROOT)), "primary_null": str(OUT_NULL.relative_to(ROOT)), "sensitivities": str(OUT_SENS.relative_to(ROOT))},
        "claim_ceiling": contract["claim_ceiling"],
    }
    OUT_RESULT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
