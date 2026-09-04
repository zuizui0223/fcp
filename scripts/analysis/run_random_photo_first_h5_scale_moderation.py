#!/usr/bin/env python3
"""Run post-H4a exploratory H5: spatial-scale moderation of within-species climate-colour coupling.

H5 does not rescue the non-significant H1/H2/H4a primary tests. It asks a new,
explicitly post-H4a question motivated by the H4a cell-threshold sensitivity:
are species-level climate-colour associations stronger in spatially localized
species and weaker as occupied geographic span increases?

The primary response is exactly the H4a species-level Spearman correlation between
pairwise colour JSD and pairwise RMS standardized macroclimate distance. The primary
moderator is log1p(mean pairwise great-circle distance among occupied H1-cell
centroids). The null preserves each species' geography and colour-distance matrix
while independently permuting colour cell labels within species.
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
CONTRACT = ROOT / "docs/supporting/random_photo_first_h5_scale_moderation_contract_v1.json"
MEASURED = ROOT / "data/derived/random_photo_first_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_measurement_result_v1.json"
CAPACITY = ROOT / "data/derived/random_photo_first_h4_within_species_capacity_v1.csv"
CLIMATE = ROOT / "data/derived/random_photo_first_h2_climate_cells_250km_v1.csv"
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"
H2_RESULT = ROOT / "docs/supporting/random_photo_first_h2_result_v1.json"
H4A_RESULT = ROOT / "docs/supporting/random_photo_first_h4a_result_v1.json"

OUT_RESULT = ROOT / "docs/supporting/random_photo_first_h5_scale_moderation_result_v1.json"
OUT_SPECIES = ROOT / "data/derived/random_photo_first_h5_species_primary_v1.csv"
OUT_NULL = ROOT / "data/derived/random_photo_first_h5_null_primary_v1.csv"
OUT_SENS = ROOT / "data/derived/random_photo_first_h5_sensitivities_v1.csv"

BIOLOGICAL_MORPHS = {"white", "yellow_orange", "red_pink", "blue_purple"}
PALETTE = [
    "flower_fraction_white",
    "flower_fraction_yellow",
    "flower_fraction_orange",
    "flower_fraction_red",
    "flower_fraction_pink",
    "flower_fraction_magenta",
    "flower_fraction_purple",
    "flower_fraction_blue",
    "flower_fraction_bronze",
]
CLIMATE_Z = ["z_bio1", "z_bio4", "z_bio12", "z_bio15"]
COLOUR_COLS = ["soft_white", "soft_yellow_orange", "soft_red_pink", "soft_blue_purple"]
N_LON = 18
N_SINLAT = 9


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def jsd_bits(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    out = 0.0
    keep = p > 0
    if keep.any():
        out += 0.5 * float(np.sum(p[keep] * np.log2(p[keep] / m[keep])))
    keep = q > 0
    if keep.any():
        out += 0.5 * float(np.sum(q[keep] * np.log2(q[keep] / m[keep])))
    return float(np.clip(out, 0.0, 1.0))


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    dx = x - float(np.mean(x))
    dy = y - float(np.mean(y))
    nx = float(np.sqrt(np.sum(dx * dx)))
    ny = float(np.sqrt(np.sum(dy * dy)))
    if nx <= 0.0 or ny <= 0.0:
        return float("nan")
    return float(np.clip(np.sum(dx * dy) / (nx * ny), -1.0, 1.0))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rankdata(np.asarray(x, dtype=float), method="average"), rankdata(np.asarray(y, dtype=float), method="average"))


def fixed_rank_unit(x: np.ndarray) -> np.ndarray:
    r = rankdata(np.asarray(x, dtype=float), method="average")
    r = r - float(np.mean(r))
    norm = float(np.sqrt(np.sum(r * r)))
    if norm <= 0.0:
        raise RuntimeError("fixed rank vector is constant")
    return r / norm


def spearman_against_fixed_unit(fixed_unit: np.ndarray, y: np.ndarray) -> float:
    r = rankdata(np.asarray(y, dtype=float), method="average")
    r = r - float(np.mean(r))
    norm = float(np.sqrt(np.sum(r * r)))
    if norm <= 0.0:
        return float("nan")
    return float(np.clip(np.sum(fixed_unit * (r / norm)), -1.0, 1.0))


def stable_species_permutation(n: int, *, seed: int, permutation: int, species: str) -> np.ndarray:
    token = f"{seed}|{permutation}|{species}".encode("utf-8")
    child = int.from_bytes(hashlib.sha256(token).digest()[:8], "big", signed=False)
    return np.random.default_rng(child).permutation(n)


def cell_centroid(cell_id: int) -> tuple[float, float]:
    cell = int(cell_id)
    if cell < 0 or cell >= N_LON * N_SINLAT:
        raise ValueError(f"invalid H1 cell_id: {cell}")
    row = cell // N_LON
    col = cell % N_LON
    lon = -180.0 + (col + 0.5) * (360.0 / N_LON)
    sin_lo = -1.0 + row * (2.0 / N_SINLAT)
    sin_hi = -1.0 + (row + 1) * (2.0 / N_SINLAT)
    lat = float(np.degrees(np.arcsin(np.clip(0.5 * (sin_lo + sin_hi), -1.0, 1.0))))
    return lat, lon


def great_circle_km(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    lat1r = np.radians(np.asarray(lat1, dtype=float))
    lon1r = np.radians(np.asarray(lon1, dtype=float))
    lat2r = np.radians(np.asarray(lat2, dtype=float))
    lon2r = np.radians(np.asarray(lon2, dtype=float))
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def build_species_cells(measured: pd.DataFrame, climate: pd.DataFrame) -> pd.DataFrame:
    work = measured.loc[measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS)].copy()
    raw = work[PALETTE].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).all():
        raise RuntimeError("classifiable palette fractions contain non-finite values")
    four = np.column_stack(
        [
            raw[:, 0],
            raw[:, 1] + raw[:, 2] + raw[:, 8],
            raw[:, 3] + raw[:, 4] + raw[:, 5],
            raw[:, 6] + raw[:, 7],
        ]
    )
    mass = four.sum(axis=1)
    if np.any(mass <= 0.0):
        raise RuntimeError("classifiable row has zero biological colour mass")
    four = four / mass[:, None]
    for idx, name in enumerate(COLOUR_COLS):
        work[name] = four[:, idx]

    agg = (
        work.groupby(["species", "cell_id"], sort=True, observed=True)
        .agg(classifiable_photos=("measurement_id", "size"), **{c: (c, "mean") for c in COLOUR_COLS})
        .reset_index()
    )
    row_mass = agg[COLOUR_COLS].sum(axis=1).to_numpy(dtype=float)
    agg.loc[:, COLOUR_COLS] = agg[COLOUR_COLS].to_numpy(dtype=float) / row_mass[:, None]

    needed = {"cell_id", "complete_macroclimate", *CLIMATE_Z}
    if not needed.issubset(climate.columns):
        raise RuntimeError("frozen H2 climate cell table lacks required H5 fields")
    joined = agg.merge(climate[["cell_id", "complete_macroclimate", *CLIMATE_Z]], on="cell_id", how="left", validate="many_to_one")
    joined = joined.loc[joined["complete_macroclimate"].fillna(False).astype(bool)].dropna(subset=CLIMATE_Z).copy()
    if not np.isfinite(joined[CLIMATE_Z].to_numpy(dtype=float)).all():
        raise RuntimeError("climate-complete H5 rows contain non-finite climate")
    return joined


def colour_distance_matrix(colour: np.ndarray) -> np.ndarray:
    n = len(colour)
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = jsd_bits(colour[i], colour[j])
            out[i, j] = d
            out[j, i] = d
    return out


def prepare_species_record(species: str, group: pd.DataFrame, capacity_row: pd.Series, *, permutations: int, seed: int) -> tuple[dict[str, Any], np.ndarray] | None:
    group = group.sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    n = len(group)
    if n < 3:
        return None
    colour = group[COLOUR_COLS].to_numpy(dtype=float)
    climate = group[CLIMATE_Z].to_numpy(dtype=float)
    ii, jj = np.triu_indices(n, k=1)
    climate_delta = climate[ii] - climate[jj]
    climate_distance = np.sqrt(np.mean(climate_delta * climate_delta, axis=1))
    cmat = colour_distance_matrix(colour)
    colour_distance = cmat[ii, jj]
    observed_rho = spearman(climate_distance, colour_distance)
    if not np.isfinite(observed_rho):
        return None

    centers = np.array([cell_centroid(int(c)) for c in group["cell_id"].to_numpy()], dtype=float)
    geo = great_circle_km(centers[ii, 0], centers[ii, 1], centers[jj, 0], centers[jj, 1])
    climate_rank_unit = fixed_rank_unit(climate_distance)

    null_rhos = np.empty(permutations, dtype=float)
    for perm_index in range(permutations):
        order = stable_species_permutation(n, seed=seed, permutation=perm_index, species=species)
        permuted_colour_distance = cmat[order[ii], order[jj]]
        rho = spearman_against_fixed_unit(climate_rank_unit, permuted_colour_distance)
        if not np.isfinite(rho):
            raise RuntimeError(f"{species} became non-evaluable under permutation {perm_index}")
        null_rhos[perm_index] = rho

    record = {
        "species": species,
        "pre_frame_classifiable_photos": int(capacity_row["classifiable_photos"]),
        "pre_frame_occupied_h1_cells": int(capacity_row["occupied_h1_cells"]),
        "pre_frame_morph_levels": int(capacity_row["morph_levels"]),
        "climate_complete_cells": int(n),
        "climate_complete_classifiable_photos": int(group["classifiable_photos"].sum()),
        "cell_pairs": int(len(ii)),
        "climate_colour_spearman_rho": float(observed_rho),
        "mean_pairwise_h1_cell_centroid_great_circle_km": float(np.mean(geo)),
        "max_pairwise_h1_cell_centroid_great_circle_km": float(np.max(geo)),
        "mean_climate_rms_distance": float(np.mean(climate_distance)),
        "mean_colour_jsd_bits": float(np.mean(colour_distance)),
    }
    return record, null_rhos


def frame_mask(table: pd.DataFrame, *, min_photos: int, min_cells: int, min_morphs: int) -> np.ndarray:
    return (
        (table["pre_frame_classifiable_photos"].to_numpy(dtype=int) >= min_photos)
        & (table["pre_frame_occupied_h1_cells"].to_numpy(dtype=int) >= min_cells)
        & (table["pre_frame_morph_levels"].to_numpy(dtype=int) >= min_morphs)
        & (table["climate_complete_cells"].to_numpy(dtype=int) >= min_cells)
    )


def moderation_test(
    table: pd.DataFrame,
    null_by_species: dict[str, np.ndarray],
    *,
    moderator: str,
    transform_log1p: bool,
    permutations: int,
) -> tuple[dict[str, Any], np.ndarray]:
    species = table["species"].astype(str).tolist()
    x = table[moderator].to_numpy(dtype=float)
    if transform_log1p:
        x = np.log1p(x)
    y = table["climate_colour_spearman_rho"].to_numpy(dtype=float)
    observed = spearman(x, y)
    if not np.isfinite(observed):
        raise RuntimeError(f"H5 moderator is not evaluable: {moderator}")
    fixed_x = fixed_rank_unit(x)
    null = np.empty(permutations, dtype=float)
    matrix = np.vstack([null_by_species[s] for s in species])
    for perm_index in range(permutations):
        value = spearman_against_fixed_unit(fixed_x, matrix[:, perm_index])
        if not np.isfinite(value):
            raise RuntimeError(f"H5 moderation null became non-evaluable: {moderator}, {perm_index}")
        null[perm_index] = value
    p_lower = float((1 + np.count_nonzero(null <= observed)) / (permutations + 1))
    return {
        "moderator": moderator,
        "transform": "log1p" if transform_log1p else "none",
        "evaluable_species": int(len(table)),
        "observed_scale_moderation_spearman_rho": float(observed),
        "p_lower": p_lower,
        "null_mean": float(np.mean(null)),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }, null


def h4a_sensitivity(h4a: dict[str, Any], name: str) -> dict[str, Any]:
    for row in h4a.get("sensitivities", []):
        if row.get("name") == name:
            return row
    raise RuntimeError(f"completed H4a lacks required sensitivity: {name}")


def main() -> int:
    contract = load_json(CONTRACT)
    measurement = load_json(MEASUREMENT_RESULT)
    h1 = load_json(H1_RESULT)
    h2 = load_json(H2_RESULT)
    h4a = load_json(H4A_RESULT)

    if contract.get("status") != "exploratory_frozen_after_h4a_before_any_h5_scale_moderation_test":
        raise RuntimeError("H5 contract is not frozen before H5 testing")
    if contract.get("introduced_after_h4a_outcomes") is not True:
        raise RuntimeError("H5 must remain explicitly post-H4a exploratory")
    if h1.get("decision") != "no_support_excess_recurrent_boundary_concentration":
        raise RuntimeError("H1 decision drifted")
    if h2.get("hierarchical_decision") != "diagnostic_only_h1_not_supported_no_climate_mechanism_claim":
        raise RuntimeError("H2 decision drifted")
    if h4a.get("decision") != "no_support_within_species_climate_colour_divergence_do_not_open_h4b":
        raise RuntimeError("H4a decision drifted")
    if float(h4a["primary"]["p_upper"]) != 0.157 or h4a.get("h4b_directional_decomposition_opened") is not False:
        raise RuntimeError("H4a primary lineage drifted")
    if measurement.get("measurement_table_sha256") != sha256(MEASURED):
        raise RuntimeError("measurement SHA256 drifted")
    if int(measurement.get("classified_rows", -1)) != 10103:
        raise RuntimeError("classified measurement denominator drifted")

    permutations = int(contract["null"]["permutations"])
    seed = int(contract["null"]["seed"])
    if permutations != 999 or seed != 20260906:
        raise RuntimeError("H5 frozen null settings drifted")

    measured = pd.read_csv(MEASURED)
    needed_measure = {"measurement_id", "species", "cell_id", "morph", *PALETTE}
    missing = sorted(needed_measure.difference(measured.columns))
    if missing:
        raise RuntimeError(f"measurement table lacks H5 columns: {missing}")
    classifiable_n = int(measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS).sum())
    if classifiable_n != 10103:
        raise RuntimeError("H5 classifiable rows do not replay measurement manifest")

    capacity = pd.read_csv(CAPACITY)
    needed_capacity = {"species", "classifiable_photos", "occupied_h1_cells", "morph_levels"}
    if not needed_capacity.issubset(capacity.columns):
        raise RuntimeError("capacity table lacks H5 fields")
    capacity = capacity.copy()
    capacity["species"] = capacity["species"].astype(str)
    capacity_index = capacity.set_index("species", drop=False)

    climate = pd.read_csv(CLIMATE)
    species_cells = build_species_cells(measured, climate)

    union_species = set(
        capacity.loc[
            (capacity["classifiable_photos"] >= 5)
            & (capacity["occupied_h1_cells"] >= 3)
            & (capacity["morph_levels"] >= 2),
            "species",
        ].astype(str)
    )
    candidate = species_cells.loc[species_cells["species"].astype(str).isin(union_species)].copy()

    records: list[dict[str, Any]] = []
    null_by_species: dict[str, np.ndarray] = {}
    for species, group in candidate.groupby("species", sort=True, observed=True):
        species_name = str(species)
        prepared = prepare_species_record(
            species_name,
            group,
            capacity_index.loc[species_name],
            permutations=permutations,
            seed=seed,
        )
        if prepared is None:
            continue
        record, null_rhos = prepared
        records.append(record)
        null_by_species[species_name] = null_rhos

    all_species = pd.DataFrame(records).sort_values("species", kind="mergesort").reset_index(drop=True)
    if len(all_species) < 30:
        raise RuntimeError("H5 has fewer than 30 climate-colour-evaluable species in union frame")

    frame_specs = {
        "primary_n10_cells3_morphs2": (10, 3, 2),
        "broader_n5_cells3_morphs2": (5, 3, 2),
        "strict_n20_cells3_morphs2": (20, 3, 2),
    }
    frames: dict[str, pd.DataFrame] = {}
    for name, (min_photos, min_cells, min_morphs) in frame_specs.items():
        frames[name] = all_species.loc[
            frame_mask(all_species, min_photos=min_photos, min_cells=min_cells, min_morphs=min_morphs)
        ].copy().reset_index(drop=True)

    # Replay the completed H4a sensitivity frames before testing H5 moderation.
    replay_targets = {
        "primary_n10_cells3_morphs2": "n10_cells3_morphs2",
        "broader_n5_cells3_morphs2": "n5_cells3_morphs2",
        "strict_n20_cells3_morphs2": "n20_cells3_morphs2",
    }
    replay: dict[str, Any] = {}
    for h5_name, h4_name in replay_targets.items():
        frame = frames[h5_name]
        h4 = h4a_sensitivity(h4a, h4_name)
        observed_mean = float(frame["climate_colour_spearman_rho"].mean())
        if int(h4["evaluable_species"]) != len(frame):
            raise RuntimeError(f"H5 frame does not replay H4a evaluable species for {h4_name}")
        if not np.isclose(observed_mean, float(h4["observed_mean_spearman_rho"]), atol=1e-12, rtol=0.0):
            raise RuntimeError(f"H5 species rhos do not replay completed H4a mean for {h4_name}")
        replay[h4_name] = {
            "evaluable_species": int(len(frame)),
            "observed_mean_spearman_rho": observed_mean,
            "completed_h4a_mean": float(h4["observed_mean_spearman_rho"]),
            "exact_replay_within_1e_12": True,
        }

    primary_frame = frames["primary_n10_cells3_morphs2"]
    primary, primary_null = moderation_test(
        primary_frame,
        null_by_species,
        moderator="mean_pairwise_h1_cell_centroid_great_circle_km",
        transform_log1p=True,
        permutations=permutations,
    )
    minimum_species = int(contract["minimum_evaluable_species"])
    alpha = float(contract["null"]["alpha"])
    support = (
        primary["evaluable_species"] >= minimum_species
        and primary["observed_scale_moderation_spearman_rho"] < 0.0
        and primary["p_lower"] < alpha
    )
    primary.update(
        {
            "name": "primary_n10_cells3_morphs2_mean_pairwise_geographic_span",
            "primary": True,
            "minimum_evaluable_species_required": minimum_species,
            "permutations": permutations,
            "seed": seed,
            "alpha": alpha,
            "support": bool(support),
        }
    )

    sensitivity_rows: list[dict[str, Any]] = []
    for name, moderator, log_transform, frame_name in [
        ("max_pairwise_geographic_span", "max_pairwise_h1_cell_centroid_great_circle_km", True, "primary_n10_cells3_morphs2"),
        ("occupied_h1_cell_count", "climate_complete_cells", False, "primary_n10_cells3_morphs2"),
        ("broader_n5_cells3_morphs2", "mean_pairwise_h1_cell_centroid_great_circle_km", True, "broader_n5_cells3_morphs2"),
        ("strict_n20_cells3_morphs2", "mean_pairwise_h1_cell_centroid_great_circle_km", True, "strict_n20_cells3_morphs2"),
    ]:
        result, _ = moderation_test(
            frames[frame_name],
            null_by_species,
            moderator=moderator,
            transform_log1p=log_transform,
            permutations=permutations,
        )
        result.update({"name": name, "primary": False, "can_rescue_primary": False})
        sensitivity_rows.append(result)

    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    primary_frame.to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    pd.DataFrame({"permutation": np.arange(permutations, dtype=int), "scale_moderation_spearman_rho": primary_null}).to_csv(
        OUT_NULL, index=False, lineterminator="\n"
    )
    pd.DataFrame(sensitivity_rows).to_csv(OUT_SENS, index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"],
        "status": "complete_exploratory_h5_evaluable",
        "introduced_after_h4a_outcomes": True,
        "claim_role": contract["claim_role"],
        "h4a_replay": replay,
        "primary": primary,
        "sensitivities": sensitivity_rows,
        "decision": "support_scale_dependent_weakening_exploratory_only" if support else "no_support_scale_dependent_weakening",
        "h4b_directional_decomposition_opened": False,
        "frozen_primary_decisions_unchanged": {
            "h1": h1["decision"],
            "h1_p_upper": float(h1["primary"]["p_upper"]),
            "h2": h2["hierarchical_decision"],
            "h4a": h4a["decision"],
            "h4a_p_upper": float(h4a["primary"]["p_upper"]),
        },
        "lineage": {
            "contract_sha256": sha256(CONTRACT),
            "measurement_table_sha256": sha256(MEASURED),
            "capacity_table_sha256": sha256(CAPACITY),
            "climate_cells_sha256": sha256(CLIMATE),
            "h4a_result_sha256": sha256(H4A_RESULT),
            "permutations": permutations,
            "seed": seed,
        },
        "files": {
            "primary_species": str(OUT_SPECIES.relative_to(ROOT)),
            "primary_null": str(OUT_NULL.relative_to(ROOT)),
            "sensitivities": str(OUT_SENS.relative_to(ROOT)),
        },
        "claim_ceiling": contract["claim_ceiling"],
    }
    OUT_RESULT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
