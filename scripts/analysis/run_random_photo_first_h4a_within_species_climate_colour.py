#!/usr/bin/env python3
"""Run the frozen post-H1/H2 exploratory H4a.

H4a asks whether, within species, greater macroclimate separation between occupied
H1 cells is associated with greater flower-colour-composition separation. The
analysis contract was committed before this script opens any H4a climate-colour
join or statistic.

The four coarse morph labels are never ordered. Colour is represented by the
continuous nine-palette fractions retained by the fresh measurement, prospectively
collapsed to four soft groups and summarized within species x H1 cell. Pairwise
colour distance is Jensen-Shannon divergence; climate distance is RMS separation
across the already-frozen H2 z-BIO1/z-BIO4/z-BIO12/z-BIO15 cell values.

Pairwise distances are dependent, so no pairwise iid p-values are used. The global
statistic is an equal-species mean of species-level Spearman correlations, and the
null independently permutes colour cell labels within each species for 999 matched
global permutations.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h4a_within_species_climate_colour_contract_v1.json"
MEASURED = ROOT / "data/derived/random_photo_first_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_measurement_result_v1.json"
CAPACITY = ROOT / "data/derived/random_photo_first_h4_within_species_capacity_v1.csv"
CAPACITY_RESULT = ROOT / "docs/supporting/random_photo_first_h4_within_species_capacity_v1.json"
CLIMATE = ROOT / "data/derived/random_photo_first_h2_climate_cells_250km_v1.csv"
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"
H2_RESULT = ROOT / "docs/supporting/random_photo_first_h2_result_v1.json"

OUT_RESULT = ROOT / "docs/supporting/random_photo_first_h4a_result_v1.json"
OUT_SPECIES = ROOT / "data/derived/random_photo_first_h4a_species_primary_v1.csv"
OUT_NULL = ROOT / "data/derived/random_photo_first_h4a_null_primary_v1.csv"
OUT_SENS = ROOT / "data/derived/random_photo_first_h4a_sensitivities_v1.csv"

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


@dataclass(frozen=True)
class FrameSpec:
    name: str
    min_photos: int
    min_cells: int
    min_morphs: int
    primary: bool = False


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
    terms = 0.0
    keep = p > 0
    if keep.any():
        terms += 0.5 * float(np.sum(p[keep] * np.log2(p[keep] / m[keep])))
    keep = q > 0
    if keep.any():
        terms += 0.5 * float(np.sum(q[keep] * np.log2(q[keep] / m[keep])))
    return float(np.clip(terms, 0.0, 1.0))


def average_ranks(x: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(x, dtype=float)).rank(method="average").to_numpy(dtype=float)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    dx = x - float(np.mean(x))
    dy = y - float(np.mean(y))
    sx = float(np.sqrt(np.sum(dx * dx)))
    sy = float(np.sqrt(np.sum(dy * dy)))
    if sx <= 0.0 or sy <= 0.0:
        return float("nan")
    return float(np.clip(np.sum(dx * dy) / (sx * sy), -1.0, 1.0))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(average_ranks(x), average_ranks(y))


def stable_species_permutation(n: int, *, seed: int, permutation: int, species: str) -> np.ndarray:
    token = f"{seed}|{permutation}|{species}".encode("utf-8")
    child = int.from_bytes(hashlib.sha256(token).digest()[:8], "big", signed=False)
    return np.random.default_rng(child).permutation(n)


def pair_vectors(colour: np.ndarray, climate: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = len(colour)
    ii, jj = np.triu_indices(n, k=1)
    climate_delta = climate[ii] - climate[jj]
    climate_distance = np.sqrt(np.mean(climate_delta * climate_delta, axis=1))
    colour_distance = np.fromiter(
        (jsd_bits(colour[i], colour[j]) for i, j in zip(ii, jj)),
        dtype=float,
        count=len(ii),
    )
    return ii, jj, climate_distance, colour_distance


def permuted_colour_distances(colour: np.ndarray, ii: np.ndarray, jj: np.ndarray, order: np.ndarray) -> np.ndarray:
    perm = colour[order]
    return np.fromiter(
        (jsd_bits(perm[i], perm[j]) for i, j in zip(ii, jj)),
        dtype=float,
        count=len(ii),
    )


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
    row_sum = four.sum(axis=1)
    if np.any(row_sum <= 0.0):
        raise RuntimeError("classifiable row has zero four-group biological colour mass")
    four = four / row_sum[:, None]
    for idx, name in enumerate(["soft_white", "soft_yellow_orange", "soft_red_pink", "soft_blue_purple"]):
        work[name] = four[:, idx]

    colour_cols = ["soft_white", "soft_yellow_orange", "soft_red_pink", "soft_blue_purple"]
    agg = (
        work.groupby(["species", "cell_id"], sort=True, observed=True)
        .agg(
            classifiable_photos=("measurement_id", "size"),
            **{c: (c, "mean") for c in colour_cols},
        )
        .reset_index()
    )
    sums = agg[colour_cols].sum(axis=1).to_numpy(dtype=float)
    if np.any(sums <= 0.0):
        raise RuntimeError("species-cell mean colour composition has zero mass")
    agg.loc[:, colour_cols] = agg[colour_cols].to_numpy(dtype=float) / sums[:, None]

    climate_cols = ["cell_id", "complete_macroclimate", *CLIMATE_Z]
    if not set(climate_cols).issubset(climate.columns):
        raise RuntimeError("frozen climate cell table lacks required H4a columns")
    joined = agg.merge(climate[climate_cols], on="cell_id", how="left", validate="many_to_one", sort=False)
    joined = joined.loc[joined["complete_macroclimate"].fillna(False).astype(bool)].copy()
    joined = joined.dropna(subset=CLIMATE_Z)
    if not np.isfinite(joined[CLIMATE_Z].to_numpy(dtype=float)).all():
        raise RuntimeError("climate-complete H4a cells contain non-finite z climate")
    return joined


def frame_pre_species(capacity: pd.DataFrame, spec: FrameSpec) -> set[str]:
    mask = (
        (capacity["classifiable_photos"] >= spec.min_photos)
        & (capacity["occupied_h1_cells"] >= spec.min_cells)
        & (capacity["morph_levels"] >= spec.min_morphs)
    )
    return set(capacity.loc[mask, "species"].astype(str))


def analyze_frame(
    species_cells: pd.DataFrame,
    capacity: pd.DataFrame,
    spec: FrameSpec,
    *,
    permutations: int,
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame, np.ndarray]:
    pre = frame_pre_species(capacity, spec)
    candidate = species_cells.loc[species_cells["species"].astype(str).isin(pre)].copy()

    colour_cols = ["soft_white", "soft_yellow_orange", "soft_red_pink", "soft_blue_purple"]
    records: list[dict[str, Any]] = []
    prepared: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []

    for species, group in candidate.groupby("species", sort=True, observed=True):
        group = group.sort_values("cell_id", kind="mergesort").reset_index(drop=True)
        if len(group) < spec.min_cells:
            continue
        colour = group[colour_cols].to_numpy(dtype=float)
        climate = group[CLIMATE_Z].to_numpy(dtype=float)
        ii, jj, climate_dist, colour_dist = pair_vectors(colour, climate)
        rho = spearman(climate_dist, colour_dist)
        if not np.isfinite(rho):
            continue
        species_name = str(species)
        prepared.append((species_name, colour, ii, jj, climate_dist, colour_dist))
        records.append(
            {
                "species": species_name,
                "pre_frame_classifiable_photos": int(capacity.loc[capacity["species"].astype(str) == species_name, "classifiable_photos"].iloc[0]),
                "climate_complete_cells": int(len(group)),
                "climate_complete_classifiable_photos": int(group["classifiable_photos"].sum()),
                "cell_pairs": int(len(ii)),
                "spearman_rho": float(rho),
                "mean_climate_rms_distance": float(np.mean(climate_dist)),
                "mean_colour_jsd_bits": float(np.mean(colour_dist)),
            }
        )

    species_table = pd.DataFrame(records)
    observed_n = len(species_table)
    if observed_n == 0:
        return (
            {
                "name": spec.name,
                "status": "not_evaluable_no_species",
                "pre_frame_species": int(len(pre)),
                "evaluable_species": 0,
            },
            species_table,
            np.array([], dtype=float),
        )

    observed_mean = float(species_table["spearman_rho"].mean())
    observed_median = float(species_table["spearman_rho"].median())
    positive_fraction = float((species_table["spearman_rho"] > 0).mean())

    null = np.empty(permutations, dtype=float)
    for perm_index in range(permutations):
        rhos: list[float] = []
        for species_name, colour, ii, jj, climate_dist, _ in prepared:
            order = stable_species_permutation(
                len(colour), seed=seed, permutation=perm_index, species=species_name
            )
            colour_dist = permuted_colour_distances(colour, ii, jj, order)
            rho = spearman(climate_dist, colour_dist)
            if not np.isfinite(rho):
                raise RuntimeError("species became non-evaluable under a label permutation")
            rhos.append(float(rho))
        null[perm_index] = float(np.mean(rhos))

    p_upper = float((1 + np.count_nonzero(null >= observed_mean)) / (permutations + 1))
    result = {
        "name": spec.name,
        "status": "complete_evaluable",
        "primary": bool(spec.primary),
        "threshold": {
            "minimum_classifiable_photos_before_climate_join": spec.min_photos,
            "minimum_occupied_h1_cells_before_climate_join": spec.min_cells,
            "minimum_coarse_morph_levels_before_climate_join": spec.min_morphs,
        },
        "pre_frame_species": int(len(pre)),
        "evaluable_species": int(observed_n),
        "evaluable_species_cell_rows": int(sum(len(species_cells.loc[(species_cells["species"].astype(str) == s)]) for s in species_table["species"].astype(str))),
        "observed_mean_spearman_rho": observed_mean,
        "observed_median_spearman_rho": observed_median,
        "positive_species_fraction": positive_fraction,
        "permutations": int(permutations),
        "seed": int(seed),
        "p_upper": p_upper,
        "null_mean": float(np.mean(null)),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }
    return result, species_table, null


def main() -> int:
    contract = load_json(CONTRACT)
    measurement = load_json(MEASUREMENT_RESULT)
    capacity_result = load_json(CAPACITY_RESULT)
    h1 = load_json(H1_RESULT)
    h2 = load_json(H2_RESULT)

    if contract.get("status") != "exploratory_frozen_after_h1_h2_outcomes_before_any_h4a_climate_colour_join_or_test":
        raise RuntimeError("H4a contract is not frozen at the required pre-join state")
    if measurement.get("measurement_table_sha256") != sha256(MEASURED):
        raise RuntimeError("measurement SHA256 drifted")
    if measurement.get("classified_rows") != 10103:
        raise RuntimeError("classifiable measurement denominator drifted")
    if capacity_result.get("status") != "complete_exploratory_capacity_audit_no_h4_model_fit":
        raise RuntimeError("pre-H4a capacity audit is missing or invalid")
    if h1.get("status") != "complete_h1_evaluable" or h1.get("decision") != "no_support_excess_recurrent_boundary_concentration":
        raise RuntimeError("completed H1 lineage/decision drifted")
    if h2.get("status") != "complete_h2_evaluable" or h2.get("hierarchical_decision") != "diagnostic_only_h1_not_supported_no_climate_mechanism_claim":
        raise RuntimeError("completed H2 lineage/decision drifted")
    if h2.get("files", {}).get("climate_cells") != str(CLIMATE.relative_to(ROOT)):
        raise RuntimeError("H4a climate table is not the completed frozen H2 primary cell frame")

    measured = pd.read_csv(MEASURED)
    required_measure = {"measurement_id", "species", "cell_id", "morph", *PALETTE}
    missing = sorted(required_measure.difference(measured.columns))
    if missing:
        raise RuntimeError(f"measurement table lacks H4a columns: {missing}")
    classifiable_n = int(measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS).sum())
    if classifiable_n != 10103:
        raise RuntimeError("H4a classifiable row count does not replay measurement manifest")

    capacity = pd.read_csv(CAPACITY)
    required_capacity = {"species", "classifiable_photos", "occupied_h1_cells", "morph_levels"}
    if not required_capacity.issubset(capacity.columns):
        raise RuntimeError("capacity table lacks required species eligibility columns")

    climate = pd.read_csv(CLIMATE)
    species_cells = build_species_cells(measured, climate)

    primary_cfg = contract["primary_species_frame"]
    specs = [
        FrameSpec(
            name="primary_n10_cells5_morphs2",
            min_photos=int(primary_cfg["minimum_classifiable_photos_before_climate_join"]),
            min_cells=int(primary_cfg["minimum_occupied_h1_cells_before_climate_join"]),
            min_morphs=int(primary_cfg["minimum_coarse_morph_levels_before_climate_join"]),
            primary=True,
        )
    ]
    specs.extend(
        FrameSpec(
            name=str(s["name"]),
            min_photos=int(s["minimum_classifiable_photos"]),
            min_cells=int(s["minimum_occupied_h1_cells"]),
            min_morphs=int(s["minimum_coarse_morph_levels"]),
            primary=False,
        )
        for s in contract["predeclared_sensitivities"]
    )

    permutations = int(contract["null"]["permutations"])
    seed = int(contract["null"]["seed"])
    if permutations != 999 or seed != 20260905:
        raise RuntimeError("H4a frozen permutation identity drifted")

    frame_results: list[dict[str, Any]] = []
    primary_species = pd.DataFrame()
    primary_null = np.array([], dtype=float)
    for spec in specs:
        result, species_table, null = analyze_frame(
            species_cells,
            capacity,
            spec,
            permutations=permutations,
            seed=seed,
        )
        frame_results.append(result)
        if spec.primary:
            primary_species = species_table
            primary_null = null

    primary = frame_results[0]
    min_species = int(primary_cfg["minimum_evaluable_species_to_run"])
    alpha = float(contract["decision"]["alpha"])
    evaluable = primary.get("status") == "complete_evaluable" and int(primary.get("evaluable_species", 0)) >= min_species
    support = bool(
        evaluable
        and float(primary["observed_mean_spearman_rho"]) > 0.0
        and float(primary["p_upper"]) < alpha
    )
    primary["minimum_evaluable_species_required"] = min_species
    primary["alpha"] = alpha
    primary["support"] = support

    for result in frame_results[1:]:
        if result.get("status") == "complete_evaluable":
            result["nominal_positive"] = bool(
                float(result["observed_mean_spearman_rho"]) > 0.0
                and float(result["p_upper"]) < alpha
            )
            result["can_rescue_primary"] = False

    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    primary_species.to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    pd.DataFrame(
        {
            "permutation": np.arange(len(primary_null), dtype=int),
            "null_equal_species_mean_spearman_rho": primary_null,
        }
    ).to_csv(OUT_NULL, index=False, lineterminator="\n")
    pd.DataFrame(frame_results[1:]).to_csv(OUT_SENS, index=False, lineterminator="\n")

    h4b_open = bool(support)
    payload = {
        "protocol": contract["protocol"],
        "status": "complete_exploratory_h4a_evaluable" if evaluable else "not_evaluable_exploratory_h4a",
        "introduced_after_h1_h2_outcomes": True,
        "claim_role": contract["claim_role"],
        "primary": primary,
        "sensitivities": frame_results[1:],
        "decision": (
            "support_within_species_climate_colour_divergence_open_h4b"
            if support
            else "no_support_within_species_climate_colour_divergence_do_not_open_h4b"
            if evaluable
            else "not_evaluable_do_not_open_h4b"
        ),
        "h4b_directional_decomposition_opened": h4b_open,
        "colour_representation": {
            "ordered_morph_axis_used": False,
            "soft_four_group_species_cell_composition": True,
            "pairwise_colour_distance": "Jensen-Shannon divergence in bits",
        },
        "climate_representation": {
            "variables": CLIMATE_Z,
            "pairwise_distance": "RMS standardized macroclimate separation",
            "source_scale_km": 250,
        },
        "frozen_primary_decisions_unchanged": {
            "h1": h1["decision"],
            "h1_p_upper": float(h1["primary"]["p_upper"]),
            "h2": h2["hierarchical_decision"],
        },
        "lineage": {
            "contract_sha256": sha256(CONTRACT),
            "measurement_table_sha256": sha256(MEASURED),
            "climate_cells_sha256": sha256(CLIMATE),
            "capacity_table_sha256": sha256(CAPACITY),
            "permutations": permutations,
            "seed": seed,
        },
        "files": {
            "primary_species": str(OUT_SPECIES.relative_to(ROOT)),
            "primary_null": str(OUT_NULL.relative_to(ROOT)),
            "sensitivities": str(OUT_SENS.relative_to(ROOT)),
        },
    }
    OUT_RESULT.parent.mkdir(parents=True, exist_ok=True)
    OUT_RESULT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
