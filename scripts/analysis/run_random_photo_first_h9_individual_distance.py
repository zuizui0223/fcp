#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h9_individual_distance_contract_v1.json"
EXECUTION = ROOT / "docs/supporting/random_photo_first_h9_inference_execution_v1.json"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_h9_measurement_result_v1.json"
MEASURED = ROOT / "data/derived/random_photo_first_h9_measured_photos_v1.csv"
SUPPORT = ROOT / "data/derived/random_photo_first_h9_measurement_species_support_v1.csv"

OUT_SELECTED = ROOT / "data/derived/random_photo_first_h9_selected_20_per_species_v1.csv"
OUT_SPECIES = ROOT / "data/derived/random_photo_first_h9_species_spatial_result_v1.csv"
OUT_NULL = ROOT / "data/derived/random_photo_first_h9_global_null_v1.csv"
OUT_RESULT = ROOT / "docs/supporting/random_photo_first_h9_individual_distance_result_v1.json"

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
FOUR = ["soft_white", "soft_yellow_orange", "soft_red_pink", "soft_blue_purple"]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selection_digest(measurement_id: str, seed: int) -> str:
    token = str(measurement_id).encode("utf-8") + str(int(seed)).encode("utf-8")
    return hashlib.sha256(token).hexdigest()


def jsd_bits(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    out = 0.0
    kp = p > 0
    if kp.any():
        out += 0.5 * float(np.sum(p[kp] * np.log2(p[kp] / m[kp])))
    kq = q > 0
    if kq.any():
        out += 0.5 * float(np.sum(q[kq] * np.log2(q[kq] / m[kq])))
    return float(np.clip(out, 0.0, 1.0))


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx = x - float(np.mean(x))
    dy = y - float(np.mean(y))
    nx = float(np.sqrt(np.sum(dx * dx)))
    ny = float(np.sqrt(np.sum(dy * dy)))
    if nx <= 0.0 or ny <= 0.0:
        return float("nan")
    return float(np.clip(np.sum(dx * dy) / (nx * ny), -1.0, 1.0))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(
        rankdata(np.asarray(x, dtype=float), method="average"),
        rankdata(np.asarray(y, dtype=float), method="average"),
    )


def fixed_rank_unit(x: np.ndarray) -> np.ndarray:
    r = rankdata(np.asarray(x, dtype=float), method="average")
    r = r - float(np.mean(r))
    n = float(np.sqrt(np.sum(r * r)))
    if n <= 0.0:
        raise RuntimeError("constant geographic rank vector")
    return r / n


def spearman_against_fixed_unit(fixed: np.ndarray, y: np.ndarray) -> float:
    r = rankdata(np.asarray(y, dtype=float), method="average")
    r = r - float(np.mean(r))
    n = float(np.sqrt(np.sum(r * r)))
    if n <= 0.0:
        return float("nan")
    return float(np.clip(np.sum(fixed * (r / n)), -1.0, 1.0))


def great_circle_km(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    lat1 = np.radians(np.asarray(lat1, dtype=float))
    lat2 = np.radians(np.asarray(lat2, dtype=float))
    lon1 = np.radians(np.asarray(lon1, dtype=float))
    lon2 = np.radians(np.asarray(lon2, dtype=float))
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def child_permutation(n: int, *, seed: int, permutation: int, taxon_id: str) -> np.ndarray:
    token = f"{seed}|{permutation}|{taxon_id}".encode("utf-8")
    child = int.from_bytes(hashlib.sha256(token).digest()[:8], "big", signed=False)
    return np.random.default_rng(child).permutation(n)


def four_group_matrix(frame: pd.DataFrame) -> np.ndarray:
    raw = frame[PALETTE].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).all():
        raise RuntimeError("non-finite palette fractions in H9 selected frame")
    four = np.column_stack([
        raw[:, 0],
        raw[:, 1] + raw[:, 2] + raw[:, 8],
        raw[:, 3] + raw[:, 4] + raw[:, 5],
        raw[:, 7] + raw[:, 6],
    ])
    mass = four.sum(axis=1)
    if np.any(mass <= 0.0):
        raise RuntimeError("zero biological four-group colour mass in H9 selected frame")
    return four / mass[:, None]


def colour_distance_matrix(colour: np.ndarray) -> np.ndarray:
    n = len(colour)
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = jsd_bits(colour[i], colour[j])
            out[i, j] = d
            out[j, i] = d
    return out


def main() -> int:
    contract = load_json(CONTRACT)
    execution = load_json(EXECUTION)
    measurement = load_json(MEASUREMENT_RESULT)

    if contract.get("status") != "prospective_h9_frozen_after_metadata_feasibility_before_any_h9_query_or_pixel":
        raise RuntimeError("H9 parent contract state drifted")
    if execution.get("status") != "frozen_after_h9_measurement_gate_before_any_h9_spatial_statistic":
        raise RuntimeError("H9 inference execution contract is not frozen pre-statistic")
    if measurement.get("status") != execution["required_measurement_status"]:
        raise RuntimeError("H9 measurement result status drifted")
    gate = measurement.get("postmeasurement_gate", {})
    if gate.get("pass") is not True or int(gate.get("evaluable_species", -1)) != 34:
        raise RuntimeError("H9 postmeasurement gate is not the frozen 34-species pass")
    if measurement.get("h9_spatial_inference_run") is not False:
        raise RuntimeError("H9 spatial inference was already marked as run")
    if sha256_file(MEASURED) != execution["required_measured_table_sha256"]:
        raise RuntimeError("H9 measured table SHA drifted")
    if sha256_file(SUPPORT) != execution["required_species_support_sha256"]:
        raise RuntimeError("H9 species-support SHA drifted")

    analysis = contract["analysis"]
    photos_per_species = int(analysis["photos_per_evaluable_species"])
    analysis_seed = int(analysis["analysis_seed"])
    permutations = int(analysis["permutations"])
    permutation_seed = int(analysis["permutation_seed"])
    if (photos_per_species, analysis_seed, permutations, permutation_seed) != (20, 20260909, 999, 20260910):
        raise RuntimeError("H9 frozen analysis settings drifted")

    measured = pd.read_csv(MEASURED, dtype={"measurement_id": str, "inat_taxon_id": str})
    required = {"measurement_id", "species", "inat_taxon_id", "latitude", "longitude", "measurement_status", "morph", *PALETTE}
    missing = sorted(required.difference(measured.columns))
    if missing:
        raise RuntimeError(f"H9 measurement table lacks required fields: {missing}")
    if len(measured) != 2280 or measured["measurement_id"].nunique() != 2280:
        raise RuntimeError("H9 measurement denominator drifted")
    if measured["inat_taxon_id"].nunique() != 38:
        raise RuntimeError("H9 raw species denominator drifted")

    classifiable = measured.loc[
        measured["measurement_status"].astype(str).eq("classified_four_state_morph")
        & measured["morph"].astype(str).isin(BIOLOGICAL_MORPHS)
    ].copy()
    if len(classifiable) != 1244:
        raise RuntimeError(f"H9 classifiable denominator drifted: {len(classifiable)}")

    support = pd.read_csv(SUPPORT, dtype={"inat_taxon_id": str})
    eligible_ids = set(support.loc[support["h9_measurement_evaluable"].astype(str).str.lower().isin({"true", "1"}), "inat_taxon_id"].astype(str))
    if len(eligible_ids) != 34:
        raise RuntimeError(f"H9 evaluable species denominator drifted: {len(eligible_ids)}")

    classifiable = classifiable.loc[classifiable["inat_taxon_id"].astype(str).isin(eligible_ids)].copy()
    classifiable["selection_hash"] = [selection_digest(x, analysis_seed) for x in classifiable["measurement_id"].astype(str)]
    selected_parts = []
    for taxon_id, group in classifiable.groupby("inat_taxon_id", sort=True, observed=True):
        ordered = group.sort_values(["selection_hash", "measurement_id"], kind="mergesort").head(photos_per_species).copy()
        if len(ordered) != photos_per_species:
            raise RuntimeError(f"H9 taxon {taxon_id} lacks exactly 20 selected classifiable photos")
        ordered["analysis_rank"] = np.arange(1, photos_per_species + 1, dtype=int)
        selected_parts.append(ordered)
    selected = pd.concat(selected_parts, ignore_index=True)
    if len(selected) != 34 * 20 or selected["inat_taxon_id"].nunique() != 34:
        raise RuntimeError("H9 exact fixed-n selected frame drifted")

    species_records: list[dict[str, Any]] = []
    null_rows = np.empty((34, permutations), dtype=float)
    species_order: list[str] = []

    for row_idx, (taxon_id, group) in enumerate(selected.groupby("inat_taxon_id", sort=True, observed=True)):
        taxon_id = str(taxon_id)
        group = group.sort_values("analysis_rank", kind="mergesort").reset_index(drop=True)
        colour = four_group_matrix(group)
        for i, name in enumerate(FOUR):
            group[name] = colour[:, i]
        ii, jj = np.triu_indices(photos_per_species, k=1)
        lat = pd.to_numeric(group["latitude"], errors="raise").to_numpy(dtype=float)
        lon = pd.to_numeric(group["longitude"], errors="raise").to_numpy(dtype=float)
        geo = great_circle_km(lat[ii], lon[ii], lat[jj], lon[jj])
        cmat = colour_distance_matrix(colour)
        colour_dist = cmat[ii, jj]
        rho = spearman(geo, colour_dist)
        if not np.isfinite(rho):
            raise RuntimeError(f"H9 species statistic non-finite for taxon {taxon_id}")
        fixed = fixed_rank_unit(geo)
        for p in range(permutations):
            order = child_permutation(photos_per_species, seed=permutation_seed, permutation=p, taxon_id=taxon_id)
            prho = spearman_against_fixed_unit(fixed, cmat[order[ii], order[jj]])
            if not np.isfinite(prho):
                raise RuntimeError(f"H9 permutation statistic non-finite for taxon {taxon_id}, p={p}")
            null_rows[row_idx, p] = prho
        species_order.append(taxon_id)
        species_records.append({
            "species": str(group["species"].iloc[0]),
            "inat_taxon_id": taxon_id,
            "analysis_photos": photos_per_species,
            "photo_pairs": int(len(ii)),
            "geographic_colour_spearman_rho": float(rho),
            "mean_pairwise_great_circle_km": float(np.mean(geo)),
            "maximum_pairwise_great_circle_km": float(np.max(geo)),
            "mean_pairwise_colour_jsd_bits": float(np.mean(colour_dist)),
        })

    species_df = pd.DataFrame(species_records).sort_values("inat_taxon_id", kind="mergesort").reset_index(drop=True)
    if species_order != species_df["inat_taxon_id"].astype(str).tolist():
        raise RuntimeError("H9 species ordering drifted during null construction")
    observed = float(species_df["geographic_colour_spearman_rho"].mean())
    global_null = null_rows.mean(axis=0)
    p_upper = float((1 + np.count_nonzero(global_null >= observed)) / (permutations + 1))
    supported = bool(observed > 0.0 and p_upper < float(analysis["alpha"]))
    decision = execution["decision"]["support_label"] if supported else execution["decision"]["no_support_label"]

    OUT_SELECTED.parent.mkdir(parents=True, exist_ok=True)
    selected_out = selected[["species", "inat_taxon_id", "measurement_id", "analysis_rank", "selection_hash"]].copy()
    selected_out.to_csv(OUT_SELECTED, index=False, lineterminator="\n")
    species_df.to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    pd.DataFrame({"permutation": np.arange(1, permutations + 1, dtype=int), "global_equal_species_mean_rho": global_null}).to_csv(
        OUT_NULL, index=False, lineterminator="\n"
    )

    result = {
        "protocol": contract["protocol"],
        "status": "complete_h9_individual_distance_inference",
        "measurement_gate_passed": True,
        "raw_species": 38,
        "evaluable_species": 34,
        "analysis_photos_per_species": photos_per_species,
        "analysis_photos": int(len(selected_out)),
        "photo_pairs_per_species": 190,
        "permutations": permutations,
        "analysis_seed": analysis_seed,
        "permutation_seed": permutation_seed,
        "primary": {
            "observed_equal_species_mean_spearman_rho": observed,
            "observed_median_species_spearman_rho": float(species_df["geographic_colour_spearman_rho"].median()),
            "positive_species_fraction": float((species_df["geographic_colour_spearman_rho"] > 0).mean()),
            "p_upper": p_upper,
            "null_mean": float(np.mean(global_null)),
            "null_median": float(np.median(global_null)),
            "null_q025": float(np.quantile(global_null, 0.025)),
            "null_q975": float(np.quantile(global_null, 0.975)),
            "support": supported,
            "decision": decision,
        },
        "claim_ceiling": contract["claim_ceiling"],
        "lineage": {
            "measured_table_sha256": sha256_file(MEASURED),
            "species_support_sha256": sha256_file(SUPPORT),
            "parent_contract_sha256": sha256_file(CONTRACT),
            "inference_execution_contract_sha256": sha256_file(EXECUTION),
            "selected_20_table_sha256": sha256_file(OUT_SELECTED),
            "species_result_sha256": sha256_file(OUT_SPECIES),
            "global_null_sha256": sha256_file(OUT_NULL),
            "new_image_pixels_opened": False,
            "new_measurement_run": False,
        },
        "files": {
            "selected_20": str(OUT_SELECTED.relative_to(ROOT)),
            "species_result": str(OUT_SPECIES.relative_to(ROOT)),
            "global_null": str(OUT_NULL.relative_to(ROOT)),
        },
    }
    OUT_RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
