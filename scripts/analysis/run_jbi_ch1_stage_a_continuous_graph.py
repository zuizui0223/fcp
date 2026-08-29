#!/usr/bin/env python3
"""Run frozen Stage-A species-conditioned continuous colour graph null."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.continuous_colour_boundaries import (
    edge_colour_discontinuity,
    weighted_graph_discontinuity,
)
from fcp_pipeline.spatial_graph import spherical_knn_edges

PROTOCOL = "jbi-ch1-stage-a-continuous-graph-v1"
EXPECTED_SPECIES = [
    "Antirrhinum majus",
    "Dactylorhiza sambucina",
    "Gentiana lutea",
    "Ipomoea purpurea",
    "Lysimachia arvensis",
    "Raphanus sativus",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def mc_summary(observed: float, null: np.ndarray) -> dict:
    null = np.asarray(null, dtype=float)
    mean = float(null.mean())
    sd = float(null.std(ddof=1))
    lower_p = float((1 + np.count_nonzero(null <= observed)) / (len(null) + 1))
    observed_dev = abs(observed - mean)
    two_p = float((1 + np.count_nonzero(np.abs(null - mean) >= observed_dev)) / (len(null) + 1))
    q = np.quantile(null, [0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995])
    return {
        "observed": float(observed),
        "null_mean": mean,
        "null_sd": sd,
        "clustering_deficit": float(mean - observed),
        "standardized_clustering_deficit": float((mean - observed) / sd) if sd > 0 else None,
        "p_lower_tail": lower_p,
        "p_two_sided_descriptive": two_p,
        "null_quantiles": {
            "p005": float(q[0]),
            "p025": float(q[1]),
            "p05": float(q[2]),
            "p50": float(q[3]),
            "p95": float(q[4]),
            "p975": float(q[5]),
            "p995": float(q[6]),
        },
    }


def run_k(
    rows: list[dict],
    coord_by_photo: dict[str, tuple[str, float, float]],
    representation: dict,
    *,
    k: int,
    n_permutations: int,
    seed: int,
    min_rows: int,
) -> tuple[dict, np.ndarray, list[str]]:
    species_order = list(EXPECTED_SPECIES)
    observed_species = np.empty(len(species_order), dtype=float)
    null_species = np.empty((n_permutations, len(species_order)), dtype=float)
    species_details = {}

    for j, species in enumerate(species_order):
        group = [
            r for r in rows
            if r.get("species") == species
            and r.get("feature_status") == "ok"
            and r.get("continuous_colour_vector_z") is not None
        ]
        group = sorted(group, key=lambda r: str(r["photo_id"]))
        if len(group) < min_rows:
            raise ValueError(f"{species}: only {len(group)} evaluable rows; minimum is {min_rows}")

        expected_dim = len(representation["per_species"][species]["feature_names"])
        values = np.asarray([r["continuous_colour_vector_z"] for r in group], dtype=float)
        if values.ndim != 2 or values.shape != (len(group), expected_dim):
            raise ValueError(f"{species}: continuous vector dimension mismatch")
        if not np.isfinite(values).all():
            raise ValueError(f"{species}: non-finite standardized colour vector")

        coords = []
        for r in group:
            photo_id = str(r["photo_id"])
            if photo_id not in coord_by_photo:
                raise ValueError(f"missing frozen coordinate for photo {photo_id}")
            coord_species, lat, lon = coord_by_photo[photo_id]
            if coord_species != species:
                raise ValueError(f"species mismatch for photo {photo_id}")
            coords.append((lat, lon))
        coords_arr = np.asarray(coords, dtype=float)
        edges, edge_distance = spherical_knn_edges(
            coords_arr[:, 0],
            coords_arr[:, 1],
            k=k,
            max_edge_km=None,
        )

        observed_scores = edge_colour_discontinuity(values, edges)
        observed_q = weighted_graph_discontinuity(observed_scores)
        observed_species[j] = observed_q

        # Independent deterministic stream per species and graph k. Complete rows are
        # shuffled; coordinates/edges remain fixed for every replicate.
        rng = np.random.default_rng(seed + 100_000 * k + 1_009 * (j + 1))
        for b in range(n_permutations):
            permuted = values[rng.permutation(len(values))]
            null_scores = edge_colour_discontinuity(permuted, edges)
            null_species[b, j] = weighted_graph_discontinuity(null_scores)

        species_details[species] = {
            "n_evaluation_total": sum(r.get("species") == species for r in rows),
            "n_measurement_evaluable": len(group),
            "n_measurement_not_evaluable": sum(r.get("species") == species for r in rows) - len(group),
            "feature_dimension": expected_dim,
            "n_graph_edges": int(len(edges)),
            "edge_distance_km": {
                "min": float(np.min(edge_distance)),
                "median": float(np.median(edge_distance)),
                "p95": float(np.quantile(edge_distance, 0.95)),
                "max": float(np.max(edge_distance)),
            },
            "q": mc_summary(observed_q, null_species[:, j]),
        }

    observed_global = float(observed_species.mean())
    null_global = null_species.mean(axis=1)
    result = {
        "k": int(k),
        "max_edge_km": None,
        "edge_weights": "uniform",
        "species_weighting": "equal arithmetic mean after within-species normalization",
        "species": species_details,
        "global_equal_species_mean_q": mc_summary(observed_global, null_global),
    }
    return result, np.column_stack([null_global, null_species]), species_order


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/evaluation/jbi_ch1_florence_evaluation_features_v1.jsonl"))
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument("--representation", type=Path, default=Path("docs/supporting/jbi_ch1_continuous_colour_representation_v1.json"))
    parser.add_argument("--contract", type=Path, default=Path("docs/supporting/jbi_ch1_stage_a_continuous_graph_contract_v1.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/supporting/jbi_ch1_stage_a_continuous_graph_v1.json"))
    parser.add_argument("--primary-null-csv", type=Path, default=Path("data/evaluation/jbi_ch1_stage_a_primary_null_v1.csv"))
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    representation = json.loads(args.representation.read_text(encoding="utf-8"))
    if contract.get("protocol") != PROTOCOL or contract.get("status") != "frozen_before_any_evaluation_spatial_colour_result":
        raise ValueError("Stage A contract is not frozen")
    if representation.get("status") != "frozen_before_evaluation_values_inspected":
        raise ValueError("continuous representation is not frozen")

    rows = load_jsonl(args.features)
    if len(rows) != 720 or len({str(r["photo_id"]) for r in rows}) != 720:
        raise ValueError("expected exactly 720 unique evaluation feature records")
    if any(r.get("evaluation_row") is not True or r.get("final_label") is not False for r in rows):
        raise ValueError("evaluation feature contract violation")
    if sorted({str(r["species"]) for r in rows}) != EXPECTED_SPECIES:
        raise ValueError("unexpected evaluation species set")

    split = pd.read_csv(args.split)
    eval_split = split.loc[split["split"].astype(str).eq("evaluation")].copy()
    if len(eval_split) != 720:
        raise ValueError("frozen split does not contain 720 evaluation rows")
    coord_by_photo = {
        str(row.photo_id): (str(row.species), float(row.latitude), float(row.longitude))
        for row in eval_split.itertuples(index=False)
    }
    if set(coord_by_photo) != {str(r["photo_id"]) for r in rows}:
        raise ValueError("evaluation feature IDs do not exactly match frozen evaluation coordinates")

    n_permutations = int(contract["null"]["permutations"])
    seed = int(contract["null"]["random_seed"])
    min_rows = int(contract["input"]["minimum_evaluable_rows_per_species"])
    primary_k = int(contract["geometry"]["primary_k"])
    sensitivity_k = [int(x) for x in contract["geometry"]["predeclared_sensitivity_k"]]

    analyses = {}
    primary_null = None
    primary_species_order = None
    for k in [primary_k, *sensitivity_k]:
        result, null_matrix, species_order = run_k(
            rows,
            coord_by_photo,
            representation,
            k=k,
            n_permutations=n_permutations,
            seed=seed,
            min_rows=min_rows,
        )
        analyses[str(k)] = result
        if k == primary_k:
            primary_null = null_matrix
            primary_species_order = species_order

    assert primary_null is not None and primary_species_order is not None
    args.primary_null_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.primary_null_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["permutation", "global_equal_species_mean_q", *primary_species_order])
        for b, values in enumerate(primary_null, start=1):
            writer.writerow([b, *[format(float(v), ".12g") for v in values]])

    primary = analyses[str(primary_k)]["global_equal_species_mean_q"]
    result = {
        "protocol": PROTOCOL,
        "status": "stage_a_evaluation_complete",
        "primary_representation": "species-specific continuous colour vector standardized from frozen calibration parameters",
        "n_evaluation_records": 720,
        "n_permutations": n_permutations,
        "random_seed": seed,
        "primary_k": primary_k,
        "sensitivity_k": sensitivity_k,
        "primary_global_result": primary,
        "primary_direction": "lower Q indicates greater local colour similarity than species-conditioned random labelling",
        "primary_rejects_random_labelling_at_0_05": bool(primary["p_lower_tail"] <= 0.05),
        "analyses_by_k": analyses,
        "evaluation_feature_sha256": sha256(args.features),
        "frozen_split_sha256": sha256(args.split),
        "frozen_representation_sha256": sha256(args.representation),
        "frozen_contract_sha256": sha256(args.contract),
        "primary_null_csv_sha256": None,
        "environment_used": False,
        "geographic_reference_library_used": False,
        "interpretation_limit": "Stage A tests within-species local colour organization only; it does not establish cross-species shared boundaries or causes.",
        "next_gate": "if spatial organization is supported, construct species-specific transition-intensity surfaces and test shared-boundary concentration under the complete within-species permutation pipeline"
    }
    result["primary_null_csv_sha256"] = sha256(args.primary_null_csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
