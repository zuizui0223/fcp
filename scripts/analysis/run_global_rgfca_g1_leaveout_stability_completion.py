#!/usr/bin/env python3
"""Complete the preoutcome G1 realm/family leave-out stability gates.

The frozen primary G1 decision is immutable. This script keeps the exact original
200-outer species/photo schedule and removes group contributions without replacing
omitted species. It uses the same 999 within-species complete-colour-vector nulls.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from fcp_pipeline.global_barrier_field import equal_area_grid_centers
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_rgfca_engine import (
    COLOUR_COLUMNS,
    _rank_edges,
    _rank_edges_matrix,
    _raw_jsd_from_source_rows,
    _raw_jsd_matrix_from_source_rows,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
    prepare_sparse_outer_geometry,
)

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_g1_leaveout_stability_completion_contract_v1.json"
AUTH = ROOT / "docs/supporting/global_rgfca_g1_leaveout_stability_completion_authorization_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
METADATA = ROOT / "data/frozen/global_rgfca_g5_species_metadata_v1.csv"
PRIMARY = ROOT / "docs/supporting/global_rgfca_g1_result_v1.json"


@dataclass(frozen=True)
class Scenario:
    kind: str
    group: str
    members: frozenset[str]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--null-batch-size", type=int, default=64)
    return p.parse_args()


def _bool_flag(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.casefold().isin({"true", "1", "yes"})


def load_pool():
    frame = pd.read_csv(MEASURED)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"measured pool missing {missing}")
    frame = frame[_bool_flag(frame["global_classifiable"])].copy()
    counts = frame.groupby("species", observed=True).size()
    eligible = set(counts[counts >= 40].index.astype(str))
    frame = frame[frame["species"].astype(str).isin(eligible)].copy()
    pool = canonical_colour_pool(frame[["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]])
    if len(pool.species_labels) != 369 or len(pool.photo_ids) != 21424:
        raise RuntimeError(f"G1 pool drift: {len(pool.species_labels)} species / {len(pool.photo_ids)} rows")
    return pool


def clean_label(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.casefold() in {"nan", "none", "na", "null"}:
        return None
    return text


def build_scenarios(species_labels: tuple[str, ...]) -> tuple[list[Scenario], dict[str, object]]:
    auth = json.loads(AUTH.read_text(encoding="utf-8"))
    expected_sha = str(auth.get("metadata_source_sha256") or "")
    actual_sha = sha256_file(METADATA)
    if not expected_sha or actual_sha != expected_sha:
        raise RuntimeError(f"G5 metadata SHA mismatch: {actual_sha} != {expected_sha}")
    meta = pd.read_csv(METADATA, usecols=["species", "accepted_family", "dominant_realm"])
    meta["species"] = meta["species"].astype(str)
    if meta["species"].duplicated().any():
        raise RuntimeError("G5 metadata has duplicate species")
    meta = meta.set_index("species").reindex(list(species_labels))
    if len(meta) != 369:
        raise RuntimeError("metadata reindex drift")

    realms: dict[str, set[str]] = {}
    families: dict[str, set[str]] = {}
    missing_realm = 0
    missing_family = 0
    for species, row in meta.iterrows():
        realm = clean_label(row["dominant_realm"])
        family = clean_label(row["accepted_family"])
        if realm is None:
            missing_realm += 1
        else:
            realms.setdefault(realm, set()).add(str(species))
        if family is None:
            missing_family += 1
        else:
            families.setdefault(family, set()).add(str(species))

    scenarios: list[Scenario] = []
    for realm in sorted(realms):
        scenarios.append(Scenario("realm", realm, frozenset(realms[realm])))
    major_threshold = math.ceil(0.05 * 369)
    major = {name: members for name, members in families.items() if len(members) >= major_threshold}
    for family in sorted(major):
        scenarios.append(Scenario("major_family", family, frozenset(major[family])))
    if not scenarios:
        raise RuntimeError("no leave-out scenarios were defined")
    audit = {
        "metadata_sha256": actual_sha,
        "realm_groups": {name: len(realms[name]) for name in sorted(realms)},
        "missing_realm_species": missing_realm,
        "major_family_threshold_species": major_threshold,
        "major_family_groups": {name: len(major[name]) for name in sorted(major)},
        "missing_family_species": missing_family,
        "scenario_count": len(scenarios),
    }
    return scenarios, audit


def weighted_concentration(field: np.ndarray, opportunity: np.ndarray) -> float:
    keep = np.isfinite(field) & np.isfinite(opportunity) & (opportunity > 0)
    if not np.any(keep):
        return float("nan")
    x = field[keep]
    w = opportunity[keep]
    mean = float(np.average(x, weights=w))
    return float(np.average(np.square(x - mean), weights=w))


def group_matrix_and_masks(geometry, outer_labels: np.ndarray, scenarios: list[Scenario], species_to_scenarios: dict[str, tuple[int, ...]]):
    n_scen = len(scenarios)
    n_cells = geometry.grid.n_cells
    support_decrement = np.zeros((n_scen, n_cells), dtype=np.int16)
    omitted_slot = np.zeros((n_scen, len(outer_labels)), dtype=bool)
    for slot, label_obj in enumerate(outer_labels):
        label = str(label_obj)
        scen_ids = species_to_scenarios.get(label, ())
        if not scen_ids:
            continue
        start, stop = geometry.edge_species_slices[slot]
        contribution = np.asarray(geometry.weighted_kernel[start:stop].sum(axis=0)).ravel() > 0
        for s in scen_ids:
            omitted_slot[s, slot] = True
            support_decrement[s] += contribution.astype(np.int16)

    blocks = []
    group_opp = np.zeros((n_scen, n_cells), dtype=float)
    for s in range(n_scen):
        edge_mask = omitted_slot[s, geometry.edge_species_index]
        if np.any(edge_mask):
            block = geometry.weighted_kernel.multiply(edge_mask[:, None]).tocsr()
            block.eliminate_zeros()
            group_opp[s] = np.asarray(block.sum(axis=0)).ravel()
        else:
            block = sparse.csr_matrix(geometry.weighted_kernel.shape, dtype=float)
        blocks.append(block)
    group_matrix = sparse.hstack(blocks, format="csr")
    retained_opp = geometry.opportunity[None, :] - group_opp
    retained_opp[np.abs(retained_opp) < 1e-15] = 0.0
    retained_support = geometry.distinct_species_support[None, :].astype(np.int32) - support_decrement.astype(np.int32)
    evaluable = (retained_support >= 5) & (retained_opp > 0)
    return group_matrix, retained_opp, evaluable


def main() -> int:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    auth = json.loads(AUTH.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    if contract.get("status") != "technical_completion_frozen_before_any_leaveout_colour_result_is_computed":
        raise RuntimeError("leave-out contract drift")
    if auth.get("status") != "authorize_exactly_one_g1_leaveout_stability_completion":
        raise RuntimeError("leave-out authorization drift")
    if primary.get("p_upper") != 0.07 or primary.get("g1_supported") is not False:
        raise RuntimeError("primary G1 reference drift")
    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 128:
        raise ValueError("null batch size must lie in [1,128]")

    pool = load_pool()
    scenarios, metadata_audit = build_scenarios(pool.species_labels)
    n_scen = len(scenarios)
    species_to_scenarios: dict[str, tuple[int, ...]] = {}
    for species in pool.species_labels:
        ids = [i for i, scenario in enumerate(scenarios) if species in scenario.members]
        species_to_scenarios[species] = tuple(ids)

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
    cache = build_pairwise_jsd_cache(pool)
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    null_indices = np.arange(999, dtype=np.int64)
    print(json.dumps({"stage":"build_null_source_rows","nulls":999,"pool_rows":len(pool.photo_ids)}), flush=True)
    null_source = null_source_row_matrix(pool, null_indices, master_seed=2026090403)
    grid = equal_area_grid_centers(36, 18)

    observed_agg_num = np.zeros((n_scen, grid.n_cells), dtype=float)
    aggregate_opp = np.zeros((n_scen, grid.n_cells), dtype=float)
    null_agg_num = np.zeros((999, n_scen, grid.n_cells), dtype=np.float64)
    max_evaluable_outer_cells = np.zeros(n_scen, dtype=np.int32)

    for outer in range(200):
        geometry = prepare_sparse_outer_geometry(
            pool,
            schedule.outer_photo_ids[outer],
            schedule.outer_species[outer],
            grid=grid,
            k=3,
            kernel_km=500.0,
            cutoff_multiplier=3.0,
            minimum_distinct_species=5,
        )
        group_matrix, retained_opp, evaluable = group_matrix_and_masks(
            geometry, schedule.outer_species[outer], scenarios, species_to_scenarios
        )
        max_evaluable_outer_cells = np.maximum(max_evaluable_outer_cells, evaluable.sum(axis=1).astype(np.int32))
        aggregate_opp += np.where(evaluable, retained_opp, 0.0)

        raw_obs = _raw_jsd_from_source_rows(cache, geometry.edge_nodes, identity)
        rank_obs = _rank_edges(raw_obs, geometry.edge_species_slices)
        full_obs_num = np.asarray(geometry.weighted_kernel.T.dot(rank_obs)).ravel()
        group_obs_flat = np.asarray(group_matrix.T.dot(rank_obs)).ravel()
        group_obs_num = group_obs_flat.reshape(n_scen, grid.n_cells)
        retained_obs_num = full_obs_num[None, :] - group_obs_num
        observed_agg_num += np.where(evaluable, retained_obs_num, 0.0)

        for start in range(0, 999, batch_size):
            stop = min(start + batch_size, 999)
            source = null_source[start:stop]
            raw_null = _raw_jsd_matrix_from_source_rows(cache, geometry.edge_nodes, source)
            rank_null = _rank_edges_matrix(raw_null, geometry.edge_species_slices)
            full_null_num = np.asarray(geometry.weighted_kernel.T.dot(rank_null.T).T)
            group_null_flat = np.asarray(group_matrix.T.dot(rank_null.T).T)
            group_null_num = group_null_flat.reshape(stop - start, n_scen, grid.n_cells)
            retained_null = full_null_num[:, None, :] - group_null_num
            null_agg_num[start:stop] += np.where(evaluable[None, :, :], retained_null, 0.0)
        if (outer + 1) % 10 == 0 or outer == 0:
            print(json.dumps({"stage":"outer_complete","outer":outer+1,"total":200,"scenarios":n_scen}), flush=True)

    rows = []
    all_positive = True
    n_estimable = 0
    for s, scenario in enumerate(scenarios):
        opp = aggregate_opp[s]
        keep = opp > 0
        estimable = bool(np.any(keep) and max_evaluable_outer_cells[s] > 0)
        result = {
            "kind": scenario.kind,
            "group": scenario.group,
            "omitted_species_count": len(scenario.members),
            "estimable": estimable,
            "evaluable_cells": int(np.count_nonzero(keep)),
            "maximum_evaluable_cells_in_any_outer": int(max_evaluable_outer_cells[s]),
        }
        if estimable:
            n_estimable += 1
            obs_field = np.full(grid.n_cells, np.nan, dtype=float)
            obs_field[keep] = observed_agg_num[s, keep] / opp[keep]
            obs_conc = weighted_concentration(obs_field, opp)
            null_conc = np.empty(999, dtype=float)
            for j in range(999):
                field = np.full(grid.n_cells, np.nan, dtype=float)
                field[keep] = null_agg_num[j, s, keep] / opp[keep]
                null_conc[j] = weighted_concentration(field, opp)
            null_mean = float(np.mean(null_conc))
            excess = float(obs_conc - null_mean)
            positive = bool(excess > 0)
            all_positive = all_positive and positive
            result.update({
                "observed_concentration": float(obs_conc),
                "null_mean_concentration": null_mean,
                "G1_excess": excess,
                "G1_excess_ratio": float(obs_conc / null_mean) if null_mean > 0 else None,
                "positive_excess_gate": positive,
                "descriptive_p_upper": float((1 + np.count_nonzero(null_conc >= obs_conc)) / 1000),
                "null_q025": float(np.quantile(null_conc, 0.025)),
                "null_q975": float(np.quantile(null_conc, 0.975)),
            })
        else:
            result.update({
                "observed_concentration": None,
                "null_mean_concentration": None,
                "G1_excess": None,
                "G1_excess_ratio": None,
                "positive_excess_gate": None,
                "descriptive_p_upper": None,
                "null_q025": None,
                "null_q975": None,
            })
        rows.append(result)

    out = pd.DataFrame(rows).sort_values(["kind", "group"], kind="mergesort").reset_index(drop=True)
    realm_estimable = out[(out["kind"] == "realm") & out["estimable"]]
    family_estimable = out[(out["kind"] == "major_family") & out["estimable"]]
    realm_gate = bool(len(realm_estimable) > 0 and realm_estimable["positive_excess_gate"].eq(True).all())
    family_gate = bool(len(family_estimable) > 0 and family_estimable["positive_excess_gate"].eq(True).all())

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "global_rgfca_g1_leaveout_stability_completion_groups_v1.csv"
    out.to_csv(csv_path, index=False)
    payload = {
        "protocol": contract["protocol"],
        "status": "complete_prespecified_g1_leaveout_stability_gates",
        "primary_g1": {"p_upper": 0.07, "supported": False, "decision_immutable": True},
        "metadata_audit": metadata_audit,
        "scenarios_total": n_scen,
        "scenarios_estimable": n_estimable,
        "realm_estimable_groups": int(len(realm_estimable)),
        "major_family_estimable_groups": int(len(family_estimable)),
        "leave_one_realm_positive_excess_every_estimable": realm_gate,
        "leave_one_major_family_positive_excess_every_estimable": family_gate,
        "both_leaveout_components_pass": bool(realm_gate and family_gate),
        "known_running_100_vs_200_gate_pass": False,
        "full_original_strong_stability_claim_pass": False,
        "reason_full_strong_stability_false": "The predeclared running-consensus 100-vs-200 field correlation was 0.8809461021778187, below the required 0.9, independently of these leave-out results.",
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "schedule_replacement_after_omission": False,
        "selective_stopping_used": False,
        "claim_boundary": contract["claim_ceiling"],
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "authorization_sha256": sha256_file(AUTH),
            "metadata_sha256": sha256_file(METADATA),
            "measured_table_sha256": sha256_file(MEASURED),
            "primary_g1_sha256": sha256_file(PRIMARY),
            "group_result_sha256": sha256_file(csv_path),
        }
    }
    result_path = args.output_dir / "global_rgfca_g1_leaveout_stability_completion_result_v1.json"
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
