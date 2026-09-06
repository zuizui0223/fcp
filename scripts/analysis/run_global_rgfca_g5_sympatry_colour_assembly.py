#!/usr/bin/env python3
"""Run the frozen primary G5 sympatry colour-assembly inference.

The matched-set file must already have been committed by the outcome-blind
metadata stage.  This script then opens the location-blind colour measurements,
reconstructs the exact RGFCA photo schedule, and evaluates the prespecified
conditional randomization test.
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_community_colour_assembly import matched_sympatry_colour_assembly_test
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule, schedule_audit
from fcp_pipeline.global_rgfca_engine import COLOUR_COLUMNS

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_g5_metadata_matching_implementation_contract_v1.json"
PARENT = ROOT / "docs/supporting/global_sympatry_colour_assembly_contract_v1.json"
MATCHED = ROOT / "data/frozen/global_rgfca_g5_matched_sets_v1.csv"
MATCH_MANIFEST = ROOT / "docs/supporting/global_rgfca_g5_matching_manifest_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def js_distance_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Jensen-Shannon distance (sqrt divergence, base 2) row-wise."""
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.shape != bb.shape or aa.ndim != 2 or aa.shape[1] != 4:
        raise ValueError("palette arrays must have matching shape (n,4)")
    m = 0.5 * (aa + bb)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(aa > 0, aa * np.log2(aa / m), 0.0).sum(axis=1)
        kb = np.where(bb > 0, bb * np.log2(bb / m), 0.0).sum(axis=1)
    divergence = np.clip(0.5 * (ka + kb), 0.0, 1.0)
    return np.sqrt(divergence)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    matching_manifest = json.loads(MATCH_MANIFEST.read_text(encoding="utf-8"))
    measurement = json.loads(MEASUREMENT.read_text(encoding="utf-8"))
    if matching_manifest.get("colour_outcome_opened") is not False or matching_manifest.get("matching_used_colour") is not False:
        raise RuntimeError("G5 outcome-blind matching firewall not satisfied")
    if sha256_file(MATCHED) != str(matching_manifest.get("output_sha256", {}).get("matched_sets") or ""):
        raise RuntimeError("G5 matched-set bytes differ from frozen matching manifest")
    if sha256_file(MEASURED) != str(measurement.get("lineage", {}).get("measured_table_sha256") or ""):
        raise RuntimeError("measured colour table lineage drift")

    matched = pd.read_csv(MATCHED)
    required_match = {"set_id", "focal_species", "sympatric_partner", "control_rank", "control_species"}
    missing = sorted(required_match - set(matched.columns))
    if missing:
        raise RuntimeError(f"frozen matched-set table lacks columns: {missing}")
    matched = matched.sort_values(["set_id", "control_rank"], kind="mergesort").reset_index(drop=True)

    frame = pd.read_csv(MEASURED)
    required_colour = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing_colour = sorted(required_colour - set(frame.columns))
    if missing_colour:
        raise RuntimeError(f"measured table lacks G5 colour inputs: {missing_colour}")
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    pool["species"] = pool["species"].astype(str)
    minimum_pool = int(contract["colour_schedule_implementation"]["frozen_parameters"]["minimum_pool_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = set(counts[counts >= minimum_pool].index.astype(str))
    pool = pool.loc[pool["species"].isin(eligible)].copy().reset_index(drop=True)
    if pool["species"].nunique() != int(contract["target_universe"]["expected_eligible_species"]):
        raise RuntimeError("G5 colour pool species count differs from frozen 369-species universe")

    colours = pool[list(COLOUR_COLUMNS)].apply(pd.to_numeric, errors="raise").to_numpy(float)
    mass = colours.sum(axis=1)
    if np.any(~np.isfinite(colours)) or np.any(colours < 0) or np.any(mass <= 0):
        raise RuntimeError("invalid four-group colour vectors")
    colours = colours / mass[:, None]
    colour_by_photo = {int(pid): colours[i] for i, pid in enumerate(pool["photo_id"].astype(np.int64).to_numpy())}

    params = contract["colour_schedule_implementation"]["frozen_parameters"]
    schedule = build_repeated_atlas_schedule(
        pool["photo_id"].astype(np.int64).to_numpy(),
        pool["species"].to_numpy(dtype=object),
        n_outer=int(params["n_outer"]),
        species_per_outer=int(params["species_per_outer"]),
        photos_per_species=int(params["photos_per_species"]),
        minimum_pool_photos_per_species=minimum_pool,
        species_seed=int(params["species_seed"]),
        photo_master_seed=int(params["photo_master_seed"]),
    )
    appearance: dict[str, list[np.ndarray]] = {label: [] for label in schedule.species_labels}
    for outer in range(schedule.n_outer):
        for slot in range(schedule.species_per_outer):
            label = str(schedule.outer_species[outer, slot])
            ids = schedule.outer_photo_ids[outer, slot]
            palette = np.mean(np.vstack([colour_by_photo[int(pid)] for pid in ids]), axis=0)
            palette = palette / palette.sum()
            appearance[label].append(palette)
    fixed_ranks = 135
    palette_by_species: dict[str, np.ndarray] = {}
    for label, rows in appearance.items():
        if len(rows) < fixed_ranks:
            raise RuntimeError(f"species {label} has fewer than 135 frozen schedule appearances")
        palette_by_species[label] = np.vstack(rows[:fixed_ranks])

    needed = set(matched["focal_species"].astype(str)) | set(matched["sympatric_partner"].astype(str)) | set(matched["control_species"].astype(str))
    absent = sorted(needed - set(palette_by_species))
    if absent:
        raise RuntimeError(f"matched species absent from colour schedule: {absent[:10]}")
    pair_cache: dict[tuple[str, str], float] = {}
    def pair_distance(a: str, b: str) -> float:
        key = tuple(sorted((str(a), str(b))))
        if key not in pair_cache:
            values = js_distance_rows(palette_by_species[key[0]], palette_by_species[key[1]])
            pair_cache[key] = float(np.mean(values))
        return pair_cache[key]

    focal_labels = sorted(matched["focal_species"].astype(str).unique().tolist())
    focal_code = {label: i for i, label in enumerate(focal_labels)}
    set_focal: list[int] = []
    set_sym: list[float] = []
    control_rows: list[list[float]] = []
    set_ids: list[str] = []
    set_details: list[dict[str, object]] = []
    for set_id, group in matched.groupby("set_id", sort=True):
        group = group.sort_values("control_rank", kind="mergesort")
        focal_values = group["focal_species"].astype(str).unique()
        partner_values = group["sympatric_partner"].astype(str).unique()
        if len(focal_values) != 1 or len(partner_values) != 1:
            raise RuntimeError(f"matched set {set_id} has inconsistent focal/partner")
        focal = str(focal_values[0])
        partner = str(partner_values[0])
        controls = group["control_species"].astype(str).tolist()
        sym_distance = pair_distance(focal, partner)
        control_distance = [pair_distance(focal, control) for control in controls]
        set_ids.append(str(set_id))
        set_focal.append(int(focal_code[focal]))
        set_sym.append(sym_distance)
        control_rows.append(control_distance)
        set_details.append(
            {
                "set_id": str(set_id),
                "focal_species": focal,
                "sympatric_partner": partner,
                "sympatric_colour_distance": sym_distance,
                "mean_allopatric_colour_distance": float(np.mean(control_distance)),
                "delta_sympatric_minus_allopatric": float(sym_distance - np.mean(control_distance)),
                "n_controls": len(control_distance),
            }
        )
    max_controls = max((len(row) for row in control_rows), default=0)
    controls_matrix = np.full((len(control_rows), max_controls), np.nan, dtype=float)
    for i, row in enumerate(control_rows):
        controls_matrix[i, : len(row)] = row

    inference = contract["inference"]
    result = matched_sympatry_colour_assembly_test(
        focal_species=np.asarray(set_focal, dtype=int),
        sympatric_colour_distance=np.asarray(set_sym, dtype=float),
        allopatric_control_colour_distance=controls_matrix,
        minimum_controls_per_set=int(inference["minimum_controls_per_set"]),
        minimum_sets_per_focal=int(inference["minimum_sets_per_focal"]),
        minimum_focal_species=int(inference["minimum_focal_species"]),
        permutations=int(inference["permutations"]),
        seed=int(inference["seed"]),
        species_fdr_alpha=float(inference["species_fdr_alpha"]),
    )
    global_result = dataclasses.asdict(result["global"])
    species_result: dict[str, dict[str, object]] = {}
    inverse_focal = {value: key for key, value in focal_code.items()}
    for code, payload in result.get("species", {}).items():
        species_result[inverse_focal[int(code)]] = payload

    status = str(global_result["status"])
    p_two = float(global_result["p_two_sided"]) if math.isfinite(float(global_result["p_two_sided"])) else None
    mean_delta = float(global_result["mean_focal_delta"]) if math.isfinite(float(global_result["mean_focal_delta"])) else None
    supported = bool(status == "evaluated" and p_two is not None and p_two < float(inference["global_alpha"]))
    if supported and mean_delta is not None and mean_delta < 0:
        direction = "convergence_in_sympatry"
    elif supported and mean_delta is not None and mean_delta > 0:
        direction = "divergence_in_sympatry"
    else:
        direction = "no_supported_direction"

    details_path = args.output_dir / "global_rgfca_g5_matched_set_colour_distances_v1.csv"
    pd.DataFrame(set_details).to_csv(details_path, index=False)
    species_path = args.output_dir / "global_rgfca_g5_species_colour_assembly_v1.csv"
    species_rows = [{"focal_species": species, **payload} for species, payload in sorted(species_result.items())]
    pd.DataFrame(species_rows).to_csv(species_path, index=False)
    result_path = args.output_dir / "global_rgfca_g5_sympatry_colour_assembly_result_v1.json"
    payload = {
        "protocol": parent["protocol"],
        "status": "complete_frozen_g5_primary" if status == "evaluated" else status,
        "parent_contract_status": parent["status"],
        "matching_manifest_status": matching_manifest.get("status"),
        "primary": global_result,
        "supported_two_sided_alpha_0_05": supported,
        "supported_direction": direction,
        "sign_convention": result["sign_convention"],
        "primary_inference": result.get("primary_inference"),
        "species_fdr_alpha": result.get("species_fdr_alpha"),
        "species_characterization": species_result,
        "colour_schedule": {
            "appearance_ranks_used_per_species": fixed_ranks,
            "photo_schedule_audit": schedule_audit(schedule),
            "pair_colour_distance": "mean of 135 appearance-rank Jensen-Shannon distances between 20-photo species mean soft palettes"
        },
        "coverage": {
            "frozen_matched_sets": int(len(set_ids)),
            "unique_focal_species_in_matching": int(len(focal_labels)),
            "unique_pair_distances_computed": int(len(pair_cache))
        },
        "lineage": {
            "implementation_contract_sha256": sha256_file(CONTRACT),
            "parent_contract_sha256": sha256_file(PARENT),
            "matching_manifest_sha256": sha256_file(MATCH_MANIFEST),
            "matched_sets_sha256": sha256_file(MATCHED),
            "measured_table_sha256": sha256_file(MEASURED),
            "matched_set_colour_distances_sha256": sha256_file(details_path),
            "species_colour_assembly_sha256": sha256_file(species_path)
        },
        "claim_ceiling": parent["claim_ceiling"],
        "cannot_rescue_or_reclassify": ["G1", "species_disjoint_commonness", "RESOLVE", "WorldClim", "expanded_environmental_panel", "environmental_interaction_family"]
    }
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "n_focal_species": global_result["n_focal_species"],
        "n_matched_sets": global_result["n_matched_sets"],
        "mean_focal_delta": mean_delta,
        "p_two_sided": p_two,
        "supported": supported,
        "direction": direction,
    }, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
