#!/usr/bin/env python3
"""Run one fixed CHELSA process block in the expanded environmental panel.

Each four-axis block is evaluated with the exact RGFCA edge schedule and same 999
species-conditioned colour nulls.  This script reports a raw block p-value only;
final support is withheld until Holm adjustment across all evaluable fixed blocks.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as warp_transform

from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_rgfca_engine import (
    _rank_edges,
    _rank_edges_matrix,
    _raw_jsd_from_source_rows,
    _raw_jsd_matrix_from_source_rows,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
)

ROOT = Path(__file__).resolve().parents[2]
PANEL_CONTRACT = ROOT / "docs/supporting/global_rgfca_expanded_environmental_process_panel_contract_v1.json"
PANEL_EXECUTION = ROOT / "docs/supporting/global_rgfca_expanded_environmental_panel_execution_v1.json"
CHELSA_SOURCE = ROOT / "docs/supporting/global_rgfca_chelsa_expanded_source_verification_v1.json"
WORLDCLIM_SCRIPT = ROOT / "scripts/analysis/run_global_rgfca_worldclim_edge_secondary.py"
ALLOWED_BLOCKS = ("thermal_regime", "water_balance", "atmospheric_energy_dryness")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_worldclim_module():
    name = "rgfca_worldclim_secondary_for_chelsa"
    spec = importlib.util.spec_from_file_location(name, WORLDCLIM_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import frozen RGFCA edge implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def stream_global_z_and_sample(path: Path, latitude: np.ndarray, longitude: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    with rasterio.open(path) as ds:
        if ds.count != 1:
            raise RuntimeError(f"expected one-band CHELSA raster: {path}")
        nodata = ds.nodata
        scale = float(ds.scales[0]) if ds.scales else 1.0
        offset = float(ds.offsets[0]) if ds.offsets else 0.0
        if not np.isfinite(scale) or scale == 0 or not np.isfinite(offset):
            raise RuntimeError(f"invalid scale/offset in {path}")

        count = 0
        total = 0.0
        total_sq = 0.0
        raw_min = float("inf")
        raw_max = float("-inf")
        for _, window in ds.block_windows(1):
            raw = ds.read(1, window=window).astype(np.float64, copy=False)
            keep = np.isfinite(raw)
            if nodata is not None:
                keep &= raw != float(nodata)
            if not np.any(keep):
                continue
            values = raw[keep] * scale + offset
            count += int(len(values))
            total += float(np.sum(values, dtype=np.float64))
            total_sq += float(np.sum(values * values, dtype=np.float64))
            raw_min = min(raw_min, float(np.min(values)))
            raw_max = max(raw_max, float(np.max(values)))
        if count < 2:
            raise RuntimeError(f"insufficient global finite support in {path}")
        mean = total / count
        variance = max(0.0, total_sq / count - mean * mean)
        sd = float(np.sqrt(variance))
        if not np.isfinite(sd) or sd <= 0:
            raise RuntimeError(f"cannot globally standardize {path}")

        xs = np.asarray(longitude, dtype=float)
        ys = np.asarray(latitude, dtype=float)
        if ds.crs is None:
            raise RuntimeError(f"CHELSA raster has no CRS: {path}")
        if str(ds.crs) != "EPSG:4326":
            x2, y2 = warp_transform("EPSG:4326", ds.crs, xs.tolist(), ys.tolist())
            coords = list(zip(x2, y2))
        else:
            coords = list(zip(xs.tolist(), ys.tolist()))
        raw_sample = np.asarray([v[0] for v in ds.sample(coords)], dtype=np.float64)
        valid = np.isfinite(raw_sample)
        if nodata is not None:
            valid &= raw_sample != float(nodata)
        sample = np.full(len(raw_sample), np.nan, dtype=np.float64)
        physical = raw_sample[valid] * scale + offset
        sample[valid] = (physical - mean) / sd

        audit = {
            "global_finite_cells": int(count),
            "global_mean_physical": float(mean),
            "global_sd_ddof0_physical": float(sd),
            "global_min_physical": raw_min,
            "global_max_physical": raw_max,
            "scale": scale,
            "offset": offset,
            "photo_complete_count": int(np.count_nonzero(np.isfinite(sample))),
            "photo_complete_fraction": float(np.mean(np.isfinite(sample))),
            "width": int(ds.width),
            "height": int(ds.height),
            "crs": str(ds.crs),
        }
        return sample, audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--block", choices=ALLOWED_BLOCKS, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-result", type=Path, required=True)
    parser.add_argument("--output-species", type=Path, required=True)
    parser.add_argument("--output-null", type=Path, required=True)
    parser.add_argument("--null-batch-size", type=int, default=64)
    args = parser.parse_args()

    block_name = str(args.block)
    panel = json.loads(PANEL_CONTRACT.read_text())
    order = json.loads(PANEL_EXECUTION.read_text())
    source = json.loads(CHELSA_SOURCE.read_text())
    if panel.get("status") != "postoutcome_secondary_panel_frozen_before_any_expanded_environmental_panel_result_or_new_source_payload_is_opened":
        raise RuntimeError("expanded environmental panel contract drift")
    if order.get("status") != "technical_execution_order_frozen_before_any_expanded_panel_block_colour_alignment_result":
        raise RuntimeError("expanded environmental execution-order drift")
    if source.get("status") != "pass_exact_expanded_chelsa_source_acquisition_before_any_expanded_chelsa_colour_alignment":
        raise RuntimeError("exact expanded CHELSA source verification missing")
    expected_order = order["execution_order"]
    if block_name not in expected_order:
        raise RuntimeError("CHELSA block not in frozen panel order")

    entries = [f for f in source["files"] if f["block"] == block_name]
    if len(entries) != 4:
        raise RuntimeError(f"expected four fixed CHELSA axes for {block_name}, got {len(entries)}")
    raster_paths: list[Path] = []
    for entry in entries:
        p = args.source_dir / entry["filename"]
        if not p.exists():
            raise RuntimeError(f"missing CHELSA source: {p}")
        if sha256_file(p) != entry["sha256"]:
            raise RuntimeError(f"CHELSA source SHA drift: {entry['id']}")
        raster_paths.append(p)

    wc = load_worldclim_module()
    frame, inference, _secondary, _primary = wc.load_frozen_pool()
    pool = canonical_colour_pool(frame)
    outer = inference["outer_schedule"]
    g1 = inference["g1_primary"]
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=int(inference["input_gate"]["minimum_classifiable_photos_per_species"]),
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
    )
    blocks = wc.build_edge_blocks(pool, schedule, k=int(g1["k"]))
    total_occurrences = int(sum(len(b.edge_nodes) for b in blocks))
    print(json.dumps({"stage":"edges_built","block":block_name,"edge_occurrences":total_occurrences}), flush=True)

    columns=[]; source_audit={}
    for entry, path in zip(entries, raster_paths):
        z, audit = stream_global_z_and_sample(path, pool.latitude, pool.longitude)
        columns.append(z); source_audit[entry["id"]] = {**audit, "sha256":entry["sha256"]}
        print(json.dumps({"stage":"axis_ready","block":block_name,"axis":entry["id"],"photo_complete":audit["photo_complete_count"]}), flush=True)
    photo_z = np.column_stack(columns)
    joint_photo_complete = np.isfinite(photo_z).all(axis=1)
    print(json.dumps({
        "stage":"chelsa_block_join_ready","block":block_name,
        "joint_photo_complete":int(np.count_nonzero(joint_photo_complete)),
        "photo_total":int(len(photo_z)),
        "joint_photo_complete_fraction":float(np.mean(joint_photo_complete)),
    }), flush=True)

    species_parts=[]; distance_parts=[]; external_parts=[]; cursor=0
    for b in blocks:
        left=photo_z[b.edge_nodes[:,0]]; right=photo_z[b.edge_nodes[:,1]]
        keep=np.isfinite(left).all(axis=1) & np.isfinite(right).all(axis=1)
        external=np.full(len(b.edge_nodes),np.nan,dtype=float)
        if np.any(keep):
            delta=left[keep]-right[keep]
            external[keep]=np.sqrt(np.mean(delta*delta,axis=1))
        b.keep_external=keep
        n_keep=int(np.count_nonzero(keep))
        species_parts.append(b.global_species_index[keep])
        distance_parts.append(b.edge_distance_km[keep])
        external_parts.append(external[keep])
        b.sorted_target_positions=np.arange(cursor,cursor+n_keep,dtype=np.int64)
        cursor += n_keep
    if cursor == 0:
        raise RuntimeError(f"no edges evaluable for {block_name}")

    species_unsorted=np.concatenate(species_parts)
    distance_unsorted=np.concatenate(distance_parts)
    external_unsorted=np.concatenate(external_parts)
    sort_order=np.argsort(species_unsorted,kind="stable")
    inverse=np.empty(len(sort_order),dtype=np.int64); inverse[sort_order]=np.arange(len(sort_order),dtype=np.int64)
    species_sorted=species_unsorted[sort_order]
    distance_sorted=distance_unsorted[sort_order]
    external_sorted=external_unsorted[sort_order]
    species_counts=np.bincount(species_sorted,minlength=len(pool.species_labels)).astype(np.int64)
    cursor=0
    for b in blocks:
        keep=b.keep_external; assert keep is not None
        n_keep=int(np.count_nonzero(keep))
        b.sorted_target_positions=inverse[np.arange(cursor,cursor+n_keep,dtype=np.int64)]
        cursor += n_keep

    cache=build_pairwise_jsd_cache(pool)
    identity=np.arange(len(pool.photo_ids),dtype=np.int64)
    observed_scores=np.empty(len(species_sorted),dtype=float)
    for b in blocks:
        keep=b.keep_external; target=b.sorted_target_positions
        assert keep is not None and target is not None
        raw=_raw_jsd_from_source_rows(cache,b.edge_nodes,identity)
        rank=_rank_edges(raw,b.edge_species_slices)
        observed_scores[target]=rank[keep]

    minimum_edges=int(panel["primary_edge_test"]["minimum_edges_per_species"])
    minimum_species=int(panel["primary_edge_test"]["minimum_evaluable_species"])
    observed_values, observed_species_ids=wc._species_partial_values(
        observed_scores,external_sorted,distance_sorted,species_counts,minimum_edges=minimum_edges
    )
    if len(observed_values) < minimum_species:
        raise RuntimeError(f"only {len(observed_values)} species evaluable for {block_name}")
    observed_mean=float(np.mean(observed_values))
    observed_median=float(np.median(observed_values))
    observed_positive=float(np.mean(observed_values>0))
    print(json.dumps({"stage":"observed_block","block":block_name,"n_species":int(len(observed_values)),"mean_partial_rho":observed_mean,"median_partial_rho":observed_median,"positive_fraction":observed_positive}), flush=True)

    batch_size=int(args.null_batch_size)
    if batch_size<1 or batch_size>256: raise ValueError("null batch size must lie in [1,256]")
    null_indices=np.arange(999,dtype=np.int64); null_stat=np.empty(999,dtype=float)
    for batch_start in range(0,999,batch_size):
        batch_indices=null_indices[batch_start:batch_start+batch_size]
        source_rows=null_source_row_matrix(pool,batch_indices,master_seed=int(g1["null_master_seed"]))
        batch_scores=np.empty((len(batch_indices),len(species_sorted)),dtype=float)
        for b in blocks:
            keep=b.keep_external; target=b.sorted_target_positions
            assert keep is not None and target is not None
            raw=_raw_jsd_matrix_from_source_rows(cache,b.edge_nodes,source_rows)
            rank=_rank_edges_matrix(raw,b.edge_species_slices)
            batch_scores[:,target]=rank[:,keep]
        for local_index,null_index in enumerate(batch_indices):
            values,ids=wc._species_partial_values(batch_scores[local_index],external_sorted,distance_sorted,species_counts,minimum_edges=minimum_edges)
            if not np.array_equal(ids,observed_species_ids):
                raise RuntimeError(f"null {int(null_index)} changed evaluable species identity")
            null_stat[int(null_index)]=float(np.mean(values))
        print(json.dumps({"stage":"null_batch_complete","block":block_name,"first_null":int(batch_indices[0]),"last_null":int(batch_indices[-1]),"completed":int(batch_indices[-1])+1,"total":999}),flush=True)

    p_upper=float((1+np.count_nonzero(null_stat>=observed_mean))/1000)
    alpha=float(panel["primary_edge_test"]["alpha"])
    raw_positive=bool(observed_mean>0 and p_upper<alpha)
    args.output_result.parent.mkdir(parents=True,exist_ok=True)
    args.output_species.parent.mkdir(parents=True,exist_ok=True)
    args.output_null.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame({
        "species":[pool.species_labels[int(i)] for i in observed_species_ids],
        "species_index":observed_species_ids,
        "retained_edge_occurrences":species_counts[observed_species_ids],
        f"partial_rho_colour_vs_{block_name}_given_distance":observed_values,
    }).to_csv(args.output_species,index=False,lineterminator="\n")
    pd.DataFrame({"null_index":null_indices,"global_mean_species_partial_rho":null_stat}).to_csv(args.output_null,index=False,lineterminator="\n")

    payload={
        "protocol":panel["protocol"],
        "status":f"complete_{block_name}_block_raw_pending_full_panel_holm",
        "block":block_name,
        "inferential_role":panel["inferential_role"],
        "final_panel_support_decision_available":False,
        "reason_final_support_pending":"Holm adjustment requires all evaluable fixed blocks in the nonselective five-block panel",
        "predictor_source":"exact SHA-verified CHELSA v2.1 1981-2010 four-axis process block",
        "axes":[e["id"] for e in entries],
        "observed_mean_species_partial_rho":observed_mean,
        "observed_median_species_partial_rho":observed_median,
        "observed_positive_species_fraction":observed_positive,
        "n_evaluable_species":int(len(observed_values)),
        "minimum_evaluable_species":minimum_species,
        "total_scheduled_edge_occurrences":total_occurrences,
        "block_evaluable_edge_occurrences":int(len(species_sorted)),
        "block_edge_occurrence_coverage_fraction":float(len(species_sorted)/total_occurrences),
        "joint_photo_complete_count":int(np.count_nonzero(joint_photo_complete)),
        "joint_photo_complete_fraction":float(np.mean(joint_photo_complete)),
        "source_audit":source_audit,
        "null_permutations":999,"null_indices_exact_0_998":True,
        "null_mean":float(np.mean(null_stat)),"null_q025":float(np.quantile(null_stat,.025)),"null_q975":float(np.quantile(null_stat,.975)),
        "raw_p_upper":p_upper,"raw_nominal_positive_before_panel_Holm":raw_positive,
        "alpha":alpha,"p_holm":None,"panel_supported":None,
        "individual_variable_decomposition_open":False,"causal_language_allowed":False,
        "execution_audit":{
            "all_200_outer_realizations_used":len(blocks)==200,
            "same_RGFCA_edge_geometry":True,"repeated_edge_occurrences_deduplicated":False,
            "missing_environment_recoded_to_zero":False,"four_fixed_axes_only":True,
            "same_999_RGFCA_colour_nulls":True,"block_result_cannot_stop_later_fixed_blocks":True,
        },
        "lineage":{
            "panel_contract_sha256":sha256_file(PANEL_CONTRACT),
            "panel_execution_sha256":sha256_file(PANEL_EXECUTION),
            "chelsa_source_verification_sha256":sha256_file(CHELSA_SOURCE),
            "measurement_manifest_sha256":sha256_file(wc.MEASUREMENT),
            "measured_table_sha256":sha256_file(wc.MEASURED),
            "primary_g1_result_sha256":sha256_file(wc.PRIMARY_G1),
            "species_result_sha256":sha256_file(args.output_species),
            "null_result_sha256":sha256_file(args.output_null),
        },
        "files":{"species":str(args.output_species),"null":str(args.output_null)},
    }
    args.output_result.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2),flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
