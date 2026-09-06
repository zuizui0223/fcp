#!/usr/bin/env python3
"""Freeze outcome-blind G5 GBIF occupancy metadata and matched sets.

This script MUST run before any G5 flower-colour outcome is opened.  It consumes
only species eligibility and metadata identifiers from the global colour project,
plus independent GBIF/CHELSA/realm metadata.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import re
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests
import rasterio

from fcp_pipeline.global_g5_matching import (
    N_CELLS,
    build_matched_sets,
    cell_center,
    equal_area_cell_id,
)

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_g5_metadata_matching_implementation_contract_v1.json"
PARENT = ROOT / "docs/supporting/global_sympatry_colour_assembly_contract_v1.json"
SUPPORT = ROOT / "data/derived/global_monte_carlo_measurement_species_support_v1.csv"
CANDIDATES = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
REALM_GRID = ROOT / "data/atlas/environment/climate_ecoregion_grid_100km.csv"
CHELSA_VERIFY = ROOT / "docs/supporting/global_rgfca_chelsa_expanded_source_verification_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def request_json(url: str, *, params: dict[str, Any] | None = None) -> dict[str, Any]:
    waits = (2, 4, 8)
    last: Exception | None = None
    headers = {"User-Agent": "fcp-rgfca-g5-metadata-freeze/1.0 (github.com/zuizui0223/fcp)"}
    for attempt in range(4):
        try:
            response = requests.get(url, params=params, headers=headers, timeout=(30, 120))
            response.raise_for_status()
            payload = response.json()
            if not isinstance(payload, dict):
                raise RuntimeError("GBIF response is not a JSON object")
            return payload
        except Exception as exc:  # fixed transport recovery only
            last = exc
            if attempt >= 3:
                break
            time.sleep(waits[attempt])
    raise RuntimeError(f"request failed after fixed retries: {url}: {last}") from last


def resolve_taxon(species: str, api_base: str) -> dict[str, Any]:
    payload = request_json(f"{api_base}/species/match", params={"name": species})
    usage = payload.get("usageKey") or payload.get("speciesKey")
    return {
        "species": species,
        "gbif_usage_key": None if usage is None else int(usage),
        "gbif_match_type": str(payload.get("matchType") or ""),
        "gbif_confidence": payload.get("confidence"),
        "accepted_scientific_name": str(payload.get("scientificName") or payload.get("canonicalName") or ""),
        "accepted_genus": str(payload.get("genus") or ""),
        "accepted_family": str(payload.get("family") or ""),
        "gbif_status": str(payload.get("status") or ""),
    }


def _inat_identifier_match(record: dict[str, Any], excluded_observation_ids: set[int]) -> bool:
    for key in ("occurrenceID", "references"):
        text = str(record.get(key) or "")
        if "inaturalist" not in text.casefold():
            continue
        ids = {int(x) for x in re.findall(r"(?<!\d)(\d{3,})(?!\d)", text)}
        if ids & excluded_observation_ids:
            return True
    return False


def fetch_species_occurrences(
    taxon: dict[str, Any],
    *,
    api_base: str,
    cap: int,
    page_limit: int,
    excluded_observation_ids: set[int],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    species = str(taxon["species"])
    taxon_key = taxon.get("gbif_usage_key")
    audit: dict[str, Any] = {
        "species": species,
        "gbif_usage_key": taxon_key,
        "status": "unresolved_taxon",
        "reported_total": None,
        "records_after_filters": 0,
        "excluded_inaturalist_colour_records": 0,
        "unmatchable_inaturalist_identifier_rows": 0,
    }
    if taxon_key is None or not taxon.get("accepted_genus") or not taxon.get("accepted_family"):
        return audit, []
    base_params = {
        "taxon_key": int(taxon_key),
        "has_coordinate": "true",
        "occurrence_status": "PRESENT",
    }
    probe = request_json(f"{api_base}/occurrence/search", params={**base_params, "limit": 1, "offset": 0})
    total = int(probe.get("count") or 0)
    audit["reported_total"] = total
    if total > int(cap):
        audit["status"] = "not_evaluable_reported_total_above_frozen_cap"
        return audit, []
    rows: dict[int, dict[str, Any]] = {}
    for offset in range(0, total, int(page_limit)):
        payload = request_json(
            f"{api_base}/occurrence/search",
            params={**base_params, "limit": int(page_limit), "offset": int(offset)},
        )
        results = payload.get("results") or []
        if not isinstance(results, list):
            raise RuntimeError(f"GBIF results malformed for {species}")
        for record in results:
            if not isinstance(record, dict):
                continue
            key = record.get("key")
            lat = record.get("decimalLatitude")
            lon = record.get("decimalLongitude")
            try:
                key_i = int(key)
                lat_f = float(lat)
                lon_f = float(lon)
            except (TypeError, ValueError):
                continue
            if not (math.isfinite(lat_f) and math.isfinite(lon_f)):
                continue
            if lat_f < -90 or lat_f > 90 or lon_f < -180 or lon_f > 180:
                continue
            if _inat_identifier_match(record, excluded_observation_ids):
                audit["excluded_inaturalist_colour_records"] += 1
                continue
            text = f"{record.get('occurrenceID') or ''} {record.get('references') or ''}"
            if "inaturalist" in text.casefold() and not re.search(r"(?<!\d)\d{3,}(?!\d)", text):
                audit["unmatchable_inaturalist_identifier_rows"] += 1
            rows[key_i] = {
                "gbif_key": key_i,
                "species": species,
                "latitude": lat_f,
                "longitude": lon_f,
                "cell_id": equal_area_cell_id(lat_f, lon_f),
            }
    records = [rows[key] for key in sorted(rows)]
    audit["records_after_filters"] = len(records)
    audit["status"] = "complete_snapshot"
    return audit, records


def normalized_column_lookup(columns: list[str], priorities: tuple[str, ...]) -> str | None:
    norm = {re.sub(r"[^a-z0-9]+", "_", str(col).strip().casefold()).strip("_"): str(col) for col in columns}
    for candidate in priorities:
        if candidate in norm:
            return norm[candidate]
    return None


def load_realm_map(path: Path) -> tuple[dict[int, str], dict[str, str]]:
    grid = pd.read_csv(path)
    realm_col = normalized_column_lookup(list(grid.columns), ("realm", "realm_name", "biogeographic_realm"))
    if realm_col is None:
        raise RuntimeError(f"realm column not found in frozen 100-km grid: {list(grid.columns)}")
    cell_col = normalized_column_lookup(list(grid.columns), ("cell_id",))
    if cell_col is not None:
        cells = pd.to_numeric(grid[cell_col], errors="coerce")
    else:
        lat_col = normalized_column_lookup(list(grid.columns), ("latitude", "lat", "lat_center", "latitude_center"))
        lon_col = normalized_column_lookup(list(grid.columns), ("longitude", "lon", "lon_center", "longitude_center"))
        if lat_col is None or lon_col is None:
            raise RuntimeError(f"cell_id or latitude/longitude columns not found: {list(grid.columns)}")
        cells = pd.Series(
            [
                equal_area_cell_id(float(lat), float(lon)) if pd.notna(lat) and pd.notna(lon) else np.nan
                for lat, lon in zip(grid[lat_col], grid[lon_col])
            ]
        )
    by_cell: dict[int, list[str]] = defaultdict(list)
    for cell, realm in zip(cells, grid[realm_col]):
        if pd.isna(cell) or pd.isna(realm):
            continue
        value = str(realm).strip()
        if not value or value.casefold() in {"nan", "none", "na"}:
            continue
        cell_i = int(cell)
        if 0 <= cell_i < N_CELLS:
            by_cell[cell_i].append(value)
    realm_map: dict[int, str] = {}
    for cell, values in by_cell.items():
        counts = Counter(values)
        maximum = max(counts.values())
        realm_map[cell] = sorted([name for name, count in counts.items() if count == maximum])[0]
    return realm_map, {"realm_column": realm_col, "cell_column": cell_col or "reconstructed_from_lat_lon"}


def download_verified_chelsa(ids: tuple[str, ...], outdir: Path) -> dict[str, Path]:
    verification = json.loads(CHELSA_VERIFY.read_text(encoding="utf-8"))
    files = {str(item["id"]): item for item in verification.get("files", [])}
    selected: dict[str, Path] = {}
    for variable in ids:
        item = files.get(variable)
        if item is None:
            raise RuntimeError(f"verified CHELSA source absent for {variable}")
        path = outdir / str(item["filename"])
        h = hashlib.sha256()
        with requests.get(str(item["url"]), stream=True, timeout=(30, 600)) as response:
            response.raise_for_status()
            with path.open("wb") as fh:
                for chunk in response.iter_content(chunk_size=1 << 20):
                    if not chunk:
                        continue
                    fh.write(chunk)
                    h.update(chunk)
        if h.hexdigest() != str(item["sha256"]):
            raise RuntimeError(f"CHELSA SHA256 mismatch for {variable}")
        selected[variable] = path
    return selected


def sample_chelsa_by_cell(rasters: dict[str, Path], cell_ids: list[int]) -> dict[str, dict[int, float]]:
    coords = []
    for cell in cell_ids:
        lat, lon = cell_center(cell)
        coords.append((lon, lat))
    out: dict[str, dict[int, float]] = {}
    for variable, path in rasters.items():
        with rasterio.open(path) as ds:
            values: dict[int, float] = {}
            nodata = ds.nodata
            for cell, sample in zip(cell_ids, ds.sample(coords)):
                value = float(sample[0])
                if not math.isfinite(value) or (nodata is not None and value == float(nodata)):
                    continue
                values[int(cell)] = value
            out[variable] = values
    return out


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    finite = np.isfinite(values.to_numpy())
    out = pd.Series(np.nan, index=series.index, dtype=float)
    if np.count_nonzero(finite) < 2:
        return out
    mean = float(values[finite].mean())
    std = float(values[finite].std(ddof=0))
    if not math.isfinite(std) or std <= 0:
        return out
    out.loc[finite] = (values.loc[finite] - mean) / std
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    if contract.get("status") != "prospective_technical_implementation_frozen_before_any_new_g5_gbif_payload_or_any_g5_colour_outcome_is_opened":
        raise RuntimeError("G5 implementation contract drift")
    if parent.get("status") != "prospective_design_frozen_during_metadata_only_global_species_discovery_before_any_global_monte_carlo_colour_field":
        raise RuntimeError("G5 parent contract drift")

    support = pd.read_csv(SUPPORT)
    eligible = support.loc[support["global_measurement_evaluable"].astype(str).str.casefold() == "true", ["species", "inat_taxon_id"]].copy()
    eligible["species"] = eligible["species"].astype(str)
    eligible = eligible.sort_values("species", kind="mergesort").drop_duplicates("species", keep="first").reset_index(drop=True)
    expected = int(contract["target_universe"]["expected_eligible_species"])
    if len(eligible) != expected:
        raise RuntimeError(f"eligible species count drift: {len(eligible)} != {expected}")

    excluded_ids = set(
        pd.read_csv(CANDIDATES, usecols=["observation_id"])["observation_id"].dropna().astype(int).tolist()
    )
    api_base = str(contract["gbif_transport"]["api"]).rstrip("/")
    cap = int(contract["gbif_transport"]["maximum_total_records_for_full_snapshot"])
    page_limit = int(contract["gbif_transport"]["occurrence_parameters"]["limit"])

    with cf.ThreadPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        taxon_rows = list(pool.map(lambda name: resolve_taxon(name, api_base), eligible["species"].tolist()))
    taxon_frame = pd.DataFrame(taxon_rows).sort_values("species", kind="mergesort").reset_index(drop=True)

    audits: list[dict[str, Any]] = []
    record_rows: list[dict[str, Any]] = []
    def one(row: dict[str, Any]):
        return fetch_species_occurrences(
            row,
            api_base=api_base,
            cap=cap,
            page_limit=page_limit,
            excluded_observation_ids=excluded_ids,
        )
    with cf.ThreadPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        futures = [pool.submit(one, row) for row in taxon_frame.to_dict(orient="records")]
        for index, future in enumerate(futures, start=1):
            audit, rows = future.result()
            audits.append(audit)
            record_rows.extend(rows)
            if index % 25 == 0 or index == len(futures):
                print(json.dumps({"gbif_species_complete": index, "total_species": len(futures)}), flush=True)

    audit_frame = pd.DataFrame(audits).sort_values("species", kind="mergesort").reset_index(drop=True)
    records = pd.DataFrame(record_rows, columns=["gbif_key", "species", "latitude", "longitude", "cell_id"])
    if len(records):
        records = records.sort_values(["species", "gbif_key"], kind="mergesort").drop_duplicates("gbif_key", keep="first")
    occupancy = (
        records.groupby(["species", "cell_id"], as_index=False).size().rename(columns={"size": "record_count"})
        if len(records)
        else pd.DataFrame(columns=["species", "cell_id", "record_count"])
    )
    occupancy = occupancy.sort_values(["species", "cell_id"], kind="mergesort").reset_index(drop=True)
    occupancy_sets = {
        species: set(group["cell_id"].astype(int).tolist())
        for species, group in occupancy.groupby("species", sort=True)
    }

    meta = eligible.merge(taxon_frame, on="species", how="left", validate="one_to_one")
    meta = meta.merge(audit_frame, on=["species", "gbif_usage_key"], how="left", validate="one_to_one")
    occupied_counts = occupancy.groupby("species")["cell_id"].nunique().rename("occupied_cells")
    meta = meta.merge(occupied_counts, on="species", how="left")
    meta["occupied_cells"] = meta["occupied_cells"].fillna(0).astype(int)
    meta["records_after_filters"] = meta["records_after_filters"].fillna(0).astype(int)
    meta["occurrence_eligible"] = (
        (meta["status"] == "complete_snapshot")
        & (meta["records_after_filters"] >= int(parent["sympatry_frame"]["species_occurrence_gate"]["minimum_georeferenced_records"]))
        & (meta["occupied_cells"] >= int(parent["sympatry_frame"]["species_occurrence_gate"]["minimum_occupied_cells"]))
    )

    realm_map, realm_columns = load_realm_map(REALM_GRID)
    dominant_realms: list[str | None] = []
    for row in meta.itertuples(index=False):
        cells = occupancy_sets.get(str(row.species), set())
        realms = [realm_map[cell] for cell in cells if cell in realm_map]
        if not realms:
            dominant_realms.append(None)
            continue
        counts = Counter(realms)
        maximum = max(counts.values())
        dominant_realms.append(sorted([name for name, count in counts.items() if count == maximum])[0])
    meta["dominant_realm"] = dominant_realms

    union_cells = sorted({cell for species, cells in occupancy_sets.items() if bool(meta.loc[meta["species"] == species, "occurrence_eligible"].iloc[0]) for cell in cells})
    chelsa_dir = outdir / "chelsa"
    chelsa_dir.mkdir(parents=True, exist_ok=True)
    chelsa_ids = ("bio01", "bio04", "bio12", "bio15")
    raster_paths = download_verified_chelsa(chelsa_ids, chelsa_dir)
    cell_climate = sample_chelsa_by_cell(raster_paths, union_cells)
    for variable in chelsa_ids:
        values: list[float] = []
        lookup = cell_climate[variable]
        for row in meta.itertuples(index=False):
            cells = occupancy_sets.get(str(row.species), set())
            finite = [lookup[cell] for cell in cells if cell in lookup]
            values.append(float(np.mean(finite)) if finite else np.nan)
        meta[variable] = values
        meta[f"z_{variable}"] = zscore(meta[variable].where(meta["occurrence_eligible"]))

    meta["log1p_occupied_cells"] = np.log1p(meta["occupied_cells"].astype(float))
    meta["log1p_gbif_occurrences"] = np.log1p(meta["records_after_filters"].astype(float))
    meta["z_log1p_occupied_cells"] = zscore(meta["log1p_occupied_cells"].where(meta["occurrence_eligible"]))
    meta["z_log1p_gbif_occurrences"] = zscore(meta["log1p_gbif_occurrences"].where(meta["occurrence_eligible"]))
    finite_cols = [f"z_{v}" for v in chelsa_ids] + ["z_log1p_occupied_cells", "z_log1p_gbif_occurrences"]
    finite_matrix = meta[finite_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    meta["matching_eligible"] = (
        meta["occurrence_eligible"].astype(bool)
        & meta["dominant_realm"].notna()
        & meta["accepted_genus"].astype(str).str.len().gt(0)
        & meta["accepted_family"].astype(str).str.len().gt(0)
        & np.isfinite(finite_matrix).all(axis=1)
    )

    matched = build_matched_sets(
        meta,
        occupancy_sets,
        maximum_distance=float(contract["matching"]["maximum_distance"]),
        controls_target=int(contract["matching"]["controls_target"]),
        controls_minimum=int(contract["matching"]["controls_minimum"]),
    )
    n_sets = int(matched["set_id"].nunique()) if len(matched) else 0
    set_focals = matched[["set_id", "focal_species"]].drop_duplicates() if len(matched) else pd.DataFrame(columns=["set_id", "focal_species"])
    focal_counts = set_focals.groupby("focal_species").size() if len(set_focals) else pd.Series(dtype=int)
    inference_ready_focals = int((focal_counts >= int(parent["primary_statistic"]["minimum_matched_sets_per_focal"])).sum())
    minimum_focals = int(parent["primary_statistic"]["minimum_evaluable_focal_species"])
    coverage_status = "pass_primary_inference_coverage" if inference_ready_focals >= minimum_focals else "not_evaluable_primary_inference_coverage"

    species_out = outdir / "global_rgfca_g5_species_metadata_v1.csv"
    occupancy_out = outdir / "global_rgfca_g5_occupancy_cells_v1.csv.gz"
    matched_out = outdir / "global_rgfca_g5_matched_sets_v1.csv"
    manifest_out = outdir / "global_rgfca_g5_matching_manifest_v1.json"
    meta.sort_values("species", kind="mergesort").to_csv(species_out, index=False)
    occupancy.to_csv(occupancy_out, index=False, compression={"method": "gzip", "mtime": 0})
    matched.to_csv(matched_out, index=False)

    manifest = {
        "protocol": contract["protocol"],
        "status": coverage_status,
        "target_species": int(len(meta)),
        "resolved_taxa": int(meta["gbif_usage_key"].notna().sum()),
        "complete_gbif_snapshots": int((meta["status"] == "complete_snapshot").sum()),
        "above_cap_species": int((meta["status"] == "not_evaluable_reported_total_above_frozen_cap").sum()),
        "occurrence_eligible_species": int(meta["occurrence_eligible"].sum()),
        "matching_eligible_species": int(meta["matching_eligible"].sum()),
        "occupancy_rows": int(len(occupancy)),
        "directed_matched_sets": n_sets,
        "focals_with_at_least_two_sets": inference_ready_focals,
        "minimum_evaluable_focal_species": minimum_focals,
        "colour_outcome_opened": False,
        "matching_used_colour": False,
        "realm_columns": realm_columns,
        "gbif": {
            "api_base": api_base,
            "record_cap": cap,
            "page_limit": page_limit,
            "total_records_after_filters": int(meta["records_after_filters"].sum()),
            "excluded_inaturalist_colour_records": int(meta["excluded_inaturalist_colour_records"].fillna(0).sum()),
        },
        "chelsa": {
            "variables": list(chelsa_ids),
            "verification_manifest_sha256": sha256_file(CHELSA_VERIFY),
            "unique_occupied_cells_sampled": int(len(union_cells)),
        },
        "source_sha256": {
            "implementation_contract": sha256_file(CONTRACT),
            "parent_contract": sha256_file(PARENT),
            "species_support": sha256_file(SUPPORT),
            "candidate_photo_metadata_source": sha256_file(CANDIDATES),
            "realm_grid": sha256_file(REALM_GRID),
        },
        "output_sha256": {
            "species_metadata": sha256_file(species_out),
            "occupancy_cells": sha256_file(occupancy_out),
            "matched_sets": sha256_file(matched_out),
        },
        "claim_ceiling": "Geographic sympatry and outcome-blind matched allopatry only; no G5 colour result exists in this manifest."
    }
    manifest_out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
