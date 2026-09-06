#!/usr/bin/env python3
"""Cache-resistant metadata-only repeated global species discovery V2."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.random_photo_pool import InaturalistObservationClient, inat_query_for_cell, parse_candidate_observation
from fcp_pipeline.shared_transition_surface import EqualAreaGrid

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_species_discovery_v2_contract.json"
V1_INDEX = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V1_SPECIES = ROOT / "data/frozen/global_monte_carlo_species_discovery_species_v1.csv"
V1_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_manifest_v1.json"

OBS_INDEX = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
SPECIES = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_species_v1.csv"
ROUND_AUDIT = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_round_audit_v1.csv"
CELL_AUDIT = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_cell_page_audit_v1.csv"
COMBINED_SPECIES = ROOT / "data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_v2_manifest_v1.json"
OUTPUTS = (OBS_INDEX, SPECIES, ROUND_AUDIT, CELL_AUDIT, COMBINED_SPECIES, MANIFEST)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def jaccard(a: set[str], b: set[str]) -> float:
    u = a | b
    return float(len(a & b) / len(u)) if u else float("nan")


def page_order(cell_id: int, max_page: int, seed: int) -> list[int]:
    if max_page <= 1:
        return [1]
    pages = list(range(2, max_page + 1))
    return sorted(
        pages,
        key=lambda p: hashlib.sha256(f"{seed}|{cell_id}|{p}".encode("utf-8")).hexdigest(),
    )


def compact(frame: pd.DataFrame, *, round_id: int, page: int) -> pd.DataFrame:
    cols = ["cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"]
    out = frame[cols].copy()
    out.insert(0, "page", int(page))
    out.insert(0, "round_id", int(round_id))
    return out


def main() -> None:
    existing = [str(p.relative_to(ROOT)) for p in OUTPUTS if p.exists()]
    if existing:
        raise RuntimeError("V2 discovery outputs already exist; refusing rerun: " + ", ".join(existing))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_after_v1_cache_diagnostic_before_any_v2_request_and_before_any_global_colour_pixel":
        raise RuntimeError("V2 discovery contract not frozen pre-query")
    v1_manifest = json.loads(V1_MANIFEST.read_text(encoding="utf-8"))
    if v1_manifest.get("status") != "complete_metadata_only_global_species_discovery":
        raise RuntimeError("V1 metadata lineage is incomplete")
    if v1_manifest.get("candidate_image_pixels_opened") is not False or v1_manifest.get("flower_colour_used") is not False:
        raise RuntimeError("V1 lineage opened colour/pixels")

    v1_index = pd.read_csv(V1_INDEX, usecols=["observation_id", "photo_id"])
    excluded_obs = set(v1_index["observation_id"].dropna().astype(int))
    excluded_photo = set(v1_index["photo_id"].dropna().astype(int))
    v1_species = pd.read_csv(V1_SPECIES)
    v1_species_set = set(v1_species["species"].astype(str))

    spec = contract["fresh_discovery"]
    grid = EqualAreaGrid(n_lon=int(spec["grid"]["n_lon"]), n_sinlat=int(spec["grid"]["n_sinlat"]))
    rounds = int(spec["rounds"])
    per_page = int(spec["per_page"])
    max_api_page = int(spec["maximum_api_page"])
    seed = int(spec["page_schedule_seed"])
    allowed = frozenset(str(x).casefold() for x in spec["allowed_photo_licenses"])
    client = InaturalistObservationClient(
        request_interval_seconds=float(spec["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(spec["request_retries"]),
        user_agent="fcp-global-monte-carlo-species-discovery-v2/1.0 (github.com/zuizui0223/fcp)",
    )

    cell_max_page: dict[int, int] = {}
    cell_orders: dict[int, list[int]] = {}
    cell_pages_used: dict[int, list[int]] = {i: [] for i in range(grid.n_cells)}
    round_sets: list[set[str]] = []
    parts: list[pd.DataFrame] = []
    round_rows: list[dict[str, object]] = []
    total_errors = 0
    total_attempts = 0
    total_returned = 0
    total_accepted = 0
    cumulative: set[str] = set(v1_species_set)
    previous: set[str] | None = None

    for round_id in range(1, rounds + 1):
        rows: list[dict[str, object]] = []
        round_errors = 0
        round_returned = 0
        for cell_id in range(grid.n_cells):
            if round_id == 1:
                page = 1
            else:
                order = cell_orders[cell_id]
                page = order[(round_id - 2) % len(order)]
            params = inat_query_for_cell(
                grid,
                cell_id,
                per_page=per_page,
                taxon_id=int(spec["taxon_id"]),
                flowering_term_id=int(spec["flowering_term_id"]),
                flowering_term_value_id=int(spec["flowering_term_value_id"]),
                maximum_positional_accuracy_m=int(spec["maximum_positional_accuracy_m"]),
                allowed_photo_licenses=tuple(spec["allowed_photo_licenses"]),
            )
            params["order_by"] = str(spec["stable_order_by"])
            params["order"] = str(spec["stable_order"])
            params["page"] = int(page)
            total_attempts += 1
            cell_pages_used[cell_id].append(int(page))
            try:
                payload = client.observations(params)
                raw = payload.get("results") or []
                if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
                    raise RuntimeError("iNaturalist results is not a sequence")
                if round_id == 1:
                    total_results = int(payload.get("total_results") or 0)
                    max_page = max(1, min(max_api_page, int(math.ceil(total_results / per_page))))
                    cell_max_page[cell_id] = max_page
                    cell_orders[cell_id] = page_order(cell_id, max_page, seed)
            except Exception:
                round_errors += 1
                total_errors += 1
                if round_id == 1:
                    cell_max_page[cell_id] = 1
                    cell_orders[cell_id] = [1]
                continue
            round_returned += len(raw)
            total_returned += len(raw)
            seen_obs: set[int] = set()
            seen_photo: set[int] = set()
            for obs in raw:
                if not isinstance(obs, Mapping):
                    continue
                parsed = parse_candidate_observation(
                    obs,
                    expected_cell_id=cell_id,
                    grid=grid,
                    maximum_positional_accuracy_m=float(spec["maximum_positional_accuracy_m"]),
                    allowed_photo_licenses=allowed,
                )
                if parsed is None:
                    continue
                oid = int(parsed["observation_id"])
                pid = int(parsed["photo_id"])
                if oid in excluded_obs or pid in excluded_photo or oid in seen_obs or pid in seen_photo:
                    continue
                seen_obs.add(oid)
                seen_photo.add(pid)
                rows.append(parsed)
        fresh = pd.DataFrame(rows)
        if len(fresh):
            # Page is cell-specific, so recover it from the fixed cell/round schedule.
            chunks = []
            for cell_id, group in fresh.groupby("cell_id", sort=False):
                chunks.append(compact(group, round_id=round_id, page=cell_pages_used[int(cell_id)][round_id - 1]))
            parts.append(pd.concat(chunks, ignore_index=True))
            fresh_species = set(fresh["species"].astype(str))
        else:
            fresh_species = set()
        new_species = fresh_species - cumulative
        cumulative |= fresh_species
        round_rows.append({
            "round_id": round_id,
            "request_attempts": grid.n_cells,
            "request_errors": round_errors,
            "api_returned": round_returned,
            "accepted_after_v1_id_exclusion": int(len(fresh)),
            "round_species": len(fresh_species),
            "new_species_beyond_v1_and_previous_v2": len(new_species),
            "cumulative_v1_plus_v2_species": len(cumulative),
            "jaccard_previous_v2_round": np.nan if previous is None else jaccard(previous, fresh_species),
        })
        round_sets.append(fresh_species)
        previous = fresh_species
        total_accepted += int(len(fresh))
        print(json.dumps(round_rows[-1], sort_keys=True), flush=True)

    expected_attempts = int(spec["total_request_attempts"])
    if total_attempts != expected_attempts:
        raise RuntimeError(f"V2 request count drift: {total_attempts} != {expected_attempts}")

    if parts:
        index = pd.concat(parts, ignore_index=True)
        index = index.sort_values(["round_id", "cell_id", "observation_id", "photo_id"], kind="mergesort").reset_index(drop=True)
    else:
        index = pd.DataFrame(columns=["round_id", "page", "cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"])
    unique = index.drop_duplicates(["observation_id", "photo_id"], keep="first")
    if len(unique):
        first = index.groupby(["species", "inat_taxon_id"], observed=True)["round_id"].min().rename("first_v2_round")
        seen = index.groupby(["species", "inat_taxon_id"], observed=True)["round_id"].nunique().rename("n_v2_rounds_seen")
        obs_n = unique.groupby(["species", "inat_taxon_id"], observed=True)["observation_id"].nunique().rename("unique_v2_observations")
        cells_n = unique.groupby(["species", "inat_taxon_id"], observed=True)["cell_id"].nunique().rename("occupied_v2_cells")
        species_v2 = pd.concat([first, seen, obs_n, cells_n], axis=1).reset_index()
    else:
        species_v2 = pd.DataFrame(columns=["species", "inat_taxon_id", "first_v2_round", "n_v2_rounds_seen", "unique_v2_observations", "occupied_v2_cells"])
    species_v2 = species_v2.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)

    v1_simple = v1_species[["species", "inat_taxon_id"]].copy()
    v1_simple["source_v1"] = True
    v2_simple = species_v2[["species", "inat_taxon_id"]].copy()
    v2_simple["source_v2"] = True
    combined = pd.merge(v1_simple, v2_simple, on="inat_taxon_id", how="outer", suffixes=("_v1", "_v2"))
    combined["species"] = combined["species_v2"].fillna(combined["species_v1"])
    combined["source_v1"] = combined["source_v1"].fillna(False).astype(bool)
    combined["source_v2"] = combined["source_v2"].fillna(False).astype(bool)
    combined = combined[["species", "inat_taxon_id", "source_v1", "source_v2"]].sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)

    odd = set().union(*(round_sets[i] for i in range(0, len(round_sets), 2))) if round_sets else set()
    even = set().union(*(round_sets[i] for i in range(1, len(round_sets), 2))) if len(round_sets) > 1 else set()
    v2_union = set(species_v2["species"].astype(str))
    cell_audit = pd.DataFrame([
        {
            "cell_id": cell_id,
            "available_pages_capped_50": cell_max_page[cell_id],
            "distinct_pages_used": len(set(cell_pages_used[cell_id])),
            "page_reuse_fraction": 1.0 - len(set(cell_pages_used[cell_id])) / float(rounds),
            "fewer_than_20_distinct_available_pages": bool(cell_max_page[cell_id] < 20),
        }
        for cell_id in range(grid.n_cells)
    ])
    error_fraction = total_errors / float(total_attempts)
    status = "complete_cache_resistant_metadata_only_global_species_discovery_v2" if error_fraction <= float(contract["failure_policy"]["maximum_error_fraction"]) else str(contract["failure_policy"]["if_exceeded"])

    OBS_INDEX.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    index.to_csv(OBS_INDEX, index=False, compression="gzip", lineterminator="\n")
    species_v2.to_csv(SPECIES, index=False, lineterminator="\n")
    pd.DataFrame(round_rows).to_csv(ROUND_AUDIT, index=False, lineterminator="\n")
    cell_audit.to_csv(CELL_AUDIT, index=False, lineterminator="\n")
    combined.to_csv(COMBINED_SPECIES, index=False, lineterminator="\n")

    manifest = {
        "protocol": contract["protocol"],
        "status": status,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "v1_independence_diagnostic_accepted": False,
        "v1_species_retained_as_metadata_discoveries": int(len(v1_species)),
        "v2": {
            "rounds": rounds,
            "request_attempts": total_attempts,
            "request_errors": total_errors,
            "request_error_fraction": error_fraction,
            "api_returned": total_returned,
            "accepted_after_v1_id_exclusion": total_accepted,
            "unique_observations": int(len(unique)),
            "unique_species": int(len(species_v2)),
            "new_species_beyond_v1": int((~combined["source_v1"] & combined["source_v2"]).sum()),
            "odd_even_species_jaccard": jaccard(odd, even),
            "v1_vs_v2_species_jaccard": jaccard(v1_species_set, v2_union),
            "cells_with_fewer_than_20_distinct_available_pages": int(cell_audit["fewer_than_20_distinct_available_pages"].sum()),
            "median_distinct_pages_used": float(cell_audit["distinct_pages_used"].median()),
        },
        "combined": {"species": int(len(combined))},
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "v1_manifest_sha256": sha256_file(V1_MANIFEST),
            "v1_observation_index_sha256": sha256_file(V1_INDEX),
            "v1_species_sha256": sha256_file(V1_SPECIES),
        },
        "files": {
            "v2_observation_index": {"path": str(OBS_INDEX.relative_to(ROOT)), "sha256": sha256_file(OBS_INDEX)},
            "v2_species": {"path": str(SPECIES.relative_to(ROOT)), "sha256": sha256_file(SPECIES)},
            "v2_round_audit": {"path": str(ROUND_AUDIT.relative_to(ROOT)), "sha256": sha256_file(ROUND_AUDIT)},
            "v2_cell_page_audit": {"path": str(CELL_AUDIT.relative_to(ROOT)), "sha256": sha256_file(CELL_AUDIT)},
            "combined_species": {"path": str(COMBINED_SPECIES.relative_to(ROOT)), "sha256": sha256_file(COMBINED_SPECIES)},
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
