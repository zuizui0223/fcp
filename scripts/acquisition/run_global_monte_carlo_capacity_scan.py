#!/usr/bin/env python3
"""Run the frozen metadata-only capacity scan over the discovered global species frame.

No candidate image pixel is opened here.  Each discovered species receives one
worldwide random metadata request.  The automatic 100/80/60 raw-photo target is
selected only from metadata support after prior-ID exclusion and observer capping.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.random_photo_h9_pool import (
    h9_query_for_species,
    observer_cap,
    parse_h9_observation,
)
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"
DISCOVERY = ROOT / "data/frozen/global_monte_carlo_species_discovery_species_v1.csv"
DISCOVERY_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_manifest_v1.json"
AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v1.csv"
SELECTED = ROOT / "data/frozen/global_monte_carlo_capacity_scan_selected_species_v1.csv"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v1.json"
OUTPUTS = (AUDIT, SELECTED, MANIFEST)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def row_id_digest(frame: pd.DataFrame) -> str:
    if frame.empty:
        return hashlib.sha256(b"").hexdigest()
    text = "\n".join(
        f"{int(o)}:{int(p)}"
        for o, p in sorted(zip(frame["observation_id"], frame["photo_id"]))
    ) + "\n"
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def maximum_span_km(frame: pd.DataFrame) -> float:
    if len(frame) < 2:
        return 0.0
    lat = np.deg2rad(frame["latitude"].to_numpy(dtype=float))
    lon = np.deg2rad(frame["longitude"].to_numpy(dtype=float))
    cos_lat = np.cos(lat)
    xyz = np.column_stack([cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return float(np.max(np.arccos(dot)) * 6371.0088)


def load_exclusions(paths: Sequence[Path]) -> tuple[set[int], set[int], dict[str, str]]:
    obs: set[int] = set()
    photo: set[int] = set()
    hashes: dict[str, str] = {}
    for path in paths:
        if not path.exists():
            raise RuntimeError(f"missing prior-experiment exclusion source: {path}")
        frame = pd.read_csv(path, usecols=lambda col: col in {"observation_id", "photo_id"})
        if "observation_id" not in frame.columns or "photo_id" not in frame.columns:
            raise RuntimeError(f"exclusion source lacks observation/photo IDs: {path}")
        obs.update(frame["observation_id"].dropna().astype(int).tolist())
        photo.update(frame["photo_id"].dropna().astype(int).tolist())
        hashes[str(path.relative_to(ROOT))] = sha256_file(path)
    return obs, photo, hashes


def main() -> None:
    existing = [str(path.relative_to(ROOT)) for path in OUTPUTS if path.exists()]
    if existing:
        raise RuntimeError("capacity-scan output already exists; refusing new random species queries: " + ", ".join(existing))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract["status"] != "frozen_before_global_species_discovery_outcome_and_before_capacity_queries":
        raise RuntimeError("capacity-scan contract is not in the frozen pre-query state")
    if contract["outcome_firewall"]["flower_colour_used"] is not False:
        raise RuntimeError("capacity-scan contract unexpectedly permits flower colour")
    discovery_manifest = json.loads(DISCOVERY_MANIFEST.read_text(encoding="utf-8"))
    if discovery_manifest.get("status") != "complete_metadata_only_global_species_discovery":
        raise RuntimeError("global species discovery did not pass its coverage gate; capacity scan is closed")
    if discovery_manifest.get("candidate_image_pixels_opened") is not False:
        raise RuntimeError("species-discovery lineage unexpectedly opened pixels")

    species = pd.read_csv(DISCOVERY)
    required_columns = {"species", "inat_taxon_id"}
    if not required_columns.issubset(species.columns):
        raise RuntimeError("global discovery species frame lacks species/taxon ID")
    species = species.drop_duplicates(["inat_taxon_id"], keep="first").copy()
    species["inat_taxon_id"] = species["inat_taxon_id"].astype(int)
    species = species.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)

    exclusion_paths = [ROOT / path for path in contract["prior_experiment_exclusion_sources"]]
    excluded_obs, excluded_photo, exclusion_hashes = load_exclusions(exclusion_paths)

    query = contract["query"]
    allowed = frozenset(str(x).casefold() for x in query["allowed_photo_licenses"])
    targets = [int(x) for x in contract["target_rule"]["candidate_raw_photos_per_species"]]
    if targets != sorted(targets, reverse=True):
        raise RuntimeError("capacity targets must be ordered largest to smallest")
    observer_cap_n = int(query["observer_cap"])
    client = InaturalistObservationClient(
        request_interval_seconds=float(query["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(query["request_retries"]),
        user_agent="fcp-global-monte-carlo-capacity-scan/1.0 (github.com/zuizui0223/fcp)",
    )

    audit_rows: list[dict[str, object]] = []
    errors = 0
    for index, row in enumerate(species.itertuples(index=False), start=1):
        taxon_id = int(row.inat_taxon_id)
        species_name = str(row.species)
        params = h9_query_for_species(
            taxon_id,
            per_page=int(query["per_page"]),
            maximum_positional_accuracy_m=int(query["maximum_positional_accuracy_m"]),
            flowering_term_id=int(query["flowering_term_id"]),
            flowering_term_value_id=int(query["flowering_term_value_id"]),
            allowed_photo_licenses=tuple(query["allowed_photo_licenses"]),
        )
        try:
            payload = client.observations(params)
            raw = payload.get("results") or []
            if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
                raise RuntimeError("iNaturalist results is not a sequence")
            request_error = ""
        except Exception as exc:
            raw = []
            request_error = f"{type(exc).__name__}:{str(exc)[:180]}"
            errors += 1

        parsed_rows: list[dict[str, object]] = []
        locally_eligible = 0
        prior_excluded = 0
        wrong_taxon = 0
        for observation in raw:
            if not isinstance(observation, Mapping):
                continue
            taxon = observation.get("taxon") or {}
            if isinstance(taxon, Mapping) and int(taxon.get("id") or -1) != taxon_id:
                wrong_taxon += 1
                continue
            parsed = parse_h9_observation(
                observation,
                expected_taxon_id=taxon_id,
                maximum_positional_accuracy_m=float(query["maximum_positional_accuracy_m"]),
                allowed_photo_licenses=allowed,
            )
            if parsed is None:
                continue
            locally_eligible += 1
            oid = int(parsed["observation_id"])
            pid = int(parsed["photo_id"])
            if oid in excluded_obs or pid in excluded_photo:
                prior_excluded += 1
                continue
            parsed_rows.append(parsed)

        fresh = pd.DataFrame(parsed_rows)
        capped = observer_cap(fresh, observer_cap_n) if len(fresh) else fresh
        capped_n = int(len(capped))
        entry: dict[str, object] = {
            "species": species_name,
            "inat_taxon_id": taxon_id,
            "raw_results": int(len(raw)),
            "locally_eligible": int(locally_eligible),
            "prior_excluded": int(prior_excluded),
            "wrong_taxon": int(wrong_taxon),
            "after_observer_cap": capped_n,
            "maximum_span_km": maximum_span_km(capped),
            "capped_row_id_sha256": row_id_digest(capped),
            "request_error": request_error,
        }
        for target in targets:
            entry[f"eligible_raw_{target}"] = bool(capped_n >= target)
        audit_rows.append(entry)
        if index % 100 == 0 or index == len(species):
            print(json.dumps({"processed_species": index, "total_species": len(species), "request_errors": errors}), flush=True)

    audit = pd.DataFrame(audit_rows)
    target_counts = {str(target): int(audit[f"eligible_raw_{target}"].sum()) for target in targets}
    required_species = int(contract["target_rule"]["minimum_metadata_eligible_species"])
    selected_target = next((target for target in targets if target_counts[str(target)] >= required_species), None)
    error_fraction = float(errors / len(species)) if len(species) else 1.0
    error_ceiling = float(contract["request_failure_policy"]["maximum_error_fraction_for_capacity_decision"])
    request_coverage_ok = bool(error_fraction <= error_ceiling)

    if not request_coverage_ok:
        status = str(contract["request_failure_policy"]["if_exceeded"])
        selected_target = None
    elif selected_target is None:
        status = str(contract["target_rule"]["if_none_pass"])
    else:
        status = "complete_metadata_only_capacity_scan_target_selected"

    if selected_target is None:
        selected = audit.iloc[0:0][["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"]].copy()
    else:
        selected = audit.loc[
            audit[f"eligible_raw_{selected_target}"],
            ["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"],
        ].copy()
        selected = selected.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
        selected["selected_raw_photo_target"] = int(selected_target)

    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(AUDIT, index=False, lineterminator="\n")
    selected.to_csv(SELECTED, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": status,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "discovered_species_scanned": int(len(species)),
        "request_errors": int(errors),
        "request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": error_ceiling,
        "request_coverage_ok": request_coverage_ok,
        "target_counts": target_counts,
        "minimum_metadata_eligible_species": required_species,
        "selected_raw_photo_target": (None if selected_target is None else int(selected_target)),
        "selected_species": int(len(selected)),
        "actual_image_acquisition_authorized": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "discovery_species_sha256": sha256_file(DISCOVERY),
            "discovery_manifest_sha256": sha256_file(DISCOVERY_MANIFEST),
            "prior_exclusion_source_sha256": exclusion_hashes,
            "unique_prior_observation_ids": int(len(excluded_obs)),
            "unique_prior_photo_ids": int(len(excluded_photo)),
        },
        "files": {
            "species_audit": {"path": str(AUDIT.relative_to(ROOT)), "sha256": sha256_file(AUDIT)},
            "selected_species_frame": {"path": str(SELECTED.relative_to(ROOT)), "sha256": sha256_file(SELECTED)},
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
