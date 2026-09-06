#!/usr/bin/env python3
"""Run one deterministic shard of the frozen metadata-only global capacity scan."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.random_photo_h9_pool import h9_query_for_species, observer_cap, parse_h9_observation
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_v2_execution_contract.json"
V2_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_v2_manifest_v1.json"
SPECIES_FRAME = ROOT / "data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def row_id_digest(frame: pd.DataFrame) -> str:
    if frame.empty:
        return hashlib.sha256(b"").hexdigest()
    text = "\n".join(f"{int(o)}:{int(p)}" for o, p in sorted(zip(frame["observation_id"], frame["photo_id"]))) + "\n"
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
            raise RuntimeError(f"missing exclusion source: {path}")
        frame = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id"})
        if not {"observation_id", "photo_id"}.issubset(frame.columns):
            raise RuntimeError(f"bad exclusion source: {path}")
        obs.update(frame["observation_id"].dropna().astype(int).tolist())
        photo.update(frame["photo_id"].dropna().astype(int).tolist())
        hashes[str(path.relative_to(ROOT))] = sha256_file(path)
    return obs, photo, hashes


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--shard-index", type=int, required=True)
    p.add_argument("--shard-count", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    execution = json.loads(EXECUTION.read_text(encoding="utf-8"))
    v2 = json.loads(V2_MANIFEST.read_text(encoding="utf-8"))
    if parent["status"] != "frozen_before_global_species_discovery_outcome_and_before_capacity_queries":
        raise RuntimeError("parent capacity contract drifted")
    if execution["status"] != "frozen_while_cache_resistant_v2_discovery_is_running_before_v2_outcome_and_before_capacity_queries":
        raise RuntimeError("capacity execution contract drifted")
    if execution["biological_rules_changed"] is not False:
        raise RuntimeError("capacity execution changed biological rules")
    if v2.get("status") != "complete_cache_resistant_metadata_only_global_species_discovery_v2":
        raise RuntimeError("cache-resistant V2 discovery did not complete successfully")
    if v2.get("candidate_image_pixels_opened") is not False or v2.get("flower_colour_used") is not False:
        raise RuntimeError("V2 discovery opened forbidden outcomes")

    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != int(execution["execution"]["deterministic_shards"]):
        raise RuntimeError("shard count differs from frozen execution")
    if shard_index < 0 or shard_index >= shard_count:
        raise RuntimeError("invalid shard index")

    species = pd.read_csv(SPECIES_FRAME)
    if not {"species", "inat_taxon_id"}.issubset(species.columns):
        raise RuntimeError("combined V2 species frame lacks species/taxon ID")
    species = species.drop_duplicates(["inat_taxon_id"], keep="first").copy()
    species["inat_taxon_id"] = species["inat_taxon_id"].astype(int)
    species["species"] = species["species"].astype(str)
    species = species.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
    species["global_row_index"] = np.arange(len(species), dtype=int)
    shard = species.loc[(species["global_row_index"] % shard_count) == shard_index].copy().reset_index(drop=True)
    if shard.empty:
        raise RuntimeError("capacity shard unexpectedly empty")

    exclusion_paths = [ROOT / path for path in parent["prior_experiment_exclusion_sources"]]
    excluded_obs, excluded_photo, exclusion_hashes = load_exclusions(exclusion_paths)
    query = parent["query"]
    targets = [int(x) for x in parent["target_rule"]["candidate_raw_photos_per_species"]]
    allowed = frozenset(str(x).casefold() for x in query["allowed_photo_licenses"])
    cap_n = int(query["observer_cap"])
    client = InaturalistObservationClient(
        request_interval_seconds=float(query["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(query["request_retries"]),
        user_agent=f"fcp-global-monte-carlo-capacity-shard-{shard_index}/1.0 (github.com/zuizui0223/fcp)",
    )

    audit_rows: list[dict[str, object]] = []
    errors = 0
    for processed, row in enumerate(shard.itertuples(index=False), start=1):
        taxon_id = int(row.inat_taxon_id)
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
        capped = observer_cap(fresh, cap_n) if len(fresh) else fresh
        capped_n = int(len(capped))
        entry: dict[str, object] = {
            "global_row_index": int(row.global_row_index),
            "shard_index": shard_index,
            "species": str(row.species),
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
        if processed % 100 == 0 or processed == len(shard):
            print(json.dumps({"shard": shard_index, "processed": processed, "shard_species": len(shard), "errors": errors}), flush=True)

    audit = pd.DataFrame(audit_rows).sort_values("global_row_index", kind="mergesort").reset_index(drop=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = args.output_dir / f"capacity_scan_shard_{shard_index:02d}.csv"
    manifest_path = args.output_dir / f"capacity_scan_shard_{shard_index:02d}.json"
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    manifest = {
        "protocol": execution["protocol"],
        "status": "complete_metadata_only_capacity_scan_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "global_species_count": int(len(species)),
        "shard_species_count": int(len(shard)),
        "request_attempts": int(len(shard)),
        "request_errors": int(errors),
        "target_counts": {str(t): int(audit[f"eligible_raw_{t}"].sum()) for t in targets},
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "lineage": {
            "parent_contract_sha256": sha256_file(PARENT),
            "execution_contract_sha256": sha256_file(EXECUTION),
            "v2_manifest_sha256": sha256_file(V2_MANIFEST),
            "combined_species_sha256": sha256_file(SPECIES_FRAME),
            "exclusion_source_sha256": exclusion_hashes,
            "audit_sha256": sha256_file(audit_path),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
