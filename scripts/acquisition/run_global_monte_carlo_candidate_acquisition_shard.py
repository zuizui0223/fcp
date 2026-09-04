#!/usr/bin/env python3
"""Run one frozen metadata-only candidate-acquisition shard for the global atlas."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.global_candidate_acquisition import (
    deterministic_candidate_pages,
    stable_candidate_query,
)
from fcp_pipeline.global_capacity_handoff import resolve_capacity_handoff
from fcp_pipeline.random_photo_h9_pool import (
    geographic_maximin,
    observer_cap,
    parse_h9_observation,
)
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_contract_v1.json"
AUTHORIZATION = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_authorization_v1.json"

PARSER_COLUMNS = [
    "species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
    "photo_license", "attribution", "latitude", "longitude", "positional_accuracy_m",
    "observed_on", "observer_id", "observer",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_exclusions(paths: Sequence[Path]) -> tuple[set[int], set[int], dict[str, str]]:
    observation_ids: set[int] = set()
    photo_ids: set[int] = set()
    hashes: dict[str, str] = {}
    for path in paths:
        if not path.exists():
            raise RuntimeError(f"missing prior-colour exclusion source: {path}")
        frame = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id"})
        if not {"observation_id", "photo_id"}.issubset(frame.columns):
            raise RuntimeError(f"bad prior-colour exclusion source: {path}")
        observation_ids.update(frame["observation_id"].dropna().astype(int).tolist())
        photo_ids.update(frame["photo_id"].dropna().astype(int).tolist())
        hashes[str(path.relative_to(ROOT))] = sha256_file(path)
    return observation_ids, photo_ids, hashes


def maximum_span_km(frame: pd.DataFrame) -> float:
    if len(frame) < 2:
        return 0.0
    lat = np.deg2rad(frame["latitude"].to_numpy(dtype=float))
    lon = np.deg2rad(frame["longitude"].to_numpy(dtype=float))
    c = np.cos(lat)
    xyz = np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return float(np.max(np.arccos(dot)) * 6371.0088)


def _results(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    raw = payload.get("results") or []
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
        raise RuntimeError("iNaturalist results is not a sequence")
    return [row for row in raw if isinstance(row, Mapping)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))

    if contract.get("status") != "frozen_after_v2_discovery_success_before_capacity_outcome_and_before_candidate_pixels":
        raise RuntimeError("candidate acquisition contract drifted")
    handoff = resolve_capacity_handoff(
        ROOT,
        AUTHORIZATION,
        minimum_species=int(contract["premeasurement_gate"]["minimum_full_target_species"]),
    )
    capacity = handoff.manifest
    target = int(handoff.target)
    selected = handoff.selected.copy()
    selected = selected.drop_duplicates("inat_taxon_id", keep="first").copy()
    selected["inat_taxon_id"] = selected["inat_taxon_id"].astype(int)
    selected["species"] = selected["species"].astype(str)
    selected = selected.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
    if len(selected) != int(capacity["selected_species"]):
        raise RuntimeError("capacity manifest/species-frame count mismatch")
    selected["global_row_index"] = np.arange(len(selected), dtype=int)

    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != int(contract["execution"]["deterministic_shards"]):
        raise RuntimeError("shard count differs from frozen candidate acquisition")
    if shard_index < 0 or shard_index >= shard_count:
        raise RuntimeError("invalid shard index")
    shard = selected.loc[(selected["global_row_index"] % shard_count) == shard_index].copy().reset_index(drop=True)
    if shard.empty:
        raise RuntimeError("candidate acquisition shard unexpectedly empty")

    exclusion_paths = [ROOT / value for value in parent["prior_experiment_exclusion_sources"]]
    excluded_obs, excluded_photo, exclusion_hashes = load_exclusions(exclusion_paths)
    q = contract["query"]
    selection = contract["selection"]
    allowed = frozenset(str(x).casefold() for x in q["allowed_photo_licenses"])
    client = InaturalistObservationClient(
        request_interval_seconds=float(q["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(q["request_retries"]),
        user_agent=f"fcp-global-candidate-acquisition-s{shard_index}/1.0 (github.com/zuizui0223/fcp)",
    )

    request_attempts = 0
    request_errors = 0
    audits: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []

    def request(taxon_id: int, page: int) -> tuple[Mapping[str, object] | None, str]:
        nonlocal request_attempts, request_errors
        request_attempts += 1
        params = stable_candidate_query(
            taxon_id,
            page=page,
            per_page=int(q["per_page"]),
            maximum_positional_accuracy_m=int(q["maximum_positional_accuracy_m"]),
            flowering_term_id=int(q["flowering_term_id"]),
            flowering_term_value_id=int(q["flowering_term_value_id"]),
            allowed_photo_licenses=tuple(q["allowed_photo_licenses"]),
        )
        try:
            return client.observations(params), ""
        except Exception as exc:
            request_errors += 1
            return None, f"{type(exc).__name__}:{str(exc)[:180]}"

    for processed, row in enumerate(shard.itertuples(index=False), start=1):
        taxon_id = int(row.inat_taxon_id)
        probe, probe_error = request(taxon_id, 1)
        total_results = 0
        pages: tuple[int, ...] = ()
        payload_by_page: dict[int, Mapping[str, object]] = {}
        page_errors: dict[int, str] = {}
        if probe is not None:
            try:
                total_results = int(probe.get("total_results") or len(_results(probe)))
                pages = deterministic_candidate_pages(
                    taxon_id,
                    total_results,
                    per_page=int(q["per_page"]),
                    maximum_api_page=int(q["maximum_api_page"]),
                    pages_per_species=int(q["candidate_pages_per_species"]),
                    seed=int(q["candidate_page_seed"]),
                )
                if 1 in pages:
                    payload_by_page[1] = probe
            except Exception as exc:
                probe_error = f"{type(exc).__name__}:{str(exc)[:180]}"
        if probe_error:
            page_errors[1] = probe_error
        if probe is not None:
            for page in pages:
                if page == 1:
                    continue
                payload, error = request(taxon_id, page)
                if payload is not None:
                    payload_by_page[int(page)] = payload
                if error:
                    page_errors[int(page)] = error

        raw_rows: list[Mapping[str, object]] = []
        for page in pages:
            payload = payload_by_page.get(int(page))
            if payload is None:
                continue
            raw_rows.extend(_results(payload))

        parsed_rows: list[dict[str, object]] = []
        locally_eligible = 0
        prior_excluded = 0
        wrong_taxon = 0
        for observation in raw_rows:
            taxon = observation.get("taxon") or {}
            if isinstance(taxon, Mapping) and int(taxon.get("id") or -1) != taxon_id:
                wrong_taxon += 1
                continue
            parsed = parse_h9_observation(
                observation,
                expected_taxon_id=taxon_id,
                maximum_positional_accuracy_m=float(q["maximum_positional_accuracy_m"]),
                allowed_photo_licenses=allowed,
            )
            if parsed is None:
                continue
            locally_eligible += 1
            if int(parsed["observation_id"]) in excluded_obs or int(parsed["photo_id"]) in excluded_photo:
                prior_excluded += 1
                continue
            parsed_rows.append(parsed)

        fresh = pd.DataFrame(parsed_rows, columns=PARSER_COLUMNS)
        if len(fresh):
            fresh = fresh.drop_duplicates("observation_id", keep="first")
            fresh = fresh.drop_duplicates("photo_id", keep="first").reset_index(drop=True)
        deduplicated = int(len(fresh))
        capped = observer_cap(fresh, int(selection["observer_cap_per_species"])) if len(fresh) else fresh
        capped_n = int(len(capped))
        full_target = capped_n >= target
        chosen = geographic_maximin(capped, target) if full_target else capped.iloc[0:0].copy()
        if full_target and len(chosen) != target:
            raise RuntimeError("geographic maximin failed to return the exact target")
        if full_target:
            chosen = chosen.copy()
            chosen["global_row_index"] = int(row.global_row_index)
            chosen["shard_index"] = shard_index
            for item in chosen.to_dict(orient="records"):
                candidate_rows.append(item)

        audits.append({
            "global_row_index": int(row.global_row_index),
            "shard_index": shard_index,
            "species": str(row.species),
            "inat_taxon_id": taxon_id,
            "selected_raw_photo_target": target,
            "total_results_probe": int(total_results),
            "candidate_pages": json.dumps(list(pages), separators=(",", ":")),
            "candidate_pages_requested": int(len(pages)),
            "candidate_page_errors": json.dumps(page_errors, sort_keys=True, separators=(",", ":")),
            "raw_results_from_candidate_pages": int(len(raw_rows)),
            "locally_eligible": int(locally_eligible),
            "prior_colour_experiment_excluded": int(prior_excluded),
            "wrong_taxon": int(wrong_taxon),
            "deduplicated_metadata_rows": deduplicated,
            "after_observer_cap": capped_n,
            "maximum_span_km_after_observer_cap": maximum_span_km(capped),
            "full_target": bool(full_target),
            "selected_candidate_rows": int(len(chosen)),
        })
        if processed % 50 == 0 or processed == len(shard):
            print(json.dumps({
                "shard": shard_index,
                "processed": processed,
                "shard_species": len(shard),
                "request_attempts": request_attempts,
                "request_errors": request_errors,
                "full_target_species": int(sum(bool(x["full_target"]) for x in audits)),
            }), flush=True)

    audit = pd.DataFrame(audits).sort_values("global_row_index", kind="mergesort").reset_index(drop=True)
    candidate_columns = PARSER_COLUMNS + ["row_hash", "global_row_index", "shard_index"]
    candidates = pd.DataFrame(candidate_rows)
    for column in candidate_columns:
        if column not in candidates.columns:
            candidates[column] = pd.Series(dtype="object")
    candidates = candidates[candidate_columns]
    if len(candidates):
        candidates = candidates.sort_values(["global_row_index", "row_hash"], kind="mergesort").reset_index(drop=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = args.output_dir / f"candidate_acquisition_shard_{shard_index:02d}_species.csv"
    candidate_path = args.output_dir / f"candidate_acquisition_shard_{shard_index:02d}_photos.csv"
    manifest_path = args.output_dir / f"candidate_acquisition_shard_{shard_index:02d}.json"
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    candidates.to_csv(candidate_path, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": "complete_metadata_only_candidate_acquisition_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "capacity_source": handoff.source,
        "global_selected_species": int(len(selected)),
        "shard_species": int(len(shard)),
        "selected_raw_photo_target": target,
        "request_attempts": int(request_attempts),
        "request_errors": int(request_errors),
        "full_target_species": int(audit["full_target"].sum()),
        "candidate_rows": int(len(candidates)),
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "lineage": {
            "candidate_contract_sha256": sha256_file(CONTRACT),
            "capacity_manifest_sha256": sha256_file(handoff.manifest_path),
            "selected_species_sha256": sha256_file(handoff.selected_path),
            "prior_colour_exclusion_sha256": exclusion_hashes,
            "species_audit_sha256": sha256_file(audit_path),
            "candidate_photos_sha256": sha256_file(candidate_path),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
