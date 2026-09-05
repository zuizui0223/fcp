#!/usr/bin/env python3
"""Retry one deterministic shard of transport-only HTTP 429 capacity failures."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.global_capacity_recovery import frozen_recovery_rows
from fcp_pipeline.random_photo_h9_pool import h9_query_for_species, observer_cap, parse_h9_observation
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_capacity_429_recovery_contract_v1.json"
PARENT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"
FAILED_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json"
FAILED_AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v2.csv"


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
    c = np.cos(lat)
    xyz = np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return float(np.max(np.arccos(dot)) * 6371.0088)


def load_exclusions(paths: Sequence[Path]) -> tuple[set[int], set[int]]:
    obs: set[int] = set()
    photo: set[int] = set()
    for path in paths:
        frame = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id"})
        if not {"observation_id", "photo_id"}.issubset(frame.columns):
            raise RuntimeError(f"bad exclusion source: {path}")
        obs.update(frame["observation_id"].dropna().astype(int).tolist())
        photo.update(frame["photo_id"].dropna().astype(int).tolist())
    return obs, photo


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--shard-index", type=int, required=True)
    p.add_argument("--shard-count", type=int, default=4)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    failed = json.loads(FAILED_MANIFEST.read_text(encoding="utf-8"))
    trigger = contract["trigger"]
    if failed.get("status") != trigger["required_parent_status"]:
        raise RuntimeError("capacity v2 is not the frozen transport-failure state")
    if float(failed.get("request_error_fraction") or 0.0) <= float(trigger["required_parent_request_error_fraction_above"]):
        raise RuntimeError("capacity v2 did not exceed the frozen request-error ceiling")
    if failed.get("candidate_image_pixels_opened") is not False or failed.get("flower_colour_used") is not False:
        raise RuntimeError("capacity v2 opened forbidden outcomes")

    audit = pd.read_csv(FAILED_AUDIT)
    retry = frozen_recovery_rows(audit)
    if retry.empty:
        raise RuntimeError("transport recovery has no frozen 429 rows")
    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != int(contract["execution"]["deterministic_shards"]):
        raise RuntimeError("recovery shard count drifted")
    retry["recovery_position"] = np.arange(len(retry), dtype=int)
    shard = retry.loc[(retry["recovery_position"] % shard_count) == shard_index].copy().reset_index(drop=True)

    q = contract["query"]
    targets = [int(x) for x in contract["target_rule"]["candidate_raw_photos_per_species"]]
    allowed = frozenset(str(x).casefold() for x in q["allowed_photo_licenses"])
    exclusion_paths = [ROOT / path for path in parent["prior_experiment_exclusion_sources"]]
    excluded_obs, excluded_photo = load_exclusions(exclusion_paths)
    client = InaturalistObservationClient(
        request_interval_seconds=float(q["request_interval_seconds"]),
        timeout_seconds=float(q["timeout_seconds"]),
        max_retries=int(q["request_retries"]),
        user_agent=f"fcp-global-capacity-429-recovery-s{shard_index}/1.0 (github.com/zuizui0223/fcp)",
    )

    rows: list[dict[str, object]] = []
    for processed, row in enumerate(shard.itertuples(index=False), start=1):
        taxon_id = int(row.inat_taxon_id)
        params = h9_query_for_species(
            taxon_id,
            per_page=int(q["per_page"]),
            maximum_positional_accuracy_m=int(q["maximum_positional_accuracy_m"]),
            flowering_term_id=int(q["flowering_term_id"]),
            flowering_term_value_id=int(q["flowering_term_value_id"]),
            allowed_photo_licenses=tuple(q["allowed_photo_licenses"]),
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

        parsed_rows: list[dict[str, object]] = []
        locally_eligible = prior_excluded = wrong_taxon = 0
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
        fresh = pd.DataFrame(parsed_rows)
        capped = observer_cap(fresh, int(q["observer_cap"])) if len(fresh) else fresh
        capped_n = int(len(capped))
        entry = {
            "global_row_index": int(row.global_row_index),
            "shard_index": int(row.shard_index),
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
        rows.append(entry)
        if processed % 100 == 0 or processed == len(shard):
            print(json.dumps({"recovery_shard": shard_index, "processed": processed, "rows": len(shard)}), flush=True)

    out = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"capacity_429_recovery_shard_{shard_index:02d}.csv"
    json_path = args.output_dir / f"capacity_429_recovery_shard_{shard_index:02d}.json"
    out.to_csv(csv_path, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": "complete_transport_429_recovery_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "frozen_retry_rows": int(len(retry)),
        "shard_rows": int(len(shard)),
        "remaining_errors": int(out["request_error"].fillna("").astype(str).str.len().gt(0).sum()) if len(out) else 0,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "failed_capacity_manifest_sha256": sha256_file(FAILED_MANIFEST),
            "failed_capacity_audit_sha256": sha256_file(FAILED_AUDIT),
            "recovery_csv_sha256": sha256_file(csv_path)
        }
    }
    json_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
