#!/usr/bin/env python3
"""Run one shard of the frozen metadata-only high-depth RGFCA capacity pilot."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import h9_query_for_species, observer_cap, parse_h9_observation
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_sharedness_v2_high_depth_capacity_pilot_contract_v1.json"
PARENT_AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rank_key(taxon_id: int) -> tuple[str, int]:
    value = int(taxon_id)
    return hashlib.sha256(str(value).encode("ascii")).hexdigest(), value


def load_exclusions(paths: Sequence[Path]) -> tuple[set[int], set[int]]:
    obs: set[int] = set()
    photos: set[int] = set()
    for path in paths:
        frame = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id"})
        if not {"observation_id", "photo_id"}.issubset(frame.columns):
            raise RuntimeError(f"bad exclusion source: {path}")
        obs.update(frame["observation_id"].dropna().astype(int).tolist())
        photos.update(frame["photo_id"].dropna().astype(int).tolist())
    return obs, photos


def frozen_pilot_frame(contract: dict[str, object]) -> pd.DataFrame:
    frame = pd.read_csv(PARENT_AUDIT)
    pool = frame.loc[pd.to_numeric(frame["raw_results"], errors="raise").astype(int).eq(200), ["species", "inat_taxon_id"]].copy()
    if len(pool) != int(contract["parent"]["required_right_censored_species"]):
        raise RuntimeError("right-censored pool differs from frozen contract")
    pool["inat_taxon_id"] = pool["inat_taxon_id"].astype(int)
    pool["rank_hash"] = [rank_key(x)[0] for x in pool["inat_taxon_id"]]
    pool = pool.sort_values(["rank_hash", "inat_taxon_id"], kind="mergesort").head(int(contract["pilot_species_selection"]["n"])).reset_index(drop=True)
    pool["pilot_rank"] = range(len(pool))
    return pool


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--shard-index", type=int, required=True)
    p.add_argument("--shard-count", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract["status"] != "frozen_before_any_high_depth_capacity_request":
        raise RuntimeError("high-depth capacity contract drifted")
    ex = contract["execution"]
    if ex["candidate_image_pixels_may_open"] is not False or ex["flower_colour_may_open"] is not False:
        raise RuntimeError("contract unexpectedly opens biological outcomes")
    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != int(ex["deterministic_shards"]) or not (0 <= shard_index < shard_count):
        raise RuntimeError("shard topology differs from frozen contract")

    pilot = frozen_pilot_frame(contract)
    shard = pilot.loc[(pilot["pilot_rank"] % shard_count) == shard_index].copy().reset_index(drop=True)
    q = contract["query"]
    targets = [int(x) for x in contract["fixed_targets_after_observer_cap"]]
    max_target = max(targets)
    allowed = frozenset(str(x).casefold() for x in q["allowed_photo_licenses"])
    exclusion_paths = [ROOT / p for p in q["prior_experiment_exclusion_sources"]]
    excluded_obs, excluded_photo = load_exclusions(exclusion_paths)
    client = InaturalistObservationClient(
        request_interval_seconds=float(q["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(q["request_retries"]),
        user_agent=f"fcp-rgfca-sharedness-v2-capacity-{shard_index}/1.0 (github.com/zuizui0223/fcp)",
    )

    rows: list[dict[str, object]] = []
    any_terminal_errors = 0
    for pos, species_row in enumerate(shard.itertuples(index=False), start=1):
        taxon_id = int(species_row.inat_taxon_id)
        parsed: list[dict[str, object]] = []
        seen_obs: set[int] = set()
        seen_photo: set[int] = set()
        pages_requested = 0
        terminal_error = ""
        api_exhausted = False
        total_results_seen: int | None = None
        raw_rows_seen = 0
        for page in range(1, int(q["maximum_pages_per_species"]) + 1):
            params = h9_query_for_species(
                taxon_id,
                per_page=int(q["per_page"]),
                maximum_positional_accuracy_m=int(q["maximum_positional_accuracy_m"]),
                flowering_term_id=int(q["flowering_term_id"]),
                flowering_term_value_id=int(q["flowering_term_value_id"]),
                allowed_photo_licenses=tuple(q["allowed_photo_licenses"]),
            )
            params["order_by"] = str(q["stable_order_by"])
            params["order"] = str(q["stable_order"])
            params["page"] = page
            pages_requested += 1
            try:
                payload = client.observations(params)
            except Exception as exc:
                terminal_error = f"{type(exc).__name__}:{str(exc)[:180]}"
                any_terminal_errors += 1
                break
            raw = payload.get("results") or []
            if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
                terminal_error = "RuntimeError:iNaturalist results is not a sequence"
                any_terminal_errors += 1
                break
            try:
                total_results_seen = int(payload.get("total_results")) if payload.get("total_results") is not None else total_results_seen
            except (TypeError, ValueError):
                pass
            raw_rows_seen += len(raw)
            for observation in raw:
                if not isinstance(observation, Mapping):
                    continue
                candidate = parse_h9_observation(
                    observation,
                    expected_taxon_id=taxon_id,
                    maximum_positional_accuracy_m=float(q["maximum_positional_accuracy_m"]),
                    allowed_photo_licenses=allowed,
                )
                if candidate is None:
                    continue
                oid = int(candidate["observation_id"])
                pid = int(candidate["photo_id"])
                if oid in excluded_obs or pid in excluded_photo or oid in seen_obs or pid in seen_photo:
                    continue
                seen_obs.add(oid)
                seen_photo.add(pid)
                parsed.append(candidate)
            frame = pd.DataFrame(parsed)
            capped = observer_cap(frame, int(q["observer_cap"])) if len(frame) else frame
            if len(capped) >= int(q["stop_when_after_observer_cap_reaches"]):
                break
            if len(raw) < int(q["per_page"]):
                api_exhausted = True
                break
            if total_results_seen is not None and page * int(q["per_page"]) >= total_results_seen:
                api_exhausted = True
                break

        frame = pd.DataFrame(parsed)
        capped = observer_cap(frame, int(q["observer_cap"])) if len(frame) else frame
        capped_n = int(len(capped))
        page_cap_hit = (
            terminal_error == ""
            and not api_exhausted
            and pages_requested >= int(q["maximum_pages_per_species"])
            and capped_n < max_target
        )
        entry: dict[str, object] = {
            "pilot_rank": int(species_row.pilot_rank),
            "rank_hash": str(species_row.rank_hash),
            "species": str(species_row.species),
            "inat_taxon_id": taxon_id,
            "pages_requested": pages_requested,
            "total_results_seen": total_results_seen,
            "raw_rows_seen": int(raw_rows_seen),
            "unique_locally_eligible_after_exclusion": int(len(frame)),
            "after_observer_cap": capped_n,
            "api_exhausted": bool(api_exhausted),
            "page_cap_hit_before_400": bool(page_cap_hit),
            "terminal_request_error": terminal_error,
        }
        for target in targets:
            entry[f"confirmed_at_least_{target}"] = bool(capped_n >= target)
            entry[f"confirmed_below_{target}"] = bool(terminal_error == "" and api_exhausted and capped_n < target)
            entry[f"indeterminate_{target}"] = bool(terminal_error != "" or (page_cap_hit and capped_n < target))
        rows.append(entry)
        if pos % 10 == 0 or pos == len(shard):
            print(json.dumps({"shard": shard_index, "processed": pos, "shard_species": len(shard), "terminal_errors": any_terminal_errors}), flush=True)

    out = pd.DataFrame(rows).sort_values("pilot_rank", kind="mergesort").reset_index(drop=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"rgfca_sharedness_v2_capacity_pilot_shard_{shard_index:02d}.csv"
    json_path = args.output_dir / f"rgfca_sharedness_v2_capacity_pilot_shard_{shard_index:02d}.json"
    out.to_csv(csv_path, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": "complete_metadata_only_high_depth_capacity_pilot_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "pilot_species_total": int(len(pilot)),
        "shard_species": int(len(out)),
        "terminal_request_error_species": int(out["terminal_request_error"].fillna("").astype(str).str.len().gt(0).sum()),
        "confirmed_at_least": {str(t): int(out[f"confirmed_at_least_{t}"].sum()) for t in targets},
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "lineage": {"contract_sha256": sha256_file(CONTRACT), "parent_audit_sha256": sha256_file(PARENT_AUDIT), "csv_sha256": sha256_file(csv_path)},
    }
    json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(manifest, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
