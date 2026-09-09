#!/usr/bin/env python3
"""Freeze one metadata-only shard of the RGFCA sharedness-v2 random geometry."""
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
CONTRACT = ROOT / "docs/supporting/rgfca_sharedness_v2_geometry_freeze_contract_v1.json"
CAPACITY_RESULT = ROOT / "docs/supporting/rgfca_sharedness_v2_high_depth_capacity_pilot_result_v1.json"
CAPACITY_ROWS = ROOT / "data/frozen/rgfca_sharedness_v2_high_depth_capacity_pilot_v1.csv"
SELECTION_RULE = ROOT / "docs/supporting/rgfca_sharedness_v2_post_capacity_selection_rule_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def row_hash(observation_id: int, photo_id: int) -> str:
    return hashlib.sha256(f"{int(observation_id)}:{int(photo_id)}".encode()).hexdigest()


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


def selected_species() -> pd.DataFrame:
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    capacity = json.loads(CAPACITY_RESULT.read_text(encoding="utf-8"))
    rule = json.loads(SELECTION_RULE.read_text(encoding="utf-8"))
    if capacity["status"] != c["activation"]["required_capacity_status"] or capacity["primary_feasibility_pass"] is not True:
        raise RuntimeError("capacity gate did not activate geometry freeze")
    if rule["status"] != "frozen_before_high_depth_capacity_pilot_final_result":
        raise RuntimeError("post-capacity selection rule drifted")
    frame = pd.read_csv(CAPACITY_ROWS)
    eligible = frame.loc[frame["confirmed_at_least_300"].astype(bool)].copy()
    eligible = eligible.sort_values("pilot_rank", kind="mergesort").head(int(rule["primary_species_sample"]["n"])).reset_index(drop=True)
    if len(eligible) != 150:
        raise RuntimeError("frozen 150-species selection is unavailable")
    eligible["species_order"] = range(150)
    eligible["split_hash"] = [hashlib.sha256(f"train-eval|{int(t)}".encode("ascii")).hexdigest() for t in eligible["inat_taxon_id"]]
    split_order = eligible.sort_values(["split_hash", "inat_taxon_id"], kind="mergesort").index.tolist()
    training = set(split_order[:75])
    eligible["sample_role"] = ["training" if i in training else "evaluation" for i in eligible.index]
    if (eligible["sample_role"] == "training").sum() != 75 or (eligible["sample_role"] == "evaluation").sum() != 75:
        raise RuntimeError("species-disjoint split failed")
    return eligible[["species_order", "sample_role", "species", "inat_taxon_id", "pilot_rank"]].copy()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--shard-index", type=int, required=True)
    p.add_argument("--shard-count", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c["status"] != "frozen_after_capacity_pass_before_any_v2_geometry_request":
        raise RuntimeError("geometry contract drifted")
    if c["image_pixels_opened"] is not False or c["flower_colour_opened"] is not False:
        raise RuntimeError("geometry contract opens prohibited outcomes")
    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != int(c["deterministic_shards"]) or not 0 <= shard_index < shard_count:
        raise RuntimeError("shard topology differs from contract")

    species = selected_species()
    shard = species.loc[(species["species_order"] % shard_count) == shard_index].copy().reset_index(drop=True)
    q = c["metadata_sampling"]
    allowed = frozenset(str(x).casefold() for x in q["allowed_photo_licenses"])
    excluded_obs, excluded_photo = load_exclusions([ROOT / p for p in q["prior_experiment_exclusion_sources"]])
    client = InaturalistObservationClient(
        request_interval_seconds=float(q["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(q["request_retries"]),
        user_agent=f"fcp-rgfca-sharedness-v2-geometry-{shard_index}/1.0 (github.com/zuizui0223/fcp)",
    )
    target = int(c["fixed_photo_frame"]["photos_per_species"])
    stop_n = int(q["stop_when_after_observer_cap_reaches"])
    max_batches = int(q["maximum_random_batches_per_species"])
    obs_frames: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []

    for pos, sr in enumerate(shard.itertuples(index=False), start=1):
        taxon_id = int(sr.inat_taxon_id)
        parsed: list[dict[str, object]] = []
        seen_obs: set[int] = set()
        seen_photo: set[int] = set()
        terminal_error = ""
        batches = 0
        raw_rows_seen = 0
        for batch in range(1, max_batches + 1):
            params = h9_query_for_species(
                taxon_id,
                per_page=int(q["per_page"]),
                maximum_positional_accuracy_m=int(q["maximum_positional_accuracy_m"]),
                flowering_term_id=int(q["flowering_term_id"]),
                flowering_term_value_id=int(q["flowering_term_value_id"]),
                allowed_photo_licenses=tuple(q["allowed_photo_licenses"]),
            )
            params["order_by"] = "random"
            params["page"] = 1
            batches += 1
            try:
                payload = client.observations(params)
            except Exception as exc:
                terminal_error = f"{type(exc).__name__}:{str(exc)[:180]}"
                break
            raw = payload.get("results") or []
            if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
                terminal_error = "RuntimeError:iNaturalist results is not a sequence"
                break
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
            if len(capped) >= stop_n:
                break
        frame = pd.DataFrame(parsed)
        capped = observer_cap(frame, int(q["observer_cap"])) if len(frame) else frame
        if len(capped):
            capped = capped.copy()
            if "row_hash" not in capped.columns:
                capped["row_hash"] = [row_hash(o, p) for o, p in zip(capped["observation_id"], capped["photo_id"])]
            selected = capped.sort_values("row_hash", kind="mergesort").head(target).copy().reset_index(drop=True)
            selected["sample_role"] = str(sr.sample_role)
            selected["species_order"] = int(sr.species_order)
            selected["photo_order"] = range(1, len(selected) + 1)
            allowed_fields = c["allowed_saved_fields"]
            for field in allowed_fields:
                if field not in selected.columns:
                    raise RuntimeError(f"required saved field unavailable: {field}")
            selected = selected[allowed_fields]
            obs_frames.append(selected)
        else:
            selected = capped
        audits.append({
            "species_order": int(sr.species_order),
            "sample_role": str(sr.sample_role),
            "species": str(sr.species),
            "inat_taxon_id": taxon_id,
            "random_batches_requested": batches,
            "raw_rows_seen": raw_rows_seen,
            "unique_locally_eligible_after_exclusion": int(len(frame)),
            "after_observer_cap": int(len(capped)),
            "retained": int(len(selected)),
            "full_fixed_300": bool(len(selected) == target),
            "terminal_request_error": terminal_error,
        })
        if pos % 5 == 0 or pos == len(shard):
            print(json.dumps({"shard": shard_index, "processed": pos, "shard_species": len(shard)}), flush=True)

    observations = pd.concat(obs_frames, ignore_index=True) if obs_frames else pd.DataFrame(columns=c["allowed_saved_fields"])
    audit = pd.DataFrame(audits).sort_values("species_order", kind="mergesort").reset_index(drop=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    obs_path = args.output_dir / f"rgfca_sharedness_v2_geometry_photos_shard_{shard_index:02d}.csv"
    audit_path = args.output_dir / f"rgfca_sharedness_v2_geometry_audit_shard_{shard_index:02d}.csv"
    manifest_path = args.output_dir / f"rgfca_sharedness_v2_geometry_shard_{shard_index:02d}.json"
    observations.to_csv(obs_path, index=False, lineterminator="\n")
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    manifest = {
        "protocol": c["protocol"],
        "status": "complete_metadata_only_geometry_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "shard_species": int(len(audit)),
        "full_fixed_300_species": int(audit["full_fixed_300"].sum()),
        "retained_photos": int(len(observations)),
        "terminal_request_error_species": int(audit["terminal_request_error"].fillna("").astype(str).str.len().gt(0).sum()),
        "image_pixels_opened": False,
        "flower_colour_opened": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "capacity_result_sha256": sha256_file(CAPACITY_RESULT),
            "capacity_rows_sha256": sha256_file(CAPACITY_ROWS),
            "selection_rule_sha256": sha256_file(SELECTION_RULE),
            "photos_csv_sha256": sha256_file(obs_path),
            "audit_csv_sha256": sha256_file(audit_path)
        }
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(manifest, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
