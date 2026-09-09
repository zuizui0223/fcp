#!/usr/bin/env python3
"""Run one metadata-only recovery-acquisition shard for fresh RGFCA global axes.

This script implements the protocol frozen in
`docs/RGFCA_INDEPENDENT_GLOBAL_AXIS_METADATA_RECOVERY_ACQUISITION_PROTOCOL_20260909.md`.
It never dereferences image URLs and never computes colour or M1-M3 outcomes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.global_candidate_acquisition import deterministic_candidate_pages, stable_candidate_query
from fcp_pipeline.random_photo_h9_pool import parse_h9_observation
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
QUERY_CONTRACT = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_contract_v1.json"
QUALIFIER = ROOT / "scripts/analysis/qualify_rgfca_independent_global_axis_metadata_20260909.py"

PROTOCOL_COMMIT = "928b8eee67e93192232ad725126600c2da524628"
EXPECTED_CAPACITY_SPECIES = 4730
EXPECTED_PRIOR_CANDIDATE_SPECIES = 1000
EXPECTED_RECOVERY_SOURCE = 3730
MIN_CAPACITY_SPAN_KM = 500.0
REQUEST_INTERVAL_SECONDS = 1.30

PARSER_COLUMNS = [
    "species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
    "photo_license", "attribution", "latitude", "longitude", "positional_accuracy_m",
    "observed_on", "observer_id", "observer", "row_hash",
]


def load_qualifier():
    spec = importlib.util.spec_from_file_location("rgfca_independent_axis_qualifier", QUALIFIER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen metadata qualifier")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ids(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True).fillna("")


def id_set(frame: pd.DataFrame, column: str) -> set[str]:
    if column not in frame.columns:
        return set()
    return {x for x in ids(frame[column]).tolist() if x}


def _results(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    raw = payload.get("results") or []
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
        raise RuntimeError("iNaturalist results is not a sequence")
    return [row for row in raw if isinstance(row, Mapping)]


def numeric_taxon_sort(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="raise").astype(np.int64)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--shard-index", type=int, required=True)
    p.add_argument("--shard-count", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--capacity-selected", type=Path, required=True)
    p.add_argument("--prior-candidate-audit", type=Path, required=True)
    p.add_argument("--reserve", type=Path, required=True)
    p.add_argument("--h9-fresh", type=Path, required=True)
    p.add_argument("--prior-ledger", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.shard_count != 8:
        raise RuntimeError("recovery protocol requires exactly 8 shards")
    if not 0 <= args.shard_index < args.shard_count:
        raise RuntimeError("invalid shard index")
    for path in (
        args.capacity_selected, args.prior_candidate_audit, args.reserve,
        args.h9_fresh, args.prior_ledger, QUERY_CONTRACT, QUALIFIER,
    ):
        if not path.exists():
            raise RuntimeError(f"missing required metadata input: {path}")

    qmod = load_qualifier()
    contract = json.loads(QUERY_CONTRACT.read_text(encoding="utf-8"))
    q = contract["query"]
    if int(q["candidate_pages_per_species"]) != 3 or int(q["candidate_page_seed"]) != 20260917:
        raise RuntimeError("upstream fixed query contract drifted")

    capacity = pd.read_csv(args.capacity_selected)
    prior = pd.read_csv(args.prior_candidate_audit, usecols=["inat_taxon_id"])
    if capacity["inat_taxon_id"].nunique() != EXPECTED_CAPACITY_SPECIES:
        raise RuntimeError("capacity selected-species denominator drifted")
    if prior["inat_taxon_id"].nunique() != EXPECTED_PRIOR_CANDIDATE_SPECIES:
        raise RuntimeError("prior candidate denominator drifted")

    capacity["inat_taxon_id"] = ids(capacity["inat_taxon_id"])
    prior["inat_taxon_id"] = ids(prior["inat_taxon_id"])
    prior_ids = set(prior["inat_taxon_id"])
    capacity_ids = set(capacity["inat_taxon_id"])
    if not prior_ids.issubset(capacity_ids):
        raise RuntimeError("prior 1000 candidate taxa are not a subset of capacity-selected taxa")
    source = capacity.loc[~capacity["inat_taxon_id"].isin(prior_ids)].copy()
    if source["inat_taxon_id"].nunique() != EXPECTED_RECOVERY_SOURCE:
        raise RuntimeError("recovery source is not exactly 3730 unqueried taxa")
    source["maximum_span_km"] = pd.to_numeric(source["maximum_span_km"], errors="raise")
    source["span_possible"] = source["maximum_span_km"].ge(MIN_CAPACITY_SPAN_KM)
    source = source.loc[source["span_possible"]].copy()
    source["taxon_sort"] = numeric_taxon_sort(source["inat_taxon_id"])
    source = source.sort_values(["taxon_sort", "species"], kind="mergesort").reset_index(drop=True)
    source["recovery_row_index"] = np.arange(len(source), dtype=int)
    shard = source.loc[(source["recovery_row_index"] % args.shard_count) == args.shard_index].copy()
    if shard.empty:
        raise RuntimeError("recovery shard unexpectedly empty")

    reserve = pd.read_csv(
        args.reserve,
        usecols=["species", "inat_taxon_id", "observation_id", "photo_id", "observer_id"],
    )
    h9 = pd.read_csv(
        args.h9_fresh,
        usecols=["species", "inat_taxon_id", "observation_id", "photo_id", "observer_id"],
    )
    ledger = pd.read_csv(args.prior_ledger, usecols=["observation_id", "photo_id", "source"])
    for frame in (reserve, h9):
        for col in ("inat_taxon_id", "observation_id", "photo_id", "observer_id"):
            frame[col] = ids(frame[col])
    for col in ("observation_id", "photo_id"):
        ledger[col] = ids(ledger[col])

    reserve_taxa = id_set(reserve, "inat_taxon_id")
    reserve_species = set(reserve["species"].dropna().astype(str))
    reserve_obs = id_set(reserve, "observation_id")
    reserve_photo = id_set(reserve, "photo_id")
    reserve_observer = id_set(reserve, "observer_id")
    h9_obs = id_set(h9, "observation_id")
    h9_photo = id_set(h9, "photo_id")
    h9_observer = id_set(h9, "observer_id")
    ledger_obs = id_set(ledger, "observation_id")
    ledger_photo = id_set(ledger, "photo_id")

    source_overlap = set(shard["inat_taxon_id"]) & reserve_taxa
    if source_overlap or set(shard["species"].astype(str)) & reserve_species:
        raise RuntimeError("unqueried recovery source unexpectedly overlaps reserve species")

    allowed = frozenset(str(x).casefold() for x in q["allowed_photo_licenses"])
    client = InaturalistObservationClient(
        request_interval_seconds=REQUEST_INTERVAL_SECONDS,
        timeout_seconds=45.0,
        max_retries=int(q["request_retries"]),
        user_agent=f"fcp-rgfca-global-axis-recovery-s{args.shard_index}/1.0 (github.com/zuizui0223/fcp)",
    )

    request_attempts = 0
    request_errors = 0
    audit_rows: list[dict[str, object]] = []
    group_rows: list[dict[str, object]] = []
    photo_parts: list[pd.DataFrame] = []

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
        species = str(row.species)
        taxon_text = str(row.inat_taxon_id)
        taxon_id = int(taxon_text)
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
                payload, error = request(taxon_id, int(page))
                if payload is not None:
                    payload_by_page[int(page)] = payload
                if error:
                    page_errors[int(page)] = error

        raw_rows: list[Mapping[str, object]] = []
        for page in pages:
            payload = payload_by_page.get(int(page))
            if payload is not None:
                raw_rows.extend(_results(payload))

        parsed_rows: list[dict[str, object]] = []
        locally_eligible = 0
        wrong_taxon = 0
        excluded_obs = 0
        excluded_photo = 0
        excluded_observer = 0
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
            obs_id = str(parsed["observation_id"])
            photo_id = str(parsed["photo_id"])
            observer_id = str(parsed["observer_id"])
            if obs_id in reserve_obs or obs_id in h9_obs or obs_id in ledger_obs:
                excluded_obs += 1
                continue
            if photo_id in reserve_photo or photo_id in h9_photo or photo_id in ledger_photo:
                excluded_photo += 1
                continue
            if observer_id in reserve_observer or observer_id in h9_observer:
                excluded_observer += 1
                continue
            parsed_rows.append(parsed)

        clean = pd.DataFrame(parsed_rows, columns=PARSER_COLUMNS)
        if len(clean):
            for col in ("inat_taxon_id", "observation_id", "photo_id", "observer_id"):
                clean[col] = ids(clean[col])
            clean = clean.drop_duplicates("photo_id", keep="first")
            clean = clean.drop_duplicates("observation_id", keep="first").reset_index(drop=True)
            x, y = qmod.epsg6933_xy(clean["longitude"].to_numpy(float), clean["latitude"].to_numpy(float))
            clean["x_m"] = x
            clean["y_m"] = y
            clean["gx"] = np.floor(clean["x_m"] / qmod.CELL_M).astype(int)
            clean["gy"] = np.floor(clean["y_m"] / qmod.CELL_M).astype(int)

        local_groups = []
        if len(clean):
            for (gx, gy), g in clean.groupby(["gx", "gy"], sort=True):
                n_photo = int(g["photo_id"].nunique())
                n_obs = int(g["observation_id"].nunique())
                n_observer = int(g["observer_id"].nunique())
                qualified = (
                    n_photo >= qmod.MIN_PHOTOS_PER_GROUP
                    and n_obs >= qmod.MIN_PHOTOS_PER_GROUP
                    and n_observer >= qmod.MIN_OBSERVERS_PER_GROUP
                )
                local_groups.append({
                    "inat_taxon_id": taxon_text,
                    "species": species,
                    "gx": int(gx),
                    "gy": int(gy),
                    "n_photos": n_photo,
                    "n_observations": n_obs,
                    "n_observers": n_observer,
                    "centroid_x_m": float(g["x_m"].mean()),
                    "centroid_y_m": float(g["y_m"].mean()),
                    "qualified_group": bool(qualified),
                })
        gframe = pd.DataFrame(local_groups)
        qgroups = gframe.loc[gframe["qualified_group"]].copy() if len(gframe) else gframe
        max_sep_m = qmod.max_pairwise_distance(
            qgroups["centroid_x_m"].to_numpy(float),
            qgroups["centroid_y_m"].to_numpy(float),
        ) if len(qgroups) else 0.0
        qualifies = bool(len(qgroups) >= qmod.MIN_GROUPS and max_sep_m >= qmod.SEPARATION_M)

        retained_groups = 0
        retained_photos = 0
        if qualifies:
            chosen_cells = qmod.choose_dispersed_groups(qgroups)
            for group_order, (gx, gy) in enumerate(chosen_cells, start=1):
                gr = qgroups.loc[(qgroups["gx"] == gx) & (qgroups["gy"] == gy)].iloc[0].to_dict()
                gr.update({
                    "selected_group": True,
                    "selected_group_order": int(group_order),
                    "recovery_row_index": int(row.recovery_row_index),
                    "shard_index": int(args.shard_index),
                })
                group_rows.append(gr)
                source_group = clean.loc[(clean["gx"] == gx) & (clean["gy"] == gy)].copy()
                six = qmod.select_six_photos(source_group, species, int(gx), int(gy))
                six["selected_group_order"] = int(group_order)
                six["recovery_row_index"] = int(row.recovery_row_index)
                six["shard_index"] = int(args.shard_index)
                six["ecological_image_pixels_opened"] = False
                six["flower_colour_used"] = False
                six["M1_M3_scores_used"] = False
                photo_parts.append(six)
                retained_groups += 1
                retained_photos += len(six)

        audit_rows.append({
            "recovery_row_index": int(row.recovery_row_index),
            "shard_index": int(args.shard_index),
            "species": species,
            "inat_taxon_id": taxon_text,
            "capacity_after_observer_cap": int(row.after_observer_cap),
            "capacity_maximum_span_km": float(row.maximum_span_km),
            "total_results_probe": int(total_results),
            "candidate_pages": json.dumps(list(pages), separators=(",", ":")),
            "candidate_pages_requested": int(len(pages)),
            "candidate_page_errors": json.dumps(page_errors, sort_keys=True, separators=(",", ":")),
            "raw_results_from_pages": int(len(raw_rows)),
            "locally_eligible": int(locally_eligible),
            "wrong_taxon": int(wrong_taxon),
            "excluded_opened_observation": int(excluded_obs),
            "excluded_opened_photo": int(excluded_photo),
            "excluded_known_opened_observer": int(excluded_observer),
            "identity_clean_rows": int(len(clean)),
            "identity_clean_distinct_observers": int(clean["observer_id"].nunique()) if len(clean) else 0,
            "qualifying_spatial_groups": int(len(qgroups)),
            "max_qualifying_group_centroid_separation_km": float(max_sep_m / 1000.0),
            "spatial_species_qualified": qualifies,
            "retained_support_groups": int(retained_groups),
            "retained_photo_rows": int(retained_photos),
        })

        if processed % 25 == 0 or processed == len(shard):
            print(json.dumps({
                "shard": args.shard_index,
                "processed": processed,
                "shard_species": len(shard),
                "request_attempts": request_attempts,
                "request_errors": request_errors,
                "spatial_species_qualified": int(sum(bool(r["spatial_species_qualified"]) for r in audit_rows)),
            }), flush=True)

    audit = pd.DataFrame(audit_rows).sort_values("recovery_row_index", kind="mergesort").reset_index(drop=True)
    groups = pd.DataFrame(group_rows)
    if len(groups):
        groups = groups.sort_values(["recovery_row_index", "selected_group_order"], kind="mergesort").reset_index(drop=True)
    photos = pd.concat(photo_parts, ignore_index=True) if photo_parts else pd.DataFrame()
    if len(photos):
        photos = photos.sort_values(
            ["recovery_row_index", "selected_group_order", "observer_id", "photo_id"], kind="mergesort"
        ).reset_index(drop=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = args.output_dir / f"recovery_shard_{args.shard_index:02d}_species.csv"
    groups_path = args.output_dir / f"recovery_shard_{args.shard_index:02d}_groups.csv"
    photos_path = args.output_dir / f"recovery_shard_{args.shard_index:02d}_photos.csv"
    manifest_path = args.output_dir / f"recovery_shard_{args.shard_index:02d}.json"
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    groups.to_csv(groups_path, index=False, lineterminator="\n")
    photos.to_csv(photos_path, index=False, lineterminator="\n")

    shard_manifest = {
        "protocol": "rgfca-independent-global-axis-metadata-recovery-acquisition-v1",
        "protocol_commit": PROTOCOL_COMMIT,
        "status": "complete_metadata_only_recovery_shard",
        "shard_index": int(args.shard_index),
        "shard_count": int(args.shard_count),
        "source_species_before_span_filter": EXPECTED_RECOVERY_SOURCE,
        "source_species_after_span_filter_global": int(len(source)),
        "shard_species": int(len(shard)),
        "request_attempts": int(request_attempts),
        "request_errors": int(request_errors),
        "spatial_species_qualified": int(audit["spatial_species_qualified"].sum()),
        "retained_support_groups": int(len(groups)),
        "retained_photo_rows": int(len(photos)),
        "candidate_image_pixels_opened": False,
        "ecological_image_pixels_opened": False,
        "flower_colour_used": False,
        "M1_M3_scores_used": False,
        "lineage": {
            "capacity_selected_sha256": sha256_file(args.capacity_selected),
            "prior_candidate_audit_sha256": sha256_file(args.prior_candidate_audit),
            "reserve_sha256": sha256_file(args.reserve),
            "h9_fresh_sha256": sha256_file(args.h9_fresh),
            "prior_ledger_sha256": sha256_file(args.prior_ledger),
            "query_contract_sha256": sha256_file(QUERY_CONTRACT),
            "species_audit_sha256": sha256_file(audit_path),
            "groups_sha256": sha256_file(groups_path),
            "photos_sha256": sha256_file(photos_path),
        },
    }
    manifest_path.write_text(json.dumps(shard_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(shard_manifest, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
