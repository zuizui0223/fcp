#!/usr/bin/env python3
"""Reduce all frozen RGFCA global-axis metadata recovery shards.

Metadata only: never dereferences image URLs and never computes colour or M1-M3.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

PROTOCOL_COMMIT = "928b8eee67e93192232ad725126600c2da524628"
EXPECTED_SHARDS = 8
EXPECTED_SOURCE_BEFORE_SPAN = 3730
TARGET_SPECIES = 300
MAX_ERROR_FRACTION = 0.05
EXPECTED_CAPACITY_SHA = "850501d7293a3968cf6566f8e9d7a396afd3927e8166ab2468e58ed22afcabde"
EXPECTED_PRIOR_AUDIT_SHA = "c9e9e76883d18bbb06355e1dec46a83aaf8f9a8680fb9b2c99d3cc811d62b152"
EXPECTED_RESERVE_SHA = "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6"
EXPECTED_H9_SHA = "111d0f964618c0c3df749a6e4bd29f834214cda5d7a2d9bda7d43cdc9dbf4c6f"
EXPECTED_LEDGER_SHA = "f9a6894740e9974399c055f92cba237be8ada41707d84e1807ba61b902c91b99"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ids(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True).fillna("")


def read_csv_allow_empty(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except EmptyDataError:
        return pd.DataFrame()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--output-root", type=Path, default=Path("."))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    manifest_paths = sorted(args.input_dir.glob("recovery_shard_*.json"))
    if len(manifest_paths) != EXPECTED_SHARDS:
        raise RuntimeError(f"expected {EXPECTED_SHARDS} shard manifests, found {len(manifest_paths)}")

    manifests = [json.loads(p.read_text(encoding="utf-8")) for p in manifest_paths]
    shard_ids = sorted(int(m["shard_index"]) for m in manifests)
    if shard_ids != list(range(EXPECTED_SHARDS)):
        raise RuntimeError(f"shard index set drifted: {shard_ids}")

    source_after_values = {int(m["source_species_after_span_filter_global"]) for m in manifests}
    if len(source_after_values) != 1:
        raise RuntimeError(f"shards disagree on post-span source denominator: {source_after_values}")
    source_after_span = source_after_values.pop()

    for m in manifests:
        if m.get("protocol") != "rgfca-independent-global-axis-metadata-recovery-acquisition-v1":
            raise RuntimeError("recovery shard protocol drifted")
        if m.get("protocol_commit") != PROTOCOL_COMMIT:
            raise RuntimeError("recovery shard protocol commit drifted")
        if int(m.get("shard_count", -1)) != EXPECTED_SHARDS:
            raise RuntimeError("recovery shard count drifted")
        if int(m.get("source_species_before_span_filter", -1)) != EXPECTED_SOURCE_BEFORE_SPAN:
            raise RuntimeError("recovery source denominator drifted")
        for flag in (
            "candidate_image_pixels_opened", "ecological_image_pixels_opened",
            "flower_colour_used", "M1_M3_scores_used",
        ):
            if m.get(flag) is not False:
                raise RuntimeError(f"sealed metadata boundary violated in shard {m['shard_index']}: {flag}")
        lineage = m["lineage"]
        expected = {
            "capacity_selected_sha256": EXPECTED_CAPACITY_SHA,
            "prior_candidate_audit_sha256": EXPECTED_PRIOR_AUDIT_SHA,
            "reserve_sha256": EXPECTED_RESERVE_SHA,
            "h9_fresh_sha256": EXPECTED_H9_SHA,
            "prior_ledger_sha256": EXPECTED_LEDGER_SHA,
        }
        for key, value in expected.items():
            if lineage.get(key) != value:
                raise RuntimeError(f"shard {m['shard_index']} lineage mismatch for {key}")

    audits, groups_parts, photo_parts = [], [], []
    for shard_id in range(EXPECTED_SHARDS):
        ap = args.input_dir / f"recovery_shard_{shard_id:02d}_species.csv"
        gp = args.input_dir / f"recovery_shard_{shard_id:02d}_groups.csv"
        pp = args.input_dir / f"recovery_shard_{shard_id:02d}_photos.csv"
        for p in (ap, gp, pp):
            if not p.exists():
                raise RuntimeError(f"missing shard output {p}")
        a = read_csv_allow_empty(ap)
        g = read_csv_allow_empty(gp)
        ph = read_csv_allow_empty(pp)
        manifest = next(m for m in manifests if int(m["shard_index"]) == shard_id)
        if len(a) != int(manifest["shard_species"]):
            raise RuntimeError(f"shard {shard_id} audit denominator mismatch")
        if len(g) != int(manifest["retained_support_groups"]):
            raise RuntimeError(f"shard {shard_id} group denominator mismatch")
        if len(ph) != int(manifest["retained_photo_rows"]):
            raise RuntimeError(f"shard {shard_id} photo denominator mismatch")
        audits.append(a)
        if len(g):
            groups_parts.append(g)
        if len(ph):
            photo_parts.append(ph)

    audit = pd.concat(audits, ignore_index=True)
    groups = pd.concat(groups_parts, ignore_index=True) if groups_parts else pd.DataFrame()
    photos = pd.concat(photo_parts, ignore_index=True) if photo_parts else pd.DataFrame()

    if len(audit) != source_after_span:
        raise RuntimeError(f"pooled audit denominator {len(audit)} != frozen post-span source {source_after_span}")
    if audit["recovery_row_index"].nunique() != len(audit):
        raise RuntimeError("recovery row index duplicated across shards")
    if audit["inat_taxon_id"].nunique() != len(audit):
        raise RuntimeError("taxon ID duplicated across recovery shards")

    for frame in (audit, groups, photos):
        if len(frame) and "inat_taxon_id" in frame:
            frame["inat_taxon_id"] = ids(frame["inat_taxon_id"])
    if len(photos):
        for col in ("photo_id", "observation_id", "observer_id"):
            photos[col] = ids(photos[col])
        if photos["photo_id"].duplicated().any() or photos["observation_id"].duplicated().any():
            raise RuntimeError("selected shard metadata contain duplicate photo/observation IDs")
        for flag in ("ecological_image_pixels_opened", "flower_colour_used", "M1_M3_scores_used"):
            if flag not in photos.columns or photos[flag].fillna(True).astype(bool).any():
                raise RuntimeError(f"photo metadata opening firewall violated: {flag}")

    qualified = audit.loc[audit["spatial_species_qualified"].astype(bool)].copy()
    if len(qualified):
        group_counts = groups.groupby(["inat_taxon_id", "species"]).size() if len(groups) else pd.Series(dtype=int)
        photo_counts = photos.groupby(["inat_taxon_id", "species"]).size() if len(photos) else pd.Series(dtype=int)
        for r in qualified.itertuples(index=False):
            key = (str(r.inat_taxon_id), str(r.species))
            expected_groups = min(5, int(r.qualifying_spatial_groups))
            got_groups = int(group_counts.get(key, 0))
            got_photos = int(photo_counts.get(key, 0))
            if got_groups != expected_groups:
                raise RuntimeError(f"retained group count mismatch for {key}: {got_groups} != {expected_groups}")
            if got_photos != 6 * expected_groups:
                raise RuntimeError(f"retained photo count mismatch for {key}: {got_photos} != {6 * expected_groups}")

    request_attempts = sum(int(m["request_attempts"]) for m in manifests)
    request_errors = sum(int(m["request_errors"]) for m in manifests)
    error_fraction = float(request_errors / request_attempts) if request_attempts else 1.0
    request_gate_pass = bool(request_attempts > 0 and error_fraction <= MAX_ERROR_FRACTION)

    ranked = qualified.sort_values(
        [
            "qualifying_spatial_groups", "identity_clean_distinct_observers",
            "max_qualifying_group_centroid_separation_km", "inat_taxon_id",
        ],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    ranked["selection_rank"] = np.arange(1, len(ranked) + 1)
    audit = audit.merge(
        ranked[["inat_taxon_id", "species", "selection_rank"]],
        on=["inat_taxon_id", "species"], how="left",
    )

    spatial_pass = bool(request_gate_pass and len(ranked) >= TARGET_SPECIES)
    selected = ranked.head(TARGET_SPECIES).copy() if spatial_pass else ranked.iloc[0:0].copy()
    selected_keys = set(zip(selected["inat_taxon_id"], selected["species"]))

    if spatial_pass:
        if len(groups) == 0 or len(photos) == 0:
            raise RuntimeError("spatial pass cannot have empty retained group/photo metadata")
        gmask = pd.Series(list(zip(groups["inat_taxon_id"], groups["species"])), index=groups.index).isin(selected_keys)
        pmask = pd.Series(list(zip(photos["inat_taxon_id"], photos["species"])), index=photos.index).isin(selected_keys)
        selected_groups = groups.loc[gmask].copy()
        selected_photos = photos.loc[pmask].copy()
    else:
        selected_groups = groups.iloc[0:0].copy()
        selected_photos = photos.iloc[0:0].copy()

    if spatial_pass:
        if selected_photos["inat_taxon_id"].nunique() != TARGET_SPECIES:
            raise RuntimeError("recovery pass did not freeze exactly 300 species")
        if not (5400 <= len(selected_photos) <= 9000):
            raise RuntimeError("selected recovery photo-row count outside frozen 300-species bounds")
        if len(selected_photos) != 6 * len(selected_groups):
            raise RuntimeError("selected recovery support groups do not have exactly six photos each")

    pre_h7_observer_identity_complete = False
    all_opened_frame_observer_identity_complete = False
    final_identity_gate_pass = False
    measurement_gate_passed = False
    opening_authorized = False
    ecological_metadata_frame_qualified = False

    if not request_gate_pass:
        status = "not_evaluable_recovery_due_request_failure"
        block_reason = "aggregate_request_error_fraction_above_0_05"
    elif not spatial_pass:
        status = "spatial_recovery_prequalification_failed"
        block_reason = "fewer_than_300_species_met_frozen_spatial_support_rule"
    else:
        status = "spatial_recovery_prequalification_passed_final_identity_gate_blocked"
        block_reason = "legacy_pre_h7_and_other_opened_model_validation_observer_identity_not_yet_fully_reconstructed"

    out = args.output_root
    audit_path = out / "data/frozen/rgfca_independent_global_axis_metadata_recovery_species_audit_20260909.csv"
    groups_path = out / "data/frozen/rgfca_independent_global_axis_metadata_recovery_support_groups_20260909.csv"
    photos_path = out / "data/frozen/rgfca_independent_global_axis_metadata_recovery_prequalified_photo_manifest_20260909.csv"
    result_path = out / "docs/supporting/rgfca_independent_global_axis_metadata_recovery_result_20260909.json"
    for p in (audit_path, groups_path, photos_path, result_path):
        p.parent.mkdir(parents=True, exist_ok=True)

    audit = audit.sort_values(["selection_rank", "recovery_row_index"], na_position="last", kind="mergesort")
    if len(selected_groups):
        selected_groups = selected_groups.sort_values(
            ["inat_taxon_id", "selected_group_order", "gx", "gy"], kind="mergesort"
        ).reset_index(drop=True)
    if len(selected_photos):
        selected_photos = selected_photos.sort_values(
            ["inat_taxon_id", "selected_group_order", "observer_id", "photo_id"], kind="mergesort"
        ).reset_index(drop=True)

    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selected_groups.to_csv(groups_path, index=False, lineterminator="\n")
    selected_photos.to_csv(photos_path, index=False, lineterminator="\n")

    upstream_lineage = manifests[0]["lineage"]
    result = {
        "protocol": "rgfca-independent-global-axis-metadata-recovery-acquisition-v1",
        "protocol_commit": PROTOCOL_COMMIT,
        "status": status,
        "candidate_image_pixels_opened": False,
        "ecological_image_pixels_opened": False,
        "flower_colour_used": False,
        "M1_M3_scores_used": False,
        "climate_used": False,
        "island_status_used": False,
        "source_species_before_span_filter": EXPECTED_SOURCE_BEFORE_SPAN,
        "source_species_after_span_filter": int(source_after_span),
        "request_attempts": int(request_attempts),
        "request_errors": int(request_errors),
        "request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": MAX_ERROR_FRACTION,
        "request_gate_pass": request_gate_pass,
        "spatial_species_qualified_count": int(len(ranked)),
        "spatial_target_species": TARGET_SPECIES,
        "spatial_recovery_prequalification_pass": spatial_pass,
        "selected_species": int(selected_photos["inat_taxon_id"].nunique()) if len(selected_photos) else 0,
        "selected_support_groups": int(len(selected_groups)),
        "selected_photo_rows": int(len(selected_photos)),
        "selected_observers": int(selected_photos["observer_id"].nunique()) if len(selected_photos) else 0,
        "pre_h7_observer_identity_complete": pre_h7_observer_identity_complete,
        "all_opened_frame_observer_identity_complete": all_opened_frame_observer_identity_complete,
        "final_identity_gate_pass": final_identity_gate_pass,
        "measurement_gate_passed": measurement_gate_passed,
        "ecological_metadata_frame_qualified": ecological_metadata_frame_qualified,
        "opening_authorized": opening_authorized,
        "block_reason": block_reason,
        "selection": {
            "rank": [
                "descending_qualifying_spatial_groups",
                "descending_identity_clean_distinct_observers",
                "descending_max_centroid_separation",
                "ascending_stable_taxon_id",
            ],
            "retained_groups_per_species": "min(5, qualifying_group_count) using deterministic farthest-point dispersion",
            "photos_per_retained_group": 6,
            "observer_rule": "deterministic round-robin, max 2 photos per observer per group",
        },
        "lattice": {
            "crs": "EPSG:6933",
            "cell_m": 300000.0,
            "origin_shift_m": [0.0, 0.0],
            "assignment": "floor(x/300000),floor(y/300000)",
            "minimum_groups_per_species": 3,
            "minimum_photos_per_group": 6,
            "minimum_observers_per_group": 3,
            "minimum_centroid_pair_separation_m": 500000.0,
        },
        "lineage": {
            "capacity_selected_sha256": upstream_lineage["capacity_selected_sha256"],
            "prior_candidate_audit_sha256": upstream_lineage["prior_candidate_audit_sha256"],
            "reserve_sha256": upstream_lineage["reserve_sha256"],
            "h9_fresh_sha256": upstream_lineage["h9_fresh_sha256"],
            "prior_ledger_sha256": upstream_lineage["prior_ledger_sha256"],
            "query_contract_sha256": upstream_lineage["query_contract_sha256"],
            "shard_manifest_sha256": [sha256_file(p) for p in manifest_paths],
            "species_audit_sha256": sha256_file(audit_path),
            "support_groups_sha256": sha256_file(groups_path),
            "photo_manifest_sha256": sha256_file(photos_path),
        },
    }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
