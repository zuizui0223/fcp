#!/usr/bin/env python3
"""Metadata-only prequalification for the independent RGFCA global-axis replication.

No image pixels, colour measurements, M1-M3 scores, island status, climate, or
previous ecological outcomes are read here. The runner uses the exact fixed
0-shift 300-km EPSG:6933 lattice from the discovery recurrence code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

PROTOCOL_COMMIT = "6ea2783970d0e0fdea959a4ce13abaafc4e9b5d9"
CELL_M = 300_000.0
SEPARATION_M = 500_000.0
TARGET_SPECIES = 300
MIN_GROUPS = 3
MIN_PHOTOS_PER_GROUP = 6
MIN_OBSERVERS_PER_GROUP = 3
MAX_PHOTOS_PER_OBSERVER_GROUP = 2
MAX_RETAINED_GROUPS = 5

# Exact WGS84 / EPSG:6933 constants copied from the frozen discovery runner.
_A = 6378137.0
_F = 1.0 / 298.257223563
_E = math.sqrt(_F * (2.0 - _F))
_PHI1 = math.radians(30.0)
_K0 = math.cos(_PHI1) / math.sqrt(1.0 - _E * _E * math.sin(_PHI1) ** 2)


def epsg6933_xy(longitude_deg: np.ndarray, latitude_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lam = np.radians(np.asarray(longitude_deg, dtype=float))
    phi = np.radians(np.asarray(latitude_deg, dtype=float))
    s = np.sin(phi)
    q = (1.0 - _E * _E) * (
        s / (1.0 - _E * _E * s * s)
        - (1.0 / (2.0 * _E)) * np.log((1.0 - _E * s) / (1.0 + _E * s))
    )
    x = _A * _K0 * lam
    y = _A * q / (2.0 * _K0)
    return x, y


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def id_strings(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)
    return out.fillna("")


def id_set(frame: pd.DataFrame, column: str) -> set[str]:
    if column not in frame.columns:
        return set()
    return {x for x in id_strings(frame[column]).tolist() if x}


def stable_key(*parts: object) -> str:
    return hashlib.sha256("|".join(str(x) for x in parts).encode("utf-8")).hexdigest()


def max_pairwise_distance(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return 0.0
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return float(np.sqrt(dx * dx + dy * dy).max())


def choose_dispersed_groups(groups: pd.DataFrame) -> list[tuple[int, int]]:
    g = groups.sort_values(["gx", "gy"], kind="mergesort").reset_index(drop=True)
    cells = [(int(r.gx), int(r.gy)) for r in g.itertuples(index=False)]
    if len(cells) <= MAX_RETAINED_GROUPS:
        return cells
    xy = g[["centroid_x_m", "centroid_y_m"]].to_numpy(float)
    dist = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    best_pair = None
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            pair_cells = tuple(sorted((cells[i], cells[j])))
            candidate = (float(dist[i, j]), pair_cells, i, j)
            if best_pair is None or candidate[0] > best_pair[0] + 1e-9 or (
                abs(candidate[0] - best_pair[0]) <= 1e-9 and candidate[1] < best_pair[1]
            ):
                best_pair = candidate
    assert best_pair is not None
    selected = [best_pair[2], best_pair[3]]
    while len(selected) < MAX_RETAINED_GROUPS:
        options = []
        for i in range(len(g)):
            if i in selected:
                continue
            min_d = float(np.min(dist[i, selected]))
            options.append((-min_d, cells[i], i))
        options.sort()
        selected.append(options[0][2])
    return [cells[i] for i in selected]


def select_six_photos(group: pd.DataFrame, species: str, gx: int, gy: int) -> pd.DataFrame:
    g = group.drop_duplicates("photo_id", keep="first").drop_duplicates("observation_id", keep="first").copy()
    g["observer_id"] = id_strings(g["observer_id"])
    g["photo_id"] = id_strings(g["photo_id"])
    g["observation_id"] = id_strings(g["observation_id"])
    buckets: dict[str, list[int]] = {}
    for observer, h in g.groupby("observer_id", sort=False):
        rows = list(h.index)
        rows.sort(key=lambda idx: stable_key(species, gx, gy, observer, g.at[idx, "photo_id"], g.at[idx, "observation_id"]))
        buckets[observer] = rows[:MAX_PHOTOS_PER_OBSERVER_GROUP]
    observers = sorted(buckets, key=lambda o: stable_key(species, gx, gy, o))
    chosen: list[int] = []
    for round_index in range(MAX_PHOTOS_PER_OBSERVER_GROUP):
        for observer in observers:
            rows = buckets[observer]
            if round_index < len(rows):
                chosen.append(rows[round_index])
                if len(chosen) == MIN_PHOTOS_PER_GROUP:
                    out = g.loc[chosen].copy()
                    if out["observer_id"].nunique() < MIN_OBSERVERS_PER_GROUP:
                        raise RuntimeError("deterministic photo selection violated observer floor")
                    return out
    raise RuntimeError("qualified support group could not supply exact six-photo round-robin sample")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--reserve", type=Path, required=True)
    p.add_argument("--h9-fresh", type=Path, required=True)
    p.add_argument("--prior-ledger", type=Path, required=True)
    p.add_argument("--output-root", type=Path, default=Path("."))
    p.add_argument("--protocol-commit", default=PROTOCOL_COMMIT)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    for path in (args.candidate, args.reserve, args.h9_fresh, args.prior_ledger):
        if not path.exists():
            raise RuntimeError(f"missing required metadata input: {path}")

    candidate_cols = [
        "species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
        "photo_license", "attribution", "latitude", "longitude", "positional_accuracy_m",
        "observed_on", "observer_id", "observer", "row_hash", "global_row_index", "shard_index",
    ]
    cand = pd.read_csv(args.candidate, usecols=candidate_cols)
    reserve = pd.read_csv(
        args.reserve,
        usecols=["species", "inat_taxon_id", "observation_id", "photo_id", "observer_id"],
    )
    h9 = pd.read_csv(
        args.h9_fresh,
        usecols=["species", "inat_taxon_id", "observation_id", "photo_id", "observer_id"],
    )
    ledger = pd.read_csv(args.prior_ledger, usecols=["observation_id", "photo_id", "source"])

    if len(cand) != 100_000:
        raise RuntimeError(f"candidate row count drifted: {len(cand)}")
    if cand["inat_taxon_id"].nunique() != 1000:
        raise RuntimeError("candidate species denominator drifted from 1000 taxon IDs")
    if len(reserve) != 50_000 or reserve["inat_taxon_id"].nunique() != 500:
        raise RuntimeError("reserve identity denominator drifted")

    for col in ("inat_taxon_id", "observation_id", "photo_id", "observer_id"):
        cand[col] = id_strings(cand[col])
        reserve[col] = id_strings(reserve[col])
        h9[col] = id_strings(h9[col])
    for col in ("observation_id", "photo_id"):
        ledger[col] = id_strings(ledger[col])

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

    raw_candidate_taxa = set(cand["inat_taxon_id"])
    raw_candidate_species = set(cand["species"].dropna().astype(str))
    taxon_overlap = raw_candidate_taxa & reserve_taxa
    species_name_overlap = raw_candidate_species & reserve_species

    cand["excluded_reserve_species"] = cand["inat_taxon_id"].isin(reserve_taxa) | cand["species"].isin(reserve_species)
    cand["excluded_opened_observation"] = cand["observation_id"].isin(reserve_obs | h9_obs | ledger_obs)
    cand["excluded_opened_photo"] = cand["photo_id"].isin(reserve_photo | h9_photo | ledger_photo)
    cand["excluded_known_opened_observer"] = cand["observer_id"].isin(reserve_observer | h9_observer)
    cand["identity_excluded"] = cand[
        [
            "excluded_reserve_species", "excluded_opened_observation", "excluded_opened_photo",
            "excluded_known_opened_observer",
        ]
    ].any(axis=1)

    work = cand.loc[~cand["identity_excluded"]].copy()
    required_nonmissing = [
        "species", "inat_taxon_id", "observation_id", "photo_id", "observer_id", "latitude", "longitude"
    ]
    for col in required_nonmissing:
        if col in {"latitude", "longitude"}:
            work = work.loc[pd.to_numeric(work[col], errors="coerce").notna()].copy()
        else:
            work = work.loc[work[col].astype(str).str.len().gt(0)].copy()
    work = work.drop_duplicates("photo_id", keep="first").drop_duplicates("observation_id", keep="first").copy()
    work["latitude"] = pd.to_numeric(work["latitude"], errors="raise")
    work["longitude"] = pd.to_numeric(work["longitude"], errors="raise")
    x, y = epsg6933_xy(work["longitude"].to_numpy(), work["latitude"].to_numpy())
    work["x_m"] = x
    work["y_m"] = y
    work["gx"] = np.floor(work["x_m"] / CELL_M).astype(int)
    work["gy"] = np.floor(work["y_m"] / CELL_M).astype(int)

    group_rows = []
    for (taxon_id, species, gx, gy), g in work.groupby(
        ["inat_taxon_id", "species", "gx", "gy"], sort=True
    ):
        n_photo = int(g["photo_id"].nunique())
        n_obs = int(g["observation_id"].nunique())
        n_observer = int(g["observer_id"].nunique())
        qualified = (
            n_photo >= MIN_PHOTOS_PER_GROUP
            and n_obs >= MIN_PHOTOS_PER_GROUP
            and n_observer >= MIN_OBSERVERS_PER_GROUP
        )
        group_rows.append({
            "inat_taxon_id": taxon_id,
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
    groups = pd.DataFrame(group_rows)
    qualified_groups = groups.loc[groups["qualified_group"]].copy()

    audit_rows = []
    all_taxon_rows = (
        cand[["inat_taxon_id", "species"]]
        .drop_duplicates()
        .sort_values(["inat_taxon_id", "species"], kind="mergesort")
    )
    raw_by = {
        (str(k[0]), str(k[1])): g
        for k, g in cand.groupby(["inat_taxon_id", "species"], sort=False)
    }
    eligible_by = {
        (str(k[0]), str(k[1])): g
        for k, g in work.groupby(["inat_taxon_id", "species"], sort=False)
    }
    qg_by = {
        (str(k[0]), str(k[1])): g
        for k, g in qualified_groups.groupby(["inat_taxon_id", "species"], sort=False)
    }
    empty_work = work.iloc[0:0]
    empty_qg = qualified_groups.iloc[0:0]
    for row in all_taxon_rows.itertuples(index=False):
        taxon_id = str(row.inat_taxon_id)
        species = str(row.species)
        key = (taxon_id, species)
        raw = raw_by[key]
        eligible = eligible_by.get(key, empty_work)
        qg = qg_by.get(key, empty_qg)
        max_sep = max_pairwise_distance(
            qg["centroid_x_m"].to_numpy(float), qg["centroid_y_m"].to_numpy(float)
        ) if len(qg) else 0.0
        qualifies = len(qg) >= MIN_GROUPS and max_sep >= SEPARATION_M
        audit_rows.append({
            "inat_taxon_id": taxon_id,
            "species": species,
            "raw_candidate_rows": int(len(raw)),
            "excluded_reserve_species": bool(raw["excluded_reserve_species"].any()),
            "n_rows_excluded_opened_observation": int(raw["excluded_opened_observation"].sum()),
            "n_rows_excluded_opened_photo": int(raw["excluded_opened_photo"].sum()),
            "n_rows_excluded_known_opened_observer": int(raw["excluded_known_opened_observer"].sum()),
            "identity_clean_candidate_rows": int(len(eligible)),
            "identity_clean_distinct_observers": int(eligible["observer_id"].nunique()),
            "qualifying_spatial_groups": int(len(qg)),
            "max_qualifying_group_centroid_separation_km": float(max_sep / 1000.0),
            "spatial_species_qualified": bool(qualifies),
        })
    audit = pd.DataFrame(audit_rows)
    ranked = audit.loc[audit["spatial_species_qualified"]].copy()
    ranked = ranked.sort_values(
        [
            "qualifying_spatial_groups", "identity_clean_distinct_observers",
            "max_qualifying_group_centroid_separation_km", "inat_taxon_id",
        ],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    ranked["selection_rank"] = np.arange(1, len(ranked) + 1)
    audit = audit.merge(ranked[["inat_taxon_id", "species", "selection_rank"]], on=["inat_taxon_id", "species"], how="left")

    spatial_prequalification_pass = len(ranked) >= TARGET_SPECIES
    selected = ranked.head(TARGET_SPECIES).copy() if spatial_prequalification_pass else ranked.iloc[0:0].copy()
    selected_keys = set(zip(selected["inat_taxon_id"], selected["species"]))

    selected_group_rows = []
    manifest_parts = []
    cell_source_by = {
        (str(k[0]), str(k[1]), int(k[2]), int(k[3])): g
        for k, g in work.groupby(["inat_taxon_id", "species", "gx", "gy"], sort=False)
    }
    if spatial_prequalification_pass:
        for taxon_id, species in sorted(selected_keys, key=lambda z: ((0, int(z[0])) if z[0].isdigit() else (1, z[0]), z[1])):
            qg = qg_by[(taxon_id, species)].copy()
            chosen_cells = choose_dispersed_groups(qg)
            for order_index, (gx, gy) in enumerate(chosen_cells, start=1):
                gr = qg.loc[(qg["gx"] == gx) & (qg["gy"] == gy)].iloc[0]
                selected_group_rows.append({
                    **gr.to_dict(),
                    "selected_group": True,
                    "selected_group_order": order_index,
                })
                source = cell_source_by[(taxon_id, species, gx, gy)].copy()
                six = select_six_photos(source, species, gx, gy)
                six["selected_group_order"] = order_index
                manifest_parts.append(six)

    selected_groups = pd.DataFrame(selected_group_rows)
    manifest = pd.concat(manifest_parts, ignore_index=True) if manifest_parts else pd.DataFrame(columns=list(work.columns) + ["selected_group_order"])
    if len(manifest):
        manifest = manifest.sort_values(
            ["inat_taxon_id", "selected_group_order", "gx", "gy", "observer_id", "photo_id"],
            kind="mergesort",
        ).reset_index(drop=True)
        manifest["ecological_image_pixels_opened"] = False
        manifest["flower_colour_used"] = False
        manifest["M1_M3_scores_used"] = False

    # Known identity sources are enforced. The legacy pre-H7 ledger does not retain
    # observer IDs, so final all-opened-frame observer clearance is intentionally false.
    pre_h7_rows = int(ledger["source"].astype(str).eq("pre_h7_experiments").sum())
    pre_h7_observer_identity_complete = False if pre_h7_rows > 0 else True
    all_opened_frame_observer_identity_complete = False
    final_identity_gate_pass = False
    ecological_metadata_frame_qualified = bool(spatial_prequalification_pass and final_identity_gate_pass)

    if len(manifest):
        overlap_checks = {
            "reserve_species_taxon_overlap": int(len(set(manifest["inat_taxon_id"]) & reserve_taxa)),
            "reserve_photo_overlap": int(len(set(manifest["photo_id"]) & reserve_photo)),
            "reserve_observation_overlap": int(len(set(manifest["observation_id"]) & reserve_obs)),
            "reserve_observer_overlap": int(len(set(manifest["observer_id"]) & reserve_observer)),
            "h9_photo_overlap": int(len(set(manifest["photo_id"]) & h9_photo)),
            "h9_observation_overlap": int(len(set(manifest["observation_id"]) & h9_obs)),
            "h9_observer_overlap": int(len(set(manifest["observer_id"]) & h9_observer)),
            "prior_ledger_photo_overlap": int(len(set(manifest["photo_id"]) & ledger_photo)),
            "prior_ledger_observation_overlap": int(len(set(manifest["observation_id"]) & ledger_obs)),
        }
    else:
        overlap_checks = {k: 0 for k in [
            "reserve_species_taxon_overlap", "reserve_photo_overlap", "reserve_observation_overlap",
            "reserve_observer_overlap", "h9_photo_overlap", "h9_observation_overlap", "h9_observer_overlap",
            "prior_ledger_photo_overlap", "prior_ledger_observation_overlap",
        ]}
    if any(overlap_checks.values()):
        raise RuntimeError(f"retained prequalification manifest has a known identity overlap: {overlap_checks}")

    root = args.output_root
    audit_path = root / "data/frozen/rgfca_independent_global_axis_species_audit_20260909.csv"
    groups_path = root / "data/frozen/rgfca_independent_global_axis_support_groups_20260909.csv"
    manifest_path = root / "data/frozen/rgfca_independent_global_axis_prequalified_photo_manifest_20260909.csv"
    result_path = root / "docs/supporting/rgfca_independent_global_axis_metadata_prequalification_result_20260909.json"
    for path in (audit_path, groups_path, manifest_path, result_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    audit = audit.sort_values(["selection_rank", "inat_taxon_id", "species"], na_position="last", kind="mergesort")
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selected_groups.to_csv(groups_path, index=False, lineterminator="\n")
    manifest.to_csv(manifest_path, index=False, lineterminator="\n")

    result = {
        "protocol": "rgfca-independent-global-signal-axis-replication-metadata-prequalification-v1",
        "protocol_commit": args.protocol_commit,
        "status": (
            "spatial_prequalification_passed_final_identity_gate_blocked"
            if spatial_prequalification_pass
            else "spatial_prequalification_failed"
        ),
        "candidate_image_pixels_opened": False,
        "ecological_image_pixels_opened": False,
        "flower_colour_used": False,
        "M1_M3_scores_used": False,
        "climate_used": False,
        "island_status_used": False,
        "measurement_gate_passed": False,
        "candidate_rows": int(len(cand)),
        "candidate_species": int(cand["inat_taxon_id"].nunique()),
        "reserve_rows": int(len(reserve)),
        "reserve_species": int(reserve["inat_taxon_id"].nunique()),
        "candidate_reserve_taxon_overlap_before_exclusion": int(len(taxon_overlap)),
        "candidate_reserve_species_name_overlap_before_exclusion": int(len(species_name_overlap)),
        "rows_after_known_identity_exclusions": int(len(work)),
        "species_after_known_identity_exclusions": int(work["inat_taxon_id"].nunique()),
        "spatial_species_qualified_count": int(len(ranked)),
        "spatial_target_species": TARGET_SPECIES,
        "spatial_prequalification_pass": bool(spatial_prequalification_pass),
        "selected_species": int(manifest["inat_taxon_id"].nunique()) if len(manifest) else 0,
        "selected_support_groups": int(len(selected_groups)),
        "selected_photo_rows": int(len(manifest)),
        "selected_observers": int(manifest["observer_id"].nunique()) if len(manifest) else 0,
        "pre_h7_observer_identity_complete": bool(pre_h7_observer_identity_complete),
        "all_opened_frame_observer_identity_complete": bool(all_opened_frame_observer_identity_complete),
        "final_identity_gate_pass": bool(final_identity_gate_pass),
        "ecological_metadata_frame_qualified": bool(ecological_metadata_frame_qualified),
        "opening_authorized": False,
        "block_reason": (
            "legacy_pre_h7_and_other_opened_model_validation_observer_identity_not_yet_fully_reconstructed"
            if spatial_prequalification_pass
            else "fewer_than_300_species_met_frozen_spatial_support_rule"
        ),
        "known_identity_overlap_checks": overlap_checks,
        "lattice": {
            "crs": "EPSG:6933",
            "cell_m": CELL_M,
            "origin_shift_m": [0.0, 0.0],
            "assignment": "floor(x/300000),floor(y/300000)",
            "minimum_groups_per_species": MIN_GROUPS,
            "minimum_photos_per_group": MIN_PHOTOS_PER_GROUP,
            "minimum_observers_per_group": MIN_OBSERVERS_PER_GROUP,
            "minimum_centroid_pair_separation_m": SEPARATION_M,
        },
        "selection": {
            "rank": [
                "descending_qualifying_spatial_groups",
                "descending_identity_clean_distinct_observers",
                "descending_max_centroid_separation",
                "ascending_stable_taxon_id",
            ],
            "retained_groups_per_species": "min(5, qualifying_group_count) using deterministic farthest-point dispersion",
            "photos_per_retained_group": MIN_PHOTOS_PER_GROUP,
            "observer_rule": "deterministic round-robin, max 2 photos per observer per group",
        },
        "lineage": {
            "candidate_path": str(args.candidate),
            "candidate_sha256": sha256_file(args.candidate),
            "reserve_path": str(args.reserve),
            "reserve_sha256": sha256_file(args.reserve),
            "h9_fresh_path": str(args.h9_fresh),
            "h9_fresh_sha256": sha256_file(args.h9_fresh),
            "prior_ledger_path": str(args.prior_ledger),
            "prior_ledger_sha256": sha256_file(args.prior_ledger),
        },
    }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
