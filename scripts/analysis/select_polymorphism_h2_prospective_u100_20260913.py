#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CENSUS = ROOT / "results" / "rgfca_42111_breadth_depth_step8b_20260911" / "species_capacity_census.csv.gz"
DISC = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RES = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
OUT = ROOT / "results" / "polymorphism_h2_prospective_u100_selection_20260913"

U0_N = 42111
U100_MIN = 100
U100_N = 4730
P500_N = 500
SALT = "FCP_H2_PROSPECTIVE_20260913"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_species(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def selection_hash(taxon_id: int, species: str) -> str:
    return hashlib.sha256(f"{SALT}|{int(taxon_id)}|{species}".encode("utf-8")).hexdigest()


def main() -> None:
    # Outcome firewall: only these explicitly named columns are read.
    cap = pd.read_csv(CENSUS, usecols=["inat_taxon_id", "species", "after_observer_cap"])
    disc = pd.read_csv(DISC, usecols=["species"])
    res = pd.read_csv(RES, usecols=["species"])

    if len(cap) != U0_N or cap["inat_taxon_id"].nunique() != U0_N:
        raise RuntimeError(f"42,111 census fingerprint drift: rows={len(cap)}, unique taxon={cap['inat_taxon_id'].nunique()}")

    cap["species"] = norm_species(cap["species"])
    cap["inat_taxon_id"] = pd.to_numeric(cap["inat_taxon_id"], errors="raise").astype(int)
    cap["after_observer_cap"] = pd.to_numeric(cap["after_observer_cap"], errors="raise").astype(int)
    if cap["species"].eq("").any():
        raise RuntimeError("empty species name in capacity census")

    u100 = cap.loc[cap["after_observer_cap"] >= U100_MIN, ["inat_taxon_id", "species", "after_observer_cap"]].copy()
    if len(u100) != U100_N:
        raise RuntimeError(f"U100 fingerprint drift: {len(u100)} != {U100_N}")

    disc_species = set(norm_species(disc["species"])) - {""}
    res_species = set(norm_species(res["species"])) - {""}
    legacy_union = disc_species | res_species
    legacy_intersection = disc_species & res_species
    u100_species = set(u100["species"])

    disc_u100 = disc_species & u100_species
    res_u100 = res_species & u100_species
    excluded_u100 = legacy_union & u100_species
    legacy_outside_u100 = legacy_union - u100_species

    p100 = u100.loc[~u100["species"].isin(legacy_union)].copy()
    if len(p100) < P500_N:
        verdict = "PROSPECTIVE_SELECTION_UNDERIDENTIFIED"
        raise RuntimeError(f"{verdict}: P100={len(p100)} < {P500_N}")

    p100["selection_hash"] = [selection_hash(t, s) for t, s in zip(p100["inat_taxon_id"], p100["species"])]
    p100 = p100.sort_values(["selection_hash", "inat_taxon_id"], kind="mergesort").reset_index(drop=True)
    p100["prospective_rank"] = range(1, len(p100) + 1)
    p500 = p100.head(P500_N).copy()

    OUT.mkdir(parents=True, exist_ok=True)
    p100_path = OUT / "p100_outcome_blind_pool.csv"
    p500_path = OUT / "p500_frozen_selection.csv"
    p100.to_csv(p100_path, index=False)
    p500.to_csv(p500_path, index=False)

    p500_sha = sha256_file(p500_path)
    p100_sha = sha256_file(p100_path)

    result = {
        "analysis": "polymorphism_h2_prospective_u100_selection",
        "date_jst": "2026-09-13",
        "status": "complete_pre_outcome_species_selection",
        "protocol": "docs/POLYMORPHISM_H2_PROSPECTIVE_U100_SELECTION_PROTOCOL_20260913.md",
        "outcome_firewall": {
            "capacity_columns_read": ["inat_taxon_id", "species", "after_observer_cap"],
            "discovery_columns_read": ["species"],
            "reserve_columns_read": ["species"],
            "colour_outcome_columns_read": False,
            "morph_columns_read": False,
            "H2_geometry_read": False,
            "H3_predictors_read": False,
        },
        "sampling_frame": {
            "U0_species": int(len(cap)),
            "U100_min_after_observer_cap": U100_MIN,
            "U100_species": int(len(u100)),
        },
        "legacy_high_depth_identity_audit": {
            "discovery_species_in_measured_file": int(len(disc_species)),
            "reserve_species_in_measured_file": int(len(res_species)),
            "discovery_reserve_identity_overlap": int(len(legacy_intersection)),
            "legacy_union_species": int(len(legacy_union)),
            "discovery_species_inside_U100": int(len(disc_u100)),
            "reserve_species_inside_U100": int(len(res_u100)),
            "legacy_union_inside_U100_excluded": int(len(excluded_u100)),
            "legacy_union_outside_U100": int(len(legacy_outside_u100)),
        },
        "prospective_pool": {
            "P100_species": int(len(p100)),
            "P500_species": int(len(p500)),
            "selection_salt": SALT,
            "ordering": "ascending SHA256(salt|inat_taxon_id|species), then inat_taxon_id",
            "replacement_after_outcome_opening_allowed": False,
            "p100_csv": str(p100_path.relative_to(ROOT)),
            "p100_sha256": p100_sha,
            "p500_csv": str(p500_path.relative_to(ROOT)),
            "p500_sha256": p500_sha,
        },
        "downstream_status": {
            "P500_colour_outcomes_opened": False,
            "P500_H1_evaluated": False,
            "P500_H2_W_evaluated": False,
        },
        "selection_verdict": "P500_FROZEN_OUTCOME_BLIND",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# Prospective H2 U100 selection receipt",
        "",
        "Verdict: **P500_FROZEN_OUTCOME_BLIND**",
        "",
        f"- U0: {len(cap):,}",
        f"- U100: {len(u100):,}",
        f"- discovery species identities in measured file: {len(disc_species):,}",
        f"- reserve species identities in measured file: {len(res_species):,}",
        f"- discovery/reserve identity overlap: {len(legacy_intersection):,}",
        f"- legacy high-depth union excluded from U100: {len(excluded_u100):,}",
        f"- P100 prospective pool: {len(p100):,}",
        f"- frozen P500: {len(p500):,}",
        f"- P500 SHA256: `{p500_sha}`",
        "",
        "Selection read species identity and U100 capacity only. No colour-outcome, morph, H2-geometry, or H3-predictor column was read.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
