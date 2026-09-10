#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CENSUS = ROOT / "results/rgfca_42111_breadth_depth_step8b_20260911/species_capacity_census.csv.gz"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
DISC_MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
RES_MEASURED = ROOT / "data/derived/rgfca_reserve_replication_measured_photos_v1.csv"
H9_LEDGER = ROOT / "data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv"
H9_META = ROOT / "data/frozen/random_photo_first_h9_fresh_metadata_v1.csv"
OUT = ROOT / "results/rgfca_42111_tiered_measurement_step8c_20260911"
SEED = 20260911
N_SPECIES = 42111


def depth_target(n: int) -> int:
    if n >= 20:
        return 20
    if n >= 10:
        return 10
    if n >= 5:
        return 5
    if n >= 2:
        return 2
    if n >= 1:
        return 1
    return 0


def read_ids(path: Path) -> tuple[set[int], set[int], set[str]]:
    if not path.exists():
        return set(), set(), set()
    df = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id", "species"})
    obs = set(pd.to_numeric(df["observation_id"], errors="coerce").dropna().astype(int)) if "observation_id" in df else set()
    photos = set(pd.to_numeric(df["photo_id"], errors="coerce").dropna().astype(int)) if "photo_id" in df else set()
    species = set(df["species"].dropna().astype(str)) if "species" in df else set()
    return obs, photos, species


def choose_hash(group: pd.DataFrame) -> pd.Series:
    keys = [
        hashlib.sha256(f"{SEED}|{int(t)}|{int(o)}|{int(p)}".encode("utf-8")).hexdigest()
        for t, o, p in zip(group["inat_taxon_id"], group["observation_id"], group["photo_id"])
    ]
    return group.iloc[int(np.argmin(np.asarray(keys, dtype=object)))]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    census = pd.read_csv(CENSUS)
    if len(census) != N_SPECIES or census["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("Step 8B census is not the exact 42,111-species universe")

    parts = []
    for source, path in (("v1", V1), ("v2", V2)):
        x = pd.read_csv(path, usecols=["observation_id", "photo_id", "species", "inat_taxon_id"])
        x["discovery_source"] = source
        parts.append(x)
    cand = pd.concat(parts, ignore_index=True)
    cand = cand.dropna(subset=["observation_id", "photo_id", "inat_taxon_id", "species"]).copy()
    for c in ["observation_id", "photo_id", "inat_taxon_id"]:
        cand[c] = pd.to_numeric(cand[c], errors="raise").astype(int)
    cand = cand.drop_duplicates(["observation_id", "photo_id"], keep="first")
    universe = set(census["inat_taxon_id"].astype(int))
    cand = cand[cand["inat_taxon_id"].isin(universe)].copy()
    if cand["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("some discovered species have no V1/V2 photo anchor metadata")

    excluded_obs: set[int] = set()
    excluded_photo: set[int] = set()
    opened_species: set[str] = set()
    for path in (DISC_MEASURED, RES_MEASURED, H9_LEDGER, H9_META):
        obs, photos, species = read_ids(path)
        excluded_obs |= obs
        excluded_photo |= photos
        if path in (DISC_MEASURED, RES_MEASURED):
            opened_species |= species

    cand["previously_opened_id"] = cand["observation_id"].isin(excluded_obs) | cand["photo_id"].isin(excluded_photo)
    counts = cand.groupby("inat_taxon_id", observed=True).agg(
        discovery_anchor_candidates=("photo_id", "nunique"),
        previously_opened_anchor_candidates=("previously_opened_id", "sum"),
    ).reset_index()

    any_rows = []
    fresh_rows = []
    for taxon_id, g in cand.groupby("inat_taxon_id", sort=False, observed=True):
        q = choose_hash(g)
        any_rows.append({
            "inat_taxon_id": int(taxon_id),
            "anchor_any_observation_id": int(q["observation_id"]),
            "anchor_any_photo_id": int(q["photo_id"]),
            "anchor_any_source": str(q["discovery_source"]),
        })
        f = g.loc[~g["previously_opened_id"]]
        if len(f):
            r = choose_hash(f)
            fresh_rows.append({
                "inat_taxon_id": int(taxon_id),
                "anchor_fresh_observation_id": int(r["observation_id"]),
                "anchor_fresh_photo_id": int(r["photo_id"]),
                "anchor_fresh_source": str(r["discovery_source"]),
            })
    any_df = pd.DataFrame(any_rows)
    fresh_df = pd.DataFrame(fresh_rows)

    out = census.merge(counts, on="inat_taxon_id", how="left", validate="one_to_one")
    out = out.merge(any_df, on="inat_taxon_id", how="left", validate="one_to_one")
    out = out.merge(fresh_df, on="inat_taxon_id", how="left", validate="one_to_one")
    out["opened_existing_species"] = out["species"].astype(str).isin(opened_species)
    out["capacity_depth_target"] = out["after_observer_cap"].astype(int).map(depth_target)
    out["breadth_anchor_available"] = out["anchor_any_photo_id"].notna()
    out["fresh_discovery_anchor_available"] = out["anchor_fresh_photo_id"].notna()
    # Every discovered species has an original discovery anchor. A capacity-zero species may therefore
    # still support a one-photo breadth observation from the frozen discovery draw.
    out["combined_target_depth"] = np.maximum(1, out["capacity_depth_target"].to_numpy(int))
    out["new_measurement_target_if_existing_reused"] = np.where(
        out["opened_existing_species"], 0, out["combined_target_depth"]
    ).astype(int)
    out["anchor_status"] = np.where(
        out["opened_existing_species"],
        "opened_existing_species",
        np.where(out["fresh_discovery_anchor_available"], "fresh_discovery_anchor_available", "fresh_requery_required"),
    )

    if not out["breadth_anchor_available"].all():
        raise RuntimeError("not all 42,111 species have an original discovery photo anchor")

    target_counts = {str(k): int(v) for k, v in out["combined_target_depth"].value_counts().sort_index().items()}
    capacity_target_counts = {str(k): int(v) for k, v in out["capacity_depth_target"].value_counts().sort_index().items()}
    result = {
        "analysis": "rgfca_42111_anchor_allocation_step8c",
        "status": "complete_metadata_only_anchor_allocation",
        "species_universe": N_SPECIES,
        "pixels_opened": False,
        "flower_colour_used": False,
        "all_species_have_original_discovery_anchor": bool(out["breadth_anchor_available"].all()),
        "species_with_fresh_discovery_anchor_after_prior_id_exclusion": int(out["fresh_discovery_anchor_available"].sum()),
        "species_requiring_fresh_requery_for_unopened_anchor": int((~out["opened_existing_species"] & ~out["fresh_discovery_anchor_available"]).sum()),
        "opened_existing_species": int(out["opened_existing_species"].sum()),
        "capacity_depth_target_counts": capacity_target_counts,
        "combined_target_depth_counts_with_discovery_anchor_floor1": target_counts,
        "capacity_tier_total_records": int(out["capacity_depth_target"].sum()),
        "combined_total_target_records": int(out["combined_target_depth"].sum()),
        "new_measurement_target_records_if_existing_reused": int(out["new_measurement_target_if_existing_reused"].sum()),
        "claim_boundary": "Metadata-only allocation. Original discovery photo IDs establish breadth-anchor availability but no image or colour outcome is opened here.",
    }
    out.to_csv(OUT / "species_anchor_allocation.csv.gz", index=False, compression="gzip")
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# RGFCA Step 8C — 42,111-species anchor allocation",
        "",
        f"- species universe: **{N_SPECIES:,}**",
        f"- all species have an original frozen discovery photo anchor: **{result['all_species_have_original_discovery_anchor']}**",
        f"- opened existing species: **{result['opened_existing_species']:,}**",
        f"- fresh discovery anchors after prior-ID exclusion: **{result['species_with_fresh_discovery_anchor_after_prior_id_exclusion']:,} species**",
        f"- unopened species requiring fresh requery for an anchor: **{result['species_requiring_fresh_requery_for_unopened_anchor']:,}**",
        f"- capacity-tier target records: **{result['capacity_tier_total_records']:,}**",
        f"- combined target after giving capacity-zero species their frozen discovery anchor: **{result['combined_total_target_records']:,}**",
        f"- new target records if the existing opened cohort is reused descriptively: **{result['new_measurement_target_records_if_existing_reused']:,}**",
        "",
        "No image pixel or flower-colour result was opened in this allocation step.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
