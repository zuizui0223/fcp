#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ALLOC = ROOT / "results/rgfca_42111_tiered_measurement_step8c_20260911/species_anchor_allocation.csv.gz"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
DISC_MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
RES_MEASURED = ROOT / "data/derived/rgfca_reserve_replication_measured_photos_v1.csv"
H9_LEDGER = ROOT / "data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv"
H9_META = ROOT / "data/frozen/random_photo_first_h9_fresh_metadata_v1.csv"
OUT = ROOT / "results/rgfca_42111_discovery_fill_step8c_20260911"
SEED = 20260911
N_SPECIES = 42111


def id_sets(path: Path) -> tuple[set[int], set[int]]:
    if not path.exists():
        return set(), set()
    q = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id"})
    obs = set(pd.to_numeric(q["observation_id"], errors="coerce").dropna().astype(int)) if "observation_id" in q else set()
    pho = set(pd.to_numeric(q["photo_id"], errors="coerce").dropna().astype(int)) if "photo_id" in q else set()
    return obs, pho


def hash_key(t: int, o: int, p: int) -> str:
    return hashlib.sha256(f"{SEED}|{t}|{o}|{p}".encode("utf-8")).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    alloc = pd.read_csv(ALLOC)
    if len(alloc) != N_SPECIES or alloc["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("allocation universe mismatch")
    target = dict(zip(alloc["inat_taxon_id"].astype(int), alloc["combined_target_depth"].astype(int)))
    opened_species = set(alloc.loc[alloc["opened_existing_species"].astype(bool), "species"].astype(str))

    parts = []
    for source, path in (("v1", V1), ("v2", V2)):
        x = pd.read_csv(path, usecols=["observation_id", "photo_id", "species", "inat_taxon_id"])
        x["source"] = source
        parts.append(x)
    cand = pd.concat(parts, ignore_index=True).dropna().copy()
    for c in ["observation_id", "photo_id", "inat_taxon_id"]:
        cand[c] = pd.to_numeric(cand[c], errors="raise").astype(int)
    cand = cand.drop_duplicates(["observation_id", "photo_id"], keep="first")
    cand = cand[cand["inat_taxon_id"].isin(target)].copy()
    cand["selection_hash"] = [hash_key(t, o, p) for t, o, p in zip(cand.inat_taxon_id, cand.observation_id, cand.photo_id)]

    excluded_obs: set[int] = set()
    excluded_pho: set[int] = set()
    for path in (DISC_MEASURED, RES_MEASURED, H9_LEDGER, H9_META):
        o, p = id_sets(path); excluded_obs |= o; excluded_pho |= p
    cand["prior_id_excluded"] = cand["observation_id"].isin(excluded_obs) | cand["photo_id"].isin(excluded_pho)

    rows = []
    selected_desc = []
    selected_fresh = []
    for taxon_id, g in cand.groupby("inat_taxon_id", sort=False, observed=True):
        k = int(target[int(taxon_id)])
        g = g.sort_values(["selection_hash", "observation_id", "photo_id"], kind="mergesort")
        gd = g.head(k)
        gf = g.loc[~g["prior_id_excluded"]].head(k)
        selected_desc.append(gd.assign(selection_surface="descriptive"))
        if len(gf): selected_fresh.append(gf.assign(selection_surface="fresh_id_excluded"))
        rows.append({
            "inat_taxon_id": int(taxon_id),
            "species": str(g["species"].iloc[0]),
            "target_depth": k,
            "discovery_candidate_count": int(len(g)),
            "discovery_fill_descriptive": int(len(gd)),
            "discovery_shortfall_descriptive": int(k - len(gd)),
            "fresh_candidate_count_after_prior_id_exclusion": int((~g["prior_id_excluded"]).sum()),
            "discovery_fill_fresh": int(len(gf)),
            "discovery_shortfall_fresh": int(k - len(gf)),
            "opened_existing_species": bool(str(g["species"].iloc[0]) in opened_species),
        })
    fill = pd.DataFrame(rows)
    if len(fill) != N_SPECIES:
        raise RuntimeError(f"candidate grouping lost species: {len(fill)}")

    desc = pd.concat(selected_desc, ignore_index=True)
    fresh = pd.concat(selected_fresh, ignore_index=True) if selected_fresh else pd.DataFrame(columns=desc.columns)
    desc.to_csv(OUT / "selected_discovery_candidates_descriptive.csv.gz", index=False, compression="gzip")
    fresh.to_csv(OUT / "selected_discovery_candidates_fresh.csv.gz", index=False, compression="gzip")
    fill.to_csv(OUT / "species_fill_audit.csv.gz", index=False, compression="gzip")

    unopened = fill.loc[~fill.opened_existing_species]
    result = {
        "analysis": "rgfca_42111_discovery_fill_step8c",
        "status": "complete_metadata_only_discovery_fill_audit",
        "species_universe": N_SPECIES,
        "target_records_total": int(fill.target_depth.sum()),
        "descriptive_records_filled_from_v1v2": int(fill.discovery_fill_descriptive.sum()),
        "descriptive_records_shortfall_after_v1v2": int(fill.discovery_shortfall_descriptive.sum()),
        "species_fully_filled_descriptive": int((fill.discovery_shortfall_descriptive == 0).sum()),
        "species_partially_filled_descriptive": int(((fill.discovery_fill_descriptive > 0) & (fill.discovery_shortfall_descriptive > 0)).sum()),
        "species_zero_fill_descriptive": int((fill.discovery_fill_descriptive == 0).sum()),
        "unopened_species": int(len(unopened)),
        "fresh_target_records_unopened": int(unopened.target_depth.sum()),
        "fresh_records_filled_from_v1v2_after_prior_id_exclusion": int(unopened.discovery_fill_fresh.sum()),
        "fresh_records_shortfall_requiring_new_query": int(unopened.discovery_shortfall_fresh.sum()),
        "unopened_species_fully_filled_from_v1v2": int((unopened.discovery_shortfall_fresh == 0).sum()),
        "unopened_species_with_any_new_query_shortfall": int((unopened.discovery_shortfall_fresh > 0).sum()),
        "pixels_opened": False,
        "flower_colour_used": False,
        "claim_boundary": "This is metadata allocation only. Fill counts are frozen candidate availability, not successful image resolution or colour measurement.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# RGFCA Step 8C — fill from frozen V1/V2 discovery IDs",
        "",
        f"- target records across 42,111 species: **{result['target_records_total']:,}**",
        f"- descriptive target records already fillable from frozen V1/V2 IDs: **{result['descriptive_records_filled_from_v1v2']:,}**",
        f"- descriptive shortfall requiring additional metadata query: **{result['descriptive_records_shortfall_after_v1v2']:,}**",
        f"- species fully filled descriptively from V1/V2: **{result['species_fully_filled_descriptive']:,} / 42,111**",
        "",
        f"- unopened species: **{result['unopened_species']:,}**",
        f"- fresh target records in unopened species: **{result['fresh_target_records_unopened']:,}**",
        f"- fresh records fillable from V1/V2 after prior-ID exclusion: **{result['fresh_records_filled_from_v1v2_after_prior_id_exclusion']:,}**",
        f"- fresh shortfall requiring a new query: **{result['fresh_records_shortfall_requiring_new_query']:,}**",
        f"- unopened species fully filled from V1/V2: **{result['unopened_species_fully_filled_from_v1v2']:,}**",
        f"- unopened species needing any new-query supplementation: **{result['unopened_species_with_any_new_query_shortfall']:,}**",
        "",
        "No image pixel or flower-colour outcome was opened.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
