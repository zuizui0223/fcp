#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
SPECIES_ANCHOR = ROOT / "results/rgfca_42111_tiered_measurement_step8c_20260911/species_anchor_allocation.csv.gz"
DIAGNOSTIC = ROOT / "results/rgfca_42111_taxon_cell_conflict_diagnostic_step8e2_20260911/result.json"
ADDENDUM = ROOT / "docs/RGFCA_42111_GEOGRAPHIC_BREADTH_STEP8E2_UNIQUENESS_ADDENDUM_20260911.md"
OUT = ROOT / "results/rgfca_42111_taxon_cell_anchor_step8e2_20260911"
SEED = 20260911
EXPECTED_PAIRS = 85337
EXPECTED_CELLS = 128
EXPECTED_SPECIES = 42111
EXPECTED_CANDIDATES = 364772


def candidate_hash(cell_id: int, taxon_id: int, observation_id: int, photo_id: int) -> str:
    return hashlib.sha256(
        f"{SEED}|{cell_id}|{taxon_id}|{observation_id}|{photo_id}".encode("utf-8")
    ).hexdigest()


def pair_hash(cell_id: int, taxon_id: int) -> str:
    return hashlib.sha256(f"{SEED}|pair|{cell_id}|{taxon_id}".encode("utf-8")).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    if diagnostic.get("status") != "complete_metadata_only_conflict_diagnostic":
        raise RuntimeError("conflict diagnostic is not complete")
    if int(diagnostic.get("taxon_cell_pairs", -1)) != EXPECTED_PAIRS:
        raise RuntimeError("diagnostic pair denominator drift")
    if int(diagnostic.get("candidate_rows", -1)) != EXPECTED_CANDIDATES:
        raise RuntimeError("diagnostic candidate denominator drift")
    if int(diagnostic.get("deterministic_greedy_unique_pairs_matched", -1)) != EXPECTED_PAIRS:
        raise RuntimeError("diagnostic did not establish complete unique matching feasibility")
    if int(diagnostic.get("deterministic_greedy_unresolved_pairs", -1)) != 0:
        raise RuntimeError("diagnostic retained unresolved unique-matching pairs")
    if not ADDENDUM.exists():
        raise RuntimeError("pre-colour uniqueness addendum is missing")

    parts: list[pd.DataFrame] = []
    for source, path in (("v1", V1), ("v2", V2)):
        x = pd.read_csv(path, usecols=["cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"])
        x["discovery_source"] = source
        parts.append(x)
    x = pd.concat(parts, ignore_index=True)
    x = x.dropna(subset=["cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"]).copy()
    for c in ["cell_id", "observation_id", "photo_id", "inat_taxon_id"]:
        x[c] = pd.to_numeric(x[c], errors="raise").astype(int)
    # Preserve candidate identities across pairs. Only exact duplicate candidate rows are removed.
    x = x.drop_duplicates(["cell_id", "inat_taxon_id", "observation_id", "photo_id"], keep="first").copy()
    if len(x) != EXPECTED_CANDIDATES:
        raise RuntimeError(f"candidate-row count drift: {len(x)} != {EXPECTED_CANDIDATES}")

    unique_pairs = x[["inat_taxon_id", "cell_id"]].drop_duplicates()
    if len(unique_pairs) != EXPECTED_PAIRS:
        raise RuntimeError(f"taxon-cell pair count drift: {len(unique_pairs)} != {EXPECTED_PAIRS}")
    if unique_pairs["cell_id"].nunique() != EXPECTED_CELLS:
        raise RuntimeError("occupied discovery-cell count drift")
    if unique_pairs["inat_taxon_id"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError("species universe drift in taxon-cell links")

    x["candidate_hash"] = [
        candidate_hash(int(c), int(t), int(o), int(p))
        for c, t, o, p in zip(x["cell_id"], x["inat_taxon_id"], x["observation_id"], x["photo_id"])
    ]
    x = x.sort_values(
        ["inat_taxon_id", "cell_id", "candidate_hash", "observation_id", "photo_id"],
        kind="mergesort",
    ).reset_index(drop=True)

    groups: dict[tuple[int, int], pd.DataFrame] = {
        (int(taxon_id), int(cell_id)): g
        for (taxon_id, cell_id), g in x.groupby(["inat_taxon_id", "cell_id"], sort=False, observed=True)
    }
    ordered_pairs = sorted(groups, key=lambda k: pair_hash(k[1], k[0]))

    used_observation_ids: set[int] = set()
    used_photo_ids: set[int] = set()
    rows: list[dict[str, object]] = []
    unresolved: list[dict[str, object]] = []
    for taxon_id, cell_id in ordered_pairs:
        g = groups[(taxon_id, cell_id)]
        chosen = None
        for r in g.itertuples(index=False):
            oid = int(r.observation_id)
            pid = int(r.photo_id)
            if oid in used_observation_ids or pid in used_photo_ids:
                continue
            chosen = r
            break
        if chosen is None:
            unresolved.append({
                "inat_taxon_id": taxon_id,
                "cell_id": cell_id,
                "candidate_rows_in_taxon_cell": int(len(g)),
            })
            continue
        used_observation_ids.add(int(chosen.observation_id))
        used_photo_ids.add(int(chosen.photo_id))
        rows.append({
            "inat_taxon_id": taxon_id,
            "species": str(chosen.species),
            "cell_id": cell_id,
            "observation_id": int(chosen.observation_id),
            "photo_id": int(chosen.photo_id),
            "discovery_source": str(chosen.discovery_source),
            "candidate_rows_in_taxon_cell": int(len(g)),
            "candidate_hash": str(chosen.candidate_hash),
            "pair_hash": pair_hash(cell_id, taxon_id),
        })

    unresolved_df = pd.DataFrame(unresolved, columns=["inat_taxon_id", "cell_id", "candidate_rows_in_taxon_cell"])
    unresolved_df.to_csv(OUT / "unresolved_taxon_cell_pairs.csv", index=False, lineterminator="\n")
    out = pd.DataFrame(rows).sort_values(["cell_id", "inat_taxon_id"], kind="mergesort").reset_index(drop=True)
    if len(out) != EXPECTED_PAIRS or len(unresolved_df) != 0:
        raise RuntimeError("constrained unique matching did not retain all 85,337 pairs")
    if out[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != EXPECTED_PAIRS:
        raise RuntimeError("frozen frame is not one row per taxon-cell pair")
    if out["observation_id"].nunique() != EXPECTED_PAIRS or out["photo_id"].nunique() != EXPECTED_PAIRS:
        raise RuntimeError("constrained frame does not have globally unique observation/photo IDs")

    s = pd.read_csv(SPECIES_ANCHOR, usecols=["inat_taxon_id", "anchor_any_observation_id", "anchor_any_photo_id"])
    s["anchor_any_observation_id"] = pd.to_numeric(s["anchor_any_observation_id"], errors="raise").astype(int)
    s["anchor_any_photo_id"] = pd.to_numeric(s["anchor_any_photo_id"], errors="raise").astype(int)
    pair_keys = set(zip(out["observation_id"].astype(int), out["photo_id"].astype(int)))
    species_keys = set(zip(s["anchor_any_observation_id"].astype(int), s["anchor_any_photo_id"].astype(int)))
    overlap = len(pair_keys & species_keys)

    richness = out.groupby("cell_id", observed=True)["inat_taxon_id"].nunique().sort_index()
    result = {
        "analysis": "rgfca_42111_taxon_cell_anchor_step8e2",
        "status": "complete_metadata_only_unique_taxon_cell_anchor_freeze",
        "candidate_rows": int(len(x)),
        "taxon_cell_anchors": int(len(out)),
        "species": int(out["inat_taxon_id"].nunique()),
        "occupied_cells": int(out["cell_id"].nunique()),
        "unresolved_pairs": int(len(unresolved_df)),
        "unique_observation_ids": int(out["observation_id"].nunique()),
        "unique_photo_ids": int(out["photo_id"].nunique()),
        "species_anchor_overlap_rows": int(overlap),
        "additional_unique_anchors_beyond_species_breadth": int(len(out) - overlap),
        "cell_species_richness": {
            "min": int(richness.min()),
            "median": float(richness.median()),
            "max": int(richness.max()),
        },
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "selection_rule": "ascending SHA256 pair order; within pair ascending SHA256 candidate order; first observation/photo identity unused globally",
        "boundary": "This frame is for cell-wise geographic breadth. It must not replace the one-anchor-per-species estimator for global species-equal composition."
    }
    out.to_csv(OUT / "taxon_cell_anchor_frame.csv.gz", index=False, compression="gzip", lineterminator="\n")
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8E2 — unique taxon×cell geographic breadth anchor frame\n\n"
        f"- candidate rows: **{len(x):,}**\n"
        f"- taxon×cell anchors: **{len(out):,} / {EXPECTED_PAIRS:,}**\n"
        f"- unresolved pairs: **{len(unresolved_df):,}**\n"
        f"- unique observation IDs: **{out['observation_id'].nunique():,}**\n"
        f"- unique photo IDs: **{out['photo_id'].nunique():,}**\n"
        f"- species represented: **{out['inat_taxon_id'].nunique():,}**\n"
        f"- occupied cells: **{out['cell_id'].nunique()}**\n"
        f"- overlap with one-per-species anchors: **{overlap:,}**\n"
        f"- additional anchors required for full cell representation: **{len(out)-overlap:,}**\n"
        f"- cell species richness min / median / max: **{int(richness.min())} / {float(richness.median()):.1f} / {int(richness.max())}**\n"
        "- image pixels opened: **false**\n"
        "- flower colour used: **false**\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
