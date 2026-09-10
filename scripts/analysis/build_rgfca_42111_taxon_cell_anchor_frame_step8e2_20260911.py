#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
SPECIES_ANCHOR = ROOT / "results/rgfca_42111_tiered_measurement_step8c_20260911/species_anchor_allocation.csv.gz"
OUT = ROOT / "results/rgfca_42111_taxon_cell_anchor_step8e2_20260911"
SEED = 20260911
EXPECTED_PAIRS = 85337
EXPECTED_CELLS = 128
EXPECTED_SPECIES = 42111


def choose(group: pd.DataFrame) -> pd.Series:
    keys = [
        hashlib.sha256(
            f"{SEED}|{int(c)}|{int(t)}|{int(o)}|{int(p)}".encode("utf-8")
        ).hexdigest()
        for c, t, o, p in zip(
            group["cell_id"], group["inat_taxon_id"], group["observation_id"], group["photo_id"]
        )
    ]
    return group.iloc[int(np.argmin(np.asarray(keys, dtype=object)))]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    parts = []
    for source, path in (("v1", V1), ("v2", V2)):
        x = pd.read_csv(path, usecols=["cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"])
        x["discovery_source"] = source
        parts.append(x)
    x = pd.concat(parts, ignore_index=True)
    x = x.dropna(subset=["cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"]).copy()
    for c in ["cell_id", "observation_id", "photo_id", "inat_taxon_id"]:
        x[c] = pd.to_numeric(x[c], errors="raise").astype(int)
    x = x.drop_duplicates(["observation_id", "photo_id"], keep="first")

    unique_pairs = x[["inat_taxon_id", "cell_id"]].drop_duplicates()
    if len(unique_pairs) != EXPECTED_PAIRS:
        raise RuntimeError(f"taxon-cell pair count drift: {len(unique_pairs)} != {EXPECTED_PAIRS}")
    if unique_pairs["cell_id"].nunique() != EXPECTED_CELLS:
        raise RuntimeError("occupied discovery-cell count drift")
    if unique_pairs["inat_taxon_id"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError("species universe drift in taxon-cell links")

    rows = []
    for (taxon_id, cell_id), g in x.groupby(["inat_taxon_id", "cell_id"], sort=False, observed=True):
        r = choose(g)
        rows.append({
            "inat_taxon_id": int(taxon_id),
            "species": str(r["species"]),
            "cell_id": int(cell_id),
            "observation_id": int(r["observation_id"]),
            "photo_id": int(r["photo_id"]),
            "discovery_source": str(r["discovery_source"]),
            "candidate_rows_in_taxon_cell": int(len(g)),
        })
    out = pd.DataFrame(rows).sort_values(["cell_id", "inat_taxon_id"], kind="mergesort").reset_index(drop=True)
    if len(out) != EXPECTED_PAIRS or out[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != EXPECTED_PAIRS:
        raise RuntimeError("frozen taxon-cell anchor frame is not one row per pair")
    if out["observation_id"].nunique() != len(out) or out["photo_id"].nunique() != len(out):
        raise RuntimeError("selected taxon-cell anchors are not globally unique observation/photo IDs")

    s = pd.read_csv(SPECIES_ANCHOR, usecols=["inat_taxon_id", "anchor_any_observation_id", "anchor_any_photo_id"])
    s["anchor_any_observation_id"] = pd.to_numeric(s["anchor_any_observation_id"], errors="raise").astype(int)
    s["anchor_any_photo_id"] = pd.to_numeric(s["anchor_any_photo_id"], errors="raise").astype(int)
    pair_keys = set(zip(out["observation_id"].astype(int), out["photo_id"].astype(int)))
    species_keys = set(zip(s["anchor_any_observation_id"].astype(int), s["anchor_any_photo_id"].astype(int)))
    overlap = len(pair_keys & species_keys)

    richness = out.groupby("cell_id", observed=True)["inat_taxon_id"].nunique().sort_index()
    result = {
        "analysis": "rgfca_42111_taxon_cell_anchor_step8e2",
        "status": "complete_metadata_only_taxon_cell_anchor_freeze",
        "taxon_cell_anchors": int(len(out)),
        "species": int(out["inat_taxon_id"].nunique()),
        "occupied_cells": int(out["cell_id"].nunique()),
        "species_anchor_overlap_rows": int(overlap),
        "additional_unique_anchors_beyond_species_breadth": int(len(out) - overlap),
        "cell_species_richness": {
            "min": int(richness.min()),
            "median": float(richness.median()),
            "max": int(richness.max()),
        },
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "selection_rule": "minimum SHA256 of 20260911|cell_id|inat_taxon_id|observation_id|photo_id within each frozen taxon-cell pair",
        "boundary": "This frame is for cell-wise geographic breadth. It must not replace the one-anchor-per-species estimator for global species-equal composition."
    }
    out.to_csv(OUT / "taxon_cell_anchor_frame.csv.gz", index=False, compression="gzip", lineterminator="\n")
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8E2 — taxon×cell geographic breadth anchor frame\n\n"
        f"- taxon×cell anchors: **{len(out):,}**\n"
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
