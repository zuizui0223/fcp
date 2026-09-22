#!/usr/bin/env python3
"""Outcome-blind species selection for FCP v2 measurement validity.

Inputs are identity/capacity manifests only. This script must never read image
pixels, colour measurements, morphs, D, q_white/W, spatial outcomes, H3
predictors, family, climate or pollinator variables.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

P_SALT = "FCP_V2_PANEL_P_20260923"
N_SALT = "FCP_V2_PANEL_N_20260923"
P_TARGET = 200
N_TARGET = 200
P_MIN_CAPACITY = 200
N_MIN_CAPACITY = 100

ALLOWED = ("inat_taxon_id", "species", "after_observer_cap")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_identity(frame: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in ALLOWED if c not in frame.columns]
    if missing:
        raise ValueError(f"identity/capacity input missing columns: {missing}")
    x = frame.loc[:, ALLOWED].copy()
    x["inat_taxon_id"] = pd.to_numeric(x["inat_taxon_id"], errors="raise").astype("int64")
    x["species"] = x["species"].astype(str).str.strip()
    x["after_observer_cap"] = pd.to_numeric(
        x["after_observer_cap"], errors="raise"
    ).astype("int64")
    if x["inat_taxon_id"].duplicated().any():
        raise ValueError("taxon identity is not unique")
    if x["species"].duplicated().any():
        raise ValueError("species identity is not unique")
    if x["species"].eq("").any():
        raise ValueError("blank species identity")
    return x


def hash_rank(frame: pd.DataFrame, *, salt: str) -> pd.DataFrame:
    x = frame.copy()
    x["selection_hash"] = [
        hashlib.sha256(f"{salt}|{taxon}|{species}".encode("utf-8")).hexdigest()
        for taxon, species in zip(x["inat_taxon_id"], x["species"])
    ]
    x = x.sort_values(
        ["selection_hash", "inat_taxon_id"],
        kind="mergesort",
    ).reset_index(drop=True)
    x["selection_rank"] = range(1, len(x) + 1)
    return x


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-high-depth-audit", type=Path, required=True)
    p.add_argument("--p100-pool", type=Path, required=True)
    p.add_argument("--p500-selected", type=Path, required=True)
    p.add_argument("--third-selected", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    legacy_raw = pd.read_csv(args.legacy_high_depth_audit)
    p100_raw = pd.read_csv(args.p100_pool)
    p500_raw = pd.read_csv(args.p500_selected)
    third_raw = pd.read_csv(args.third_selected, sep="\t")

    # Explicit firewall: only identity/capacity columns are read into selection.
    legacy = canonical_identity(legacy_raw)
    p100 = canonical_identity(p100_raw)
    p500 = canonical_identity(p500_raw)

    if not {"inat_taxon_id", "species"}.issubset(third_raw.columns):
        raise ValueError("third-cohort identity manifest schema drift")
    third = third_raw[["inat_taxon_id", "species"]].copy()
    third["inat_taxon_id"] = pd.to_numeric(
        third["inat_taxon_id"], errors="raise"
    ).astype("int64")
    third["species"] = third["species"].astype(str).str.strip()
    if len(third) != 500:
        raise ValueError(f"third selected denominator drift: {len(third)}")
    if third["inat_taxon_id"].nunique() != 500 or third["species"].nunique() != 500:
        raise ValueError("third selected identities are not unique")

    if len(legacy) != 1000:
        raise ValueError(f"legacy high-depth denominator drift: {len(legacy)}")
    if len(p100) != 3730:
        raise ValueError(f"P100 denominator drift: {len(p100)}")
    if len(p500) != 500:
        raise ValueError(f"P500 denominator drift: {len(p500)}")

    p_candidates = legacy.loc[
        legacy["after_observer_cap"] >= P_MIN_CAPACITY
    ].copy()
    if len(p_candidates) < P_TARGET:
        raise RuntimeError("Panel P has insufficient metadata-only fresh-photo capacity")
    p_ranked = hash_rank(p_candidates, salt=P_SALT)
    p_selected = p_ranked.iloc[:P_TARGET].copy()
    p_selected.insert(0, "panel", "P")

    p500_taxa = set(p500["inat_taxon_id"].astype(int))
    third_taxa = set(third["inat_taxon_id"].astype(int))
    n_candidates = p100.loc[
        (~p100["inat_taxon_id"].isin(p500_taxa))
        & (~p100["inat_taxon_id"].isin(third_taxa))
        & (p100["after_observer_cap"] >= N_MIN_CAPACITY)
    ].copy()

    if len(n_candidates) != 2730:
        raise RuntimeError(
            f"Panel N expected 2730 unused P100 species, observed {len(n_candidates)}"
        )
    n_ranked = hash_rank(n_candidates, salt=N_SALT)
    n_selected = n_ranked.iloc[:N_TARGET].copy()
    n_selected.insert(0, "panel", "N")

    overlap_taxa = set(p_selected["inat_taxon_id"]) & set(n_selected["inat_taxon_id"])
    overlap_species = set(p_selected["species"]) & set(n_selected["species"])
    if overlap_taxa or overlap_species:
        raise RuntimeError("Panel P/N overlap detected")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    p_ranked.to_csv(out / "panel_p_candidate_pool.csv", index=False, lineterminator="\n")
    n_ranked.to_csv(out / "panel_n_candidate_pool.csv", index=False, lineterminator="\n")
    p_selected.to_csv(out / "panel_p_selected_species.csv", index=False, lineterminator="\n")
    n_selected.to_csv(out / "panel_n_selected_species.csv", index=False, lineterminator="\n")
    combined = pd.concat([p_selected, n_selected], ignore_index=True)
    combined.to_csv(out / "selected_species_manifest.tsv", sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema": "fcp_v2_species_selection_v1",
        "date_jst": "2026-09-23",
        "status": "COMPLETE_OUTCOME_BLIND_METADATA_ONLY_SELECTION",
        "panel_p": {
            "candidate_species": int(len(p_ranked)),
            "minimum_after_observer_cap": P_MIN_CAPACITY,
            "selected_species": P_TARGET,
            "salt": P_SALT,
            "role": "paired_species_fresh_image_transport",
        },
        "panel_n": {
            "candidate_species": int(len(n_ranked)),
            "minimum_after_observer_cap": N_MIN_CAPACITY,
            "selected_species": N_TARGET,
            "salt": N_SALT,
            "role": "novel_species_measurement_generalization",
        },
        "selected_species_total": int(len(combined)),
        "panel_overlap_taxa": 0,
        "panel_overlap_species": 0,
        "selection_columns_read": list(ALLOWED),
        "third_manifest_columns_read": ["inat_taxon_id", "species"],
        "flower_colour_outcomes_read": False,
        "image_pixels_opened": False,
        "D_read": False,
        "H2_read": False,
        "spatial_outcomes_read": False,
        "H3_predictors_read": False,
        "family_or_climate_used": False,
        "replacement_after_pixel_opening_allowed": False,
        "input_sha256": {
            "legacy_high_depth_audit": file_sha256(args.legacy_high_depth_audit),
            "p100_pool": file_sha256(args.p100_pool),
            "p500_selected": file_sha256(args.p500_selected),
            "third_selected": file_sha256(args.third_selected),
        },
        "output_sha256": {
            name: file_sha256(out / name)
            for name in (
                "panel_p_candidate_pool.csv",
                "panel_n_candidate_pool.csv",
                "panel_p_selected_species.csv",
                "panel_n_selected_species.csv",
                "selected_species_manifest.tsv",
            )
        },
    }
    (out / "result.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
