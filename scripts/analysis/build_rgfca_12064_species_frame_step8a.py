#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv"
OUTDIR = ROOT / "analysis_outputs/rgfca_12064_step8a"
TARGET_N = 12064
DEPTHS = [100, 80, 60, 50, 40, 30, 20, 10]
SEED = "20260910-step8a"


def stable_key(taxon_id: object, species: object) -> str:
    return hashlib.sha256(f"{SEED}|{taxon_id}|{species}".encode()).hexdigest()


def main() -> int:
    df = pd.read_csv(AUDIT)
    required = {"inat_taxon_id", "species", "after_observer_cap"}
    missing = required.difference(df.columns)
    if missing:
        raise RuntimeError(f"missing required columns: {sorted(missing)}; columns={list(df.columns)}")

    df["after_observer_cap"] = pd.to_numeric(df["after_observer_cap"], errors="coerce").fillna(0).astype(int)
    counts = {str(d): int((df["after_observer_cap"] >= d).sum()) for d in DEPTHS}
    selected_depth = next((d for d in DEPTHS if counts[str(d)] >= TARGET_N), None)

    result = {
        "protocol": "rgfca-12064-species-step8a-v1",
        "source_rows": int(len(df)),
        "target_species": TARGET_N,
        "candidate_depth_counts": counts,
        "selected_common_raw_photo_depth": selected_depth,
        "pixels_opened": False,
        "flower_colour_used": False,
        "selection_status": "eligible" if selected_depth is not None else "not_evaluable_no_fixed_depth_supports_12064_species",
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    if selected_depth is None:
        (OUTDIR / "capacity_result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
        return 2

    eligible = df.loc[df["after_observer_cap"] >= selected_depth].copy()
    eligible["selection_hash"] = [stable_key(t, s) for t, s in zip(eligible["inat_taxon_id"], eligible["species"])]
    eligible = eligible.sort_values(["selection_hash", "inat_taxon_id", "species"], kind="mergesort")
    frame = eligible.head(TARGET_N).copy()
    if len(frame) != TARGET_N or frame["inat_taxon_id"].nunique() != TARGET_N:
        raise RuntimeError("failed to construct exact unique 12,064-species frame")

    keep = [c for c in ["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km", "selection_hash"] if c in frame.columns]
    frame[keep].to_csv(OUTDIR / "rgfca_12064_species_frame.csv", index=False)

    result.update({
        "eligible_species_at_selected_depth": int(len(eligible)),
        "selected_species": int(len(frame)),
        "minimum_after_observer_cap_in_selected": int(frame["after_observer_cap"].min()),
        "maximum_after_observer_cap_in_selected": int(frame["after_observer_cap"].max()),
        "nested_checkpoints": [1000, 2000, 4000, 8000, 12064],
    })
    (OUTDIR / "capacity_result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
