#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SELECTION_SCRIPT = ROOT / "scripts" / "analysis" / "select_polymorphism_h2_prospective_u100_20260913.py"
SELECTION_OUT = ROOT / "results" / "polymorphism_h2_prospective_u100_selection_20260913"
P500 = SELECTION_OUT / "p500_frozen_selection.csv"
OUT = ROOT / "results" / "polymorphism_h2_p500_acquisition_preflight_20260913"
EXPECTED_P500_SHA = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"

EXCLUSION_SOURCES = [
    ROOT / "data" / "frozen" / "random_photo_first_h9_exclusion_ledger_v1.csv",
    ROOT / "data" / "frozen" / "random_photo_first_h9_fresh_metadata_v1.csv",
    ROOT / "results" / "rgfca_42111_breadth_measurement_step8f_20260911" / "rgfca_42111_species_breadth_measured.csv.gz",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    # Recreate the frozen P500 from metadata only and require exact byte identity.
    subprocess.run(["python", str(SELECTION_SCRIPT)], cwd=ROOT, check=True)
    observed_p500_sha = sha256_file(P500)
    if observed_p500_sha != EXPECTED_P500_SHA:
        raise RuntimeError(f"P500 selection drift: {observed_p500_sha}")
    p500 = pd.read_csv(P500, usecols=["inat_taxon_id", "species", "prospective_rank", "selection_hash"])
    if len(p500) != 500 or p500["inat_taxon_id"].nunique() != 500:
        raise RuntimeError("P500 identity fingerprint drift")
    p500_species = set(p500["species"].fillna("").astype(str).str.strip())

    source_audits = []
    exclusion_obs: set[int] = set()
    exclusion_photo: set[int] = set()
    breadth_p500_species: set[str] = set()

    for path in EXCLUSION_SOURCES:
        if not path.exists():
            raise RuntimeError(f"missing exclusion source: {path}")
        header = pd.read_csv(path, nrows=0).columns.tolist()
        required = {"observation_id", "photo_id"}
        if not required.issubset(header):
            raise RuntimeError(f"exclusion source lacks observation/photo ids: {path}; columns={header}")
        usecols = ["observation_id", "photo_id"] + (["species"] if "species" in header else [])
        frame = pd.read_csv(path, usecols=usecols, low_memory=False)
        obs = pd.to_numeric(frame["observation_id"], errors="coerce").dropna().astype("int64")
        photo = pd.to_numeric(frame["photo_id"], errors="coerce").dropna().astype("int64")
        exclusion_obs.update(obs.tolist())
        exclusion_photo.update(photo.tolist())
        overlap_n = None
        if "species" in frame.columns:
            species = set(frame["species"].fillna("").astype(str).str.strip()) - {""}
            overlap = species & p500_species
            overlap_n = len(overlap)
            if "breadth_measurement_step8f" in str(path):
                breadth_p500_species |= overlap
        source_audits.append({
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256_file(path),
            "rows": int(len(frame)),
            "header_has_species": "species" in header,
            "P500_species_overlap": overlap_n,
            "unique_observation_ids": int(obs.nunique()),
            "unique_photo_ids": int(photo.nunique()),
            "colour_columns_read": False,
        })

    result = {
        "analysis": "polymorphism_h2_p500_acquisition_preflight",
        "date_jst": "2026-09-13",
        "status": "pass_metadata_only_acquisition_preflight",
        "P500": {
            "species": 500,
            "sha256": observed_p500_sha,
            "colour_outcomes_opened": False,
        },
        "exclusion_sources": source_audits,
        "combined_exclusion_ids": {
            "unique_observation_ids": int(len(exclusion_obs)),
            "unique_photo_ids": int(len(exclusion_photo)),
            "P500_species_with_prior_breadth_anchor": int(len(breadth_p500_species)),
        },
        "outcome_firewall": {
            "only_identity_capacity_and_observation_photo_ids_read": True,
            "morph_read": False,
            "palette_read": False,
            "D_read": False,
            "H2_W_read": False,
            "coordinates_read": False,
        },
        "next_gate": "candidate_metadata_acquisition_may_run; image pixels remain closed",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
