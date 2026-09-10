#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CENSUS = ROOT / "results/rgfca_42111_breadth_depth_step8b_20260911/species_capacity_census.csv.gz"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
OUT = ROOT / "results/rgfca_42111_identifiability_ladder_step8b3_20260911"
N = 42111


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(CENSUS)
    if len(x) != N or x["inat_taxon_id"].nunique() != N:
        raise RuntimeError("census denominator drift")
    links = []
    for path in (V1, V2):
        q = pd.read_csv(path, usecols=["inat_taxon_id", "cell_id"])
        q["inat_taxon_id"] = pd.to_numeric(q["inat_taxon_id"], errors="raise").astype(int)
        q["cell_id"] = pd.to_numeric(q["cell_id"], errors="raise").astype(int)
        links.append(q)
    cell = pd.concat(links, ignore_index=True).drop_duplicates(["inat_taxon_id", "cell_id"])
    ncell = cell.groupby("inat_taxon_id", observed=True)["cell_id"].nunique().rename("n_discovery_cells")
    x = x.merge(ncell, on="inat_taxon_id", how="left", validate="one_to_one")
    x["n_discovery_cells"] = x["n_discovery_cells"].fillna(0).astype(int)
    cap = x["after_observer_cap"].astype(int)
    span = pd.to_numeric(x["maximum_span_km"], errors="coerce").fillna(0.0)

    masks = {
        "breadth_anchor_available": pd.Series(True, index=x.index),
        "at_least_2_raw_photos": cap >= 2,
        "at_least_5_raw_photos": cap >= 5,
        "at_least_10_raw_photos": cap >= 10,
        "at_least_20_raw_photos": cap >= 20,
        "depth20_and_ge2_discovery_cells": (cap >= 20) & (x["n_discovery_cells"] >= 2),
        "depth20_and_ge5_discovery_cells": (cap >= 20) & (x["n_discovery_cells"] >= 5),
        "depth20_and_span_ge100km": (cap >= 20) & (span >= 100.0),
        "depth20_ge5cells_span_ge100km": (cap >= 20) & (x["n_discovery_cells"] >= 5) & (span >= 100.0),
    }
    rows = []
    for name, mask in masks.items():
        n = int(mask.sum())
        rows.append({"tier": name, "species": n, "fraction_of_42111": n / float(N)})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "identifiability_ladder.csv", index=False, lineterminator="\n")

    result = {
        "analysis": "rgfca_42111_identifiability_ladder_step8b3",
        "species_universe": N,
        "tiers": {r["tier"]: {"species": int(r["species"]), "fraction": float(r["fraction_of_42111"])} for r in rows},
        "pixels_opened": False,
        "flower_colour_used": False,
        "interpretation": {
            "breadth_anchor_available": "Supports one-photo species-equal breadth summaries only; not species polymorphism.",
            "at_least_2_raw_photos": "Only a repeated-observation opportunity bound; two photos are not a reliable D estimate.",
            "at_least_5_raw_photos": "Coarse within-species colour-diversity opportunity; uncertainty must be explicit.",
            "at_least_10_raw_photos": "Moderate within-species colour-diversity opportunity; still not a C*/S* guarantee.",
            "at_least_20_raw_photos": "Step-8 tiered depth ceiling for new large-scale polymorphism summaries.",
            "depth20_and_ge2_discovery_cells": "Upper-bound geometry for spatial comparison, not proof of estimable C*/S*.",
            "depth20_and_ge5_discovery_cells": "Stronger pre-colour spatial-support geometry; still requires measured/classifiable local support.",
        },
        "boundary": "These are pre-colour opportunity/identifiability counts, not positive biological findings."
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = ["# RGFCA Step 8B3 — 42,111-species identifiability ladder", ""]
    for r in rows:
        lines.append(f"- {r['tier']}: **{int(r['species']):,}** ({100*float(r['fraction_of_42111']):.1f}%)")
    lines += ["", "No image pixel or flower-colour outcome was opened.", "These are opportunity bounds, not biological positives."]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
