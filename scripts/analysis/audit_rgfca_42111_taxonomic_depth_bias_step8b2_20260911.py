#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CENSUS = ROOT / "results/rgfca_42111_breadth_depth_step8b_20260911/species_capacity_census.csv.gz"
COMBINED = ROOT / "data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
OUT = ROOT / "results/rgfca_42111_taxonomic_depth_bias_step8b2_20260911"
N = 42111


def genus_of(name: object) -> str:
    parts = str(name or "").strip().split()
    return parts[0] if parts else ""


def effective_number(labels: pd.Series) -> float:
    p = labels.value_counts(normalize=True)
    return float(1.0 / np.square(p.to_numpy(float)).sum()) if len(p) else float("nan")


def jsd_from_counts(a: pd.Series, b: pd.Series) -> float:
    keys = sorted(set(a.index.astype(str)) | set(b.index.astype(str)))
    pa = np.asarray([float(a.get(k, 0.0)) for k in keys], dtype=float)
    pb = np.asarray([float(b.get(k, 0.0)) for k in keys], dtype=float)
    pa = pa / pa.sum()
    pb = pb / pb.sum()
    m = 0.5 * (pa + pb)
    def kl(p: np.ndarray, q: np.ndarray) -> float:
        keep = p > 0
        return float(np.sum(p[keep] * np.log2(p[keep] / q[keep])))
    return 0.5 * kl(pa, m) + 0.5 * kl(pb, m)


def qstats(s: pd.Series) -> dict[str, float]:
    x = pd.to_numeric(s, errors="coerce").dropna().to_numpy(float)
    if len(x) == 0:
        return {"min": float("nan"), "median": float("nan"), "p90": float("nan"), "max": float("nan")}
    return {
        "min": float(np.min(x)),
        "median": float(np.median(x)),
        "p90": float(np.quantile(x, 0.9)),
        "max": float(np.max(x)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(CENSUS)
    if len(x) != N or x["inat_taxon_id"].nunique() != N:
        raise RuntimeError("census is not the exact 42,111-species universe")
    x["genus_label"] = x["species"].map(genus_of)
    if x["genus_label"].eq("").any():
        raise RuntimeError("one or more species names lack a first-token genus label")

    combined = pd.read_csv(COMBINED)
    if len(combined) != N or combined["inat_taxon_id"].nunique() != N:
        raise RuntimeError("combined V1/V2 species frame is not 42,111 unique taxa")
    combined["source_v1"] = combined["source_v1"].astype(str).str.lower().isin(["true", "1"])
    combined["source_v2"] = combined["source_v2"].astype(str).str.lower().isin(["true", "1"])
    source_counts = {
        "v1_only": int((combined["source_v1"] & ~combined["source_v2"]).sum()),
        "v2_only": int((~combined["source_v1"] & combined["source_v2"]).sum()),
        "both": int((combined["source_v1"] & combined["source_v2"]).sum()),
    }

    full_genus_counts = x["genus_label"].value_counts().sort_index()
    singleton_genera = int((full_genus_counts == 1).sum())
    depth_rows = []
    depth_summary: dict[str, dict[str, object]] = {}
    for depth in [1, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100]:
        sub = x.loc[x["after_observer_cap"].astype(int) >= depth].copy()
        counts = sub["genus_label"].value_counts().sort_index()
        depth_rows.append({
            "minimum_raw_photos": depth,
            "species": int(len(sub)),
            "genera": int(sub["genus_label"].nunique()),
            "genus_retention": float(sub["genus_label"].nunique() / x["genus_label"].nunique()),
            "effective_genera_simpson": effective_number(sub["genus_label"]),
            "genus_composition_jsd_vs_full": jsd_from_counts(full_genus_counts, counts),
        })
        depth_summary[str(depth)] = dict(depth_rows[-1])
    depth_table = pd.DataFrame(depth_rows)
    depth_table.to_csv(OUT / "depth_taxonomic_retention.csv", index=False, lineterminator="\n")

    links = []
    for source, path in (("v1", V1), ("v2", V2)):
        q = pd.read_csv(path, usecols=["cell_id", "inat_taxon_id"])
        q["source"] = source
        links.append(q)
    cell = pd.concat(links, ignore_index=True)
    cell["cell_id"] = pd.to_numeric(cell["cell_id"], errors="raise").astype(int)
    cell["inat_taxon_id"] = pd.to_numeric(cell["inat_taxon_id"], errors="raise").astype(int)
    cell = cell.drop_duplicates(["inat_taxon_id", "cell_id"])
    occupied_per_species = cell.groupby("inat_taxon_id", observed=True)["cell_id"].nunique()
    richness = cell.groupby("cell_id", observed=True)["inat_taxon_id"].nunique().sort_values(ascending=False)
    total_links = int(len(cell))
    top10_share = float(richness.head(10).sum() / total_links) if total_links else float("nan")

    result = {
        "analysis": "rgfca_42111_taxonomic_depth_bias_step8b2",
        "species_universe": N,
        "genus_labels_first_token": int(x["genus_label"].nunique()),
        "singleton_genus_labels": singleton_genera,
        "singleton_genus_fraction": float(singleton_genera / x["genus_label"].nunique()),
        "effective_genera_simpson_full": effective_number(x["genus_label"]),
        "species_per_genus": qstats(full_genus_counts),
        "v1_v2_source_species": source_counts,
        "taxon_cell_links": total_links,
        "occupied_equal_area_cells": int(richness.size),
        "cell_species_richness": qstats(richness),
        "top10_cells_taxon_cell_link_share": top10_share,
        "species_occupied_cell_counts": qstats(occupied_per_species),
        "species_with_one_discovery_cell": int((occupied_per_species == 1).sum()),
        "species_with_ge_2_discovery_cells": int((occupied_per_species >= 2).sum()),
        "species_with_ge_5_discovery_cells": int((occupied_per_species >= 5).sum()),
        "species_with_ge_10_discovery_cells": int((occupied_per_species >= 10).sum()),
        "depth20": depth_summary["20"],
        "pixels_opened": False,
        "flower_colour_used": False,
        "boundary": "Genus labels are the first token of the frozen scientific name and are descriptive labels, not an independently resolved phylogeny or family-level audit."
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    d20 = depth_summary["20"]
    (OUT / "RESULT.md").write_text(
        "# RGFCA Step 8B2 — 42,111-species taxonomic/depth audit\n\n"
        f"- full species: **{N:,}**\n"
        f"- first-token genus labels: **{result['genus_labels_first_token']:,}**\n"
        f"- singleton genus labels: **{singleton_genera:,}**\n"
        f"- V1 only / V2 only / both: **{source_counts['v1_only']:,} / {source_counts['v2_only']:,} / {source_counts['both']:,}**\n"
        f"- unique taxon-cell links: **{total_links:,}**\n"
        f"- occupied equal-area cells: **{int(richness.size)}**\n"
        f"- species seen in >=2 / >=5 / >=10 cells: **{result['species_with_ge_2_discovery_cells']:,} / {result['species_with_ge_5_discovery_cells']:,} / {result['species_with_ge_10_discovery_cells']:,}**\n\n"
        "## Depth >=20 subset\n\n"
        f"- species: **{int(d20['species']):,}**\n"
        f"- genus labels: **{int(d20['genera']):,}**\n"
        f"- genus retention vs full: **{float(d20['genus_retention']):.3f}**\n"
        f"- genus-composition JSD vs full: **{float(d20['genus_composition_jsd_vs_full']):.4f}**\n\n"
        "No image pixel or flower-colour outcome was opened.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
