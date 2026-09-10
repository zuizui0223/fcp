#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPECIES = ROOT / "data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv"
CAPACITY = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv"
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
DISC_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_v2_manifest_v1.json"
CAP_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json"
OUT = ROOT / "results/rgfca_42111_breadth_depth_step8b_20260911"

N_SPECIES = 42111
DEPTHS = [1, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100]
BREADTH_N = [500, 1000, 2000, 5000, 10000, 20000, 30000, 42111]
N_REP = 200
SEED0 = 20260911
N_CELLS = 162


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def availability_bin(n: int) -> str:
    if n <= 0:
        return "0_structural_no_current_eligible_photo"
    if n == 1:
        return "1"
    if n <= 4:
        return "2_4"
    if n <= 9:
        return "5_9"
    if n <= 19:
        return "10_19"
    if n <= 29:
        return "20_29"
    if n <= 39:
        return "30_39"
    if n <= 49:
        return "40_49"
    if n <= 59:
        return "50_59"
    if n <= 79:
        return "60_79"
    if n <= 99:
        return "80_99"
    if n <= 149:
        return "100_149"
    return "150_plus"


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def load_and_validate() -> tuple[pd.DataFrame, pd.DataFrame]:
    dm = json.loads(DISC_MANIFEST.read_text(encoding="utf-8"))
    cm = json.loads(CAP_MANIFEST.read_text(encoding="utf-8"))
    if int(dm["combined"]["species"]) != N_SPECIES:
        raise RuntimeError("discovery universe is not 42,111 species")
    if dm.get("candidate_image_pixels_opened") is not False or dm.get("flower_colour_used") is not False:
        raise RuntimeError("discovery frame opened forbidden outcomes")
    if int(cm["discovered_species_scanned"]) != N_SPECIES:
        raise RuntimeError("capacity scan did not cover all 42,111 species")
    if int(cm.get("remaining_request_errors", 1)) != 0:
        raise RuntimeError("capacity recovery still has unresolved request errors")
    if cm.get("candidate_image_pixels_opened") is not False or cm.get("flower_colour_used") is not False:
        raise RuntimeError("capacity scan opened forbidden outcomes")

    sp = pd.read_csv(SPECIES)
    cap = pd.read_csv(CAPACITY)
    if len(sp) != N_SPECIES or sp["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError(f"species frame mismatch: {len(sp)} rows")
    if len(cap) != N_SPECIES or cap["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError(f"capacity frame mismatch: {len(cap)} rows")

    keep = ["inat_taxon_id", "after_observer_cap", "maximum_span_km", "request_error"]
    missing = [c for c in keep if c not in cap.columns]
    if missing:
        raise RuntimeError(f"capacity columns missing: {missing}")
    if cap["request_error"].fillna("").astype(str).str.len().gt(0).any():
        raise RuntimeError("capacity v3 contains non-empty request errors")

    census = sp.merge(cap[keep], on="inat_taxon_id", how="left", validate="one_to_one")
    if census["after_observer_cap"].isna().any():
        raise RuntimeError("capacity merge left missing species")
    census["after_observer_cap"] = pd.to_numeric(census["after_observer_cap"], errors="raise").astype(int)
    census["maximum_span_km"] = pd.to_numeric(census["maximum_span_km"], errors="coerce")
    census["availability_bin"] = census["after_observer_cap"].map(availability_bin)
    for d in DEPTHS:
        census[f"eligible_ge_{d}"] = census["after_observer_cap"] >= d
    census["breadth_hash"] = [
        hashlib.sha256(f"{SEED0}|{int(t)}|{s}".encode("utf-8")).hexdigest()
        for t, s in zip(census["inat_taxon_id"], census["species"])
    ]
    census = census.sort_values(["breadth_hash", "inat_taxon_id"], kind="mergesort").reset_index(drop=True)
    census["breadth_rank"] = np.arange(1, len(census) + 1)

    links = []
    for path in (V1, V2):
        q = pd.read_csv(path, usecols=["inat_taxon_id", "cell_id"])
        links.append(q)
    lk = pd.concat(links, ignore_index=True).dropna()
    lk["inat_taxon_id"] = lk["inat_taxon_id"].astype(int)
    lk["cell_id"] = lk["cell_id"].astype(int)
    lk = lk.drop_duplicates(["inat_taxon_id", "cell_id"], keep="first")
    universe = set(census["inat_taxon_id"].astype(int))
    lk = lk[lk["inat_taxon_id"].isin(universe)].copy()
    if lk["inat_taxon_id"].nunique() != N_SPECIES:
        raise RuntimeError("not every discovered species has a retained discovery taxon-cell link")
    if lk["cell_id"].min() < 0 or lk["cell_id"].max() >= N_CELLS:
        raise RuntimeError("cell ids outside frozen 162-cell grid")
    return census, lk


def breadth_saturation(census: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    taxon_to_idx = {int(t): i for i, t in enumerate(census["inat_taxon_id"].astype(int))}
    link_species_idx = links["inat_taxon_id"].map(taxon_to_idx).to_numpy(dtype=np.int32)
    link_cell = links["cell_id"].to_numpy(dtype=np.int16)
    full_counts = np.bincount(link_cell, minlength=N_CELLS).astype(float)
    full_occ = int(np.sum(full_counts > 0))
    full_links = int(len(links))
    if full_occ <= 0 or full_links <= 0:
        raise RuntimeError("empty full-frame geography")

    rows: list[dict[str, object]] = []
    nsp = len(census)
    for rep in range(N_REP):
        rng = np.random.default_rng(SEED0 + rep)
        perm = rng.permutation(nsp)
        rank = np.empty(nsp, dtype=np.int32)
        rank[perm] = np.arange(nsp, dtype=np.int32)
        link_rank = rank[link_species_idx]
        for n in BREADTH_N:
            mask = link_rank < n
            counts = np.bincount(link_cell[mask], minlength=N_CELLS).astype(float)
            rows.append({
                "replicate": rep + 1,
                "seed": SEED0 + rep,
                "species_n": n,
                "occupied_cells": int(np.sum(counts > 0)),
                "occupied_cell_retention": float(np.sum(counts > 0) / full_occ),
                "cell_species_count_pearson": pearson(counts, full_counts),
                "taxon_cell_link_retention": float(np.sum(mask) / full_links),
            })
    return pd.DataFrame(rows)


def qdict(x: pd.Series) -> dict[str, float | None]:
    z = pd.to_numeric(x, errors="coerce").dropna().to_numpy(float)
    if not len(z):
        return {"q10": None, "median": None, "q90": None}
    return {
        "q10": float(np.quantile(z, 0.10)),
        "median": float(np.quantile(z, 0.50)),
        "q90": float(np.quantile(z, 0.90)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    census, links = load_and_validate()
    sat = breadth_saturation(census, links)

    census_out = OUT / "species_capacity_census.csv.gz"
    sat_out = OUT / "breadth_saturation_200x.csv"
    census.to_csv(census_out, index=False, compression="gzip")
    sat.to_csv(sat_out, index=False)

    threshold_counts = {str(d): int((census["after_observer_cap"] >= d).sum()) for d in DEPTHS}
    bin_counts = {str(k): int(v) for k, v in census["availability_bin"].value_counts().sort_index().items()}
    source_counts = {
        "v1": int(census["source_v1"].astype(bool).sum()),
        "v2": int(census["source_v2"].astype(bool).sum()),
        "both": int((census["source_v1"].astype(bool) & census["source_v2"].astype(bool)).sum()),
        "v1_only": int((census["source_v1"].astype(bool) & ~census["source_v2"].astype(bool)).sum()),
        "v2_only": int((~census["source_v1"].astype(bool) & census["source_v2"].astype(bool)).sum()),
    }

    sat_summary: dict[str, dict[str, float]] = {}
    for n, g in sat.groupby("species_n", sort=True):
        sat_summary[str(int(n))] = {
            "median_occupied_cell_retention": float(g["occupied_cell_retention"].median()),
            "q05_occupied_cell_retention": float(g["occupied_cell_retention"].quantile(.05)),
            "median_cell_species_count_pearson": float(g["cell_species_count_pearson"].median()),
            "q05_cell_species_count_pearson": float(g["cell_species_count_pearson"].quantile(.05)),
            "median_taxon_cell_link_retention": float(g["taxon_cell_link_retention"].median()),
        }

    result = {
        "analysis": "rgfca_42111_breadth_depth_step8b",
        "status": "complete_metadata_only_breadth_depth_census",
        "species_universe": N_SPECIES,
        "pixels_opened": False,
        "flower_colour_used": False,
        "threshold_counts_after_observer_cap": threshold_counts,
        "availability_bin_counts": bin_counts,
        "source_membership": source_counts,
        "maximum_span_km_quantiles_all_species": qdict(census["maximum_span_km"]),
        "maximum_span_km_quantiles_ge1": qdict(census.loc[census["after_observer_cap"] >= 1, "maximum_span_km"]),
        "discovery_taxon_cell_links": int(len(links)),
        "occupied_equal_area_cells_full": int(links["cell_id"].nunique()),
        "breadth_saturation": sat_summary,
        "lineage": {
            "species_frame_sha256": sha256_file(SPECIES),
            "capacity_audit_sha256": sha256_file(CAPACITY),
            "v1_index_sha256": sha256_file(V1),
            "v2_index_sha256": sha256_file(V2),
            "census_sha256": sha256_file(census_out),
            "saturation_sha256": sha256_file(sat_out),
        },
        "claim_boundary": "Metadata-only sampling/capacity result. No flower-colour prevalence, polymorphism, C*/S*, or biological mechanism is inferred here.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# RGFCA Step 8B — 42,111-species breadth–depth census",
        "",
        f"- full species universe: **{N_SPECIES:,}**",
        f"- discovery taxon–cell links: **{len(links):,}**",
        f"- occupied equal-area cells: **{links['cell_id'].nunique()} / {N_CELLS}**",
        "- image pixels opened: **false**",
        "- flower colour used: **false**",
        "",
        "## Available raw-photo depth after observer cap",
        "",
    ]
    for d in DEPTHS:
        lines.append(f"- >= {d}: **{threshold_counts[str(d)]:,} species**")
    lines += ["", "## Repeated geographic breadth saturation", ""]
    for n in BREADTH_N:
        z = sat_summary[str(n)]
        lines.append(
            f"- N={n:,}: occupied-cell retention median **{z['median_occupied_cell_retention']:.3f}** "
            f"(5th pct **{z['q05_occupied_cell_retention']:.3f}**); cell-count r median **{z['median_cell_species_count_pearson']:.3f}** "
            f"(5th pct **{z['q05_cell_species_count_pearson']:.3f}**); taxon-cell-link retention median **{z['median_taxon_cell_link_retention']:.3f}**"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "All 42,111 discovered species remain in the sampling universe. Species with zero current eligible capacity are retained as structural missingness. This result diagnoses sampling breadth and depth only; it does not infer flower-colour frequencies or polymorphism.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
