#!/usr/bin/env python3
"""Post hoc photo-density balanced 50-km flower-colour spatial analysis.

All analyses use the immutable FCP measured source and preserve each species'
complete colour-state composition under whole-vertex permutations. Geographic
cells are photographic sampling blocks, not biological populations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    MORPH_CODE, pairwise_geo_km, load,
)

RADIUS_KM = 50.0
CELL_SIZE_DEGREES = 0.25
MIN_PHOTO_PAIRS = 30
MIN_CELL_PAIR_PHOTOS = 4
MIN_PHOTOS_PER_CELL = 3
MIN_CELL_PAIRS = 2
PERMUTATIONS = 199
COHORTS = ("discovery", "validation", "third")
MIN_COHORT_SPECIES = 30
SEED_PREFIX = "fcp-posthoc-site-pair-balance-20261017-v1"


def seed_for(cohort: str, taxon: int) -> int:
    b = f"{SEED_PREFIX}|{cohort}|{taxon}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(b).digest()[:8], "little")


def photo_cells(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if lat.shape != lon.shape or not (np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError("nonfinite or mismatched coordinates")
    if np.any(np.abs(lat) > 90) or np.any(np.abs(lon) > 180):
        raise ValueError("coordinate outside latitude/longitude domain")
    ilat = np.minimum(np.floor((lat + 90) / CELL_SIZE_DEGREES).astype(np.int64), 719)
    ilon = np.floor(((lon + 180) % 360) / CELL_SIZE_DEGREES).astype(np.int64)
    return ilat * 1440 + ilon


def analyse_sitebalanced_species(g: pd.DataFrame, cohort: str) -> tuple[dict | None, np.ndarray | None]:
    g = g.dropna(subset=["latitude", "longitude"]).sort_values("photo_id", kind="stable")
    n = len(g)
    if n < 40:
        return None, None
    lat = g.latitude.to_numpy(float)
    lon = g.longitude.to_numpy(float)
    cells = photo_cells(lat, lon)
    cell_ids, cell_counts = np.unique(cells, return_counts=True)
    supported_cells = cell_ids[cell_counts >= MIN_PHOTOS_PER_CELL]
    if len(supported_cells) < 2:
        return None, None
    labels = np.array([MORPH_CODE[str(x)] for x in g.morph], dtype=np.int8)
    u, v = np.triu_indices(n, 1)
    D_pair = float(np.mean(labels[u] != labels[v]))
    distances = pairwise_geo_km(lat, lon)[u, v]
    original_within = distances <= RADIUS_KM
    within = original_within & np.isin(cells[u], supported_cells) & np.isin(cells[v], supported_cells)
    if np.count_nonzero(within) < MIN_PHOTO_PAIRS:
        return None, None
    uu, vv = u[within], v[within]
    aa, bb = np.minimum(cells[uu], cells[vv]), np.maximum(cells[uu], cells[vv])
    keys = aa * (720 * 1440) + bb
    unique, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    kept_group = counts >= MIN_CELL_PAIR_PHOTOS
    if int(np.sum(kept_group)) < MIN_CELL_PAIRS:
        return None, None
    keep = kept_group[inverse]
    uu, vv, inverse = uu[keep], vv[keep], inverse[keep]
    kept_keys = unique[kept_group]
    left = kept_keys // (720 * 1440)
    right = kept_keys % (720 * 1440)
    if len(np.unique(np.concatenate((left, right)))) < 2:
        return None, None
    mapping = np.full(len(counts), -1, dtype=int)
    mapping[np.flatnonzero(kept_group)] = np.arange(np.sum(kept_group))
    idx = mapping[inverse]
    count_per_group = np.bincount(idx)
    edge_weight = 1.0 / (len(count_per_group) * count_per_group[idx])
    assert np.isclose(edge_weight.sum(), 1.0, rtol=0, atol=1e-12)
    observed = float(np.dot(edge_weight, (labels[uu] != labels[vv]).astype(float)))
    raw_local = float(np.mean(labels[u[original_within]] != labels[v[original_within]]))
    rng = np.random.default_rng(seed_for(cohort, int(g.inat_taxon_id.iloc[0])))
    permutations = np.stack([rng.permutation(labels) for _ in range(PERMUTATIONS)])
    local_null = (permutations[:, uu] != permutations[:, vv]).astype(float) @ edge_weight
    null_depletion = D_pair - local_null
    row = {
        "cohort": cohort,
        "inat_taxon_id": int(g.inat_taxon_id.iloc[0]),
        "species": str(g.species.iloc[0]),
        "n_classifiable": n,
        "n_geocells": int(len(supported_cells)),
        "n_all_photo_geocells": int(len(cell_ids)),
        "n_cell_pair_blocks": int(len(count_per_group)),
        "n_retained_50km_photo_pairs": int(len(uu)),
        "n_all_50km_photo_pairs": int(np.count_nonzero(original_within)),
        "n_geocell_supported_50km_photo_pairs": int(np.count_nonzero(within)),
        "max_to_median_block_depth": float(np.max(count_per_group) / np.median(count_per_group)),
        "D_specieswide": D_pair,
        "D_local_unweighted": raw_local,
        "D_local_sitebalanced": observed,
        "depletion_unweighted": D_pair - raw_local,
        "depletion_sitebalanced": D_pair - observed,
        "null_mean_depletion_sitebalanced": float(np.mean(null_depletion)),
    }
    return row, null_depletion


def summarise(rows: list[dict], nulls: list[np.ndarray]) -> dict:
    if not rows:
        return {"evaluable": False, "n_species": 0, "interpretation": "No eligible taxa after density-control gates"}
    df = pd.DataFrame(rows)
    null_means = np.vstack(nulls).mean(axis=0)
    obs = float(df.depletion_sitebalanced.mean())
    rng = np.random.default_rng(20261017)
    vals = df.depletion_sitebalanced.to_numpy(float)
    boot = vals[rng.integers(0, len(vals), size=(1999, len(vals)))].mean(axis=1)
    return {
        "evaluable": True,
        "n_species": len(rows),
        "mean_sitebalanced_depletion": obs,
        "mean_unweighted_depletion_same_species": float(df.depletion_unweighted.mean()),
        "median_sitebalanced_depletion": float(df.depletion_sitebalanced.median()),
        "positive_species_fraction": float(np.mean(vals > 0)),
        "total_retained_cell_pair_blocks": int(df.n_cell_pair_blocks.sum()),
        "median_geocells_per_species": float(df.n_geocells.median()),
        "null_mean": float(np.mean(null_means)),
        "null_q025": float(np.quantile(null_means, .025)),
        "null_q975": float(np.quantile(null_means, .975)),
        "p_upper": float((1 + np.count_nonzero(null_means >= obs)) / (PERMUTATIONS + 1)),
        "species_bootstrap_mean_95": [float(z) for z in np.quantile(boot, [.025, .975])],
    }


def analyse_cohort(path: Path, cohort: str) -> tuple[dict, pd.DataFrame]:
    d = load(path, cohort)
    rows, nulls = [], []
    for _, g in d.groupby("inat_taxon_id", sort=True):
        row, null = analyse_sitebalanced_species(g, cohort)
        if row is not None:
            rows.append(row)
            nulls.append(null)
    return summarise(rows, nulls), pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    for c in COHORTS:
        parser.add_argument("--" + c, type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    a = parser.parse_args()
    outdir = a.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    stats, frames = {}, []
    for c in COHORTS:
        stat, detail = analyse_cohort(getattr(a, c), c)
        stats[c] = stat
        frames.append(detail)
    result = {
        "schema": "fcp_site_pair_balanced_50km_posthoc_v1",
        "role": "posthoc_spatial_photo_density_sensitivity",
        "original_frozen_results_changed": False,
        "independent_source_replication": False,
        "design": {
            "primary_radius_km": RADIUS_KM,
            "cell_size_degrees": CELL_SIZE_DEGREES,
            "minimum_original_local_photo_pairs": MIN_PHOTO_PAIRS,
            "minimum_cell_pair_photo_pairs": MIN_CELL_PAIR_PHOTOS,
            "minimum_photos_per_geocell": MIN_PHOTOS_PER_CELL,
            "minimum_distinct_cell_pair_blocks": MIN_CELL_PAIRS,
            "permutations": PERMUTATIONS,
            "minimum_cohort_species_for_crosscohort_support": MIN_COHORT_SPECIES,
            "weights": "Each retained unordered geographic cell-pair block contributes 1/K, regardless of photographed pair count; exact species morph multiset remains fixed by vertex permutation",
        },
        "cohorts": stats,
        "positive_all_three": bool(all(s.get("evaluable") and s["mean_sitebalanced_depletion"] > 0 for s in stats.values())),
        "p_under_0_05_all_three": bool(all(s.get("evaluable") and s["n_species"] >= MIN_COHORT_SPECIES and s["mean_sitebalanced_depletion"] > 0 and s["p_upper"] < .05 for s in stats.values())),
        "hard_nonclaims": [
            "post hoc study of photographed phenotypes, not genetic FCP",
            "coarse site cells are not actual biological populations",
            "not adaptation, selection, pigment fitness or seasonal climate",
            "not independent photo-provider replication",
        ],
    }
    pd.concat(frames, ignore_index=True).to_csv(outdir / "species_sitebalanced_50km.csv", index=False)
    (outdir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print("SITE_BALANCED_50KM_COMPLETE", json.dumps({c: stats[c] for c in COHORTS}, sort_keys=True))


if __name__ == "__main__":
    main()
