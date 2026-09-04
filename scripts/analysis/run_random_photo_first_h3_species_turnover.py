#!/usr/bin/env python3
"""Exploratory post-H1/H2 decomposition of colour persistence by species turnover.

This analysis was introduced after the confirmatory H1/H2 outcomes were known.
It therefore cannot rescue, replace, or reinterpret the frozen H1/H2 decisions.
It uses only already-frozen inputs and the exact 999 H1 null persistence maps.

Question: do edges with larger species-composition turnover also show larger
observed colour-boundary persistence, and is that alignment stronger than expected
under the H1 null that already preserves species geography and each species' morph
marginals?

Species composition is outcome-blind: within each frozen 18x9 cell, every species
is weighted by min(candidate_photo_count, 2), matching the H1 species cap while
using all 20,845 frozen candidate records regardless of colour-classification
success. Edge turnover is Jensen-Shannon divergence in bits on those normalized
species weights.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = ROOT / "data/frozen/random_photo_first_candidate_pool_v1.csv"
CANDIDATE_MANIFEST = ROOT / "docs/supporting/random_photo_first_candidate_pool_manifest_v1.json"
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"
H2_RESULT = ROOT / "docs/supporting/random_photo_first_h2_result_v1.json"
H1_EDGES = ROOT / "data/derived/random_photo_first_h1_edges_v1.csv"
H1_NULL_MAPS = ROOT / "data/derived/random_photo_first_h1_null_maps_v1.npz"
OUT_EDGE = ROOT / "data/derived/random_photo_first_h3_species_turnover_edges_v1.csv"
OUT_NULL = ROOT / "data/derived/random_photo_first_h3_species_turnover_null_correlations_v1.csv"
OUT_JSON = ROOT / "docs/supporting/random_photo_first_h3_species_turnover_v1.json"

SPECIES_CAP = 2


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sparse_jsd_bits(
    a: Mapping[str, float],
    b: Mapping[str, float],
) -> float:
    total_a = float(sum(a.values()))
    total_b = float(sum(b.values()))
    if total_a <= 0.0 or total_b <= 0.0:
        return float("nan")
    out = 0.0
    for species in set(a) | set(b):
        p = float(a.get(species, 0.0)) / total_a
        q = float(b.get(species, 0.0)) / total_b
        m = 0.5 * (p + q)
        if p > 0.0:
            out += 0.5 * p * np.log2(p / m)
        if q > 0.0:
            out += 0.5 * q * np.log2(q / m)
    return float(np.clip(out, 0.0, 1.0))


def weighted_pearson(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x = x[valid]
    y = y[valid]
    w = w[valid]
    if len(x) < 3:
        raise ValueError("weighted correlation requires at least three supported edges")
    sw = float(w.sum())
    if sw <= 0.0:
        raise ValueError("weighted correlation has zero total weight")
    mx = float(np.sum(w * x) / sw)
    my = float(np.sum(w * y) / sw)
    dx = x - mx
    dy = y - my
    vx = float(np.sum(w * dx * dx) / sw)
    vy = float(np.sum(w * dy * dy) / sw)
    if vx <= 0.0 or vy <= 0.0:
        raise ValueError("weighted correlation has zero variance")
    cov = float(np.sum(w * dx * dy) / sw)
    return float(np.clip(cov / np.sqrt(vx * vy), -1.0, 1.0))


def main() -> int:
    candidate_manifest = load_json(CANDIDATE_MANIFEST)
    h1 = load_json(H1_RESULT)
    h2 = load_json(H2_RESULT)
    if candidate_manifest.get("status") != "metadata_pool_frozen_before_candidate_image_pixels":
        raise RuntimeError("candidate pool is not the one-shot frozen metadata frame")
    if h1.get("status") != "complete_h1_evaluable":
        raise RuntimeError("completed evaluable H1 required")
    if h2.get("status") != "complete_h2_evaluable":
        raise RuntimeError("completed H2 required so this analysis is unambiguously post-H1/H2")
    if int(h1["primary"]["permutations"]) != 999:
        raise RuntimeError("expected exact 999 H1 null maps")
    if int(h1["primary"]["sampling_seed"]) != 20260903 or int(h1["primary"]["permutation_seed"]) != 20260904:
        raise RuntimeError("H1 seed lineage drifted")

    # Exact frozen-input identity checks.
    if sha256(CANDIDATE) != candidate_manifest["candidate_table_sha256"]:
        raise RuntimeError("candidate table SHA256 drifted")
    if sha256(H1_EDGES) != h1["file_sha256"]["edge_table"]:
        raise RuntimeError("H1 edge table SHA256 drifted")
    if sha256(H1_NULL_MAPS) != h1["file_sha256"]["null_persistence_maps"]:
        raise RuntimeError("H1 null-map SHA256 drifted")
    if h1["candidate_table_sha256"] != candidate_manifest["candidate_table_sha256"]:
        raise RuntimeError("H1 and candidate manifest disagree on frozen candidate identity")

    candidates = pd.read_csv(CANDIDATE, usecols=["cell_id", "species"])
    if len(candidates) != int(candidate_manifest["counts"]["observations"]):
        raise RuntimeError("candidate row count drifted")
    if candidates["species"].astype(str).nunique() != int(candidate_manifest["counts"]["species"]):
        raise RuntimeError("candidate species count drifted")

    counts = (
        candidates.groupby(["cell_id", "species"], observed=True)
        .size()
        .rename("candidate_n")
        .reset_index()
    )
    counts["capped_weight"] = counts["candidate_n"].clip(upper=SPECIES_CAP).astype(float)
    cell_species: dict[int, dict[str, float]] = {}
    cell_candidate_n: dict[int, int] = candidates.groupby("cell_id", observed=True).size().astype(int).to_dict()
    cell_species_n: dict[int, int] = counts.groupby("cell_id", observed=True).size().astype(int).to_dict()
    for cell_id, group in counts.groupby("cell_id", sort=False, observed=True):
        cell_species[int(cell_id)] = {
            str(species): float(weight)
            for species, weight in zip(group["species"], group["capped_weight"])
        }

    edges = pd.read_csv(H1_EDGES)
    required = {"edge_id", "cell_i", "cell_j", "opportunities", "persistence"}
    if not required.issubset(edges.columns):
        raise RuntimeError("H1 edge table lacks required fields")

    turnover: list[float] = []
    shared_species: list[int] = []
    union_species: list[int] = []
    for row in edges.itertuples(index=False):
        a = cell_species.get(int(row.cell_i), {})
        b = cell_species.get(int(row.cell_j), {})
        turnover.append(sparse_jsd_bits(a, b))
        sa = set(a)
        sb = set(b)
        shared_species.append(len(sa & sb))
        union_species.append(len(sa | sb))

    audit = edges.copy()
    audit["species_turnover_jsd_bits"] = turnover
    audit["candidate_n_i"] = audit["cell_i"].map(cell_candidate_n).fillna(0).astype(int)
    audit["candidate_n_j"] = audit["cell_j"].map(cell_candidate_n).fillna(0).astype(int)
    audit["species_n_i"] = audit["cell_i"].map(cell_species_n).fillna(0).astype(int)
    audit["species_n_j"] = audit["cell_j"].map(cell_species_n).fillna(0).astype(int)
    audit["shared_species"] = shared_species
    audit["union_species"] = union_species
    audit["species_jaccard_similarity"] = np.where(
        audit["union_species"] > 0,
        audit["shared_species"] / audit["union_species"],
        np.nan,
    )

    with np.load(H1_NULL_MAPS, allow_pickle=False) as archive:
        edge_ids = archive["edge_ids"].astype(str)
        null_persistence = archive["null_persistence"].astype(float)
    if null_persistence.shape != (999, len(audit)):
        raise RuntimeError("H1 null persistence map dimensions drifted")
    if not np.array_equal(edge_ids, audit["edge_id"].astype(str).to_numpy()):
        raise RuntimeError("H1 null edge order drifted")

    supported = (
        (audit["opportunities"].to_numpy(dtype=float) > 0)
        & np.isfinite(audit["persistence"].to_numpy(dtype=float))
        & np.isfinite(audit["species_turnover_jsd_bits"].to_numpy(dtype=float))
    )
    if int(np.count_nonzero(supported)) != int(h1["primary"]["supported_edges"]):
        raise RuntimeError("species-turnover audit does not cover exactly the H1-supported edge set")

    x = audit.loc[supported, "species_turnover_jsd_bits"].to_numpy(dtype=float)
    y = audit.loc[supported, "persistence"].to_numpy(dtype=float)
    w = audit.loc[supported, "opportunities"].to_numpy(dtype=float)
    observed_r = weighted_pearson(x, y, w)

    null_r = np.empty(999, dtype=float)
    for i in range(999):
        null_y = null_persistence[i, supported]
        if not np.isfinite(null_y).all():
            raise RuntimeError(f"non-finite H1 null persistence on supported edge set: permutation {i}")
        null_r[i] = weighted_pearson(x, null_y, w)

    p_upper = float((1 + np.count_nonzero(null_r >= observed_r)) / (len(null_r) + 1))
    p_two_sided_abs = float((1 + np.count_nonzero(np.abs(null_r) >= abs(observed_r))) / (len(null_r) + 1))
    percentile = float((np.count_nonzero(null_r < observed_r) + 0.5 * np.count_nonzero(null_r == observed_r)) / len(null_r))

    null_frame = pd.DataFrame({"permutation": np.arange(999, dtype=int), "weighted_r": null_r})
    OUT_NULL.parent.mkdir(parents=True, exist_ok=True)
    null_frame.to_csv(OUT_NULL, index=False, lineterminator="\n")
    audit.to_csv(OUT_EDGE, index=False, lineterminator="\n")

    # This is explicitly descriptive/exploratory. r^2 is not promoted to a
    # confirmatory causal fraction explained.
    interpretation = (
        "species_turnover_alignment_above_h1_null"
        if p_upper < 0.05 and observed_r > 0
        else "no_species_turnover_alignment_above_h1_null"
    )
    payload = {
        "protocol": "random-photo-first-h3-species-turnover-exploratory-v1",
        "status": "complete_exploratory_post_h1_h2_species_turnover_audit",
        "introduced_after_h1_h2_outcomes": True,
        "claim_role": "exploratory_decomposition_only_cannot_rescue_or_replace_h1_h2",
        "question": (
            "Do H1-supported edges with greater outcome-blind species-composition turnover show greater "
            "colour-boundary persistence than expected under the exact H1 null that preserves species geography?"
        ),
        "species_composition": {
            "source": "one-shot frozen candidate pool before image pixels",
            "candidate_rows": int(len(candidates)),
            "species": int(candidates["species"].astype(str).nunique()),
            "occupied_cells": int(candidates["cell_id"].nunique()),
            "cell_species_weight": "min(candidate_photo_count, 2)",
            "species_cap": SPECIES_CAP,
            "edge_metric": "Jensen-Shannon divergence in bits",
            "uses_colour_or_measurement_success": False,
        },
        "result": {
            "supported_edges": int(np.count_nonzero(supported)),
            "observed_weighted_r": observed_r,
            "observed_weighted_r_squared_descriptive": float(observed_r ** 2),
            "p_upper_against_exact_h1_null": p_upper,
            "p_two_sided_abs_against_exact_h1_null": p_two_sided_abs,
            "observed_r_percentile_in_h1_null": percentile,
            "null_r_mean": float(np.mean(null_r)),
            "null_r_median": float(np.median(null_r)),
            "null_r_q025": float(np.quantile(null_r, 0.025)),
            "null_r_q975": float(np.quantile(null_r, 0.975)),
            "species_turnover_jsd_mean": float(np.average(x, weights=w)),
            "species_turnover_jsd_median_unweighted": float(np.median(x)),
            "species_turnover_jsd_min": float(np.min(x)),
            "species_turnover_jsd_max": float(np.max(x)),
            "interpretation": interpretation,
        },
        "null_logic": (
            "The exact 999 H1 null persistence maps preserve each species' geography and species-specific morph "
            "marginals. Therefore the null correlation is the alignment expected from that preserved species "
            "structure; only excess positive alignment beyond this distribution would indicate additional "
            "turnover-associated colour geography."
        ),
        "frozen_primary_decisions_unchanged": {
            "h1": h1["decision"],
            "h1_p_upper": float(h1["primary"]["p_upper"]),
            "h2": h2["hierarchical_decision"],
        },
        "lineage": {
            "candidate_table_sha256": sha256(CANDIDATE),
            "h1_edge_table_sha256": sha256(H1_EDGES),
            "h1_null_maps_sha256": sha256(H1_NULL_MAPS),
            "h1_permutations": 999,
            "h1_sampling_seed": 20260903,
            "h1_permutation_seed": 20260904,
        },
        "files": {
            "edge_audit": str(OUT_EDGE.relative_to(ROOT)),
            "null_correlations": str(OUT_NULL.relative_to(ROOT)),
        },
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
